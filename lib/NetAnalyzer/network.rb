require 'rubystats'
#require 'nmatrix'
#require 'nmatrix/lapacke'
require 'numo/narray'
require 'numo/linalg'
require 'npy'
require 'parallel'

#require 'pp'
require 'benchmark'
#require 'nmatrix_expansion'

require 'semtools'
require 'expcalc'


class Network

	attr_accessor :association_values, :control_connections, :kernels, :reference_nodes, :group_nodes, :threads

	## BASIC METHODS
	############################################################
	def initialize(layers)
		@threads = 0
		@nodes = {}
		@edges = {}
		@reference_nodes = []
		@group_nodes = {}
		@adjacency_matrices = {}
		@kernels = {}
		@layers = layers
		@association_values = {}
		@control_connections = {}
		@compute_pairs = :conn
		@compute_autorelations = true
		@loaded_obos = []
		@ontologies = []
		@layer_ontologies = {}
	end

	def set_compute_pairs(use_pairs, get_autorelations)
		@compute_pairs = use_pairs
		@compute_autorelations = get_autorelations
	end

	def add_node(nodeID, nodeType = 0)
		@nodes[nodeID] = Node.new(nodeID, nodeType)
	end

	def add_edge(nodeID1, nodeID2)
		add_edge2hash(nodeID1, nodeID2)
		add_edge2hash(nodeID2, nodeID1)
	end

	def add_edge2hash(nodeA, nodeB)
		query = @edges[nodeA]
		if query.nil?
			@edges[nodeA] = [nodeB]
		else
			query << nodeB
		end
	end

	def set_layer(layer_definitions, node_name)
		layer = nil
		if layer_definitions.length > 1
			layer_definitions.each do |layer_name, regexp|
				if node_name =~ regexp
					layer = layer_name
					break
				end
			end
			raise("The node '#{node_name}' not match with any layer regex") if layer.nil?			
		else
			layer = layer_definitions.first.first
		end
		@layers << layer if !@layers.include?(layer)
		return layer
	end

	def generate_adjacency_matrix(layerA, layerB)
		layerAidNodes = @nodes.select{|id, node| node.type == layerA}.keys
		layerBidNodes = @nodes.select{|id, node| node.type == layerB}.keys
		matrix = Numo::DFloat.zeros(layerAidNodes.length, layerBidNodes.length)
		layerAidNodes.each_with_index do |nodeA, i|
			layerBidNodes.each_with_index do |nodeB, j|
				if @edges[nodeB].include?(nodeA)
					matrix[i, j] = 1
				else
					matrix[i, j] = 0
				end
			end
		end
		all_info_matrix = [matrix, layerAidNodes, layerBidNodes]

		if layerA == layerB
			@adjacency_matrices[[layerA]] = all_info_matrix
		else
			@adjacency_matrices[[layerA, layerB]] = all_info_matrix
		end
		return all_info_matrix
	end


	def delete_nodes(node_list, mode='d')
		if mode == 'd'
			@nodes.reject!{|n| node_list.include?(n)}
			@edges.reject!{|n, connections| node_list.include?(n)}
			@edges.each do |n, connections|
				connections.reject!{|c| node_list.include?(c)}
			end
		elsif mode == 'r'
			@nodes.select!{|n| node_list.include?(n)}
			@edges.select!{|n, connections| node_list.include?(n)}
			@edges.each do |n, connections|
				connections.select!{|c| node_list.include?(c)}
			end
		end
		@edges.reject!{|n, connections| connections.empty?}
	end

	def get_connected_nodes(node_id, from_layer)
		return @edges[node_id].map{|id| @nodes[id]}.select{|node| node.type == from_layer}.map{|node| node.id}
	end

	def get_nodes_from_layer(from_layer)
		return @nodes.values.select{|node| node.type == from_layer}.map{|node| node.id}
	end

	def get_bipartite_subgraph(from_layer_node_ids, from_layer, to_layer)
		bipartite_subgraph = {}
		from_layer_node_ids.each do |from_layer_node_id| 
			connected_nodes = @edges[from_layer_node_id]
			connected_nodes.each do |connected_node| 
				if @nodes[connected_node].type == to_layer
					query = bipartite_subgraph[connected_node]
					if query.nil?
						bipartite_subgraph[connected_node] = get_connected_nodes(connected_node, from_layer)
					end
				end
			end
		end
		return bipartite_subgraph
	end


	def get_edge_number
		node_connections = get_degree(zscore = false).values.inject(0){|sum, n| sum + n}
		return node_connections/2
	end

	def get_degree(zscore=true)
		degree = {}
		@edges.each do |id, nodes|
			degree[id] = nodes.length
		end
		if zscore
			degree_values = degree.values
			mean_degree = degree_values.mean
			std_degree = degree_values.standard_deviation
			degree.transform_values!{|v| (v - mean_degree).fdiv(std_degree)}
		end
		return degree
	end

	def get_node_attributes(attr_names)
		attrs = []
		attr_names.each do |attr_name|
			if attr_name == 'get_degree'
				attrs << get_degree
			elsif attr_name == 'get_degreeZ'
				attrs << get_degree(zscore=true)
			end
		end
		node_ids = attrs.first.keys
		node_attrs = []
		node_ids.each do |n|
			node_attrs << [n].concat(attrs.map{|at| at[n]})
		end
		return node_attrs
	end

	def plot_network(options = {})
		net_data = {
			group_nodes: @group_nodes,
			reference_nodes: @reference_nodes,
			nodes: @nodes,
			edges: @edges,
			layers: @layers
		}
		Net_plotter.new(net_data, options)
	end

	def compute_group_metrics(output_filename)
		metrics = []
		header = ['group']
		@group_nodes.keys.each do |k|
			metrics << [k]
		end
		header << 'comparative_degree'
		comparative_degree = communities_comparative_degree(@group_nodes)
		comparative_degree.each_with_index{|val,i| metrics[i] << replace_nil_vals(val)}
		header << 'avg_sht_path'
		avg_sht_path = communities_avg_sht_path(@group_nodes)
		avg_sht_path.each_with_index{|val,i| metrics[i] << replace_nil_vals(val)}
		if !@reference_nodes.empty?
			header.concat(%w[node_com_assoc_by_edge node_com_assoc_by_node])
			node_com_assoc = compute_node_com_assoc_in_precomputed_communities(@group_nodes, @reference_nodes.first)
			node_com_assoc.each_with_index{|val,i| metrics[i].concat(val)}
		end
		File.open(output_filename, 'w') do |f|
			f.puts header.join("\t")
			metrics.each do |gr|
				f. puts gr.join("\t")
			end
		end
	end

	def replace_nil_vals(val)
		return val.nil? ? 'NULL' : val
	end

	def communities_comparative_degree(coms) 
		comparative_degrees = []
		coms.each do |com_id, com|
			comparative_degrees << compute_comparative_degree(com)
		end
		return comparative_degrees
	end

	def communities_avg_sht_path(coms) 
		avg_sht_path = []
		coms.each do |com_id, com|
			dist, paths = compute_avg_sht_path(com)
			avg_sht_path << dist
		end
		return avg_sht_path
	end

	def compute_node_com_assoc_in_precomputed_communities(coms, ref_node)
		node_com_assoc = []
		coms.each do |com_id, com|
			node_com_assoc << [compute_node_com_assoc(com, ref_node)]
		end
		return node_com_assoc
	end

	def compute_comparative_degree(com) # see Girvan-Newman Benchmark control parameter in http://networksciencebook.com/chapter/9#testing (communities chapter)
		internal_degree = 0
		external_degree = 0
		com.each do |nodeID|
			nodeIDneigh = @edges[nodeID]
			next if nodeIDneigh.nil?
			internal_degree += (nodeIDneigh & com).length
			external_degree += (nodeIDneigh - com).length
		end
		comparative_degree = external_degree.fdiv(external_degree + internal_degree)
		return comparative_degree
	end

	def compute_avg_sht_path(com, paths=false)
		path_lengths = []
		all_paths = []
		group = com.dup
		while !group.empty?
			node_start = group.shift
			sht_paths = Parallel.map(group, in_processes: @threads) do |node_stop|
			#group.each do |node_stop|
				dist, path = shortest_path(node_start, node_stop, paths)
				[dist, path]
				#path_lengths << dist if !dist.nil?
				#all_paths << path if !path.empty?
			end
			sht_paths.each do |dist, path|
				path_lengths << dist
				all_paths << path				
			end
		end
		if path_lengths.include?(nil)
			avg_sht_path = nil
		else
			avg_sht_path = path_lengths.inject(0){|sum,l| sum + l}.fdiv(path_lengths.length)
		end
		return avg_sht_path, all_paths
	end

	# https://pythoninwonderland.wordpress.com/2017/03/18/how-to-implement-breadth-first-search-in-python/
	# finds shortest path between 2 nodes of a graph using BFS
	def bfs_shortest_path(start, goal, paths=false)
		dist = nil
	    explored = {} # keep track of explored nodes
	    previous = {}
	    queue = [[start, 0]] # keep track of all the paths to be checked
	    is_goal = false
	    while !queue.empty? && !is_goal # keeps looping until all possible paths have been checked
	        node, dist = queue.pop # pop the first path from the queue
	        if !explored.include?(node) # get the last node from the path
	            neighbours = @edges[node] 
	            explored[node] = true # mark node as explored
	            next if neighbours.nil?
	            dist += 1 
	            neighbours.each do |neighbour| # go through all neighbour nodes, construct a new path
	            	next if explored.include?(neighbour)
	                queue.unshift([neighbour, dist]) # push it into the queue
	                previous[neighbour] = node if paths
	                if neighbour == goal # return path if neighbour is goal
	                	is_goal = true
	                	break
	                end
	            end
	        end
	    end
		if is_goal
			path = build_path(previous, start, goal) if paths
		else
			dist = nil 
			path = []
		end
	    return dist, path
	end

	def build_path(previous, startNode, stopNode)
		path = []
		currentNode = stopNode
		path << currentNode
		while currentNode != startNode
		    currentNode = previous[currentNode]
			path << currentNode
		end
		return path
	end

	def shortest_path(node_start, node_stop, paths=false)
		#https://betterprogramming.pub/5-ways-to-find-the-shortest-path-in-a-graph-88cfefd0030f
		#return bidirectionalSearch(node_start, node_stop)
		#https://efficientcodeblog.wordpress.com/2017/12/13/bidirectional-search-two-end-bfs/
		dist, all_paths = bfs_shortest_path(node_start, node_stop, paths)
		return  dist, all_paths
	end

	def expand_clusters(expand_method)
		clusters = {}
		@group_nodes.each do |id, nodes|
			if expand_method == 'sht_path'
				dist, paths = compute_avg_sht_path(nodes, paths=true) # this uses bfs, maybe Dijkstra is the best one
				new_nodes = paths.flatten.uniq
				clusters[id] = nodes | new_nodes # If some node pair are not connected, recover them
			end
		end
		return clusters
	end

	def compute_node_com_assoc(com, ref_node)
		ref_cons = 0
		ref_secondary_cons = 0
		secondary_nodes = {}
		other_cons = 0
		other_nodes = {}

		refNneigh = @edges[ref_node]
		com.each do |nodeID|
			nodeIDneigh = @edges[nodeID]
			next if nodeIDneigh.nil?
			ref_cons += 1 if nodeIDneigh.include?(ref_node)
			if !refNneigh.nil?
				common_nodes = nodeIDneigh & refNneigh
				common_nodes.each {|id| secondary_nodes[id] = true}
				ref_secondary_cons += common_nodes.length 
			end
			specific_nodes = nodeIDneigh - refNneigh - [ref_node]
			specific_nodes.each {|id| other_nodes[id] = true}
			other_cons += specific_nodes.length
		end
		by_edge = (ref_cons + ref_secondary_cons).fdiv(other_cons)
		by_node = (ref_cons + secondary_nodes.length).fdiv(other_nodes.length)
		return by_edge, by_node
	end

	def get_all_intersections
		intersection_lengths = get_all_pairs do |node1, node2|
			intersection(node1, node2).length
		end
		return intersection_lengths
	end

	def get_all_pairs(args = {})
		all_pairs = []
		default = {:layers => :all}
		args = default.merge(args)
		nodeIDsA, nodeIDsB = collect_nodes(args)
		if @compute_autorelations
			if @compute_pairs == :all
				while !nodeIDsA.empty?
					node1 = nodeIDsA.shift
					pairs = Parallel.map(nodeIDsA, in_processes: @threads) do |node2|
						yield(node1, node2)
					end
					all_pairs.concat(pairs)
				end
			elsif @compute_pairs == :conn # TODO: Review this case to avoid return nil values
				while !nodeIDsA.empty?
					node1 = nodeIDsA.shift
					ids_connected_to_n1 = @edges[node1]
					pairs = Parallel.map(nodeIDsA, in_processes: @threads) do |node2|
						result = nil
						ids_connected_to_n2 = @edges[node2]
						if exist_connections?(ids_connected_to_n1, ids_connected_to_n2)
							result = yield(node1, node2)
						end
						result
					end
					pairs.compact!
					all_pairs.concat(pairs)
				end
			end
		else
			#MAIN METHOD
			if @compute_pairs == :conn
				all_pairs = Parallel.map(nodeIDsA, in_processes: @threads) do |node1|
					ids_connected_to_n1 = @edges[node1]
					node1_pairs = []
					nodeIDsB.each do |node2|
						ids_connected_to_n2 = @edges[node2]
						if exist_connections?(ids_connected_to_n1, ids_connected_to_n2)
							node1_pairs << yield(node1, node2)
						end
					end
					node1_pairs
				end
				all_pairs.flatten!(1)
			elsif @compute_pairs == :all
				raise 'Not implemented'
			end
		end

		return all_pairs
	end

	def collect_nodes(args)
		nodeIDsA = nil
		nodeIDsB = nil
		if @compute_autorelations
			if args[:layers] == :all
				nodeIDsA = @nodes.keys
			else
				nodeIDsA = []
				args[:layers].each do |layer|
					nodeIDsA.concat(@nodes.select{|id, node| node.type == layer}.keys)
				end
			end
		else
			if args[:layers] != :all
				nodeIDsA = @nodes.select{|id, node| node.type == args[:layers][0]}.keys
				nodeIDsB = @nodes.select{|id, node| node.type == args[:layers][1]}.keys
			end
		end
		return nodeIDsA, nodeIDsB
	end


	def get_nodes_layer(layers)
		#for creating ny value in hypergeometric and pcc index
		nodes = []
		layers.each do |layer|
			nodes.concat(@nodes.select{|nodeId, node| node.type == layer}.values)
		end
		return nodes
	end

	def intersection(node1, node2)
		shared_nodes = []
		associatedIDs_node1 = @edges[node1]
		associatedIDs_node2 = @edges[node2]
		intersectedIDs = associatedIDs_node1 & associatedIDs_node2
		intersectedIDs.each do |id|
			shared_nodes << @nodes[id]
		end
		return shared_nodes
	end

	def clean_autorelations_on_association_values
		@association_values.each do |meth, values|
			values.select!{|relation| @nodes[relation[0]].type != @nodes[relation[1]].type}
		end
	end

	## ASSOCIATION METHODS
	############################################################
	def get_association_values(layers, base_layer, meth)
		relations = [] #node A, node B, val
		if meth == :counts
			relations = get_counts_association(layers, base_layer)
		elsif meth == :jaccard #all networks
			relations = get_jaccard_association(layers, base_layer)
		elsif meth == :simpson #all networks
			relations = get_simpson_association(layers, base_layer)
		elsif meth == :geometric #all networks
			relations = get_geometric_associations(layers, base_layer)
		elsif meth == :cosine #all networks
			relations = get_cosine_associations(layers, base_layer)
		elsif meth == :pcc #all networks
			relations = get_pcc_associations(layers, base_layer)
		elsif meth == :hypergeometric #all networks
			relations = get_hypergeometric_associations(layers, base_layer)
		elsif meth == :hypergeometric_bf #all networks
			relations = get_hypergeometric_associations(layers, base_layer, :bonferroni)
		elsif meth == :hypergeometric_bh #all networks
			relations = get_hypergeometric_associations(layers, base_layer, :benjamini_hochberg)
		elsif meth == :hypergeometric_elim #tripartite networks?
			relations = get_hypergeometric_associations_with_topology(layers, base_layer, :elim)
		elsif meth == :hypergeometric_weight #tripartite networks?
			relations = get_hypergeometric_associations_with_topology(layers, base_layer, :weight)			
		elsif meth == :csi #all networks
			relations = get_csi_associations(layers, base_layer)
		elsif meth == :transference #tripartite networks
			relations = get_association_by_transference_resources(layers, base_layer)
		end
		return relations
	end

	## association methods adjacency matrix based
	#---------------------------------------------------------
	def get_association_by_transference_resources(firstPairLayers, secondPairLayers, lambda_value1 = 0.5, lambda_value2 = 0.5)
		relations = []
		matrix1 = @adjacency_matrices[firstPairLayers].first
		matrix2 = @adjacency_matrices[secondPairLayers].first
		finalMatrix = Adv_mat_calc.tranference_resources(matrix1, matrix2, lambda_value1 = lambda_value1, lambda_value2 = lambda_value2)
		rowIds = @adjacency_matrices[firstPairLayers][1]
		colIds =  @adjacency_matrices[secondPairLayers][2]
		relations = matrix2relations(finalMatrix, rowIds, colIds)
		@association_values[:transference] = relations
		return relations
	end

	## association methods node pairs based
	#---------------------------------------------------------
	# Bass 2013, doi:10.1038/nmeth.2728
	def get_associations(layers, base_layer) # BASE METHOD
		associations = get_all_pairs(layers: layers) do |node1, node2|
			associatedIDs_node1 = @edges[node1].map{|id| @nodes[id]}.select{|node| node.type == base_layer}.map{|node| node.id}
			associatedIDs_node2 = @edges[node2].map{|id| @nodes[id]}.select{|node| node.type == base_layer}.map{|node| node.id}
			intersectedIDs = associatedIDs_node1 & associatedIDs_node2
			associationValue = yield(associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2)
			[node1, node2, associationValue]  
		end
		return associations
	end

	def get_counts_association(layers, base_layer)
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|
			countValue = intersectedIDs.length	
		end
		@association_values[:counts] = relations
		return relations
	end

	def get_jaccard_association(layers, base_layer)
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|
			unionIDS = associatedIDs_node1 | associatedIDs_node2
			jaccValue = intersectedIDs.length.to_f/unionIDS.length		
		end
		@association_values[:jaccard] = relations
		return relations
	end

	def get_simpson_association(layers, base_layer)
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|
			minLength = [associatedIDs_node1.length, associatedIDs_node2.length].min
			simpsonValue = intersectedIDs.length.to_f/minLength
		end
		@association_values[:simpson] = relations
		return relations
	end

	def get_geometric_associations(layers, base_layer)
		#wang 2016 method
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|	
			intersectedIDs = intersectedIDs.length**2
			productLength = Math.sqrt(associatedIDs_node1.length * associatedIDs_node2.length)
			geometricValue = intersectedIDs.to_f/productLength
		end
		@association_values[:geometric] = relations
		return relations
	end

	def get_cosine_associations(layers, base_layer)
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|
			productLength = Math.sqrt(associatedIDs_node1.length * associatedIDs_node2.length)
			cosineValue = intersectedIDs.length/productLength
		end
		@association_values[:cosine] = relations
		return relations
	end

	def get_pcc_associations(layers, base_layer)
		#for Ny calcule use get_nodes_layer
		base_layer_nodes = get_nodes_layer([base_layer])
		ny = base_layer_nodes.length
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|
			intersProd = intersectedIDs.length * ny
			nodesProd = associatedIDs_node1.length * associatedIDs_node2.length
			nodesSubs = intersProd - nodesProd
			nodesAInNetwork = ny - associatedIDs_node1.length
			nodesBInNetwork = ny - associatedIDs_node2.length
			pccValue = nodesSubs.to_f / Math.sqrt(nodesProd * nodesAInNetwork * nodesBInNetwork)
		end
		@association_values[:pcc] = relations
		return relations
	end

	def get_hypergeometric_associations(layers, base_layer, pvalue_adj_method= nil)
		ny = get_nodes_layer([base_layer]).length
		fet = Rubystats::FishersExactTest.new
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|
			fisher = 0
			intersection_lengths = intersectedIDs.length
			if intersection_lengths > 0
				n1_items = associatedIDs_node1.length
				n2_items = associatedIDs_node2.length
				fisher = fet.calculate(
						intersection_lengths,
						n1_items - intersection_lengths,
						n2_items - intersection_lengths, 
						ny - (n1_items + n2_items - intersection_lengths)
				)
				fisher = fisher[:right]
			end
			fisher
		end
		if pvalue_adj_method == :bonferroni
			meth = :hypergeometric_bf
			compute_adjusted_pvalue_bonferroni(relations)
		elsif pvalue_adj_method == :benjamini_hochberg
			meth = :hypergeometric_bh
			compute_adjusted_pvalue_benjaminiHochberg(relations)
		else
			meth = :hypergeometric
			compute_log_transformation(relations)
		end
		@association_values[meth] = relations
		return relations
	end

	def get_hypergeometric_associations_with_topology(layers, base_layer, mode, thresold = 0.01)
		relations = []
		reference_layer = (layers - @layer_ontologies.keys).first
		ontology_layer = (layers - [reference_layer]).first
		ref_nodes = get_nodes_from_layer(reference_layer) # get nodes from NOT ontology layer
		ontology = @layer_ontologies[ontology_layer]
		base_layer_length = @nodes.values.count{|n| n.type == base_layer}
		ref_nodes.each do |ref_node|
			base_nodes = get_connected_nodes(ref_node, base_layer)
			ontology_base_subgraph = get_bipartite_subgraph(base_nodes, base_layer, ontology_layer) # get shared nodes between nodes from NOT ontology layer and ONTOLOGY layer. Also get the conections between shared nodes and ontology nodes.
			next if ontology_base_subgraph.empty?
			ontology_base_subgraph.transform_keys!{|k| k.to_sym}
			ontology.load_item_relations_to_terms(ontology_base_subgraph, remove_old_relations = true)
			term_pvals = ontology.compute_relations_to_items(base_nodes, base_layer_length, mode, thresold)
			relations.concat(term_pvals.map{|term| [ref_node, term[0], term[1]]})
		end
		compute_log_transformation(relations)
		if mode == :elim
			meth = :hypergeometric_elim
		elsif mode == :weight
			meth = :hypergeometric_weight
		end
		@association_values[meth] = relations
		return relations
	end

	def compute_adjusted_pvalue(relations, log_val=true)
		relations.each_with_index do |data, i| #p1, p2, pval
			pval_adj = yield(data.last, i)		
			pval_adj = -Math.log10(pval_adj) if log_val && pval_adj > 0
			data[2] = pval_adj 
		end
	end

	def compute_log_transformation(relations) #Only perform log transform whitout adjust pvalue. Called when adjusted method is not defined 
		compute_adjusted_pvalue(relations) do |pval, index| 
			pval
		end
	end

	def compute_adjusted_pvalue_bonferroni(relations)
		n_comparations = relations.length
		compute_adjusted_pvalue(relations) do |pval, index|
			adj = pval * n_comparations
			adj = 1 if adj > 1
			adj
		end
	end

	def compute_adjusted_pvalue_benjaminiHochberg(relations)
		adj_pvalues = get_benjaminiHochberg_pvalues(relations.map{|rel| rel.last})
		compute_adjusted_pvalue(relations) do |pval, index|
			adj_pvalues[index]
		end
	end

	def get_csi_associations(layers, base_layer)
		pcc_relations = get_pcc_associations(layers, base_layer)
		clean_autorelations_on_association_values if layers.length > 1
		nx = get_nodes_layer(layers).length
		pcc_vals = {}
		node_rels = {}
		pcc_relations.each do |node1, node2, assoc_index|
			add_nested_record(pcc_vals, node1, node2, assoc_index.abs)
			add_nested_record(pcc_vals, node2, node1, assoc_index.abs)
			add_record(node_rels, node1, node2)
			add_record(node_rels, node2, node1)
		end
		relations = []
		pcc_relations.each do |node1, node2 ,assoc_index|
			pccAB = assoc_index - 0.05
			valid_nodes = 0
			node_rels[node1].each do |node|
				valid_nodes += 1 if pcc_vals[node1][node] >= pccAB
			end
			node_rels[node2].each do |node|
				valid_nodes += 1 if pcc_vals[node2][node] >= pccAB
			end
			csiValue = 1 - (valid_nodes-1).fdiv(nx) 
			# valid_nodes-1 is done due to the connection node1-node2 is counted twice (one for each loop)
			relations << [node1, node2, csiValue]
		end
		@association_values[:csi] = relations
		return relations
	end

	def get_kernel(layer2kernel, kernel, normalization=false)
		matrix, node_names = @adjacency_matrices[layer2kernel]
		matrix_result = Adv_mat_calc.get_kernel(matrix, node_names, kernel, normalization=normalization)
		@kernels[layer2kernel] = matrix_result
	end

	def write_kernel(layer2kernel, output_file)
		@kernels[layer2kernel].save(output_file)
	end

	def link_ontology(ontology_file_path, layer_name)
		if !@loaded_obos.include?(ontology_file_path) #Load new ontology
			ontology = Ontology.new(file: ontology_file_path, load_file: true)
			@loaded_obos << ontology_file_path
			@ontologies << ontology
		else #Link loaded ontology to current layer
			ontology = @ontologies[@loaded_obos.index(ontology_file_path)]
		end
		@layer_ontologies[layer_name] = ontology
	end

	## RAMDOMIZATION METHODS
	############################################################
	def randomize_monopartite_net_by_nodes(layer)
		nodeIds = @adjacency_matrices[layer][1]
		nodeIds.shuffle!
		@adjacency_matrices[layer][1] = nodeIds
		@adjacency_matrices[layer][2] = nodeIds
	end	

	def randomize_bipartite_net_by_nodes(layers)
		layerA, layerB = layers
		rowIds = @adjacency_matrices[[layerA, layerB]][1]
		rowIds.shuffle!
		@adjacency_matrices[[layerA, layerB]][1] = rowIds
	end

	def randomize_monopartite_net_by_links(layers)
		nodesA = []
		nodesB = []
		relations = diagonal2relations(@adjacency_matrices[layers].first, @adjacency_matrices[layers][1], @adjacency_matrices[layers][2])
		relations.each do |relation|
			nodesA << relation[0]
			nodesB << relation[1]
		end
		nodesB.shuffle!
		@edges = {}
		nodesA.each do |nodeA|
			index_nodeB = 0
			while nodeA == nodesB[index_nodeB] 
				index_nodeB += 1
			end
			nodeB = nodesB.delete_at(index_nodeB)
			# if nodeB.nil? randomize_monopartite_net_by_links(layers)
			add_edge(nodeA, nodeB)
		end
		generate_adjacency_matrix(layers[0], layers[0])
	end


	def randomize_bipartite_net_by_links(layers) 
		nodesA = []
		nodesB = []
		relations = matrix2relations(@adjacency_matrices[layers].first, @adjacency_matrices[layers][1], @adjacency_matrices[layers][2])
		relations.each do |relation|
			nodesA << relation[0]
			nodesB << relation[1]
		end
		nodesB.shuffle!
		@edges = {}

		nodesA.each_with_index do |nodeA, i|
			add_edge(nodeA, nodesB[i])
		end
		generate_adjacency_matrix(layers[0], layers[1])

	end

	def save_adjacency_matrix(layerA, layerB, output_file)
		if layerA == layerB
			layers = [layerA]
		else 
			layers = [layerA, layerB]
		end
		Npy.save(output_file, @adjacency_matrices[layer].first)
		node_names = @nodes.values.map{|node| node.id}
		File.open(output_file+'.lst', 'w'){|f| f.print node_names.join("\n")}
	end

	def build_nodes_from_adjacency_matrix(layers_network, layers_adjacency_matrix)
		nodes_ids = @adjacency_matrices[layers_adjacency_matrix][1].concat(@adjacency_matrices[layers_adjacency_matrix][2]).uniq
		nodes_ids.each do |node_id|
			add_node(node_id, set_layer(layers_network, node_id))
		end
	end

	def build_edges_from_adjacency_matrix(layer)
		@edges = {}
		relations = matrix2relations(@adjacency_matrices[layer].first, @adjacency_matrices[layer][1], @adjacency_matrices[layer][2])
		relations.each do |relation|
			add_edge(relation[0], relation[1])
		end
	end


	## AUXILIAR METHODS
	#######################################################################################
	private

	def add_record(hash, node1, node2)
		query = hash[node1]
		if query.nil?
			hash[node1] = [node2]
		else
			query << node2
		end
	end

	def add_nested_record(hash, node1, node2, val)
		query_node1 = hash[node1]
		if query_node1.nil?
			hash[node1] = {node2 => val}
		else
			query_node1[node2] = val
		end
	end

 	def exist_connections?(ids_connected_to_n1, ids_connected_to_n2)
		res = false
		if !ids_connected_to_n1.nil? && 
			!ids_connected_to_n2.nil? && 
			!(ids_connected_to_n1 & ids_connected_to_n2).empty? # check that at least exists one node that connect to n1 and n2
			res = true
		end
		return res
	end

	def matrix2relations(finalMatrix, rowIds, colIds)
		relations = []
		rowIds.each_with_index do |rowId, rowPos|
			colIds.each_with_index do |colId, colPos|
				associationValue = finalMatrix[rowPos, colPos]
				relations << [rowId, colId, associationValue] if associationValue > 0
			end
		end
		return relations
	end

	def diagonal2relations(finalMatrix, rowIds, colIds)
		relations = []
		rowIds.each_with_index do |rowId, rowPos|
			colIds.each_with_index do |colId, colPos| 
				colMatrix = rowPos + colPos + 1
				if colMatrix < colIds.length
					associationValue = finalMatrix[rowPos, colMatrix]
					relations << [rowId, colIds[colMatrix], associationValue] if associationValue > 0
				end
			end
		end
		return relations
	end

end
