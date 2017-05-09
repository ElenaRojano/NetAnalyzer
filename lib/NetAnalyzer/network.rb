require 'nodes'
require 'nmatrix'
require 'pp'
require 'bigdecimal'

class Network
	attr_reader :association_values

	## BASIC METHODS
	############################################################
	def initialize(layers)
		@nodes = {} 
		@edges = {}
		@adjacency_matrices = {}
		@layers = layers
		@association_values = {}
		@control_connections = {}
	end

	def add_node(nodeID, nodeType = 0)
		@nodes[nodeID] = Node.new(nodeID, nodeType)
	end

	def add_edge(nodeID1, nodeID2)
		query_edge(nodeID1, nodeID2)
		query_edge(nodeID2, nodeID1)
	end

	def query_edge(nodeA, nodeB)
		query = @edges[nodeA]
		if query.nil?
			@edges[nodeA] = [nodeB]
		else
			query << nodeB
		end
	end

	def load_network_by_pairs(file, layers, split_character="\t")
		File.open(file).each("\n") do |line|
			line.chomp!
			pair = line.split(split_character)
			node1 = pair[0]
			node2 = pair[1]
			add_node(node1, set_layer(layers, node1))
			add_node(node2, set_layer(layers, node2))
			add_edge(node1, node2)	
		end
	end

	def get_edge_number
		node_connections = @edges.values.map{|connections| connections.length}.inject(0){|sum, n| sum + n}
		return node_connections/2
	end

	def plot(output_filename, layout="dot")		
		roboWrite = File.open(output_filename, 'w')
		roboWrite.puts "digraph g {"
		@edges.each do |nodeID, associatedIDs|
			associatedIDs.each do |associatedID|
				roboWrite.puts "\"#{nodeID}\"->\"#{associatedID}\";"
			end
		end
		roboWrite.puts "}"
		roboWrite.close
		cmd = "#{layout} -Tpng #{output_filename} -o #{output_filename}.png"
		system(cmd)
	end

	def get_all_intersections
		intersection_lengths = []
		get_all_pairs do |node1, node2|
			intersection_lengths << intersection(node1, node2).length
		end
		return intersection_lengths
	end

	def get_all_pairs(args = {})
		default = {:meth => :all, :layers => :all}
		args = default.merge(args)
		if args[:layers] == :all
			nodeIDs = @nodes.keys
		else
			nodeIDs = []
			args[:layers].each do |layer|
				nodeIDs.concat(@nodes.select{|id, node| node.type == layer}.keys)
			end
		end

		if args[:meth] == :all
			while !nodeIDs.empty?
				node1 = nodeIDs.shift
				nodeIDs.each do |node2|
					yield(node1, node2)
				end
			end
		#elsif args[:meth] == :conn
			
		end
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

	def generate_adjacency_matrix(layerA, layerB)
		layerAidNodes = @nodes.select{|id, node| node.type == layerA}.keys
		layerBidNodes = @nodes.select{|id, node| node.type == layerB}.keys
		adjacency_matrix = []
		layerAidNodes.each do |nodeA|
			layerBidNodes.each do |nodeB|
				if @edges[nodeB].include?(nodeA)
					adjacency_matrix << 1
				else
					adjacency_matrix << 0
				end
			end
		end
		matrix = NMatrix.new([layerAidNodes.length, layerBidNodes.length], adjacency_matrix)
		all_info_matrix = [matrix, layerAidNodes, layerBidNodes]
		@adjacency_matrices[[layerA, layerB]] = all_info_matrix
		return all_info_matrix
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
		if meth == :jaccard #all networks
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
		elsif meth == :csi #all networks
			relations = get_csi_associations(layers, base_layer)
		elsif meth == :transference #tripartite networks
			relations = get_association_by_transference_resources(layers, base_layer)
		end
		return relations
	end

	## association methods adjacency matrix based
	#---------------------------------------------------------
	# Alaimo 2014, doi: 10.3389/fbioe.2014.00071
	def get_association_by_transference_resources(firstPairLayers, secondPairLayers, lambda_value1 = 0.5, lambda_value2 = 0.5)
		matrix1 = @adjacency_matrices[firstPairLayers].first
		rowIds = @adjacency_matrices[firstPairLayers][1]
		matrix2 = @adjacency_matrices[secondPairLayers].first
		colIds =  @adjacency_matrices[secondPairLayers][2]
		m1rowNumber = matrix1.rows
		m1colNumber = matrix1.cols
		m2rowNumber = matrix2.rows
		m2colNumber = matrix2.cols
		#puts m1rowNumber, m1colNumber, m2rowNumber, m2colNumber
		matrix1Weight = graphWeights(m1colNumber, m1rowNumber, matrix1.transpose, lambda_value1)
		matrix2Weight = graphWeights(m2colNumber, m2rowNumber, matrix2.transpose, lambda_value2)
		matrixWeightProduct = matrix1Weight.dot(matrix2.dot(matrix2Weight))
		finalMatrix = matrix1.dot(matrixWeightProduct)
		relations = nmatrix2relations(finalMatrix, rowIds, colIds)
		@association_values[:transference] = relations
		return relations
	end

	## association methods node pairs based
	#---------------------------------------------------------
	# Bass 2013, doi:10.1038/nmeth.2728
	def get_associations(layers, base_layer) # BASE METHOD
		relations = []
		get_all_pairs(layers: layers) do |node1, node2|
			associatedIDs_node1 = @edges[node1].map{|id| @nodes[id]}.select{|node| node.type == base_layer}.map{|node| node.id}
			associatedIDs_node2 = @edges[node2].map{|id| @nodes[id]}.select{|node| node.type == base_layer}.map{|node| node.id}
			intersectedIDs = associatedIDs_node1 & associatedIDs_node2
			associationValue = yield(associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2)
			relations << [node1, node2, associationValue]  
		end
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
		ny = get_nodes_layer([base_layer]).length
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

	def get_hypergeometric_associations(layers, base_layer)
		ny = get_nodes_layer([base_layer]).length
		relations = get_associations(layers, base_layer) do |associatedIDs_node1, associatedIDs_node2, intersectedIDs, node1, node2|
			minLength = [associatedIDs_node1.length, associatedIDs_node2.length].min
			intersection_lengths = intersectedIDs.length
			sum = 0
			if intersection_lengths > 0
				nA = associatedIDs_node1.length
				nB = associatedIDs_node2.length
				#Using index from A layer proyected to B
				hyper_denom = binom(ny, nB)
				(intersection_lengths..minLength).each do |i|
					binom_product = binom(nA, i) * binom(ny - nA, nB - i)
					sum += binom_product.fdiv(hyper_denom)
				end
			end
			if sum == 0
				hypergeometricValue = 0
			else
				hypergeometricValue = -Math.log10(sum)
			end
			hypergeometricValue
		end
		@association_values[:hypergeometric] = relations
		return relations
	end

	def add_record(hash, node1, node2)
		query = hash[node1]
		if query.nil?
			hash[node1] = [node2]
		else
			query << node2
		end
	end


	def get_csi_associations(layers, base_layer)
		pcc_relations = get_pcc_associations(layers, base_layer)
		ny = get_nodes_layer([base_layer]).length
		pcc_vals = {}
		node_rels = {}
		pcc_relations.each do |node1, node2, assoc_index|
			pcc_vals[[node1, node2].sort] = assoc_index
			add_record(node_rels, node1, node2)
			add_record(node_rels, node2, node1)
		end
		relations = []
		pcc_relations.each do |node1, node2 ,assoc_index|
			pccAB = assoc_index - 0.05
			layer_nodes_conn2A = node_rels[node1]
			layer_nodes_conn2B = node_rels[node2]
			layer_intersectedIDs = layer_nodes_conn2B & layer_nodes_conn2A
			valid_connections = 0
			layer_intersectedIDs.each do |common_node|
				pccAcommon_value = pcc_vals[[node1, common_node].sort]
				pccBcommon_value = pcc_vals[[node2, common_node].sort]
				if pccBcommon_value < pccAB && pccAcommon_value < pccAB
					valid_connections += 1
				end
			end
			csiValue = valid_connections / ny.to_f

			relations << [node1, node2, csiValue]
		end

		@association_values[:csi] = relations
		return relations
	end


	## PERFORMANCE METHODS
	############################################################
	def load_control(ref_array)
		control = {}
		ref_array.each do |node1, node2|
			if node2 != '-'
				query = control[node1]
				if query.nil?
					control[node1] = [node2]
				else
					query << node2
				end
			end
		end
		@control_connections = control
		return control
	end

	def load_prediction(pairs_array)
		pred = {}
		min = nil
		max = nil
		pairs_array.each do |key, label, score|
			query = pred[key]
			if !min.nil? && !max.nil?
				min = score if score < min
				max = score if score > max
			else
				min = score; max = score
			end
			if query.nil?
				pred[key] = [[label], [score]]
			else
				query.first << label
				query.last << score
			end
		end
		return pred, [min, max]
	end


	# Pandey 2007, Association Analysis-based Transformations for Protein Interaction Networks: A Function Prediction Case Study
	def get_pred_rec(meth, cut_number = 100, top_number = 10000)
		performance = [] #cut, pred, rec
		preds, limits = load_prediction(@association_values[meth])
		cuts = get_cuts(limits, cut_number)
		cuts.each do |cut|
			prec, rec = pred_rec(preds, cut, top_number)
			performance << [cut, prec, rec]
		end
		return performance
	end

	def pred_rec(preds, cut, top)
		predicted_labels = 0 #m
		true_labels = 0 #n
		common_labels = 0 # k
		@control_connections.each do |key, c_labels|
			true_labels += c_labels.length #n
			pred_info = preds[key]
			if !pred_info.nil?
				labels, scores = pred_info
				reliable_labels = get_reliable_labels(labels, scores, cut, top)

				predicted_labels += reliable_labels.length #m
				common_labels += (c_labels & reliable_labels).length #k
			end
		end
		prec = common_labels.to_f/predicted_labels
		rec = common_labels.to_f/true_labels
		prec = 0.0 if prec.nan?
		rec = 0.0 if rec.nan?
		return prec, rec
	end



	## AUXILIAR METHODS
	#######################################################################################
	private

	def set_layer(layer_definitions, node_name)
		layer = nil
		if layer_definitions.length > 1
			layer_definitions.each do |layer_name, regexp|
				if node_name =~ regexp
					layer = layer_name
					break
				end
			end
		else
			layer = layer_definitions.first.first
		end
		return layer
	end

	def get_cuts(limits, n_cuts)
		cuts = []
		range = (limits.last - limits.first).to_f/n_cuts
		cut = limits.first
		n_cuts.times do
			cuts << cut
			cut += range
		end
		return cuts
	end

	def get_reliable_labels(labels, scores, cut, top)
		reliable_labels = []
		scores.each_with_index do |score, i|
			reliable_labels << [labels[i], score] if score >= cut
		end
		reliable_labels = reliable_labels.sort!{|l1,l2| l2.last <=> l1.last}[0..top-1].map{|pred| pred.first}
		return reliable_labels
	end
	
	def graphWeights (rowsNumber, colsNumber, inputMatrix, lambdaValue = 0.5)
		invMatrix = inputMatrix.sum(0).map{|e| 1.0/ e}
	 	diagonalColSums = NMatrix.diag(invMatrix)
	 	rowsSums = inputMatrix.sum(1).to_flat_a
	 	ky = NMatrix.new([rowsNumber, rowsNumber], rowsSums).map{|e| e ** lambdaValue } 	
	 	invertLambdaVal = (1 - lambdaValue)
	 	kx = NMatrix.new([rowsNumber, rowsNumber], rowsSums).transpose.map{|e| e ** invertLambdaVal } 
	 	nx = (ky * kx).map{|e| 1.0/ e}
	 	weigth = (inputMatrix.dot(diagonalColSums)).transpose
	 	weigth = inputMatrix.dot(weigth)
	 	weigth = nx * weigth
	 	return weigth
	end

	def nmatrix2relations(finalMatrix, rowIds, colIds)
		relations = []
		rowIds.each_with_index do |rowId, rowPos|
			colIds.each_with_index do |colId, colPos|
				associationValue = finalMatrix[rowPos, colPos]
				relations << [rowId, colId, associationValue]
			end
		end
		return relations
	end

	def binom(n,k)
		if k > 0 && k < n
  			res = (1+n-k..n).inject(:*)/(1..k).inject(:*)
  		else
  			res = 1
  		end
	end
end
