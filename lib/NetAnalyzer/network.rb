require 'nodes'
require 'nmatrix'
require 'nmatrix/lapacke'
#require 'pp'
require 'bigdecimal'
require 'benchmark'
require 'nmatrix_expansion'

class Network

	attr_accessor :association_values, :control_connections, :kernels

	## BASIC METHODS
	############################################################
	def initialize(layers)
		@nodes = {} 
		@edges = {}
		@adjacency_matrices = {}
		@kernels = {}
		@layers = layers
		@association_values = {}
		@control_connections = {}
		@compute_pairs = :conn
		@compute_autorelations = true
		@matrix_byte_format = :float64
	end

	def set_compute_pairs(use_pairs, get_autorelations)
		@compute_pairs = use_pairs
		@compute_autorelations = get_autorelations
	end

	def set_matrix_byte_format(matrix_format)
		@matrix_byte_format = matrix_format
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
		#File.open(file).each("\n") do |line| # old bug from io patching from lapacke gem
		File.open(file).each do |line|
			line.chomp!
			pair = line.split(split_character)
			node1 = pair[0]
			node2 = pair[1]
			add_node(node1, set_layer(layers, node1))
			add_node(node2, set_layer(layers, node2))
			add_edge(node1, node2)	
		end
	end

	def load_network_by_bin_matrix(input_file, node_file, layers)
		node_names = load_input_list(node_file)
		@adjacency_matrices[layers.map{|l| l.first}] = [NMatrix.read(input_file), node_names, node_names]
	end

	def load_network_by_plain_matrix(input_file, node_file, layers, splitChar)
		node_names = load_input_list(node_file)
		@adjacency_matrices[layers.map{|l| l.first}] = [load_matrix_file(input_file, splitChar), node_names, node_names]
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
		default = {:layers => :all}
		args = default.merge(args)
		nodeIDsA, nodeIDsB = collect_nodes(args)
		if @compute_autorelations
			if @compute_pairs == :all
				while !nodeIDsA.empty?
					node1 = nodeIDsA.shift
					nodeIDsA.each do |node2|
						yield(node1, node2)
					end
				end
			elsif @compute_pairs == :conn
				processed_node_ids = {}
				while !nodeIDsA.empty?
					node1 = nodeIDsA.shift
					ids_connected_to_n1 = @edges[node1]
					nodeIDsA.each do |node2|
						if processed_node_ids[node2].nil?
							ids_connected_to_n2 = @edges[node2]
							if exist_connections?(ids_connected_to_n1, ids_connected_to_n2)
								yield(node1, node2)
							end
						end
					end
					processed_node_ids[node1] = true
				end
			end
		else
			if @compute_pairs == :conn
				processed_node_ids = {}
				nodeIDsA.each do |node1|
					ids_connected_to_n1 = @edges[node1]
					nodeIDsB.each do |node2|
						if processed_node_ids[node2].nil?
							ids_connected_to_n2 = @edges[node2]
							if exist_connections?(ids_connected_to_n1, ids_connected_to_n2)
								yield(node1, node2)
							end
						end
					end
					processed_node_ids[node1] = true
				end
			end

		end
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

	def add_nested_record(hash, node1, node2, val)
		query_node1 = hash[node1]
		if query_node1.nil?
			hash[node1] = {node2 => val}
		else
			query_node1[node2] = val
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
		#puts "cut: #{cut} trueL: #{true_labels} predL: #{predicted_labels} commL: #{common_labels}"
		prec = common_labels.to_f/predicted_labels
		rec = common_labels.to_f/true_labels
		prec = 0.0 if prec.nan?
		rec = 0.0 if rec.nan?
		return prec, rec
	end

	## KERNEL METHODS
	#######################################################################################
	def get_kernel(layer2kernel, kernel, normalization=false)
		#matrix = NMatrix.new([3, 3],[1, 1, 0, 0, 0, 2, 0, 5, -1], dtype: @matrix_byte_format)
		matrix, node_names = @adjacency_matrices[layer2kernel]
		#I = identity matrix
		#D = Diagonal matrix
		#A = adjacency matrix
		#L = laplacian matrix = D − A
		matrix_result = nil
		dimension_elements = matrix.cols
		# In scuba code, the diagonal values of A is set to 0. In weighted matrix the kernel result is the same with or without this operation. Maybe increases the computing performance?
		# In the md kernel this operation affects the values of the final kernel
		#dimension_elements.times do |n|
		#	matrix[n,n] = 0.0
		#end
		if kernel == 'el' || kernel == 'ct' || kernel == 'rf' || 
			kernel.include?('vn') || kernel.include?('rl') || kernel == 'me'
			row_sum = matrix.sum(1) # get the total sum for each row, for this reason the sum method takes the 1 value. If sum colums is desired, use 0
			diagonal_matrix = NMatrix.diag(row_sum, dtype: @matrix_byte_format) # Make a matrix whose diagonal is row_sum
			matrix_L = diagonal_matrix - matrix
			if kernel == 'el' #Exponential Laplacian diffusion kernel(active). F Fouss 2012 | doi: 10.1016/j.neunet.2012.03.001
			    beta = 0.02
			    beta_product = matrix_L * beta
			    matrix_result = beta_product.expm
			elsif kernel == 'ct' # Commute time kernel (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			    matrix_result = matrix_L.pinv # Anibal saids that this kernel was normalized. Why?. Paper do not seem to describe this operation for ct, it describes for Kvn or for all kernels, it is not clear.
			elsif kernel == 'rf' # Random forest kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			    matrix_result = (NMatrix.eye(dimension_elements, dtype: @matrix_byte_format) + matrix_L).invert! #Krf = (I +L ) ^ −1
			elsif kernel.include?('vn') # von Neumann diffusion kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			    alpha = kernel.gsub('vn', '').to_f * matrix.max_eigenvalue ** -1  # alpha = impact_of_penalization (1, 0.5 or 0.1) * spectral radius of A. spectral radius of A = absolute value of max eigenvalue of A 
			    matrix_result = (NMatrix.eye(dimension_elements, dtype: @matrix_byte_format) - matrix * alpha ).invert! #  (I -alphaA ) ^ −1
			elsif kernel.include?('rl') # Regularized Laplacian kernel matrix (active)
			    alpha = kernel.gsub('rl', '').to_f * matrix.max_eigenvalue ** -1  # alpha = impact_of_penalization (1, 0.5 or 0.1) * spectral radius of A. spectral radius of A = absolute value of max eigenvalue of A
			    matrix_result = (NMatrix.eye(dimension_elements, dtype: @matrix_byte_format) + matrix_L * alpha ).invert! #  (I + alphaL ) ^ −1
			elsif kernel == 'me' # Markov exponential diffusion kernel (active). G Zampieri 2018 | doi.org/10.1186/s12859-018-2025-5 . Taken from compute_kernel script
				beta=0.04
				#(beta/N)*(N*I - D + A)
				id_mat = NMatrix.eye(dimension_elements, dtype: @matrix_byte_format)
				m_matrix = (id_mat * dimension_elements - diagonal_matrix + matrix ) * (beta/dimension_elements)
				matrix_result = expm(m_matrix)
			end
		elsif kernel == 'ka' # Kernelized adjacency matrix (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			lambda_value = matrix.min_eigenvalue
			matrix_result = matrix + NMatrix.eye(dimension_elements, dtype: @matrix_byte_format) * lambda_value.abs # Ka = A + lambda*I # lambda = the absolute value of the smallest eigenvalue of A
		elsif kernel.include?('md') # Markov diffusion kernel matrix. G Zampieri 2018 | doi.org/10.1186/s12859-018-2025-5 . Taken from compute_kernel script
			t = kernel.gsub('md', '').to_i
			col_sum = matrix.sum(1)
			p_mat = matrix.div_by_vector(col_sum)
			p_temp_mat = p_mat.clone
			zt_mat = p_mat.clone
			(t-1).times do
				p_temp_mat = p_temp_mat.dot(p_mat)
				zt_mat = zt_mat + p_temp_mat
			end
			zt_mat = zt_mat * (1.0/t)
			matrix_result = zt_mat.dot(zt_mat.transpose())
		else
			matrix_result = matrix
			warn('Warning: The kernel method was not specified or not exists. The adjacency matrix will be given as result')
			# This allows process a previous kernel and perform the normalization in a separated step.
		end
		matrix_result = matrix_result.cosine_normalization if normalization
		@kernels[layer2kernel] = matrix_result
	end

	def write_kernel(layer2kernel, output_file)
		@kernels[layer2kernel].write(output_file)
	end

	## AUXILIAR METHODS
	#######################################################################################
	private

	def load_input_list(file)
		return File.open(file).readlines.map!{|line| line.chomp}
	end

	def load_matrix_file(input_file, splitChar = "\t")
		dimension_elements = 0
		adjacency_vector = []
		File.open(input_file).each do |line|
		    	line.chomp!
	    		adjacency_vector.concat(line.split(splitChar).map{|c| c.to_f })
		    	dimension_elements += 1
		end
		matrix = NMatrix.new([dimension_elements, dimension_elements], adjacency_vector, @matrix_byte_format) # Create working matrix
		return matrix
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
		range = (limits.last - limits.first).abs.fdiv(n_cuts)
		range = BigDecimal(range, 10)
		cut = limits.first
		(n_cuts + 1).times do |n|
			cuts << (cut + n * range).to_f
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
	 	diagonalColSums = NMatrix.diag(invMatrix, @matrix_byte_format)
	 	rowsSums = inputMatrix.sum(1).to_flat_a
	 	ky = NMatrix.new([rowsNumber, rowsNumber], rowsSums, @matrix_byte_format).map{|e| e ** lambdaValue } 	
	 	invertLambdaVal = (1 - lambdaValue)
	 	kx = NMatrix.new([rowsNumber, rowsNumber], rowsSums, @matrix_byte_format).transpose.map{|e| e ** invertLambdaVal } 
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
