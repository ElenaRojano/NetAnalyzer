ROOT_PATH = File.dirname(__FILE__)
require File.join(ROOT_PATH, 'test_helper.rb')




class NodeTest < Minitest::Test
	def setup 
		@node_to_test = Node.new("M1", :main)
	end

	def test_node_clone
		node_clone = @node_to_test.clone
		assert_equal node_clone, @node_to_test
	end
end

class NetworkTest < Minitest::Test

	def setup
		@tripartite_layers = [[:main, /M[0-9]+/], [:projection, /P[0-9]+/], [:salient, /S[0-9]+/]]
		@tripartite_network = Net_parser.load_network_by_pairs(File.join(DATA_TEST_PATH, 'tripartite_network_for_validating.txt'), @tripartite_layers)

		@bipartite_layers = [[:main, /M[0-9]+/], [:projection, /P[0-9]+/]]
		@network_obj = Net_parser.load_network_by_pairs(File.join(DATA_TEST_PATH, 'bipartite_network_for_validating.txt'), @bipartite_layers)
		@network_obj.generate_adjacency_matrix(@bipartite_layers[0].first, @bipartite_layers[1].first)
		@monopartite_layers = [[:main, /\w/], [:main, /\w/]]
		@monopartite_network = Net_parser.load_network_by_pairs(File.join(DATA_TEST_PATH, 'monopartite_network_for_validating.txt'), @monopartite_layers)
		@monopartite_network.generate_adjacency_matrix(@monopartite_layers[0].first, @monopartite_layers[0].first)
	end

	def test_clone
		network_clone = @network_obj.clone
		assert_equal @network_obj, network_clone	
		network_clone.compute_autorelations = false
	end

	def test_clone_change
		network_clone = @network_obj.clone
		network_clone.add_node('M8', network_clone.set_layer(@bipartite_layers, 'M8'))
		refute_equal @network_obj, network_clone	
	end

	def test_add_node
		network_clone = @network_obj.clone
		new_node = 'M8'
		node_test = Node.new(new_node, :main)
		network_clone.add_node(new_node, :main)
		test_add_node = network_clone.nodes[new_node]
		
		assert_equal node_test, test_add_node
	end

	def test_add_edge
		network_clone = @network_obj.clone
		node_1 = 'M8'
		node_2 = 'M9'
		network_clone.add_node(node_1, network_clone.set_layer(@bipartite_layers, node_1))
		network_clone.add_node(node_2, network_clone.set_layer(@bipartite_layers, node_2))
		network_clone.add_edge(node_1, node_2)
		test_node_1 = network_clone.get_connected_nodes('M9', :main).include?('M8')
		assert_equal true, test_node_1
		test_node_1 = network_clone.get_connected_nodes('M8', :main).include?('M9')
		assert_equal true, test_node_1
	end

	def test_add_edge2hash
		network_clone = @network_obj.clone
		node_1 = 'M8'
		node_2 = 'M9'
		network_clone.add_node(node_2, network_clone.set_layer(@bipartite_layers, node_2))
		network_clone.add_edge2hash(node_1, node_2)
		result_test = network_clone.get_connected_nodes('M8', :main).include?('M9')
		assert_equal true, test_add_edge
	end

	def test_set_layer
		network_clone = @network_obj.clone
		node_name = "M8"
		layer_test = network_clone.set_layer(@bipartite_layers, node_name)
		expected_result = :main
		assert_equal expected_result, layer_test
	end

	def test_set_layer_excep
		network_clone = @network_obj.clone
		node_name = "Z8"
		begin
			network_clone.set_layer(@bipartite_layers, node_name)
		rescue Exception => e
			test_result = e.message
		end
		expected_result = "The node '#{node_name}' not match with any layer regex"
		assert_equal expected_result, test_result
	end

	def test_generate_adjacency_matrix_monopartite
		test_result = @monopartite_network.adjacency_matrices
		matrix_values = Numo::DFloat[[0, 1, 1, 0, 0],[1, 0, 0, 0, 0], [1, 0, 0, 0, 1], [0, 0, 0, 0, 1], [0, 0, 1, 1, 0]]
		expected_result = {[:main] => [matrix_values, ['A', 'C', 'E', 'B', 'D'], ['A', 'C', 'E', 'B', 'D']]} 
		assert_equal expected_result, test_result	
	end

	def test_generate_adjacency_matrix_bipartite
		test_result = @network_obj.adjacency_matrices
		matrix_values = Numo::DFloat[[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 1, 1, 0, 0, 0, 0],[1, 1, 1, 1, 0, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]]
		expected_result = {[:main, :projection] => [matrix_values, ['M1', 'M2', 'M3', 'M4', 'M5', 'M6'], ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10']]} 
		assert_equal expected_result, test_result	
	end

	def test_delete_nodes_d_mono
		network_clone = @monopartite_network.clone
		network_clone.delete_nodes(['E'])
		nodes_test_result = network_clone.get_nodes_layer([:main]).length
		assert_equal 4, nodes_test_result
		edges_test_result = network_clone.get_edge_number
		assert_equal 2, edges_test_result
	end

	def test_delete_nodes_d_bi
		network_clone = @network_obj.clone
		network_clone.delete_nodes(['M1', 'M2'])
		nodes_test_result = network_clone.get_nodes_layer([:main]).length
		assert_equal 4, nodes_test_result
		edges_test_result = network_clone.get_edge_number
		assert_equal 20, edges_test_result
	end

	def test_delete_nodes_r_mono
		network_clone = @monopartite_network.clone
		network_clone.delete_nodes(['E'])
		nodes_test_result = network_clone.get_nodes_layer([:main]).length
		assert_equal 4, nodes_test_result
		edges_test_result = network_clone.get_edge_number
		assert_equal 2, edges_test_result
	end

	def test_delete_nodes_r_bi
		network_clone = @network_obj.clone
		network_clone.delete_nodes(['M1', 'M2'])
		nodes_test_result = network_clone.get_nodes_layer([:main]).length
		assert_equal 4, nodes_test_result
		edges_test_result = network_clone.get_edge_number
		assert_equal 20, edges_test_result
	end

	def test_get_connected_nodes
		test_result = @monopartite_network.get_connected_nodes('A', :main)
		expected_result = ['C', 'E']
		assert_equal expected_result, test_result
	end

	def test_get_nodes_from_layer
		test_result = @network_obj.get_nodes_from_layer(:main)
		expected_result = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6']
		assert_equal expected_result, test_result
	end

	def test_get_bipartite_subgraph
		bipartirte_test = @tripartite_network.get_bipartite_subgraph(['M1', 'M2', 'M3', 'M4', 'M5', 'M6'], :salient, :projection)
		expected_result = {'P1' =>["S2"], "P2"=>["S2"], "P3"=>["S1", "S6"], "P4"=>["S2"], "P5"=>["S1", "S4"], "P6"=>[], "P7"=>["S3"], "P8"=>["S5"], "P9"=>["S6"], "P10"=>["S4", "S6"]}
		assert_equal expected_result, bipartirte_test
	end

	def test_get_edge_number
		edge_number_test = @monopartite_network.get_edge_number
		assert_equal 4, edge_number_test
	end

	def test_get_degree_zscore
		degree_test = @monopartite_network.get_degree
		expected_result = {'A' => 0.8164965809277259, 'C' => -1.2247448713915894, 'E' => 0.8164965809277259, 'B' => -1.2247448713915894, 'D' => 0.8164965809277259}
		assert_equal expected_result, degree_test
	end

	def test_get_degree
		degree_test = @monopartite_network.get_degree(zscore=false)
		expected_result = {'A' => 2, 'C' => 1, 'E' => 2, 'B' => 1, 'D' => 2}
		assert_equal expected_result, degree_test
	end

	def test_get_all_intersections_autorr_all_layers_conn
		network_clone = @tripartite_network.clone
		test_result = network_clone.get_all_intersections
		expected_result = [10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 7, 5, 6, 4, 4, 3, 3, 2, 2, 5, 6, 4, 4, 3, 3, 2, 2, 
			5, 5, 4, 3, 3, 3, 3, 4, 4, 3, 3, 2, 2, 4, 3, 3, 2, 3, 3, 3, 2, 2, 3, 2, 2, 2, 2, 3, 8, 6, 4, 2, 2,
			3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 1, 1, 1, 1, 4, 2, 2, 3, 1, 1, 2, 1, 3, 1, 2, 1, 1, 1]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_autorr_all_layers_all
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:all, true)
		test_result = network_clone.get_all_intersections
		expected_result = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 7, 5, 6, 4, 4, 3, 3, 2, 2, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 6, 4, 4, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 4, 3, 3, 3, 3, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 3, 2, 3, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 1, 1, 1, 1, 4, 2, 2, 3, 0, 1, 0, 1, 2, 1, 3, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_autorr_some_layers_conn
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:conn, true)
		test_result = network_clone.get_all_intersections({:layers =>[:main, :salient]})
		expected_result = [10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 
			1, 1, 1, 1, 4, 2, 2, 3, 1, 1, 2, 1, 3, 1, 2, 1, 1, 1]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_autorr_some_layers_all
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:all, true)
		test_result = network_clone.get_all_intersections({:layers =>[:main, :salient]})
		expected_result = [10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 
			1, 1, 1, 1, 4, 2, 2, 3, 0, 1, 0, 1, 2, 1, 3, 0, 0, 0, 1, 0, 2, 
			0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_no_autorr_all_layers_conn
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:conn, false)
		test_result = network_clone.get_all_intersections
		expected_result = []
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_no_autorr_some_layers_conn
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:conn, false)
		test_result = network_clone.get_all_intersections({:layers =>[:main, :salient]})
		expected_result = [2, 3, 1, 2, 1, 3, 2, 3, 1, 2, 1, 3, 2, 3, 1, 1, 1, 1, 2, 3, 1, 1, 1, 3, 1, 2]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_no_autorr_all
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:all, false)
		begin
			network_clone.get_all_intersections()
		rescue Exception => e
			test_result = e.message
		end
		expected_result = 'Not implemented'
		assert_equal expected_result, test_result		
	end

	def test_get_all_intersections_autorr_all_layers_conn
		network_clone = @tripartite_network.clone
		test_result = network_clone.get_all_pairs() do |node1, node2|
			network_clone.intersection(node1, node2).length
		end
		expected_result = [10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 7, 5, 6, 4, 4, 3, 3, 2, 2, 5, 6, 4, 4, 3, 3, 2, 2, 
			5, 5, 4, 3, 3, 3, 3, 4, 4, 3, 3, 2, 2, 4, 3, 3, 2, 3, 3, 3, 2, 2, 3, 2, 2, 2, 2, 3, 8, 6, 4, 2, 2,
			3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 1, 1, 1, 1, 4, 2, 2, 3, 1, 1, 2, 1, 3, 1, 2, 1, 1, 1]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_autorr_all_layers_all
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:all, true)
		test_result = network_clone.get_all_pairs() do |node1, node2|
			network_clone.intersection(node1, node2).length
		end
		expected_result = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 7, 5, 6, 4, 4, 3, 3, 2, 2, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 6, 4, 4, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 4, 3, 3, 3, 3, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 3, 2, 3, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 1, 1, 1, 1, 4, 2, 2, 3, 0, 1, 0, 1, 2, 1, 3, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_autorr_some_layers_conn
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:conn, true)
		test_result = network_clone.get_all_pairs({:layers =>[:main, :salient]}) do |node1, node2|
			network_clone.intersection(node1, node2).length
		end
		expected_result = [10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 
			1, 1, 1, 1, 4, 2, 2, 3, 1, 1, 2, 1, 3, 1, 2, 1, 1, 1]
		assert_equal expected_result, test_result
	end

	def test_get_all_intersections_autorr_some_layers_all
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:all, true)
		test_result = network_clone.get_all_pairs({:layers =>[:main, :salient]}) do |node1, node2|
			network_clone.intersection(node1, node2).length
		end
		expected_result = [10, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 8, 6, 4, 2, 2, 3, 1, 2, 1, 3, 6, 4, 2, 2, 3, 
			1, 1, 1, 1, 4, 2, 2, 3, 0, 1, 0, 1, 2, 1, 3, 0, 0, 0, 1, 0, 2, 
			0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
		assert_equal expected_result, test_result
	end

	def test_get_all_pairs_no_autorr_all_layers_conn
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:conn, false)
		test_result = network_clone.get_all_pairs
		expected_result = []
		assert_equal expected_result, test_result
	end

	def test_get_all_pairs_no_autorr_some_layers_conn
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:conn, false)
		test_result = network_clone.get_all_pairs({:layers =>[:main, :salient]}) do |node1, node2|
			network_clone.intersection(node1, node2).length
		end
		expected_result = [2, 3, 1, 2, 1, 3, 2, 3, 1, 2, 1, 3, 2, 3, 1, 1, 1, 1, 2, 3, 1, 1, 1, 3, 1, 2]
		assert_equal expected_result, test_result
	end

	def test_get_all_pairs_no_autorr_all
		network_clone = @network_obj.clone
		network_clone.set_compute_pairs(:all, false)
		begin
			network_clone.get_all_pairs()
		rescue Exception => e
			test_result = e.message
		end
		expected_result = 'Not implemented'
		assert_equal expected_result, test_result		
	end

	def test_collect_nodes_autorr_some_layers
		nodesA_test, nodesB_test = @tripartite_network.collect_nodes({:layers =>[:main, :salient]})
		expected_result_nodesA = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6']
		assert_equal expected_result_nodesA, nodesA_test
		assert_nil nodesB_test

	end

	def test_collect_nodes_autorr_all_layers
		nodesA_test, nodesB_test = @network_obj.collect_nodes({:layers => :all})
		expected_result_nodesA = ['M1', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'M2', 'M3', 'M4', 'M5', 'M6']
		assert_equal expected_result_nodesA, nodesA_test
		assert_nil nodesB_test
	end

	def test_collect_nodes_no_autorr_some_layers
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:all, false)
		nodesA_test, nodesB_test = network_clone.collect_nodes({:layers =>[:main, :salient]})
		expected_result_nodesA = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6']
		expected_result_nodesB = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
		assert_equal expected_result_nodesA, nodesA_test
		assert_equal expected_result_nodesB, nodesB_test
	end

	def test_collect_nodes_no_autorr_all_layers
		network_clone = @tripartite_network.clone
		network_clone.set_compute_pairs(:conn, false)
		nodesA_test, nodesB_test = network_clone.collect_nodes({:layers => :all})
		assert_nil nodesA_test
		assert_nil nodesB_test
	end

	def test_get_nodes_layer
		nodes_from_layers_test = @tripartite_network.get_nodes_layer([:main, :salient]).map{|node| node.id}
		expected_result = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 
			'S1', 'S2', 'S3', 'S4', 'S5', 'S6']
		assert_equal expected_result, nodes_from_layers_test
	end

	def test_intersection
		test_result = @network_obj.intersection('M3', 'M6').map{|node| node.id}
		expected_result = ['P1', 'P2']
		assert_equal expected_result, test_result
	end

	def test_get_node_attributes
		node_attribute_test = @monopartite_network.get_node_attributes(['get_degree'])
		expected_result = [['A', 2], ['C', 1], ['E', 2],['B', 1], ['D', 2]]
		assert_equal expected_result, node_attribute_test
	end

	def test_get_node_attributes_zscore
		node_attribute_test = @monopartite_network.get_node_attributes(['get_degreeZ', 'get_degree'])
		expected_result = [['A', 0.8164965809277259, 2], ['C', -1.2247448713915894, 1], ['E', 0.8164965809277259, 2],['B', -1.2247448713915894, 1], ['D', 0.8164965809277259, 2]]
		node_attribute_test = @monopartite_network.get_node_attributes(['get_degreeZ', "get_degree"])
	end

	def test_get_counts_association
		test_association = @network_obj.get_counts_association([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2]]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'counts_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end

	def test_get_jaccard_association
		test_association = @network_obj.get_jaccard_association([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'jaccard_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value.round(6)]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end

	def test_get_simpson_association
		test_association = @network_obj.get_simpson_association([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'simpson_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value.round(6)]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end

	def test_get_geometric_associations
		test_association = @network_obj.get_geometric_associations([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'geometric_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value.round(6)]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end


	def test_get_cosine_associations
		test_association = @network_obj.get_cosine_associations([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'cosine_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value.round(6)]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end

=begin
	def test_get_pcc_associations
		test_association = @network_obj.get_pcc_associations([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'pcc_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value.round(6)]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end

	def test_get_hypergeometric_associations
		test_association = @network_obj.get_hypergeometric_associations([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'hyi_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value.round(6)]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end

	def test_get_csi_associations
		test_association = @network_obj.get_csi_associations([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
		all_association_values = []
		File.open(File.join(ROOT_PATH, 'csi_results.txt')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			association_value = fields.pop.to_f
			all_association_values << [fields[0], fields[1], association_value.round(6)]
		end
		assert_equal(all_association_values.sort, test_association.sort)
	end
=end
	#def test_randomize_monopartite_net_by_nodes
	#	nodes =  @monopartite_network.get_nodes_from_layer(@monopartite_layers[0].first).length	
	#	edges = @monopartite_network.get_edge_number
	#	
	#	@monopartite_network.randomize_monopartite_net_by_nodes([@monopartite_layers[0].first])
	#	
	#	random_nodes =  @monopartite_network.get_nodes_from_layer(@monopartite_layers[0].first).length	
	#	random_edges = @monopartite_network.get_edge_number
	#	assert_equal([nodes, edges], [random_nodes, random_edges])
#
	#end
	
#	def test_randomize_bipartite_net_by_nodes
#		layerA_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[0].first).length
#		layerB_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[1].first).length
#		edges = @network_obj.get_edge_number
#		
#		@network_obj.randomize_bipartite_net_by_nodes(@bipartite_layers.map(&:first))
#	
#		random_layerA_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[0].first).length
#		random_layerB_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[1].first).length
#		random_edges = @network_obj.get_edge_number
#		assert_equal([layerA_nodes, layerB_nodes, edges], [random_layerA_nodes, random_layerB_nodes, random_edges])
#	end
	
	#def test_randomize_monopartite_net_by_links
	#	previous_degree = @monopartite_network.get_degree(zscore = false)
	#
	#	@monopartite_network.randomize_monopartite_net_by_links([@monopartite_layers[0].first])
	#
	#	random_degree = @monopartite_network.get_degree(zscore = false)
	#	assert_equal(previous_degree, random_degree)
	#end
	
	def test_randomize_bipartite_net_by_links
		previous_degree = @network_obj.get_degree(zscore = false)
	
		@network_obj.randomize_bipartite_net_by_links(@bipartite_layers.map(&:first))
	
		random_degree = @network_obj.get_degree(zscore = false)
		assert_equal(previous_degree, random_degree)
	end

end
