ROOT_PATH = File.dirname(__FILE__)
require File.join(ROOT_PATH, 'test_helper.rb')

class Net_parserTest < Minitest::Test

	def setup
		@bipartite_layers = [[:main, /M[0-9]+/], [:projection, /P[0-9]+/]]

		@monopartite_layers = [[:main, /M[0-9]+/], [:main, /M[0-9]+/]]
		@monopartite_network_node_names = File.join(DATA_TEST_PATH, 'monopartite_network_node_names.txt')


	end

	def test_load_pairs
		options = {:input_format => 'pair', :input_file => File.join(DATA_TEST_PATH, 'bipartite_network_for_validating.txt'), :layers => @bipartite_layers, :split_char => "\t"}
		network_obj = Net_parser.load(options)
		test_main_layer = network_obj.get_nodes_layer([:main]).length
		assert_equal(6, test_main_layer)
		test_projection_layer = network_obj.get_nodes_layer([:projection]).length
		assert_equal(10, test_projection_layer)
		test_connections = network_obj.get_edge_number
		assert_equal(40, test_connections)	
	end

	def test_load_bin_matrix
		options = {:input_format => 'bin', :input_file => File.join(DATA_TEST_PATH, 'monopartite_network_bin_matrix.npy'), :layers => @monopartite_layers, :node_file => @monopartite_network_node_names}
		monopartite_network_by_bin_matrix = Net_parser.load(options)
		test_adjacency_matrices = monopartite_network_by_bin_matrix.adjacency_matrices
		adjacency_matrices_values = Numo::DFloat[[0, 1, 1, 0, 0],[1, 0, 0, 0, 0], [1, 0, 0, 0, 1], [0, 0, 0, 0, 1], [0, 0, 1, 1, 0]]
		expected_adjacency_matrices = {[:main, :main] => [adjacency_matrices_values, ['A', 'B', 'C', 'D', 'E'], ['A', 'B', 'C', 'D', 'E']]} 
		assert_equal expected_adjacency_matrices, test_adjacency_matrices
	end

	def test_load_plain_matrix
		options = {:input_format => 'matrix', :input_file => File.join(DATA_TEST_PATH, 'monopartite_network_matrix'), :layers => @monopartite_layers, :node_file => @monopartite_network_node_names}
		monopartite_network_by_plain_matrix = Net_parser.load(options)
		test_adjacency_matrices = monopartite_network_by_plain_matrix.adjacency_matrices
		adjacency_matrices_values = Numo::DFloat[[0, 0, 1, 0, 1],[0, 0, 0, 1, 0], [1, 0, 0, 0, 0], [0, 1, 0, 0, 1], [1, 0, 0, 1, 0]]
		expected_adjacency_matrices = {[:main, :main] => [adjacency_matrices_values, ['A', 'B', 'C', 'D', 'E'], ['A', 'B', 'C', 'D', 'E']]}
		assert_equal expected_adjacency_matrices, test_adjacency_matrices
	end

	def test_load_network_by_pairs
		network_obj = Net_parser.load_network_by_pairs(File.join(DATA_TEST_PATH, 'bipartite_network_for_validating.txt'), @bipartite_layers)
		test_main_layer = network_obj.get_nodes_layer([:main]).length
		assert_equal(6, test_main_layer)
		test_projection_layer = network_obj.get_nodes_layer([:projection]).length
		assert_equal(10, test_projection_layer)
		test_connections = network_obj.get_edge_number
		assert_equal(40, test_connections)	
	end

	def test_load_network_by_bin_matrix
		monopartite_network_by_bin_matrix = Net_parser.load_network_by_bin_matrix(File.join(DATA_TEST_PATH, 'monopartite_network_bin_matrix.npy'), @monopartite_network_node_names, @monopartite_layers)
		test_adjacency_matrices = monopartite_network_by_bin_matrix.adjacency_matrices
		adjacency_matrices_values = Numo::DFloat[[0, 1, 1, 0, 0],[1, 0, 0, 0, 0], [1, 0, 0, 0, 1], [0, 0, 0, 0, 1], [0, 0, 1, 1, 0]]
		expected_adjacency_matrices = {[:main, :main] => [adjacency_matrices_values, ['A', 'B', 'C', 'D', 'E'], ['A', 'B', 'C', 'D', 'E']]} 
		assert_equal expected_adjacency_matrices, test_adjacency_matrices
	end

	def test_load_network_by_plain_matrix
		monopartite_network_by_plain_matrix = Net_parser.load_network_by_plain_matrix(File.join(DATA_TEST_PATH, 'monopartite_network_matrix'), @monopartite_network_node_names, @monopartite_layers)
		test_adjacency_matrices = monopartite_network_by_plain_matrix.adjacency_matrices
		adjacency_matrices_values = Numo::DFloat[[0, 0, 1, 0, 1],[0, 0, 0, 1, 0], [1, 0, 0, 0, 0], [0, 1, 0, 0, 1], [1, 0, 0, 1, 0]]
		expected_adjacency_matrices = {[:main, :main] => [adjacency_matrices_values, ['A', 'B', 'C', 'D', 'E'], ['A', 'B', 'C', 'D', 'E']]}
		assert_equal expected_adjacency_matrices, test_adjacency_matrices
	end

end