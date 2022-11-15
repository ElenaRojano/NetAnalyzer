ROOT_PATH = File.dirname(__FILE__)
require File.join(ROOT_PATH, 'test_helper.rb')

class Net_parserTest < Minitest::Test

	def setup
		@bipartite_layers = [[:main, /M[0-9]+/], [:projection, /P[0-9]+/]]
		@network_obj = Net_parser.load_network_by_pairs(File.join(DATA_TEST_PATH, 'bipartite_network_for_validating.txt'), @bipartite_layers)
		@network_obj.generate_adjacency_matrix(@bipartite_layers[0].first, @bipartite_layers[1].first) 

		@monopartite_layers = [[:main, /M[0-9]+/], [:main, /M[0-9]+/]]
		@monopartite_network_node_names = File.join(DATA_TEST_PATH, 'monopartite_network_node_names.txt')
	end

	def test_load_network_by_pairs
		test_main_layer = @network_obj.get_nodes_layer([:main]).length
		assert_equal(6, test_main_layer)
		test_projection_layer = @network_obj.get_nodes_layer([:projection]).length
		assert_equal(10, test_projection_layer)
		test_connections = @network_obj.get_edge_number
		assert_equal(40, test_connections)	
	end

	def test_load_network_by_bin_matrix
		monopartite_network_by_bin_matrix = Net_parser.load_network_by_bin_matrix(File.join(DATA_TEST_PATH, 'monopartite_network_bin_matrix.npy'), @monopartite_network_node_names, @monopartite_layers)
		test_adjacency_matrices = monopartite_network_by_bin_matrix.adjacency_matrices
		adjacency_matrices_values = Numo::DFloat[[0, 1, 1, 0, 0],[1, 0, 0, 0, 0], [1, 0, 0, 0, 1], [0, 0, 0, 0, 1], [0, 0, 1, 1, 0]]
		adjacency_matrices = {[:main, :main] => [adjacency_matrices_values, ['A', 'B', 'C', 'D', 'E'], ['A', 'B', 'C', 'D', 'E']]} # Numo::DFloat.new(5, 5).seq#.[[0, 0, 0, 0, 0][0, 0, 0, 0, 0][0, 0, 0, 0, 0][0, 0, 0, 0, 0][0, 0, 0, 0, 0]]
		assert_equal adjacency_matrices, test_adjacency_matrices
	end

	def test_load_network_by_plain_matrix
		monopartite_network_by_plain_matrix = Net_parser.load_network_by_plain_matrix(File.join(DATA_TEST_PATH, 'monopartite_network_matrix'), @monopartite_network_node_names, @monopartite_layers)
		test_adjacency_matrices = monopartite_network_by_plain_matrix.adjacency_matrices
		adjacency_matrices_values = Numo::DFloat[[0, 0, 1, 0, 1],[0, 0, 0, 1, 0], [1, 0, 0, 0, 0], [0, 1, 0, 0, 1], [1, 0, 0, 1, 0]]
		adjacency_matrices = {[:main, :main] => [adjacency_matrices_values, ['A', 'B', 'C', 'D', 'E'], ['A', 'B', 'C', 'D', 'E']]} # Numo::DFloat.new(5, 5).seq#.[[0, 0, 0, 0, 0][0, 0, 0, 0, 0][0, 0, 0, 0, 0][0, 0, 0, 0, 0][0, 0, 0, 0, 0]]
		assert_equal adjacency_matrices, test_adjacency_matrices
	end

end