ROOT_PATH = File.dirname(__FILE__)
require File.join(ROOT_PATH, 'test_helper.rb')

class NetworkTest < Test::Unit::TestCase

	def setup
		layers = [[:main, /M[0-9]+/], [:projection, /P[0-9]+/]]
		@network_obj = Network.new(layers.map{|layer| layer.first})
		@network_obj.load_network_by_pairs(File.join(ROOT_PATH, 'bipartite_network_for_validating.txt'), layers)
	end

	def test_load_network_by_pairs
		test_main_layer = @network_obj.get_nodes_layer([:main]).length
		assert_equal(20, test_main_layer)
		test_projection_layer = @network_obj.get_nodes_layer([:projection]).length
		assert_equal(10, test_projection_layer)
		test_connections = @network_obj.get_edge_number
		assert_equal(110, test_connections)	
	end
end