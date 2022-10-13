ROOT_PATH = File.dirname(__FILE__)
require File.join(ROOT_PATH, 'test_helper.rb')

class NetworkTest < Minitest::Test

	def setup
		@bipartite_layers = [[:main, /M[0-9]+/], [:projection, /P[0-9]+/]]
		@network_obj = Net_parser.load_network_by_pairs(File.join(ROOT_PATH, 'bipartite_network_for_validating.txt'), @bipartite_layers)
		@network_obj.generate_adjacency_matrix(@bipartite_layers[0].first, @bipartite_layers[1].first)

		@monopartite_layers = [[:main, /\w/], [:main, /\w/]]
		@monopartite_network = Net_parser.load_network_by_pairs(File.join(ROOT_PATH, 'monopartite_network_for_validating.txt'), @monopartite_layers)
		@monopartite_network.generate_adjacency_matrix(@monopartite_layers[0].first, @monopartite_layers[0].first)
	end

	def test_load_network_by_pairs
		test_main_layer = @network_obj.get_nodes_layer([:main]).length
		assert_equal(6, test_main_layer)
		test_projection_layer = @network_obj.get_nodes_layer([:projection]).length
		assert_equal(10, test_projection_layer)
		test_connections = @network_obj.get_edge_number
		assert_equal(40, test_connections)	
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

	def test_randomize_monopartite_net_by_nodes
		nodes =  @monopartite_network.get_nodes_from_layer(@monopartite_layers[0].first).length	
		edges = @monopartite_network.get_edge_number
		
		@monopartite_network.randomize_monopartite_net_by_nodes([@monopartite_layers[0].first])
		
		random_nodes =  @monopartite_network.get_nodes_from_layer(@monopartite_layers[0].first).length	
		random_edges = @monopartite_network.get_edge_number
		assert_equal([nodes, edges], [random_nodes, random_edges])

	end
	
	def test_randomize_bipartite_net_by_nodes
		layerA_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[0].first).length
		layerB_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[1].first).length
		edges = @network_obj.get_edge_number
		
		@network_obj.randomize_bipartite_net_by_nodes(@bipartite_layers.map(&:first))
	
		random_layerA_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[0].first).length
		random_layerB_nodes = @network_obj.get_nodes_from_layer(@bipartite_layers[1].first).length
		random_edges = @network_obj.get_edge_number
		assert_equal([layerA_nodes, layerB_nodes, edges], [random_layerA_nodes, random_layerB_nodes, random_edges])

	end
	
	def test_randomize_monopartite_net_by_links
		previous_degree = @monopartite_network.get_degree(zscore = false)
	
		@monopartite_network.randomize_monopartite_net_by_links([@monopartite_layers[0].first])
	
		random_degree = @monopartite_network.get_degree(zscore = false)
		assert_equal(previous_degree, random_degree)
	end
	
	def test_randomize_bipartite_net_by_links
		previous_degree = @network_obj.get_degree(zscore = false)
	
		@network_obj.randomize_bipartite_net_by_links(@bipartite_layers.map(&:first))
	
		random_degree = @network_obj.get_degree(zscore = false)
		assert_equal(previous_degree, random_degree)
	end
end