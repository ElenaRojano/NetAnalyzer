ROOT_PATH = File.dirname(__FILE__)
require File.join(ROOT_PATH, 'test_helper.rb')

class NetworkTest < Minitest::Test

	def setup
		layers = [[:main, /M[0-9]+/], [:projection, /P[0-9]+/]]
		@network_obj = Network.new(layers.map{|layer| layer.first})
		@network_obj.load_network_by_pairs(File.join(ROOT_PATH, 'bipartite_network_for_validating.txt'), layers)
	end

	def test_load_network_by_pairs
		test_main_layer = @network_obj.get_nodes_layer([:main]).length
		assert_equal(9, test_main_layer)
		test_projection_layer = @network_obj.get_nodes_layer([:projection]).length
		assert_equal(16, test_projection_layer)
		test_connections = @network_obj.get_edge_number
		assert_equal(64, test_connections)	
	end

	def test_get_counts_association
		test_association = @network_obj.get_counts_association([:main], :projection) 
		test_association.map!{|a| [a[0], a[1], a[2].round(6)]}
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
end