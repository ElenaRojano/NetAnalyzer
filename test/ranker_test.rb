ROOT_PATH = File.dirname(__FILE__)
FILE_PATH = File.join(ROOT_PATH,'ranker_validation')
require File.join(ROOT_PATH, 'test_helper.rb')

class RankerTest < Minitest::Test

	def setup
		@matrix = Npy.load(File.join(FILE_PATH, 'kernel_for_validating'))
		@kernel_nodes = lst2arr(File.join(FILE_PATH, 'kernel_for_validating.lst'))
		@genes_seed = load_genes_by_group(File.join(FILE_PATH, 'seed_genes_for_validating'), ",")
		@seed_indexes = get_seed_indexes(@kernel_nodes, @genes_seed.values.flatten)
	end

	def ranked_genes2array(ranked_genes)
		test_ranked_genes = []
  		ranked_genes.each do |seed_name, rank_list|
  			rank_list.each do |ranked_gene| 
            	test_ranked_genes << ranked_gene.map{|el| el.to_s} + [seed_name]
          	end
  		end
  		return test_ranked_genes
	end

	def load_results(file_name)
		validate_ranked_genes = []
  		File.open(File.join(FILE_PATH, file_name)).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			validate_ranked_genes << fields
		end
		return validate_ranked_genes
	end

	def test_lst2arr
		test_nodes_loaded = @kernel_nodes
		validate_nodes_loaded = ["A","B","C","D","E"]
		assert_equal(validate_nodes_loaded, test_nodes_loaded)	
	end

	def test_load_genes_by_group
		test_seed_genes_loaded = @genes_seed
		validate_seed_genes_loaded = {'toy_group1'=> ["A","B"],'toy_group2'=> ["C","D","E"],'toy_group3'=> ["A","D"]}
		assert_equal(validate_seed_genes_loaded, test_seed_genes_loaded)	
	end

	def test_get_seed_indexes
		test_seed_indexes = @seed_indexes
		validate_seed_indexes={"A"=>0, "B"=>1, "C"=>2, "D"=>3, "E"=>4}
		assert_equal(validate_seed_indexes, test_seed_indexes)
	end

	def test_leave_one_out
		ranked_genes = {}
		@genes_seed.each do |seed_name, seed|  
        	ranked_genes[seed_name] = leave_one_out_validation(@matrix, @kernel_nodes, seed)
        end
  		test_ranked_genes = ranked_genes2array(ranked_genes)
  		validate_ranked_genes = load_results('leave_one_out_by_seedgene_results')
  		assert_equal(validate_ranked_genes, test_ranked_genes)
	end

	def test_rank_by_seed_gene
		ranked_genes = {}
		@genes_seed.each do |seed_name, seed|
			ranked_genes[seed_name] = rank_by_seedgen(@matrix, @seed_indexes, seed, @kernel_nodes)
  		end
  		test_ranked_genes = ranked_genes2array(ranked_genes)
  		validate_ranked_genes = []
  		validate_ranked_genes = load_results('rank_by_seedgene_results')
  		assert_equal(validate_ranked_genes, test_ranked_genes)
	end

	def test_get_filtered
		genes_to_keep = load_genes_by_group(File.join(FILE_PATH, 'genes2filter_for_validating'),",")
    	ranked_genes = {}
    	File.open(File.join(FILE_PATH, 'ranked_genes')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			seed_name = fields.shift
			values = fields[0].split(";").map{|row| row.split(",")}
			ranked_genes[seed_name] = values
		end

		filtered_ranked_genes = get_filtered(genes_to_keep, ranked_genes)
		test_filtered_ranked_genes = ranked_genes2array(filtered_ranked_genes)
    	validate_filtered_ranked_genes = load_results('filter_results')
  		assert_equal(validate_filtered_ranked_genes, test_filtered_ranked_genes)
  	end

	def test_get_top
		ranked_genes = {}
    	File.open(File.join(FILE_PATH, 'ranked_genes')).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			seed_name = fields.shift
			values = fields[0].split(";").map{|row| row.split(",")}
			ranked_genes[seed_name] = values
		end

		top_ranked_genes = get_top(2, ranked_genes)
		test_top_ranked_genes = ranked_genes2array(top_ranked_genes)
		validate_top_ranked_genes = load_results('top_results')
		assert_equal(validate_top_ranked_genes,test_top_ranked_genes)
	end

end