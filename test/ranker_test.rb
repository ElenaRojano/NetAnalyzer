ROOT_PATH = File.dirname(__FILE__)
FILE_PATH = File.join(ROOT_PATH,'ranker_validation')
require File.join(ROOT_PATH, 'test_helper.rb')

class RankerTest < Minitest::Test

	def load_results(file_name)
		validate_ranked_genes = []
  		File.open(File.join(FILE_PATH, file_name)).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			validate_ranked_genes << fields
		end
		return validate_ranked_genes
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

	def load_ranking(file)
		ranked_genes = {}
    	File.open(file).each("\n") do |line|
			line.chomp!
			fields = line.split("\t")
			seed_name = fields.shift
			values = fields[0].split(";").map{|row| row.split(",")}
			values.each do |row|
				row.map!.with_index do |val,idx| 
					if idx == 0
						val 
					elsif idx >= 3 
						val = val.to_i
					else
						val = val.to_f
					end
				end
			end
			ranked_genes[seed_name] = values
		end
		return ranked_genes
	end

	def setup
		@ranker = Ranker.new()
		@ranker.matrix = Npy.load(File.join(FILE_PATH, 'kernel_for_validating'))
		@ranker.load_seeds(File.join(FILE_PATH, 'seed_genes_for_validating'), sep: ",") # Should be in its test method but modifications are deleted from one test to another
		@ranker.load_nodes_from_file(File.join(FILE_PATH, 'kernel_for_validating.lst')) # Should be in its test method but modifications are deleted from one test to another

		@ranker_with_ranking = Ranker.new()
		@ranker_with_ranking.load_references(File.join(FILE_PATH, 'genes2filter_for_validating'))
		@ranker_with_ranking.ranking = load_ranking(File.join(FILE_PATH, 'ranked_genes'))
	end

	def test_load_nodes_from_file
		assert_equal(["A","B","C","D","E"], @ranker.nodes)	
	end

	def test_load_seeds
		validate_seed_genes_loaded = {'toy_group1'=> ["A","B"],'toy_group2'=> ["C","D","E"],'toy_group3'=> ["A","D"]}
		assert_equal(validate_seed_genes_loaded, @ranker.seeds)	
	end

	def test_get_seed_indexes
		validate_seed_indexes={"A"=>0, "B"=>1, "C"=>2, "D"=>3, "E"=>4}
		assert_equal(validate_seed_indexes, @ranker.get_seed_indexes)
	end

	def test_leave_one_out
		@ranker.do_ranking(leave_one_out: true)
  		test_ranked_genes = ranked_genes2array(@ranker.ranking)
  		assert_equal(load_results('leave_one_out_by_seedgene_results'), test_ranked_genes)
	end

	def test_rank_by_seed
		@ranker.do_ranking()
  		test_ranked_genes = ranked_genes2array(@ranker.ranking)
  		assert_equal(load_results('rank_by_seedgene_results'), test_ranked_genes)
	end

	def test_get_filtered
		test_filtered_ranked_genes = ranked_genes2array(@ranker_with_ranking.get_reference_ranks)
  		assert_equal(load_results('filter_results'), test_filtered_ranked_genes)
  	end

	def test_get_top
		test_top_ranked_genes = ranked_genes2array(@ranker_with_ranking.get_top(2))
		assert_equal(load_results('top_results'),test_top_ranked_genes)
	end

	def test_filtered_top_compatibility
		filtered_ranked_genes = ranked_genes2array(@ranker_with_ranking.get_reference_ranks)
		test_top_ranked_genes = ranked_genes2array(@ranker_with_ranking.get_top(2))
		assert_equal(load_results('top_results'),test_top_ranked_genes)
	end

end