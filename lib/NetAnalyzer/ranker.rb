#!/usr/bin/env ruby

require 'numo/narray'
require "numo/linalg"

class Ranker
  attr_accessor :matrix, :nodes, :seeds, :ranking

  def initialize()
    @matrix = nil
    @nodes = [] # kernel_nodes
    @seeds = {} # genes_seed
    @reference_nodes = {}
    @ranking = {} # ranked_genes
  end

  def load_seeds(node_groups, sep: ',')
    @seeds = load_nodes_by_group(node_groups, sep: sep)
  end

  def load_references(node_groups, sep: ',')
    @reference_nodes = load_nodes_by_group(node_groups, sep: sep)
  end

  def load_nodes_by_group(node_groups, sep: ',')
    if File.exist?(node_groups) 
      group_nodes = load_node_groups_from_file(node_groups, sep: sep)
    else
      group_nodes = {"seed_genes" => node_groups.split(sep)}
    end
    return group_nodes
  end
          
  def load_node_groups_from_file(file, sep: ',')
    group_nodes = {}
    File.open(file).each do |line|
      set_name, nodes = line.chomp.split("\t")
      group_nodes[set_name] = nodes.split(sep)
    end
    return group_nodes
  end

  def load_nodes_from_file(file)
    File.open(file).each do |line|
      @nodes << line.chomp!
    end
  end

  def do_ranking(leave_one_out: false, threads: 0)
    if leave_one_out
      @seeds.each do |seed_name, seeds|  
        @ranking[seed_name] = leave_one_out_validation(seeds) # Benchmarking mode
      end
    else
      seed_indexes = get_seed_indexes
      seed_groups = @seeds.to_a # Array conversion needed for parallelization
      ranked_lists = Parallel.map(seed_groups, in_processes: threads) do |seed_name, seed|
        # The code in this block CANNOT modify nothing outside
        rank_list = rank_by_seed(seed_indexes, seed) # Production mode
        [seed_name, rank_list]
      end
      ranked_lists.each do |seed_name, rank_list| # Transfer resuls to hash
        @ranking[seed_name] = rank_list
      end
    end
  end

  def rank_by_seed(seed_indexes, seeds)
    ordered_gene_score = []
    genes_pos = seeds.map{|s| seed_indexes[s]}.compact
    number_of_seed_genes = genes_pos.length
    number_of_all_nodes = @nodes.length
    
    if number_of_seed_genes > 0
      subsets_gen_values = @matrix[genes_pos,true]
      integrated_gen_values = subsets_gen_values.sum(0)
      gen_list = 1.fdiv(number_of_seed_genes) * integrated_gen_values.inplace

      ordered_indexes = gen_list.sort_index # from smallest to largest
      
      last_val = nil
      n_elements = ordered_indexes.shape.first 
      n_elements.times do |pos|
        order_index = ordered_indexes[pos]    
        val = gen_list[order_index]
        node_name = @nodes[order_index]

        rank = get_position_for_items_with_same_score(pos, val, last_val, gen_list, n_elements, ordered_gene_score) # number of items behind
        rank = n_elements - rank # number of nodes below or equal
        rank_percentage = rank.fdiv(number_of_all_nodes)

        ordered_gene_score << [node_name, val, rank_percentage, rank]
        last_val = val
      end

      ordered_gene_score = ordered_gene_score.reverse # from largest to smallest
      ordered_gene_score = add_absolute_rank_column(ordered_gene_score)
    end
    
    return ordered_gene_score
  end

  def get_position_for_items_with_same_score(pos, val, prev_val, gen_list, n_elements, ordered_gene_score)
      members_behind = 0
      if !prev_val.nil?
        if prev_val < val
          members_behind = pos
        else 
          members_behind = n_elements - ordered_gene_score.last[3]
        end
      end
      return members_behind
  end

  def add_absolute_rank_column(ranking)
    ranking_with_new_column = []
    absolute_rank = 1
    n_rows = ranking.length
    n_rows.times do |row_pos|
      if row_pos == 0
        new_row = ranking[row_pos] << absolute_rank
        ranking_with_new_column << new_row
      else
        prev_val = ranking[row_pos-1][2]
        val = ranking[row_pos][2]
        if val > prev_val
          absolute_rank +=1
        end
        new_row = ranking[row_pos] << absolute_rank
        ranking_with_new_column << new_row
      end
    end
    return ranking_with_new_column
  end

  def leave_one_out_validation(seed_genes)
    group_number = seed_genes.length - 1
    one_out_seeds = seed_genes.combination(group_number).to_a
    
    out_genes_score = []
    one_out_seeds.each do |one_out_seed|
      gene_to_predict = seed_genes - one_out_seed
      gene_to_predict = gene_to_predict[0]
      ranked_one_out = get_individual_rank(one_out_seed, gene_to_predict)
      out_genes_score << ranked_one_out if !ranked_one_out.nil?
    end
    return out_genes_score
  end

  def get_individual_rank(seed_genes, node_of_interest)
    genes_pos = get_nodes_indexes(seed_genes)
    node_of_interest_pos = @nodes.find_index(node_of_interest)
    ordered_gene_score = nil

    if !genes_pos.empty? && !node_of_interest_pos.nil?

      subsets_gen_values = @matrix[genes_pos,true]
      integrated_gen_values = subsets_gen_values.sum(0)
      integrated_gen_values = 1.fdiv(genes_pos.length) * integrated_gen_values.inplace

      ref_value = integrated_gen_values[node_of_interest_pos]

      members_below_test = 0
      integrated_gen_values.each do |gen_value|
        members_below_test += 1 if gen_value >= ref_value
      end

      rank_percentage = members_below_test.fdiv(@nodes.length)
      rank = members_below_test
      rank_absolute = get_individual_absolute_rank(integrated_gen_values.to_a,ref_value)

      ordered_gene_score = [node_of_interest, ref_value, rank_percentage, rank, rank_absolute]
    end

    return ordered_gene_score

  end

  def get_individual_absolute_rank(values_list,ref_value)
    ref_pos = nil
    values_list = values_list.sort.reverse.uniq
    values_list.each_with_index do |value,pos|
      if value == ref_value
        ref_pos = pos+1
        break
      end
    end
    return ref_pos
  end


  def get_reference_ranks
    filtered_ranked_genes = {}
    
  	@ranking.each do |seed_name, ranking|
  	  filtered_ranked_genes[seed_name] = []
      next if @reference_nodes[seed_name].nil?
  	  ranking.each do |rank|
  	    next if !@reference_nodes[seed_name].include?(rank[0])
  	    filtered_ranked_genes[seed_name] << rank
  	  end
  	end
      
    return filtered_ranked_genes
  end

  def get_top(top_n)
  	top_ranked_genes = {}
  	@ranking.each do |seed_name, ranking|
        top_ranked_genes[seed_name] = ranking[0..top_n-1] if !ranking.nil?
  	end
  	return top_ranked_genes
  end

  def get_nodes_indexes(nodes)
    node_indxs = []
    nodes.each do |node|
      index_node = @nodes.find_index(node)
      node_indxs << index_node if !index_node.nil?
    end
    return node_indxs
  end

  def get_seed_indexes
    indexes = {}
    @seeds.values.flatten.each do |node|
      if !indexes.include?(node)
        indx = @nodes.index(node)
        indexes[node] = indx if !indx.nil?
      end
    end
    return indexes
  end
end