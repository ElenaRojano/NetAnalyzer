#!/usr/bin/env ruby

require 'numo/narray'
require "numo/linalg"

def load_genes_by_group(group_genes, group_genes_sep)
  if File.exist?(group_genes) 
    group_genes = load_genes_from_file(group_genes, group_genes_sep)
  else
    group_genes = {"seed_genes" => group_genes.split(",")}
  end
  return group_genes
end
        
def load_genes_from_file(file, group_genes_sep)
  group_genes = {}
  File.open(file).each do |line|
    line = line.chomp.split("\t")
    set_name = line[0]
    genes = line[1].split(group_genes_sep) 
    group_genes[set_name] = genes
  end
  return group_genes
end

def lst2arr(lst_file)
  nodes = []

  File.open(lst_file,"r").each do |line|
  line.chomp!
  nodes.append(line)
  end
  
  return nodes
end

def rank_by_seedgen(kernel_matrix, seed_indexes, seed_genes, kernels_nodes)
  ordered_gene_score = []
  genes_pos = seed_genes.map{|s| seed_indexes[s]}.compact
  number_of_seed_genes = genes_pos.length
  number_of_all_nodes = kernels_nodes.length
  
  if number_of_seed_genes > 0
    subsets_gen_values = kernel_matrix[genes_pos,true]
    integrated_gen_values = subsets_gen_values.sum(0)
    gen_list = 1.fdiv(number_of_seed_genes) * integrated_gen_values.inplace

    ordered_indexes = gen_list.sort_index
    
    last_val = nil
    n_elements = ordered_indexes.shape.first
    n_elements.times do |pos|
      order_index = ordered_indexes[pos]    
      val = gen_list[order_index]
      node_name = kernels_nodes[order_index]
      rank = get_position_for_items_with_same_score(pos, val, last_val, gen_list, n_elements, ordered_gene_score)
      rank = n_elements - rank
      rank_percentage = rank.fdiv(number_of_all_nodes)
      ordered_gene_score << [node_name, val, rank_percentage, rank]
      last_val = val
    end
  end
  return ordered_gene_score.reverse
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

def leave_one_out_validation(kernel_matrix, kernels_nodes, seed_genes)
  group_number = seed_genes.length - 1
  one_out_seeds = seed_genes.combination(group_number).to_a
  
  out_genes_score = []
  one_out_seeds.each do |one_out_seed|
    gene_to_predict = seed_genes - one_out_seed
    gene_to_predict = gene_to_predict[0]
    ranked_one_out = get_individual_rank(kernel_matrix, kernels_nodes, one_out_seed, gene_to_predict)
    out_genes_score.append(ranked_one_out) if !ranked_one_out.nil?
  end

  return out_genes_score
end

def get_individual_rank(kernel_matrix, kernels_nodes, seed_genes, node_of_interest)
  genes_pos = get_nodes_indexes(kernels_nodes, seed_genes)
  node_of_interest_pos = kernels_nodes.find_index(node_of_interest)

  return nil if genes_pos.empty? || node_of_interest_pos.nil?

  subsets_gen_values = kernel_matrix[genes_pos,true]
  integrated_gen_values = subsets_gen_values.sum(0)
  integrated_gen_values = 1.fdiv(genes_pos.length) * integrated_gen_values.inplace


  ref_value = integrated_gen_values[node_of_interest_pos]

  members_below_test = 0
  integrated_gen_values.each do |gen_value|
    members_below_test += 1 if gen_value >= ref_value
  end

  rank_percentage = members_below_test.fdiv(kernels_nodes.length)
  rank = members_below_test

  return [node_of_interest, ref_value, rank_percentage, rank]
end

def get_filtered(genes_to_keep, ranked_genes)
  ranked_genes.reject!{|seed_name, ranking| genes_to_keep[seed_name].nil?}

  filtered_ranked_genes = {}

	ranked_genes.each do |seed_name, ranking|
	  filtered_ranked_genes[seed_name] = []
	  ranking.each do |rank|
	    next if !genes_to_keep[seed_name].include?(rank[0])
	    filtered_ranked_genes[seed_name] << rank
	  end
	end
    
    return filtered_ranked_genes
end

def get_top(top_n, ranked_genes)
	top_ranked_genes = {}
	ranked_genes.each do |seed_name, ranking|
      top_ranked_genes[seed_name] = ranking[0..top_n-1] if !ranking.nil?
	end
	return top_ranked_genes
end

def get_nodes_indexes(all_nodes, nodes_of_interest)
  nodes_of_interest_pos = []
  nodes_of_interest.each do |node_of_interest|
    index_node = all_nodes.find_index(node_of_interest)
    nodes_of_interest_pos.append(index_node) if !index_node.nil?
  end
  return nodes_of_interest_pos
end

def get_seed_indexes(all_nodes, all_nodes_of_interest)
  indexes = {}
  all_nodes_of_interest.each do |node|
    if !indexes.include?(node)
      indx = all_nodes.index(node)
      indexes[node] = indx if !indx.nil?
    end
  end
  return indexes
end
