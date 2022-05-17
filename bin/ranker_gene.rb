#! /usr/bin/env ruby
require 'optparse'
require 'npy'
require 'numo/narray'

########################### FUNCTIONS #######################
#############################################################


def load_seed_genes(seed_genes, seed_genes_sep)
  if File.exist?(seed_genes) 
    seed_genes = load_genes_from_file(seed_genes, seed_genes_sep)
  else
    seed_genes = {"seed_genes" => seed_genes.split(",")}
  end
  return seed_genes
end
        

def load_genes_from_file(file, seed_genes_sep)
  seed_genes = {}
  File.open(file).each do |line|
    line = line.chomp.split("\t")
    set_name = line[0]
    genes = line[1].split(seed_genes_sep) 
    seed_genes[set_name] = genes
  end
  return seed_genes
end


def load_filter(file)
  filter = {}
  if File.exist?(file)
    File.open(file).each do |line|
      line = line.chomp.split("\t")
      filter[line[0]] = line[1].split(",") 
    end
  end
  return filter
end


def lst2arr(lst_file)
  nodes = []

  File.open(lst_file,"r").each do |line|
  line.chomp!
  nodes.append(line)
  end
  
  return nodes
end


def rank_by_seedgen(kernel_matrix, kernels_nodes, seed_genes)
  genes_pos = get_nodes_indexes(kernels_nodes, seed_genes)
  number_of_seed_genes = genes_pos.length
  number_of_all_nodes = kernels_nodes.length
  
  return nil if number_of_seed_genes == 0
  
  subsets_gen_values = kernel_matrix[genes_pos,true]
  integrated_gen_values = subsets_gen_values.sum(0)
  gen_list = 1.fdiv(number_of_seed_genes) * integrated_gen_values

  ordered_indexes = gen_list.sort_index.to_a.reverse
  
  ordered_gene_score = []
  ordered_indexes.each_with_index do |order_index, pos| 
    val = gen_list[order_index]
    node_name = kernels_nodes[order_index]
    if pos == 0
      members_below = 0
    elsif pos > 0
      prev_val = gen_list[ordered_indexes[pos-1]]
      if prev_val > val
        members_below = pos
      else 
        members_below = ordered_gene_score.last[3]
      end
    end
    rank = members_below
    rank_percentage = members_below.fdiv(number_of_all_nodes)
    ordered_gene_score.append([kernels_nodes[order_index], val, rank_percentage, rank])
  end

  return ordered_gene_score
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


  ref_value = integrated_gen_values[node_of_interest_pos]

  members_below_test = 0
  integrated_gen_values.each do |gen_value|
    members_below_test += 1 if gen_value > ref_value
  end

  rank_percentage = members_below_test.fdiv(kernels_nodes.length)
  rank = members_below_test

  return [node_of_interest, ref_value, rank_percentage, rank]
end


def get_nodes_indexes(all_nodes, nodes_of_interest)
  nodes_of_interest_pos = []
  nodes_of_interest.each do |node_of_interest|
    index_node = all_nodes.find_index(node_of_interest)
    nodes_of_interest_pos.append(index_node) if !index_node.nil?
  end
  return nodes_of_interest_pos
end

########################### OPTPARSE ########################
#############################################################

options = {}
OptionParser.new do  |opts|

  options[:kernel_file] = nil
  opts.on("-k","--input_kernels KER", "The roots from each kernel to integrate") do |ker|
    options[:kernel_file] = ker
  end

  options[:node_file] = nil
  opts.on("-n","--input_nodes NODE", "The list of node for each kernel in lst format") do |node_file|
    options[:node_file] = node_file
  end

  options[:genes_seed] = nil
  opts.on("-s","--genes_seed SEED", "The name of the gene to look for backups") do |genes_seed|
    options[:genes_seed] = genes_seed
  end

  options[:seed_genes_sep] = ","
  opts.on("-S","--genes_seed_sep SEP", "Separator of seed genes. Only use when -s point to a file") do |genes_seed|
    options[:genes_seed] = genes_seed
  end

  options[:filter] = nil   
  opts.on("-f","--filter NAME", "PATH to file with seed_name and genes to keep in output") do |file|
    options[:filter] = file
  end

  options[:leave_one_out] = false 
  opts.on("-l","--leave_one_out", "Perform leave one out from a seed genes group") do 
    options[:leave_one_out] = true
  end
 
  options[:top_n] = nil   
  opts.on("-t","--top_n INT", "Top N genes to print in output") do |str|
    options[:top_n] = str.to_i
  end

  options[:output_name] = "ranked_genes"
  opts.on("-o","--output_name NAME", "The name of the ranked file") do |output_name|
    options[:output_name] = output_name
  end

end.parse!

########################### MAIN ############################
#############################################################

matrix = Npy.load(options[:kernel_file])
kernel_nodes = lst2arr(options[:node_file])
genes_seed = load_seed_genes(options[:genes_seed], options[:seed_genes_sep])
output_name = options[:output_name]

ranked_genes = {}
if options[:leave_one_out]
  genes_seed.each do |seed_name, seed|  
    ranked_genes[seed_name] = leave_one_out_validation(matrix, kernel_nodes, seed) # Benchmarking mode
  end
else
  genes_seed.each do |seed_name, seed|  
    ranked_genes[seed_name] = rank_by_seedgen(matrix, kernel_nodes, seed) # Production mode
  end
end

genes_to_keep = []
if !options[:filter].nil? && !options[:leave_one_out]
  genes_to_keep = load_filter(options[:filter])
  ranked_genes.reject!{|seed_name, seed| genes_to_keep[seed_name].nil?}
end

if genes_to_keep.empty?
  filtered_ranked_genes = ranked_genes
else
  filtered_ranked_genes = {}
  ranked_genes.each do |seed_name, ranking|
    filtered_ranked_genes[seed_name] = []
    next if ranking.nil?
    ranking.each do |line|
      next if !genes_to_keep[seed_name].include?(line[0])
      filtered_ranked_genes[seed_name] << line
    end
  end
end

top_n = options[:top_n]
top_n = matrix.length if top_n.nil?

discarded_seeds = []
File.open("#{output_name}_all_candidates" ,'w') do |f|
  filtered_ranked_genes.each do |seed_name, ranking|
    top_counter = 0
    if !ranking.nil?
      ranking.each do |ranked_gene| 
        break if top_counter > top_n
        f.puts "#{ranked_gene.join("\t")}\t#{seed_name}"
        top_counter += 1
      end
      
    else 
      discarded_seeds << seed_name
    end
  end
end

if !discarded_seeds.empty?
  File.open("#{output_name}_discarded",'w') do |f|
      discarded_seeds.each do |seed_name|
        f.puts "#{seed_name}\t#{genes_seed[seed_name].join(options[:seed_genes_sep])}"
      end
  end
end


