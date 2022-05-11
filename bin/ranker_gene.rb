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

def rank_by_seedgen(kernel_matrix, kernels_nodes, seed_genes)
  #gen_pos = seed_genes.map{|gen| kernels_nodes.find_index(gen)}
  genes_pos = []
  seed_genes.each do |gene|
    index_gene = kernels_nodes.find_index(gene)
    genes_pos.append(index_gene) if !index_gene.nil?
  end

  number_of_seed_genes = genes_pos.length
  
  return nil if number_of_seed_genes == 0
  
  subsets_gen_values = []
  genes_pos.each do |gen_pos|
    subsets_gen_values.append(kernel_matrix[gen_pos,true])
  end

  gen_list = 1.fdiv(number_of_seed_genes) * subsets_gen_values.sum
  percentiles = (1..gen_list.length).to_a
  percentiles.map!{|percentile| percentile/percentiles.length.to_f}

  ordered_indexes = gen_list.sort_index.to_a.reverse
  ordered_gene_score = []
  ordered_indexes.each.with_index{|order_index, pos| ordered_gene_score.append([kernels_nodes[order_index], gen_list[order_index], percentiles[pos], pos])}

  return ordered_gene_score
end

def lst2arr(lst_file)
	nodes = []

	File.open(lst_file,"r").each do |line|
	line.chomp!
	nodes.append(line)
	end
  
	return nodes
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

genes_seed.each do |seed_name, seed|  
    ranked_genes[seed_name] = rank_by_seedgen(matrix, kernel_nodes, seed)
end

genes_to_keep = []
if !options[:filter].nil?
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


