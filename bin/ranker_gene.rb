#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$LOAD_PATH.unshift(File.expand_path(File.join(ROOT_PATH, '..', 'lib')))

require 'optparse'
require 'npy'
require 'parallel'
require 'NetAnalyzer'

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

  options[:threads] = 0
  opts.on( '-T', '--threads INTEGER', 'Number of threads to use in computation, one thread will be reserved as manager.' ) do |opt|
      options[:threads] = opt.to_i - 1
  end
end.parse!

########################### MAIN ############################
#############################################################

matrix = Npy.load(options[:kernel_file])
kernel_nodes = lst2arr(options[:node_file])
genes_seed = load_genes_by_group(options[:genes_seed], options[:seed_genes_sep])
output_name = options[:output_name]
ranked_genes = {}

if options[:leave_one_out]
  genes_seed.each do |seed_name, seed|  
    ranked_genes[seed_name] = leave_one_out_validation(matrix, kernel_nodes, seed) # Benchmarking mode
  end
else
  seed_indexes = get_seed_indexes(kernel_nodes, genes_seed.values.flatten)
  seed_groups = genes_seed.to_a # Array conversion needed for parallelization
  ranked_lists = Parallel.map(seed_groups, in_processes: options[:threads]) do |seed_name, seed|
    # The code in this block CANNOT modify nothing outside
    rank_list = rank_by_seedgen(matrix, seed_indexes, seed, kernel_nodes) # Production mode
    [seed_name, rank_list]
  end
  ranked_lists.each do |seed_name, rank_list| # Transfer resuls to hash
    ranked_genes[seed_name] = rank_list
  end
end

discarded_seeds = ranked_genes.select{|seed_name, ranks| ranks.empty?}.keys
if !discarded_seeds.empty?
  File.open("#{output_name}_discarded",'w') do |f|
    discarded_seeds.each do |seed_name|
      f.puts "#{seed_name}\t#{genes_seed[seed_name].join(options[:seed_genes_sep])}"
    end
  end
end

if !options[:filter].nil? && !options[:leave_one_out]
  genes_to_keep = load_genes_by_group(options[:filter],",")
  filtered_ranked_genes = get_filtered(genes_to_keep, ranked_genes)
else
  filtered_ranked_genes = ranked_genes
end

if !options[:top_n].nil?
  filtered_ranked_genes = get_top(options[:top_n], filtered_ranked_genes)
end

if !filtered_ranked_genes.empty?
  File.open("#{output_name}_all_candidates" ,'w') do |f|
    filtered_ranked_genes.each do |seed_name, ranking|
      ranking.each do |ranked_gene| 
        f.puts "#{ranked_gene.join("\t")}\t#{seed_name}"
      end
    end
  end
end



