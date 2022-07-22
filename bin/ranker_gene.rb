#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$LOAD_PATH.unshift(File.expand_path(File.join(ROOT_PATH, '..', 'lib')))

require 'optparse'
require 'npy'
require 'parallel'
require 'NetAnalyzer'


########################### METHODS ########################
#############################################################

def write_ranking(file, ranking_list)
  File.open(file ,'w') do |f|
    ranking_list.each do |seed_name, ranking|
      ranking.each do |ranked_gene| 
        f.puts "#{ranked_gene.join("\t")}\t#{seed_name}"
      end
    end
  end
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

  options[:output_top] = nil   
  opts.on("--output_top PATH", "File to save Top N genes") do |path|
    options[:output_top] = path
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

ranker = Ranker.new()
ranker.matrix = Npy.load(options[:kernel_file])
ranker.load_nodes_from_file(options[:node_file])
ranker.load_seeds(options[:genes_seed], sep: options[:seed_genes_sep])
ranker.do_ranking(leave_one_out: options[:leave_one_out], threads: options[:threads])
rankings = ranker.ranking

discarded_seeds = rankings.select{|seed_name, ranks| ranks.empty?}.keys
if !discarded_seeds.empty?
  File.open("#{options[:output_name]}_discarded",'w') do |f|
    discarded_seeds.each do |seed_name|
      f.puts "#{seed_name}\t#{ranker.seeds[seed_name].join(options[:seed_genes_sep])}"
    end
  end
end

if !options[:top_n].nil?
  top_n = ranker.get_top(options[:top_n])
  if options[:output_top].nil?
    rankings = top_n
  else
    write_ranking(options[:output_top], top_n)
  end
end

if !options[:filter].nil? && !options[:leave_one_out]
  ranker.load_references(options[:filter], sep: ",")
  rankings = ranker.get_reference_ranks
end

if !rankings.empty?
  write_ranking("#{options[:output_name]}_all_candidates", rankings)
end



