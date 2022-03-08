#!/usr/bin/env ruby

require 'optparse'

##############################
#FUNCTIONS
##############################


def load_clusters(options)
	clusters = {}
	File.open(options[:input_file]).each do |line|
		line = line.chomp.split(options[:column_sep])
		cluster = line[options[:cluster_index]]
		clusters[cluster] = [] if clusters[cluster].nil?
		node = line[options[:node_index]]
		node = node.split(options[:node_sep]) if !options[:node_sep].nil?
		clusters[cluster] << node
		clusters[cluster].flatten!
	end
	return clusters
end


def random_sample(nodes, replacement, all_sizes, seed)
	random_clusters = {}
	nodes_list = nodes.dup
	all_sizes.each_with_index do |cluster_size, counter|
		abort("Not enough nodes to generate clusters. Please activate replacement or change random mode") if cluster_size > nodes_list.size
		random_nodes = nodes_list.uniq.sample(cluster_size, random: Random.new(seed))
		if !replacement 
			nodes_list = nodes_list - random_nodes
		end
		random_clusters["#{counter}_random"] = random_nodes
		seed += 1
	end
	return random_clusters
end

def write_clusters(clusters, output_file, sep)
	File.open(output_file, 'w') do |outfile|
		clusters.each do |cluster, nodes|
			nodes = [nodes.join(sep)] if !sep.nil?
			nodes.each do |node|
				outfile.puts [cluster, node].flatten.join("\t")
			end
		end		
	end	
end

##############################
#OPTPARSE
##############################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file") do |input_file|
    options[:input_file] = input_file
  end

  options[:node_index] = 1
  opts.on("-N", "--node_column INTEGER", "Number of the nodes column. Default = #{options[:node_index]}") do |node_i|
    options[:node_index] = node_i.to_i - 1
  end

  options[:cluster_index] = 0
  opts.on("-C", "--cluster_column INTEGER", "Number of the clusters column. Default = #{options[:cluster_index]}" ) do |cluster_i|
    options[:cluster_index] = cluster_i.to_i - 1
  end

  options[:column_sep] = "\t"
  opts.on("-S", "--split_char CHARACTER", "Character for splitting input file. Default: tab") do |split_char|
    options[:column_sep] = split_char
  end

  options[:node_sep] = nil
  opts.on("-s", "--node_sep CHARACTER", "Node split character. This option must to be used when input file is aggregated.") do |split_char|
    options[:node_sep] = split_char
  end

  options[:random_type] = ["size"]
  opts.on("-r", "--random_type STRING", "Indicate random mode. 'size' for radomize clusters with the same size as input clusters. 'full_size' same as 'size' but all nodes are repaeted as same as input. 'fixed:n:s' for generate 'n' clusters of 's' nodes. Default = #{options[:output_file]}") do |random_type|
    options[:random_type] = random_type.split(":")
  end

  options[:replacement] = false
  opts.on("-R", "--replacement", "Boolean. Activates ramdom sampling with replacement. Sampling witout replacement will be executed instead.") do 
    options[:replacement] = true
  end

  options[:output_file] = "./random_clusters.txt"
  opts.on("-o", "--output_file FILEPATH", "Output file") do |output_file|
    options[:output_file] = output_file
  end

  options[:aggregate_sep] = nil
  opts.on("-a", "--aggregate_sep CHARACTER", "This option activates aggregation in output. Separator character must be provided") do |split_char|
    options[:aggregate_sep] = split_char
  end

end.parse!
##########################
#MAIN
##########################

clusters = load_clusters(options)

nodes = clusters.values.flatten
nodes = nodes.uniq if !options[:random_type][0] == "full_size"

if options[:random_type][0].include?("size") && options[:random_type].size == 1
	all_sizes = clusters.map{|cluster, nodes| nodes.size}
elsif options[:random_type][0] == "fixed" && options[:random_type].size == 3
	all_sizes = Array.new(options[:random_type][1].to_i, options[:random_type][2].to_i)
end

random_clusters = random_sample(nodes, options[:replacement], all_sizes, 123)
write_clusters(random_clusters, options[:output_file], options[:aggregate_sep])
