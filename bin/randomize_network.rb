#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$LOAD_PATH.unshift(File.expand_path(File.join(ROOT_PATH, '..', 'lib')))
require 'optparse'
require 'benchmark'
require 'NetAnalyzer'



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

  options[:node_file] = nil
  opts.on("-n", "--node_names_file PATH", "File with node names corresponding to the input matrix, only use when -i is set to bin or matrix.") do |node_file|
    options[:node_file] = node_file
  end

  options[:input_format] = 'pair'
  opts.on("-f", "--input_format STRING", "Input file format: pair (default), bin, matrix") do |input_format|
    options[:input_format] = input_format
  end

  options[:split_char] = "\t"
  opts.on("-s", "--split_char CHARACTER", "Character for splitting input file. Default: tab") do |split_char|
    options[:split_char] = split_char
  end

  options[:layers] = [:layer, '-']
  opts.on("-l", "--layers STRING", "Layer definition on network: layer1name,regexp1;layer2name,regexp2...") do |layers|
  	layers_definition = layers.split(";").map{|layer_attr| layer_attr.split(',')}
	  layers_definition.map!{|layer_attr| [layer_attr.first.to_sym, /#{layer_attr.last}/]}
    options[:layers] = layers_definition
  end

  options[:type_random] = nil
  opts.on("-r", "--type_random network", "Randomized basis. 'nodes' Node-baseds randomize or 'links' Links-baseds randomize") do |type_random|
    options[:type_random] = type_random
  end

  options[:output_file] = nil
  opts.on("-o", "--output_file FILEPATH", "Output file") do |output_file|
    options[:output_file] = output_file
  end

end.parse!


##########################
#MAIN
##########################
fullNet = Network.new(options[:layers].map{|layer| layer.first})
#fullNet.threads = options[:threads]
#fullNet.group_nodes = options[:group_nodes]
#puts options[:layers].map{|layer| layer.first}.inspect
puts "Loading network data"

if options[:layers].length == 1
  layerA = layerB = options[:layers][0].first
elsif  options[:layers].length == 2
  layerA = options[:layers][0].first
  layerB = options[:layers][1].first
end

if options[:input_format] == 'pair'
  fullNet.load_network_by_pairs(options[:input_file], options[:layers], options[:split_char])
  fullNet.generate_adjacency_matrix(layerA, layerB)
elsif options[:input_format] == 'bin' && !options[:node_file].nil?
  fullNet.load_network_by_bin_matrix(options[:input_file], options[:node_file], options[:layers])
elsif options[:input_format] == 'matrix' && !options[:node_file].nil?
  fullNet.load_network_by_plain_matrix(options[:input_file], options[:node_file], options[:layers], options[:splitChar])
else
  raise("ERROR: The format #{options[:input_format]} is not defined")
  exit
end


if options[:type_random] == 'nodes'
  if options[:layers].length == 1
    fullNet.randomize_monopartite_net_by_nodes(layer)
    #fullNet.build_edges_from_adjacency_matrix([layer])
  elsif options[:layers].length == 2
    fullNet.randomize_bipartite_net_by_nodes(layerA, layerB)
    #fullNet.build_edges_from_adjacency_matrix([layerA, layerB])

  else 
    raise("ERROR: The randomization is not available for #{options[:layers].length} types of nodes")
    exit
  end
elsif options[:type_random] == 'links'
  if options[:layers].length == 1
    fullNet.randomize_monopartite_net_by_links([layer])
    #fullNet.build_nodes_from_adjacency_matrix(options[:layers], [layer])
  elsif options[:layers].length == 2
    fullNet.randomize_bipartite_net_by_links([layerA, layerB])
    #fullNet.build_nodes_from_adjacency_matrix(options[:layers], [layerA, layerB])
  else
    raise("ERROR: The randomization is not available for #{options[:layers].length} types of nodes")
    exit
  end
else
  raise("ERROR: The #{options[:type_random]}-based random is not defined")
  exit
end

fullNet.save_adjacency_matrix(layerA, layerB, options[:output_file])


