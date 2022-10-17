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
puts "Loading network data"

if options[:layers].length == 1
  layerA = layerB = options[:layers][0].first
elsif  options[:layers].length == 2
  layerA = options[:layers][0].first
  layerB = options[:layers][1].first
end

if options[:input_format] == 'pair'
  fullNet.load_network_by_pairs(options[:input_file], options[:layers], options[:split_char])
elsif options[:input_format] == 'bin' && !options[:node_file].nil?
  fullNet.load_network_by_bin_matrix(options[:input_file], options[:node_file], options[:layers])
elsif options[:input_format] == 'matrix' && !options[:node_file].nil?
  fullNet.load_network_by_plain_matrix(options[:input_file], options[:node_file], options[:layers], options[:splitChar])
else
  raise("ERROR: The format #{options[:input_format]} is not defined")
  exit
end


fullNet.randomize_network(options[:type_random])


#fullNet.save_adjacency_matrix(layerA, layerB, options[:output_file])


