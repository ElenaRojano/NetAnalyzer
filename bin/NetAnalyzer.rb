#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'NetAnalyzer'))
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'NetAnalyzer', 'methods'))

require 'network'
require 'optparse'

##############################
# MAIN METHODS
##############################

def set_layer(layer_definitions, node_name)
	layer = nil
	if layer_definitions.length > 1
		layer_definitions.each do |layer_name, regexp|
			if node_name =~ regexp
				layer = layer_name
				break
			end
		end
	else
		layer = layer_definitions.first.first
	end
	return layer
end

##############################
#OPTPARSE
##############################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file to create bipartite networks for further analysis") do |input_file|
    options[:input_file] = input_file
  end

  options[:split_char] = "\t"
  opts.on("-s", "--split_char STRING", "Character for splitting input file. Default: tab") do |split_char|
    options[:split_char] = split_char
  end

  options[:output_file] = "network2plot"
  opts.on("-o", "--output_file PATH", "Output file name") do |output_file|
    options[:output_file] = output_file
  end

  options[:assoc_file] = "assoc_values.txt"
  opts.on("-a", "--assoc_file PATH", "Output file name for association values") do |output_file|
    options[:assoc_file] = output_file
  end

  options[:performance_file] = "perf_values.txt"
  opts.on("-p", "--performance_file PATH", "Output file name for performance values") do |output_file|
    options[:performance_file] = output_file
  end

  options[:layers] = [:layer, '-']
  opts.on("-l", "--layers STRING", "Layer definition on network: layer1name,regexp1;layer2name,regexp2...") do |layers|
	layers_definition = layers.split(";").map{|layer_attr| layer_attr.split(',')}
	layers_definition.map!{|layer_attr| [layer_attr.first.to_sym, /#{layer_attr.last}/]}
    options[:layers] = layers_definition
  end

  options[:use_layers] = []
  opts.on("-u", "--use_layers STRING", "Set which layers must be used on association methods: layer1,layer2;layerA,layerB") do |string|
    options[:use_layers] = string.split(";").map{|layer_attr| layer_attr.split(',').map{|layer| layer.to_sym}}
  end

  options[:control_file] = nil
  opts.on("-c", "--control_file PATH", "Control file name") do |file|
    options[:control_file] = file
  end

  options[:output_style] = "neato"
  opts.on("-t", "--output_style STRING", "Style to plot output file") do |output_style|
    options[:output_style] = output_style
  end

  options[:meth] = nil
  opts.on("-m", "--association_method STRING", "Association method to use on network") do |meth|
    options[:meth] = meth.to_sym
  end

  options[:no_autorelations] = FALSE
  opts.on("-N", "--no_autorelations", "Remove association values between nodes os same type") do
    options[:no_autorelations] = TRUE
  end

end.parse!

##########################
#MAIN
##########################

fullNet = Network.new(options[:layers].map{|layer| layer.first})
puts "Loading network data"
File.open(options[:input_file]).each("\n") do |line|
	line.chomp!
	pair = line.split(options[:splitChar])
	node1 = pair[0]
	node2 = pair[1]
	fullNet.add_node(node1, set_layer(options[:layers], node1))
	fullNet.add_node(node2, set_layer(options[:layers], node2))
	fullNet.add_edge(node1, node2)	
end
#fullNet.plot(options[:output_file], options[:output_style])


if !options[:meth].nil?
	puts "Performing association method #{options[:meth]} on network"
	if options[:meth] == :transference
		fullNet.generate_adjacency_matrix(options[:use_layers][0][0], options[:use_layers][0][1])
		fullNet.generate_adjacency_matrix(options[:use_layers][1][0], options[:use_layers][1][1])
		fullNet.get_association_values(
			[options[:use_layers][0][0], options[:use_layers][0][1]], 
			[options[:use_layers][1][0], options[:use_layers][1][1]],
			:transference)
	else
		fullNet.get_association_values(
			options[:use_layers][0], 
			options[:use_layers][1].first, 
			options[:meth])
	end
	puts 'Clean autorelations'  if options[:no_autorelations]
	fullNet.clean_autorelations_on_association_values if options[:no_autorelations]
	File.open(options[:assoc_file], 'w') do |f|
		fullNet.association_values[options[:meth]].each do |val|
			f.puts val.join("\t")
		end
	end
end

if !options[:meth].nil? && !options[:control_file].nil?
	puts "Doing validation on association values obtained from method #{options[:meth]}"
	control = []
	File.open(options[:control_file]).each("\n") do |line|
		line.chomp!
		control << line.split("\t")
	end
	fullNet.load_control(control)
	performance = fullNet.get_pred_rec(options[:meth])
	File.open(options[:performance_file], 'w') do |f|
		f.puts %w[cut prec rec meth].join("\t")
		performance.each do |item|
			item << options[:meth].to_s
			f.puts item.join("\t")
		end
	end
end

puts "End of analysis: #{options[:meth]}"
