#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'NetAnalyzer'))
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'NetAnalyzer', 'methods'))

require 'network'
require 'optparse'

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

  options[:node_file] = nil
  opts.on("-n", "--node_names_file PATH", "File with node names corresponding to the input matrix, only use when -i is set to bin or matrix.") do |node_file|
    options[:node_file] = node_file
  end

  options[:input_format] = 'pair'
  opts.on("-f", "--input_format STRING", "Input file format: pair (default), bin, matrix") do |input_format|
    options[:input_format] = input_format
  end

  options[:split_char] = "\t"
  opts.on("-s", "--split_char STRING", "Character for splitting input file. Default: tab") do |split_char|
    options[:split_char] = split_char
  end

  options[:use_pairs] = :conn
  opts.on("-P", "--use_pairs STRING", "Which pairs must be computed. 'all' means all posible pair node combinations and 'conn' means the pair are truly connected in the network. Default 'conn' ") do |use_pairs|
    options[:use_pairs] = use_pairs.to_sym
  end

  options[:output_file] = "network2plot"
  opts.on("-o", "--output_file PATH", "Output file name") do |output_file|
    options[:output_file] = output_file
  end

  options[:assoc_file] = "assoc_values.txt"
  opts.on("-a", "--assoc_file PATH", "Output file name for association values") do |output_file|
    options[:assoc_file] = output_file
  end

  options[:kernel_file] = "kernel_values"
  opts.on("-K", "--kernel_file PATH", "Output file name for kernel values") do |output_file|
    options[:kernel_file] = output_file
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

  options[:kernel] = nil
  opts.on("-k", "--kernel_method STRING", "Kernel operation to perform with the adjacency matrix") do |kernel|
    options[:kernel] = kernel
  end

  options[:no_autorelations] = false
  opts.on("-N", "--no_autorelations", "Remove association values between nodes os same type") do
    options[:no_autorelations] = true
  end

  options[:normalize_kernel] = false
  opts.on("-z", "--normalize_kernel_values", "Apply cosine normalization to the obtained kernel") do
    options[:normalize_kernel] = true
  end

  options[:byte_format] = :float64
  opts.on( '-b', '--byte_format STRING', 'Format of the numeric values stored in matrix. Default: float64, warning set this to less precission can modify computation results using this matrix.' ) do |opt|
      options[:byte_format] = opt.to_sym
  end
end.parse!

##########################
#MAIN
##########################

fullNet = Network.new(options[:layers].map{|layer| layer.first})
fullNet.set_compute_pairs(options[:use_pairs])
fullNet.set_matrix_byte_format(options[:byte_format])
#puts options[:layers].map{|layer| layer.first}.inspect
puts "Loading network data"
if options[:input_format] == 'pair'
  fullNet.load_network_by_pairs(options[:input_file], options[:layers], options[:splitChar])
elsif options[:input_format] == 'bin'
  fullNet.load_network_by_bin_matrix(options[:input_file], options[:node_file], options[:layers])
elsif options[:input_format] == 'matrix'
  fullNet.load_network_by_plain_matrix(options[:input_file], options[:node_file], options[:layers], options[:splitChar])
else
  raise("ERROR: The format #{options[:input_format]} is not defined")
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
  if options[:no_autorelations]
  	puts 'Clean autorelations'  
  	fullNet.clean_autorelations_on_association_values
  end
	File.open(options[:assoc_file], 'w') do |f|
		fullNet.association_values[options[:meth]].each do |val|
			f.puts val.join("\t")
		end
	end
  if !options[:control_file].nil?
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
end

if !options[:kernel].nil?
  layer2kernel = options[:use_layers].first # we use only a layer to perform the kernel, so only one item it is selected.
  fullNet.get_kernel(layer2kernel, options[:kernel], options[:normalize_kernel])
  fullNet.write_kernel(layer2kernel, options[:kernel_file])
end
