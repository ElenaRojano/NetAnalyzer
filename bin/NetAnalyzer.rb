#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$LOAD_PATH.unshift(File.expand_path(File.join(ROOT_PATH, '..', 'lib')))
require 'optparse'
require 'benchmark'
require 'NetAnalyzer'

######################################
## METHODS
######################################
def load_file(path)
  data = []
  File.open(path).each do |line|
    data << line.chomp.split("\t")
  end
  return data
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

  options[:ontologies] = []
  opts.on("-O", "--ontology STRING", "String that define which ontologies must be used with each layer. String definition:'layer_name1:path_to_obo_file1;layer_name2:path_to_obo_file2'") do |ontologies|
    options[:ontologies] = ontologies.split(';').map{|pair| pair.split(':')}
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

  options[:graph_file] = nil
  opts.on("-g", "--graph_file PATH", "Build a graphic representation of the network") do |item|
    options[:graph_file] = item
  end

  options[:graph_options] = {method: 'el_grapho', layout: 'forcedir', steps: '30'}
  opts.on("--graph_options STRING", "Set graph parameters as 'NAME1=value1,NAME2=value2,...") do |item|
    options[:graph_options] = {}
    item.split(',').each do |pair|
      fields = pair.split('=')
      options[:graph_options][fields.first.to_sym] = fields.last 
    end
  end

  options[:threads] = 0
  opts.on( '-T', '--threads INTEGER', 'Number of threads to use in computation, one thread will be reserved as manager.' ) do |opt|
      options[:threads] = opt.to_i - 1
  end

  options[:reference_nodes] = []
  opts.on("-r", "--reference_nodes STRING", "Node ids comma separared") do |item|
    options[:reference_nodes] = item.split(',')
  end

  options[:group_nodes] = {}
  opts.on("-G", "--group_nodes STRING", "File path or groups separated by ';' and group node ids comma separared") do |item|
    if File.exists?(item)
      File.open(item).each do |line|
        groupID, nodeID = line.chomp.split("\t")
        query = options[:group_nodes][groupID]
        query.nil? ? options[:group_nodes][groupID] = [nodeID] : query << nodeID
      end
    else
      item.split(';').each_with_index do |group, i| 
        options[:group_nodes][i] = group.split(',')
      end
    end
  end

  options[:group_metrics] = false
  opts.on("-M", "--group_metrics", "Perform group group_metrics") do 
    options[:group_metrics] = true
  end

  options[:expand_clusters] = nil
  opts.on("-x", "--expand_clusters STRING", "Method to expand clusters Available methods: sht_path") do |item|
    options[:expand_clusters] = item
  end

  options[:get_attributes] = []
  opts.on("-A", "--attributes STRING", "String separadted by commas with the name of network attribute") do |item|
    options[:get_attributes] = item.split(',')
  end

  options[:delete_nodes] = []
  opts.on("-d", "--delete PATH", "Remove nodes from file. If PATH;r then nodes not included in file are removed") do |item|
    options[:delete_nodes] = item.split(';')
  end
end.parse!
##########################
#MAIN
##########################
puts "Loading network data"
fullNet = Net_parser.load(options)
fullNet.reference_nodes = options[:reference_nodes]
fullNet.threads = options[:threads]
fullNet.group_nodes = options[:group_nodes]
fullNet.set_compute_pairs(options[:use_pairs], !options[:no_autorelations])

if !options[:delete_nodes].empty?
  node_list = load_file(options[:delete_nodes].first).flatten
  options[:delete_nodes].length > 1 ? mode = options[:delete_nodes][1] : 'd'
  fullNet.delete_nodes(node_list, mode)
end

options[:ontologies].each do |layer_name, ontology_file_path|
  fullNet.link_ontology(ontology_file_path, layer_name.to_sym)
end

if !options[:get_attributes].empty?
  node_attributes = fullNet.get_node_attributes(options[:get_attributes])
  File.open(File.join(File.dirname(options[:output_file]), 'node_attributes.txt'), 'w' ) do |f|
    node_attributes.each do |attributes|
      f.puts(attributes.join("\t"))
    end
  end
end

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
  	Performancer.load_control(control)
    predictions = fullNet.association_values[options[:meth]]
  	performance = Performancer.get_pred_rec(predictions)
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

if !options[:graph_file].nil?
  options[:graph_options][:output_file] = options[:graph_file]
  fullNet.plot_network(options[:graph_options]) 
end

if options[:group_metrics]
  fullNet.compute_group_metrics(File.join(File.dirname(options[:output_file]), 'group_metrics.txt'))
end

if !options[:expand_clusters].nil?
  expanded_clusters = fullNet.expand_clusters(options[:expand_clusters])
  File.open(File.join(File.dirname(options[:output_file]), 'expand_clusters.txt'), 'w' ) do |f|
    expanded_clusters.each do |cl_id, nodes|
      nodes.each do |node|
        f.puts "#{cl_id}\t#{node}"
      end
    end
  end
end