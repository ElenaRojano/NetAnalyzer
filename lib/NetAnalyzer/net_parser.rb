require 'expcalc'

class Net_parser
	def self.load(options)
		net = nil
		if options[:input_format] == 'pair'
		  net = load_network_by_pairs(options[:input_file], options[:layers], options[:split_char])
		elsif options[:input_format] == 'bin'
		  net = load_network_by_bin_matrix(options[:input_file], options[:node_file], options[:layers])
		elsif options[:input_format] == 'matrix'
		  net = load_network_by_plain_matrix(options[:input_file], options[:node_file], options[:layers], options[:splitChar])
		else
		  raise("ERROR: The format #{options[:input_format]} is not defined")
		end
		return net
	end

	def self.load_network_by_pairs(file, layers, split_character="\t")
		net = Network.new(layers.map{|layer| layer.first})
		File.open(file).each do |line|
			line.chomp!
			pair = line.split(split_character)
			node1 = pair[0]
			node2 = pair[1]
			net.add_node(node1, net.set_layer(layers, node1))
			net.add_node(node2, net.set_layer(layers, node2))
			net.add_edge(node1, node2)	
		end
		return net
	end

	def self.load_network_by_bin_matrix(input_file, node_file, layers)
		net = Network.new(layers.map{|layer| layer.first})
		node_names = load_input_list(node_file)
		net.adjacency_matrices[layers.map{|l| l.first}] = [Numo::NArray.load(input_file, type='npy'), node_names, node_names]
		return net
	end

	def self.load_network_by_plain_matrix(input_file, node_file, layers, splitChar="\t")
		net = Network.new(layers.map{|layer| layer.first})
		node_names = load_input_list(node_file)
		net.adjacency_matrices[layers.map{|l| l.first}] = [Numo::NArray.load(input_file, type='txt', splitChar=splitChar), node_names, node_names]
		return net
	end

	private
	def self.load_input_list(file)
		return File.open(file).readlines.map!{|line| line.chomp}
	end
end