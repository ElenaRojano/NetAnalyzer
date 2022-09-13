require 'gv'

#For javascrip plotting
require 'erb'
require 'base64'
require 'json'
require 'zlib'

TEMPLATES = File.join(File.dirname(__FILE__), 'templates')


class Net_plotter
	def initialize(net_data, options = {})
		@group_nodes = net_data[:group_nodes]
		@reference_nodes = net_data[:reference_nodes]
		@nodes = net_data[:nodes]
		@edges = net_data[:edges]
		@layers = net_data[:layers]

		if options[:method] == 'graphviz'
			plot_dot(options)
		elsif options[:method] == 'cyt_app'
			plot_cyt_app(options)
		else
			if options[:method] == 'elgrapho'
				template = 'el_grapho'
			elsif options[:method] == 'cytoscape'
				template = 'cytoscape'
			elsif options[:method] == 'sigma'
				template = 'sigma'
			end
			renderered_template = ERB.new(File.open(File.join(TEMPLATES, template + '.erb')).read).result(binding)
			File.open(options[:output_file] + '.html', 'w'){|f| f.puts renderered_template}
		end	
	end

	def plot_cyt_app(user_options = {})
		options = {}
		options = options.merge(user_options)

		node_cyt_ids = {}
		nodes = []
		count = 0
		group_nodes = {}
		@group_nodes.each do |groupID, gNodes|
			gNodes.each do |gNode|
				group_nodes[gNode] = groupID
			end
		end
		@nodes.each do |id, node|
			cyt_app_add_node(nodes, count, node, group_nodes)
			node_cyt_ids[id] = count.to_s
			count += 1
		end
		edges = cyt_app_add_edges(node_cyt_ids, count)
		cys_net = {
			'elements' => {
				'nodes' => nodes,
				'edges' => edges
			}
		}
		File.open(options[:output_file]+ '.cyjs', 'w'){|f| f.print JSON.pretty_generate(cys_net)}
	end

	def cyt_app_add_node(nodes, count, node, group_nodes)
		id = node.id
		cyt_node = {
			'data' => {
				'id' => count.to_s,
				'name' => id
			}
		}
		cyt_node['data']['type'] = node.type
		if !@reference_nodes.empty?
			ref = @reference_nodes.include?(id) ? 'y' : 'n'
			cyt_node['data']['ref'] = ref
		end
		if !group_nodes.empty?
			query = group_nodes[id]
			cyt_node['data']['group'] = query if !query.nil?			
		end
		nodes << cyt_node
	end

	def cyt_app_add_edges(node_ids, count)
		edges = []
		plotted_edges = {}
		@edges.each do |nodeID, associatedIDs|
			associatedIDs.each do |associatedID|
				pair = [nodeID, associatedID].sort.join('_').to_sym
				if !plotted_edges[pair]
					edges << {
						'data' => {
							'id' => count.to_s,
							'source' => node_ids[nodeID],
							'target' => node_ids[associatedID],
							"interaction" => "-",
							"weight" => 1.0
						}
					}
					count +=1
					plotted_edges[pair] = true
				end
			end
		end
		return edges
	end

	def plot_dot(user_options = {}) # input keys: layout
		options = {layout: "sfdp"}
		options = options.merge(user_options)
		graphviz_colors = %w[lightsteelblue1 lightyellow1 lightgray orchid2]
		palette = {}
		@layers.each do |layer|
			palette[layer] = graphviz_colors.shift
		end
		graph = GV::Graph.open('g', type = :undirected)
		plotted_edges = {}
		@edges.each do |nodeID, associatedIDs|
			associatedIDs.each do |associatedID|
				pair = [nodeID, associatedID].sort.join('_').to_sym
				if !plotted_edges[pair]
					graph.edge 'e', 
						graph.node(nodeID, label: '', style: 'filled', fillcolor: palette[@nodes[nodeID].type]), 
						graph.node(associatedID, label: '', style: 'filled' , fillcolor: palette[@nodes[associatedID].type])
					plotted_edges[pair] = true
				end
			end
		end
		@reference_nodes.each do |nodeID|
			graph.node(nodeID, style: 'filled', fillcolor: 'firebrick1', label: '')
		end
		graphviz_border_colors = %w[blue darkorange red olivedrab4]
		@group_nodes.each do |groupID, gNodes|
			border_color = graphviz_border_colors.shift
			gNodes.each do |nodeID|
					graph.node(nodeID, color: border_color, penwidth: '10', label: '')
			end
		end
		graph[:overlap] = false
		STDERR.puts 'Save graph'
		graph.save(options[:output_file] + '.png', format='png', layout=options[:layout])
	end

end