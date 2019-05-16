#! /usr/bin/env ruby

require 'benchmark'
require 'optparse'
require 'nmatrix'
require 'pp'
#############################################################################
## METHODS
##############################################################################

def load_matrix_file(source)
        dimension_elements = 0
        adjacency_vector = []
        source.each do |line|
                line.chomp!
                adjacency_vector.concat(line.split("\t").map{|c| c.to_f })
                dimension_elements += 1
        end
        matrix = NMatrix.new([dimension_elements, dimension_elements], adjacency_vector, dtype: :float32) # Create working matrix
        return matrix
end

def load_pair_file(source)
	connections = {}
	source.each do |line|
		node_a, node_b, weight = line.chomp.split("\t")
		weight.nil? ? weight = 1.0 : weight = weight.to_f
		add_pair(node_a, node_b, weight, connections)
		add_pair(node_b, node_a, weight, connections)
	end
	names = connections.keys
	matrix = NMatrix.new( names.length, 0.0, dtype: :float32)
	count = 0
	connections.each do |nodeA, subhash|
		subhash.each do |nodeB, weight|
			index_B = names.index(nodeB)
			matrix[count, index_B] = weight
		end
		count += 1
	end	
	return matrix, names
end

def add_pair(node_a, node_b, weight, connections)
	query = connections[node_a]
	if !query.nil?
		query[node_b] = weight
	else
		subhash = Hash.new(0.0)
		subhash[node_b] = weight
		connections[node_a] = subhash
	end
end

#############################################################################
## OPTPARSE
##############################################################################
options = {}

optparse = OptionParser.new do |opts|
    options[:input_file] = nil
    opts.on( '-i', '--input_file PATH', 'Input file' ) do |opt|
        options[:input_file] = opt
    end

    options[:output_matrix_file] = nil
    opts.on( '-o', '--output_matrix_file PATH', 'Output matrix file' ) do |opt|
        options[:output_matrix_file] = opt
    end

    options[:input_type] = 'pair'
    opts.on( '-t', '--input_type PATH', 'Set input format file. "pair" or "matrix"' ) do |opt|
        options[:input_type] = opt
    end

    options[:output_type] = 'bin'
    opts.on( '-O', '--output_type PATH', 'Set output format file. "bin" for binary (default) or "mat" for tabulated text file matrix' ) do |opt|
        options[:output_type] = opt
    end

   opts.banner = "Usage: #{File.basename(__FILE__)} [options] \n\n"

    opts.on( '-h', '--help', 'Display this screen' ) do
            puts opts
            exit
    end
end

optparse.parse!

################################################################################
## MAIN
###############################################################################
if options[:input_file] == '-'
	source = STDIN
else
	source = File.open(options[:input_file])
end

if options[:input_type] == 'matrix'
	matrix = load_matrix_file(source)
elsif options[:input_type] == 'pair'
	matrix, names = load_pair_file(source)
	File.open(options[:output_matrix_file]+'.lst', 'w'){|f| f.print names.join("\n")}
end
if options[:output_type] == 'bin'	
	matrix.write(options[:output_matrix_file])
elsif options[:output_type] == 'mat'
	File.open(options[:output_matrix_file], 'w') do |f|
		matrix.each_row do |r|
			f.puts r.to_a.join("\t")
		end
	end
end
