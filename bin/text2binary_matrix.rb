#! /usr/bin/env ruby

require 'benchmark'
require 'optparse'
require 'nmatrix'
require 'pp'
#############################################################################
## METHODS
##############################################################################

def load_matrix_file(source, byte_format = :float32)
    dimension_elements = 0
    adjacency_vector = []
    source.each do |line|
            line.chomp!
            adjacency_vector.concat(line.split("\t").map{|c| c.to_f })
            dimension_elements += 1
    end
    matrix = NMatrix.new([dimension_elements, dimension_elements], adjacency_vector, dtype: byte_format) # Create working matrix
    return matrix
end

def load_pair_file(source, byte_format = :float32)
	connections = {}
	source.each do |line|
		node_a, node_b, weight = line.chomp.split("\t")
		weight.nil? ? weight = 1.0 : weight = weight.to_f
		add_pair(node_a, node_b, weight, connections)
		add_pair(node_b, node_a, weight, connections)
	end
	names = connections.keys
	matrix = NMatrix.new( names.length, 0.0, dtype: byte_format)
	count = 0
	connections.each do |nodeA, subhash|
		index_A = names.index(nodeA)
		subhash.each do |nodeB, weight|
			index_B = names.index(nodeB)
			matrix[index_A, index_B] = weight
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

def get_stats(matrix)
	stats = []
	primary_stats = get_primary_stats(matrix)
	stats << ['Matrix - Symmetric?', matrix.symmetric?]
	stats << ['Matrix - Dimensions', matrix.shape.join('x')]
	stats << ['Matrix - Elements', primary_stats[:count]]
	stats << ['Matrix - Elements Non Zero', primary_stats[:countNonZero]]
	stats << ['Matrix - Non Zero Density', primary_stats[:countNonZero].fdiv(primary_stats[:count])]
	stats << ['Weigth - Max', primary_stats[:max]]
	stats << ['Weigth - Min', primary_stats[:min]]
	stats << ['Weigth - Average', primary_stats[:average]]
	stats << ['Weigth - Variance', primary_stats[:variance]]
	stats << ['Weigth - Standard Deviation', primary_stats[:standardDeviation]]
	stats << ['Weigth - Q1', primary_stats[:q1]]
	stats << ['Weigth - Median', primary_stats[:median]]
	stats << ['Weigth - Q3', primary_stats[:q3]]
	stats << ['Weigth - Min Non Zero', primary_stats[:minNonZero]]
	stats << ['Weigth - Average Non Zero', primary_stats[:averageNonZero]]
	stats << ['Weigth - Variance Non Zero', primary_stats[:varianceNonZero]]
	stats << ['Weigth - Standard Deviation Non Zero', primary_stats[:standardDeviationNonZero]]
	stats << ['Weigth - Q1 Non Zero', primary_stats[:q1NonZero]]
	stats << ['Weigth - Median Non Zero', primary_stats[:medianNonZero]]
	stats << ['Weigth - Q3 Non Zero', primary_stats[:q3NonZero]]
	connections = get_connection_number(matrix)
	connection_stats = get_primary_stats(connections)
	stats << ['Node - Elements', connection_stats[:count]]
	stats << ['Node - Elements Non Zero', connection_stats[:countNonZero]]
	stats << ['Node - Non Zero Density', connection_stats[:countNonZero].fdiv(connection_stats[:count])]
	stats << ['Edges - Max', connection_stats[:max]]
	stats << ['Edges - Min', connection_stats[:min]]
	stats << ['Edges - Average', connection_stats[:average]]
	stats << ['Edges - Variance', connection_stats[:variance]]
	stats << ['Edges - Standard Deviation', connection_stats[:standardDeviation]]
	stats << ['Edges - Q1', connection_stats[:q1]]
	stats << ['Edges - Median', connection_stats[:median]]
	stats << ['Edges - Q3', connection_stats[:q3]]
	stats << ['Edges - Min Non Zero', primary_stats[:minNonZero]]
	stats << ['Edges - Average Non Zero', connection_stats[:averageNonZero]]
	stats << ['Edges - Variance Non Zero', connection_stats[:varianceNonZero]]
	stats << ['Edges - Standard Deviation Non Zero', connection_stats[:standardDeviationNonZero]]
	stats << ['Edges - Q1 Non Zero', connection_stats[:q1NonZero]]
	stats << ['Edges - Median Non Zero', connection_stats[:medianNonZero]]
	stats << ['Edges - Q3 Non Zero', connection_stats[:q3NonZero]]	
	return stats
end

def get_connection_number(matrix)
	connections = NMatrix.new([1, matrix.cols], 0, dtype: :float32)
	i = 0
	matrix.each_column do |column|
		count = 0
		column.each do |value|
			count += 1 if value != 0
		end
		connections[0, i] = count - 1
		i += 1
	end
	return connections
end

def transform_keys(hash)
	new_hash = {}
	hash.each do |key, val|
		new_key = yield(key)
		new_hash[new_key] = val
	end
	return new_hash
end

def get_primary_stats(matrix)
	stats = Hash.new(0)
	max = matrix[0, 0] # Initialize max value
	min = matrix[0, 0] # Initialize min value
	min_non_zero = matrix[0, 0] # Initialize min value
	matrix.each do |value|
		stats[:count] += 1
		stats[:countNonZero] += 1 if value != 0
		stats[:sum] += value
		max = value if value > max
		min = value if value < min
		min_non_zero = value if value != 0 && value < min
	end
	stats[:max] = max
	stats[:min] = min
	stats[:minNonZero] = min_non_zero
	values = matrix.to_a
	values.flatten! if values.first.class == Array
	values.sort!
	quartile_stats = get_quartiles(values, stats[:count])
	stats.merge!(transform_keys(quartile_stats){|k| k.to_sym})
	values.select!{|v| v != 0}
	quartile_stats_non_zero = get_quartiles(values, stats[:countNonZero])
	stats.merge!(transform_keys(quartile_stats_non_zero){|k| (k + 'NonZero').to_sym})
	get_composed_stats(stats, matrix)
	return stats
end

def get_quartiles(values, n_items)
	stats = {}
	q1_coor = n_items * 0.25 - 1
	median = n_items * 0.5 - 1
	q3_coor = n_items * 0.75 - 1
	if n_items % 2 == 0
		stats['q1'] = (values[q1_coor.to_i] + values[q1_coor.to_i + 1]).fdiv(2)
		stats['median'] = (values[median.to_i] + values[median.to_i + 1]).fdiv(2)
		stats['q3'] = (values[q3_coor.to_i] + values[q3_coor.to_i + 1]).fdiv(2)		
	else
		stats['q1'] = values[q1_coor.ceil]
		stats['median'] = values[median.ceil]
		stats['q3'] = values[q3_coor.ceil]
	end
	return stats
end

def get_composed_stats(stats, matrix)
	average = stats[:sum].fdiv(stats[:count])
	average_non_zero = stats[:sum].fdiv(stats[:countNonZero])
	stats[:average] = average
	stats[:averageNonZero] = average_non_zero
	matrix.each do |value|
		stats[:sumDevs] = (value - average) ** 2
		stats[:sumDevsNonZero] = (value - average_non_zero) ** 2 if value != 0
	end
	stats[:variance] = stats[:sumDevs].fdiv(stats[:count])
	stats[:varianceNonZero] = stats[:sumDevsNonZero].fdiv(stats[:countNonZero])
	stats[:standardDeviation] = stats[:variance] ** 0.5
	stats[:standardDeviationNonZero] = stats[:varianceNonZero] ** 0.5
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

    options[:byte_format] = :float64
    opts.on( '-b', '--byte_format STRING', 'Format of the numeric values stored in matrix. Default: float64, warning set this to less precission can modify computation results using this matrix.' ) do |opt|
        options[:byte_format] = opt.to_sym
    end

    options[:input_type] = 'pair'
    opts.on( '-t', '--input_type STRING', 'Set input format file. "pair" or "matrix"' ) do |opt|
        options[:input_type] = opt
    end

    options[:set_diagonal] = false
    opts.on( '-d', '--set_diagonal', 'Set to 1.0 the main diagonal' ) do 
        options[:set_diagonal] = true
    end

    options[:stats] = false
    opts.on( '-s', '--get_stats', 'Get stats from the processed matrix' ) do
        options[:stats] = true
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

if options[:input_type] == 'bin'
	matrix = NMatrix.read(options[:input_file]) # the method needs a path not a IO object 
elsif options[:input_type] == 'matrix'
	matrix = load_matrix_file(source, options[:byte_format])
elsif options[:input_type] == 'pair'
	matrix, names = load_pair_file(source, options[:byte_format])
	File.open(options[:output_matrix_file]+'.lst', 'w'){|f| f.print names.join("\n")}
end

if options[:set_diagonal]
	elements = matrix.cols
	elements.times do |n|
		matrix[n, n] = 1.0
	end
end

if options[:stats]
	stats = get_stats(matrix)
	stats.each do |stat|
		puts stat.join("\t")
	end
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
