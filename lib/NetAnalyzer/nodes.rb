class Node
	attr_reader :type, :id
	def initialize(id, type)
		@id = id
		@type = type
	end

	def clone
		node_clone = Node.new(@id.clone, @type.clone)
		return node_clone
	end

	def ==(other)
		are_equal = true
		if self.id != other.id ||
		   self.type != other.type	
		   are_equal = false
		end   
		return are_equal
	end

end