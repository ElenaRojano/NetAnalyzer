require 'numo/narray'
require 'numo/linalg'

class Adv_mat_calc # Advanced matrix calculations
	## KERNEL METHODS
	#######################################################################################
	def self.get_kernel(matrix, node_names, kernel, normalization=false)
		#I = identity matrix
		#D = Diagonal matrix
		#A = adjacency matrix
		#L = laplacian matrix = D − A
		matrix_result = nil
		dimension_elements = matrix.shape.last
		# In scuba code, the diagonal values of A is set to 0. In weighted matrix the kernel result is the same with or without this operation. Maybe increases the computing performance?
		# In the md kernel this operation affects the values of the final kernel
		#dimension_elements.times do |n|
		#	matrix[n,n] = 0.0
		#end
		if kernel == 'el' || kernel == 'ct' || kernel == 'rf' || 
			kernel.include?('vn') || kernel.include?('rl') || kernel == 'me'
			diagonal_matrix = matrix.sum(1).diag 	# get the total sum for each row, for this reason the sum method takes the 1 value. If sum colums is desired, use 0
													# Make a matrix whose diagonal is row_sum
			matrix_L = diagonal_matrix - matrix
			if kernel == 'el' #Exponential Laplacian diffusion kernel(active). F Fouss 2012 | doi: 10.1016/j.neunet.2012.03.001
			    beta = 0.02
			    beta_product = matrix_L * -beta
			    #matrix_result = beta_product.expm
			    matrix_result = Numo::Linalg.expm(beta_product, 14)
			elsif kernel == 'ct' # Commute time kernel (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			    matrix_result = Numo::Linalg.pinv(matrix_L) # Anibal saids that this kernel was normalized. Why?. Paper do not seem to describe this operation for ct, it describes for Kvn or for all kernels, it is not clear.
			elsif kernel == 'rf' # Random forest kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			    matrix_result = Numo::Linalg.inv(Numo::DFloat.eye(dimension_elements) + matrix_L) #Krf = (I +L ) ^ −1
			elsif kernel.include?('vn') # von Neumann diffusion kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			    alpha = kernel.gsub('vn', '').to_f * matrix.max_eigenvalue ** -1  # alpha = impact_of_penalization (1, 0.5 or 0.1) * spectral radius of A. spectral radius of A = absolute value of max eigenvalue of A 
			    matrix_result = Numo::Linalg.inv(Numo::DFloat.eye(dimension_elements) - matrix * alpha ) #  (I -alphaA ) ^ −1
			elsif kernel.include?('rl') # Regularized Laplacian kernel matrix (active)
			    alpha = kernel.gsub('rl', '').to_f * matrix.max_eigenvalue ** -1  # alpha = impact_of_penalization (1, 0.5 or 0.1) * spectral radius of A. spectral radius of A = absolute value of max eigenvalue of A
			    matrix_result = Numo::Linalg.inv(Numo::DFloat.eye(dimension_elements) + matrix_L * alpha ) #  (I + alphaL ) ^ −1
			elsif kernel == 'me' # Markov exponential diffusion kernel (active). G Zampieri 2018 | doi.org/10.1186/s12859-018-2025-5 . Taken from compute_kernel script
				beta=0.04
				#(beta/N)*(N*I - D + A)
				id_mat = Numo::DFloat.eye(dimension_elements)
				m_matrix = (id_mat * dimension_elements - diagonal_matrix + matrix ) * (beta/dimension_elements)
				#matrix_result = m_matrix.expm
			    matrix_result = Numo::Linalg.expm(m_matrix, 16)
			end
		elsif kernel == 'ka' # Kernelized adjacency matrix (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
			lambda_value = matrix.min_eigenvalue
			matrix_result = matrix + Numo::DFloat.eye(dimension_elements) * lambda_value.abs # Ka = A + lambda*I # lambda = the absolute value of the smallest eigenvalue of A
		elsif kernel.include?('md') # Markov diffusion kernel matrix. G Zampieri 2018 | doi.org/10.1186/s12859-018-2025-5 . Taken from compute_kernel script
			t = kernel.gsub('md', '').to_i
			#TODO: check implementation with Numo::array
			col_sum = matrix.sum(1)
			p_mat = matrix.div_by_vector(col_sum)
			p_temp_mat = p_mat.clone
			zt_mat = p_mat.clone
			(t-1).times do
				p_temp_mat = p_temp_mat.dot(p_mat)
				zt_mat = zt_mat + p_temp_mat
			end
			zt_mat = zt_mat * (1.0/t)
			matrix_result = zt_mat.dot(zt_mat.transpose)
		else
			matrix_result = matrix
			warn('Warning: The kernel method was not specified or not exists. The adjacency matrix will be given as result')
			# This allows process a previous kernel and perform the normalization in a separated step.
		end
		matrix_result = matrix_result.cosine_normalization if normalization #TODO: check implementation with Numo::array
		return matrix_result
	end

	# Alaimo 2014, doi: 10.3389/fbioe.2014.00071
	def self.tranference_resources(matrix1, matrix2, lambda_value1 = 0.5, lambda_value2 = 0.5)
		m1rowNumber, m1colNumber = matrix1.shape
		m2rowNumber, m2colNumber = matrix2.shape
		#puts m1rowNumber, m1colNumber, m2rowNumber, m2colNumber
		matrix1Weight = self.graphWeights(m1colNumber, m1rowNumber, matrix1.transpose, lambda_value1)
		matrix2Weight = self.graphWeights(m2colNumber, m2rowNumber, matrix2.transpose, lambda_value2)
		matrixWeightProduct = Numo::Linalg.dot(matrix1Weight, Numo::Linalg.dot(matrix2, matrix2Weight))
		finalMatrix = Numo::Linalg.dot(matrix1, matrixWeightProduct)
		return finalMatrix
	end

	def self.graphWeights(rowsNumber, colsNumber, inputMatrix, lambdaValue = 0.5)
	 	ky = (1.0 / inputMatrix.sum(0)).diag #sum cols
	 	weigth = Numo::Linalg.dot(inputMatrix, ky).transpose
	 	ky = nil #free memory
	 	weigth = Numo::Linalg.dot(inputMatrix, weigth)

	 	kx = inputMatrix.sum(1) #sum rows
	 	
	 	kx_lamb = kx ** lambdaValue
	 	kx_lamb_mat = Numo::DFloat.zeros(rowsNumber, rowsNumber)
	 	rowsNumber.times do |j|
	 		rowsNumber.times do |i|	 			
	 			kx_lamb_mat[j,i] = kx_lamb[i]
	 		end
	 	end
	 	kx_lamb = nil #free memory

	 	kx_inv_lamb = kx ** (1 - lambdaValue)
	 	kx_inv_lamb_mat = Numo::DFloat.zeros(rowsNumber, rowsNumber)
	 	rowsNumber.times do |j|
	 		rowsNumber.times do |i|
	 			kx_inv_lamb_mat[i, j] = kx_inv_lamb[i]
	 		end
	 	end
	 	kx_inv_lamb = nil #free memory

	 	nx = 1.0/(kx_lamb_mat.inplace * kx_inv_lamb_mat).inplace # inplace marks a matrix to be used by reference, not for value
	 	kx_lamb_mat = nil #free memory
	 	kx_inv_lamb_mat = nil #free memory
	 	weigth.inplace * nx
	 	return weigth
	end

end
