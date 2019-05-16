require 'nmatrix'
require 'nmatrix/lapacke'
require 'cmath'

class NMatrix
	def div(second_mat) #Matrix division A/B => A.dot(B.pinv) #https://stackoverflow.com/questions/49225693/matlab-matrix-division-into-python
		return self.dot(second_mat.pinv)
	end

	def div_by_vector(vector, by=:col)
		new_matrix =  NMatrix.zeros(self.shape, dtype: :float32)
		if by == :col
			self.cols.times do |n|
				vector.each_with_indices do |val, i, j|
					new_matrix[i, n] = self[i, n].fdiv(val)
				end
			end
		elsif by == :row

		end
		return new_matrix
	end

	def frobenius_norm
		fro = 0.0
		self.each do |value|
			fro += value ** 2
		end
		return fro ** 0.5
	end

	def vector_product(vec_b)
		product = 0.0
		self.each_with_indices do |val, i, j|
			product += val * vec_b[i, j]
		end
		return product
	end

	def vector_self_product
		product = 0.0
		self.each_stored_with_indices do |val, i, j|
			product +=  val ** 2
		end
		return product
	end

	def max_eigenvalue(n=100, error = 10e-12) # do not set error too low or the eigenvalue cannot stabilised around the real one
		max_eigenvalue = 0.0
		length = self.cols
		v = NMatrix.random([self.cols, 1])
		# http://web.mit.edu/18.06/www/Spring17/Power-Method.pdf 
		#IMPLEMENTATION PROBLEM: RESULTS ARE TOO VARIABLE
		last_max_eigenvalue = nil
		n.times do
			v = self.dot(v) # calculate the matrix-by-vector product Mv
			v = v / v.frobenius_norm # calculate the norm and normalize the vector
			max_eigenvalue = v.vector_product(self.dot(v)) / v.vector_self_product #Rayleigh quotient 
			# Rayleigh quotient: lambda = vMv/vv
			# v is a vector so vv is inner product of one vector with self (use vector_self_product); 
			# Mv gives a vector, so vMv is the inner product of two different vectors (use vector_product)
			break if !last_max_eigenvalue.nil? && last_max_eigenvalue - max_eigenvalue <= error
			last_max_eigenvalue = max_eigenvalue
		end
		return max_eigenvalue
	end

	def min_eigenvalue(n=100, error = 10e-12)
		return  self.invert.max_eigenvalue(n, error)
	end

	def expm
		#expm(matrix) = V*diag(exp(diag(D)))/V; V => eigenvectors(right), D => eigenvalues (right). # https://es.mathworks.com/help/matlab/ref/expm.html
		eigenvalues, eigenvectors = NMatrix::LAPACK.geev(self, :right)
		eigenvalues.map!{|val| Math.exp(val)}
		numerator = eigenvectors.dot(NMatrix.diagonal(eigenvalues, dtype: :float32))
		matrix_exp = numerator.div(eigenvectors)
		return matrix_exp
	end

	def cosine_normalization
		normalized_matrix =  NMatrix.zeros(self.shape, dtype: :float32)
		#normalized_matrix =  NMatrix.zeros(self.shape, dtype: :complex64)
		self.each_with_indices do |val, i, j|
			norm = val/CMath.sqrt(self[i, i] * self[j,j])
			#abort("#{norm} has non zero imaginary part" ) if norm.imag != 0
			normalized_matrix[i, j] = norm#.real
		end
		return normalized_matrix
	end
end



