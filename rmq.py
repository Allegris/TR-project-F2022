########################################################
# Range Minimum Query (RMQ)
########################################################

'''
Computes the range minimum query (RMQ) of interval [i,j), i.e. (index and value of)
leftmost occurrence of min value in range [i,j)
Returns RMQ of the form (index, min_value)
'''
def RMQ(RMQ_matrix, array, L, R):
	interval = array[L:R] # The interval to do RMQ on
	j = len(interval).bit_length() - 1 # log(interval_len) floor
	# Min of left and right interval of exponents 2^j
	# I.e. interval starting at pos L with length 2^j
	# and interval starting at pos R - 2**j of length 2^j
	# There may be an overlap in the two intervals, but this is OK, result will not change
	right_idx = int(R - (2**j))
	if RMQ_matrix[L][j][1] <= RMQ_matrix[right_idx][j][1]:
		return RMQ_matrix[L][j]
	else:
		return RMQ_matrix[right_idx][j]


'''
Preprocess matrix for RMQ in time O(n*log(n))
Returns RMQ matrix of the form [[(idx, val), (idx, val),...],[(idx, val),...],...]
Where rows are LCP indices and
cols are j = 1, 2, 4, 8, 16,... where we have calculated RMQ for intervals of lengths 2^j
So col 1 is j = 0 (interval length 2^j = 2^0 = 1),
col 2 is j = 1 (interval length 2^1=2), etc.
'''
def RMQ_preprocess(array):
	n = len(array)
	# M matrix to fill: n x log(n),
	# where entry M[i][j] is the RMQ for interval starting at idx i of length 2^j
	log_n = n.bit_length() # this is log(n) ceil, so eg. for 15 => 4 (2^4 = 16)
	#M = [[(sys.maxsize, sys.maxsize)]*(log_n) for _ in range(n)]
	M = [[(None, None)]*(log_n) for _ in range(n)]
	# Intervals of length 1 first:
	for i in range(n):
		M[i][0] = (i, array[i])

	# Preprocessing/filling out M for intervals of length 2^k (2, 4, 8, 16...) next:
	# Run through all exponents, in increasing order: j = 1,2,4,8,..., log n
	for j in range(1, log_n):
		# Run through all intervals of length 2**j (stop if interval exceeds list length)
		last_idx = (n - (2**j)) + 1
		for i in range(last_idx): # note: last_index excluded
			# Every interval of length 2^j can be seen as two intervals of size 2^(j-1)
			# So for an interval of size 2^j, take the minimum of the two 2^(j-1) values
			# These two are the intervals:
			# 1) Starting at pos i with length 2^(j-1)
			# 2) Starting at pos i + 2^(j-1) with length 2^(j-1)
			left_min = M[i][j-1]
			right_min = M[i+(2**(j-1))][j-1]
			M[i][j] = left_min if left_min[1] <= right_min[1] else right_min
	return M
