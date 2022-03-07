
import skew
import numpy as np

########################################################
# Class representing a suffix array for a string
########################################################

class suffix_array:
	def __init__(self, string):
		self.string = string #+ "$"
		self.length = len(self.string)
		self.array = skew.skew_rec(self.string)
		self.isa = None
		self.lcp = None

	'''
	# OLD: Computes SA, but inefficiently
	def construct_array(self):
		self.array = [s[1] for s in sorted((self.string[i:], i) for i in range(self.length))]
		return self.array
	'''

	# Constructs the inverse suffix array (for given suffix index, returns suffix rank)
	def construct_isa(self):
		return {self.array[i]: i for i in range(self.length)}

	# Constructs the lcp array for this suffix array
	# (i.e. prefix shared between entry in suffix array and previous entry)
	def construct_lcp_array(self):
		self.lcp = lcp_array(self)
		return self.lcp


########################################################
# Class representing an lcp array for a string
########################################################

class lcp_array:
	def __init__(self, sa):
		self.sa = sa
		self.string = sa.string
		self.length = sa.length
		self.array = self.construct_lcp()

	# Check length of shared prefix between suffix i and j
	def compare_lcp(self, string, i, j):
		m = min(len(string) - i, len(string) - j)
		for k in range(m):
			if string[i + k] != string[j + k]:
				return k
		return m

	# Constructs lcp array from suffix array
	def construct_lcp(self):
		lcp = [None] * self.length
		# Inverse suffix array (contains ranks)
		isa = self.sa.construct_isa()

		offset = 0
		for i in range(self.length):
			offset = max(0, offset - 1)
			# ii <- rank of suffix at index i
			ii = isa[i]
			# If suffix i is the first entry in suffix array, set lcp = 0
			if ii == 0:
				lcp[ii] = 0
				continue
			# Index of the suffix above suffix ii in suffix array
			j = self.sa.array[ii - 1]
			# Prefix that suffix ii and suffix ii-1 share
			offset += self.compare_lcp(self.string, i + offset, j + offset)
			lcp[ii] = offset
		return lcp



	# Get child intervald for an L-lcp interval [i,j)
	def get_child_intervals(self, i, j):
		res = []
		(prev_i, L) = self.RMQ(i, j)
		res.append((i, prev_i)) # first interval
		(ii, LL) = self.RMQ(prev_i + 1, j)
		while LL == L:
			res.append((prev_i, ii))
			prev_i = ii
			(ii, LL) = self.RMQ(prev_i + 1, j)
		res.append((prev_i, j))
		return res

########################################################
# Range Minimum Query (RMQ)
########################################################
	'''
	Computes the range minimum query (RMQ) of interval [i,j), i.e. leftmost occurrence of min value in range [i,j)
	Returns RMQ of the form (index, min_value)
	'''
	def RMQ(self, lcp, i, j):
		M = RMQ_preprocess(len(lcp))
		return i, 0

	def RMQ_preprocess(n):

		# M matrix to fill: n x log(n) (floored),
		# where entry M[i][j] is the RMQ for interval starting at idx i of length 2^j
		log_n = n.bit_length() - 1 # this is log(n) floored, so eg. for 15 => 3 (2^3 = 8)
		M = np.empty((n, log_n))

		# Intervals of length 1
		for i in range(n):
			M[i][0] = lcp[i]

		# Preprocessing/filling out M for intervals of length 2^k (2, 4, 8, 16...)


		# Run through all smaller exponents
		# So if we want to know RMQ of a 2^k, then we run through j = 1,2,3,...,k-1
		for j in range(1, log_n - 1):
			# Run through all intervals of length 2**j (stop if interval exceeds list length)
			last_idx = (n - 2**j) + 1 # last starting idx in list where interval does not exceed list length
			for i in range(last_idx):
				# Every interval of length 2^j can be seen as two intervals of size 2^(j-1)
				# So for an interval of size 2^j, take the minimum of the two 2^(j-1) values
				# These two are the intervals:
				# 1) Starting at pos i with length 2^(j-1)
				# 2) Starting at pos i + 2^(j-1) with length 2^(j-1)
				M[i][j] = min(M[i][j-1], M[i+2**(j-1)][j-1])




########################################################
# TEST CODE
########################################################

s = "mississippi$"
print("str: ", s)

sa1 = suffix_array(s)
print("sa:  ", sa1.array)

sa1.construct_lcp_array()
print("lcp: ", sa1.lcp.array)




