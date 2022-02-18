import skew

# Class representing a suffix array for a string
class suffix_array:
	def __init__(self, string):
		self.string = string #+ "$"
		self.length = len(self.string)
		self.array = skew.skew_rec(self.string)
		self.isa = None
		self.lcp = None

	'''
	# Computes SA, but inefficiently
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


# Class representing an lcp array for a string
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


	# TO DO
	def RMQ(self, i, j):
		return i, 0

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



### TEST CODE ###

s = "mississippi$"
print("str: ", s)

sa1 = suffix_array(s)
print("sa:  ", sa1.array)

sa1.construct_lcp_array()
print("lcp: ", sa1.lcp.array)




