import skew


class suffix_array:
	def __init__(self, string):
		self.string = string #+ "$"
		self.length = len(self.string)
		self.array = skew.skew_rec(self.string)
		self.isa = None
		self.lcp = None

	# Computes SA, but inefficiently
	def construct_array(self):
		self.array = [s[1] for s in sorted((self.string[i:], i) for i in range(self.length))]
		return self.array

	def construct_isa(self):
		return {self.array[i]: i for i in range(self.length)}

	def construct_lcp_array(self):
		self.lcp = lcp_array(self)
		return self.lcp




class lcp_array:
	def __init__(self, sa):
		self.string = sa.string
		self.length = sa.length
		self.array = [None] * self.length

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




sa1 = suffix_array("mississippi")
print(sa1.array)







