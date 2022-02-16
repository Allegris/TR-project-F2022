
class suffix_array:
	def __init__(self, string):
		self.string = string + "$"
		self.length = len(self.string)
		self.array = [None] * self.length # self.construct_array_skew()
		self.isa = None
		self.lcp = None

	# Computes SA, but inefficiently
	def construct_array(self):
		self.array = [s[1] for s in sorted((self.string[i:], i) for i in range(self.length))]
		return self.array


	##### SKEW ALGORITHM FOR COMPUTING SUFFIX ARRAY #####
	def construct_array_skew(self):
		padded_string = self.string + "$$"
		sa_12 = self.construct_sa_12(padded_string)
		self.array = sa_12

	def construct_sa_12(self, string):

		print(string)

		# Create triples and sort them
		# [('$$$', 11), ('i$$', 10), ('ipp', 7), ...]
		suffixes_12 = [s for s in sorted((string[i:i+3], i) for i in range(self.length) if i%3 != 0)]

		# [11, 10, 7, ...]
		sa_12 = [s[1] for s in suffixes_12]
		# ['$$$', 'i$$', 'ipp', ...]
		sa_12_triples =  [s[0] for s in suffixes_12]
		print(sa_12_triples)

		# If no duplicates in triple list, we are done
		if len(sa_12_triples) == len(set(sa_12_triples)):
			print("DONE")
			return sa_12
		else:
			# Assign new lex names to triples
			# [(11, '$$$', 0), (10, 'i$$', 1), (7, 'ipp', 2), ...]
			index_to_lexname_dict = dict([(suffixes_12[i][1], (suffixes_12[i][0], i)) for i in range(len(suffixes_12))])

			# Dict with key: index of suffix into string, and val: (triple, lex name)
			# {11: ('$$$', 0), 10: ('i$$', 1), 7: ('ipp', 2), ...}
			lexname = 0
			index_to_lexname_dict = dict()
			index_to_lexname_dict[suffixes_12[0][1]] = (suffixes_12[0][0], lexname)
			for i in range(1, len(suffixes_12)):
				if suffixes_12[i][0] != suffixes_12[i-1][0]:
					lexname += 1
				index_to_lexname_dict[suffixes_12[i][1]] = (suffixes_12[i][0], lexname)


			# All the mod 1 suffixes (non-sorted)
			suffixes_1 = [i for i in range(self.length) if i%3 == 1]
			# All the mod 2 suffixes (non-sorted)
			suffixes_2 = [i for i in range(self.length) if i%3 == 2]

			# Create u string to recurse on
			# 5540#3321  (lex names for mod 2 suffixes, '#', and lex names for mod 1 suffixes)
			u = ""
			for i in suffixes_2:
				u += str(index_to_lexname_dict[i][1])
			u += "#"
			for i in suffixes_1:
				u += str(index_to_lexname_dict[i][1])

			self.construct_sa_12(u)


	######################################################


	def construct_isa(self):
		self.isa = [None] * self.length
		for i in range(0, self.length):
			self.isa[self.array[i]] = i
		return self.isa

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
lcp1 = sa1.construct_lcp_array()
#print(lcp1.string)

sa1.construct_array()
#print(sa1.array)

sa1.construct_array_skew()
#print(sa1.array)







