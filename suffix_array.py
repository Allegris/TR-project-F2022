
import skew
import sys
from pandas import *

########################################################
# Class representing a suffix array for a string
########################################################

class suffix_array:
	def __init__(self, string):
		self.string = string #+ "$"
		self.length = len(self.string)
		self.array = skew.skew_rec(self.string)
		self.isa = {self.array[i]: i for i in range(self.length)}
		self.lcp = None


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
		self.sa = sa # Suffix array
		self.isa = sa.isa # Inverse suffix array (contains ranks)
		self.string = sa.string
		self.length = sa.length
		self.array = self.construct_lcp()
		self.RMQ_matrix = self.RMQ_preprocess(len(self.array))

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

		offset = 0
		for i in range(self.length):
			offset = max(0, offset - 1)
			# ii <- rank of suffix at index i
			ii = self.isa[i]
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


	########################################################
	# Range Minimum Query (RMQ) for LCP array
	########################################################

	'''
	Computes the range minimum query (RMQ) of interval [i,j), i.e. (index and value of)
	leftmost occurrence of min value in range [i,j)
	Returns RMQ of the form (index, min_value)
	'''
	def RMQ(self, L, R):
		interval = self.array[L:R] # The interval to do RMQ on
		j = len(interval).bit_length() - 1 # log(interval_len) floor
		# Min of left and right interval of exponents 2^j
		# I.e. interval starting at pos L with length 2^j
		# and interval starting at pos R - 2**j of length 2^j
		# There may be an overlap in the two intervals, but this is OK, result will not change
		right_idx = int(R - (2**j))
		if self.RMQ_matrix[L][j][1] <= self.RMQ_matrix[right_idx][j][1]:
			return self.RMQ_matrix[L][j]
		else:
			return self.RMQ_matrix[right_idx][j]


	'''
	Preprocess matrix for RMQ in time O(n*log(n))
	Returns RMQ matrix of the form [[(idx, val), (idx, val),...],[(idx, val),...],...]
	Where rows are LCP indices and
	Cols are j=1,2,4,8,16,... where we have calculated RMQ for intervals of lengths 2^j
	So col 1 is j=0 (interval length 2^j = 2^0 = 1), col 2 is j=1 (interval length 2^1=2), etc.
	'''
	def RMQ_preprocess(self, n):
		# M matrix to fill: n x log(n),
		# where entry M[i][j] is the RMQ for interval starting at idx i of length 2^j
		log_n = n.bit_length() # this is log(n) ceil, so eg. for 15 => 4 (2^4 = 16)
		#M = [[(sys.maxsize, sys.maxsize)]*(log_n) for _ in range(n)]
		M = [[(None, None)]*(log_n) for _ in range(n)]
		# Intervals of length 1 first:
		for i in range(n):
			M[i][0] = (i, self.array[i])

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


	########################################################
	# Finding branching tandem repeats
	########################################################

	'''
	Finds all tandem repeats in a string, given a list on branching tandem repeats by using left-rotations
	Eg. if we have a branching TR at pos i+1 consisting of string w*a, then we can check if we also have a
	tandem repeat by left rotation, i.e. by checking if the symbol string[i] == a.
	'''
	def find_all_tandem_repeats(string, branching_TRs):
		res = []
		for (idx, length) in branching_TRs:
			i = idx - 1
			while i >= 0:
				btr_symbol = string[idx + length - 1]
				if string[i] == btr_symbol:
					res.append((i, length))
				else:
					break



	'''
	Finds all branching tandem repeats using the "smaller half trick" to get a
	running time of O(n log(n)) for a string of length n
	Returns all branching tandem repeats in a list [(string_idx, L)]
	Where L is the length of the repeated substring, so the tandem repeat "AA" has length 2*L
	'''
	def branching_TR_smaller_half(self, i, j):
		res = []
		# Suffix array and inverse suffix array
		sa = self.sa.array
		isa = self.isa
		# All inner nodes in suffix tree
		inner_nodes = [(i, j) for (i, j) in self.child_intervals_rec(i, j) if j - i > 1]
		# For each inner node v, check if each leaf below is a tandem repeat:
		# A leaf is a tandem repeat if isa[sa[leaf] + depth(v)] is below v,
		# but in a different subtree than leaf.
		for (node_i, node_j) in inner_nodes:
			child_intervals = self.child_intervals_rec(node_i, node_j)
			# The widest subtree/interval
			(w_i, w_j) = self.widest(child_intervals)
			# L is the shared node depth for all children in this subtree (i.e. child interval):
			(_, L) = self.RMQ(node_i+1, node_j)
			# Run through each subtree, EXCEPT THE WIDEST SUBTREE
			for (ii, jj) in child_intervals:
				if (ii, jj) == (w_i, w_j):
					continue
				# Run through each leaf in subtree
				for q in range(ii, jj):
					# *** Check to the right ***
					# "Leaf sa[q]+L" that may potentially form a branching TR with "leaf sa[q]"
					r = isa[sa[q] + L]
					# If "leaf sa[q]+L" is in a different subtree
					# than "leaf sa[q]", then we have a branching TR:
					if r in range(node_i, ii) or r in range(jj, node_j):
						res.append((sa[q], L))
					# *** Check to the left (widest subtree) ***
					r = isa[sa[q] - L]
					# If "leaf sa[q]-L" is in the widest subtree, then we would have missed it,
					# since we don't run through the widest subtree, so we add it here
					if r in range(w_i, w_j):
						res.append((sa[r], L))
		return list(set(res)) # remove duplicates

	'''
	Returns the widest interval (i, j) from a list of intervals
	'''
	def widest(self, intervals):
		max_size = 0
		max_interval = (None, None)
		for (i, j) in intervals:
			if j - i > max_size:
				max_size = j - i
				max_interval = (i, j)
		return max_interval


	'''
	Finds all branching tandem repeats in running time O(n^2) for a string of length n
	Returns all branching tandem repeats in a list [(string_idx, L)]
	Where L is the length of the repeated substring, so the tandem repeat "AA" has length 2*L
	'''
	def branching_TR(self, i, j):
		res = []
		# Suffix array and inverse suffix array
		sa = self.sa.array
		isa = self.isa
		# All inner nodes in suffix tree
		inner_nodes = [(i, j) for (i, j) in self.child_intervals_rec(i, j) if j - i > 1]
		# For each inner node v, check if each leaf below is a tandem repeat:
		# A leaf is a tandem repeat if isa[sa[leaf] + depth(v)] is below v,
		# but in a different subtree than leaf.
		for (node_i, node_j) in inner_nodes:
			child_intervals = self.child_intervals_rec(node_i, node_j)
			# L is the shared node depth for all children in this subtree (i.e. child interval):
			(_, L) = self.RMQ(node_i+1, node_j)
			# Run through each subtree
			for (ii, jj) in child_intervals:
				# Run through each leaf in subtree
				for q in range(ii, jj):
					# "Leaf sa[q]+L" that may potentially form a branching TR with "leaf sa[q]"
					r = isa[sa[q] + L]
					# If "Leaf sa[q]+L" is in a different subtree
					# than "leaf sa[q]", then we have a branching TR:
					if r in range(node_i, ii) or r in range(jj, node_j):
						res.append((sa[q], L))
		return list(set(res)) # remove duplicates


	'''
	Recursive funtion that finds ALL the child intervals
	(so also the child intervals of the child intervals, all the way down to the leaves, (i, i+1))
	'''
	def child_intervals_rec(self, i, j):
		res = []
		child_intervals = self.get_child_intervals(i, j)
		res += child_intervals
		if len(child_intervals) > 1:
			for (ii, jj) in child_intervals:
				if jj-ii > 1: # if non-empty, non-singleton interval
					res += self.child_intervals_rec(ii, jj)
		return list(set(res)) # remove duplicates


	'''
	Get (direct) child intervals for an L-lcp interval [i,j)
	'''
	def get_child_intervals(self, i, j):
		# Empty interval, eg. (5, 5)
		if i == j:
			return []
		# Singleton interval, eg. (4, 5)
		elif i + 1 == j:
			return [(i, j)]
		# Non-empty, non-singleton interval of length 2, eg. (3, 5)
		# => We want both (3, 4) and (4, 5) as child intervals
		elif i + 2 == j:
			return [(i, i + 1), (i + 1, j)]
		# Non-empty, non-singleton interval of length > 2, eg. (2, 5)
		else:
			res = []
			(prev_i, L) = self.RMQ(i + 1, j) # L is min val in interval, ie. split point
			res.append((i, prev_i)) # first interval (i.e. we split at index prev_i)
			(ii, LL) = self.RMQ(prev_i + 1, j) # potential second split point
			while LL == L: # If LL is in fact a split point
				res.append((prev_i, ii)) # Add interval to list
				prev_i = ii
				if prev_i + 1 >= j: # If we have no more interval left
					break
				else:
					(ii, LL) = self.RMQ(prev_i + 1, j)
			res.append((prev_i, j))
			return res



########################################################
# TEST CODE
########################################################

s = "mississippi$"
print("str: ", s)

sa1 = suffix_array(s)
print("sa:  ", sa1.array)

lcp = sa1.construct_lcp_array()
print("lcp: ", lcp.array)

#print(lcp.RMQ(0, 12))
#print(DataFrame(sa1.lcp.RMQ_matrix))

#child_intervals = lcp.get_child_intervals(0, len(lcp.array))
#print("Child intervals: ", child_intervals)

#print("Child intervals: ", lcp.get_child_intervals(0, 9))
#print("All child intervals (recursive): ", lcp.child_intervals_rec(0, 12))
print("Branding tandem repeats: ", lcp.branching_TR(0, 12))
print("Branding tandem repeats SMALLER HALF: ", lcp.branching_TR_smaller_half(0, 12))




