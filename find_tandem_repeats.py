import skew
from collections import Counter
from pandas import *

########################################################
# Class representing a suffix array (SA) for a string
########################################################

class suffix_array:
	def __init__(self, string):
		self.string = string
		self.length = len(self.string)
		self.array = self.construct_array()
		self.isa = {self.array[i]: i for i in range(self.length)}
		self.lcp = None

	'''
	Constructs the suffix array using the Skew algorithm in time O(n)
	'''
	def construct_array(self):
		alpha, ints = skew.map_string_to_ints(self.string)
		return skew.skew_rec(ints, len(alpha))

	'''
	Constructs the LCP array for this suffix array
	(i.e. prefix shared between entry in suffix array and previous entry)
	'''
	def construct_lcp_array(self):
		self.lcp = lcp_array(self)
		return self.lcp

	'''
	Time consuming (but simple) function for finding SA
	Used to check if we get correct result from Skew Algo
	'''
	def slow_SA(self):
		return [s[1] for s in sorted((self.string[i:],i) for i in range(self.length))]


########################################################
# Class representing an LCP array for a string
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
		isa = self.isa
		sa = self.sa.array
		for i in range(self.length):
			offset = max(0, offset - 1)
			# ii <- rank of suffix at index i
			ii = isa[i]
			# If suffix i is the first entry in suffix array, set lcp = 0
			if ii == 0:
				lcp[ii] = 0
				continue
			# Index of the suffix above suffix ii in suffix array
			j = sa[ii - 1]
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
	cols are j = 1, 2, 4, 8, 16,... where we have calculated RMQ for intervals of lengths 2^j
	So col 1 is j = 0 (interval length 2^j = 2^0 = 1),
	col 2 is j = 1 (interval length 2^1=2), etc.
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
	def find_all_tandem_repeats(self, string, branching_TRs):
		res = branching_TRs.copy()
		stack = branching_TRs.copy()
		while len(stack) > 0:
			(idx, length) = stack.pop()
			last_symbol = string[idx + length - 1]
			if string[idx - 1] == last_symbol:
				res.append((idx - 1, length))
				stack.append((idx - 1, length))
		return list(set(res)) # Remove duplicates


	'''
	Finds all branching tandem repeats using the "smaller half trick" to get a
	running time of O(n log(n)) for a string of length n
	Returns all branching tandem repeats in a list [(string_idx, L)]
	Where L is the length of the repeated substring, so the tandem repeat "AA" has length 2*L
	'''
	def old_branching_TR_smaller_half(self, string):
		res = []
		# Suffix array and inverse suffix array
		sa = self.sa.array
		isa = self.isa
		# All inner nodes in suffix tree
		inner_nodes = [(a, b) for (a, b) in self.child_intervals_rec(0, len(string)) if b - a > 1]
		# For each inner node v, check if each leaf below is a tandem repeat:
		# A leaf is a tandem repeat if isa[sa[leaf] + depth(v)] is below v,
		# but in a different subtree than leaf
		for (node_i, node_j) in inner_nodes:
			child_intervals = self.child_intervals_rec(node_i, node_j)
			# The widest subtree/interval
			(w_i, w_j) = self.widest(child_intervals)
			# L is the shared node depth for all children in this subtree (i.e. child interval):
			(_, L) = self.RMQ(node_i + 1, node_j)
			# Run through each subtree, EXCEPT THE WIDEST SUBTREE
			for (ii, jj) in child_intervals:
				if (ii, jj) == (w_i, w_j): # ignore widest subtree
					continue
				# Run through each leaf in subtree
				for q in range(ii, jj):
					# *** Check to the right ***
					# "Leaf sa[q]+L" that may potentially form a branching TR with "leaf sa[q]"
					if sa[q] + L in range(0, len(string)):
						r = isa[sa[q] + L]
						# If "leaf sa[q]+L" is in a different subtree
						# than "leaf sa[q]", then we have a branching TR:
						if r in range(node_i, ii) or r in range(jj, node_j):
							res.append((sa[q], L))
					# *** Check to the left (widest subtree) ***
					if sa[q] - L in range(0, len(string)):
						r = isa[sa[q] - L]
						# If "leaf sa[q]-L" is in the widest subtree, then we would have missed it,
						# since we don't run through the widest subtree, so we add it here
						if r in range(w_i, w_j):
							res.append((sa[r], L))
		return list(set(res)) # remove duplicates

	'''
	Finds all branching tandem repeats using the "smaller half trick" to get a
	running time of O(n log(n)) for a string of length n
	Returns all branching tandem repeats in a list [(string_idx, L)]
	Where L is the length of the repeated substring, so the tandem repeat "AA" has length 2*L
	'''
	def branching_TR_smaller_half(self, x):
		res = []
		# Suffix array and inverse suffix array
		sa = self.sa.array
		isa = self.isa
		# All inner nodes in suffix tree
		L_intervals = [(i, j) for (i, j) in self.child_intervals_rec(0, len(x)) if j - i > 1]
		# For each inner node v, check if each leaf below is a tandem repeat:
		# A leaf is a tandem repeat if isa[sa[leaf] + depth(v)] is below v,
		# but in a different subtree than leaf
		for (i, j) in L_intervals:
			child_intervals = self.child_intervals_rec(i, j)
			# The widest subtree/interval
			(w_i, w_j) = self.widest(child_intervals)
			# L is the shared node depth for all children in this subtree (i.e. child interval):
			(_, L) = self.RMQ(i + 1, j)
			# Run through each subtree, EXCEPT THE WIDEST SUBTREE
			for (ii, jj) in child_intervals:
				if (ii, jj) == (w_i, w_j): # ignore widest subtree
					continue
				# Run through each leaf in subtree
				for q in range(ii, jj):
					# *** Check to the right ***
					# "Leaf sa[q]+L" that may potentially form a branching TR with "leaf sa[q]"
					if sa[q] + L in range(0, len(x)):
						r = isa[sa[q] + L]
						# If "leaf sa[q]+L" is in a different subtree
						# than "leaf sa[q]", then we have a branching TR:
						if r in range(i, ii) or r in range(jj, j):
							res.append((sa[q], L))
					# *** Check to the left (widest subtree) ***
					if sa[q] - L in range(0, len(x)):
						r = isa[sa[q] - L]
						# If "leaf sa[q]-L" is in the widest subtree, then we would have missed it,
						# since we don't run through the widest subtree, so we add it here
						if r in range(w_i, w_j):
							res.append((sa[r], L))
		return list(set(res)) # remove duplicates


	'''
	Returns the widest interval [i, j) from a list of intervals
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
	def branching_TR(self, string):
		res = []
		# Suffix array and inverse suffix array
		sa = self.sa.array
		isa = self.isa
		# All inner nodes in suffix tree
		inner_nodes = [(i, j) for (i, j) in self.child_intervals_rec(0, len(string)) if j - i > 1]
		# For each inner node v, check if each leaf below is a tandem repeat:
		# A leaf is a tandem repeat if isa[sa[leaf] + depth(v)] is below v,
		# but in a different subtree than leaf.
		for (node_i, node_j) in inner_nodes:
			child_intervals = self.child_intervals_rec(node_i, node_j)
			# L is the shared node depth for all children in this subtree (i.e. child interval):
			(_, L) = self.RMQ(node_i + 1, node_j)
			# Run through each subtree
			for (ii, jj) in child_intervals:
				# Run through each leaf in subtree
				for q in range(ii, jj):
					if sa[q] + L in range(0, len(string)):
						# "Leaf sa[q]+L" that may potentially form a branching TR with "leaf sa[q]"
						r = isa[sa[q] + L]
						# If "Leaf sa[q]+L" is in a different subtree
						# than "leaf sa[q]", then we have a branching TR:
						if r in range(node_i, ii) or r in range(jj, node_j):
							res.append((sa[q], L))
		return list(set(res)) # remove duplicates


	'''
	Recursive funtion that finds ALL the child intervals
	(so also the child intervals of the child intervals, all the way down to the leaves, [i, i+1))
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
			if prev_i + 1 < j:
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

	'''
	Prints the tandem repeats
	'''
	def print_TRs(self, string, TRs):
		print("**********************")
		print("TANDEM REPEATS:")
		for tr in TRs:
			first_end =  tr[0] + tr[1]
			if tr[1] <= 20:
				print("Idx", tr[0], ": ", string[tr[0]: first_end], ",", string[first_end:first_end + tr[1]])
			else:
				print("Idx", tr[0], ": repeat of length", tr[1])
		print("**********************")



########################################################
# TEST CODE
########################################################

# Test to compare two lists, i.e. do they contain the same elements
def compare_lists(s, t):
    return Counter(s) == Counter(t)

# Test strings
#s = "mississippi0"
#s = "abcabcabc0"
#s = "banana0"
#s = "aaaaa0"
#s = "abababaccccccccccccccccaaaaaaabeaa0"
s = "ACCACCAGTGT$"

print("Input string of len", len(s), ":", s)

sa1 = suffix_array(s)
print("sa:  ", sa1.array)
print("Is Skew SA correct: ", compare_lists(sa1.array, sa1.slow_SA()))

lcp = sa1.construct_lcp_array()
print("LCP: ", lcp.array)

# The two methods for finding branching tandem repeats
branching_TRs_1 = lcp.branching_TR(s)
branching_TRs = lcp.branching_TR_smaller_half(s)
print("Are branching TRs found by the two methods the same: ", compare_lists(branching_TRs_1, branching_TRs))

# Finding all tandem repeats
TRs = lcp.find_all_tandem_repeats(s, branching_TRs)
print("Tandem repeats (idx, len): ", TRs)
#lcp.print_TRs(s, TRs)



# OLD TESTS:
#print(lcp.RMQ(0, 12))
#print(DataFrame(sa1.lcp.RMQ_matrix))
#print("All child intervals (recursive): ", lcp.child_intervals_rec(0, 12))
#print("Branching tandem repeats: ", branching_TRs_1)
#print("SMALLER HALF TRICK Branching tandem repeats: ", branching_TRs)