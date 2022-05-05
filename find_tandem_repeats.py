import time
import matplotlib.pyplot as plt
from math import log2
import skew


########################################################
# Class representing a suffix array (SA) for a string
########################################################

class suffix_array:
	def __init__(self, string):
		self.string = string
		self.length = len(self.string)
		self.array = self.construct_array()
		self.isa = {self.array[i]: i for i in range(self.length)}

	'''
	Constructs the suffix array using the Skew algorithm in time O(n)
	'''
	def construct_array(self):
		alpha, ints = skew.map_string_to_ints(self.string)
		return skew.skew_rec(ints, len(alpha))

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
	def __init__(self, x, sa):
		self.sa = sa # Suffix array
		self.string = x
		self.length = len(x)
		self.array = self.construct_lcp()
		self.RMQ_matrix = self.RMQ_preprocess(len(self.array))

	# Check length of shared prefix between suffix i and j
	def compare_lcp(self, x, i, j):
		m = min(len(x) - i, len(x) - j)
		for k in range(m):
			if x[i + k] != x[j + k]:
				return k
		return m

	# Constructs lcp array from suffix array
	def construct_lcp(self):
		lcp = [None] * self.length
		offset = 0
		sa = self.sa
		isa = {sa[i]: i for i in range(self.length)}
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
	# Find child intervals for LCP array
	########################################################

	'''
	Get direct child intervals for an L-interval [i,j)
	I.e. get direct child nodes of a node [i, j)
	'''
	def get_child_intervals(self, i, j):
		if i + 1 == j:
			return [(i, j)]
		else:
			res = []
			(prev_i, L) = lcp.RMQ(i + 1, j) # L is min val in interval, ie. split point
			res.append((i, prev_i)) # first interval (i.e. we split at index prev_i)
			if prev_i + 1 < j:
				(ii, LL) = lcp.RMQ(prev_i + 1, j) # potential second split point
				while LL == L: # If LL is in fact a split point
					res.append((prev_i, ii)) # Add interval to list
					prev_i = ii
					if prev_i + 1 >= j: # If we have no more interval left
						break
					else:
						(ii, LL) = lcp.RMQ(prev_i + 1, j)
			res.append((prev_i, j))
			return res


	def alt_get_child_intervals(self, i, j):
		if i <= j < i+2:
			return
		(prev_i, L) = lcp.RMQ(i + 1, j) # L is min val in interval, ie. split point
		yield (i, prev_i) # first interval (i.e. we split at index prev_i)
		if prev_i + 1 < j:
			(ii, LL) = lcp.RMQ(prev_i + 1, j) # potential second split point
			while LL == L: # If LL is in fact a split point
				yield (prev_i, ii) # Add interval to list
				prev_i = ii
				if prev_i + 1 >= j: # If we have no more interval left
					break
				else:
					(ii, LL) = lcp.RMQ(prev_i + 1, j)
		yield (prev_i, j)



	def get_inner_nodes(self, a, b):
		if b - a <= 2:
			return # No inner nodes for interval of length b-a
		stack = [(a, b)]
		while stack:
			(i, j) = stack.pop()
			yield (i, j)
			stack.extend((ii, jj)
				for (ii, jj) in self.get_child_intervals(i, j)
				if jj > ii + 1)




########################################################
# Finding branching tandem repeats
########################################################

'''
Returns the widest interval [i, j) from a list of intervals
'''
def widest(intervals):
	max_size = 0
	max_interval = (None, None)
	for (i, j) in intervals:
		if j - i > max_size:
			max_size = j - i
			max_interval = (i, j)
	return max_interval



def valid_index(sa, i, j, offset):
	for q in range(i, j):
		if 0 <= sa[q] + offset < len(sa):
			yield q

'''
Finds all branching tandem repeats using the "smaller half trick" to get a
running time of O(n log(n)) for a string of length n
Returns all branching tandem repeats in a list [(string_idx, length)]
(eg. ABCABC has length 6)
'''
def branching_TR_smaller_half(x, sa, lcp):
	isa = [None]*len(sa)
	for i, a in enumerate(sa):
		isa[a] = i
	for (i, j) in lcp.get_inner_nodes(0, len(x)):
		child_nodes = lcp.get_child_intervals(i, j)
		(w_i, w_j) = widest(child_nodes)
		(_, L) = lcp.RMQ(i + 1, j)
		for (ii, jj) in child_nodes:
			if (ii, jj) == (w_i, w_j):
				continue
			for q in valid_index(sa, ii, jj, +L):
				r = isa[sa[q] + L]
				if (i <= r < j) and not (ii <= r < jj):
					yield (sa[q], 2*L)
			for q in valid_index(sa, ii, jj, -L):
				r = isa[sa[q] - L]
				if w_i <= r < w_j:
					yield (sa[r], 2*L)



'''
Finds all tandem repeats in a string, given a list of branching tandem repeats by using left-rotations
Eg. if we have a branching TR at pos i consisting of string w*a, then we can check if we also have a
tandem repeat by left rotation, i.e. by checking if the symbol string[i-1] == a.
'''
def old_find_all_tandem_repeats(x, branching_TRs):
	res = list(branching_TRs)
	stack = res.copy()
	while len(stack) > 0:
		(idx, length) = stack.pop()
		last_symbol = x[idx + length - 1]
		if x[idx - 1] == last_symbol:
			res.append((idx - 1, length))
			stack.append((idx - 1, length))
	return res


def find_all_tandem_repeats(x, branching_TRs):
	for (i, L) in branching_TRs:
		yield (i, L)
		while can_rotate(x, i, L):
			yield (i-1, L)
			i -= 1


def can_rotate(x, i, L):
	return i > 0 and x[i - 1] == x[i + L - 1]

'''
Prints the tandem repeats
'''
def print_TRs(x, TRs):
	print("**********************")
	print("TANDEM REPEATS, (index, length), for string:")
	if len(x) <= 100:
		print("x:", x)
	else:
		print("x of length", str(len(x)), "is too long, so not printed")
	for tr in sorted(TRs):
		idx = tr[0]
		L = tr[1]
		if L <= 40:
			print(str(tr) + ": " + x[idx: idx+L//2], x[idx+L//2: idx+L])
		else:
			print(str(tr) + ": too long, not printed")
	print("**********************")



########################################################
# TEST CODE
########################################################

# Input string
#x = "ACCACCAGTGT$"
#x = "ACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGT$"
#x = "mississippi$"
#x = "banana$"
x = "aaaaaa$"
#x = "aaaa$"

# Suffix array and LCP array
sa = suffix_array(x).array
lcp = lcp_array(x, sa)

#print(lcp.get_inner_nodes(0, len(x)))


# Find all branching TRs
branching_TRs = branching_TR_smaller_half(x, sa, lcp)

# Continuously left-rotate to get all TRs
TRs = find_all_tandem_repeats(x, branching_TRs)

# Print the TRs
print_TRs(x, TRs)

########################################################
# TEST CODE: Running time
########################################################
'''
from numpy.random import choice

alpha = ["A", "C", "G", "T"]
#probs = [0.25] * len(alpha)
probs = [0.1, 0.3, 0.1, 0.5]

def random_string(n):
	return "".join(choice(alpha, n, p=probs))

N = 1000000
lens = range(1, N, 100000) #10

xs = []
for i in lens:
	x = random_string(i)
	#x = "A" * i
	xs.append(x + "$")

times = []
exp_times = []
exp_times2 = []

for x in xs:
	ts = []
	for i in range(5):
		start = time.time() # Start timer
		sa = suffix_array(x).array
		lcp = lcp_array(x, sa)
		branching_TRs = branching_TR_smaller_half(x, sa, lcp)
		tr = find_all_tandem_repeats(x, branching_TRs)
		end = time.time() # Stop timer
		ts.append(end - start)
	# Average running time
	t = sum(ts)/(len(ts))
	times.append(t)
	n = len(x)
	print(n)
	#exp_times.append(t/((n*log2(n))+len(tr)))
	exp_times.append(t/(n**2))
	exp_times2.append(t/((n*log2(n))+len(tr)))


# Time plot
plt.scatter(list(lens), times, color = "blue")
#plt.xticks(range(1,21))
#plt.title("n = " + str(n))
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time (sec)", fontsize = 13)
plt.savefig("time_plot_" + str(N))
plt.show()
plt.clf() # Clear plot

plt.scatter(list(lens), exp_times, color = "blue")
#plt.xticks(range(1,21))
#plt.title("n = " + str(n))
plt.ylim(-0.0001, 0.0001)
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time (sec) / n^2", fontsize = 13)
plt.savefig("time_plot_exp_" + str(N))
plt.show()
plt.clf()

plt.scatter(list(lens), exp_times2, color = "blue")
#plt.xticks(range(1,21))
#plt.title("n = " + str(n))
plt.ylim(-0.0001, 0.0001)
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time (sec) / nlogn + z", fontsize = 13)
plt.savefig("time_plot_exp2_" + str(N))
plt.show()
'''


