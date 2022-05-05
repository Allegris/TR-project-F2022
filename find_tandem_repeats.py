import time
import matplotlib.pyplot as plt
from math import log2
import skew
import rmq


########################################################
# Finding branching tandem repeats
########################################################

'''
Finds all branching tandem repeats using the "smaller half trick" to get a
running time of O(n log(n)) for a string of length n
Returns all branching tandem repeats in a list [(string_idx, length)]
(eg. ABCABC has length 6)
'''
def branching_TR_smaller_half(x, sa, lcp):
	isa = construct_isa(sa)
	for (i, j) in get_inner_nodes(lcp, 0, len(x)):
		child_nodes = list(get_child_intervals(lcp, i, j))
		(w_i, w_j) = widest(child_nodes)
		(_, L) = rmq.RMQ(lcp, i + 1, j)
		for (ii, jj) in child_nodes:
			if (ii, jj) == (w_i, w_j):
				continue
			for q in valid_isa_index(sa, ii, jj, +L):
				r = isa[sa[q] + L]
				if (i <= r < j) and not (ii <= r < jj):
					yield (sa[q], 2*L)
			for q in valid_isa_index(sa, ii, jj, -L):
				r = isa[sa[q] - L]
				if w_i <= r < w_j:
					yield (sa[r], 2*L)

'''
Yields the "inner nodes" in the "suffix tree"
as L-intervals; i.e. only the non-singleton L-intervals
'''
def get_inner_nodes(lcp, a, b):
	if b - a <= 2:
		return # No inner nodes for interval of length < 2
	stack = [(a, b)]
	while stack:
		(i, j) = stack.pop()
		yield (i, j)
		stack.extend((ii, jj)
			for (ii, jj) in get_child_intervals(lcp, i, j)
			if jj > ii + 1)

'''
Get direct child intervals for an L-interval [i,j)
I.e. get direct child nodes of a node [i, j)
'''
def get_child_intervals(lcp, i, j):
	# Singleton
	if i + 1 == j:
		yield (i, j)
	# Non-singleton
	else:
		(prev_i, L) = rmq.RMQ(lcp, i + 1, j)
		yield (i, prev_i)
		if prev_i + 1 < j:
			(ii, LL) = rmq.RMQ(lcp, prev_i + 1, j)
			while LL == L:
				yield (prev_i, ii)
				prev_i = ii
				if prev_i + 1 >= j:
					break
				else:
					(ii, LL) = rmq.RMQ(lcp, prev_i + 1, j)
		yield (prev_i, j)

'''
HELPER FUNCTION
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

'''
HELPER FUNCTION
Yields the q's in the interval [i, j) for
which sa[q]+offset is within the interval [0, len(sa)).
Used for checking if a given index is valid for isa.
'''
def valid_isa_index(sa, i, j, offset):
	for q in range(i, j):
		if 0 <= sa[q] + offset < len(sa):
			yield q

########################################################
# Finding all tandem repeats from the branching ones
########################################################

'''
Finds all tandem repeats in a string, given a list of branching tandem repeats by using left-rotations
Eg. if we have a branching TR at pos i consisting of string w*a, then we can check if we also have a
tandem repeat by left rotation, i.e. by checking if the symbol string[i-1] == a.
'''
def find_all_tandem_repeats(x, branching_TRs):
	for (i, L) in branching_TRs:
		yield (i, L)
		while can_rotate(x, i, L):
			yield (i-1, L)
			i -= 1

'''
Checks if we can left-rotate x[i:n) and get a tandem repeat
'''
def can_rotate(x, i, L):
	return i > 0 and x[i - 1] == x[i + L - 1]


########################################################
# Arrays: SA, ISA, LCP
########################################################

'''
Construct suffix array using the Skew algorithm in time O(n)
'''
def construct_sa(x):
	alpha, ints = skew.map_string_to_ints(x)
	return skew.skew_rec(ints, len(alpha))

'''
Slow, but simple, function for constructing suffix array.
Only used to test against Skew algo results.
'''
def construct_sa_slow(x):
	return [s[1] for s in sorted((x[i:], i) for i in range(len(x)))]

'''
Constructs the inverse suffix array from the suffix array
'''
def construct_isa(sa):
	isa = [None]*len(sa)
	for i, a in enumerate(sa):
		isa[a] = i
	return isa

'''
Construct LCP array from suffix array
'''
def construct_lcp(x, sa):
	isa = construct_isa(sa)
	lcp = [None] * len(x)
	offset = 0
	for i in range(len(x)):
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
		offset += shared_prefix_len(x, i + offset, j + offset)
		lcp[ii] = offset
	return lcp

'''
Returns the length of longest shared prefix between suffix i and j
'''
def shared_prefix_len(x, i, j):
	m = min(len(x) - i, len(x) - j)
	for k in range(m):
		if x[i + k] != x[j + k]:
			return k
	return m



########################################################
# Printing TRs
########################################################

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
# TEST CODE: Correctness
########################################################

# Input string
#x = "ACCACCAGTGT$"
#x = "ACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGT$"
x = "mississippi$"
#x = "banana$"
#x = "aaaaaa$"
#x = "aaaa$"

# Suffix array and LCP array
sa = construct_sa(x)
#lcp = lcp_array(x, sa)
lcp = construct_lcp(x, sa)

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


