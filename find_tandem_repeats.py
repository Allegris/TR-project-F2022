import time
import matplotlib.pyplot as plt
from math import log2
import skew
from rmq import RMQ
from rmq import RMQ_preprocess
from line_profiler import LineProfiler
lprofiler = LineProfiler()



########################################################
# Finding branching tandem repeats
########################################################

'''
Finds all branching tandem repeats using the "smaller half trick" to get a
running time of O(n log(n)) for a string of length n.
Yields all branching tandem repeats as tuplets (string_idx, length)
(eg. ABCABC has length 6).
'''
def branching_TR_smaller_half(x, sa, lcp):
	isa = construct_isa(sa)
	M = RMQ_preprocess(lcp)
	for (i, j) in get_inner_nodes(lcp, M, 0, len(x)):
		child_nodes = list(get_child_nodes(lcp, M, i, j))
		(w_i, w_j) = widest(child_nodes)
		(_, L) = RMQ(lcp, M, i + 1, j)
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
Yields the "inner nodes" in the "suffix tree";
i.e. all non-singleton as L-intervals below [a,b).
'''
def get_inner_nodes(lcp, M, a, b):
	if b - a <= 2:
		return # No inner nodes for interval of length < 2
	stack = [(a, b)]
	while stack:
		(i, j) = stack.pop()
		yield (i, j)
		stack.extend((ii, jj)
			for (ii, jj) in get_child_nodes(lcp, M, i, j)
			if jj > ii + 1)

'''
Yields direct "child nodes" of a "node" [i, j) in "suffix tree";
i.e. direct child intervals for an L-interval [i,j).
'''
def get_child_nodes(lcp, M, i, j):
	# Singleton
	if i + 1 == j:
		yield (i, j)
	# Non-singleton
	else:
		(prev_i, L) = RMQ(lcp, M, i + 1, j)
		yield (i, prev_i)
		if prev_i + 1 < j:
			(ii, LL) = RMQ(lcp, M, prev_i + 1, j)
			while LL == L:
				yield (prev_i, ii)
				prev_i = ii
				if prev_i + 1 >= j:
					break
				else:
					(ii, LL) = RMQ(lcp, M, prev_i + 1, j)
		yield (prev_i, j)

'''
HELPER FUNCTION
Returns the widest interval [i, j) from a list of intervals.
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
Yields all tandem repeats in a string, given a list
of branching tandem repeats by using left-rotations.
Eg. if we have a branching TR at pos i consisting of
string w*a, then we can check if we also have a
tandem repeat by left rotation, i.e. by checking if
the symbol string[i-1] == a.
'''
def find_all_tandem_repeats(x, branching_TRs):
	for (i, L) in branching_TRs:
		yield (i, L)
		while can_rotate(x, i, L):
			yield (i-1, L)
			i -= 1

'''
HELPER FUNCTION
Checks if we can left-rotate x[i:L) and get a tandem repeat
'''
def can_rotate(x, i, L):
	return i > 0 and x[i - 1] == x[i + L - 1]


########################################################
# Arrays: SA, ISA, LCP
########################################################

'''
Constructs suffix array using the Skew algorithm in time O(n)
'''
def construct_sa(x):
	alpha, indices = skew.map_string_to_ints(x)
	return skew.skew_rec(indices, len(alpha))

'''
Slow (but simple) function for constructing suffix array.
Only used to test against Skew algo results.
'''
def construct_sa_slow(x):
	return [s[1] for s in sorted((x[i:], i) for i in range(len(x)))]

'''
Constructs the inverse suffix array from a suffix array.
'''
def construct_isa(sa):
	isa = [None]*len(sa)
	for i, a in enumerate(sa):
		isa[a] = i
	return isa

'''
Constructs the LCP array from a suffix array.
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
HELPER FUNCTION
Returns the length of longest shared prefix between suffix i and j.
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
# TEST CODE: Correctness + profiling
########################################################
'''
def main_func():
	# Input string
	x = "AAAAAAAAAAAAAAAAAAAAAAAAAA$"
	#x = "ACCACCAGTGT$"
	#x = "ACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGTACCACCAGTGT$"
	#x = "mississippi$"
	#x = "banana$"
	#x = "aaaaaa$"
	#x = "aaaa$"

	# Suffix array and LCP array
	sa = construct_sa(x)
	#lcp = lcp_array(x, sa)
	lcp = construct_lcp(x, sa)
	global M;
	M = RMQ_preprocess(lcp)

	# Find all branching TRs
	branching_TRs = branching_TR_smaller_half(x, sa, lcp)

	# Continuously left-rotate to get all TRs
	TRs = find_all_tandem_repeats(x, branching_TRs)

	# Print the TRs
	print_TRs(x, TRs)

#main_func()

# Profiling
def prof():
	lprofiler = LineProfiler()
	lprofiler.add_function(find_all_tandem_repeats)
	lp_wrapper = lprofiler(main_func)
	lp_wrapper()
	lprofiler.print_stats()

#prof()
'''
########################################################
# TEST CODE: Running time
########################################################

from numpy.random import choice

alpha = ["A", "C", "G", "T"]
probs = [0.25] * len(alpha)
#probs = [0.1, 0.3, 0.1, 0.5]

def random_string(n):
	return "".join(choice(alpha, n, p=probs))


N = 10000
lens = range(1, N, 500)

xs = []
for n in lens:
	x = random_string(n) # random data
	#x = "A" * n # worst case data
	xs.append(x + "$")


## Split into two algos ##
# Algo 1: Find BRANCHING TRs
ns = []
times = []
exp_times = []
# Algo 2: Find aLL TRs (from branching ones)
ns2 = []
zs = []
times2 = []
exp_times2 = []

for x in xs:
	ts = []
	ts2 = []
	tr = []
	for i in range(5):
		sa = construct_sa(x)
		lcp = construct_lcp(x, sa)

		# Algo 1
		start = time.perf_counter_ns() # Start timer
		branching_TRs = branching_TR_smaller_half(x, sa, lcp)
		end = time.perf_counter_ns() # Stop timer
		# Append time and convert from nanosec to sec
		ts.append(end*10**(-9) - start*10**(-9))

		#Algo 2
		start2 = time.perf_counter_ns()
		tr = find_all_tandem_repeats(x, branching_TRs)
		end2 = time.perf_counter_ns()
		# Append time and convert from nanosec to sec
		ts2.append(end2*10**(-9) - start2*10**(-9))

	# Algo 1
	# Average running time
	n = len(x)
	ns.append(n)
	t = sum(ts)/len(ts)
	times.append(t)
	exp_times.append(t/((n*log2(n))))

	# Algo 2
	# If expected time should be O(n^2)
	'''
	z = len(list(tr))
	zs.append(z)
	t2 =  sum(ts2)/len(ts2)
	times2.append(t2)
	exp_times2.append(t2/(n**2))
	'''

	# If expected time should be O(z) - make sure we do not divide by z if z = 0
	z = len(list(tr))
	if z != 0:
		zs.append(z)
		ns2.append(n)
		t2 =  sum(ts2)/len(ts2)
		times2.append(t2)
		exp_times2.append(t2/z)

######## PLOTS ########

### Algo 1 ###
# Time plot
plt.scatter(ns, times, color = "blue")
plt.ylim(0, 6*(10**(-6)))
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time (sec)", fontsize = 13)
plt.savefig("btr_random_time_plot_" + str(N))
plt.show()
plt.clf() # Clear plot

# Exp time plot
plt.scatter(ns, exp_times, color = "blue")
plt.ylim(0, 5*(10**(-10)))
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time / (n log n)", fontsize = 13)
plt.savefig("btr_random_time_plot_exp_" + str(N))
plt.show()
plt.clf() # Clear plot

### Algo 2 ###
# Z plot
plt.scatter(ns2, zs, color = "blue")     #ns
#plt.ylim(0, 4*(10**(-4)))
plt.xlabel("n", fontsize = 13)
plt.ylabel("z", fontsize = 13)
plt.savefig("z_alltr_random_time_plot_" + str(N))
plt.show()
plt.clf() # Clear plot

# Time plot
plt.scatter(ns2, times2, color = "blue")   #ns
plt.ylim(0, 2*(10**(-4)))
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time (sec)", fontsize = 13)
plt.savefig("alltr_random_time_plot_" + str(N))
plt.show()
plt.clf() # Clear plot

# Exp time plot
'''
plt.scatter(ns, exp_times2, color = "red")
plt.ylim(0, 4*(10**(-11)))
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time / n^2", fontsize = 13)
plt.savefig("alltr_wc_time_plot_exp_" + str(N))
plt.show()
plt.clf() # Clear plot
'''

# Exp time plot against Z
plt.scatter(ns2, exp_times2, color = "blue")
#plt.ylim(0, 4*(10**(-10))) #wc
plt.ylim(0, 2*(10**(-7))) #random
plt.xlabel("n", fontsize = 13)
plt.ylabel("Time / z", fontsize = 13)
plt.savefig("ZZZalltr_random_time_plot_exp_" + str(N))
plt.show()
plt.clf() # Clear plot
