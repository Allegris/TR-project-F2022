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








# Test to compare two lists, i.e. do they contain the same elements
def compare_lists(s, t):
    return Counter(s) == Counter(t)

# Test strings
#s = "mississippi0"
#s = "abcabcabc0"
#s = "banana0"
#s = "aaaaa0"
#s = "abababaccccccccccccccccaaaaaaabeaa0"

'''
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
'''