
import collections

SENTINEL_terminal = "$"
SENTINEL_central = "#"


########################################################
# SKEW ALGORITHM - COMPUTE SUFFIX ARRAY IN TIME O(N)
########################################################

'''
Returns the suffix array of a string
'''
def skew_rec(string):
	### STEP 1: COMPUTE SA_12 ###
	sa_12 = compute_sa_12(string)
	alphabet = get_alphabet(string, sa_12)
	# If all unique in sa_12, don't recurse
	# If the two have equal length, then we have a duplicate in sa_12, because alphabet contains sentinel
	if len(sa_12) >= len(alphabet):
		# Compute u string and recurse to find its suffix array
		u = construct_u(string, alphabet)
		sa_u = skew_rec(u)
		# The index of the central sentinel
		m = len(u) // 2
		# Map the indices in sa_u back to indices in the original string
		sa_12 = [map_u_to_string_index(i, m) for i in sa_u if i != m] # i == m is central sentinel (#)

	### STEP 2: COMPUTE SA_3 from SA_12 ###
	sa_3 = []
	# Special case: if last index in string is 0 mod 3, then this suffix should be first in sa_3, since it is the shortest suffix
	# It will be sorted based on the first (and in this case only) character afterwards
	if len(string) % 3 == 1:
		sa_3.append(len(string) - 1)
	# All strings mod 3 = 0 are just strings mod 3 = 1, with one char in front, so:
	# Use the order found in sa_12, and then stable sort on the first character afterwards
	sa_3 += [i - 1 for i in sa_12 if i % 3 == 1]
	# Stable sort on first character
	sa_3 = bucket_sort_first_char(string, sa_3) # OLD: Works, but O(n*lg n) # sa_3 = sorted(sa_3, key = lambda i : safe_get_char(string, i)[:1])

	### STEP 3: MERGE SA_12 and SA_3 into the final suffix array ###
	return merge(string, sa_12, sa_3)


'''
COMPUTES SA_12
'''
# Computes the suffix array for the suffixes starting at indices != 0 mod 3
def compute_sa_12(string):
	triplets_12 = []
	# Create triplets for all indices != 0 mod 3 in the string
	for i in range(len(string)):
		# Only use mod 3 = 1 or 2
		if i % 3 != 0:
			triplet = get_triplet(string, i) # Triplet of form ("A", "B", "A")
			triplets_12.append((i,) + triplet) # Add idx to triplet: (idx, "A", "B", "A")
	# Sort the triplets (using stable bucket sort)
	triplets_12 = radix_3(string, triplets_12) # OLD: Works, but O(n*lg n) # triplets_12 = sorted(tuple(triplets_12), key = lambda tup: (tup[1], tup[2], tup[3]))
	# Extract the suffix indices
	sa_12 = [t[0] for t in triplets_12] # Use only the idx from triplets of form (idx, "A", "B", "A") in sa_12
	return sa_12

'''
MERGES SA_12 and SA_3 into the final suffix array
'''
def merge(string, sa_12, sa_3):
	# Create inverse suffix array for sa_12 (contains sa-ranks for given string idx)
    isa_12 = {sa_12[i]: i for i in range(len(sa_12))}
    sa = []
	# Let i be an index into sa_12 and j be an index into sa_3
    i = j = 0
	# Check which suffix is smaller: sa_12[i] or sa_3[j]
    while i < len(sa_12) and j < len(sa_3):
        if smaller(string, sa_12[i], sa_3[j], isa_12):
			# sa_12[i] is smaller
            sa.append(sa_12[i])
            i += 1
        else:
			# sa_3[j] is smaller
            sa.append(sa_3[j])
            j += 1
	# If there are suffixes left in one of the suffix arrays after the other is empty,
	# add these suffixes to sa
    sa.extend(sa_12[i:])
    sa.extend(sa_3[j:])
    return sa



########################################################
# HELPER FUNCTIONS
########################################################


# Returns the char in string at given index
# If we reach beyond string, return terminal sentinel ("pad" string with terminal sentinel)
def safe_get_char(string, index):
	if index >= len(string):
		return SENTINEL_terminal
	else:
		return string[index]


# Get triplet starting at pos i in string as a tuple ("A", "B", "A")
def get_triplet(string, i):
	return (safe_get_char(string, i), safe_get_char(string, i + 1), safe_get_char(string, i + 2))


# Gets alphabet of lex-names for the triples starting at given indices
# Input should be a string and suffix array (sa_12)
# Returns dictionary of form {('$', '$', '$'): 0, ('A', '$', '$'): 1, ...}
def get_alphabet(string, indices):
	alphabet = {('$', '$', '$'): 0} # never used: but gives alphabet correct size
	letter = 1
	for i in indices:
		triplet = get_triplet(string, i)
		# If triplet not seen before, add it to the alphabet with new, incremented lex name
		if triplet not in alphabet:
			alphabet[triplet] = letter
			letter += 1
	return alphabet


# Construct the u-string used by Skew Algo
# The u string consists of:
# Lex names for i mod 3 = 1, CENTRAL SENTINEL, lex names for i mod 3 = 2
def construct_u(string, alphabet):
	u = ""
	# Lex names for i mod 3 = 1
	for  i in range(1, len(string), 3):
		u += str(alphabet[get_triplet(string, i)])
	# Add sentinel
	u += SENTINEL_central
	# Lex names for i mod 3 = 2
	for  i in range(2, len(string), 3):
		u += str(alphabet[get_triplet(string, i)])
	return u


# Maps indices in u back to indices in the original input string
def map_u_to_string_index(i, m):
	# Left of the central sentinel
	if i < m:
		return 1 + 3 * i
	# Right of the central sentinel
	else:
		return 2 + 3 * (i - m - 1)


# Determines which suffix is smaller (i or j), using char comparisons and isa
def smaller(string, i, j, isa):
	# Char comparison first
	a = safe_get_char(string, i)
	b = safe_get_char(string, j)
	if a < b:
		return True
	if a > b:
		return False
	# If true, both suffix i and j will be in sa_12:
	# We can compare them using isa!
	if i % 3 != 0 and j % 3 != 0:
		return isa[i] < isa[j]
	# Else, we compare i + 1 and j + 1 recursively
	# We will max get a rec depth of 2, because of the skewed construction
	return smaller(string, i + 1, j + 1, isa)



########################################################
# STABLE BUCKET SORT AND RADIX SORT FUNCTIONS
########################################################

# Creates ordered dictionary of letters in alphabet from a string (always contains $ as key also):
# {A: ["ATC", "AAA", "ACG"], C: ["CTC", "CGG"], $: ["$AB"]}
def create_alphabet_dict(string):
	d = {key:[] for key in list(string) + [SENTINEL_terminal]}
	ordered_d = collections.OrderedDict(sorted(d.items()))
	return ordered_d

# Puts triplets into buckets, depending on idx (idx of 1, 2 or 3 is first, second and third char)
# Returns a list of [bucket 1] + [bucket 2] + ... + [bucket k]
def bucket_sort_triplet(string, triplets, idx):
	alpha_dict = create_alphabet_dict(string)
	for t in triplets:
		alpha_dict[t[idx]].append(t)
	res_list = sum(alpha_dict.values(), [])
	return res_list

# Sorts triplets using stable bucket sort on the three letters (last letter first, then second, then first)
def radix_3(string, triplets):
	triplets = bucket_sort_triplet(string, triplets, 3)
	triplets = bucket_sort_triplet(string, triplets, 2)
	triplets = bucket_sort_triplet(string, triplets, 1)
	return triplets

# Sorts suffixes of string with indices as in suffix_index according to their first character
def bucket_sort_first_char(string, suffix_indices):
	alpha_dict = create_alphabet_dict(string)
	for idx in suffix_indices:
		char = safe_get_char(string, idx)
		alpha_dict[char].append(idx)
	res_list = sum(alpha_dict.values(), [])
	return res_list



########################################################
# TEST CODE
########################################################


s = "mississippi$"
print(skew_rec(s))


