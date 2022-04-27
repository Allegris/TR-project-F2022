
SENTINEL = 0

################################################################################
# SKEW ALGORITHM without central sentinel - COMPUTE SUFFIX ARRAY IN TIME O(N)
################################################################################

'''
Computes the suffix array for a string in time O(n)

Input is a list of integers:
Our input is indirectly the string, e.g. "ACGTAA0", but in the form of a integer list,
as we get from "map_string_to_ints(string)[1]", e.g. "[1, 2, 3, 4, 1, 1, 0]"
Remember to add the sentinel, 0, to the input

Returns the suffix array as a sorted list of indices into the string
'''
def skew_rec(x, alpha_size):
	### STEP 1: COMPUTE SA_12 ###
	sa_12 = compute_sa_12(x, alpha_size)
	# Alphabet dict, {triplet: lex-name}
	new_alphabet = get_triplet_alphabet(x, sa_12)

	# If we have ducplicate(s) in sa_12, recurse:
	if len(sa_12) > len(new_alphabet):
		# Compute u string and recurse to find its suffix array
		u = construct_u(x, new_alphabet)
		sa_u = skew_rec(u, len(new_alphabet))
		# m is where suffixes mod 1 and mod 2 separate in u
		m = (len(sa_u) + 1) // 2
		sa_12 = [map_u_to_string_index(i, m) for i in sa_u]

	### STEP 2: COMPUTE SA_3 from SA_12 ###
	sa_3 = []
	# Special case: if last index in string is 0 mod 3,
	# then this suffix should be first in sa_3, since it is the shortest suffix
	# It will be sorted based on the first (and in this case only) character afterwards
	if len(x) % 3 == 1:
		sa_3.append(len(x) - 1)
	# All strings mod 3 = 0 are just strings mod 3 = 1, with one char in front, so:
	# Use the order found in sa_12, and then stable sort on the first character afterwards
	sa_3 += [i - 1 for i in sa_12 if i % 3 == 1]
	# Stable sort on first character
	sa_3 = bucket_sort_first_char(x, sa_3, alpha_size)

	### STEP 3: MERGE SA_12 and SA_3 into the final suffix array ###
	return merge(x, sa_12, sa_3)



################################################################################
# Functions for computing SA_12, the u-string and merging SA_12 and SA_3
################################################################################

'''
Computes SA_12, i.e. the suffix array for the suffixes starting at indices != 0 mod 3
'''
def compute_sa_12(x, alpha_size):
	triplets_12 = []
	# Create triplets for all indices != 0 mod 3 in the string
	for i in range(len(x)):
		# Only use mod 3 = 1 or 2
		if i % 3 != 0:
			triplet = get_triplet(x, i) # Triplet of form (1, 2, 3)
			triplets_12.append((i,) + triplet) # Add idx to triplet: (idx, 1, 2, 3)
	# Sort the triplets (using stable bucket sort)
	triplets_12 = radix_3(x, triplets_12, alpha_size)
	# Extract the suffix indices
	sa_12 = [t[0] for t in triplets_12] # Use only the idx from triplets of form (idx, 1, 2, 3) in sa_12
	return sa_12


'''
Constructs the u-string (as an integer list)
The u string consists of:
Lex-names for i mod 3 = 1, lex-names for i mod 3 = 2
We do NOT add a central (#) sentinel, since it is superfluous when we use a terminal sentinel (0)
'''
def construct_u(x, alphabet):
	u = []
	# Lex names for i mod 3 = 1
	for  i in range(1, len(x), 3):
		u.append(alphabet[get_triplet(x, i)])
	# Add sentinel
	# u.append(SENTINEL)
	# Lex names for i mod 3 = 2
	for  i in range(2, len(x), 3):
		u.append(alphabet[get_triplet(x, i)])
	return u

'''
Maps u indices back to indices in the original string
'''
def map_u_to_string_index(i, m):
	return 1 + 3 * i if i < m else 2 + 3 * (i - m)


'''
Merges SA_12 and SA_3 into the final suffix array
'''
def merge(x, sa_12, sa_3):
	# Create inverse suffix array for sa_12 (contains sa-ranks for given x idx)
    isa_12 = {sa_12[i]: i for i in range(len(sa_12))}
    sa = []
	# Let i be an index into sa_12 and j be an index into sa_3
    i = j = 0
	# Check which suffix is smaller: sa_12[i] or sa_3[j]
    while i < len(sa_12) and j < len(sa_3):
        if smaller(x, sa_12[i], sa_3[j], isa_12):
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


################################################################################
# HELPER FUNCTIONS
################################################################################

'''
Determines which suffix is lexicographically smaller, i or j, using char comparisons and isa
'''
def smaller(string, i, j, isa):
	# Char comparison first
	a = safe_string_idx(string, i)
	b = safe_string_idx(string, j)

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


'''
Returns a mapping between numbers and letter, eg. {0: "0", 1: "A", 2: "C", 3: "G", 4: "T"}, 0 is the sentinel
and a list of the string, eg. for "ACGT" we get [1, 2, 3, 4]
'''
def map_string_to_ints(string):
	letters = ''.join(set(string))
	letters = sorted(letters)
	num = 0
	num_to_letter_dict = {}
	letter_to_num_dict = {}
	for letter in letters:
		num_to_letter_dict[num] = letter
		letter_to_num_dict[letter] = num
		num += 1
	num_ls = []
	for s in string:
		num_ls.append(letter_to_num_dict[s])
	return num_to_letter_dict, num_ls


'''
Returns the integer in x at given index
If we reach beyond x, return sentinel, 0 ("pad" string with sentinel)
'''
def safe_string_idx(x, index):
	if index >= len(x):
		return SENTINEL
	else:
		return x[index]


'''
Returns triplet starting at pos i in x as a tuple, e.g. (1, 2, 3)
Pads with sentinel (0), if we reach beyond x
'''
def get_triplet(x, i):
	return (safe_string_idx(x, i), safe_string_idx(x, i + 1), safe_string_idx(x, i + 2))


'''
Gets alphabet of lex-names for the triples starting at given indices
Returns dictionary of form {triplet: lex-name}, e.g. {(1, 2, 1): 1, (1, 3, 1): 2, (2, 1, 1): 3...}
'''
def get_triplet_alphabet(x, indices):
	alphabet = {}
	letter = 1
	for i in indices:
		triplet = get_triplet(x, i)
		# If triplet not seen before, add it to the alphabet with new, incremented lex name
		if triplet not in alphabet:
			alphabet[triplet] = letter
			letter += 1
	return alphabet



########################################################
# STABLE BUCKET SORT AND RADIX SORT FUNCTIONS
########################################################

'''
Puts triplets into buckets, depending on idx (idx of 1, 2 or 3 is first, second and third char)
Returns a list of combined buckets in lex-order: [bucket 1, bucket 2, ..., bucket k]
'''
def bucket_sort_triplet(x, triplets, alpha_size, idx):
	# Bucket for
	buckets = [[]]
	buckets += [[] for i in range(alpha_size)]
	for t in triplets:
		buckets[t[idx]].append(t)
	#combined_buckets = [list(itertools.chain(*sub)) for sub in buckets]
	combined_buckets = [item for sublist in buckets for item in sublist]
	return combined_buckets


'''
Sorts triplets using stable bucket sort on the three symbols
(last symbol first, then second, then first)
Returns a list of the sorted triplets
'''
#
def radix_3(x, triplets, alpha_size):
	triplets = bucket_sort_triplet(x, triplets, alpha_size, 3)
	triplets = bucket_sort_triplet(x, triplets, alpha_size, 2)
	triplets = bucket_sort_triplet(x, triplets, alpha_size, 1)
	return triplets


'''
Sorts suffixes of x with indices as in suffix_index according to their first character
Returns a list of the sorted suffix indices
'''
def bucket_sort_first_char(x, suffix_indices, alpha_size):
	buckets = [[]]
	buckets += [[] for i in range(alpha_size)]
	for idx in suffix_indices:
		char = safe_string_idx(x, idx)
		buckets[char].append(idx)
	combined_buckets = [item for sublist in buckets for item in sublist]
	return combined_buckets



########################################################
# TEST CODE
########################################################

'''
# TEST STRINGS
#string = "ACGTAA0"
#string = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa0"
#string = "mississippi0"
#string = "ABCABC0"
string = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa0"

alpha, ints = map_string_to_ints(string)
#print(string)
#print(alpha)
#print(ints) # <-- ints is our input ("indirect input string") to skew_rec

sa = skew_rec(ints, len(alpha))
print("SA:", sa)
# Print all the suffixes in order of SA
for s in sa:
	print(string[s:])

'''
