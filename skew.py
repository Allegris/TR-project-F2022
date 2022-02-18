'''
	Skew algorithm for computing the suffix array for a string in time O(n)
'''


SENTINEL_terminal = "$"
SENTINEL_central = "#"


# If we reach beyond string, return terminal sentinel ("pad" string with terminal sentinel)
def safe_get_char(string, index):
	if index >= len(string):
		return SENTINEL_terminal
	else:
		return string[index]


# Get triplet starting at pos i in string
def get_triplet(string, i):
	return (safe_get_char(string, i), safe_get_char(string, i + 1), safe_get_char(string, i + 2))


# Get alphabet of lex-names for the triples starting at given indices
# Input should be a string and suffix array (sa_12)
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


# Map indices in u back to indices in the original input string
def map_u_to_string_index(i, m):
	# Left of the central sentinel
	if i < m:
		return 1 + 3 * i
	# Right of the central sentinel
	else:
		return 2 + 3 * (i - m - 1)

# Computes the suffix array for the suffixes starting at indices != 0 mod 3
def compute_sa_12(string):
	triplets_12 = []
	# Create triplets for all indices != 0 mod 3 in the string
	for i in range(len(string)):
		if i % 3 != 0:
			triplet = get_triplet(string, i)
			triplets_12.append((i,) + triplet)
	# Sort the triplets
	triplets_12 = sorted(tuple(triplets_12), key = lambda tup: (tup[1], tup[2], tup[3]))
	# Extract the suffix indices
	sa_12 = [t[0] for t in triplets_12]
	return sa_12


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


# Merges sa_12 and sa_3 into the final suffix array
def merge(string, sa_12, sa_3):
	# Create inverse suffix array for sa_12
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


# Returns the suffix array of a string
def skew_rec(string):
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

	# Construct sa_3 from sa_12
	sa_3 = []
	# Special case: if last index in string is 0 mod 3, then this suffix should be first in sa_3, since it is the shortest suffix
	# It will be sorted based on the first (and in this case only) character afterwards
	if len(string) % 3 == 1:
		sa_3.append(len(string) - 1)
	# All strings mod 3 = 0 are just strings mod 3 = 1, with one char in front, so:
	# Use the order found in sa_12, and then stable sort on the first character afterwards
	sa_3 += [i - 1 for i in sa_12 if i % 3 == 1]
	# Stable sort on first character
	sa_3 = sorted(sa_3, key = lambda i : safe_get_char(string, i)[:1])

	#print(string, "sa_12: ", sa_12)
	#print(string, "SA_3: ", sa_3)
	return merge(string, sa_12, sa_3)



### TEST CODE ###

s = "mississippi$"
print(skew_rec(s))


