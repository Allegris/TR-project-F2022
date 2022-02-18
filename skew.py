s = "mississippi$" # must not include sentinel
SENTINEL_terminal = "$"
SENTINEL_central = "#"

# Make sure we don't reach beyond string length when getting a char from the string
def safe_get_char(string, index):
	if index >= len(string):
		return SENTINEL_terminal
	else:
		return string[index]

# Get triplet starting at pos i in string
def get_triplet(string, i):
	return (safe_get_char(string, i), safe_get_char(string, i + 1), safe_get_char(string, i + 2))



# Get alphabet of lex-names for the triples starting at given indices
# Input should be string and suffix array (sa_12)
def get_alphabet(string, indices):
	alphabet = {('$', '$', '$'): 0}
	letter = 1
	for i in indices:
		triplet = get_triplet(string, i)
		if triplet not in alphabet:
			alphabet[triplet] = letter
			letter += 1
	return alphabet


# The u string consists of:
# Lex names for i mod 3 = 1, SENTINEL, lex names for i mod 3 = 2
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
	if i < m:
		return 1 + 3 * i
	else:
		return 2 + 3 * (i - m - 1)

def compute_sa_12(string):
	triplets_12 = []
	for i in range(len(string)):
		if i % 3 != 0:
			triplet = get_triplet(string, i)
			triplets_12.append((i,) + triplet)
	triplets_12 = sorted(tuple(triplets_12), key = lambda tup: (tup[1], tup[2], tup[3]))
	#print("Triplets: ", triplets_12)
	sa_12 = [t[0] for t in triplets_12]
	return sa_12

def construct_isa(sa):
	isa = {}#[None] * len(sa)
	for i in range(0, len(sa)):
		isa[sa[i]] = i
	return isa


def less(string, i, j, isa):
	a = safe_get_char(string, i)
	b = safe_get_char(string, j)
	if a < b:
		return True
	if a > b:
		return False
	if i % 3 != 0 and j % 3 != 0:
		return isa[i] < isa[j]
	return less(string, i + 1, j + 1, isa)



def merge(string, sa_12, sa_3):
    isa_12 = {sa_12[i]: i for i in range(len(sa_12))}
    sa = []
    i = j = 0
    while i < len(sa_12) and j < len(sa_3):
        if less(string, sa_12[i], sa_3[j], isa_12):
            sa.append(sa_12[i])
            i += 1
        else:
            sa.append(sa_3[j])
            j += 1
    sa.extend(sa_12[i:])
    sa.extend(sa_3[j:])
    return sa



def skew_rec(string):
	#print("BEGIN: ", string)
	sa_12 = compute_sa_12(string)
	alphabet = get_alphabet(string, sa_12)
	#print("alphabet: ", alphabet)
	#print("sa_12: ", sa_12)

	# If all unique in sa_12, don't recurse (if they have equal length, then we have a duplicate in sa_12, because alphabet contains sentinel)
	if len(sa_12) >= len(alphabet):
		u = construct_u(string, alphabet)
		#print("u string: ", u)
		sa_u = skew_rec(u)
		m = len(u) // 2
		#print(string, "sa_u: ", sa_u)
		sa_12 = [map_u_to_string_index(i, m) for i in sa_u if i != m] # i == m is central sentinel (#)
		#print(string, m, "sa_12: ", sa_12)


	# Construct sa_3 from sa_12
	sa_3 = []
	# Special case if last index in string is 0 mod 3, then this suffix should be first in sa_3, since it is the shortest suffix
	# (it will be sorted based on the first (and in this case only) character later)
	if len(string) % 3 == 1:
		sa_3.append(len(string) - 1)
	# Use the order found in sa_12
	sa_3 += [i - 1 for i in sa_12 if i % 3 == 1]
	# Stable sort on first character
	sa_3 = sorted(sa_3, key = lambda x : string[x][:1])

	print(string, "sa_12: ", sa_12)
	print(string, "SA_3: ", sa_3)
	return merge(string, sa_12, sa_3)


print(skew_rec(s))



#print(sorted(["aaa", "aba", "aaa"], key=lambda x:x[0:1]))


