s = "mississippi"
SENTINEL = "$"

# Make sure we don't reach beyond string length when getting a char from the string
def safe_get_char(string, index):
	if index >= len(string):
		return SENTINEL
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
		u += alphabet[get_triplet(string, i)]
	# Add sentinel
	u += SENTINEL
	# Lex names for i mod 3 = 2
	for  i in range(2, len(string), 3):
		u += alphabet[get_triplet(string, i)]
	return u


# Map indices in u back to indices in the original input string
def map_u_to_string(i, m):
	if i < m:
		return 1 + 3 * i
	else:
		return 2 + 3 * (i - m - 1)

def compute_sa_12(string):
	triples_12 = []
	for i in range(len(string)):
		if i % 3 != 0:
			triplet = get_triplet(string, i)
			triples_12.append((i,) + triplet)
	triples_12 = sorted(tuple(triples_12), key = lambda tup: (tup[1], tup[2], tup[3]))
	sa_12 = [t[0] for t in triples_12]
	return sa_12


def skew_rec(string):
	sa_12 = compute_sa_12(string)
	alphabet = get_alphabet(string, sa_12)
	return alphabet



print(skew_rec(s))
#print(get_alphabet(s, [11,10,7,1,4,8,2,5]))



#def radix3(triple_list):









