s = "mississippi$"
SENTINEL = "#"

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
# Input should be string and suffix array
def get_alphabet(string, indices):
	alphabet = {('#', '#', '#', 0)}
	letter = 0
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





