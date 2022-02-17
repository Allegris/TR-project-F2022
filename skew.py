s = "mississippi$"
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
def get_alphabet(string, indices):
	alphabet = {}
	letter = 0
	for i in indices:
		triplet = get_triplet(string, i)
		if triplet not in alphabet:
			alphabet[triplet] = letter
			letter += 1
	return alphabet
