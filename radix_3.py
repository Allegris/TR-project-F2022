#################################################
# Radix sort on first 3 symbols in strings
#################################################

import collections


# Create ordered dictionary of letters in alphabet from a string:
# {A: ["ATC", "AAA", "ACG"], C: ["CTC", "CGG"]}
def create_alphabet_dict(string):
	d = {key:[] for key in list(string)}
	ordered_d = collections.OrderedDict(sorted(d.items()))
	return ordered_d

# Put triplets into buckets, depending on idx (idx of 1, 2 or 3 is first, second and third char)
# Returns a list of [bucket 1] + [bucket 2] + ... + [bucket k]
def bucket_sort_triplet(string, triplets, idx):
	alpha_dict = create_alphabet_dict(string)
	for t in triplets:
		alpha_dict[t[idx]].append(t)
	res_list = sum(alpha_dict.values(), [])
	return res_list

# Sort triplets using stable bucket sort on the three letters (last letter first, then second, then first)
def radix_3(string, triplets):
	triplets = bucket_sort_triplet(string, triplets, 3)
	triplets = bucket_sort_triplet(string, triplets, 2)
	triplets = bucket_sort_triplet(string, triplets, 1)
	return triplets

def bucket_sort_string(string, suffix_indices, idx):
	alpha_dict = create_alphabet_dict(string)
	for idx in suffix_indices:
		char = safe_get_char(string, idx)
		alpha_dict[char].append(idx)
	res_list = sum(alpha_dict.values(), [])
	return res_list



#### TEST CODE

s = "mississippi$"

triplets = [(0, "m", "i", "s"), (1, "i", "s", "s"), (2, "s", "s", "i"), (3, "s", "i", "s"), (4, "s", "s", "i"), (5, "s", "i", "p"), (6, "i", "p", "p"), (7, "p", "p", "i"), (8, "p", "i", "$"), (9, "i", "$", "$"), (10, "$", "$", "$")]

print(radix_3(s, triplets))
