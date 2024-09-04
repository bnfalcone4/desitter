import itertools

from collections import Counter

N = 500

# Define the ranges for n, l, and r
n_range = range(0, N)  # Replace with your range for n
l_range = range(0, N)  # Replace with your range for l
r_range = range(0, N)  # Replace with your range for r

# Generate all combinations
combinations = list(itertools.product(n_range, l_range, r_range))

sums = []

for n, l, r in combinations:

    sums.append(n**2 + l**2 + r**2)

print(len(sums))

sums_set = set(sums)

print(len(sums_set))

sums_dict = dict(Counter(sums))

print(len(sums_dict))