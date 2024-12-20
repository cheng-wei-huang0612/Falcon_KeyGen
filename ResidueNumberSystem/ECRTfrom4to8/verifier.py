import sys


def crt(moduli, residues):
    # A simple implementation of CRT for pairwise coprime moduli
    # Assumes moduli are pairwise coprime
    from functools import reduce
    M = 1
    for m in moduli:
        M *= m
    x = 0
    for m, r in zip(moduli, residues):
        M_i = M // m
        # Compute inverse of M_i mod m
        inv = pow(M_i, -1, m)  # Python 3.8+ syntax
        x += r * M_i * inv
    return x % M, M



# Define n, m, and p as in data.h
n = 4
m = 8
p = [
1073707009,1073698817,1073692673,1073682433,1073668097,1073655809,1073651713,1073643521
]


input_data = """
Partial: -199638979 14365605 -149137902 -199514714
Full: -199638979 14365605 -149137902 -199514714 -98845824 230738384 198024729 -156150512
""".strip()


# Parse the input_data
partial_line = None
full_line = None
for line in input_data.splitlines():
    line = line.strip()
    if line.startswith("Partial:"):
        partial_line = line
    elif line.startswith("Full:"):
        full_line = line

if partial_line is None or full_line is None:
    print("Error: Could not find partial or full lines in input_data")
    sys.exit(1)

# Parse arrays
partial_vals = partial_line.split()[1:] # skip the word "Partial:"
full_vals = full_line.split()[1:]       # skip the word "Full:"

if len(partial_vals) != n:
    print("Error: partial array length mismatch")
    sys.exit(1)

if len(full_vals) != m:
    print("Error: full array length mismatch")
    sys.exit(1)

partial = list(map(int, partial_vals))
full_rep = list(map(int, full_vals))

# Adjust partial residues to [0, p[i]-1] for CRT
residues_n_nonneg = []
for i in range(n):
    val = partial[i]
    if val < 0:
        val += p[i]
    residues_n_nonneg.append(val)

# Compute u mod P using CRT from the partial residues
u_mod_P, mod_product = crt(p[:n], residues_n_nonneg)
u_mod_P_int = int(u_mod_P)

# Verify full representation
for j in range(m):
    pj = p[j]
    val = u_mod_P_int % pj
    if val > pj//2:
        val -= pj
    if val != full_rep[j]:
        print(f"Mismatch at modulus p[{j}]={pj}: expected {val}, got {full_rep[j]}")
        sys.exit(1)

print("Verification passed: full representation matches computed u.")
