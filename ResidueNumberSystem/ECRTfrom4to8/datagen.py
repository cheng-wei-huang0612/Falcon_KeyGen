from math import gcd

# Given moduli
p = [
    1073707009,
    1073698817,
    1073692673,
    1073682433,
    1073668097,
    1073655809,
    1073651713,
    1073643521
]

n = 4
m = 8

def extended_gcd(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = extended_gcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)

def mod_inverse(a, m):
    g, x, _ = extended_gcd(a, m)
    if g != 1:
        raise ValueError("No modular inverse")
    return x % m

# Compute P = p_1 * p_2 * ... * p_n
P = 1
for i in range(n):
    P *= p[i]

# Compute P_i and q_i
P_i_list = [P // p[i] for i in range(n)]
q_list = [mod_inverse(P_i_list[i] % p[i], p[i]) for i in range(n)]

# Precompute P mod p_j and P_i mod p_j
P_mod = [P % p_j for p_j in p]

P_i_mod = []
for i in range(n):
    row = [(P_i_list[i] % p_j) for p_j in p]
    P_i_mod.append(row)

# Write data.h and data.c
with open("data.h", "w") as fh:
    fh.write("#ifndef DATA_H\n")
    fh.write("#define DATA_H\n\n")
    fh.write("#include <stdint.h>\n\n")

    fh.write("extern const int n;\n")
    fh.write("extern const int m;\n")
    fh.write("extern const uint32_t p[];\n")
    fh.write("extern const uint32_t P_mod[];\n")
    fh.write("extern const uint32_t P_i_mod[];\n")
    fh.write("extern const uint32_t q[];\n")

    fh.write("#endif // DATA_H\n")

with open("data.c", "w") as fc:
    fc.write("#include \"data.h\"\n\n")
    fc.write("const int n = %d;\n" % n)
    fc.write("const int m = %d;\n" % m)

    # p array
    fc.write("const uint32_t p[] = {")
    fc.write(",".join(str(x) for x in p))
    fc.write("};\n")

    # P_mod array
    fc.write("const uint32_t P_mod[] = {")
    fc.write(",".join(str(x) for x in P_mod))
    fc.write("};\n")

    # q array
    fc.write("const uint32_t q[] = {")
    fc.write(",".join(str(q_list[i]) for i in range(n)))
    fc.write("};\n")

    # Flatten P_i_mod for storage
    fc.write("const uint32_t P_i_mod[] = {")
    for i in range(n):
        fc.write(",".join(str(x) for x in P_i_mod[i]))
        if i < n-1:
            fc.write(",")
    fc.write("};\n")

# Write data.py to store the same data in Python format
with open("data.py", "w") as fpy:
    fpy.write("# This file was auto-generated.\n\n")
    fpy.write(f"n = {n}\n")
    fpy.write(f"m = {m}\n")

    fpy.write("p = [")
    fpy.write(",".join(str(x) for x in p))
    fpy.write("]\n")

    fpy.write("P_mod = [")
    fpy.write(",".join(str(x) for x in P_mod))
    fpy.write("]\n")

    fpy.write("q = [")
    fpy.write(",".join(str(x) for x in q_list))
    fpy.write("]\n")

    # Flatten P_i_mod for Python
    # We'll store it as a single list, but you can reshape as needed
    flat_P_i_mod = []
    for i in range(n):
        flat_P_i_mod.extend(P_i_mod[i])

    fpy.write("P_i_mod = [")
    fpy.write(",".join(str(x) for x in flat_P_i_mod))
    fpy.write("]\n")
