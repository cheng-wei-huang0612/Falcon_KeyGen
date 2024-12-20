from math import gcd
from fractions import Fraction

# Given moduli
p = [1073707009, 1073698817, 1073692673, 1073682433]
n = 2  # number of moduli for the partial representation
m = len(p)

# Extended Euclidean Algorithm to find inverses
def extended_gcd(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = extended_gcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)

def mod_inverse(a, m):
    g, x, _ = extended_gcd(a, m)
    if g != 1:
        raise ValueError("Modular inverse does not exist.")
    return x % m

# Precompute P = p_1 * p_2 * ... * p_n
P = 1
for i in range(n):
    P *= p[i]

# Compute P_i and q_i
P_i_list = [P // p[i] for i in range(n)]
q_list = [mod_inverse(P_i_list[i], p[i]) for i in range(n)]

# We will also store the moduli array and related constants in data.c
# We'll write to data.c and data.h
# data.h will contain forward declarations, data.c will have definitions.

with open("data.h", "w") as fh:
    fh.write("#ifndef DATA_H\n")
    fh.write("#define DATA_H\n\n")
    fh.write("#include <stdint.h>\n\n")

    fh.write("extern const int n;\n")
    fh.write("extern const int m;\n")
    fh.write("extern const uint32_t p[];\n")
    fh.write("extern const uint64_t P;\n")
    fh.write("extern const uint64_t P_i[];\n")
    fh.write("extern const uint64_t q[];\n")
    fh.write("#endif // DATA_H\n")

with open("data.c", "w") as fc:
    fc.write("#include \"data.h\"\n\n")
    fc.write("const int n = %d;\n" % n)
    fc.write("const int m = %d;\n" % m)

    # p array
    fc.write("const uint32_t p[] = {")
    fc.write(",".join(str(x) for x in p))
    fc.write("};\n")

    # P
    fc.write("const uint64_t P = %dULL;\n" % P)

    # P_i
    fc.write("const uint64_t P_i[] = {")
    fc.write(",".join(str(P_i_list[i]) + "ULL" for i in range(n)))
    fc.write("};\n")

    # q
    fc.write("const uint64_t q[] = {")
    fc.write(",".join(str(q_list[i]) + "ULL" for i in range(n)))
    fc.write("};\n")
