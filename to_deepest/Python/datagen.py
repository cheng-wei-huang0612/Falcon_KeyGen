"""
This script performs the following steps:
1. Generate a list of prime moduli P.
2. For each modulus p in P, find a primitive 2048-th root of unity.
3. Generate required NTT-related arrays (omegas, mu_omegas).
4. Generate random polynomials f and g following a Gaussian distribution.
5. Write all generated data into data.py for later usage.

"""



logn = 10
bit_length = 30
num_of_moduli = 304
num_of_used_moduli_fg = [1, 1, 1, 2, 4, 8, 16, 28, 56, 108, 216]
num_of_used_moduli_FG = [1, 2, 4, 8, 12, 20, 40, 80, 152, 304, 216]



num_of_test = 10



# --- Generate prime moduli P ---

# Initialization function: compute the contents of primetest_pp[].
def mk_primetest():
    primetest_pp = []

    max_p = 65536   # ceil(sqrt(2^31))
    max_q = 256     # ceil(sqrt(max_p))
    tab = [True] * max_p
    tab[0] = False
    tab[1] = False
    for i in range(2, max_q):
        if tab[i]:
            for j in range(2*i, max_p, i):
                tab[j] = False
    for i in range(0, max_p):
        if tab[i]:
            primetest_pp.append(i)
    
    return primetest_pp

# Function to test if a number is prime.
def is_prime(p,primetest_pp):
    if p < 2:
        return p == 2
    if (p & 1) == 0: # p is even
        return False
    for j in range(0, len(primetest_pp)):
        t = primetest_pp[j]
        if t*t > p:
            return True
        if (p % t) == 0:
            return False
    return True

def mk_moduli(bit_length = 31, num_of_moduli = 241):

    primetest_pp = mk_primetest()


    moduli = []
    p = 2**bit_length + 1
    while len(moduli) < num_of_moduli:
        if is_prime(p,primetest_pp):
            moduli.append(p)
        p -= 2048

    return moduli




# --- Generate 2048 roots ---

def prime_factors(n):
    """Return the prime factorization of n as a list of prime factors (not necessarily unique)."""
    factors = []
    # Trial division (works for smaller n)
    # For large n, consider a more efficient method.
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1 if d == 2 else 2  # increment by 1 if 2, else by 2 for optimization
    if n > 1:
        factors.append(n)
    return factors

def is_primitive_root(g, p, factors):
    """Check if g is a primitive root modulo p, given the prime factors of p-1."""
    # p-1 = product of prime factors in 'factors'
    order = p - 1
    for f in set(factors):
        # If g^{(p-1)/f} ≡ 1 (mod p), g is not a primitive root
        if pow(g, order // f, p) == 1:
            return False
    return True

def find_primitive_root(p):
    """Find a primitive root modulo p."""
    # Factor p-1
    factors = prime_factors(p - 1)
    
    # Try small candidates for g
    for g in range(2, p):
        if is_primitive_root(g, p, factors):
            return g
    # In theory, we should always find one for prime p, but return None if not found.
    return None

def primitive_2048th_root_of_unity(p):
    # Check if 2048 divides p-1
    if (p - 1) % 2048 != 0:
        return 0  # No such root exists
    
    # Find a primitive root modulo p
    g = find_primitive_root(p)
    if g is None:
        # This should not happen for a prime p, but just in case:
        return 0
    
    # Compute g^{(p-1)/2048} mod p
    exponent = (p - 1) // 2048
    root = pow(g, exponent, p)
    return root





# --- Copmute P and omega_2048 (using the above functions) ---

P = mk_moduli(bit_length, num_of_moduli)
logn_top = 10
omega_2048 = [primitive_2048th_root_of_unity(p) for p in P]



# Here begin the computation
def MOD(a, q):
    result = a % q
    if result >= (q>>1):
        result = result - q
    return result


# Bit-reverse index
def bitrev(i):
    rev = 0
    for j in range(10):
        rev = (rev << 1) | (i & 1)
        i >>= 1
    return rev



n = 1 << logn_top


# polynomial sampler
import numpy as np
import math
def prng(sigma):
    return int(abs((np.random.normal(0, sigma, 1))[0]))

def Falcon_poly_generator(n):
    q = 12289
    sigma = math.sqrt(q/(2*n)) * 1.17

    f = []
    for i in range(n):
        f.append(prng(sigma))

    # Check for irreducibility

    # Check the inner product is small enough
    
    return f

def Mu(b,q):
    mu = 2*int(((b << 31)/(2*q))+0.5)
    return mu

def omegas_gen(omega,n,mod):
    """ 
    returns omegas[n] for each given omega and n
    """
    val = 1
    omegas = [1] * (n)

    for i in range(n):
        omegas[bitrev(i)] = val
        val = MOD((val * omega) , mod)

    return omegas

def mu_omegas_gen(omegas,mod):
    """
    returns mu_omegas[n] for given omegas and modulo mod
    """
    mu_omegas = [Mu(omega,mod) for omega in omegas]
    return mu_omegas


f_poly = Falcon_poly_generator(n)
g_poly = Falcon_poly_generator(n)


# 重新產生 omegas 與 mu_omegas，使用二維結構
omegas = []
mu_omegas = []
for i in range(num_of_moduli):
    current_omegas = omegas_gen(omega_2048[i], n, P[i])
    current_mu_omegas = mu_omegas_gen(current_omegas, P[i])
    omegas.append(current_omegas)
    mu_omegas.append(current_mu_omegas)

# 假設 omegas 與 mu_omegas 已經是二維陣列
# omegas[i]：對應第 i 個模數的 NTT omegas，長度為 n
# mu_omegas[i]：對應第 i 個模數的 mu_omegas，長度為 n

# 將二維陣列攤平成一維陣列
flat_omegas = []
flat_mu_omegas = []

for i in range(num_of_moduli):
    flat_omegas.extend(omegas[i])
    flat_mu_omegas.extend(mu_omegas[i])

with open("data.py", "w") as file:
    file.write("# This file consists of: logn_top, P, omega_2048, omegas, mu_omegas, f, g\n\n")
    file.write("logn_top = %d\n" % logn_top)

    # 輸出 P
    elems_per_line = 10
    file.write("# Prime moduli (P)\n")
    file.write("P = [\n")
    for i, val in enumerate(P):
        if i % elems_per_line == 0:
            file.write("    ")
        file.write("%d," % val)
        if (i + 1) % elems_per_line == 0:
            file.write("\n")
    if len(P) % elems_per_line != 0:
        file.write("\n")
    file.write("]\n\n")

    # 輸出 omega_2048
    file.write("# 2048-th roots of unity for each modulus\n")
    file.write("omega_2048 = [\n")
    for i, val in enumerate(omega_2048):
        if i % elems_per_line == 0:
            file.write("    ")
        file.write("%d," % val)
        if (i + 1) % elems_per_line == 0:
            file.write("\n")
    if len(omega_2048) % elems_per_line != 0:
        file.write("\n")
    file.write("]\n\n")

    # 輸出平坦化後的一維 omegas 陣列
    # 註解中標示出對應模數區間
    file.write("# Omegas array (bit-reversed order), flattened.\n")
    file.write("# Index ranges: For modulus i, omegas are in flat_omegas[i*n : (i+1)*n]\n")
    file.write("omegas = [\n")
    for i, val in enumerate(flat_omegas):
        if i % elems_per_line == 0:
            file.write("    ")
        file.write("%d," % val)
        if (i + 1) % elems_per_line == 0:
            file.write("\n")
    if len(flat_omegas) % elems_per_line != 0:
        file.write("\n")
    file.write("]\n\n")

    # 輸出平坦化後的一維 mu_omegas 陣列
    file.write("# mu_omegas array, flattened.\n")
    file.write("# Index ranges: For modulus i, mu_omegas are in flat_mu_omegas[i*n : (i+1)*n]\n")
    file.write("mu_omegas = [\n")
    for i, val in enumerate(flat_mu_omegas):
        if i % elems_per_line == 0:
            file.write("    ")
        file.write("%d," % val)
        if (i + 1) % elems_per_line == 0:
            file.write("\n")
    if len(flat_mu_omegas) % elems_per_line != 0:
        file.write("\n")
    file.write("]\n\n")

    # 輸出 f
    file.write("# Polynomial f\n")
    file.write("f = [\n")
    for i, val in enumerate(f_poly):
        if i % elems_per_line == 0:
            file.write("    ")
        file.write("%d," % val)
        if (i + 1) % elems_per_line == 0:
            file.write("\n")
    if len(f_poly) % elems_per_line != 0:
        file.write("\n")
    file.write("]\n\n")

    # 輸出 g
    file.write("# Polynomial g\n")
    file.write("g = [\n")
    for i, val in enumerate(g_poly):
        if i % elems_per_line == 0:
            file.write("    ")
        file.write("%d," % val)
        if (i + 1) % elems_per_line == 0:
            file.write("\n")
    if len(g_poly) % elems_per_line != 0:
        file.write("\n")
    file.write("]\n")



    # Write the same data to data.c and data.h

    # Write data.h
    with open("../C/data.h", "w") as file:
        file.write("#ifndef DATA_H\n")
        file.write("#define DATA_H\n\n")
        file.write("#include <stdint.h>\n\n")
        file.write("extern const int logn_top;\n")
        file.write("extern const uint32_t P[%d];\n" % len(P))
        file.write("extern const uint32_t omega_2048[%d];\n" % len(omega_2048))
        file.write("extern const uint32_t omegas[%d];\n" % len(flat_omegas))
        file.write("extern const uint32_t mu_omegas[%d];\n" % len(flat_mu_omegas))
        file.write("extern const int32_t f[%d];\n" % len(f_poly))
        file.write("extern const int32_t g[%d];\n" % len(g_poly))
        file.write("\n#endif // DATA_H\n")

    # Write data.c
    with open("../C/data.c", "w") as file:
        file.write("#include \"data.h\"\n\n")
        file.write("const int logn_top = %d;\n\n" % logn_top)

        # Write P
        file.write("const uint32_t P[%d] = {\n" % len(P))
        for i, val in enumerate(P):
            if i % elems_per_line == 0:
                file.write("    ")
            file.write("%d," % val)
            if (i + 1) % elems_per_line == 0:
                file.write("\n")
        if len(P) % elems_per_line != 0:
            file.write("\n")
        file.write("};\n\n")

        # Write omega_2048
        file.write("const uint32_t omega_2048[%d] = {\n" % len(omega_2048))
        for i, val in enumerate(omega_2048):
            if i % elems_per_line == 0:
                file.write("    ")
            file.write("%d," % val)
            if (i + 1) % elems_per_line == 0:
                file.write("\n")
        if len(omega_2048) % elems_per_line != 0:
            file.write("\n")
        file.write("};\n\n")

        # Write flat_omegas
        file.write("const uint32_t omegas[%d] = {\n" % len(flat_omegas))
        for i, val in enumerate(flat_omegas):
            if i % elems_per_line == 0:
                file.write("    ")
            file.write("%d," % val)
            if (i + 1) % elems_per_line == 0:
                file.write("\n")
        if len(flat_omegas) % elems_per_line != 0:
            file.write("\n")
        file.write("};\n\n")

        # Write flat_mu_omegas
        file.write("const uint32_t mu_omegas[%d] = {\n" % len(flat_mu_omegas))
        for i, val in enumerate(flat_mu_omegas):
            if i % elems_per_line == 0:
                file.write("    ")
            file.write("%d," % val)
            if (i + 1) % elems_per_line == 0:
                file.write("\n")
        if len(flat_mu_omegas) % elems_per_line != 0:
            file.write("\n")
        file.write("};\n\n")

        # Write f
        file.write("const int32_t f[%d] = {\n" % len(f_poly))
        for i, val in enumerate(f_poly):
            if i % elems_per_line == 0:
                file.write("    ")
            file.write("%d," % val)
            if (i + 1) % elems_per_line == 0:
                file.write("\n")
        if len(f_poly) % elems_per_line != 0:
            file.write("\n")
        file.write("};\n\n")

        # Write g
        file.write("const int32_t g[%d] = {\n" % len(g_poly))
        for i, val in enumerate(g_poly):
            if i % elems_per_line == 0:
                file.write("    ")
            file.write("%d," % val)
            if (i + 1) % elems_per_line == 0:
                file.write("\n")
        if len(g_poly) % elems_per_line != 0:
            file.write("\n")
        file.write("};\n")