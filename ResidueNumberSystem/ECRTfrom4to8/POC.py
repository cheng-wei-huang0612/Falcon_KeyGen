from math import gcd
from fractions import Fraction
import random
from data import n, m, p, P_mod, P_i_mod, q

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

def signed_repr(val, mod):
    # Convert val mod mod into range [-mod/2, mod/2)
    val = val % mod
    if val > mod // 2:
        val -= mod
    return val


def crt(moduli, residues):
    # A simple implementation of CRT for pairwise coprime moduli
    from functools import reduce
    M = 1
    for mm in moduli:
        M *= mm
    x = 0
    for mm, r in zip(moduli, residues):
        M_i = M // mm
        inv = pow(M_i, -1, mm)  # Python 3.8+ syntax for modular inverse
        x += r * M_i * inv

    return signed_repr(x, M)


def ecrt(u_partial, p, n, m, q, P_mod, P_i_mod):
    # Steps:
    # 1. Compute t_i
    t = []
    for i in range(n):
        pi = p[i]
        ui = u_partial[i]
        # t_i = (u_i * q_i) mod p_i, in signed form
        tt = (ui * q[i]) % pi
        tt = signed_repr(tt, pi)
        t.append(tt)

    # 2. Compute alpha and r
    # alpha = sum(t_i / p_i)
    alpha = Fraction(0, 1)
    for i in range(n):
        alpha += Fraction(t[i], p[i])
    r = round(alpha)  # exact rounding using Fraction

    # 3. Compute full RNS
    u_full = []
    for j in range(m):
        pj = p[j]
        # sum_tiPi = sum(t_i * (P_i mod p_j))
        sum_tiPi = 0
        for i in range(n):
            val = (t[i] * P_i_mod[i*m + j]) % pj
            sum_tiPi = (sum_tiPi + val) % pj

        Pr = (P_mod[j] * r) % pj

        val = (sum_tiPi - Pr) % pj
        val = signed_repr(val, pj)
        u_full.append(val)

    return u_full






# Pick a random u with |u| < P by selecting partial residues randomly
# in [-p[i]/2, p[i]/2)
u_partial = []
for i in range(n):
    ui = random.randint(-p[i]//2, p[i]//2 - 1)
    u_partial.append(ui)

# Compute u from partial representation using CRT
u= crt(p[:n], u_partial)

print(f"u: {u}")
print(f"Partial: {u_partial}")

# Use ECRT to get full representation
u_full = ecrt(u_partial, p, n, m, q, P_mod, P_i_mod)
print(f"Full: {u_full}")

# Verification:
# Check all residues match u:
for j in range(m):
    expected = u % p[j]
    if expected > p[j]//2:
        expected -= p[j]
    if expected != u_full[j]:
        print(f"Mismatch at modulus p[{j}]={p[j]}: got {u_full[j]}, expected {expected}")
        print("M"*100)
        break
else:
    print("All residues match.")

# Additionally, verify partial stability:
for i in range(n):
    if u_full[i] != u_partial[i]:
        print(f"Partial mismatch at i={i}, expected {u_partial[i]}, got {u_full[i]}")
        break
else:
    print("Partial representation stable.")
