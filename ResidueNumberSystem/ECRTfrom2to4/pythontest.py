from math import gcd
from fractions import Fraction

def extended_gcd(a, b):
    """Return (g, x, y) such that a*x + b*y = g = gcd(a, b)."""
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = extended_gcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)

def mod_inverse(a, m):
    """Compute the modular inverse of a modulo m, if it exists."""
    g, x, _ = extended_gcd(a, m)
    if g != 1:
        raise ValueError("No modular inverse exists for given inputs.")
    return x % m

def explicit_crt_full_representation(u_list, p_list):
    """
    Given:
      - u_list = [u_1, u_2, ..., u_n]: partial RNS representation of u modulo p_1, p_2, ..., p_n
      - p_list = [p_1, p_2, ..., p_m]: full list of moduli

    It is assumed that |u| < P = p_1 * p_2 * ... * p_n.
    We produce u modulo each p_i, i=1..m, using an explicit-CRT inspired approach.
    """
    n = len(u_list)
    m = len(p_list)

    # Step 1: Compute P and partial products P_i
    P = 1
    for i in range(n):
        P *= p_list[i]

    P_i_list = [P // p_list[i] for i in range(n)]

    # Step 2: Compute q_i such that q_i * (P_i) ≡ 1 (mod p_i)
    q_list = []
    for i in range(n):
        q_i = mod_inverse(P_i_list[i], p_list[i])
        q_list.append(q_i)

    # Step 3: Compute t_i = (u_i * q_i) mod p_i, chosen as a signed representative
    # We'll convert to a signed representative in [-p_i/2, p_i/2).
    t_list = []
    for i in range(n):
        t = (u_list[i] * q_list[i]) % p_list[i]
        # Convert to signed representation
        if t > p_list[i] // 2:
            t = t - p_list[i]
        t_list.append(t)

    # Step 4: Compute alpha = sum(t_i / p_i) exactly using Fractions
    # alpha approximates u/P
    alpha = Fraction(0, 1)
    for i in range(n):
        alpha += Fraction(t_list[i], p_list[i])

    # Step 5: Compute r = round(alpha)
    # Since alpha is a Fraction, convert to float or implement rounding on fraction directly
    # We'll do a rational to nearest integer rounding:
    r_as_float = alpha.numerator / alpha.denominator
    r = int(round(r_as_float))

    # Step 6: Reconstruct u = sum(t_i * P_i) - P * r
    # This gives us an integer representative for u in [-P/2, P/2)
    # Note: u might be large, but Python can handle big integers.
    u_reconstructed = 0
    for i in range(n):
        u_reconstructed += t_list[i] * P_i_list[i]
    u_reconstructed -= P * r

    # Now we have a representative of u. Since |u| < P, this is the unique representative.
    # Step 7: Compute full RNS representation by reducing u modulo each p_j, j=1..m
    full_representation = []
    for j in range(m):
        # Compute u mod p_j and convert to signed representation in [-p_j/2, p_j/2)
        mod_val = u_reconstructed % p_list[j]
        if mod_val > p_list[j] // 2:
            mod_val -= p_list[j]
        full_representation.append(mod_val)

    return full_representation


if __name__ == "__main__":
    # Example usage:
    # Suppose we have moduli:
    p_list = [2147473409, 2147389441, 2147387393, 2147377153]  # p_1=7, p_2=11, p_3=13, p_4=17
    # We will use only the first n=2 moduli for the partial representation
    n = 2

    # Let's pick a random integer u such that |u|<P, with P=p_1*p_2=7*11=77
    u = 2147473410  # arbitrary choice, -38 < u=23 < 77 is true

    # Compute its partial RNS: u_1 mod p_1, u_2 mod p_2
    # Signed representatives:
    def signed_mod(x, m):
        val = x % m
        if val > m // 2:
            val -= m
        return val

    u_list = [signed_mod(u, p_list[i]) for i in range(n)]

    # Now we have (u_1 mod p_1, u_2 mod p_2) but we want the full representation under p_3=13 and p_4=17
    full_rep = explicit_crt_full_representation(u_list, p_list)

    # Check correctness:
    # full_rep[i] should be equal to u (mod p_list[i]), in signed form.
    print("Computed full RNS representation:", full_rep)
    # Verify correctness:
    for i, p in enumerate(p_list):
        expected = signed_mod(u, p)
        assert full_rep[i] == expected, f"Mismatch on modulus {p}: got {full_rep[i]}, expected {expected}"

    print("Verification passed! The computed representation matches the known integer u.")
