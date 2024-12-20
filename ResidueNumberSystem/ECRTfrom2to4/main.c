#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "data.h"


// Compute signed representative modulo p: a value in [-p/2, p/2)
static int64_t signed_mod(int64_t x, uint32_t m) {
    int64_t val = (int64_t)(x % m);
    if (val < 0) val += m; // ensure nonnegative
    if (val > (int64_t)m/2) val -= (int64_t)m;
    return val;
}

// signed_mod_mul: Compute (a*b) mod m, returning a signed representative in [-m/2, m/2)
static int64_t signed_mod_mul(int64_t a, int64_t b, int64_t m) {
    // First reduce a and b mod m, ensuring nonnegative
    int64_t am = a % m;
    if (am < 0) am += m;
    int64_t bm = b % m;
    if (bm < 0) bm += m;

    __int128 prod = (__int128)am * (__int128)bm;
    __int128 res = prod % m;
    int64_t r = (int64_t)res;

    // Convert to signed representative
    if (r > m/2) r -= m;
    return r;
}

// Given the partial representation (u_1 mod p_1, u_2 mod p_2, ..., u_n mod p_n)
// contained in the array u (int32_t), we will reconstruct the full RNS representation
// and store it back into the same array u, extending it to (u_1 mod p_1, ..., u_m mod p_m).
static void ecrt_full_representation(int32_t *u) {
    // u has length at least m, the first n entries are the partial RNS representation.

    // Compute t_i = (u_i * q_i) mod p_i in signed form
    // Here, u[i] is int32_t. We'll cast to 64-bit for arithmetic.
    int32_t t[2]; // for n=2 as in the original code snippet
    for (int i = 0; i < n; i++) {
        int64_t ui = (int64_t)u[i];
        int64_t pi = (int64_t)p[i];
        int64_t qi = (int64_t)q[i];

        // If ui < 0, add p[i] to make it nonnegative for modular multiplication
        //if (ui < 0) ui += pi;

        // Use signed_mod_mul to compute (ui * qi) mod p[i], in signed form
        int64_t tt = signed_mod_mul(ui, qi, pi);
        t[i] = (int32_t) tt;
    }

    // Compute alpha = sum(t_i / p_i)
    double alpha = 0.0;
    for (int i = 0; i < n; i++) {
        alpha += (double)t[i] / (double)p[i];
    }

    int64_t r = (int64_t)llround(alpha);

    // Compute u = sum(t_i * P_i) - P * r
    // Use __int128 for safe intermediate
    __int128 U = 0;
    for (int i = 0; i < n; i++) {
        __int128 part = (__int128)t[i] * (__int128)P_i[i];
        U += part;
    }
    U -= (__int128)P * (__int128)r;

    // Now we have U as the reconstructed integer
    // Convert to int64_t (safe if |U| < 2^63)
    int64_t u_reconstructed = (int64_t)U;

    // Compute full RNS representation
    for (int j = 0; j < m; j++) {
        u[j] = (int32_t)signed_mod(u_reconstructed, p[j]);
    }
}



int main() {
    // Example usage:
    // Suppose we know u is represented partially by u_1 and u_2:
    // Let's pick an example u and check:
    int64_t u_example = 123456789876543; // just an example

    // partial representation:
    int32_t u[4];
    for (int i = 0; i < n; i++) {
        // signed representative mod p[i]
        int64_t val = u_example % p[i];
        if (val < 0) val += p[i];
        if (val > (int64_t)p[i]/2) val -= (int64_t)p[i];
        u[i] = val;
    }

    printf("Partial representation: u_1 = %d, u_2 = %d\n", u[0], u[1]);

    ecrt_full_representation(u);

    // Verify correctness:
    // u_full_res[i] should be u_example mod p[i] in signed representation.
    for (int i = 0; i < m; i++) {
        int64_t expected = (int64_t)((u_example % p[i] + p[i]) % p[i]);
        if (expected > (int64_t)p[i]/2) expected -= (int64_t)p[i];
        if (u[i] != expected) {
            fprintf(stderr, "Mismatch at modulus %u: got %lld, expected %lld\n",
                    p[i], (long long)u[i], (long long)expected);
            return 1;
        }
    }

    printf("All residues match.\n");

    
    return 0;
}
