#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "data.h" // Contains n, m, p[], q[], etc.

// Compute signed representative modulo p: a value in [-p/2, p/2)
static int64_t signed_mod(int64_t x, uint32_t m) {
    int64_t val = (int64_t)(x % m);
    if (val < 0) val += m; // ensure nonnegative
    if (val > (int64_t)m/2) val -= (int64_t)m;
    return val;
}

// signed_mod_mul: Compute (a*b) mod m, returning signed representative in [-m/2, m/2)
static int64_t signed_mod_mul(int64_t a, int64_t b, uint32_t m) {
    int64_t am = a % m;
    if (am < 0) am += m;
    int64_t bm = b % m;
    if (bm < 0) bm += m;

    // 64-bit multiplication is safe: (am, bm < 2^31)
    __int128 prod = (__int128)am * (__int128)bm;
    __int128 res = prod % m;
    int64_t r = (int64_t)res;

    if (r > (int64_t)m/2) r -= m;
    return r;
}

// Given partial representation (u_1,...,u_n) in u,
// produce full representation (u_1,...,u_m) in u.
// We have precomputed q[], P_i_mod[], P_mod[].
static void ecrt_full_representation(int32_t *u) {
    // Compute t_i = (u_i * q_i) mod p_i (signed)
    // n=4 in this example
    int32_t t[4];
    for (int i = 0; i < n; i++) {
        int64_t ui = (int64_t)u[i];
        uint32_t pi = p[i];
        uint32_t qi = q[i];
        int64_t tt64 = signed_mod_mul(ui, qi, pi);
        t[i] = (int32_t)tt64;
    }

    // Compute alpha = sum(t_i / p_i)
    double alpha = 0.0;
    for (int i = 0; i < n; i++) {
        alpha += (double)t[i] / (double)p[i];
    }

    int64_t r = (int64_t)llround(alpha);

    // Now compute u mod each p_j:
    // u mod p_j = [Σ t_i (P_i mod p_j) - (P mod p_j)*r] mod p_j, in signed form
    // Retrieve P_i_mod[i][j] and P_mod[j] from data.c arrays
    // P_i_mod is stored row-wise: P_i_mod[i*m + j]
    for (int j = 0; j < m; j++) {
        uint32_t pj = p[j];
        // sum_tiPi = Σ t_i*(P_i mod p_j)
        // Handle i < n, P_i only for the first n moduli
        // If j < n, (P_i mod p_j) from the same array
        int64_t sum_tiPi = 0;
        for (int i = 0; i < n; i++) {
            uint32_t Pi_mod_pj = P_i_mod[i*m + j];
            // signed multiplication mod pj:
            // careful: t[i] fits in 32 bits, Pi_mod_pj < pj < 2^31, no overflow in 64 bits.
            int64_t part = (int64_t)t[i] * (int64_t)Pi_mod_pj;
            int64_t val = part % pj;
            if (val < 0) val += pj;
            sum_tiPi += val;
        }
        int64_t Pmod_pj = P_mod[j];
        // sum_tiPi mod pj
        sum_tiPi = sum_tiPi % pj;

        // subtract (P mod p_j)*r mod p_j
        // (Pmod_pj * r) % pj
        int64_t Pr = ((int64_t)Pmod_pj * r) % pj;
        if (Pr < 0) Pr += pj;

        int64_t u_mod_pj = sum_tiPi - Pr;
        u_mod_pj %= pj;
        if (u_mod_pj < 0) u_mod_pj += pj;

        // convert to signed representative
        if (u_mod_pj > (int64_t)pj/2) u_mod_pj -= pj;

        u[j] = (int32_t)u_mod_pj;
    }
}

// A function to generate a random residue in the range [-p/2, p/2)
static int32_t random_signed_residue(uint32_t p_i) {
    uint32_t r = (uint32_t)rand();
    uint32_t mod_val = r % p_i;
    int64_t val = (int64_t)mod_val;
    if (val > (int64_t)p_i/2) val -= (int64_t)p_i;
    return (int32_t)val;
}

int main() {
    srand((unsigned)time(NULL));

    // For demonstration, just do one test:
    // You can loop over multiple tests if desired.
    int TESTS = 1;
    for (int test = 1; test <= TESTS; test++) {
        int32_t u_array[m];
        int32_t original_partial[n];

        // Generate random partial representation
        // We know that these partial residues define some integer u mod P,
        // but we do not explicitly know u. We only check consistency.
        for (int i = 0; i < n; i++) {
            int32_t residue = random_signed_residue(p[i]);
            u_array[i] = residue;
            original_partial[i] = residue;
        }

        // At this point, u_array[0..n-1] is our random partial representation.

        // Print partial representation
        printf("Partial:");
        for (int i = 0; i < n; i++) {
            printf(" %d", u_array[i]);
        }
        printf("\n");
        
        // Compute the full RNS representation
        ecrt_full_representation(u_array);


        // Print full representation
        printf("Full:");
        for (int i = 0; i < m; i++) {
            printf(" %d", u_array[i]);
        }
        printf("\n");

        
        // Verify that the partial representation didn't change.
        // If the ECRT is consistent, these should remain correct.
        for (int i = 0; i < n; i++) {
            if (u_array[i] != original_partial[i]) {
                fprintf(stderr, "Test %d failed: mismatch at index %d\n", test, i);
                fprintf(stderr, "Expected %d, got %d\n", original_partial[i], u_array[i]);
                return 1;
            }
        }

        printf("Test %d passed: partial representation stable.\n", test);

        // Optional: If you implement a CRT here using p[0..n-1] and u_array[0..n-1],
        // you could reconstruct u (mod P) and then verify full representation.
        // But that requires additional big-integer and CRT logic.
    }

    printf("All tests passed.\n");
    return 0;
}