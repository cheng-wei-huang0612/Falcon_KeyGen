#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "data.h"

// Compute signed representative modulo p: a value in [-p/2, p/2)
static int64_t signed_mod(int64_t x, uint32_t m) {
    int64_t val = (int64_t)(x % m);
    if (val < 0) val += m; 
    if (val > (int64_t)m/2) val -= (int64_t)m;
    return val;
}

// signed_mod_mul: Compute (a*b) mod m, returning signed representative
static int64_t signed_mod_mul(int64_t a, int64_t b, uint32_t m) {
    int64_t am = a % m;
    if (am < 0) am += m;
    int64_t bm = b % m;
    if (bm < 0) bm += m;

    __int128 prod = (__int128)am * (__int128)bm;
    __int128 res = prod % m;
    int64_t r = (int64_t)res;

    if (r > (int64_t)m/2) r -= m;
    return r;
}

// ECRT for n=216, m=304
// Input: u[0..n-1] partial representation
// Output: u[0..m-1] full representation
static void ecrt_full_representation(int32_t *u) {
    int32_t t[216];

    // Compute t_i = (u_i * q_i) mod p_i in signed form
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

    // Reconstruct u mod each p_j
    // u mod p_j = (Σ t_i * (P_i mod p_j) - (P mod p_j)*r) mod p_j
    for (int j = 0; j < m; j++) {
        uint32_t pj = p[j];
        int64_t sum_tiPi = 0;
        for (int i = 0; i < n; i++) {
            // Pi_mod_pj = P_i_mod[i*m + j];
            // t[i] fits in 32-bit, Pi_mod_pj < pj
            int32_t val = signed_mod_mul(t[i], P_i_mod[i*m + j], pj);
            sum_tiPi += val;
        }

        sum_tiPi = signed_mod(sum_tiPi, pj);

        //int64_t Pmod_pj = P_mod[j];
        int64_t Pr = signed_mod_mul(P_mod[j], r, pj);

        //int64_t u_mod_pj = sum_tiPi - Pr;
        //u_mod_pj %= pj;
        //if (u_mod_pj < 0) u_mod_pj += pj;
        int64_t u_mod_pj = signed_mod(sum_tiPi - Pr, pj);


        u[j] = (int32_t)u_mod_pj;
    }
}

// Function to generate a random signed residue in [-p_i/2, p_i/2].
static int32_t random_signed_residue(uint32_t p_i) {
    // Generate a uniform random integer in [0, p_i-1]
    // Then shift to signed range.
    uint32_t r = (uint32_t)rand();
    // Use modulo to get a number in [0, p_i-1]
    // Though not perfectly uniform for large p_i, good enough for a demo
    uint32_t mod_val = r % p_i;

    // Convert to signed: if mod_val > p_i/2, subtract p_i
    int64_t val = (int64_t)mod_val;
    if (val > (int64_t)p_i/2) val -= (int64_t)p_i;
    return (int32_t)val;
}

int main() {
    srand((unsigned)time(NULL));

    // Number of tests
    int TESTS = 1; // You can increase this count as needed.

    for (int test = 1; test <= TESTS; test++) {
        // Generate a random partial representation:
        int32_t u_array[m];
        int32_t original_partial[n]; // Store original partial residues for verification

        for (int i = 0; i < n; i++) {
            int32_t residue = random_signed_residue(p[i]);
            u_array[i] = residue;
            original_partial[i] = residue;
        }

        // Print partial representation
        printf("Partial:");
        for (int i = 0; i < n; i++) {
            printf(" %d", u_array[i]);
        }
        printf("\n");

        // Compute full RNS
        ecrt_full_representation(u_array);

        // Print full representation
        printf("Full:");
        for (int i = 0; i < m; i++) {
            printf(" %d", u_array[i]);
        }
        printf("\n");

        // Verify that the first n residues match what we started with
        for (int i = 0; i < n; i++) {
            if (u_array[i] != original_partial[i]) {
                fprintf(stderr, "Test %d failed: mismatch at index %d\n", test, i);
                fprintf(stderr, "Expected %d, got %d\n", original_partial[i], u_array[i]);
                return 1;
            }
        }

        printf("Test %d passed.\n", test);
    }

    printf("All tests passed.\n");
    return 0;
}
