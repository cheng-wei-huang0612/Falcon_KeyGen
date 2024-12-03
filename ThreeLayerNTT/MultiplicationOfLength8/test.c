// Goal: Test the three layer NTT algorithm on length 8 poly mul
// With Barrett, NEON 


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "poly_mul.h"
#include "data.h"





void school_book(int32_t *p, int32_t *q, int32_t *r) {
    for (int i = 0; i < 2*N; i++) {
        r[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            r[i+j] += p[i] * q[j];
            r[i+j] %= 12289;
        }
    }

    for (int i = 0; i < N; i++) {
        r[i] = (r[i] - r[i+N]) % 12289;
    }
}

int main() {

    printf("Polynomial f: ");
    for (int i = 0; i < N; i++) {
        printf("%d,", f[i]);
    }

    printf("\n\nPolynomial g: ");
    for (int i = 0; i < N; i++) {
        printf("%d,", g[i]);
    }
    

    




    return 0;
}
