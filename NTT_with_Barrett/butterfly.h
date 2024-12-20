#ifndef BUTTERFLY_H
#define BUTTERFLY_H

#include <stdint.h>

// Function prototypes for modular arithmetic
int32_t mod_mul_standard(int32_t a, int32_t b, int32_t mod);
int32_t Barrett_mul(int32_t a, int32_t b, int32_t mod, int32_t mu_b);
int32_t mod_add(int32_t a, int32_t b, int32_t mod);
int32_t mod_sub(int32_t a, int32_t b, int32_t mod);

// Function prototypes for butterfly operations
void CT_4pts(int32_t P, 
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t omega_1, int32_t omega2, int32_t omega3,
             int32_t mu_omega_1, int32_t mu_omega2, int32_t mu_omega3);

void GS_4pts(int32_t P, 
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t omega_1, int32_t omega2, int32_t omega3,
             int32_t mu_omega_1, int32_t mu_omega2, int32_t mu_omega3);

#endif // BUTTERFLY_H
