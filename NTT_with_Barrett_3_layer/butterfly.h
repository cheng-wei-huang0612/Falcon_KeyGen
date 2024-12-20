#ifndef BUTTERFLY_H
#define BUTTERFLY_H

#include <stdint.h>

// Function prototypes for modular arithmetic
int32_t mod_mul_standard(int32_t a, int32_t b, int32_t mod);
int32_t Barrett_mul(int32_t a, int32_t b, int32_t mod, int32_t mu_b);
int32_t mod_add(int32_t a, int32_t b, int32_t mod);
int32_t mod_sub(int32_t a, int32_t b, int32_t mod);

// Function prototypes for butterfly operations

void CT_8pts(int32_t P,
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t *t4, int32_t *t5, int32_t *t6, int32_t *t7,
             int32_t omega_1, int32_t mu_omega_1,
             int32_t omega_2, int32_t mu_omega_2,
             int32_t omega_3, int32_t mu_omega_3,
             int32_t omega_4, int32_t mu_omega_4,
             int32_t omega_5, int32_t mu_omega_5,
             int32_t omega_6, int32_t mu_omega_6,
             int32_t omega_7, int32_t mu_omega_7);


void GS_8pts(int32_t P,
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t *t4, int32_t *t5, int32_t *t6, int32_t *t7,
             int32_t omega_1, int32_t mu_omega_1,
             int32_t omega_2, int32_t mu_omega_2,
             int32_t omega_3, int32_t mu_omega_3,
             int32_t omega_4, int32_t mu_omega_4,
             int32_t omega_5, int32_t mu_omega_5,
             int32_t omega_6, int32_t mu_omega_6,
             int32_t omega_7, int32_t mu_omega_7);  

void CT_4pts(int32_t P, 
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t omega_1, int32_t omega2, int32_t omega3,
             int32_t mu_omega_1, int32_t mu_omega2, int32_t mu_omega3);

void GS_4pts(int32_t P, 
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t omega_1, int32_t omega2, int32_t omega3,
             int32_t mu_omega_1, int32_t mu_omega2, int32_t mu_omega3);

#endif // BUTTERFLY_H
