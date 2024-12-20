#include "butterfly.h"
#include <stdint.h>




int32_t mod_mul_standard(int32_t a, int32_t b, int32_t mod) {
    int32_t result = (int32_t)(((int64_t)a * b) % mod);
    if (result > (mod>>1)) {
        result -= mod;
    }
    if (result < -(mod>>1)) {
        result += mod;
    }
    return result;
    
}

int32_t Barrett_mul(int32_t a, int32_t b, int32_t mod, int32_t mu_b) {
    int32_t z = a * b;
    int32_t t_high = (int32_t)(((int64_t)a * (int64_t)(mu_b<<1)) >> 32);
    int32_t result = (z - t_high * mod) ;

    if (result > (mod>>1)) {
        result = result - mod;
    }
    if (result < -(mod>>1)) {
        result = result + mod;
    }
    return result;
}

int32_t mod_add(int32_t a, int32_t b, int32_t mod) {
    int32_t result = a + b;
    if (result > (mod>>1)) {
        return result - mod;
    }
    return result;
}

int32_t mod_sub(int32_t a, int32_t b, int32_t mod) {
    int32_t result = a - b + mod;
    if (result > (mod>>1)) {
        return result - mod;
    }
    return result;
}



void CT_4pts(int32_t P, 
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t omega_1, int32_t omega2, int32_t omega3,
             int32_t mu_omega_1, int32_t mu_omega2, int32_t mu_omega3) {
    
    int32_t temp0;



    //printf("t0: %d, t1: %d, t2: %d, t3: %d\n", *t0, *t1, *t2, *t3);

    // first layer
    *t2 = Barrett_mul(*t2, omega_1, P, mu_omega_1);
    *t3 = Barrett_mul(*t3, omega_1, P, mu_omega_1);



    temp0 = mod_add(*t0, *t2, P);
    *t2 = mod_sub(*t0, *t2, P);
    *t0 = temp0;

    temp0 = mod_add(*t1, *t3, P);
    *t3 = mod_sub(*t1, *t3, P);
    *t1 = temp0;


    // second layer
    *t1 = Barrett_mul(*t1, omega2, P, mu_omega2);
    *t3 = Barrett_mul(*t3, omega3, P, mu_omega3);



    temp0 = mod_add(*t0, *t1, P);
    *t1 = mod_sub(*t0, *t1, P);
    *t0 = temp0;

    temp0 = mod_add(*t2, *t3, P);
    *t3 = mod_sub(*t2, *t3, P);
    *t2 = temp0;
}

void GS_4pts(int32_t P, 
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t omega_1, int32_t omega2, int32_t omega3,
             int32_t mu_omega_1, int32_t mu_omega2, int32_t mu_omega3) {

    int32_t temp0;
    int32_t temp1;

    temp0 = mod_add(*t0, *t1, P);
    *t1 = mod_sub(*t0, *t1, P);
    *t0 = temp0;

    temp1 = mod_add(*t2, *t3, P);
    *t3 = mod_sub(*t2, *t3, P);
    *t2 = temp1;

    *t1 = Barrett_mul(*t1, omega2, P, mu_omega2);
    *t3 = Barrett_mul(*t3, omega3, P, mu_omega3);


    temp0 = mod_add(*t0, *t2, P);
    *t2 = mod_sub(*t0, *t2, P);
    *t0 = temp0;

    temp1 = mod_add(*t1, *t3, P);
    *t3 = mod_sub(*t1, *t3, P);
    *t1 = temp1;

    *t2 = Barrett_mul(*t2, omega_1, P, mu_omega_1);
    *t3 = Barrett_mul(*t3, omega_1, P, mu_omega_1);

    }

#include <stdint.h>

// Assume mod_add, mod_sub, Barrett_mul, etc., are the same as in your code.

// CT_8pts performs one 8-point butterfly.
// Inputs/outputs: t0, t1, t2, t3, t4, t5, t6, t7
// P: modulus
// Twiddle factors: You must provide appropriate 8th-root-of-unity powers.
// For example, for each layer you’ll have a set of twiddles and their Barrett constants.
//
// Naming convention suggestion (you can adapt this):
// omega_1: W^(N/2)   (e.g. W^(4) for an 8-point NTT, used in first layer)
// omega_2, omega_3: W^(N/4), W^(3N/4) or similar for the second layer
// omega_4, omega_5, ... for the final layer.
// The exact mapping depends on your chosen NTT factorization and how you index twiddles.
//
// Here we’ll just use omega_1, omega_2, omega_3, omega_4, omega_5 for demonstration.
// You may need to verify which twiddles correspond to which pairs.

void CT_8pts(int32_t P,
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t *t4, int32_t *t5, int32_t *t6, int32_t *t7,
             int32_t omega_1, int32_t mu_omega_1,
             int32_t omega_2, int32_t mu_omega_2,
             int32_t omega_3, int32_t mu_omega_3,
             int32_t omega_4, int32_t mu_omega_4,
             int32_t omega_5, int32_t mu_omega_5,
             int32_t omega_6, int32_t mu_omega_6,
             int32_t omega_7, int32_t mu_omega_7) {

    int32_t temp0;

    // ---- Layer 1 (butterfly distance = 4) ----
    // Multiply second half (t4, t5, t6, t7) by omega_1
    *t4 = Barrett_mul(*t4, omega_1, P, mu_omega_1);
    *t5 = Barrett_mul(*t5, omega_1, P, mu_omega_1);
    *t6 = Barrett_mul(*t6, omega_1, P, mu_omega_1);
    *t7 = Barrett_mul(*t7, omega_1, P, mu_omega_1);

    // Combine (t0,t4), (t1,t5), (t2,t6), (t3,t7)
    temp0 = mod_add(*t0, *t4, P);
    *t4 = mod_sub(*t0, *t4, P);
    *t0 = temp0;

    temp0 = mod_add(*t1, *t5, P);
    *t5 = mod_sub(*t1, *t5, P);
    *t1 = temp0;

    temp0 = mod_add(*t2, *t6, P);
    *t6 = mod_sub(*t2, *t6, P);
    *t2 = temp0;

    temp0 = mod_add(*t3, *t7, P);
    *t7 = mod_sub(*t3, *t7, P);
    *t3 = temp0;


    // ---- Layer 2 (butterfly distance = 2) ----
    // Now we do pairs: (t0,t2), (t1,t3), (t4,t6), (t5,t7)
    // Assign appropriate twiddles for these pairs’ second elements:
    *t2 = Barrett_mul(*t2, omega_2, P, mu_omega_2); // for example
    *t3 = Barrett_mul(*t3, omega_2, P, mu_omega_2);
    *t6 = Barrett_mul(*t6, omega_3, P, mu_omega_3);
    *t7 = Barrett_mul(*t7, omega_3, P, mu_omega_3);

    // Combine
    temp0 = mod_add(*t0, *t2, P);
    *t2 = mod_sub(*t0, *t2, P);
    *t0 = temp0;

    temp0 = mod_add(*t1, *t3, P);
    *t3 = mod_sub(*t1, *t3, P);
    *t1 = temp0;

    temp0 = mod_add(*t4, *t6, P);
    *t6 = mod_sub(*t4, *t6, P);
    *t4 = temp0;

    temp0 = mod_add(*t5, *t7, P);
    *t7 = mod_sub(*t5, *t7, P);
    *t5 = temp0;


    // ---- Layer 3 (butterfly distance = 1) ----
    // Pairs: (t0,t1), (t2,t3), (t4,t5), (t6,t7)
    // Multiply second half of each pair by appropriate twiddles:
    *t1 = Barrett_mul(*t1, omega_4, P, mu_omega_4);
    *t3 = Barrett_mul(*t3, omega_5, P, mu_omega_5);
    *t5 = Barrett_mul(*t5, omega_6, P, mu_omega_6);
    *t7 = Barrett_mul(*t7, omega_7, P, mu_omega_7);

    // Combine
    temp0 = mod_add(*t0, *t1, P);
    *t1 = mod_sub(*t0, *t1, P);
    *t0 = temp0;

    temp0 = mod_add(*t2, *t3, P);
    *t3 = mod_sub(*t2, *t3, P);
    *t2 = temp0;

    temp0 = mod_add(*t4, *t5, P);
    *t5 = mod_sub(*t4, *t5, P);
    *t4 = temp0;

    temp0 = mod_add(*t6, *t7, P);
    *t7 = mod_sub(*t6, *t7, P);
    *t6 = temp0;
}


void GS_8pts(int32_t P,
             int32_t *t0, int32_t *t1, int32_t *t2, int32_t *t3,
             int32_t *t4, int32_t *t5, int32_t *t6, int32_t *t7,
             int32_t omega_1, int32_t mu_omega_1,
             int32_t omega_2, int32_t mu_omega_2,
             int32_t omega_3, int32_t mu_omega_3,
             int32_t omega_4, int32_t mu_omega_4,
             int32_t omega_5, int32_t mu_omega_5,
             int32_t omega_6, int32_t mu_omega_6,
             int32_t omega_7, int32_t mu_omega_7) {

    int32_t temp0;

    // ---- Layer 3 (distance = 1) ----
    // Combine pairs: (t0,t1), (t2,t3), (t4,t5), (t6,t7)
    temp0 = mod_add(*t0, *t1, P);
    *t1 = mod_sub(*t0, *t1, P);
    *t0 = temp0;

    temp0 = mod_add(*t2, *t3, P);
    *t3 = mod_sub(*t2, *t3, P);
    *t2 = temp0;

    temp0 = mod_add(*t4, *t5, P);
    *t5 = mod_sub(*t4, *t5, P);
    *t4 = temp0;

    temp0 = mod_add(*t6, *t7, P);
    *t7 = mod_sub(*t6, *t7, P);
    *t6 = temp0;

    // Multiply the differences by the appropriate twiddles:
    // (These twiddles correspond to the third layer of CT, which in GS are applied here.)
    *t1 = Barrett_mul(*t1, omega_4, P, mu_omega_4);
    *t3 = Barrett_mul(*t3, omega_5, P, mu_omega_5);
    *t5 = Barrett_mul(*t5, omega_6, P, mu_omega_6);
    *t7 = Barrett_mul(*t7, omega_7, P, mu_omega_7);


    // ---- Layer 2 (distance = 2) ----
    // Combine pairs: (t0,t2), (t1,t3), (t4,t6), (t5,t7)
    temp0 = mod_add(*t0, *t2, P);
    *t2 = mod_sub(*t0, *t2, P);
    *t0 = temp0;

    temp0 = mod_add(*t1, *t3, P);
    *t3 = mod_sub(*t1, *t3, P);
    *t1 = temp0;

    temp0 = mod_add(*t4, *t6, P);
    *t6 = mod_sub(*t4, *t6, P);
    *t4 = temp0;

    temp0 = mod_add(*t5, *t7, P);
    *t7 = mod_sub(*t5, *t7, P);
    *t5 = temp0;

    // Multiply the differences by the appropriate twiddles:
    // (Matches what was done in CT layer 2, but done after combining in GS.)
    *t2 = Barrett_mul(*t2, omega_2, P, mu_omega_2);
    *t3 = Barrett_mul(*t3, omega_2, P, mu_omega_2);
    *t6 = Barrett_mul(*t6, omega_3, P, mu_omega_3);
    *t7 = Barrett_mul(*t7, omega_3, P, mu_omega_3);


    // ---- Layer 1 (distance = 4) ----
    // Combine pairs: (t0,t4), (t1,t5), (t2,t6), (t3,t7)
    temp0 = mod_add(*t0, *t4, P);
    *t4 = mod_sub(*t0, *t4, P);
    *t0 = temp0;

    temp0 = mod_add(*t1, *t5, P);
    *t5 = mod_sub(*t1, *t5, P);
    *t1 = temp0;

    temp0 = mod_add(*t2, *t6, P);
    *t6 = mod_sub(*t2, *t6, P);
    *t2 = temp0;

    temp0 = mod_add(*t3, *t7, P);
    *t7 = mod_sub(*t3, *t7, P);
    *t3 = temp0;

    // Multiply the differences by the appropriate twiddle (omega_1 for the largest distance)
    *t4 = Barrett_mul(*t4, omega_1, P, mu_omega_1);
    *t5 = Barrett_mul(*t5, omega_1, P, mu_omega_1);
    *t6 = Barrett_mul(*t6, omega_1, P, mu_omega_1);
    *t7 = Barrett_mul(*t7, omega_1, P, mu_omega_1);

}
