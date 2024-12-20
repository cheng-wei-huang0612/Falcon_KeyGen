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
    int32_t temp1;


    //printf("t0: %d, t1: %d, t2: %d, t3: %d\n", *t0, *t1, *t2, *t3);

    // first layer
    *t2 = Barrett_mul(*t2, omega_1, P, mu_omega_1);
    *t3 = Barrett_mul(*t3, omega_1, P, mu_omega_1);



    temp0 = mod_add(*t0, *t2, P);
    *t2 = mod_sub(*t0, *t2, P);
    *t0 = temp0;

    temp1 = mod_add(*t1, *t3, P);
    *t3 = mod_sub(*t1, *t3, P);
    *t1 = temp1;


    // second layer
    *t1 = Barrett_mul(*t1, omega2, P, mu_omega2);
    *t3 = Barrett_mul(*t3, omega3, P, mu_omega3);



    temp0 = mod_add(*t0, *t1, P);
    *t1 = mod_sub(*t0, *t1, P);
    *t0 = temp0;

    temp1 = mod_add(*t2, *t3, P);
    *t3 = mod_sub(*t2, *t3, P);
    *t2 = temp1;
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
