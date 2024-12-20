#include <stdio.h>
#include <stdint.h>
#include <arm_neon.h>
#include <assert.h>
#include "data.h"
#include "utils.h"
#include "vec_butterfly.h"
#include "butterfly.h"

// We include the static inline functions for testing directly if allowed.
// Otherwise, you can move them to a separate header or make them non-static for testing.

static void test_vmod_adjust() {
    // Test with a simple P and some values.
    // P = 17, half_mod = P>>1 = 8. Range: [-8,8]
    int32_t P_arr[4]       = {17,17,17,17};
    int32_t half_mod_arr[4]= {8,8,8,8};
    int32x4_t vP       = vld1q_s32(P_arr);
    int32x4_t vHalfMod = vld1q_s32(half_mod_arr);

    // Test values
    int32_t x_arr[4] = {9, -9, 8, -8}; 
    // 9 > 8 => 9-17=-8
    // -9 < -8 => -9+17=8
    // 8 and -8 are on boundary, should remain the same
    int32x4_t x = vld1q_s32(x_arr);

    int32x4_t res = vmod_adjust(x, vHalfMod, vP);

    int32_t out[4];
    vst1q_s32(out, res);

    // Expected: { -8, 8, 8, -8 }
    assert(out[0] == -8);
    assert(out[1] == 8);
    assert(out[2] == 8);
    assert(out[3] == -8);
}

static void test_vmod_add() {
    // Same P = 17
    int32_t P_arr[4] = {17,17,17,17};
    int32x4_t vP     = vld1q_s32(P_arr);
    int32x4_t vHalf  = vshrq_n_s32(vP, 1); // half = 8

    int32_t a_arr[4] = {1, -8, 7, 8};
    int32_t b_arr[4] = {10, 9, -7, -8};
    int32x4_t a = vld1q_s32(a_arr);
    int32x4_t b = vld1q_s32(b_arr);

    // a+b mod 17, adjusted to [-8,8]
    // 1+10=11 in [-8,8]? 11>8 => 11-17=-6
    // -8+9=1 fits in [-8,8]
    // 7+(-7)=0 fits
    // 8+(-8)=0 fits
    int32x4_t res = vmod_add(a, b, vP, vHalf);
    int32_t out[4]; vst1q_s32(out, res);

    // Expected: { -6, 1, 0, 0 }
    assert(out[0] == -6);
    assert(out[1] == 1);
    assert(out[2] == 0);
    assert(out[3] == 0);
}

static void test_vmod_sub() {
    // P=17
    int32_t P_arr[4] = {17,17,17,17};
    int32x4_t vP     = vld1q_s32(P_arr);
    int32x4_t vHalf  = vshrq_n_s32(vP, 1); // 8

    int32_t a_arr[4] = {1, -8, 7, 8};
    int32_t b_arr[4] = {10, 9, -7, -8};
    int32x4_t a = vld1q_s32(a_arr);
    int32x4_t b = vld1q_s32(b_arr);

    // a-b+P mod, adjusted to [-8,8]
    // (1-10+17)=8 fits
    // (-8-9+17)=0 fits
    // (7-(-7)+17)=7+7+17=31,31>8 =>31-17=14>8 =>14-17=-3
    // (8-(-8)+17)=8+8+17=33 =>33-17=16>8 =>16-17=-1
    int32x4_t res = vmod_sub(a, b, vP, vHalf);
    int32_t out[4]; vst1q_s32(out, res);

    // Expected: {8, 0, -3, -1}
    assert(out[0] == 8);
    assert(out[1] == 0);
    assert(out[2] == -3);
    assert(out[3] == -1);
}

static void test_vBarrett_mul_const() {
    // For a simple test, let's pick:
    // N=17, mu_c representing (bR/N)/2 in Q1.31 is tricky.
    // We'll just do a scenario where mu_c = 0 (no reduction) and see if we get a*c result as expected mod N.
    // Without the final adjustment step in vBarrett_mul_const, we must rely on conditions given.
    // Here we just test the arithmetic: z = a*c - t*N, where t= sqrdmulh(a,mu_c).

    // Let's pick mu_c=0 for a trivial test: then t=0, z = a*c.
    // With mu_c=0: vqrdmulhq_s32(a,0)=0, so z=z-t*N=z.
    // This just tests if vmulq_s32 works as expected in isolation.
    int32_t N_arr[4] = {17,17,17,17};
    int32x4_t vN = vld1q_s32(N_arr);

    int32_t mu_arr[4] = {0,0,0,0};
    int32x4_t vMu = vld1q_s32(mu_arr);

    int32_t a_arr[4] = {1,2,3,4};
    int32_t c_arr[4] = {10,10,10,10};
    int32x4_t a = vld1q_s32(a_arr);
    int32x4_t c = vld1q_s32(c_arr);

    int32x4_t z = vBarrett_mul_const(a, c, vN, vMu);
    int32_t out[4]; vst1q_s32(out, z);

    // Since mu=0, no reduction, z = a*c:
    // {1*10=10, 2*10=20, 3*10=30, 4*10=40}
    assert(out[0] == 10);
    assert(out[1] == 20);
    assert(out[2] == 30);
    assert(out[3] == 40);
}


int main() {
    printf("Running unit tests for vec_butterfly helpers...\n");
    test_vmod_adjust();
    test_vmod_add();
    test_vmod_sub();
    test_vBarrett_mul_const();
    printf("All unit tests passed!\n\n");

    printf("Testing vec_butterfly...\n");

    int32_t t0[4] = {1, 10, 100, 1000};
    int32_t t1[4] = {2, 20, 200, 2000};
    int32_t t2[4] = {3, 30, 300, 3000};
    int32_t t3[4] = {4, 40, 400, 4000};

    // Use a consistent set of omega and mu_omega for both forward and inverse
    // For a single 4-pt transform, let's pick power_index = 1 for simplicity.
    int32_t omega_1[4], omega_2[4], omega_3[4];
    int32_t mu_omega_1[4], mu_omega_2[4], mu_omega_3[4];
    
    int p = 1;
    for (int i = 0; i < 4; i++) {
        omega_1[i]    = omegas[1024*i + p];
        omega_2[i]    = omegas[1024*i + (p << 1)];
        omega_3[i]    = omegas[1024*i + (p << 1) + 1];

        mu_omega_1[i] = mu_omegas[1024*i + p];
        mu_omega_2[i] = mu_omegas[1024*i + (p << 1)];
        mu_omega_3[i] = mu_omegas[1024*i + (p << 1) + 1];
    }

    printf("Before CT_4pts_vec:\n");
    PRINT_SPEC(t0, 0, 4, 1);
    PRINT_SPEC(t1, 0, 4, 1);
    PRINT_SPEC(t2, 0, 4, 1);
    PRINT_SPEC(t3, 0, 4, 1);

    // Forward transform
    CT_4pts_vec(P, t0, t1, t2, t3, omega_1, omega_2, omega_3, mu_omega_1, mu_omega_2, mu_omega_3);

    printf("After CT_4pts_vec:\n");
    PRINT_SPEC(t0, 0, 4, 1);
    PRINT_SPEC(t1, 0, 4, 1);
    PRINT_SPEC(t2, 0, 4, 1);
    PRINT_SPEC(t3, 0, 4, 1);

    // Now for the inverse transform, we supply the inverse omegas and mu_omegas.
    for (int i = 0; i < 4; i++) {
        omega_1[i]    = -omegas[1024*i + p];
        omega_2[i]    = -omegas[1024*i + (p << 1) + 1];
        omega_3[i]    = -omegas[1024*i + (p << 1)];

        mu_omega_1[i] = -(mu_omegas[1024*i + p] - 2);
        mu_omega_2[i] = -(mu_omegas[1024*i + (p << 1) + 1] - 2);
        mu_omega_3[i] = -(mu_omegas[1024*i + (p << 1)] - 2);
    }

    printf("Before GS_4pts_vec:\n");
    PRINT_SPEC(t0, 0, 4, 1);
    PRINT_SPEC(t1, 0, 4, 1);
    PRINT_SPEC(t2, 0, 4, 1);
    PRINT_SPEC(t3, 0, 4, 1);

    // Inverse transform
    GS_4pts_vec(P, t0, t1, t2, t3, omega_1, omega_2, omega_3, mu_omega_1, mu_omega_2, mu_omega_3);

    int32_t inv4[4] = {-8387584, -8385024, -8384512, -8383488};
    int32_t mu_inv4[4] = {-536870894, -536870894, -536870894, -536870894};

    for (int i = 0; i < 4; i++) {
        t0[i] = Barrett_mul(t0[i], inv4[0], P[0], mu_inv4[0]);
        t1[i] = Barrett_mul(t1[i], inv4[1], P[1], mu_inv4[1]);
        t2[i] = Barrett_mul(t2[i], inv4[2], P[2], mu_inv4[2]);
        t3[i] = Barrett_mul(t3[i], inv4[3], P[3], mu_inv4[3]);
    }

    printf("After GS_4pts_vec (and division by 4):\n");
    PRINT_SPEC(t0, 0, 4, 1);
    PRINT_SPEC(t1, 0, 4, 1);
    PRINT_SPEC(t2, 0, 4, 1);
    PRINT_SPEC(t3, 0, 4, 1);

    printf("Test completed.\n");
    return 0;
}
