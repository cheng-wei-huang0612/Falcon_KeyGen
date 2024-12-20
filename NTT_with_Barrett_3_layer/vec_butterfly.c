// vec_butterfly.c

#include "vec_butterfly.h"
#include <arm_neon.h>
#include <stdint.h>
#include <stdio.h>

// Vectorized Cooley-Turkey butterfly for four parallel FFT/NTT instances
void CT_4pts_vec(
    const int32_t P[4],
    int32_t t0[4], int32_t t1[4], int32_t t2[4], int32_t t3[4],
    const int32_t omega_1[4], const int32_t omega_2[4], const int32_t omega_3[4],
    const int32_t mu_omega_1[4], const int32_t mu_omega_2[4], const int32_t mu_omega_3[4])
{
    // Load data into NEON vectors
    int32x4_t vT0 = vld1q_s32(t0);
    int32x4_t vT1 = vld1q_s32(t1);
    int32x4_t vT2 = vld1q_s32(t2);
    int32x4_t vT3 = vld1q_s32(t3);

    // Load moduli and compute half_mod = P >> 1
    int32x4_t vP = vld1q_s32(P);
    int32x4_t vHalf = vshrq_n_s32(vP, 1);

    // Load omega and mu_omega constants
    int32x4_t vOmega1    = vld1q_s32(omega_1);
    int32x4_t vOmega2    = vld1q_s32(omega_2);
    int32x4_t vOmega3    = vld1q_s32(omega_3);
    int32x4_t vMu_Omega1 = vld1q_s32(mu_omega_1);
    int32x4_t vMu_Omega2 = vld1q_s32(mu_omega_2);
    int32x4_t vMu_Omega3 = vld1q_s32(mu_omega_3);




    // First layer: Multiply t2 and t3 by omega_1
    vT2 = vBarrett_mul_const(vT2, vOmega1, vP, vMu_Omega1);
    vT3 = vBarrett_mul_const(vT3, vOmega1, vP, vMu_Omega1);


    int32x4_t temp0 = vmod_add(vT0, vT2, vP, vHalf);
    vT2 = vmod_sub(vT0, vT2, vP, vHalf);
    vT0 = temp0;

    int32x4_t temp1 = vmod_add(vT1, vT3, vP, vHalf);
    vT3 = vmod_sub(vT1, vT3, vP, vHalf);
    vT1 = temp1;

    // Second layer: Multiply t1 by omega2 and t3 by omega3
    vT1 = vBarrett_mul_const(vT1, vOmega2, vP, vMu_Omega2);
    vT3 = vBarrett_mul_const(vT3, vOmega3, vP, vMu_Omega3);

    temp0 = vmod_add(vT0, vT1, vP, vHalf);
    vT1 = vmod_sub(vT0, vT1, vP, vHalf);
    vT0 = temp0;

    temp1 = vmod_add(vT2, vT3, vP, vHalf);
    vT3 = vmod_sub(vT2, vT3, vP, vHalf);
    vT2 = temp1;

    // Store results back to memory
    vst1q_s32(t0, vT0);
    vst1q_s32(t1, vT1);
    vst1q_s32(t2, vT2);
    vst1q_s32(t3, vT3);
}

// Vectorized Gentleman-Sande butterfly for four parallel FFT/NTT instances
void GS_4pts_vec(
    const int32_t P[4],
    int32_t t0[4], int32_t t1[4], int32_t t2[4], int32_t t3[4],
    const int32_t omega_1[4], const int32_t omega2[4], const int32_t omega_3[4],
    const int32_t mu_omega_1[4], const int32_t mu_omega_2[4], const int32_t mu_omega_3[4])
{
    // Load data into NEON vectors
    int32x4_t vT0 = vld1q_s32(t0);
    int32x4_t vT1 = vld1q_s32(t1);
    int32x4_t vT2 = vld1q_s32(t2);
    int32x4_t vT3 = vld1q_s32(t3);

    // Load moduli and compute half_mod = P >> 1
    int32x4_t vP = vld1q_s32(P);
    int32x4_t vHalf = vshrq_n_s32(vP, 1);

    // Load omega and mu_omega constants
    int32x4_t vOmega1    = vld1q_s32(omega_1);
    int32x4_t vOmega2    = vld1q_s32(omega2);
    int32x4_t vOmega3    = vld1q_s32(omega_3);
    int32x4_t vMu_Omega1 = vld1q_s32(mu_omega_1);
    int32x4_t vMu_Omega2 = vld1q_s32(mu_omega_2);
    int32x4_t vMu_Omega3 = vld1q_s32(mu_omega_3);

    int32x4_t temp0 = vmod_add(vT0, vT1, vP, vHalf);
    vT1 = vmod_sub(vT0, vT1, vP, vHalf);
    vT0 = temp0;

    int32x4_t temp1 = vmod_add(vT2, vT3, vP, vHalf);
    vT3 = vmod_sub(vT2, vT3, vP, vHalf);
    vT2 = temp1;

    // Multiply t1 by omega2 and t3 by omega3
    vT1 = vBarrett_mul_const(vT1, vOmega2, vP, vMu_Omega2);
    vT3 = vBarrett_mul_const(vT3, vOmega3, vP, vMu_Omega3);

    temp0 = vmod_add(vT0, vT2, vP, vHalf);
    vT2 = vmod_sub(vT0, vT2, vP, vHalf);
    vT0 = temp0;

    temp1 = vmod_add(vT1, vT3, vP, vHalf);
    vT3 = vmod_sub(vT1, vT3, vP, vHalf);
    vT1 = temp1;

    // Multiply t2 and t3 by omega_1
    vT2 = vBarrett_mul_const(vT2, vOmega1, vP, vMu_Omega1);
    vT3 = vBarrett_mul_const(vT3, vOmega1, vP, vMu_Omega1);

    // Store results back to memory
    vst1q_s32(t0, vT0);
    vst1q_s32(t1, vT1);
    vst1q_s32(t2, vT2);
    vst1q_s32(t3, vT3);
}
