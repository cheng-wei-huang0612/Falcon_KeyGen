#ifndef VEC_BUTTERFLY_H
#define VEC_BUTTERFLY_H

#include <stdint.h>
#include <arm_neon.h>

// Declare your vectorized functions here as static inline

// Adjust results to be within [-P/2, P/2]
static inline int32x4_t vmod_adjust(int32x4_t x, int32x4_t half_mod, int32x4_t mod) {
    // We loop until all lanes are in range
    while (1) {
        // Condition: x > half_mod
        uint32x4_t greater_mask = vcgtq_s32(x, half_mod);

        // Condition: x < -half_mod
        int32x4_t neg_half_mod = vnegq_s32(half_mod);
        uint32x4_t less_mask = vcltq_s32(x, neg_half_mod);

        // If no lane is greater and no lane is less, we are done
        // vcgtq_s32 and vcltq_s32 return a mask of all ones (0xFFFFFFFF) for true, zero for false.
        // If all false, then all lanes are in range.
        // Combine masks: if all lanes are in range, both masks should be zero.
        uint32x4_t combined = vorrq_u32(greater_mask, less_mask);
        // Check if combined is all zeros:
        if (vgetq_lane_u32(combined, 0) == 0 &&
            vgetq_lane_u32(combined, 1) == 0 &&
            vgetq_lane_u32(combined, 2) == 0 &&
            vgetq_lane_u32(combined, 3) == 0) {
            break; // all lanes are now in range
        }

        // Adjust lanes
        // Where result > half_mod, subtract mod
        int32x4_t x_sub_mod = vsubq_s32(x, mod);
        x = vbslq_s32(greater_mask, x_sub_mod, x);

        // Where result < -half_mod, add mod
        int32x4_t x_add_mod = vaddq_s32(x, mod);
        x = vbslq_s32(less_mask, x_add_mod, x);
    }

    return x;
}



// Modular addition using vmod_adjust
static inline int32x4_t vmod_add(int32x4_t a, int32x4_t b, int32x4_t P, int32x4_t half_mod) {
    int32x4_t result = vaddq_s32(a, b);
    return vmod_adjust(result, half_mod, P);
}

// Modular subtraction using vmod_adjust
static inline int32x4_t vmod_sub(int32x4_t a, int32x4_t b, int32x4_t P, int32x4_t half_mod) {
    int32x4_t result = vaddq_s32(vsubq_s32(a, b), P);
    return vmod_adjust(result, half_mod, P);
}

// Vectorized Barrett multiplication (Q1.31 fixed-point approach)
static inline int32x4_t vBarrett_mul_const(int32x4_t a, int32x4_t c, int32x4_t N, int32x4_t mu_c) {
    int32x4_t z = vmulq_s32(a, c);
    int32x4_t t = vqrdmulhq_s32(a, mu_c);
    z = vmlsq_s32(z, t, N);
    return z;
}

// Declare the external CT_4pts_vec and GS_4pts_vec
void CT_4pts_vec(
    const int32_t P[4],
    int32_t t0[4], int32_t t1[4], int32_t t2[4], int32_t t3[4],
    const int32_t omega_1[4], const int32_t omega_2[4], const int32_t omega_3[4],
    const int32_t mu_omega_1[4], const int32_t mu_omega_2[4], const int32_t mu_omega_3[4]);

void GS_4pts_vec(
    const int32_t P[4],
    int32_t t0[4], int32_t t1[4], int32_t t2[4], int32_t t3[4],
    const int32_t omega_1[4], const int32_t omega2[4], const int32_t omega_3[4],
    const int32_t mu_omega_1[4], const int32_t mu_omega_2[4], const int32_t mu_omega_3[4]);

#endif // VEC_BUTTERFLY_H
