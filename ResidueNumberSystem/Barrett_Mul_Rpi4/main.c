#include <stdio.h>
#include <stdint.h>
#include <arm_neon.h>
#include <math.h>

void Signed_ModMul(int32_t *a, int32_t *b, int32_t *q, int32_t *r) {
    for (int i = 0; i < 4; i++) {
        int64_t temp = (int64_t)a[i] * (int64_t)b[i];
        temp = temp % q[i];

        temp = (int32_t)temp;
        if (temp > (q[i] >> 2)) {
            temp = temp - q[i];
        }

        r[i] = (int32_t)temp;
    }
}



void NEON_Signed_Barrett_ModMul(int32_t *a, int32_t *b, int32_t *mu, int32_t *q, int32_t *r) {
    int32x4_t va = vld1q_s32(a);
    int32x4_t vb = vld1q_s32(b);
    int32x4_t vmu = vld1q_s32(mu);
    int32x4_t vq = vld1q_s32(q);

    int32x4_t vz = vmulq_s32(va, vb);
    int32x4_t vt = vqrdmulhq_s32(va, vmu);
    int32x4_t vtemp = vmlsq_s32(vz, vt, vq);

    //print vz vt and vtemp:
    //int32_t z[4];
    //int32_t t[4];
    //int32_t temp[4];
    //vst1q_s32(z, vz);
    //vst1q_s32(t, vt);
    //vst1q_s32(temp, vtemp);
    //printf("vz = %d %d %d %d\n", z[0], z[1], z[2], z[3]);
    //printf("vt = %d %d %d %d\n", t[0], t[1], t[2], t[3]);
    //printf("vtemp = %d %d %d %d\n", temp[0], temp[1], temp[2], temp[3]);



    vst1q_s32(r, vtemp);
}

void Mu(int32_t *b, int32_t *q, int32_t *mu) {
    for (int i = 0; i < 4; i++) {

        // Perform the calculation
        double shifted_value = (double)b[i] * (1LL << 31); 
        double divided_value = shifted_value / (2.0 * q[i]);
        double floored_value = floor(divided_value + 0.5);
        int result = (int)(2.0 * floored_value); // Cast to integer type

        mu[i] = result;
    }
}


int main () {
    // Operands
    int32_t a[4] = {123,7234,7345,745622};
    int32_t b[4] = {7234,7345,7456,756722};

    int32_t q[4] = {33550337, 33540097, 33538049, 33533953};
    int32_t mu[4] = {463032,470280,477418,48459784};
    int32_t r[4] = {0};

    // Calculate mu
    Mu(b, q, mu);

    printf("The mu(multiplier) values are\n");
    printf("Mu = %d %d %d %d\n\n", mu[0], mu[1], mu[2], mu[3]);





    Signed_ModMul(a, b, q, r);
    printf("The standard result is\n");
    printf("r = %d %d %d %d\n\n", r[0], r[1], r[2], r[3]);

    printf("Barrett with NEON intrinsics\n");
    NEON_Signed_Barrett_ModMul(a, b, mu, q, r);
    printf("r = %d %d %d %d\n", r[0], r[1], r[2], r[3]);



    return 0;
}
