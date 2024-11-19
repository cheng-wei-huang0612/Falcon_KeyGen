#include <stdint.h>
#include <stdio.h>

void barrett_mul_neon(uint32_t *a, uint32_t *b, uint32_t *N, uint32_t *mu, uint32_t *z);

int main() {
    uint32_t a[4] = {123456789, 987654321, 192837465, 564738291};
    uint32_t b[4] = {112233445, 556677889, 998877665, 443322110};
    uint32_t N[4] = {1000000007, 1000000009, 1000000021, 1000000033};
    uint32_t mu[4];
    uint32_t z[4];

    // Precompute mu[i] = floor((b[i] * R) / N[i])
    // Assuming R = 2^31 (for Q31 fixed-point representation)
    uint64_t R = (uint64_t)1 << 31;
    for (int i = 0; i < 4; i++) {
        mu[i] = (uint32_t)(((uint64_t)b[i] * R) / N[i]);
    }

    barrett_mul_neon(a, b, N, mu, z);

    printf("Result z:\n");
    for (int i = 0; i < 4; i++) {
        printf("z[%d] = %u\n", i, z[i]);
    }

    return 0;
}
