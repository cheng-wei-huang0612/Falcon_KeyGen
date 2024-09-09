#include <stdio.h>
#include <stdint.h>

uint16_t logn = 8;
uint16_t n = 256;

uint32_t p = 16760833;
uint32_t p0i = 4043292671; // p0i = -1/p mod 2^32
uint32_t r2 = 13696128; // r squared in modulo p, r = 2^32.
uint32_t inv2 = 2097024; // in Montgomery form


uint32_t f[256] = {1, 0, 2, 4, 1, 0, 0, 3, 0, 2, 5, 3, 0, 1, 2, 5, 3, 0, 3, 2, 1, 0, 1, 2, 2, 1, 3, 0, 0, 1, 1, 0, 2, 4, 1, 0, 1, 2, 3, 1, 0, 2, 0, 1, 0, 5, 2, 3, 2, 3, 1, 2, 2, 2, 3, 0, 4, 2, 1, 1, 1, 1, 1, 5, 4, 0, 0, 3, 0, 6, 0, 7, 4, 7, 2, 1, 1, 3, 1, 3, 4, 2, 4, 2, 3, 2, 4, 6, 3, 6, 3, 1, 0, 0, 1, 0, 3, 4, 3, 0, 0, 2, 7, 2, 6, 5, 4, 1, 1, 0, 1, 1, 3, 0, 0, 0, 2, 0, 1, 5, 1, 1, 1, 2, 2, 0, 4, 1, 0, 1, 0, 2, 5, 2, 2, 3, 2, 1, 4, 2, 3, 4, 8, 3, 1, 0, 3, 0, 5, 3, 1, 1, 1, 6, 7, 0, 1, 2, 2, 5, 1, 6, 6, 0, 3, 0, 1, 4, 1, 0, 6, 0, 0, 5, 0, 2, 3, 0, 2, 0, 1, 3, 1, 0, 1, 1, 0, 7, 0, 1, 2, 2, 0, 3, 0, 1, 0, 0, 1, 3, 0, 1, 5, 4, 1, 0, 0, 4, 1, 4, 5, 3, 2, 3, 4, 3, 0, 0, 1, 1, 1, 0, 4, 2, 3, 4, 2, 0, 1, 2, 0, 3, 3, 2, 4, 2, 1, 2, 1, 2, 3, 0, 5, 5, 2, 2, 4, 3, 1, 1, 5, 1, 1, 0, 4, 3};
uint32_t g[256] = {3, 0, 1, 0, 3, 3, 2, 2, 1, 2, 3, 0, 1, 1, 4, 4, 2, 3, 4, 2, 2, 1, 3, 3, 3, 1, 2, 0, 0, 0, 0, 1, 1, 3, 5, 2, 0, 6, 3, 2, 2, 6, 2, 5, 1, 3, 0, 1, 6, 1, 3, 1, 7, 0, 3, 6, 0, 1, 0, 2, 2, 0, 0, 1, 1, 0, 0, 2, 5, 4, 6, 6, 0, 2, 4, 1, 0, 0, 3, 3, 5, 1, 1, 3, 0, 0, 3, 5, 1, 1, 1, 6, 1, 0, 3, 2, 4, 1, 0, 0, 3, 1, 1, 5, 0, 2, 0, 0, 0, 0, 2, 1, 1, 0, 3, 4, 2, 0, 4, 3, 1, 5, 1, 5, 5, 2, 8, 1, 4, 6, 3, 2, 0, 1, 1, 2, 2, 3, 0, 1, 1, 0, 2, 0, 0, 2, 0, 5, 3, 2, 1, 0, 3, 1, 3, 1, 1, 0, 2, 6, 0, 0, 0, 2, 0, 5, 3, 2, 1, 8, 6, 4, 4, 2, 3, 0, 4, 2, 1, 4, 0, 2, 0, 0, 1, 3, 0, 1, 1, 1, 5, 4, 6, 0, 1, 0, 2, 8, 0, 4, 1, 1, 4, 1, 0, 0, 2, 2, 0, 1, 2, 4, 0, 0, 0, 1, 0, 3, 1, 3, 2, 4, 2, 6, 2, 1, 2, 0, 2, 0, 0, 4, 0, 3, 0, 6, 2, 4, 1, 0, 7, 0, 0, 6, 4, 3, 1, 0, 4, 2, 2, 0, 1, 1, 0, 0};
const uint32_t omegas[64] = {4194048, 14494886, 15152879, 3034156, 16101275, 9421437, 9406032, 3096948, 14719561, 3387958, 11631221, 11686617, 11930405, 1716427, 7304379, 12000572, 6994530, 5999281, 14509801, 13656077, 15369345, 11959950, 10601440, 1099823, 14899773, 2179276, 4598303, 16218269, 10059753, 365533, 9943678, 5884202, 12566785, 2265947, 1607954, 13726677, 659558, 7339396, 7354801, 13663885, 2041272, 13372875, 5129612, 5074216, 4830428, 15044406, 9456454, 4760261, 9766303, 10761552, 2251032, 3104756, 1391488, 4800883, 6159393, 15661010, 1861060, 14581557, 12162530, 542564, 6701080, 16395300, 6817155, 10876631};

uint8_t total_power[32] = {1, 33, 17, 49, 9, 41, 25, 57, 5, 37, 21, 53, 13, 45, 29, 61, 3, 35, 19, 51, 11, 43, 27, 59, 7, 39, 23, 55, 15, 47, 31, 63};

static inline uint32_t tbmask(uint32_t x) {
	return (uint32_t)(*(int32_t *)&x >> 31);
}
static inline uint32_t mp_montymul(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i) {
	uint64_t z = (uint64_t)a * (uint64_t)b;
	uint32_t w = (uint32_t)z * p0i;
	uint32_t d = (uint32_t)((z + (uint64_t)w * (uint64_t)p) >> 32) - p;
	return d + (p & tbmask(d));
}
static inline uint32_t mp_montypow(uint32_t base_mont, uint32_t exp, uint32_t p, uint32_t p0i) {
    uint32_t result_mont = 1;    // Initialize result to 1 in Montgomery form

    while (exp > 0) {
        // If the exponent is odd, multiply the result with base
        if (exp & 1) {
            result_mont = mp_montymul(result_mont, base_mont, p, p0i);
        }

        // Square the base: base = base^2 mod p (Montgomery multiplication)
        base_mont = mp_montymul(base_mont, base_mont, p, p0i);

        // Right shift the exponent by 1
        exp >>= 1;
    }

    // Since the input and result are in Montgomery form, return result_mont without additional conversion.
    return result_mont;
}





// static inline uint32_t mp_add(uint32_t a, uint32_t b, uint32_t p) {
// 	uint32_t d = a + b - p;
// 	return d + (p & tbmask(d));
// }

static inline uint32_t mp_add(uint32_t a, uint32_t b, uint32_t p) {
	uint32_t d = a + b;
    if (d>=p)
    {
        d -= p;
        return d;
    }
    return d;
}

// static inline uint32_t mp_sub(uint32_t a, uint32_t b, uint32_t p) {
// 	uint32_t d = a - b;
// 	return d + (p & tbmask(d));
// }

static inline uint32_t mp_sub(uint32_t a, uint32_t b, uint32_t p) {
	uint32_t d = a + p;
    d = d - b;

    if (d>=p)
    {
        d -= p;
        return d;
    }
    return d;
}


void to_mont_form(uint32_t *f, size_t n, uint32_t r2, uint32_t p, uint32_t p0i) {
    for (size_t i = 0; i < n; i++) {
        f[i] = mp_montymul(f[i], r2, p, p0i);
    }
}

void de_mont_form(uint32_t *f, size_t n, uint32_t r2, uint32_t p, uint32_t p0i) {
    for (size_t i = 0; i < n; i++) {
        f[i] = mp_montymul(f[i], 1, p, p0i);
    }
}


void StepDown(int layer, uint32_t *f, int n, uint8_t *total_power, int logn, uint32_t p, uint32_t p0i) {
    // Twiddle factors (apply Montgomery multiplication)
    int step_size = 1 << (logn - 3 - layer);

    int ptr = n >> (layer + 1);
    for (int j = 0; j < (1 << layer); j++) {
        for (int i = 0; i < (n >> (layer + 1)); i++) {
            // Apply the twiddle factors using Montgomery multiplication
            f[ptr + i] = mp_montymul(f[ptr + i], omegas[(n>>(layer+4))*total_power[j*(step_size)]], p, p0i);

        }
        ptr += (n >> layer);
    }

    // Cooley-Tukey butterfly using Montgomery addition and subtraction
    uint32_t f_temp[n]; // Temporary storage for f
    for (int i = 0; i < n; i++) {
        f_temp[i] = f[i];
    }

    int d = (n >> (layer + 1)); // Distance
    ptr = 0;
    for (int j = 0; j < (1 << layer); j++) {
        for (int i = 0; i < (n >> (layer + 1)); i++) {
            f[ptr + i] = mp_add(f_temp[ptr + i], f_temp[ptr + d + i], p); // Addition (butterfly operation)
            f[ptr + d + i] = mp_sub(f_temp[ptr + i], f_temp[ptr + d + i], p); // Subtraction (butterfly operation)
        }
        ptr += (n >> layer);
    }
}

void StepUp(int layer, uint32_t *f, int n, uint8_t *total_power, int logn, uint32_t p, uint32_t p0i) {
    // Cooley-Tukey butterfly using Montgomery addition and subtraction
    uint32_t f_temp[n]; // Temporary storage for f
    for (int i = 0; i < n; i++) {
        f_temp[i] = f[i];
    }

    int d = (n >> (layer + 1)); // Distance
    int ptr = 0;
    for (int j = 0; j < (1 << layer); j++) {
        for (int i = 0; i < (n >> (layer + 1)); i++) {
            f[ptr + i] = mp_add(f_temp[ptr + i], f_temp[ptr + d + i], p); // Addition (butterfly operation)
            f[ptr + d + i] = mp_sub(f_temp[ptr + i], f_temp[ptr + d + i], p); // Subtraction (butterfly operation)
        }
        ptr += (n >> layer);
    }

    // Twiddle factors (inverse twiddle)
    int step_size = 1 << (logn - 3 - layer);
    ptr = n >> (layer + 1);
    
    for (int j = 0; j < (1 << layer); j++) {
        for (int i = 0; i < (n >> (layer + 1)); i++) {
            // Apply the twiddle factors using Montgomery multiplication (inverse twiddle factors)
            f[ptr + i] = mp_montymul(f[ptr + i], omegas[(n>>2)-((n>>(layer+4))*total_power[j*(step_size)])], p, p0i);
        }
        ptr += (n >> layer);
    }

    // Normalization (dividing by 2 in Montgomery space)
    for (int i = 0; i < n; i++) {
        f[i] = mp_montymul(f[i], inv2, p, p0i); // Montgomery reduction by 2
    }
}


void TMVP_base(uint32_t *f_small_8,uint32_t *g_small_8, uint32_t twiddle, uint32_t *result) {
    int toep[15];  // The Toeplitz array will have 15 elements (8 original + 7 more)

    // Generate the Toeplitz array by multiplying f_small_8 by twiddle
    for (int i = 1; i < 8; i++) {
        toep[i - 1] = mp_montymul(f_small_8[i], twiddle, p, p0i);
    }

    for (int i = 0; i < 8; i++) {
        toep[i + 7] = f_small_8[i];
    }

    // Perform TMVP (Toeplitz Matrix Vector Product)
    for (int i = 0; i < 8; i++) {
        result[i] = 0;
        for (int j = 0; j < 8; j++) {
            result[i] += mp_montymul(g_small_8[j],toep[i + 7 - j],p,p0i);
        }
    }
}


//extern void TMVP_base_int(int *f_small_8, int *g_small_8, int *twiddle, int *result);

int main() {
    uint32_t h[256] = {0};  // Initialize h with zeros
    

    

    // mont_form of f and g
    to_mont_form(f,n,r2,p,p0i);
    to_mont_form(g,n,r2,p,p0i);

    StepDown(0,f,n,total_power,logn,p,p0i);
    StepDown(1,f,n,total_power,logn,p,p0i);
    StepDown(2,f,n,total_power,logn,p,p0i);
    StepDown(3,f,n,total_power,logn,p,p0i);
    StepDown(4,f,n,total_power,logn,p,p0i);

    StepDown(0,g,n,total_power,logn,p,p0i);
    StepDown(1,g,n,total_power,logn,p,p0i);
    StepDown(2,g,n,total_power,logn,p,p0i);
    StepDown(3,g,n,total_power,logn,p,p0i);
    StepDown(4,g,n,total_power,logn,p,p0i);


    // TMVP operation on blocks of 8
    for (int i = 0; i < n; i += 8) {
        uint32_t result[8] = {0};  // Temporary result for TMVP_base
        uint32_t twiddle = omegas[i / 8];  // Compute omega^total_power
        TMVP_base(&f[i], &g[i], omegas[total_power[i>>3]], result);
        
        // Copy the result back to h
        for (int j = 0; j < 8; j++) {
            h[i + j] = result[j];
        }
    }




    StepUp(4,h,n,total_power,logn,p,p0i);
    printf("\n\n\n");

    de_mont_form(h,n,r2,p,p0i);

    // Print the result array h
    printf("[");
    for (int i = 0; i < n; i++) {
        printf("%u, ", h[i]);
    }
    printf("\b\b]");

    // StepUp(3,h,n,total_power,logn,p,p0i);
    // StepUp(2,h,n,total_power,logn,p,p0i);
    // StepUp(1,h,n,total_power,logn,p,p0i);
    // StepUp(0,h,n,total_power,logn,p,p0i);



    // StepUp(4,f,n,total_power,logn,p,p0i);
    // StepUp(3,f,n,total_power,logn,p,p0i);
    // StepUp(2,f,n,total_power,logn,p,p0i);
    // StepUp(1,f,n,total_power,logn,p,p0i);
    // StepUp(0,f,n,total_power,logn,p,p0i);


    // StepUp(4,g,n,total_power,logn,p,p0i);
    // StepUp(3,g,n,total_power,logn,p,p0i);
    // StepUp(2,g,n,total_power,logn,p,p0i);
    // StepUp(1,g,n,total_power,logn,p,p0i);
    // StepUp(0,g,n,total_power,logn,p,p0i);




    // de_mont_form(f,n,r2,p,p0i);
    // de_mont_form(g,n,r2,p,p0i);

    // printf("[");
    // for (size_t i = 0; i < 256; i++)
    // {
    //     printf("%d, ", g[i]);
    // }
    // printf("\b\b]");
    

    return 0;
}