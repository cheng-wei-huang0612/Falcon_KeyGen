#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <arm_neon.h>
#include "data.h"
#include "utils.h"




int32_t mod_mul_standard(int32_t a, int32_t b, int32_t mod) {
    return (int32_t)(((int64_t)a * b) % mod);
}

int32_t mod_add(int32_t a, int32_t b, int32_t mod) {
    return (a + b) % mod;
}

int32_t mod_sub(int32_t a, int32_t b, int32_t mod) {
    return (a - b + mod) % mod;
}


int32_t mod_exp(int32_t base, int32_t exp, int32_t mod, int64_t R) {
    int32_t result = 4093;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = mod_mul_standard(result, base, mod);
        }

        base = mod_mul_standard(base, base, mod);
        exp /= 2;
    }
    return result;
}




void Cooley_Tuckey(int32_t *f, const uint32_t logn_top, uint32_t logn, uint32_t i, const int32_t *P, const int32_t *omegas){

    uint32_t d = 1 << (logn-2);    // the distance between each point

    for (int logj=0; logj<4; logj += 1){
    int j = logj << logn_top;
    
    
    f[j + i + 2*d] = mod_mul_standard(f[j + i + 2*d], omegas[2*j + 512], P[j]);
    f[j + i + 3*d] = mod_mul_standard(f[j + i + 3*d], omegas[2*j + 512], P[j]);  


    f[j + i      ] = (f[j + i     ] +     f[j + i + 2*d]) % P[j];
    f[j + i + 2*d] = (f[j + i     ] - 2 * f[j + i + 2*d]) % P[j];

    f[j + i +   d] = (f[j + i +  d] +     f[j + i + 3*d]) % P[j];
    f[j + i + 3*d] = (f[j + i +  d] - 2 * f[j + i + 3*d]) % P[j];



    f[j + i +   d] = mod_mul_standard(f[j + i +   d], omegas[2*j + 256], P[j]);
    f[j + i + 3*d] = mod_mul_standard(f[j + i + 3*d], omegas[2*j + 768], P[j]);


    f[j + i      ] = (f[j + i      ] +     f[j + i +   d]) % P[j];
    f[j + i +   d] = (f[j + i      ] - 2 * f[j + i +   d]) % P[j];

    f[j + i + 2*d] = (f[j + i + 2*d] +     f[j + i + 3*d]) % P[j];
    f[j + i + 3*d] = (f[j + i + 2*d] - 2 * f[j + i + 3*d]) % P[j];
    }
}





int main() {
    // Parameters
    // {'P': 33550337, 'R': 5764325841882320897, 'Omega_2048': 33547300, 'Omega_2048_Plan': 12430441, 'Inv2': 16775169, 'Inv2_Plan': 16773122}
    // {'P': 33540097, 'R': 274879891390871553, 'Omega_2048': 33533628, 'Omega_2048_Plan': 33049065, 'Inv2': 16770049, 'Inv2_Plan': 14838770} 
    // {'P': 33538049, 'R': 46165194950328321, 'Omega_2048': 33510748, 'Omega_2048_Plan': 7340307, 'Inv2': 16769025, 'Inv2_Plan': 335848} 
    // {'P': 33533953, 'R': 10511127271236980737, 'Omega_2048': 33448611, 'Omega_2048_Plan': 4349711, 'Inv2': 16766977, 'Inv2_Plan': 17747916} 

    // uint8_t logn = 10;
    // uint16_t n = 1024;
    
    
    printf("Hello, World!\n");
    for (int i=0; i<4;i++){
        printf("P[%d]: %d\n", i, P[i]);
    }
    

    int32_t temp[4096] = {0};
    for (int i=0; i<1024; i++){
        temp[i] = f[i];
        temp[i+1024] = f[i];
        temp[i+2048] = f[i];
        temp[i+3072] = f[i];
    }
    

    for (int j =0; j<4096; j+=1024){
        PRINT_SPEC(temp, j, 1024, 256);
    }

    Cooley_Tuckey(temp, logn_top, 10, 0, P, omegas);
    
    for (int j =0; j<4096; j+=1024){
        PRINT_SPEC(temp, j, 1024, 256);
    }


    
    



    
    printf("Byebye, World!\n");
    return 0;
}
