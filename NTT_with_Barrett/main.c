#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <arm_neon.h>
#include "data.h"
#include "utils.h"




int32_t mod_mul_standard(int32_t a, int32_t b, int32_t mod) {
    return (int32_t)(((int64_t)a * b) % mod);
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
    

    int32_t temp[2048] = {0};
    for (int i=0; i<2048; i++){
        temp[i] = f[i];
    }

    printf("temp = ");
    PRINT_SPEC(temp, 2048,1);

    
    



    
    printf("Byebye, World!\n");
    return 0;
}
