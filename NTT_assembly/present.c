#include <stdio.h>
#include <stdint.h>
#include <string.h>


#include "data.c"


// Print the array
void PRINT_SPEC(int32_t *f, uint16_t n, uint16_t d)
{
    printf("\n[");
	for (size_t i = 0; i < n; i+=d)
	{
		printf("%d, ", f[i]);
	}
    printf("\b\b]\n");
}


//void omegas_gen(int32_t omega_2048, int32_t *omegas_2048){
//    for (size_t i = 0; i < 1024; i++)
//    {
//        omegas_2048[i] = mod_exp(omega_2048, i, 33550337);
//    }
//}



/* Signed_Plantard
 *
 * The function receives 32 bits signed A and 32 bits signed B, they are the operands of the multiplication
 * The function receives 32 bits signed P as the modulo number
 * The function receives 64 bits R, which is centered representative of P^(-1) under modulo 2^64.
 *
 * The function outputs the centered representative of A*B*(- 2^(-64)) in modulo P
 */

// Original, this requires 3 multiplications 
int32_t Signed_Plantard(int32_t A, int32_t B, int32_t P, int64_t R){

    int64_t d = (int64_t) A * (int64_t) B;
    int32_t e = ((d * R)>>32) + 16;
    int32_t result = ((int64_t) e * (int64_t) P)>>32;

    return result;
}




// 模冪 (快速冪，適用於有符號數)
int32_t mod_exp(int32_t base, int32_t exp, int32_t mod, int64_t R) {
    int32_t result = 4093;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = Signed_Plantard(result, base, mod, R);
        }

        base = Signed_Plantard(base, base, mod, R);
        exp /= 2;
    }
    return result;
}

extern void butterfly_two_layers_ntt(int32_t *f, uint32_t i, uint32_t d, 
                              int32_t w1, int32_t w2, int32_t w3, 
                              int32_t mod, int64_t R);



void inv_butterfly_two_layers_ntt(int32_t *f, uint32_t i, uint32_t d, int32_t w1, int32_t w2, int32_t w3, int32_t inv2, int32_t mod, int64_t R) {
    // add and sub

    int32_t t1 = f[i] + f[i + d];
    int32_t t2 = f[i] - f[i + d];
    int32_t t3 = f[i + 2*d] + f[i + 3*d];
    int32_t t4 = f[i + 2*d] - f[i + 3*d];

    f[i] = t1;
    f[i + d] = t2;
    f[i + 2*d] = t3;
    f[i + 3*d] = t4;
    

    // inv-twiddle and divided all by 2 
    // w1 = omega_4, w2 = omega_8, omega_4 = omega_4 
    // f[i + d] should be multiplied by -omega_8^3
    // f[i + 3*d] should be multiplied by -omega_8

    f[i+d] = Signed_Plantard(f[i+d], -w2, mod, R);
    f[i+3*d] = Signed_Plantard(f[i+3*d], -w3, mod, R);

    // Divide all by 2

    f[i] = Signed_Plantard(f[i], inv2, mod, R);
    f[i+d] = Signed_Plantard(f[i+d], inv2, mod, R);
    f[i+2*d] = Signed_Plantard(f[i+2*d], inv2, mod, R);
    f[i+3*d] = Signed_Plantard(f[i+3*d], inv2, mod, R);
    

    // add and sub

    t1 = f[i]+ f[i + 2*d];
    t2 = f[i + d] + f[i + 3*d];
    t3 = f[i] - f[i + 2*d];
    t4 = f[i + d] - f[i + 3*d];

    f[i] = t1;
    f[i + d] = t2;
    f[i + 2*d] = t3;
    f[i + 3*d] = t4;

    
    // inv-twiddle and divided all by 2 
    // w1 = omega_4, w2 = omega_8, omega_4 = omega_4 
    // f[i + d] should be multiplied by -omega_4
    // f[i + 3*d] should be multiplied by -omega_4

    f[i + 2*d] = Signed_Plantard(f[i + 2*d], -w1, mod, R);
    f[i + 3*d] = Signed_Plantard(f[i + 3*d], -w1, mod, R);

    // Divide all by 2

    f[i] = Signed_Plantard(f[i], inv2, mod, R);
    f[i+d] = Signed_Plantard(f[i+d], inv2, mod, R);
    f[i+2*d] = Signed_Plantard(f[i+2*d], inv2, mod, R);
    f[i+3*d] = Signed_Plantard(f[i+3*d], inv2, mod, R);

}

void Forward_ntt(int32_t *omegas_2048, int32_t *f, int32_t P, int64_t R){
    uint16_t power_number[256] = {1, 257, 129, 385, 65, 321, 193, 449, 33, 289, 161, 417, 97, 353, 225, 481, 17, 273, 145, 401, 81, 337, 209, 465, 49, 305, 177, 433, 113, 369, 241, 497, 9, 265, 137, 393, 73, 329, 201, 457, 41, 297, 169, 425, 105, 361, 233, 489, 25, 281, 153, 409, 89, 345, 217, 473, 57, 313, 185, 441, 121, 377, 249, 505, 5, 261, 133, 389, 69, 325, 197, 453, 37, 293, 165, 421, 101, 357, 229, 485, 21, 277, 149, 405, 85, 341, 213, 469, 53, 309, 181, 437, 117, 373, 245, 501, 13, 269, 141, 397, 77, 333, 205, 461, 45, 301, 173, 429, 109, 365, 237, 493, 29, 285, 157, 413, 93, 349, 221, 477, 61, 317, 189, 445, 125, 381, 253, 509, 3, 259, 131, 387, 67, 323, 195, 451, 35, 291, 163, 419, 99, 355, 227, 483, 19, 275, 147, 403, 83, 339, 211, 467, 51, 307, 179, 435, 115, 371, 243, 499, 11, 267, 139, 395, 75, 331, 203, 459, 43, 299, 171, 427, 107, 363, 235, 491, 27, 283, 155, 411, 91, 347, 219, 475, 59, 315, 187, 443, 123, 379, 251, 507, 7, 263, 135, 391, 71, 327, 199, 455, 39, 295, 167, 423, 103, 359, 231, 487, 23, 279, 151, 407, 87, 343, 215, 471, 55, 311, 183, 439, 119, 375, 247, 503, 15, 271, 143, 399, 79, 335, 207, 463, 47, 303, 175, 431, 111, 367, 239, 495, 31, 287, 159, 415, 95, 351, 223, 479, 63, 319, 191, 447, 127, 383, 255, 511};

    uint16_t power;
    int32_t w1, w2, w3;




    for (size_t i = 0; i < 256; i+=1)
    {
        butterfly_two_layers_ntt(f, i, 256, omegas_2048[512], omegas_2048[256], omegas_2048[256*3], P,R);
    }
    



    for (size_t i = 0; i < 64; i+=1)
    {
        for (size_t j = 0; j < 4; j+=1)
        {
            power = power_number[64*j];
            w1 = omegas_2048[128 * (power  )];
            w2 = omegas_2048[ 64 * (power  )];
            w3 = omegas_2048[ 64 * (power+8)];
            butterfly_two_layers_ntt(f     , i + j*256, 64, w1, w2, w3, P,R);

        }

    }

    

    for (size_t i=0; i<16; i+=1) 
    {
        for (size_t j = 0; j < 16; j+=1)
        {
            power = power_number[16*j];
            w1 = omegas_2048[32 * (power   )];
            w2 = omegas_2048[16 * (power   )];
            w3 = omegas_2048[16 * (power+32)];
            
            butterfly_two_layers_ntt(f, i + j*64, 16, w1,w2,w3, P,R);
        }
    }


    

    for (size_t i=0; i<4; i+=1) 
    {
        for (size_t j = 0; j < 64; j+=1)
        {
            power = power_number[4*j];
            w1 = omegas_2048[8 * (power    )];
            w2 = omegas_2048[4 * (power    )];
            w3 = omegas_2048[4 * (power+128)];
            
            butterfly_two_layers_ntt(f, i + j*16, 4, w1,w2,w3, P,R);
        }
    }


    for (size_t i=0; i<1; i+=1) 
    {
        for (size_t j = 0; j < 256; j+=1)
        {
            power = power_number[j];
            w1 = omegas_2048[2 * (power    )];
            w2 = omegas_2048[1 * (power    )];
            w3 = omegas_2048[1 * (power+512)];
            
            butterfly_two_layers_ntt(f, i + j*4, 1, w1,w2,w3, P,R);
        }
    }
}


void backward_ntt(int32_t *omegas_2048, int32_t *f, int32_t P, int64_t R, int32_t inv2){

    uint16_t power_number[256] = {1, 257, 129, 385, 65, 321, 193, 449, 33, 289, 161, 417, 97, 353, 225, 481, 17, 273, 145, 401, 81, 337, 209, 465, 49, 305, 177, 433, 113, 369, 241, 497, 9, 265, 137, 393, 73, 329, 201, 457, 41, 297, 169, 425, 105, 361, 233, 489, 25, 281, 153, 409, 89, 345, 217, 473, 57, 313, 185, 441, 121, 377, 249, 505, 5, 261, 133, 389, 69, 325, 197, 453, 37, 293, 165, 421, 101, 357, 229, 485, 21, 277, 149, 405, 85, 341, 213, 469, 53, 309, 181, 437, 117, 373, 245, 501, 13, 269, 141, 397, 77, 333, 205, 461, 45, 301, 173, 429, 109, 365, 237, 493, 29, 285, 157, 413, 93, 349, 221, 477, 61, 317, 189, 445, 125, 381, 253, 509, 3, 259, 131, 387, 67, 323, 195, 451, 35, 291, 163, 419, 99, 355, 227, 483, 19, 275, 147, 403, 83, 339, 211, 467, 51, 307, 179, 435, 115, 371, 243, 499, 11, 267, 139, 395, 75, 331, 203, 459, 43, 299, 171, 427, 107, 363, 235, 491, 27, 283, 155, 411, 91, 347, 219, 475, 59, 315, 187, 443, 123, 379, 251, 507, 7, 263, 135, 391, 71, 327, 199, 455, 39, 295, 167, 423, 103, 359, 231, 487, 23, 279, 151, 407, 87, 343, 215, 471, 55, 311, 183, 439, 119, 375, 247, 503, 15, 271, 143, 399, 79, 335, 207, 463, 47, 303, 175, 431, 111, 367, 239, 495, 31, 287, 159, 415, 95, 351, 223, 479, 63, 319, 191, 447, 127, 383, 255, 511};

    uint16_t power;
    int32_t w1, w2, w3;


    for (size_t i=0; i<1; i+=1) 
    {
        for (size_t j = 0; j < 256; j+=1)
        {
            power = power_number[255 - j];
            w1 = omegas_2048[2 * (power    )];
            w2 = omegas_2048[1 * (power+512)];
            w3 = omegas_2048[1 * (power    )];
            
            inv_butterfly_two_layers_ntt(f, i + j*4, 1, w1,w2,w3,inv2, P,R);
        }
    }

    for (size_t i=0; i<4; i+=1) 
    {
        for (size_t j = 0; j < 64; j+=1)
        {
            power = power_number[252 - 4*j];
            w1 = omegas_2048[8 * (power    )];
            w2 = omegas_2048[4 * (power+128)];
            w3 = omegas_2048[4 * (power    )];
            
            inv_butterfly_two_layers_ntt(f, i + j*16, 4, w1,w2,w3,inv2, P,R);
        }
    }

    for (size_t i=0; i<16; i+=1) 
    {
        for (size_t j = 0; j < 16; j+=1)
        {
            power = power_number[240 -16*j];
            w1 = omegas_2048[32 * (power   )];
            w2 = omegas_2048[16 * (power+32)];
            w3 = omegas_2048[16 * (power   )];
            
            inv_butterfly_two_layers_ntt(f, i + j*64, 16, w1,w2,w3,inv2, P,R);
        }
    }



    for (size_t i = 0; i < 64; i+=1)
    {
        for (size_t j = 0; j < 4; j+=1)
        {
            power = power_number[192 - 64*j];
            w1 = omegas_2048[128 * (power  )];
            w2 = omegas_2048[ 64 * (power+8)];
            w3 = omegas_2048[ 64 * (power  )];
            inv_butterfly_two_layers_ntt(f, i + j*256, 64, w1,w2,w3,inv2, P,R);
        }

    }





    for (size_t i = 0; i < 256; i+=1)
    {
        inv_butterfly_two_layers_ntt(f, i, 256, omegas_2048[512], omegas_2048[256*3], omegas_2048[256], inv2, P, R);
    }
}

int main() {
    // Parameters
    //{'p': 33550337, 'p0i': 16773119, 'R2': 33546244, 'omega_2048': 33547300, 'inv2': 16775169}

    //uint8_t logn = 10;
    //uint16_t n = 1024;
    int32_t P = 33550337; 
    int64_t R = 5764325841882320897; // R = (P^-1) mod 2^64 = 5764325841882320897
    // Plantard_adjust = (-2^64) mod P = 4093
    
    int32_t inv2 = 16773122; // Already in Plantard form.
    

    int32_t omega_2048 = 12430441;


    // Compute the omega_powers : omega_2048_1 omega_2048_2 ,... omega_2048_1023
    int32_t omegas_2048[1024] = {0};

    omegas_2048[0] = 4093;
    for (size_t i = 1; i < 1024; i++)
    {
        omegas_2048[i] = Signed_Plantard(omegas_2048[i-1], omega_2048, P, R);
    }
    

    int32_t f_star[1024] = {0};
    for (uint16_t i = 0; i < (n>>1); i+=1)
    {
        f_star[2*i  ] =  f[2*i  ]; 
        f_star[2*i+1] = -f[2*i+1];
    }
    
    Forward_ntt(omegas_2048, f, P, R);
    Forward_ntt(omegas_2048, f_star, P, R);

    printf("NTT(f) is: \n");
    PRINT_SPEC(f,1024,1);
    printf("NTT(f_star) is: \n");
    PRINT_SPEC(f_star,1024,1);

    printf("The front 8 terms of NTT(f) and NTT(f_star) are: \n");
    PRINT_SPEC(f,8,1);
    PRINT_SPEC(f_star,8,1);


    int32_t N_f[1024] = {0};

    for (size_t i = 0; i < 1024; i+= 2)
    {
        // adjust the f[i] into Plantard form
        f[i+1] = Signed_Plantard(f[i+1], 16752649, P, R);//16797688 = -2^(128) mod P

        N_f[i  ] = Signed_Plantard(f[i], f[i+1], P, R);
        N_f[i+1] = N_f[i];
    }

    printf("NTT(N(f)) is: \n");
    PRINT_SPEC(N_f,1024,1);

    backward_ntt(omegas_2048, N_f, P, R, inv2);
    PRINT_SPEC(N_f,1024,1);
    forward_ntt(omegas_2048, N_f, P, R);




    for (size_t i = 0; i < 1024; i+= 4)
    {
        // adjust the f[i] into Plantard form
        f[i+1] = Signed_Plantard(f[i+1], 16752649, P, R);//16797688 = -2^(128) mod P

        N_f[i  ] = Signed_Plantard(f[i], f[i+1], P, R);
        N_f[i+1] = N_f[i];
    }



    

    


    




    printf("Hello, World!\n");




    
    return 0;
}
