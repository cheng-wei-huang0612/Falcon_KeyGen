#include <stdint.h>
#include <stdio.h>

#include "butterfly.h"

/* Cooley-Tuckey FFT algorithm 
 * f: the pointer to the working memory sequence, 
 * this function will perform four independent NTTs on f[0:n], f[n:2n], etc, resp.
 *
 * i: the index of the first operand of the sequence
 *
 * d: 
 *
 * P: the moduli of four independent NTTs
 *
 * omegas: the twiddle factors:
 * the i th NTT uses omegas[2048*i + j], for j from 0 to 2047
 *
 * mu_omegas: the mu value (used in Barrett mul) of the twiddle factors
 *
 */

void Cooley_Tuckey(int32_t *f, uint32_t i, uint32_t d, int32_t *P, int32_t *omegas){
    

    f[i + 2*d] = (f[i + 2*d] * omegas[512]) % P[0];
    f[i + 3*d] = (f[i + 3*d] * omegas[512]) % P[0];

    f[i      ] = (f[i     ] +     f[i + 2*d]) % P[0];
    f[i + 2*d] = (f[i     ] - 2 * f[i + 2*d]) % P[0];

    f[i +   d] = (f[i +  d] +     f[i + 3*d]) % P[0];
    f[i + 3*d] = (f[i +  d] - 2 * f[i + 3*d]) % P[0];


    
    f[i +   d] = (f[i +   d] * omegas[256]) % P[0];
    f[i + 3*d] = (f[i + 3*d] * omegas[768]) % P[0];


    f[i      ] = (f[i      ] +     f[i +   d]) % P[0];
    f[i +   d] = (f[i      ] - 2 * f[i +   d]) % P[0];

    f[i + 2*d] = (f[i + 2*d] +     f[i + 3*d]) % P[0];
    f[i + 3*d] = (f[i + 2*d] - 2 * f[i + 3*d]) % P[0];
}



