#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <arm_neon.h>
#include <stdbool.h>
#include "data.h"
#include "utils.h"
#include "butterfly.h"




int main() {
    // Parameters
    // {'P': 33550337, 'Omega_2048': 33547300, 'Inv2': 16775169}
    // {'P': 33540097, 'Omega_2048': 33533628, 'Inv2': 16770049} 
    // {'P': 33538049, 'Omega_2048': 33510748, 'Inv2': 16769025} 
    // {'P': 33533953, 'Omega_2048': 33448611, 'Inv2': 16766977} 

    // uint8_t logn = 10;
    // uint16_t n = 1024;
    
    
    printf("Hello, World!\n");


    
    
    uint16_t n = 1<<logn_top;
    int32_t temp[8192] = {0};
    for (int i=0; i<1024; i++){
        temp[i] = f[i];
        temp[i+n] = g[i];
    }


    // Perform the (complete) NTT on *f and *g and store the result in the temp[0-1023] and temp[1024-2047]
    



    int32_t *f;
    int32_t *g;
    f = &temp[0];
    g = &temp[n];


    uint16_t depth;
    uint16_t power_index;
    uint16_t hhn;
    uint16_t hhhn;


    depth = 0;
    n = 1 << (logn_top - depth);
    hhn = n>>2;
    power_index = 1;

    for (int j=0; j<1; j++){
        for (int i=0; i<hhn ; i++){
            CT_4pts(P[0], &f[j*n+i], &f[j*n+i+hhn], &f[j*n+i+2*hhn], &f[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
            CT_4pts(P[0], &g[j*n+i], &g[j*n+i+hhn], &g[j*n+i+2*hhn], &g[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
        }
        power_index++;
    }

    depth = 2;
    n = 1 << (logn_top - depth);
    hhn = n>>2;
    power_index = 4;

    for (int j=0; j<4; j++){
        for (int i=0; i<hhn ; i++){
            CT_4pts(P[0], &f[j*n+i], &f[j*n+i+hhn], &f[j*n+i+2*hhn], &f[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
            CT_4pts(P[0], &g[j*n+i], &g[j*n+i+hhn], &g[j*n+i+2*hhn], &g[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
        }
        power_index++;
    }


    depth = 4;
    n = 1 << (logn_top - depth);
    hhhn = n>>3;
    power_index = 16;
    int32_t omega_1;
    int32_t omega_2;
    int32_t omega_3;
    int32_t omega_4;
    int32_t omega_5;
    int32_t omega_6;
    int32_t omega_7;

    int32_t mu_omega_1;
    int32_t mu_omega_2;
    int32_t mu_omega_3;
    int32_t mu_omega_4;
    int32_t mu_omega_5;
    int32_t mu_omega_6;
    int32_t mu_omega_7;
    
    for (int j=0; j<16; j++){
        for (int i=0; i<hhhn ; i++){
            omega_1 = omegas[power_index];
            omega_2 = omegas[power_index<<1];
            omega_3 = omegas[(power_index<<1)+1];
            omega_4 = omegas[power_index<<2];
            omega_5 = omegas[(power_index<<2)+1];
            omega_6 = omegas[(power_index<<2)+2];
            omega_7 = omegas[(power_index<<2)+3];

            mu_omega_1 = mu_omegas[power_index];
            mu_omega_2 = mu_omegas[power_index<<1];
            mu_omega_3 = mu_omegas[(power_index<<1)+1];
            mu_omega_4 = mu_omegas[power_index<<2];
            mu_omega_5 = mu_omegas[(power_index<<2)+1];
            mu_omega_6 = mu_omegas[(power_index<<2)+2];
            mu_omega_7 = mu_omegas[(power_index<<2)+3];


            CT_8pts(P[0], &f[j*n+i], &f[j*n+i+hhhn], &f[j*n+i+2*hhhn], &f[j*n+i+3*hhhn], &f[j*n+i+4*hhhn], &f[j*n+i+5*hhhn], &f[j*n+i+6*hhhn], &f[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);
            CT_8pts(P[0], &g[j*n+i], &g[j*n+i+hhhn], &g[j*n+i+2*hhhn], &g[j*n+i+3*hhhn], &g[j*n+i+4*hhhn], &g[j*n+i+5*hhhn], &g[j*n+i+6*hhhn], &g[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);



        }
        power_index++;
    }




    depth = 7;
    n = 1 << (logn_top - depth);
    hhhn = n>>3;
    power_index = 128;
    for (int j=0; j<128; j++){
        for (int i=0; i<hhhn ; i++){
            omega_1 = omegas[power_index];
            omega_2 = omegas[power_index<<1];
            omega_3 = omegas[(power_index<<1)+1];
            omega_4 = omegas[power_index<<2];
            omega_5 = omegas[(power_index<<2)+1];
            omega_6 = omegas[(power_index<<2)+2];
            omega_7 = omegas[(power_index<<2)+3];

            mu_omega_1 = mu_omegas[power_index];
            mu_omega_2 = mu_omegas[power_index<<1];
            mu_omega_3 = mu_omegas[(power_index<<1)+1];
            mu_omega_4 = mu_omegas[power_index<<2];
            mu_omega_5 = mu_omegas[(power_index<<2)+1];
            mu_omega_6 = mu_omegas[(power_index<<2)+2];
            mu_omega_7 = mu_omegas[(power_index<<2)+3];


            CT_8pts(P[0], &f[j*n+i], &f[j*n+i+hhhn], &f[j*n+i+2*hhhn], &f[j*n+i+3*hhhn], &f[j*n+i+4*hhhn], &f[j*n+i+5*hhhn], &f[j*n+i+6*hhhn], &f[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);
            CT_8pts(P[0], &g[j*n+i], &g[j*n+i+hhhn], &g[j*n+i+2*hhhn], &g[j*n+i+3*hhhn], &g[j*n+i+4*hhhn], &g[j*n+i+5*hhhn], &g[j*n+i+6*hhhn], &g[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);



        }
        power_index++;
    }

    PRINT_SPEC(f, 0, 1024, 1);


    // Now temp is:
    // temp: |--- NTT(f) ---|--- NTT(g) ---|
    //       | 1024         | 1024         |



    n = 1<<logn_top;
    //uint16_t depth;
    depth = 1;
    for (size_t i = 0; i < (n>>depth); i++)
    {
        temp[2*n+i] = mod_mul_standard(f[2*i], f[2*i+1], P[0]);
        temp[2*n+(n>>depth)+i] = mod_mul_standard(g[2*i], g[2*i+1], P[0]);
    }

    // Now temp is organized as follows:
    // temp : |--- NTT(f) ---|--- NTT(g) ---|--- NTT(N(f)) ---|--- NTT(N(g)) ---|
    //        |<----- n ---->|<----- n ---->|<------ hn ----->|<------ hn ----->|





    n = 1<<logn_top;
    uint16_t hn = n>>1;
    depth = 2;
    f = &temp[2*n];
    g = f + hn;
    int32_t *Nf = g + hn;
    int32_t *Ng = Nf + 2*(n>>depth); // 2 is the number of prime modulus that will be used in the next layer.

    for (size_t i = 0; i < (n>>depth); i++)
    {
        Nf[i] = mod_mul_standard(f[2*i], f[2*i+1], P[0]);
        Ng[i] = mod_mul_standard(g[2*i], g[2*i+1], P[0]);
    }



    // Now temp is organized as follows:
    //
    //
    // temp : |--- NTT(f) ---|--- NTT(g) ---|--- NTT(N(f)) ---|--- NTT(N(g)) ---|--- NTT(N^2(f)) under P[0] ---|------------------------------|--- NTT(N^2(g)) under P[0] ---|------------------------------|






    n = 1 << (logn_top-2);
    f = Ng + 2*n;
    g = f + n;
    

    // temp:  | // |--- NTT(N^2(f)) under P[0] ---|------------------------------|--- NTT(N^2(g)) under P[0] ---|------------------------------|------------------------------|------------------------------|
    //             |<-- n ----------------------->|<-- n ----------------------->|<-- n ----------------------->|<-- n ----------------------->|<-- n ----------------------->|<-- n ----------------------->|
    //             ^                                                             ^                                                             ^f                             ^g
    //             Nf                                                            Ng

    // copy

    for (size_t i = 0; i < n; i++)
    {
        f[i] = Nf[i];
        g[i] = Ng[i];
    }

    // Now temp is:
    // temp:  | // |--- NTT(N^2(f)) under P[0] ---|------------------------------|--- NTT(N^2(g)) under P[0] ---|------------------------------|--- NTT(N^2(f)) under P[0] ---|--- NTT(N^2(g)) under P[0] ---|
    // INTT on the last two

    printf("Line 238");
    PRINT_SPEC(f, 0, 256, 1);
    depth = 5;
    n = 1 << (8 - depth);
    hhhn = n>>3;
    power_index = 63;
    for (int j=0; j<32; j++){
        for (int i=0; i<hhhn ; i++){
            omega_1 = -omegas[power_index];
            omega_2 = -omegas[(power_index<<1)+1];
            omega_3 = -omegas[(power_index<<1)];
            omega_4 = -omegas[(power_index<<2)+3];
            omega_5 = -omegas[(power_index<<2)+2];
            omega_6 = -omegas[(power_index<<2)+1];
            omega_7 = -omegas[(power_index<<2)];

            mu_omega_1 = -(mu_omegas[power_index]-2);
            mu_omega_2 = -(mu_omegas[(power_index<<1)+1]-2);
            mu_omega_3 = -(mu_omegas[(power_index<<1)]-2);
            mu_omega_4 = -(mu_omegas[(power_index<<2)+3]-2);
            mu_omega_5 = -(mu_omegas[(power_index<<2)+2]-2);
            mu_omega_6 = -(mu_omegas[(power_index<<2)+1]-2);
            mu_omega_7 = -(mu_omegas[(power_index<<2)]-2);


            GS_8pts(P[0], &f[j*n+i], &f[j*n+i+hhhn], &f[j*n+i+2*hhhn], &f[j*n+i+3*hhhn], &f[j*n+i+4*hhhn], &f[j*n+i+5*hhhn], &f[j*n+i+6*hhhn], &f[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);
            GS_8pts(P[0], &g[j*n+i], &g[j*n+i+hhhn], &g[j*n+i+2*hhhn], &g[j*n+i+3*hhhn], &g[j*n+i+4*hhhn], &g[j*n+i+5*hhhn], &g[j*n+i+6*hhhn], &g[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);



        }
        power_index--;
    }
    
    depth = 2;
    n = 1 << (8 - depth);
    hhhn = n>>3;
    power_index = 7;
    for (int j=0; j<4; j++){
        for (int i=0; i<hhhn ; i++){
            omega_1 = -omegas[power_index];
            omega_2 = -omegas[(power_index<<1)+1];
            omega_3 = -omegas[(power_index<<1)];
            omega_4 = -omegas[(power_index<<2)+3];
            omega_5 = -omegas[(power_index<<2)+2];
            omega_6 = -omegas[(power_index<<2)+1];
            omega_7 = -omegas[(power_index<<2)];

            mu_omega_1 = -(mu_omegas[power_index]-2);
            mu_omega_2 = -(mu_omegas[(power_index<<1)+1]-2);
            mu_omega_3 = -(mu_omegas[(power_index<<1)]-2);
            mu_omega_4 = -(mu_omegas[(power_index<<2)+3]-2);
            mu_omega_5 = -(mu_omegas[(power_index<<2)+2]-2);
            mu_omega_6 = -(mu_omegas[(power_index<<2)+1]-2);
            mu_omega_7 = -(mu_omegas[(power_index<<2)]-2);


            GS_8pts(P[0], &f[j*n+i], &f[j*n+i+hhhn], &f[j*n+i+2*hhhn], &f[j*n+i+3*hhhn], &f[j*n+i+4*hhhn], &f[j*n+i+5*hhhn], &f[j*n+i+6*hhhn], &f[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);
            GS_8pts(P[0], &g[j*n+i], &g[j*n+i+hhhn], &g[j*n+i+2*hhhn], &g[j*n+i+3*hhhn], &g[j*n+i+4*hhhn], &g[j*n+i+5*hhhn], &g[j*n+i+6*hhhn], &g[j*n+i+7*hhhn], omega_1, mu_omega_1, omega_2, mu_omega_2, omega_3, mu_omega_3, omega_4, mu_omega_4, omega_5, mu_omega_5, omega_6, mu_omega_6, omega_7, mu_omega_7);



        }
        power_index--;
    }
    printf("Line 302");
    PRINT_SPEC(f, 0, 256, 1);


    depth = 0;
    n = 1 << (8 - depth);
    hhn = n>>2;
    power_index = 1;
    for (int j=0; j<1; j++){
        for (int i=0; i<hhn ; i++){
            GS_4pts(P[0], &f[j*n+i], &f[j*n+i+hhn], &f[j*n+i+2*hhn], &f[j*n+i+3*hhn], -omegas[power_index], -omegas[(power_index<<1)+1], -omegas[power_index<<1], -mu_omegas[power_index], -mu_omegas[(power_index<<1)+1], -mu_omegas[power_index<<1]);
            GS_4pts(P[0], &g[j*n+i], &g[j*n+i+hhn], &g[j*n+i+2*hhn], &g[j*n+i+3*hhn], -omegas[power_index], -omegas[(power_index<<1)+1], -omegas[power_index<<1], -mu_omegas[power_index], -mu_omegas[(power_index<<1)+1], -mu_omegas[power_index<<1]);
        }
        power_index--;
    }




    int32_t inv256 = -131056;
    int32_t mu_inv256 = -8388606;

    for (size_t i = 0; i < (256); i++)
    {
        f[i] = Barrett_mul(f[i], inv256, P[0], mu_inv256);
        g[i] = Barrett_mul(g[i], inv256, P[0], mu_inv256);
    }

    PRINT_SPEC(f, 0, 256, 1);

    
    // Now temp is:
    // temp:  | // |--- NTT(N^2(f)) under P[0] ---|------------------------------|--- NTT(N^2(g)) under P[0] ---|------------------------------|--- N^2(f) under P[0] --------|--- N^2(g)) under P[0] -------|
    //             ^                                                             ^                                                             ^f                             ^g
    //             Nf                                                            Ng

    // ECRT (trivial)

    n = 1 << (logn_top-2);
    for (size_t i = 0; i < n; i++)
    {
        Nf[n+i] = mod_add(f[i], 0, P[1]);
        Ng[n+i] = mod_add(g[i], 0, P[1]);
    }

    // Now temp is:
    // temp:  | // |--- NTT(N^2(f)) under P[0] ---|--- N^2(f) under P[1] --------|--- NTT(N^2(g)) under P[0] ---|--- N^2(g)) under P[1] -------|--- N^2(f) under P[0] --------|--- N^2(g)) under P[0] -------|



    n = 1 << (logn_top-2);
    hhn = n>>2;
    power_index = 1;
    for (int j=0; j<1; j++){
        for (int i=0; i<hhn ; i++){
            CT_4pts(P[0], &f[j*n+i], &f[j*n+i+hhn], &f[j*n+i+2*hhn], &f[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
            CT_4pts(P[0], &g[j*n+i], &g[j*n+i+hhn], &g[j*n+i+2*hhn], &g[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
        }
        power_index++;
    }


    n = hhn;
    hhn = n>>2;

    power_index = 4;
    for (int j=0; j<4; j++){
        for (int i=0; i<hhn ; i++){
            CT_4pts(P[0], &f[j*n+i], &f[j*n+i+hhn], &f[j*n+i+2*hhn], &f[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
            CT_4pts(P[0], &g[j*n+i], &g[j*n+i+hhn], &g[j*n+i+2*hhn], &g[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
        }
        power_index++;
    }



    
    n = hhn;
    hhn = n>>2;
    power_index = 16;
    for (int j=0; j<16; j++){
        for (int i=0; i<hhn ; i++){
            CT_4pts(P[0], &f[j*n+i], &f[j*n+i+hhn], &f[j*n+i+2*hhn], &f[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
            CT_4pts(P[0], &g[j*n+i], &g[j*n+i+hhn], &g[j*n+i+2*hhn], &g[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);

        }
        power_index++;
    }



    n = hhn;
    hhn = n>>2;
    power_index = 64;
    for (int j=0; j<64; j++){
        for (int i=0; i<hhn ; i++){
            CT_4pts(P[0], &f[j*n+i], &f[j*n+i+hhn], &f[j*n+i+2*hhn], &f[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
            CT_4pts(P[0], &g[j*n+i], &g[j*n+i+hhn], &g[j*n+i+2*hhn], &g[j*n+i+3*hhn], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1]);
        }
        power_index++;
    }


    // temp:  | // |--- NTT(N^2(f)) under P[0] ---|--- NTT(N^2(f)) under P[1] ---|--- NTT(N^2(g)) under P[0] ---|--- NTT(N^2(g)) under P[1] ---|--- N^2(f) under P[0] --------|--- N^2(g)) under P[0] -------|








    
    printf("\n\nByebye, World!\n");
    return 0;
}
