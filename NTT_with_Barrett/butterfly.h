// This .h file defined CT butterfly algorithm
// The modular multiplication is 3-instruction Barrett
#include <stdint.h>

void Cooley_Tuckey(int32_t *f, uint32_t i, uint32_t d, int32_t *P, int32_t *omegas, int32_t *mu_omegas);

