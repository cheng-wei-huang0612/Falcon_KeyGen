#include <stdio.h>
#include <stdint.h>



// Print the array
void PRINT_SPEC(int32_t *f, uint16_t start, uint16_t n, uint16_t d)
{
    printf("\n[");
	for (size_t i = 0; i < n; i+=d)
	{
		printf("%d, ", f[start+i]);
	}
    printf("\b\b]\n");
}
