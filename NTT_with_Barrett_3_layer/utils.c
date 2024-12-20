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


int number_order(int32_t *temp, uint16_t length, int32_t p) {
    int32_t maximum = 0;
    for (uint16_t i = 0; i < length; i++) {
        int32_t abs_val = temp[i] >= 0 ? temp[i] : -temp[i];
        if (abs_val > maximum) {
            maximum = abs_val;
        }
    }

    int degree = 0;
    while (1) {
        if (degree * (p >> 1) < maximum && maximum < (degree + 1) * (p >> 1)) {
            return degree;
        }
        degree++;
    }
}