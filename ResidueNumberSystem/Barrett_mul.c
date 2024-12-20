
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h> // 提供 atoi 函數

// Barrett Reduction Multiplication
int32_t Barrett_mul(int32_t a, int32_t b, int32_t mod, int32_t mu_b) {
    int32_t z = a * b;
    int32_t t_high = (int32_t)(((int64_t)a * (int64_t)(mu_b << 1)) >> 32);
    return (z - t_high * mod);
}

int main(int argc, char *argv[]) {
    if (argc != 5) { // 檢查參數數量
        printf("Usage: ./Barrett_mul <a> <b> <mod> <mu_b>\n");
        return 1; // 錯誤返回值
    }

    // 將命令列參數轉為整數
    int32_t a = atoi(argv[1]);
    int32_t b = atoi(argv[2]);
    int32_t mod = atoi(argv[3]);
    int32_t mu_b = atoi(argv[4]);

    // 計算結果
    int32_t result = Barrett_mul(a, b, mod, mu_b);
    printf("Result: %d\n", result);

    return 0;
}
