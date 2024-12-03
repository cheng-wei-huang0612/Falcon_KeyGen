#include <stdio.h>
#include <stdlib.h>
#include <omp.h>  // 用於多核心並行

#define ARRAY_SIZE 1000000
#define PREFETCH_DISTANCE 64  // 預取距離（字節）

// 單核心處理函數
void single_core_compute(int *data, int size) {
    for (int i = 0; i < size; i++) {
        data[i] = data[i] * data[i];
    }
}

// 單核心處理函數（帶預取）
void single_core_compute_with_prefetch(int *data, int size) {
    for (int i = 0; i < size; i++) {
        // 預取未來的數據
        if (i + PREFETCH_DISTANCE / sizeof(int) < size) {
            __builtin_prefetch(&data[i + PREFETCH_DISTANCE / sizeof(int)], 0, 1);
        }
        data[i] = data[i] * data[i];
    }
}

// 多核心處理函數（無預取）
void multi_core_compute(int *data, int size) {
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        data[i] = data[i] * data[i];
    }
}

// 多核心處理函數（帶預取）
void multi_core_compute_with_prefetch(int *data, int size) {
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        // 預取未來的數據
        if (i + PREFETCH_DISTANCE / sizeof(int) < size) {
            __builtin_prefetch(&data[i + PREFETCH_DISTANCE / sizeof(int)], 0, 1);
        }
        data[i] = data[i] * data[i];
    }
}

int main() {
    int *data = (int *)malloc(ARRAY_SIZE * sizeof(int));
    int *data_copy = (int *)malloc(ARRAY_SIZE * sizeof(int));
    
    // 初始化數據
    for (int i = 0; i < ARRAY_SIZE; i++) {
        data[i] = i + 1;
        data_copy[i] = i + 1;
    }

    // 測試單核心無預取
    double start_time = omp_get_wtime();
    single_core_compute(data, ARRAY_SIZE);
    double end_time = omp_get_wtime();
    printf("Single core without prefetch: %f seconds\n", end_time - start_time);

    // 測試單核心帶預取
    start_time = omp_get_wtime();
    single_core_compute_with_prefetch(data_copy, ARRAY_SIZE);
    end_time = omp_get_wtime();
    printf("Single core with prefetch: %f seconds\n", end_time - start_time);

    // 測試多核心無預取
    start_time = omp_get_wtime();
    multi_core_compute(data, ARRAY_SIZE);
    end_time = omp_get_wtime();
    printf("Multi-core without prefetch: %f seconds\n", end_time - start_time);

    // 測試多核心帶預取
    start_time = omp_get_wtime();
    multi_core_compute_with_prefetch(data_copy, ARRAY_SIZE);
    end_time = omp_get_wtime();
    printf("Multi-core with prefetch: %f seconds\n", end_time - start_time);

    free(data);
    free(data_copy);

    return 0;
}