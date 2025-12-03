// GetTopPeriods_cuda.cu - CUDA-accelerated version with auto-tuned block size
// Fixed: Avoids "goto bypasses initialization" errors by separating declaration and initialization
// GetTopPeriods_cuda.cu - CUDA-accelerated version with auto-tuned block size
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_DIST 6000
#define NUMBER_OF_PERIODS 5

__global__ void updateCountsKernel(
    const int *d_flat_history,
    const int *d_offsets,
    const int *d_num,
    int *d_counts,
    int valid_length
) {
    int mer2 = blockIdx.x;
    if (mer2 >= 16 || d_num[mer2] == 0) return;

    extern __shared__ int s_counts[];
    int tid = threadIdx.x;
    int num_threads = blockDim.x;

    int num = d_num[mer2];
    int base = d_offsets[mer2];
    const int *history = d_flat_history + base;

    // shared memory
    for (int i = tid; i < MAX_DIST; i += num_threads) {
        s_counts[i] = 0;
    }
    __syncthreads();

    // count
    for (int idx = tid; idx < num; idx += num_threads) {
        int current_pos = history[idx];
        int j = idx - 1;
        while (j >= 0) {
            int dist = current_pos - history[j];
            if (dist >= MAX_DIST || dist <= 0) break;
            atomicAdd(&s_counts[dist], 1);
            j--;
        }
    }
    __syncthreads();

    // global memory
    for (int i = tid; i < MAX_DIST; i += num_threads) {
        if (s_counts[i] > 0) {
            atomicAdd(&d_counts[i], s_counts[i]);
        }
    }
}

extern "C" int has_actual_gpu() {
    int count = 0;
    return (cudaGetDeviceCount(&count) == cudaSuccess) && (count > 0);
}

extern "C" int GetTopPeriods_cuda(
    unsigned char *h_pattern,
    int length,
    int *toparray,
    int h_Index[256],
    int MAXDISTANCE
) {
    // ===  ===
    int *h_counts = NULL;
    double *h_counts2 = NULL;
    int *h_num = NULL;
    int *h_offsets = NULL;
    int *h_flat_history = NULL;
    int *d_flat_history = NULL;
    int *d_offsets = NULL;
    int *d_num = NULL;
    int *d_counts = NULL;

    int total_size;
    cudaError_t err = cudaSuccess;

    // --- CUDA occupancy tuning variables ---
    int minGridSize = 0;
    int blockSize = 0;
    size_t dynamicSMemSize = MAX_DIST * sizeof(int);
    int gridSize = 16; 

    const int valid_length = length - 2;

    // === detrend vars ===
    double xysum = 0.0, xsum = 0.0, ysum = 0.0, x2sum = 0.0;
    double n = 0.0, denom = 0.0, slope = 0.0;
    int end = 0;

    // === Allocate host memory ===
    h_counts = (int*)calloc(length, sizeof(int));
    if (!h_counts) goto cleanup;

    h_counts2 = (double*)calloc(length, sizeof(double));
    if (!h_counts2) goto cleanup;

    h_num = (int*)calloc(16, sizeof(int));
    if (!h_num) goto cleanup;

    h_offsets = (int*)calloc(16, sizeof(int));
    if (!h_offsets) goto cleanup;

    // === Count mer2 occurrences ===
    for (int i = 0; i <= valid_length; ++i) {
        int c1 = h_Index[h_pattern[i]];
        int c2 = h_Index[h_pattern[i+1]];
        int tupid = (c1 << 2) | c2;
        if (tupid >= 0 && tupid < 16) {
            h_num[tupid]++;
        }
    }

    // === Compute offsets ===
    total_size = 0;
    for (int i = 0; i < 16; ++i) {
        h_offsets[i] = total_size;
        total_size += h_num[i];
    }

    h_flat_history = (int*)malloc(total_size * sizeof(int));
    if (!h_flat_history) goto cleanup;

    // === Fill flat history ===
    memset(h_num, 0, 16 * sizeof(int)); // 

    for (int i = 0; i <= valid_length; ++i) {
        int c1 = h_Index[h_pattern[i]];
        int c2 = h_Index[h_pattern[i+1]];
        int tupid = (c1 << 2) | c2;
        if (tupid >= 0 && tupid < 16) {
            int offset = h_offsets[tupid] + h_num[tupid]++;
            h_flat_history[offset] = i;
        }
    }

    // === Device memory allocation ===
    err = cudaMalloc(&d_flat_history, total_size * sizeof(int));
    if (err != cudaSuccess) goto cleanup;

    err = cudaMalloc(&d_offsets, 16 * sizeof(int));
    if (err != cudaSuccess) goto cleanup;

    err = cudaMalloc(&d_num, 16 * sizeof(int));
    if (err != cudaSuccess) goto cleanup;

    err = cudaMalloc(&d_counts, length * sizeof(int));
    if (err != cudaSuccess) goto cleanup;

    err = cudaMemset(d_counts, 0, length * sizeof(int));
    if (err != cudaSuccess) goto cleanup;

    // === Copy data to device ===
    err = cudaMemcpy(d_flat_history, h_flat_history, total_size * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) goto cleanup;

    err = cudaMemcpy(d_offsets, h_offsets, 16 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) goto cleanup;

    err = cudaMemcpy(d_num, h_num, 16 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) goto cleanup;

    // === block size ===
    err = cudaOccupancyMaxPotentialBlockSize(
        &minGridSize,           // grid size
        &blockSize,             //  block size
        updateCountsKernel,     // kernel 
        dynamicSMemSize,        // shared memory 
        0                       // block size 
    );
    if (err != cudaSuccess) goto cleanup;

    // block size 为1024（根据kernel的设计）
    if (blockSize > 1024) {
        blockSize = 1024;
    }

    // // printf("CUDA: block=%d, grid=%d, smem=%zu bytes\n", blockSize, gridSize, dynamicSMemSize);

    // === Launch kernel ===
    updateCountsKernel<<<gridSize, blockSize, dynamicSMemSize>>>(
        d_flat_history, d_offsets, d_num, d_counts, valid_length
    );

    err = cudaGetLastError();
    if (err != cudaSuccess) goto cleanup;

    cudaDeviceSynchronize();

    // === Copy result back ===
    err = cudaMemcpy(h_counts, d_counts, length * sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) goto cleanup;

    // === Detrend ===
    for (int i = 1; i <= valid_length; ++i) {
        xysum += i * h_counts[i];
        xsum += i;
        ysum += h_counts[i];
        x2sum += i * i;
    }

    n = valid_length;
    denom = n * x2sum - xsum * xsum;
    slope = (denom != 0.0) ? (n * xysum - xsum * ysum) / denom : 0.0;

    for (int i = 1; i <= valid_length; ++i) {
        h_counts2[i] = h_counts[i] - i * slope;
    }

    end = (valid_length < MAXDISTANCE) ? valid_length : MAXDISTANCE;

    // === Find top 5 periods ===
    for (int t = 0; t < NUMBER_OF_PERIODS; ++t) {
        double topval = -1.0;
        int topind = 0;
        for (int i = 1; i <= end; ++i) {
            if (h_counts2[i] > topval) {
                topval = h_counts2[i];
                topind = i;
            }
        }
        toparray[t] = (topind > 0) ? topind : 0;
        if (topind > 0) h_counts2[topind] = -1.0;
    }

    err = cudaSuccess;

cleanup:
    // === free  host mem ===
    if (h_counts) free(h_counts);
    if (h_counts2) free(h_counts2);
    if (h_num) free(h_num);
    if (h_offsets) free(h_offsets);
    if (h_flat_history) free(h_flat_history);

    // === free device mem ===
    cudaFree(d_flat_history);
    cudaFree(d_offsets);
    cudaFree(d_num);
    cudaFree(d_counts);

    return (err == cudaSuccess) ? 0 : 1;
}

