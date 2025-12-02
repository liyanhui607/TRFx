#ifndef MAP_C
#define MAP_C

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bseq.h"
#include "trfx.h"
#include "trfrun.h"
#include "tr30dat.h"
#include "tr30dat.c"

// 修复函数签名与头文件一致
void init_readonly(readonly_vars_struct *ro_vars, const trf_opt *opt)
{
    ro_vars->maxwraplength = opt->maxwraplength;
    ro_vars->alpha = opt->match;
    ro_vars->beta = opt->mismatch;
    ro_vars->delta = opt->indel;
    ro_vars->PM = opt->PM;
    ro_vars->PI = opt->PI;
    ro_vars->minscore = opt->minscore;
    ro_vars->maxperiod = opt->maxperiod;
    ro_vars->Min_Distance_Entries = 20;
    ro_vars->Min_Distance_Window = 20;
    ro_vars->SM = NULL;
    ro_vars->waitdata = NULL;
    ro_vars->sumdata = NULL;
    ro_vars->NTS = 3;
    ro_vars->Index = NULL;

    for (int i = 0; i < 10; i++) {
        ro_vars->Tuplesize[i] = 0;
    }
}

/**************************
 * Multi-threaded mapping *
 **************************/

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void *, int, void *), void *shared_data, int n_steps);

typedef struct {
    int batch_size, n_processed, n_threads;
    const trf_opt *opt;
    bseq_file_t *fp;
    readonly_vars_struct *ro_vars;
} pipeline_t;

typedef struct {
    const pipeline_t *p;
    int n_seq;
    bseq1_t *seq;
    Result *myResult;
    int *sorted_indices;
} step_t;

static void worker_for(void *_data, long i, int tid) {
    step_t *step = (step_t *)_data;
    int original_idx = step->sorted_indices[i];
    step->myResult[original_idx].GlobalIndexList = NULL;
    step->myResult[original_idx].GlobalIndexListTail = NULL;
    step->myResult[original_idx] = TRF(&step->seq[original_idx], step->myResult[original_idx], step->p->opt, step->p->ro_vars);
    
    if (NULL == step->myResult[original_idx].GlobalIndexList) {
        return;
    }
}

// 简单的插入排序（用于小数组）
static void insertion_sort(int *indices, bseq1_t *seq, int n) {
    for (int i = 1; i < n; i++) {
        int key = indices[i];
        int j = i - 1;
        
        // 按序列长度降序排列
        while (j >= 0 && seq[indices[j]].l_seq < seq[key].l_seq) {
            indices[j + 1] = indices[j];
            j--;
        }
        indices[j + 1] = key;
    }
}

// 快速排序分区函数
static int partition(int *indices, bseq1_t *seq, int low, int high) {
    int pivot_idx = indices[high];
    int pivot_length = seq[pivot_idx].l_seq;
    int i = low - 1;
    
    for (int j = low; j < high; j++) {
        // 降序排列：当前元素长度大于基准时交换
        if (seq[indices[j]].l_seq > pivot_length) {
            i++;
            int temp = indices[i];
            indices[i] = indices[j];
            indices[j] = temp;
        }
    }
    
    int temp = indices[i + 1];
    indices[i + 1] = indices[high];
    indices[high] = temp;
    
    return i + 1;
}

// 堆排序调整函数
static void heapify(int *indices, bseq1_t *seq, int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    
    if (left < n && seq[indices[left]].l_seq > seq[indices[largest]].l_seq) {
        largest = left;
    }
    
    if (right < n && seq[indices[right]].l_seq > seq[indices[largest]].l_seq) {
        largest = right;
    }
    
    if (largest != i) {
        int temp = indices[i];
        indices[i] = indices[largest];
        indices[largest] = temp;
        
        heapify(indices, seq, n, largest);
    }
}

// 内省排序（结合快速排序、堆排序和插入排序）
static void introsort(int *indices, bseq1_t *seq, int n, int max_depth) {
    if (n <= 16) {
        insertion_sort(indices, seq, n);
        return;
    }
    
    if (max_depth == 0) {
        // 使用堆排序避免快速排序的最坏情况
        for (int i = n / 2 - 1; i >= 0; i--) {
            heapify(indices, seq, n, i);
        }
        
        for (int i = n - 1; i > 0; i--) {
            int temp = indices[0];
            indices[0] = indices[i];
            indices[i] = temp;
            
            heapify(indices, seq, i, 0);
        }
        return;
    }
    
    int pivot = partition(indices, seq, 0, n - 1);
    
    introsort(indices, seq, pivot, max_depth - 1);
    introsort(indices + pivot + 1, seq, n - pivot - 1, max_depth - 1);
}

// 计算最大递归深度（2*log2(n)）
static int calculate_max_depth(int n) {
    int depth = 0;
    while (n > 1) {
        n >>= 1;
        depth++;
    }
    return depth * 2;
}

// 基于内省排序的排序函数
static void sort_by_length(step_t *s) {
    int n = s->n_seq;
    
    if (n <= 1) return;
    
    // 计算最大递归深度
    int max_depth = calculate_max_depth(n);
    
    // 使用内省排序
    introsort(s->sorted_indices, s->seq, n, max_depth);
}

static void *worker_pipeline(void *shared, int step, void *in) {
    pipeline_t *p = (pipeline_t *)shared;
    
    if (step == 0) {
        step_t *s;
        s = (step_t *)calloc(1, sizeof(step_t));
        s->seq = bseq_read(p->fp, p->batch_size, &s->n_seq);

        if (s->seq) {
            s->p = p;
            s->sorted_indices = (int*)malloc(s->n_seq * sizeof(int));
            s->myResult = (Result *)calloc(s->n_seq, sizeof(Result));
            
            if (!s->myResult || !s->sorted_indices) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                if (s->seq) free(s->seq);
                if (s->sorted_indices) free(s->sorted_indices);
                free(s);
                return NULL;
            }
            
            // 初始化索引
            for (int i = 0; i < s->n_seq; ++i) {
                s->sorted_indices[i] = i;
            }
            return s;
        } else {
            free(s);
        }
    }
    else if (step == 1) {
        step_t *s = (step_t *)in;
        
        // 使用内省排序进行排序
        sort_by_length(s);
        
        kt_for(p->n_threads, worker_for, in, s->n_seq);
        return in;
    }
    else if (step == 2) {
        step_t *s = (step_t *)in;
        for (int i = 0; i < s->n_seq; ++i) {
            printResult(stdout, &s->seq[i], s->myResult[i]);
            freeResult(s->myResult[i]);
            free(s->seq[i].seq);
            free(s->seq[i].name);
        }
        free(s->myResult);
        free(s->seq);
        free(s->sorted_indices);
        free(s);
    }
    return 0;
}

// 修复函数签名与头文件一致
int trf_search_file(const char *fn, const trf_opt *opt, int n_threads, long tbatch_size, readonly_vars_struct *ro_vars) {
    pipeline_t pl;
    memset(&pl, 0, sizeof(pipeline_t));
    pl.fp = bseq_open(fn);
    if (pl.fp == 0) return -1;
    pl.opt = opt;
    pl.n_threads = n_threads;
    pl.batch_size = (int)tbatch_size;  // 转换为int
    pl.ro_vars = ro_vars;
    kt_pipeline(n_threads == 1 ? 1 : 2, worker_pipeline, &pl, 3);
    bseq_close(pl.fp);
    return 0;
}

#endif
