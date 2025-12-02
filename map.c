#ifndef MAP_C
#define MAP_C

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "bseq.h"
#include "trfx.h"
#include "trfrun.h"
#include "tr30dat.h"
#include "tr30dat.c"
#include "ksort.h"  // 包含ksort头文件

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
    size_t *seq_lengths; // 存储序列长度的数组
    int *index_array;    // 索引数组用于排序
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

// 简单的基数排序实现，专门用于size_t类型的排序
static void radix_sort_size_t(size_t *array, int *index, int n) {
    if (n <= 1) return;
    
    // 找到最大值来确定位数
    size_t max_val = array[0];
    for (int i = 1; i < n; i++) {
        if (array[i] > max_val) max_val = array[i];
    }
    
    // 临时数组
    int *temp_index = (int*)malloc(n * sizeof(int));
    size_t *temp_array = (size_t*)malloc(n * sizeof(size_t));
    
    // 对每个字节进行排序
    for (int shift = 0; shift < 64; shift += 8) {
        int count[256] = {0};
        
        // 统计每个桶的数量
        for (int i = 0; i < n; i++) {
            count[(array[i] >> shift) & 0xFF]++;
        }
        
        // 计算累积分布
        for (int i = 1; i < 256; i++) {
            count[i] += count[i-1];
        }
        
        // 从后往前放入桶中（保持稳定性）
        for (int i = n-1; i >= 0; i--) {
            int bucket = (array[i] >> shift) & 0xFF;
            temp_index[count[bucket] - 1] = index[i];
            temp_array[count[bucket] - 1] = array[i];
            count[bucket]--;
        }
        
        // 复制回原数组
        memcpy(index, temp_index, n * sizeof(int));
        memcpy(array, temp_array, n * sizeof(size_t));
        
        // 如果已经排序完成，提前退出
        if (shift + 8 >= 64 && (max_val >> (shift + 8)) == 0) {
            break;
        }
    }
    
    free(temp_index);
    free(temp_array);
}

// 基于序列长度的排序函数
static void sort_by_length(step_t *s) {
    int n = s->n_seq;
    
    if (n <= 1) return;
    
    // 提取序列长度
    for (int i = 0; i < n; i++) {
        s->seq_lengths[i] = (size_t)s->seq[i].l_seq;
        s->index_array[i] = i;
    }
    
    // 使用基数排序对索引数组进行排序（按长度降序）
    radix_sort_size_t(s->seq_lengths, s->index_array, n);
    
    // 反转顺序（基数排序是升序，我们需要降序）
    for (int i = 0; i < n/2; i++) {
        int temp = s->index_array[i];
        s->index_array[i] = s->index_array[n-1-i];
        s->index_array[n-1-i] = temp;
    }
    
    // 复制到sorted_indices
    memcpy(s->sorted_indices, s->index_array, n * sizeof(int));
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
            s->seq_lengths = (size_t*)malloc(s->n_seq * sizeof(size_t));
            s->index_array = (int*)malloc(s->n_seq * sizeof(int));
            s->myResult = (Result *)calloc(s->n_seq, sizeof(Result));
            
            if (!s->myResult || !s->sorted_indices || !s->seq_lengths || !s->index_array) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                if (s->seq) free(s->seq);
                if (s->sorted_indices) free(s->sorted_indices);
                if (s->seq_lengths) free(s->seq_lengths);
                if (s->index_array) free(s->index_array);
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
        
        // 使用基数排序进行排序
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
        free(s->seq_lengths);
        free(s->index_array);
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
