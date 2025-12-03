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


void init_readonly(readonly_vars_struct *ro_vars, const trf_opt *opt)
{
    // 使用结构体赋值一次性初始化所有字段
    *ro_vars = (readonly_vars_struct){
        .maxwraplength = opt->maxwraplength,
        .alpha = opt->match,
        .beta = opt->mismatch,
        .delta = opt->indel,
        .PM = opt->PM,
        .PI = opt->PI,
        .minscore = opt->minscore,
        .maxperiod = opt->maxperiod,
        .Min_Distance_Entries = 20,  // 最小存储位置数
        .Min_Distance_Window = 20,   // 最小距离窗口大小
        .SM = NULL,
        .waitdata = NULL,
        .sumdata = NULL,
        .NTS = 3,
        .Index = NULL
    };
    
    // 明确初始化 Tuplesize 数组为零
    memset(ro_vars->Tuplesize, 0, sizeof(ro_vars->Tuplesize));
}

/**************************
 * Multi-threaded mapping *
 **************************/

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void *, int, void *), void *shared_data, int n_steps);

typedef struct
{
    long batch_size;
	int  n_processed, n_threads;
	const trf_opt *opt;
	bseq_file_t *fp;
	readonly_vars_struct *ro_vars;
} pipeline_t;

typedef struct
{
	const pipeline_t *p;
	int n_seq;
	bseq1_t *seq;
	Result *myResult; // region
	///int *sorted_indices;
} step_t;



static void worker_for(void *_data, long i, int tid)
{
    step_t *step = (step_t *)_data;
    
    // 初始化结果结构
    step->myResult[i] = (Result){
        .GlobalIndexList = NULL,
        .GlobalIndexListTail = NULL
    };
    
    // 执行TRF处理
    step->myResult[i] = TRF(&step->seq[i], step->myResult[i], step->p->opt, step->p->ro_vars);
    
    // 如果结果为NULL则直接返回
    if (step->myResult[i].GlobalIndexList == NULL) {
        return;
    }
}



static void *worker_pipeline(void *shared, int step, void *in)
{
    pipeline_t *p = (pipeline_t *)shared;
    
    switch (step) {
        case 0: { // Step 0: 读取序列
            step_t *s = (step_t *)calloc(1, sizeof(step_t));
            if (!s) {
                fprintf(stderr, "Error: Memory allocation failed for step_t\n");
                return NULL;
            }
            
           // fprintf(stderr, "p->batch_size: %ld\n", p->batch_size);
            s->seq = bseq_read(p->fp, p->batch_size, &s->n_seq);
            if (!s->seq) {
                free(s);
                return NULL;
            }
            
            s->p = p;
            s->myResult = (Result *)calloc(s->n_seq, sizeof(Result));
            if (!s->myResult) {
                fprintf(stderr, "Error: Memory allocation failed for myResult\n");
                free(s->seq);
                free(s);
                return NULL;
            }
            
            return s;
        }
        
        case 1: { // Step 1: 处理序列
            kt_for(p->n_threads, worker_for, in, ((step_t *)in)->n_seq);
            return in;
        }
        
        case 2: { // Step 2: 输出结果并清理
            step_t *s = (step_t *)in;
            for (int i = 0; i < s->n_seq; ++i) {
                printResult(stdout, &s->seq[i], s->myResult[i]);
                freeResult(s->myResult[i]);
                free(s->seq[i].seq);
                free(s->seq[i].name);
            }
            free(s->myResult);
            free(s->seq);
            free(s);
            break;
        }
    }
    
    return NULL;
}


int trf_search_file(const char *fn, const trf_opt *opt, int n_threads, long tbatch_size, readonly_vars_struct *ro_vars)
{
    // 使用复合字面量初始化 pipeline_t 结构
    pipeline_t pl = {
        .fp = bseq_open(fn),
        .opt = opt,
        .n_threads = n_threads,
        .batch_size = tbatch_size,
        .ro_vars = ro_vars
    };

    // 检查文件打开是否成功
    if (!pl.fp) {
        fprintf(stderr, "Error: Failed to open file '%s'\n", fn);
        return -1;
    }

    // 确定流水线阶段数
    const int n_stages = (n_threads == 1) ? 1 : 2;
    
    // 执行处理流水线
    kt_pipeline(n_stages, worker_pipeline, &pl, 3);
    
    // 关闭文件
    bseq_close(pl.fp);
    
    return 0;
}


#endif
