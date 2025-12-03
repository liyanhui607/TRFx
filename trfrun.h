/*
Tandem Repeats Finder
Copyright (C) 1999-2020 Gary Benson

This file is part of the Tandem Repeats Finder (TRF) program.

TRF is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

TRF is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with TRF.  If not, see <https://www.gnu.org/licenses/>.
*/

/**************************************************************
 *   trfrun.h :  This file contains the code that calls the TRF
 *               algorithm.
 *               It declares all the variables that control how
 *               the program runs and where the output is saved.
 *               If the input file contains more than one
 *               sequence then the input in broken into single
 *               sequences by the control routine and fed to the
 *               algorithm one sequence at a time.
 *               The output is assembled by the control routine
 *               as described in the file readme.txt
 *
 *                                           December 10, 2001
 *


 Last updated Dec 14,2004
 ***************************************************************/

#ifndef TRFRUN_H
#define TRFRUN_H
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "trfx.h"
#include "tr30dat.c"

void PrintError(char *errortext)
{
	fprintf(stderr, "Error: %s\n", errortext);
	return;
}

void FreeList(Index_List *headptr)
{
	Index_List *nextptr;
	Index_List *holdptr;
	nextptr = headptr;
	while (nextptr != NULL)
	{
		holdptr = nextptr;
		nextptr = nextptr->next;
		free(holdptr->pattern);
		free(holdptr);
	}
	return;
}

int *newTags(int length)
{
	int *p = (int *)calloc(length, sizeof(int));
	if (p == NULL)
	{
		PrintError("out of mem");
		return NULL;
	}
	return (p);
}

char *newAlignPairtext(int length)
{
	char *p = (char *)calloc(length, sizeof(char));
	if (p == NULL)
	{
		PrintError("out of mem");
		return NULL;
	}
	return (p);
}

char *newLine(int length)
{
	char *p = (char *)calloc(length, sizeof(char));
	if (p == NULL)
	{
		PrintError("out of mem");
		return NULL;
	}
	return (p);
}

int *newAlignPairindex(int length)
{
	int *p = (int *)calloc(length, sizeof(int));
	if (p == NULL)
	{
		PrintError("out of mem");
		return NULL;
	}
	return (p);
}




void printResult(FILE *destdfp, bseq1_t *seq, Result myResult) 
{
    if (myResult.GlobalIndexList == NULL) {
        return;
    }

    fprintf(destdfp, "@%s\n", seq->name);
    const unsigned char *seq_data = seq->seq - 1;  // 保持无符号类型

    for (Index_List *lpointer = myResult.GlobalIndexList; 
         lpointer != NULL; 
         lpointer = lpointer->next) 
    {
        // 打印主要信息
        fprintf(destdfp, 
               "%d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s ",
               lpointer->first, lpointer->last, lpointer->period,
               lpointer->copies, lpointer->size, lpointer->matches,
               lpointer->indels, lpointer->score, lpointer->acount,
               lpointer->ccount, lpointer->gcount, lpointer->tcount,
               lpointer->entropy, lpointer->pattern);

        // 打印核心序列
        for (int i = lpointer->first; i <= lpointer->last; i++) {
            fputc(seq_data[i], destdfp);
        }

        // 处理侧翼序列
        const int flank_start = max(1, lpointer->first - 50);
        const int flank_end = min(seq->l_seq, lpointer->last + 50);

        // 左侧翼
        fputc(' ', destdfp);
        if (lpointer->first == 1) {
            fputc('.', destdfp);
        } else {
            for (int i = flank_start; i < lpointer->first; i++) {
                fputc(seq_data[i], destdfp);
            }
        }

        // 右侧翼
        fputc(' ', destdfp);
        if (lpointer->last == seq->l_seq) {
            fputc('.', destdfp);
        } else {
            for (int i = lpointer->last + 1; i <= flank_end; i++) {
                fputc(seq_data[i], destdfp);
            }
        }

        fputc('\n', destdfp);
    }
}


void freeResult(Result myResult)
{
	FreeList(myResult.GlobalIndexList);
	myResult.GlobalIndexList = NULL;
	myResult.GlobalIndexListTail = NULL;
}



int IntervalOverlap(Index_List *iptr, Index_List *jptr)
{
    const int beg = iptr->first > jptr->first ? iptr->first : jptr->first;
    const int end = iptr->last < jptr->last ? iptr->last : jptr->last;
    const int overlap = end - beg + 1;
    
    return overlap > 0 ? overlap : 0;
}


int IsRedundant(Index_List *iptr, Index_List *jptr)
{
    const int i_period = iptr->period;
    const int j_period = jptr->period;
    const float i_score = iptr->score;
    const float j_score = jptr->score;

    // 检查两种情况：1) i的周期是j的倍数 2) 周期相同
    return (i_period > j_period && 
            i_period % j_period == 0 && 
            i_score <= 1.1f * j_score) ||
           (i_period == j_period && 
            i_score <= j_score);
}

Index_List *RemoveRedundancy(Index_List *headptr)
{
    if (headptr == NULL || headptr->next == NULL) {
        return headptr; /* 处理空链表或单节点链表 */
    }

    Index_List **iptr_ptr = &headptr;  // 使用二级指针简化删除操作
    
    while (*iptr_ptr != NULL) {
        Index_List *iptr = *iptr_ptr;
        Index_List **jptr_ptr = &(iptr->next);
        int iinterval = iptr->last - iptr->first + 1;
        
        while (*jptr_ptr != NULL) {
            Index_List *jptr = *jptr_ptr;
            int jinterval = jptr->last - jptr->first + 1;
            int overlap = IntervalOverlap(iptr, jptr);
            
            if (overlap == 0) {
                break; /* 没有重叠，提前退出内层循环 */
            }
            
            // 检查iptr相对于jptr是否冗余
            if (overlap >= 0.9 * iinterval && IsRedundant(iptr, jptr)) {
                // 删除iptr并跳出内层循环
                *iptr_ptr = iptr->next;
                free(iptr->pattern);
                free(iptr);
                break; // 跳出内层循环，外层循环会继续处理新的*iptr_ptr
            }
            
            // 检查jptr相对于iptr是否冗余
            if (overlap >= 0.9 * jinterval && IsRedundant(jptr, iptr)) {
                // 删除jptr并继续内层循环
                *jptr_ptr = jptr->next;
                free(jptr->pattern);
                free(jptr);
                continue; // 继续检查下一个jptr
            }
            
            // 移动到下一个jptr
            jptr_ptr = &((*jptr_ptr)->next);
        }
        
        // 移动到下一个iptr（如果iptr没有被删除）
        if (*iptr_ptr == iptr) {
            iptr_ptr = &((*iptr_ptr)->next);
        }
    }
    
    return headptr;
}


Index_List *RemoveBySize(Index_List *headptr, int maxsize)
{
    Index_List **curr = &headptr;  // 使用二级指针简化删除操作
    
    while (*curr != NULL) {
        if ((*curr)->period > maxsize) {
            // 需要删除的节点
            Index_List *to_delete = *curr;
            *curr = to_delete->next;  // 直接从链表中移除节点
            
            free(to_delete->pattern);
            free(to_delete);
        } else {
            // 移动到下一个节点
            curr = &((*curr)->next);
        }
    }
    
    return headptr;
}


Index_List *SortByIndex(Index_List *headptr) 
{
    /* 处理空链表或单节点链表的边界情况 */
    if (headptr == NULL || headptr->next == NULL) {
        return headptr;
    }

    int swapped;
    Index_List **pptr;  // 使用二级指针简化交换操作
    Index_List *curr = NULL;
    Index_List *end = NULL;  // 优化：记录已排序部分的结束位置

    do {
        swapped = 0;
        pptr = &headptr;
        
        while ((curr = *pptr)->next != end) {
            if (curr->first > curr->next->first) {
                /* 交换相邻节点 */
                Index_List *temp = curr->next;
                curr->next = temp->next;
                temp->next = curr;
                *pptr = temp;
                
                swapped = 1;
            }
            pptr = &(*pptr)->next;
        }
        end = curr;  // 记录本轮排序的最后一个节点
        
    } while (swapped);

    return headptr;
}


Index_List *SortByCount(Index_List *headptr) 
{
    /* 处理空链表或单节点链表的边界情况 */
    if (headptr == NULL || headptr->next == NULL) {
        return headptr;
    }

    int swapped;
    Index_List **pptr;  // 使用二级指针简化交换操作
    Index_List *curr = NULL;
    Index_List *end = NULL;  // 优化：记录已排序部分的结束位置

    do {
        swapped = 0;
        pptr = &headptr;
        
        while ((curr = *pptr)->next != end) {
            if (curr->count > curr->next->count) {
                /* 交换相邻节点 */
                Index_List *temp = curr->next;
                curr->next = temp->next;
                temp->next = curr;
                *pptr = temp;
                
                swapped = 1;
            }
            pptr = &(*pptr)->next;
        }
        end = curr;  // 记录本轮排序的最后一个节点
        
    } while (swapped);

    return headptr;
}


Result TRF(bseq1_t *pseq, Result myResult, const trf_opt *opt, readonly_vars_struct *ro_vars)
{
    __attribute__((aligned(64))) thread_local_var_struct thread_local_var = {0};
    thread_local_var_struct *pthread_local_var = &thread_local_var;
    void *base_ptr = NULL;
    int **S = NULL;
    unsigned int i;
    Result ret = myResult;

    // 初始化线程局部变量
    *pthread_local_var = (thread_local_var_struct){
        .OUTPUTcount = 0,
        .counterInSeq = 0,
        .Minsize = 1,
        .allocated_S_cols = 2
    };

    init_bestperiodlist(pthread_local_var);
    const int maxwraplength = min(opt->maxwraplength, pseq->l_seq);
    pthread_local_var->maxwraplength = maxwraplength;
    const int row = maxwraplength;

    // 缓存频繁访问的指针
    unsigned char *Sequence = pseq->seq - 1;
    pthread_local_var->Sequence = Sequence;
    pthread_local_var->Length = pseq->l_seq;
    const int seq_length = pseq->l_seq;
    
    cons_data *pConsensus = &pthread_local_var->Consensus;

    // 优化：使用单个calloc分配共识数据内存
    const size_t consensus_size = 2 * (MAXPATTERNSIZECONSTANT + 1);
    size_t total_consensus_bytes = 9 * consensus_size * sizeof(int) + consensus_size * sizeof(unsigned char);
    unsigned char *consensus_buffer = (unsigned char*)calloc(1, total_consensus_bytes);
    if (!consensus_buffer) {
        PrintError("Unable to allocate consensus memory");
        goto cleanup;
    }

    // 设置共识数据指针
    pConsensus->pattern = consensus_buffer;
    pConsensus->A = (int*)(consensus_buffer + consensus_size * sizeof(unsigned char));
    pConsensus->C = pConsensus->A + consensus_size;
    pConsensus->G = pConsensus->C + consensus_size;
    pConsensus->T = pConsensus->G + consensus_size;
    pConsensus->dash = pConsensus->T + consensus_size;
    pConsensus->insert = pConsensus->dash + consensus_size;
    pConsensus->letters = pConsensus->insert + consensus_size;
    pConsensus->total = pConsensus->letters + consensus_size;

    // 分配S的行指针数组
    S = (int**)malloc((row + 1) * sizeof(int*));
    if (!S) {
        PrintError("Unable to allocate S row pointers");
        goto cleanup_consensus;
    }
    pthread_local_var->S = S;

    // 直接分配对齐的S数据内存 - 移除内存重用逻辑
    const size_t aligned_row_bytes = (MAXBANDWIDTH + 1) * sizeof(int);
    const size_t total_bytes = (row + 1) * aligned_row_bytes;
    
    // 直接分配，不尝试重用
    if (posix_memalign(&base_ptr, 64, total_bytes) != 0) {
        char errmsg[255];
        snprintf(errmsg, sizeof(errmsg), "Failed to allocate %zu bytes (aligned). "
                  "Reduce maxwraplength or MAXBANDWIDTH. (%s:%d)",
                  total_bytes, __FILE__, __LINE__);
        PrintError(errmsg);
        goto cleanup_S_rows;
    }

    // 设置S的行指针 - 使用实际列数 MAXBANDWIDTH + 1
    int *S_data = (int*)base_ptr;
    const int stride = pthread_local_var->allocated_S_cols;
    for (int i = 0; i <= row; i++) {
        S[i] = S_data + i * stride;
    }
    S[0][0] = 1;  // 初始化

    // 分配对齐对内存
    char *textprime = newAlignPairtext(2 * row);
    char *textsecnd = newAlignPairtext(2 * row);
    int *indexprime = newAlignPairindex(2 * row);
    int *indexsecnd = newAlignPairindex(2 * row);
    
    if (!textprime || !textsecnd || !indexprime || !indexsecnd) {
        PrintError("Unable to allocate AlignPair memory");
        goto cleanup_base_ptr;
    }
    
    pairalign *pAlignPair = &pthread_local_var->AlignPair;
    pAlignPair->textprime = textprime;
    pAlignPair->textsecnd = textsecnd;
    pAlignPair->indexprime = indexprime;
    pAlignPair->indexsecnd = indexsecnd;

    // 计算最大距离
    int maxdistance = ro_vars->maxperiod;
    if (maxdistance < 500) maxdistance = 500;
    maxdistance = min(maxdistance, (int)(seq_length * 0.6));
    maxdistance = max(maxdistance, 200);
    pthread_local_var->maxdistance = pthread_local_var->maxpatternsize = maxdistance;

    // 分配距离相关数据结构
    distancelist *Distance = new_distancelist(pthread_local_var, ro_vars);
    int *Tag = newTags(maxdistance / TAGSEP + 1);
    int Toptag = (int)ceil(maxdistance / TAGSEP);
    
    if (!Distance || !Tag) {
        PrintError("Unable to allocate distance structures");
        goto cleanup_alignpair;
    }
    
    pthread_local_var->Distance = Distance;
    pthread_local_var->Tag = Tag;
    pthread_local_var->Toptag = Toptag;

    // 初始化数据结构
    init_links(pthread_local_var);
    init_distanceseenarray(pthread_local_var);
    init_distance(pthread_local_var, ro_vars);

    // 分配统计内存
    int *Statistics_Distance = (int*)calloc(4 * maxdistance, sizeof(int));
    if (!Statistics_Distance) {
        PrintError("Unable to allocate statistics memory");
        goto cleanup_distance;
    }
    pthread_local_var->Statistics_Distance = Statistics_Distance;

    clear_distancelist(Distance, pthread_local_var, ro_vars);

    // 主算法执行
    ret = newtupbo(pthread_local_var, myResult, ro_vars);

    // 结果后处理
    Index_List *GlobalIndexList = ret.GlobalIndexList;
    GlobalIndexList = RemoveBySize(GlobalIndexList, ro_vars->maxperiod);
    GlobalIndexList = SortByIndex(GlobalIndexList);
    GlobalIndexList = RemoveRedundancy(GlobalIndexList);
    GlobalIndexList = SortByCount(GlobalIndexList);
    ret.GlobalIndexList = GlobalIndexList;

cleanup_distance:
    free(Statistics_Distance);
    free(Tag);
    free_distanceseenarray(pthread_local_var);
    free(pthread_local_var->_DistanceEntries);
    free(Distance);

cleanup_alignpair:
    free(textprime);
    free(textsecnd);
    free(indexprime);
    free(indexsecnd);

cleanup_base_ptr:
    // 直接释放分配的内存
    free(base_ptr);

cleanup_S_rows:
    free(S);

cleanup_consensus:
    free(consensus_buffer);

cleanup:
    for (i = 1; i <= ro_vars->NTS; i++) {
        free(pthread_local_var->Tuplehash[i]);
        free(pthread_local_var->History[i]);
    }
    free(pthread_local_var->Sortmultiples);
    free_bestperiodlist(pthread_local_var);

    return ret;
}

#endif
