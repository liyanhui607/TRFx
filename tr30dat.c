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

#ifndef TR30DAT_C
#define TR30DAT_C

#include <stdlib.h> /* has calloc definition */
#include <math.h>
#include "tr30dat.h"

void PrintError(char *errortext); // 正确签名

#ifdef __cplusplus
extern "C"
{
#endif
	int has_actual_gpu();
#ifdef __cplusplus
}
#endif

extern int GetTopPeriods_cuda(unsigned char *h_pattern, int length, int *toparray, int Index_cuda[256], int MAXDISTANCE_cuda);

void FreeList(Index_List *headptr);

int d_range(int d, const readonly_vars_struct *restrict ro_vars)
{
	return ((int)floor(2.3 * sqrt((float)ro_vars->PI / 100 * d)));
}

static inline int isPowerOfTwo(int n)
{
	return (n > 0) && ((n & (n - 1)) == 0);
}

/* Jan 27, 2006, Gelfand, changed to use Similarity Matrix to avoid N matching itself */
/* This function may be called multiple times (for different match/mismatch scores) */

/**************************  init_index()  *************************/

/*******************************************************************/

void init_index(readonly_vars_struct *restrict ro_vars)

{
	ro_vars->Index = (int *)calloc(256, sizeof(int));
	if (ro_vars->Index == NULL)
	{
		trf_message("\nInit_index: Out of memory!");
		exit(-1);
	}
	ro_vars->Index['A'] = 0;
	ro_vars->Index['C'] = 1;
	ro_vars->Index['G'] = 2;
	ro_vars->Index['T'] = 3;
}

// 2 4 8 16 32 48 64 80 96 112 128 144 160

void resize_S_matrix(thread_local_var_struct *restrict pthread_local_var, const int required_cols)
{
	pthread_local_var->allocated_S_cols = required_cols;
	const int rows = pthread_local_var->maxwraplength + 1;
	const size_t aligned_row_bytes = (required_cols * sizeof(int)); // ((required_cols * sizeof(int) + 3 & ~3)); //

	char *base_ptr = (char *)pthread_local_var->S[0];
	int **S = pthread_local_var->S;

	for (int i = 0; i < rows; i++)
	{
		S[i] = (int *)(base_ptr + i * aligned_row_bytes);
	}
}

/***************************** narrowbandwrap()  macros **************************/

/* *pcurr>=maxscore added 2.17.05 gary benson -- to extend alignment as far as possible */
#define test_trace_and_backwards_maxscore                          \
	if ((realr <= start - size) && (*pcurr == 0))                  \
		pleft = (*pcurr = -1000);                                  \
	else                                                           \
		end_of_trace = FALSE;                                      \
	if (*pcurr >= maxscore)                                        \
	{                                                              \
		maxscore = *pcurr;                                         \
		minrealrow = realr;                                        \
		/*mincol = c; */                                           \
		/***6/9/05 G. Benson***/ mincolbandcenter = Bandcenter[r]; \
		/*** 6/14/05 G. Benson ***/ mincolposition = i;            \
	}

#define test_maxrowscore_with_match                                 \
	if (*pcurr > maxrowscore)                                       \
	{                                                               \
		maxrowscore = *pcurr;                                       \
		if ((*pcurr == *pdiag) && (match_yes_no == ro_vars->alpha)) \
			matchatmax_col = c;                                     \
		else                                                        \
			matchatmax_col = -2;                                    \
	}

#define test_maxrowscore_without_match \
	if (*pcurr > maxrowscore)          \
		matchatmax_col = -2;

/* *pcurr>=maxscore added 2.17.05 gary benson -- to extend alignment as far as possible */
#define test_trace_and_forward_maxscore    \
	if ((realr >= start) && (*pcurr == 0)) \
		pleft = (*pcurr = -1000);          \
	else                                   \
		end_of_trace = FALSE;              \
	if (*pcurr >= maxscore)                \
	{                                      \
		maxscore = *pcurr;                 \
		maxrealrow = realr;                \
		maxrow = r;                        \
		maxcol = c;                        \
	}





__attribute__((hot, aligned(64))) void narrowbandwrap(const int start, const int size, const  int bandradius, const int bandradiusforward,
													  int option, const int tuplesize, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int g;
	register int *pup, *pdiag, *pcurr, pleft;
	int c, realr, end_of_trace,
		maxscore = 0, maxrow = 0, maxcol = 0; // mincol = 0;
	int maxrealrow = 0, minrealrow = 0;
	//unsigned char currchar;
	int matches_in_diagonal, matchatmax_col, i, k, maxrowscore,
		lastmatchatmax_col, match_yes_no;
	int mincolbandcenter = 0, zeroat, mincolposition = 0;
	unsigned int r;

	const int required_cols = max(2 * bandradius + 1, 2 * bandradiusforward + 1); // +5作为安全边界

	// fprintf(stderr, "required_cols: %d\n", required_cols);
	// FILE *log_file = fopen("pacbio_required_cols.txt", "a");
    // fprintf(log_file, "%d\n", required_cols);
    // fclose(log_file);


	if (pthread_local_var->allocated_S_cols < required_cols)
	{
		resize_S_matrix(pthread_local_var, required_cols);
	}

	unsigned char *EC = pthread_local_var->EC;
	const int delta = ro_vars->delta;					// 缓存delta值
	int *const restrict Diag = pthread_local_var->Diag; // Diag数组指针
	int *const restrict Up = pthread_local_var->Up;		// Up数组指针
	int *Bandcenter = pthread_local_var->Bandcenter;
	int **S = pthread_local_var->S; // __attribute__((aligned(64)))
	unsigned char *Sequence = pthread_local_var->Sequence;
	int Rows = pthread_local_var->Rows;
	const int match = ro_vars->alpha;
	const int mismatch = ro_vars->beta;
	 int w = bandradius;
	 int w2 = 2 * w;
	const int Length = pthread_local_var->Length;

	// fprintf(stderr, "start: %d, size: %d, bandradius: %d, bandradiusforward: %d \n",  start, size, bandradius, bandradiusforward);

	if (MAXBANDWIDTH < w2 + 1)
	{
		trf_message("\nIn narrowbandwrap, MAXBANDWIDTH: %d exceeded by 2*w+1: %d\n",
					MAXBANDWIDTH, w2 + 1);
		exit(-1);
	}

	/* fill EC */
	if (option == WITHCONSENSUS)
		for (g = 0; g < size; g++)
			EC[g] = pthread_local_var->Consensus.pattern[g];
	else
	{
		const unsigned char *src = &Sequence[start - size + 1];
		memcpy(EC, src, size);
	}

	/* backward wdp */
	// int row = pthread_local_var->maxwraplength;
	maxscore = 0;
	realr = start + 1;
	r = pthread_local_var->maxwraplength;
	Bandcenter[r] = 0;
	matches_in_diagonal = 0;
	matchatmax_col = -2;

	int *S_curr = S[r];
	pcurr = &S_curr[0];
	pdiag = &Diag[0];
	pup = &Up[0];

	//__builtin_prefetch(&Sequence[realr + 16], 0, 0);

	/* 3/14/05 gary benson -- reverse direction */
	/* change zeroth row values to put in -1000 in unreachable cells
	and gap penalty in cells beyond start location */

	for (i = 0; i <= w; i++)
	{
		*pup = (*pdiag = (*pcurr = 0 + delta * (w - i))) + delta;
		pup++;
		pdiag++;
		pcurr++;
	}

	for (i = w + 1; i <= w2; i++)
	{
		*pup = (*pdiag = (*pcurr = -1000)) + delta;
		pup++;
		pdiag++;
		pcurr++;
	}

	/* backwards */

	end_of_trace = FALSE;

	while ((!end_of_trace) && (realr > 1) && (r > 0))
	{
		r--;
		realr--;
		Rows++;
		maxrowscore = -1;
		lastmatchatmax_col = matchatmax_col;
		end_of_trace = TRUE;
		const unsigned char currchar = Sequence[realr];
		S_curr = S[r];
		pcurr = &S_curr[w2];
		pleft = -1000;
       // int a = 0;
		if (matches_in_diagonal >= tuplesize)
		{
			int new_col = matchatmax_col - 1;
			if (new_col < 0)
				new_col += size;
			Bandcenter[r] = new_col;
			//int a = Bandcenter[r];
			//a = new_col;
		}
		else
		{
			int new_center = Bandcenter[r + 1] - 1;
			if (new_center < 0)
				new_center += size;
			Bandcenter[r] = new_center;
			//a = new_center;
		}

	//	{
			const int a = Bandcenter[r];
			const int b = Bandcenter[r + 1];
			k = a - b;
			k += (k >> (32 - 1)) & size; // 无分支版本 32 = (sizeof(int)<<3)
		//}

		if (size - k <= k)
			k = -(size - k);
		//c = (Bandcenter[r] + w) % size;
		c = (a + w) % size;

		if (k <= -1) /* band shifts left */
		{
			k = -k;
			pdiag = &Diag[w2 - k + 1];
			pup = &Up[w2 - k];
			for (i = w2; i >= k; i--)
			{

				*pdiag += (match_yes_no = (currchar == EC[c]) ? match : mismatch);
				pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_with_match;
				pcurr--;
				pdiag--;
				pup--;
				// c=(c-1+size)%size;
				c--;
				if (c < 0)
					c += size; // 避免除法，仅用条件判断和加法
			}
			i = k - 1;
			*pdiag += (match_yes_no = (currchar == EC[c]) ? match : mismatch);
			pleft = (*pcurr = max3(0, *pdiag, pleft)) + delta;
			test_trace_and_backwards_maxscore;
			test_maxrowscore_with_match;
			pcurr--;
			// c=(c-1+size)%size;
			c--;
			if (c < 0)
				c += size; // 避免除法，仅用条件判断和加法

			// 优化后实现（无分支）
			// c = (c > 0) ? c - 1 : size - 1;

			for (i = k - 2; i >= 0; i--)
			{
				pleft = (*pcurr = max2(0, pleft)) + delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_without_match;
				pcurr--;
				// c=(c-1+size)%size;
				c--;
				if (c < 0)
					c += size; // 避免除法，仅用条件判断和加法
			}
		}
		else /* band shifts right */
		{
			pdiag = &Diag[w2];
			pup = &Up[w2];
			for (i = w2; i >= w2 - k + 1; i--)
			{
				pleft = (*pcurr = max2(0, pleft)) + delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_without_match;
				pcurr--;
				// c=(c-1+size)%size;
				c--;
				if (c < 0)
					c += size; // 避免除法，仅用条件判断和加法
			}
			i = w2 - k;
			pleft = (*pcurr = max3(0, *pup, pleft)) + delta;
			test_trace_and_backwards_maxscore;
			test_maxrowscore_without_match;
			pcurr--;
			pup--;

			// c=(c-1+size)%size;
			c--;
			if (c < 0)
				c += size; // 避免除法，仅用条件判断和加法
			for (i = w2 - k - 1; i >= 0; i--)
			{
				*pdiag += (match_yes_no = (currchar == EC[c]) ? match : mismatch);
				pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta;
				test_trace_and_backwards_maxscore;
				test_maxrowscore_with_match;
				pcurr--;
				pdiag--;
				pup--;

				// c=(c-1+size)%size;
				c--;
				if (c < 0)
					c += size; // 避免除法，仅用条件判断和加法
			}
		}
		pcurr = &S_curr[0];
		pdiag = &Diag[0];
		pup = &Up[0];
		for (i = 0; i <= w2; i++)
		{
			*pup = (*pdiag = *pcurr) + delta;
			pcurr++;
			pdiag++;
			pup++;
		}
		if ((matchatmax_col - lastmatchatmax_col + size) % size == size - 1)
			matches_in_diagonal++;
		else
			matches_in_diagonal = 0;
	}

	/**********************************************************************/
	/* forward */

	{

		r = 0;
		realr = minrealrow - 1;
		// Bandcenter[0] = (mincolbandcenter - 1 + size) % size;
		{ // 优化版本（条件判断代替模运算）
			int new_center = mincolbandcenter - 1;
			if (new_center < 0)
				new_center += size;
			Bandcenter[0] = new_center;
		}

		S_curr = S[r];
		matches_in_diagonal = 0;
		matchatmax_col = -2;
		maxscore = 0;
		pup = &Up[0];
		pdiag = &Diag[0];
		pcurr = &S_curr[0];
		zeroat = mincolposition - w;
		w = bandradiusforward;
		w2 = 2 * w;
		if (MAXBANDWIDTH < w2 + 1)
		{
			trf_message("\nIn narrowbandwrap, MAXBANDWIDTH: %d exceeded by bandradiusforward 2*w+1: %d\n",
						MAXBANDWIDTH, w2 + 1);
			exit(-1);
		}

		for (i = 0; i < w + zeroat; i++)
		{
			*pup = (*pdiag = (*pcurr = -1000)) + delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		for (i = w + zeroat; i <= w2; i++)
		{
			*pup = (*pdiag = (*pcurr = 0 + delta * (i - (w + zeroat)))) + delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		/* compute until end of trace */
		end_of_trace = FALSE;
		while ((!end_of_trace) && (realr < Length) && (r < pthread_local_var->maxwraplength))
		{
			r++;
			realr++;
			Rows++;
			end_of_trace = TRUE;
			maxrowscore = -1;
			lastmatchatmax_col = matchatmax_col;
			const unsigned char  currchar = Sequence[realr];

			S_curr = S[r];
			pcurr = &S_curr[0];
			pleft = -1000; /* don't use pleft for first entry */

			if (matches_in_diagonal >= tuplesize)
			{
				int new_val = matchatmax_col + 1;
				Bandcenter[r] = (new_val == size) ? 0 : new_val;
			}
			else
			{

				int prev = Bandcenter[r - 1]; // 缓存前值减少内存访问
				int next = prev + 1;
				Bandcenter[r] = (next >= size) ? 0 : next;
			}

			// k = (Bandcenter[r] - Bandcenter[r - 1] + size) % size;
			//{									  // 优化版本（无分支位运算）
				const int current = Bandcenter[r];	  // 缓存当前值减少内存访问
				const int previous = Bandcenter[r - 1]; // 缓存前值
				const int diff = current - previous;
				const int sign_mask = diff >> 31;	   //(sizeof(int) * 8 - 1); // 负数时为全1，非负为0
				k = diff + (sign_mask & size); // 等价于 diff <0 ? diff+size : diff
			//}

			if (size - k <= k)
				k = -(size - k);
			if (k >= 1) /* band shifts right */
			{
				pdiag = &Diag[k - 1];
				pup = &Up[k];
				//c = (Bandcenter[r] - w + size) % size;
				c = (current - w + size) % size;

				for (i = 0; i <= w2 - k; i++)
				{
					*pdiag += (match_yes_no = (currchar == EC[c]) ? match : mismatch);
					pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_with_match;
					pcurr++;
					pdiag++;
					pup++;

					// c=(c+1)%size;
					c++;
					if (c >= size)
						c -= size;
				}
				i = w2 - k + 1;
				*pdiag += (match_yes_no = (currchar == EC[c]) ? match : mismatch); // match_ab(currchar, EC[c], match,  mismatch));
				pleft = (*pcurr = max3(0, *pdiag, pleft)) + delta;
				test_trace_and_forward_maxscore;
				test_maxrowscore_with_match;
				pcurr++;

				// c=(c+1)%size;
				c++;
				if (c >= size)
					c -= size;

				for (i = w2 - k + 2; i <= w2; i++)
				{
					pleft = (*pcurr = max2(0, pleft)) + delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_without_match;
					pcurr++;

					// c=(c+1)%size;
					c++;
					if (c >= size)
						c -= size;
				}
			}
			else /* band shifts left */
			{
				k = -k;
				//c = (Bandcenter[r] - w + size) % size;
				c = (current - w + size) % size;

				pup = &Up[0];
				pdiag = &Diag[0];
				for (i = 0; i <= k - 1; i++)
				{
					pleft = (*pcurr = max2(0, pleft)) + delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_without_match;
					pcurr++;

					// c=(c+1)%size;
					c++;
					if (c >= size)
						c -= size;
				}
				i = k;
				pleft = (*pcurr = max3(0, *pup, pleft)) + delta;
				test_trace_and_forward_maxscore;
				test_maxrowscore_without_match;
				pcurr++;
				pup++;

				// c=(c+1)%size;
				c++;
				if (c >= size)
					c -= size;

				for (i = k + 1; i <= w2; i++)
				{
					*pdiag += (match_yes_no = (currchar == EC[c]) ? match : mismatch);
					pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta;
					test_trace_and_forward_maxscore;
					test_maxrowscore_with_match;
					pcurr++;
					pdiag++;
					pup++;
					// c=(c+1)%size;
					c++;
					if (c >= size)
						c -= size;
				}
			}
			S_curr = S[r];
			pcurr = &S_curr[0];
			pdiag = &Diag[0];
			pup = &Up[0];
			for (i = 0; i <= w2; i++)
			{
				*pup = (*pdiag = *pcurr) + delta;
				pcurr++;
				pdiag++;
				pup++;
			}

			if ((matchatmax_col - lastmatchatmax_col + size) % size == 1)
				matches_in_diagonal++;
			else
				matches_in_diagonal = 0;
		}
		if (maxscore >= 0)
		{
			pthread_local_var->Maxrealrow = maxrealrow;
			pthread_local_var->Maxrow = maxrow;
			pthread_local_var->Maxcol = maxcol;
			pthread_local_var->Maxscore = maxscore;
		}
	}
}

/*******************************************************************/

void newwrap(const int start, const int size, int consensuspresent, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)

{

	int g;
	register int *pup, *pdiag, *pcurr, pleft;
	int adjlength, adjmone, c, realr, end_of_trace,
		maxscore, minrow = 0, maxrow = 0, maxcol = 0, modstart, maxrealrow = 0;
	char currchar;
	/* Feb 16, 2016 Yozen */
	unsigned int r = 0;

	if (pthread_local_var->allocated_S_cols < size)
	{
		resize_S_matrix(pthread_local_var, size);
	}

	unsigned char *EC = pthread_local_var->EC;
	const int delta = ro_vars->delta;		   // 缓存delta值
	int *const Diag = pthread_local_var->Diag; // Diag数组指针
	int *const Up = pthread_local_var->Up;	   // Up数组指针
	int **S = pthread_local_var->S;
	unsigned char *Sequence = pthread_local_var->Sequence;
	int Rows = pthread_local_var->Rows;
	const int match = ro_vars->alpha;
	const int mismatch = ro_vars->beta;
	// int* const S_row = pthread_local_var->S[r]; // 当前行S[r]的指针
	const int Length = pthread_local_var->Length;

	const int maxwraplength = pthread_local_var->maxwraplength;

	/* fill pthread_local_var->EC */
	if (consensuspresent)
		for (g = 0; g < size; g++)
			EC[g] = pthread_local_var->Consensus.pattern[g];
	else
	{
		unsigned const char *src = &Sequence[start - size + 1];
		memcpy(EC, src, size * sizeof(char));
	}

	/* backward wdp */
	maxscore = 0;
	realr = start + 1;
	r = maxwraplength;

	adjlength = size - 1;
	adjmone = adjlength - 1;
	pup = Up;
	pdiag = Diag;
	int *S_curr = S[r];
	pcurr = &S_curr[0];

	for (c = 0; c < size; c++)
	{
		*pup = (*pdiag = (*pcurr = maxscore)) + delta;
		pup++;
		pdiag++;
		pcurr++;
	}
	/* backwards */
	end_of_trace = FALSE;
	while ((!end_of_trace) && (realr > 1) && (r > 0))
	{
		r--;
		realr--;
		Rows++;
		currchar = Sequence[realr];
		S_curr = S[r];
		pcurr = &S_curr[adjlength];
		pleft = delta; /* first pass. set S[r][0]=0 */
		pdiag = &Diag[adjlength];
		pup = &Up[adjlength];
		for (c = adjlength; c >= 0; c--)
		{
			//*pdiag += match_ab(currchar, pthread_local_var->EC[c]);
			*pdiag += (currchar == EC[c]) ? match : mismatch;

			pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta;
			pcurr--;
			pdiag--;
			pup--;
		}

		end_of_trace = TRUE; /* setup for encountering a break in the trace */
		/*second pass */
		pcurr = &S_curr[adjlength];
		/* pleft set from first pass */
		pdiag = &Diag[adjlength];
		pup = &Up[adjlength];
		for (c = adjlength; c >= 0; c--)
		{
			*pup = (pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta);
			if ((realr <= start - max(size, ro_vars->Min_Distance_Window)) && (*pcurr == 0))
				/* if (*pcurr==0) doesn't work with consensus */
				*pup = pleft = *pcurr = -1000;
			else
				end_of_trace = FALSE;

			/* 3/14/05 gary benson -- made >= to extend alignments as far as possible */
			if (*pcurr >= maxscore) /* test for maximum */
			{
				maxscore = *pcurr;
				minrow = realr;
			}
			pcurr--;
			pdiag--;
			pup--;
		}
		pcurr = &S_curr[adjlength];
		pdiag = &Diag[adjmone];
		for (c = adjmone; c >= 0; c--)
		{
			*pdiag = *pcurr;
			pcurr--;
			pdiag--;
		}
		Diag[adjlength] = *pcurr;
	}
	/***********************************/
	r = 0;
	realr = start - 1;

	pup = Up;
	S_curr = S[r];
	pdiag = &Diag[0];
	pcurr = &S_curr[0];
	/* initialize_scoring_array top row (*pcurr)*/
	/* initialize diagonal branch (*pdiag) and up branch (*pup) */
	for (c = 0; c < size; c++)
	{
		*pup = (*pdiag = (*pcurr = maxscore)) + delta;
		pup++;
		pdiag++;
		pcurr++;
	}

	end_of_trace = FALSE;
	while ((!end_of_trace) && (realr < Length) && (r < maxwraplength))
	{
		r++;
		realr++;
		Rows++;
		currchar = Sequence[realr];
		S_curr = S[r];

		pcurr = &S_curr[0];
		pleft = delta; /* first pass. S[r][adjlength]=0  */
		pdiag = &Diag[0];
		pup = &Up[0];
		for (c = 0; c < size; c++)
		{
			//*pdiag += match_ab(currchar, pthread_local_var->EC[c]);
			*pdiag += (currchar == EC[c]) ? match : mismatch;
			pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta;
			pcurr++;
			pdiag++;
			pup++;
		}
		end_of_trace = TRUE; /* setup for encountering a break in the trace */
		/*second pass */
		pcurr = &S_curr[0];
		/* pleft set from first pass */
		pdiag = &Diag[0];
		pup = &Up[0];
		for (c = 0; c < size; c++)
		{
			*pup = (pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta);
			if ((realr >= start) && (*pcurr == 0))
				*pup = pleft = *pcurr = -1000;
			else
				end_of_trace = FALSE;

			if (*pcurr >= maxscore) /* test for maximum */
			{
				maxscore = *pcurr;
				maxrow = realr; /* ??? do we need this? */
			}
			pcurr++;
			pdiag++;
			pup++;
		}
		pcurr = &S_curr[0];
		pdiag = &Diag[1];
		for (c = 0; c < size; c++)
		{
			*pdiag = *pcurr;
			pcurr++;
			pdiag++;
		}
		Diag[0] = Diag[size];
	}

	if (maxscore >= 0)
	{
		r = 0;
		modstart = start + 1;
		realr = minrow - size - 1; /* go back at least one additional pattern */
		if (realr < 0)
			realr = 0; /* length in case of misalignment in */
		/* backwards scan due to starting with */
		/* maxscore in every column */
		maxscore = 0;
		/****/
		pup = Up;
		S_curr = S[r];
		pdiag = &Diag[0];
		pcurr = &S_curr[0];

		/* initialize_scoring_array top row (*pcurr)*/
		/* initialize diagonal branch (*pdiag) and up branch (*pup) */
		for (c = 0; c < size; c++)
		{
			*pup = (*pdiag = (*pcurr = maxscore)) + delta;
			pup++;
			pdiag++;
			pcurr++;
		}

		end_of_trace = FALSE;
		while ((!end_of_trace) && (realr < Length) && (r < maxwraplength))
		{
			r++;
			realr++;
			Rows++;
			currchar = Sequence[realr];
			S_curr = S[r];
			pcurr = &S_curr[0];
			pleft = delta; /* first pass. S[r][adjlength]=0  */
			pdiag = &Diag[0];
			pup = &Up[0];
			for (c = 0; c < size; c++)
			{
				*pdiag += (currchar == EC[c]) ? match : mismatch;
				pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta;
				pcurr++;
				pdiag++;
				pup++;
			}
			end_of_trace = TRUE; /* setup for encountering a break in the trace */
			/*second pass */
			pcurr = &S_curr[0];
			pdiag = &Diag[0];
			pup = &Up[0];
			for (c = 0; c < size; c++)
			{
				*pup = (pleft = (*pcurr = max4(0, *pdiag, *pup, pleft)) + delta);
				if ((realr >= modstart) && (*pcurr == 0))
				{
					*pup = pleft = *pcurr = -1000;
				}
				else
					end_of_trace = FALSE;

				if (*pcurr > maxscore) /* test for maximum */
				{
					maxscore = *pcurr;
					maxrealrow = realr;
					maxrow = r;
					maxcol = c;
				}
				pcurr++;
				pdiag++;
				pup++;
			}
			pcurr = &S_curr[0];
			pdiag = &Diag[1];
			for (c = 0; c < size; c++)
			{
				*pdiag = *pcurr;
				pcurr++;
				pdiag++;
			}
			Diag[0] = Diag[size];
		}
	}
	if (maxscore >= 0)
	{
		pthread_local_var->Maxrealrow = maxrealrow;
		pthread_local_var->Maxrow = maxrow;
		pthread_local_var->Maxcol = maxcol;
		pthread_local_var->Maxscore = maxscore;
	}
}

/*******************************************************************/

void init_bestperiodlist(thread_local_var_struct *pthread_local_var)
{
	pthread_local_var->Bestperiodlist->next = NULL;
}

/**********************************************************/
/**********************************************************/

void free_bestperiodlist(thread_local_var_struct *pthread_local_var)
{
	bestperiodlistelement *entry, *entrylast;
	entry = pthread_local_var->Bestperiodlist->next;
	pthread_local_var->Bestperiodlist->next = NULL;
	while (entry != NULL)
	{
		entrylast = entry;
		entry = entry->next;
		free(entrylast);
	}
}

/**********************************************************/
/**********************************************************/

/*** modified 6/2/05 G. Benson ***/

void add_to_bestperiodlist(const int d, thread_local_var_struct *pthread_local_var)
{
	if (d == 1)
		return; // 提前返回简化逻辑

	int *Sortmultiples = pthread_local_var->Sortmultiples;
	pairalign *pAlignPair = &pthread_local_var->AlignPair;

	bestperiodlistelement *ptr = (bestperiodlistelement *)calloc(1, sizeof(bestperiodlistelement));
	if (ptr == NULL)
	{
		trf_message("\nAdd_to_bestperiodlist: Out of memory!");
		exit(-1);
	}

	ptr->indexlow = pAlignPair->indexprime[pAlignPair->length];
	ptr->indexhigh = pAlignPair->indexprime[1];
	memcpy(&ptr->best1, Sortmultiples, 5 * sizeof(int)); // 连续复制5个整型值
	ptr->next = pthread_local_var->Bestperiodlist->next;
	pthread_local_var->Bestperiodlist->next = ptr;
}

/*** modified 6/2/05 G. Benson ***/
void adjust_bestperiod_entry(const int d, thread_local_var_struct *pthread_local_var)
/* intens length of best period entry if the consensus alignment
   turns out to be inter than the first alignment; similar
   to distanceseen */
{
	bestperiodlistelement *ptr;
	pairalign *pAlignPair = &pthread_local_var->AlignPair;

	/*** added 6/2/05 G. Benson***/
	/* when distance is 1, no pthread_local_var->Sortmultiples are defined due to change in multiples_criteria_4*/
	if (d != 1)
	{
		ptr = pthread_local_var->Bestperiodlist->next;
		if (ptr->indexhigh > pAlignPair->indexprime[1])
			ptr->indexhigh = pAlignPair->indexprime[1];
	}
}
/**********************************************************/
/**********************************************************/

int search_for_range_in_bestperiodlist(const int start, const int distance, thread_local_var_struct *pthread_local_var)
{
	bestperiodlistelement *entry = pthread_local_var->Bestperiodlist->next;
	bestperiodlistelement *entrylast = pthread_local_var->Bestperiodlist;
	int range_covered = FALSE;

	// 预计算常用值
	const int maxdistance = pthread_local_var->maxdistance;
	const int maxdistance_threshold = start - 2 * maxdistance;
	const int waiting_time_criteria = pthread_local_var->Distance[distance].waiting_time_criteria;
	const int low_threshold = start - 2 * distance + 1 + waiting_time_criteria;

	while (entry != NULL)
	{
		bestperiodlistelement *next_entry = entry->next; // 预先获取下一个节点

		if (entry->indexhigh < maxdistance_threshold)
		{
			// 移除太远的条目
			entrylast->next = next_entry;
			free(entry);
			entry = next_entry;
			continue;
		}

		// 检查是否覆盖指定范围
		if (entry->indexlow <= low_threshold && entry->indexhigh >= start)
		{
			range_covered = TRUE;

			// 检查距离是否匹配最佳周期
			if (entry->best1 == distance || entry->best2 == distance ||
				entry->best3 == distance || entry->best4 == distance ||
				entry->best5 == distance)
			{
				return TRUE;
			}
		}

		// 移动到下一个条目
		entrylast = entry;
		entry = next_entry;
	}

	return !range_covered;
}
/*******************************************************************/
/*******************************************************************/
void init_distanceseenarray(thread_local_var_struct *pthread_local_var)
/* created 5/23/05 G. Benson */
{
	pthread_local_var->Distanceseenarray = (distanceseenarrayelement *)calloc(MAXDISTANCECONSTANT + 1, sizeof(distanceseenarrayelement));
	if (pthread_local_var->Distanceseenarray == NULL)
	{
		trf_message("\nInit pthread_local_var->Distanceseenarray: Out of memory!");
		exit(-1);
	}
}
/*******************************************************************/
/*******************************************************************/
void free_distanceseenarray(thread_local_var_struct *pthread_local_var)
/* created 5/23/05 G. Benson */
{
	free(pthread_local_var->Distanceseenarray);
}
/*******************************************************************/
/*******************************************************************/
void add_to_distanceseenarray(const int location, const int distance, const int end, const int score, thread_local_var_struct *pthread_local_var)
/* created 5/23/05 G. Benson */
/* adds the extent of an alignment (end) at the given distance
   location and score are not currently used */
{
	distanceseenarrayelement *ptr;

	ptr = &pthread_local_var->Distanceseenarray[distance];
	ptr->index = location;
	ptr->end = end;
	ptr->score = score;
}
/*******************************************************************/
/*******************************************************************/
int search_for_distance_match_in_distanceseenarray(const int distance, const int start, thread_local_var_struct *pthread_local_var)
/* created 5/23/05 G. Benson */
/* searches for an alignment with patternsize of distance in the region
   including start.  True means found and alignment should be blocked.*/
{
	distanceseenarrayelement ptr;
	ptr = pthread_local_var->Distanceseenarray[distance];

	if (ptr.end >= start)
		return (TRUE);
	else
		return (FALSE);
}

/*******************************************************************/

void init_distanceseenlist(thread_local_var_struct *pthread_local_var)
{
	pthread_local_var->Distanceseenlist->next = NULL;
}

/*******************************************************************/

void free_distanceseenlist(thread_local_var_struct *pthread_local_var)
{
	distancelistelement *entry, *entrylast;
	entry = pthread_local_var->Distanceseenlist->next;
	pthread_local_var->Distanceseenlist->next = NULL;
	while (entry != NULL)
	{
		entrylast = entry;
		entry = entry->next;
		free(entrylast);
	}
}

/*******************************************************************/

void add_to_distanceseenlist(const int location, const int distance, const int end, const int score,
	const int acceptstatus, thread_local_var_struct *pthread_local_var)
{
	distancelistelement *ptr = (distancelistelement *)calloc(1, sizeof(distancelistelement));
	if (ptr == NULL)
	{
		trf_message("\nAdd_to_distanceseenlist: Out of memory!");
		exit(-1);
	}

	distancelistelement *entry = pthread_local_var->Distanceseenlist->next;

	// 检查是否需要修改现有条目
	if (acceptstatus == WITHCONSENSUS && entry != NULL &&
		location == entry->index && end <= entry->end &&
		(distance <= 250 || distance == entry->distance))
	{
		entry->end = 0; // 标记为需要移除
	}

	// 设置新条目的字段
	ptr->index = location;
	ptr->distance = distance;
	ptr->changed_from_distance = (acceptstatus == WITHCONSENSUS && entry != NULL) ? entry->distance : distance;
	ptr->end = end;
	ptr->score = score;
	ptr->accepted = acceptstatus;

	// 插入到链表头部
	ptr->next = pthread_local_var->Distanceseenlist->next;
	pthread_local_var->Distanceseenlist->next = ptr;
}

int search_for_distance_match_in_distanceseenlist(const int distance, const int start,
												  thread_local_var_struct *pthread_local_var,
												  const readonly_vars_struct *restrict ro_vars)
{
	distancelistelement *entry = pthread_local_var->Distanceseenlist->next;
	distancelistelement *entrylast = pthread_local_var->Distanceseenlist;

	// 预计算阈值
	const int threshold = min(start - distance,
							  max(start - ro_vars->minscore / ro_vars->alpha,
								  start - ro_vars->Min_Distance_Window)) +
						  1;

	while (entry != NULL)
	{
		distancelistelement *next_entry = entry->next; // 预先获取下一个节点

		if (entry->end < threshold)
		{
			// 移除过期条目
			entrylast->next = next_entry;
			free(entry);
			entry = next_entry;
			continue;
		}

		// 计算距离差
		int absdiff = abs(distance - entry->distance);

		// 检查是否匹配
		if (absdiff == 0 || (absdiff <= 5 && distance == entry->changed_from_distance))
		{
			return TRUE;
		}

		// 移动到下一个条目
		entrylast = entry;
		entry = next_entry;
	}

	return FALSE;
}





// 辅助函数：调整索引
static inline int adjust_index(int idx, int size, int size_half) {
    if (idx < 0) idx += size;
    if (idx > size_half) idx -= size;
    return idx;
}

void get_narrowband_pair_alignment_with_copynumber(const int size,
    const int bandradius, int option,
    thread_local_var_struct *pthread_local_var,
    const readonly_vars_struct *restrict ro_vars) {
    
    // 预计算常量（保持不变）
    const int delta = ro_vars->delta;
    const int match = ro_vars->alpha;
    const int mismatch = ro_vars->beta;
    const int w = bandradius;
    const int w2 = 2 * w;
    const int size_half = size / 2;

    // 局部指针缓存
    unsigned char *x = pthread_local_var->Sequence;
    unsigned char *y = pthread_local_var->EC;
    int *Bandcenter = pthread_local_var->Bandcenter;
    int **S = pthread_local_var->S;
    pairalign *pAlignPair = &pthread_local_var->AlignPair;

    // 初始化变量
    int realr = pthread_local_var->Maxrealrow;
    int r = pthread_local_var->Maxrow;
    int c = pthread_local_var->Maxcol;
    int length = 0;
    double copy_count = 0.0;
    pthread_local_var->Copynumber = 0;

    // 预计算 fullcopy
    const int fullcopy = (pthread_local_var->Maxcol + 1) % size;

    // k 和 i 初始值
    int k = adjust_index(c - Bandcenter[r], size, size_half);
    int i = w + k;
    pAlignPair->score = S[r][i];

    // 主回溯循环
    while (r > 0 || (option == GLOBAL && c != -1)) {
        int *S_curr = S[r];
        int *S_prev = (r > 0) ? S[r - 1] : NULL;

        // 提前终止：LOCAL 模式下 S == 0 且非合法路径
        if (option == LOCAL && S_curr[i] == 0) {
            // 检查是否是合法的 0（来自匹配/插入/删除）
            int prev_k = adjust_index(c - Bandcenter[r - 1], size, size_half);
            int upi = w + prev_k;
            int band_shift_val = Bandcenter[r] - Bandcenter[r - 1];
            int band_shift = adjust_index(band_shift_val, size, size_half);
            
            int match_cost = (x[realr] == y[c]) ? match : mismatch;
            int valid = 0;

            if (band_shift >= 1) {
                if (i <= w2 - band_shift) {
                    valid = (S_curr[i] == S_prev[upi - 1] + match_cost) ||
                            (S_curr[i] == S_prev[upi] + delta) ||
                            (S_curr[i] == S_curr[i - 1] + delta);
                } else if (i == w2 - band_shift + 1) {
                    valid = (S_curr[i] == S_prev[upi - 1] + match_cost) ||
                            (S_curr[i] == S_curr[i - 1] + delta);
                } else {
                    valid = (S_curr[i] == S_curr[i - 1] + delta);
                }
            } else {
                band_shift = -band_shift;
                if (i <= band_shift - 1) {
                    valid = (S_curr[i] == S_curr[i - 1] + delta);
                } else if (i == band_shift) {
                    valid = S_prev && ((S_curr[i] == S_prev[upi] + delta) ||
                                       (S_curr[i] == S_curr[i - 1] + delta));
                } else {
                    valid = S_prev && ((S_curr[i] == S_prev[upi - 1] + match_cost) ||
                                       (S_curr[i] == S_prev[upi] + delta) ||
                                       (S_curr[i] == S_curr[i - 1] + delta));
                }
            }
            
            if (!valid) break;
        }

        // GLOBAL 模式 c == -1 特殊处理
        if (option == GLOBAL && c == -1) {
            int prev_k = adjust_index(c - Bandcenter[r - 1], size, size_half);
            int upi = w + prev_k;
            
            if (S_curr[i] == S_prev[upi] + delta) {
                length++;
                fill_align_pair(x[realr], '-', length, realr, c + 1 + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                realr--;
                r--;
                i = upi;
                continue;
            } else {
                trf_message("\nget_pair_alignment_with_copynumber: error in trace back");
                trf_message("\nfailed to go up when c=-1");
                break;
            }
        }

        // 计算上一行位置
        int prev_k = adjust_index(c - Bandcenter[r - 1], size, size_half);
        int upi = w + prev_k;
        int band_shift_val = Bandcenter[r] - Bandcenter[r - 1];
        int band_shift = adjust_index(band_shift_val, size, size_half);
        int matched = (x[realr] == y[c]) ? match : mismatch;

        // 根据 band_shift 处理
        int found = 0;
        if (band_shift >= 1) {
            if (i <= w2 - band_shift) {
                if (S_curr[i] == S_prev[upi - 1] + matched) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair(x[realr], y[c], length, realr, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair(x[realr], y[c], length, realr, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    realr--;
                    r--;
                    i = upi - 1;
                    found = 1;
                } else if (S_curr[i] == S_prev[upi] + delta) {
                    length++;
                    if (option == LOCAL) {
                        fill_align_pair(x[realr], '-', length, realr, (c + 1) % size);
                    } else {
                        fill_align_pair(x[realr], '-', length, realr, c + 1 + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                    }
                    realr--;
                    r--;
                    i = upi;
                    found = 1;
                } else if (i > 0 && S_curr[i] == S_curr[i - 1] + delta) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair('-', y[c], length, realr + 1, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair('-', y[c], length, realr + 1, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    i--;
                    found = 1;
                }
            } else if (i == w2 - band_shift + 1) {
                if (S_curr[i] == S_prev[upi - 1] + matched) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair(x[realr], y[c], length, realr, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair(x[realr], y[c], length, realr, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    realr--;
                    r--;
                    i = upi - 1;
                    found = 1;
                } else if (i > 0 && S_curr[i] == S_curr[i - 1] + delta) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair('-', y[c], length, realr + 1, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair('-', y[c], length, realr + 1, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    i--;
                    found = 1;
                }
            } else {
                if (i > 0 && S_curr[i] == S_curr[i - 1] + delta) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair('-', y[c], length, realr + 1, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair('-', y[c], length, realr + 1, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    i--;
                    found = 1;
                }
            }
        } else {
            band_shift = -band_shift;
            if (i <= band_shift - 1) {
                if (i > 0 && S_curr[i] == S_curr[i - 1] + delta) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair('-', y[c], length, realr + 1, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair('-', y[c], length, realr + 1, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    i--;
                    found = 1;
                }
            } else if (i == band_shift) {
                if (S_curr[i] == S_prev[upi] + delta) {
                    length++;
                    if (option == LOCAL) {
                        fill_align_pair(x[realr], '-', length, realr, (c + 1) % size);
                    } else {
                        fill_align_pair(x[realr], '-', length, realr, c + 1 + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                    }
                    realr--;
                    r--;
                    i = upi;
                    found = 1;
                } else if (i > 0 && S_curr[i] == S_curr[i - 1] + delta) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair('-', y[c], length, realr + 1, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair('-', y[c], length, realr + 1, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    i--;
                    found = 1;
                }
            } else {
                if (S_curr[i] == S_prev[upi - 1] + matched) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair(x[realr], y[c], length, realr, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair(x[realr], y[c], length, realr, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    realr--;
                    r--;
                    i = upi - 1;
                    found = 1;
                } else if (S_curr[i] == S_prev[upi] + delta) {
                    length++;
                    if (option == LOCAL) {
                        fill_align_pair(x[realr], '-', length, realr, (c + 1) % size);
                    } else {
                        fill_align_pair(x[realr], '-', length, realr, c + 1 + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                    }
                    realr--;
                    r--;
                    i = upi;
                    found = 1;
                } else if (i > 0 && S_curr[i] == S_curr[i - 1] + delta) {
                    length++;
                    if (c == fullcopy) copy_count += 1.0;
                    if (option == LOCAL) {
                        fill_align_pair('-', y[c], length, realr + 1, c);
                        c--;
                        if (c < 0) c += size;
                    } else {
                        fill_align_pair('-', y[c], length, realr + 1, c + pthread_local_var->Maxrealcol - pthread_local_var->Maxcol);
                        c--;
                    }
                    i--;
                    found = 1;
                }
            }
        }

        if (!found) {
            trf_message("\nget_pair_alignment_with_copynumber: error in trace back");
            break;
        }
    }

    // 设置最终结果
    pAlignPair->length = length;
    int col_diff = pthread_local_var->Maxcol - c;
    if (col_diff < 0) col_diff += size;
    pthread_local_var->Copynumber = copy_count + (double)col_diff / size;
}




void get_pair_alignment_with_copynumber(const int size, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int i, si, j, length;
	unsigned char *x, *y;
	int **S = pthread_local_var->S;
	const int delta = ro_vars->delta;
	pairalign *pAlignPair = &pthread_local_var->AlignPair;

	// 初始化主变量
	x = pthread_local_var->Sequence;
	y = pthread_local_var->EC;
	i = pthread_local_var->Maxrealrow;
	si = pthread_local_var->Maxrow;
	j = pthread_local_var->Maxcol;

	const int size_minus_1 = size - 1;
	const int fullcopy = (j + 1) % size; // 固定值：完整拷贝边界
	pAlignPair->score = S[si][j];
	length = 0;
	pthread_local_var->Copynumber = 0.0;

	const int match = ro_vars->alpha;
	const int mismatch = ro_vars->beta;

	while (1)
	{
		int *S_curr = S[si];

		// 终止条件：当前得分 <= 0
		if (S_curr[j] <= 0)
		{
			pAlignPair->length = length;

			int maxcol = pthread_local_var->Maxcol;
			if (maxcol >= j)
			{
				pthread_local_var->Copynumber += (double)(maxcol - j) / size;
			}
			else
			{
				pthread_local_var->Copynumber += (double)(maxcol + size - j) / size;
			}
			return;
		}

		// 预计算常用值
		int prev_j = (j + size_minus_1) % size; // j-1 mod size
		int *S_prev = (si > 0) ? S[si - 1] : NULL;

		// 当前字符和匹配得分（仅在需要时计算）
		unsigned char current_x = x[i];
		unsigned char current_y = y[j];
		int match_score = (current_x == current_y) ? match : mismatch;

		// 分支1: 匹配/错配（对角线）
		if (S_prev != NULL && S_curr[j] == S_prev[prev_j] + match_score)
		{
			length++;
			if (j == fullcopy)
			{
				pthread_local_var->Copynumber += 1.0;
			}
			fill_align_pair(current_x, current_y, length, i, j);
			i--;
			si--;
			j = prev_j;
		}
		// 分支2: 删除（来自上方）
		else if (S_prev != NULL && S_curr[j] == S_prev[j] + delta)
		{
			length++;
			fill_align_pair(current_x, '-', length, i, (j + 1) % size); // 注意：此处 j 不变，上一行同列
			i--;
			si--;
			// j 不变
		}
		// 分支3: 插入（来自左方）
		else if (S_curr[j] == S_curr[prev_j] + delta)
		{
			length++;
			if (j == fullcopy)
			{
				pthread_local_var->Copynumber += 1.0;
			}
			fill_align_pair('-', current_y, length, i + 1, j);
			j = prev_j;
			// i 不变, si 不变
		}
		// 错误路径
		else
		{
			trf_message("\nget_pair_alignment_with_copynumber: error in trace back");
			trf_message("\nrow=%d  column=%d", i, j);
			trf_message("\nS=%d  Sleft=%d  Sup=%d  Sdiag=%d  match_ab=%d",
						S_curr[j],
						S_curr[prev_j],
						(S_prev != NULL) ? S_prev[j] : -1,
						(S_prev != NULL) ? S_prev[prev_j] : -1,
						match_score);
			break; // 跳出循环，但未设置 length？原函数如此，保持一致
		}
	}
}

/*******************************************************************/

void reverse(thread_local_var_struct *pthread_local_var)
{
	pairalign *pAlignPair = &pthread_local_var->AlignPair;
	const int ml = pAlignPair->length;
	const int half_len = ml / 2;

	// 提前取出指针，避免每次访问结构体
	char *textprime = pAlignPair->textprime;
	char *textsecnd = pAlignPair->textsecnd;
	int *indexprime = pAlignPair->indexprime;
	int *indexsecnd = pAlignPair->indexsecnd;

	for (int j = 1; j <= half_len; j++)
	{
		int opposite_j = ml - j + 1;

		// 交换 textprime[j] 和 textprime[opposite_j]
		char temp_c = textprime[j];
		textprime[j] = textprime[opposite_j];
		textprime[opposite_j] = temp_c;

		// 交换 textsecnd[j] 和 textsecnd[opposite_j]
		temp_c = textsecnd[j];
		textsecnd[j] = textsecnd[opposite_j];
		textsecnd[opposite_j] = temp_c;

		// 交换 indexprime[j] 和 indexprime[opposite_j]
		int temp_i = indexprime[j];
		indexprime[j] = indexprime[opposite_j];
		indexprime[opposite_j] = temp_i;

		// 交换 indexsecnd[j] 和 indexsecnd[opposite_j]
		temp_i = indexsecnd[j];
		indexsecnd[j] = indexsecnd[opposite_j];
		indexsecnd[opposite_j] = temp_i;
	}
}

static inline void find_max_base(int a, int c, int g, int t, int *max, char *maxchar)
{
	*max = a;
	*maxchar = 'A';
	if (c > a)
	{
		*max = c;
		*maxchar = 'C';
	}
	if (g > *max)
	{
		*max = g;
		*maxchar = 'G';
	}
	if (t > *max)
	{
		*max = t;
		*maxchar = 'T';
	}
}

void get_consensus(int patternlength, thread_local_var_struct *pthread_local_var)
{
	pairalign *pAlignPair = &pthread_local_var->AlignPair;
	cons_data *Consensus = &pthread_local_var->Consensus;
	const int array_size = 2 * MAXPATTERNSIZECONSTANT + 1;

	// 直接引用数组（不再加 restrict）
	int *A_arr = Consensus->A;
	int *C_arr = Consensus->C;
	int *G_arr = Consensus->G;
	int *T_arr = Consensus->T;
	int *dash_arr = Consensus->dash;
	int *insert_arr = Consensus->insert;
	int *total_arr = Consensus->total;
	unsigned char *pattern_arr = Consensus->pattern;

	// 批量清零
	const size_t bytes = array_size * sizeof(int);
	memset(A_arr, 0, bytes);
	memset(C_arr, 0, bytes);
	memset(G_arr, 0, bytes);
	memset(T_arr, 0, bytes);
	memset(dash_arr, 0, bytes);
	memset(insert_arr, 0, bytes);
	memset(total_arr, 0, bytes);

	// 初始化 pattern 为 DASH
	for (int c = 0; c < array_size; c++)
	{
		pattern_arr[c] = DASH;
	}

	if (pAlignPair->length == 0)
	{
		pthread_local_var->ConsClasslength = 0;
		return;
	}

	int lastindex = -1;
	int i = 1;
	const int len = pAlignPair->length;

	// 单独处理第一个元素
	int idx_secnd = pAlignPair->indexsecnd[i];
	unsigned char text_prime = pAlignPair->textprime[i];
	int pos = 2 * idx_secnd + 1;

	switch (text_prime)
	{
	case 'A':
		A_arr[pos]++;
		break;
	case 'C':
		C_arr[pos]++;
		break;
	case 'G':
		G_arr[pos]++;
		break;
	case 'T':
		T_arr[pos]++;
		break;
	default:
		dash_arr[pos]++;
		break;
	}
	lastindex = idx_secnd;
	i++;

	// 主循环
	while (i <= len)
	{
		idx_secnd = pAlignPair->indexsecnd[i];
		text_prime = pAlignPair->textprime[i];

		if (idx_secnd != lastindex)
		{
			pos = 2 * idx_secnd + 1;
			switch (text_prime)
			{
			case 'A':
				A_arr[pos]++;
				break;
			case 'C':
				C_arr[pos]++;
				break;
			case 'G':
				G_arr[pos]++;
				break;
			case 'T':
				T_arr[pos]++;
				break;
			default:
				dash_arr[pos]++;
				break;
			}

			int total_idx = (idx_secnd == patternlength - 1) ? 0 : (2 * idx_secnd + 2);
			total_arr[total_idx]++;
			lastindex = idx_secnd;
			i++;
		}
		else
		{
			int base_idx = 2 * idx_secnd;
			insert_arr[base_idx]++;

			// ✅ 修复：使用 while 而非 do-while，避免越界
			while (i <= len && pAlignPair->indexsecnd[i] == lastindex)
			{
				text_prime = pAlignPair->textprime[i];
				switch (text_prime)
				{
				case 'A':
					A_arr[base_idx]++;
					break;
				case 'C':
					C_arr[base_idx]++;
					break;
				case 'G':
					G_arr[base_idx]++;
					break;
				case 'T':
					T_arr[base_idx]++;
					break;
				default:
					dash_arr[base_idx]++;
					break;
				}
				i++;
			}
		}
	}

	// 构建共识序列
	int max;
	char maxchar;

	// 奇数索引
	for (int k = 1; k <= 2 * patternlength; k += 2)
	{
		max = dash_arr[k];
		maxchar = DASH;
		if (A_arr[k] > max)
		{
			max = A_arr[k];
			maxchar = 'A';
		}
		if (C_arr[k] > max)
		{
			max = C_arr[k];
			maxchar = 'C';
		}
		if (G_arr[k] > max)
		{
			max = G_arr[k];
			maxchar = 'G';
		}
		if (T_arr[k] > max)
		{
			max = T_arr[k];
			maxchar = 'T';
		}
		pattern_arr[k] = maxchar;
	}

	// 偶数索引
	for (int k = 0; k <= 2 * patternlength; k += 2)
	{
		if (total_arr[k] > 0 && 2 * insert_arr[k] >= total_arr[k])
		{
			find_max_base(A_arr[k], C_arr[k], G_arr[k], T_arr[k], &max, &maxchar);
			pattern_arr[k] = maxchar;
		}
		else
		{
			pattern_arr[k] = DASH;
		}
	}

	// 压缩
	int j = 0;
	for (int k = 0; k <= 2 * patternlength; k++)
	{
		if (pattern_arr[k] != DASH)
		{
			pattern_arr[j++] = pattern_arr[k];
		}
	}
	pthread_local_var->ConsClasslength = j;
}

distancelist *new_distancelist(thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int g, N, K;
	distanceentry *ptr;
	const int maxdistance = pthread_local_var->maxdistance;
	const int Min_Distance_Entries = ro_vars->Min_Distance_Entries;
	distancelist *objptr = (distancelist *)calloc(maxdistance + 1, sizeof(distancelist));
	K = Min_Distance_Entries + 1;
	N = maxdistance + 1;
	ptr = pthread_local_var->_DistanceEntries = (distanceentry *)malloc(((K + N) * (N - K + 1) / 2 + K * (K - 1)) * sizeof(distanceentry));

	for (g = 1; g <= maxdistance; g++)
	{
		int count = max(g, Min_Distance_Entries) + 1;
		objptr[g].entry = ptr;
		memset(objptr[g].entry, 0, count * sizeof(distanceentry));
		ptr += count;
	}

	return (objptr);
}

void clear_distancelist(distancelist *objptr, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int g;
	const int maxdistance = pthread_local_var->maxdistance;
	const int min_distance_entries = ro_vars->Min_Distance_Entries;

	for (g = 1; g <= maxdistance; g++)
	{
		objptr[g].lowindex = 0;
		objptr[g].highindex = max(g, min_distance_entries);
		objptr[g].numentries = 0;
		objptr[g].nummatches = 0;
	}
}

void init_links(thread_local_var_struct *pthread_local_var)
{
	pthread_local_var->Distance[0].linkup = pthread_local_var->maxdistance + 1;
}

void add_tuple_match_to_Distance_entry(int location, int size, int d, distancelist *objptr, const readonly_vars_struct *restrict ro_vars)
{
	// 获取只读参数
	const int min_dist_win = ro_vars->Min_Distance_Window;
	const int min_dist_ent = ro_vars->Min_Distance_Entries;

	// 计算缓冲区大小：max(d, min_dist_ent) + 1
	const int windowsize = (d > min_dist_ent) ? d + 1 : min_dist_ent + 1;

	// 获取目标桶
	distancelist *list = &objptr[d];
	distanceentry *entries = list->entry;
	int *lo = &list->lowindex;
	int *hi = &list->highindex;
	int *z = &list->numentries;
	int *m = &list->nummatches;

	// === 阶段1：清理过期条目（滑动窗口左侧） ===
	if (*z > 0)
	{
		// 正确计算：min(location - d, location - min_dist_win) + 1
		int effective_offset = (d > min_dist_win) ? d : min_dist_win;
		int windowleftend = location - effective_offset + 1;

		int l = *lo;
		int num = *z;
		int matches = *m;

		while (num > 0 && entries[l].location < windowleftend)
		{
			matches -= entries[l].size;
			num--;
			l++;
			if (l >= windowsize)
				l = 0;
		}

		*lo = l;
		*z = num;
		*m = matches;
	}

	// === 阶段2：尝试合并到最后一个条目 ===
	if (*z > 0 && entries[*hi].location == location - 1)
	{
		entries[*hi].location = location; // 移动到当前 location
		entries[*hi].size++;			  // ← 只增加 1
		(*m)++;
		return;
	}

	// === 阶段3：插入新条目 ===
	int new_hi = *hi + 1;
	if (new_hi >= windowsize)
	{
		new_hi = 0;
	}

	entries[new_hi].location = location;
	entries[new_hi].size = size; // 使用传入的 size

	*hi = new_hi;
	(*z)++;
	*m += size;
}



void link_Distance_window(int d, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
    int t, f, h;
    int *Tag = pthread_local_var->Tag;                      // 缓存 Tag 数组指针
    int Toptag = pthread_local_var->Toptag;                // 缓存 Toptag 值
    distancelist *Distance = pthread_local_var->Distance;  // 缓存 Distance 数组指针

    // 计算 t 的索引
    t = (int)ceil(d / TAGSEP);

    // 获取当前 Tag[t] 的值
    int current_tag = Tag[t];

    // 处理 Tag[t] < d 的情况
    if (current_tag < d)
    {
        f = current_tag;
        while (t <= Toptag && Tag[t] < d)
        {
            Tag[t] = d;
            t++;
            current_tag = Tag[t];  // 更新 current_tag 避免重复访问数组
        }
    }
    // 处理 Tag[t] > d 的情况
    else if (current_tag > d)
    {
        f = current_tag;
        while (f > d)
        {
            f = Distance[f].linkdown;
        }
        if (f == d)
        {
            trf_message("\nTag error following links.  f==d=%d", d);
        }
    }
    // 处理 Tag[t] == d 的情况
    else
    {
        trf_message("\nTag error Tag[%d]=%d", t, d);
        exit(-2);
    }

    // 使用临时变量减少结构体访问
    distancelist *d_entry = &Distance[d];
    distancelist *f_entry = &Distance[f];

    // 设置 linkdown 和 linkup
    d_entry->linkdown = f;
    d_entry->linkup = f_entry->linkup;
    f_entry->linkup = d;

    // 更新 h 的 linkdown（如果存在）
    h = d_entry->linkup;
    if (h <= pthread_local_var->maxdistance)
    {
        Distance[h].linkdown = d;
    }
    d_entry->linked = TRUE;
}


void untag_Distance_window(int d, int linkdown, thread_local_var_struct *pthread_local_var)
{
	int t;
	int *Tag = pthread_local_var->Tag;		// Cache Tag array pointer
	int Toptag = pthread_local_var->Toptag; // Cache Toptag value

	/* get next highest tag */
	t = (int)ceil(d / TAGSEP);
	if (Tag[t] != d)
		return; /* Tag[t] is the largest index less or */
				/* equal to (t)x(TAGSEP) that is linked */
	else
	{
		while ((t <= Toptag) && (Tag[t] == d)) /* check higher tags and replace d */
											   /* with its linkdown */
		{
			Tag[t] = linkdown;
			t++;
		}
	}
}

int no_matches_so_unlink_Distance(int d, int location,
								  distancelist *objptr, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	const int min_dist_win = ro_vars->Min_Distance_Window;
	const int min_dist_ent = ro_vars->Min_Distance_Entries;
	const int maxdistance = pthread_local_var->maxdistance;

	distancelist *list = &objptr[d];
	distancelist *Distance = pthread_local_var->Distance;

	distanceentry *entries = list->entry;
	int *lo = &list->lowindex;
	int *z = &list->numentries;
	int *m = &list->nummatches;

	if (*z == 0)
		goto do_unlink;

	int effective_d = (d > min_dist_win) ? d : min_dist_win;
	int windowleftend = location - effective_d + 1;
	int win_size = (d > min_dist_ent) ? d + 1 : min_dist_ent + 1;

	int num_remaining = *z;
	int low_idx = *lo;

	while (num_remaining > 0 && entries[low_idx].location < windowleftend)
	{
		*m -= entries[low_idx].size;
		num_remaining--;
		if (++low_idx >= win_size)
			low_idx = 0;
	}

	*z = num_remaining;
	*lo = low_idx;

	if (*z == 0)
	{
	do_unlink:
	{
		int g = Distance[d].linkdown;
		int h = Distance[d].linkup;

		Distance[g].linkup = h;
		if (h <= maxdistance)
			Distance[h].linkdown = g;
		Distance[d].linked = FALSE;

		untag_Distance_window(d, g, pthread_local_var);
		return 1;
	}
	}

	return 0;
}

int GetTopPeriods(unsigned char *pattern, int length, int *toparray, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int topind;
	double topval;
	int heads[16];
	int *history;
	double* counts;
	int i,t,end,tupid;
	int curr,dist;
	double n,xysum,xsum,ysum,x2sum,s;

	/* allocate an array of counts */
	counts = (double*) calloc(length,sizeof(double));
	if(counts==NULL) return 1;

	/* allocate history array */
	history = (int*) malloc(length*sizeof(int));
	if(history==NULL)
	{
		free(counts);
		return 1;
	}

	/* clear the heads array which point into history array */
	for(i=0;i<16;i++) heads[i]=-1;

	/* scan pattern for tuples of size 2 */
	for(i=0,end=length-2;i<=end;i++)
	{
		/* figure out tuple id */
		tupid = ro_vars->Index[pattern[i]]*4+ro_vars->Index[pattern[i+1]];

		/* record last occurence into history and update heads[] pointer */
		history[i] = heads[tupid];
		heads[tupid]=i;

		/* loop into history and add distances */
		/* 11/17/15 G. Benson */
		/* limit maximum length of distance recorded between tuples to MAXDISTANCECONSTANT*3 = 6,000*/
		/* this should be long enough to deter finding periods that are not the most frequent */
		/* Without this change, this procudure is quadratic in the length, which could be several million */
		/* and caused the program to hang with long centromeric repeats */
		/* for(curr=i;history[curr]!=-1;curr=history[curr])*/
		dist = 0;
		for(curr=i;((history[curr]!=-1)&&(dist<(MAXDISTANCECONSTANT*3)));curr=history[curr])
		{
			dist = i-history[curr];
			counts[dist]+=1.0;
		}
	}

	/* compute slope using least-square regression */
	xysum=xsum=ysum=x2sum=0.0;
	end = length-2;
	for(i=1;i<=end;i++)
	{
		xysum += (i*counts[i]);
		xsum += (i);
		ysum += (counts[i]);
		x2sum += (i*i);
	}
	n = end;
	s = (n*xysum-xsum*ysum)/(n*x2sum-xsum*xsum);

	/* flatten trend by adding -s per increment */
	end = length-2;
	for(i=1;i<=end;i++)
	{
		counts[i] = counts[i] - i*s;
	}

	/* pick highest values */
	end = length-2;
	if(end>pthread_local_var->maxdistance) end = pthread_local_var->maxdistance; /* 3/14/05 accepts smaller multiples is best ones are too large */
	for(t=0;t<NUMBER_OF_PERIODS;t++)
	{
		/* do t passes to find t highes counts */
		topind=0;
		topval=0.0;
		for(i=1;i<=end;i++)
		{
			if(counts[i]>topval)
			{
				topind=i;
				topval = counts[i];
			}
		}

		/* copy to array passed as parameter */
		toparray[t] = topind;
		counts[topind]=0.0;
	}

	/* free memory */
	free(counts);
	free(history);

	return 0;
}

// this funtion
int GetTopPeriods2(unsigned char *pattern, int length, int *toparray, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int topind;
	double topval;
	int num[16];
	double *counts; //=  pthread_local_var->counts;
	int i, t, end, tupid;
	int dist;
	double n, xysum, xsum, ysum, x2sum, s;

	/* allocate an array of counts */
	counts = (double *)calloc(length, sizeof(double));
	if (counts == NULL)
		return 1;

	int (*history)[length] = (int (*)[length])malloc(16 * sizeof(int[length]));
	if (history)
		memset(history, 0, 16 * sizeof(int[length]));
	else
	{
		free(counts);
		return 1;
	}

	/* clear the heads array which point into history array */
	for (i = 0; i < 16; i++)
		num[i] = 0;

	/* scan pattern for tuples of size 2 */
	end = length - 2;
	for (i = 0; i <= end; i++)
	{
		tupid = (ro_vars->Index[pattern[i]] << 2) | ro_vars->Index[pattern[i + 1]];
		int *curr_history = history[tupid]; // 缓存行指针
		int curr_num = num[tupid];			// 缓存当前计数
		curr_history[curr_num] = i;			// 写入数据
		num[tupid] = curr_num + 1;			// 更新计数（注意顺序）
		dist = 0;
		for (int j = curr_num; j >= 0 && dist < (MAXDISTANCECONSTANT * 3); j--)
		{
			dist = i - curr_history[j]; // 直接使用缓存的指针
			counts[dist]++;
		}
	}

	/* compute slope using least-square regression */
	xysum = xsum = ysum = x2sum = 0.0;
	end = length - 2;
	for (i = 1; i <= end; i++)
	{
		xysum += (i * counts[i]);
		xsum += (i);
		ysum += (counts[i]);
		x2sum += (i * i);
	}
	n = end;
	s = (n * xysum - xsum * ysum) / (n * x2sum - xsum * xsum);

	// printf("alg2: length:%d, s:%lf, pattern_len: %d, pattern: %s\n", length, s, strlen(pattern), pattern);

	end = length - 2;
	for (i = 1; i <= end; i++)
	{
		counts[i] = counts[i] - i * s;
	}

	/* pick highest values */
	end = length - 2;
	if (end > pthread_local_var->maxdistance)
		end = pthread_local_var->maxdistance; /* 3/14/05 accepts smaller multiples is best ones are too large */
	for (t = 0; t < NUMBER_OF_PERIODS; t++)
	{
		/* do t passes to find t highes counts */
		topind = 0;
		topval = 0.0;
		for (i = 1; i <= end; i++)
		{
			if (counts[i] > topval)
			{
				topind = i;
				topval = counts[i];
			}
		}
		/* copy to array passed as parameter */
		toparray[t] = topind; // here the ind, is infact the distance between the same nucletide
							  // printf("alg2: topind:%d,counts[%d]:%lf\n", topind, topind, counts[topind]);

		counts[topind] = 0.0;
	}
	/* free memory */
	free(counts);
	free(history);
	return 0;
}

int multiples_criteria_4(int found_d, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int g;
	int topperiods[NUMBER_OF_PERIODS];
	pairalign *pAlignPair = &pthread_local_var->AlignPair;

	// 缓存结构体成员和常量
	unsigned char *Sequence = pthread_local_var->Sequence;
	int *Index = ro_vars->Index;
	const int maxdistance = pthread_local_var->maxdistance;

	// 计算通用的 lowerindex 和 upperindex
	int lowerindex = pAlignPair->indexprime[pAlignPair->length];
	int upperindex = pAlignPair->indexprime[1];
	unsigned char *pattern = Sequence + lowerindex;
	const int length = upperindex - lowerindex + 1;

	// 处理 size == 1 的情况
	if (found_d == 1)
	{
		int comps[4] = {0};
		int total = length;
		int maxind = 0;

		for (g = 0; g < length; g++)
		{
			int idx = Index[pattern[g]]; // 预计算 Index[pattern[g]]
			comps[idx]++;
		}

		for (g = 1; g < 4; g++)
		{
			if (comps[g] > comps[maxind])
				maxind = g;
		}

		float percentmatch = (float)comps[maxind] * 100.0f / total;
		return (percentmatch >= 80.0f) ? TRUE : FALSE;
	}

// GPU 加速
#if defined(HAVE_GPU)
	if (has_actual_gpu() && length >= 4000)
	{
		if (GetTopPeriods_cuda(pattern, length, topperiods, Index, maxdistance))
		{
			fprintf(stderr, "\nGPU加速失败，回退到CPU版本...");
		}
		else
		{
			goto process_results;
		}
	}
#endif

	// CPU 版本（主流程或回退）
	if (GetTopPeriods2(pattern, length, topperiods, pthread_local_var, ro_vars))
	{
		fprintf(stderr, "\n无法分配计数数组!");
		exit(-1);
	}

#if defined(HAVE_GPU)
process_results:
#endif

	// 将 topperiods 复制到 Sortmultiples
	// for (g = 0; g < NUMBER_OF_PERIODS_INTO_SORTMULTIPLES; g++)
	// 	pthread_local_var->Sortmultiples[g] = topperiods[g];

		memcpy(pthread_local_var->Sortmultiples, topperiods, 
			NUMBER_OF_PERIODS_INTO_SORTMULTIPLES * sizeof(int));
	// 检查 found_d 是否在 topperiods 中
	for (g = 0; g < NUMBER_OF_PERIODS_TO_TEST; g++)
		if (found_d == topperiods[g])
			return TRUE;

	return FALSE;
}

int new_meet_criteria_3(int d, int location, int tuplesize, thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	distancelist *Distance = pthread_local_var->Distance;
	distancelist *main_d_info = &Distance[d];

	const int min_krun_matches = main_d_info->k_run_sums_criteria;
	const int min_distance_window = ro_vars->Min_Distance_Window;
	const int max_d_min_window = (d > min_distance_window) ? d : min_distance_window;
	const int max_first_match_location = (location - max_d_min_window > 0 ? location - max_d_min_window : 0) + main_d_info->waiting_time_criteria;
	const int low_end_of_range = main_d_info->lo_d_range;
	const int high_end_of_range = main_d_info->hi_d_range;
	const int main_d_matches = main_d_info->nummatches;

	// Precompute threshold for waiting time test
	const int range_d_min_for_waiting_time_test = (int)(0.35 * (main_d_matches < min_krun_matches ? main_d_matches : min_krun_matches));

	// Compute first match location for main d
	const int main_lowindex = main_d_info->lowindex;
	distanceentry *main_entry = &main_d_info->entry[main_lowindex];
	const int main_d_first_match_location = main_entry->location - main_entry->size + tuplesize;

	int waiting_time_ok = (main_d_first_match_location <= max_first_match_location);
	int waiting_time_d = waiting_time_ok ? d : 0;

	// Early return if main d already satisfies
	if (main_d_matches >= min_krun_matches && waiting_time_ok)
	{
		return TRUE;
	}

	int m = main_d_matches;
	int d_still_best = TRUE;
	int lopointer = d;
	int hipointer = d;

	// Check lower range (d-1, d-2, ...)
	int t = main_d_info->linkdown;
	while (t >= low_end_of_range && d_still_best)
	{
		distancelist *range_d_info = &Distance[t];
		int s = range_d_info->linkdown; // Save next before potential unlink

		if (no_matches_so_unlink_Distance(t, location, Distance, pthread_local_var, ro_vars))
		{
			t = s;
			continue;
		}

		const int range_d_matches = range_d_info->nummatches;
		if (range_d_matches > main_d_matches)
		{
			d_still_best = FALSE;
		}
		else
		{
			lopointer = t;
			m += range_d_matches;

			// Only check waiting time if not already satisfied
			if (!waiting_time_ok && range_d_matches >= range_d_min_for_waiting_time_test)
			{
				const int range_lowindex = range_d_info->lowindex;
				distanceentry *range_entry = &range_d_info->entry[range_lowindex];
				const int range_d_first_match_location = range_entry->location - range_entry->size + tuplesize;
				if (range_d_first_match_location <= max_first_match_location)
				{
					waiting_time_ok = 1;
					waiting_time_d = t;
				}
			}
		}
		t = s;
	}

	if (!d_still_best)
	{
		return FALSE;
	}

	// Check upper range for dominance (only check if main d is still best)
	t = main_d_info->linkup;
	while (t <= high_end_of_range && d_still_best)
	{
		distancelist *range_d_info = &Distance[t];
		int s = range_d_info->linkup;

		if (!no_matches_so_unlink_Distance(t, location, Distance, pthread_local_var, ro_vars))
		{
			if (range_d_info->nummatches > main_d_matches)
			{
				d_still_best = FALSE;
			}
		}
		t = s;
	}

	if (!d_still_best)
	{
		return FALSE;
	}

	// If accumulated matches in lower+main satisfy criteria
	if (m >= min_krun_matches && waiting_time_ok)
	{
		return TRUE;
	}

	// Sliding window over upper range
	t = main_d_info->linkup;
	while (t <= high_end_of_range)
	{
		distancelist *range_d_info = &Distance[t];
		int s = range_d_info->linkup;

		hipointer = t;
		const int range_d_matches = range_d_info->nummatches;
		m += range_d_matches;

		// Slide left pointer: remove distances that are too far behind
		const int window_min_index = hipointer - d + low_end_of_range;
		while (lopointer < window_min_index)
		{
			m -= Distance[lopointer].nummatches;
			if (m <= 0)
			{
				trf_message("\n*** error in meet criteria, m<=0");
				exit(-16);
			}
			if (waiting_time_ok && lopointer == waiting_time_d)
			{
				waiting_time_ok = 0;
				waiting_time_d = 0;
			}
			lopointer = Distance[lopointer].linkup;
		}

		// Check if this range entry can satisfy waiting time
		if (!waiting_time_ok && range_d_matches >= range_d_min_for_waiting_time_test)
		{
			const int range_lowindex = range_d_info->lowindex;
			distanceentry *range_entry = &range_d_info->entry[range_lowindex];
			const int range_d_first_match_location = range_entry->location - range_entry->size + tuplesize;
			if (range_d_first_match_location <= max_first_match_location)
			{
				waiting_time_ok = 1;
				waiting_time_d = t;
			}
		}

		// Check if current window satisfies criteria
		if (m >= min_krun_matches && waiting_time_ok)
		{
			return TRUE;
		}

		t = s;
	}

	return FALSE;
}

void printECtoBuffer(char *trg, int start, int width, thread_local_var_struct *pthread_local_var)
{
	unsigned char *EC = pthread_local_var->EC; // Cache EC array pointer with correct type

	// Copy from start to end of EC array
	const int first_part_length = width - start;
	memcpy(trg, EC + start, first_part_length);

	// Copy from beginning of EC array to start position
	memcpy(trg + first_part_length, EC, start);

	// Add null terminator
	trg[width] = '\0';
	return;
}

Result get_statistics(thread_local_var_struct *pthread_local_var, Result myResult, const readonly_vars_struct *restrict ro_vars, int index)
{
	// 缓存频繁访问的指针和变量
	const int consensussize = pthread_local_var->Classlength;
	pairalign *pAlignPair = &pthread_local_var->AlignPair;
	const int pAlignPair_length = pAlignPair->length;
	int *Statistics_Distance = pthread_local_var->Statistics_Distance;
	int maxdistance_val = pthread_local_var->maxdistance;

	char *textsecnd = pAlignPair->textsecnd;
	char *textprime = pAlignPair->textprime;
	int *indexsecnd = pAlignPair->indexsecnd;
	int *indexprime = pAlignPair->indexprime;

	int match = 0, mismatch = 0, indel = 0;
	int mindistance = consensussize, maxdistance = consensussize;
	int d, best_match_distance = 0, best_match_count = 0;
	int i, count;
	double diversity[4] = {0}, entropy = 0.0;
	int ACGTcount[4] = {0}; // A, C, G, T
	const int startECpos = indexsecnd[pAlignPair_length];

	// 计算Statistics_Distance数组的有效长度
	int row = maxdistance_val + d_range(maxdistance_val, ro_vars);
	if (row <= 0)
	{
		trf_message("Invalid row value in get_statistics2");
		return myResult;
	}
	memset(Statistics_Distance + 1, 0, row * sizeof(int));

	// 初始化lp和rp指针
	int lp = 1, rp = 2;
	while (lp <= pAlignPair_length && textsecnd[lp] == '-')
		lp++;
	if (lp > pAlignPair_length)
	{
		trf_message("Initial left pointer exceeds pAlignPair_length");
		return myResult;
	}
	while (rp <= pAlignPair_length && indexsecnd[rp] != indexsecnd[lp])
		rp++;
	if (rp > pAlignPair_length)
	{
		trf_message("Initial right pointer exceeds pAlignPair_length");
		return myResult;
	}
	while (rp <= pAlignPair_length && textsecnd[rp] == '-')
		rp++;
	if (rp > pAlignPair_length || indexsecnd[lp] != indexsecnd[rp])
	{
		trf_message("Index mismatch in initial pointers");
		return myResult;
	}

	// 主循环处理比对结果
	while (rp <= pAlignPair_length && lp < rp)
	{
		char c_lp = textsecnd[lp], c_rp = textsecnd[rp];
		char p_lp = textprime[lp], p_rp = textprime[rp];

		if (c_lp != '-' && c_rp != '-')
		{
			if (p_lp != '-' && p_rp != '-')
			{
				if (indexsecnd[lp] != indexsecnd[rp])
				{
					trf_message("Index mismatch during match/mismatch check");
					return myResult;
				}
				if (p_lp == p_rp)
				{
					match++;
					d = abs(indexprime[rp] - indexprime[lp]);
					if (d < 4 * maxdistance_val)
					{
						Statistics_Distance[d]++;
						mindistance = (d < mindistance) ? d : mindistance;
						maxdistance = (d > maxdistance) ? d : maxdistance;
					}
				}
				else
				{
					mismatch++;
				}
				lp++, rp++;
			}
			else if (p_lp == '-' && p_rp == '-')
			{
				if (indexsecnd[lp] != indexsecnd[rp])
				{
					trf_message("Index mismatch during do nothing");
					return myResult;
				}
				lp++, rp++;
			}
			else
			{
				if (indexsecnd[lp] != indexsecnd[rp])
				{
					trf_message("Index mismatch during indel");
					return myResult;
				}
				indel++;
				lp++, rp++;
			}
		}
		else if (c_lp == '-' && c_rp == '-')
		{
			if (p_lp == p_rp)
			{
				match++;
				d = abs(indexprime[rp] - indexprime[lp]);
				if (d < 4 * maxdistance_val)
				{
					Statistics_Distance[d]++;
					mindistance = (d < mindistance) ? d : mindistance;
					maxdistance = (d > maxdistance) ? d : maxdistance;
				}
			}
			else
			{
				mismatch++;
			}
			lp++, rp++;
		}
		else if (c_lp == '-')
		{
			indel++;
			lp++;
		}
		else if (c_rp == '-')
		{
			indel++;
			rp++;
		}
	}

	if (lp >= rp)
	{
		trf_message("Left pointer >= right pointer");
		return myResult;
	}

	int total = match + mismatch + indel;

	// 寻找最佳匹配距离
	for (int g = mindistance; g <= maxdistance; g++)
	{
		if (Statistics_Distance[g] > best_match_count)
		{
			best_match_count = Statistics_Distance[g];
			best_match_distance = g;
		}
	}

	// 统计ACGT分布
	for (i = 1; i <= pAlignPair_length; i++)
	{
		char c = textprime[i];
		if (c != '-')
		{
			switch (c)
			{
			case 'A':
				ACGTcount[0]++;
				break;
			case 'C':
				ACGTcount[1]++;
				break;
			case 'G':
				ACGTcount[2]++;
				break;
			case 'T':
				ACGTcount[3]++;
				break;
			}
		}
	}
	count = ACGTcount[0] + ACGTcount[1] + ACGTcount[2] + ACGTcount[3];
	if (count > 0)
	{
		for (int idx = 0; idx < 4; idx++)
		{
			diversity[idx] = (double)ACGTcount[idx] / count;
		}
		const double log2 = log(2.0);
		for (int idx = 0; idx < 4; idx++)
		{
			if (diversity[idx] > 0)
			{
				entropy -= diversity[idx] * log(diversity[idx]) / log2;
			}
		}
	}

	// 创建链表节点
	Index_List *newptr = (Index_List *)aligned_alloc(8, sizeof(Index_List));
	if (!newptr)
	{
		FreeList(myResult.GlobalIndexList);
		myResult.GlobalIndexList = NULL;
		return myResult;
	}
	pthread_local_var->counterInSeq++;

	newptr->count = pthread_local_var->counterInSeq;
	sprintf(newptr->ref, "%d--%d,%d,%3.1f,%d,%d",
			indexprime[pAlignPair_length], indexprime[1], best_match_distance,
			pthread_local_var->Copynumber, consensussize, (int)pthread_local_var->OUTPUTcount);

	newptr->first = indexprime[pAlignPair_length];
	newptr->last = indexprime[1];
	newptr->period = best_match_distance;
	newptr->copies = pthread_local_var->Copynumber;
	newptr->size = consensussize;
	newptr->matches = (int)(100.0 * match / total);
	newptr->indels = (int)(100.0 * indel / total);
	newptr->score = pAlignPair->score;
	newptr->acount = (int)(100.0 * ACGTcount[0] / count);
	newptr->ccount = (int)(100.0 * ACGTcount[1] / count);
	newptr->gcount = (int)(100.0 * ACGTcount[2] / count);
	newptr->tcount = (int)(100.0 * ACGTcount[3] / count);
	newptr->entropy = entropy;

	// 分配模式字符串
	newptr->pattern = (char *)malloc((consensussize + 1) * sizeof(char));
	if (!newptr->pattern)
	{
		free(newptr);
		FreeList(myResult.GlobalIndexList);
		myResult.GlobalIndexList = NULL;
		return myResult;
	}
	printECtoBuffer(newptr->pattern, startECpos, consensussize, pthread_local_var);

	// 无需加锁，因为每个线程的链表是私有的
	if (!myResult.GlobalIndexList)
	{
		myResult.GlobalIndexList = myResult.GlobalIndexListTail = newptr;
	}
	else
	{
		myResult.GlobalIndexListTail->next = newptr;
		myResult.GlobalIndexListTail = newptr;
	}
	newptr->next = NULL;

	return myResult;
}

/*******************************************************************/

void init_and_fill_coin_toss_stats2000_with_4tuplesizes(readonly_vars_struct *restrict ro_vars)
{
	int *Tuplesize = ro_vars->Tuplesize;

	if (ro_vars->PM == 80)
	{
		ro_vars->NTS = 3; /* ro_vars->Tuplesize[pthread_local_var->ro_vars->NTS+1]={0,4,5,7}; */
		Tuplesize[0] = 0;
		Tuplesize[1] = 4;
		Tuplesize[2] = 5;
		Tuplesize[3] = 7;
		ro_vars->waitdata = waitdata80;
		ro_vars->sumdata = sumdata80;
	}
	else if (ro_vars->PM == 75)
	{
		ro_vars->NTS = 4;
		Tuplesize[0] = 0;
		Tuplesize[1] = 3;
		Tuplesize[2] = 4;
		Tuplesize[3] = 5;
		Tuplesize[4] = 7;
		ro_vars->waitdata = waitdata75;
		ro_vars->sumdata = sumdata75;
	}
	else
	{
		exit(-13);
	}
}

void init_distance(thread_local_var_struct *pthread_local_var, const readonly_vars_struct *restrict ro_vars)
{
	int g, d;
	distancelist *Distance = pthread_local_var->Distance;
	const int maxdist = pthread_local_var->maxdistance;
	const int PM = ro_vars->PM; // Cache PM value

	// Initialize Tuplemaxdistance based on PM value
	if (PM == 80)
	{
		pthread_local_var->Tuplemaxdistance[0] = 0;
		pthread_local_var->Tuplemaxdistance[1] = 29;
		pthread_local_var->Tuplemaxdistance[2] = 159;
		pthread_local_var->Tuplemaxdistance[3] = maxdist;
	}
	else if (PM == 75)
	{
		pthread_local_var->Tuplemaxdistance[0] = 0;
		pthread_local_var->Tuplemaxdistance[1] = 29;
		pthread_local_var->Tuplemaxdistance[2] = 43;
		pthread_local_var->Tuplemaxdistance[3] = 159;
		pthread_local_var->Tuplemaxdistance[4] = maxdist;
	}
	else
	{
		exit(-13);
	}

	// Process small distances first without function calls
	int small_end = SMALLDISTANCE < maxdist ? SMALLDISTANCE : maxdist;
	for (g = 1; g <= small_end; g++)
	{
		Distance[g].lo_d_range = g;
		Distance[g].hi_d_range = g;
	}

	// Process larger distances with d_range calculation
	for (g = small_end + 1; g <= maxdist; g++)
	{
		int dr = d_range(g, ro_vars);
		Distance[g].lo_d_range = (g - dr) > 1 ? (g - dr) : 1;
		Distance[g].hi_d_range = (g + dr) < maxdist ? (g + dr) : maxdist;
	}

	// Initialize waiting_time_criteria and k_run_sums_criteria
	int data_end = maxdist < 2000 ? maxdist : 2000;
	for (d = 1; d <= data_end; d++)
	{
		Distance[d].waiting_time_criteria = ro_vars->waitdata[d];
		Distance[d].k_run_sums_criteria = ro_vars->sumdata[d];
	}

	// For distances beyond 2000, use the value at index 2000
	if (maxdist > 2000)
	{
		int wait_val_2000 = ro_vars->waitdata[2000];
		int sum_val_2000 = ro_vars->sumdata[2000];

		for (d = 2001; d <= maxdist; d++)
		{
			Distance[d].waiting_time_criteria = wait_val_2000;
			Distance[d].k_run_sums_criteria = sum_val_2000;
		}
	}
}

Result newtupbo(thread_local_var_struct *restrict pthread_local_var, Result myResult,  readonly_vars_struct *restrict ro_vars)
{
	// 缓存频繁访问的变量
	const int Length = pthread_local_var->Length;
	unsigned char *Sequence = pthread_local_var->Sequence;
	int *Tuplesize = ro_vars->Tuplesize;
	const int NTS = ro_vars->NTS;
	int **Tuplehash = pthread_local_var->Tuplehash;
	int *Tuplecode = pthread_local_var->Tuplecode;
	historyentry **History = pthread_local_var->History;

	// 预计算常量
	const int four_to_the[11] = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576};
	const int mintuplesize = Tuplesize[1];
	const int maxtuplesize = Tuplesize[NTS];

	// 创建ACGT查找表
	static const char ACGT[256] = {
		['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1};

	// 分配对齐的内存
	size_t size = (pthread_local_var->maxwraplength + 1) * sizeof(int);
	int *aligned_bandcenter;
	if (posix_memalign((void **)&aligned_bandcenter, 64, size) != 0)
	{
		char errmsg[255];
		snprintf(errmsg, 255, "Allocation failed for %zu bytes (aligned)", size);
		PrintError(errmsg);
		exit(-1);
	}
	memset(aligned_bandcenter, 0, size);
	pthread_local_var->Bandcenter = aligned_bandcenter;

	// 初始化历史记录和哈希表
	for (int g = 1; g <= NTS; g++)
	{
		const int tuple_size = Tuplesize[g];
		const int array_size = four_to_the[tuple_size];

		Tuplehash[g] = (int *)calloc(array_size, sizeof(int));
		pthread_local_var->Historysize[g] = 2 * (pthread_local_var->Tuplemaxdistance[g] + 1) + 2;

		History[g] = (historyentry *)calloc(pthread_local_var->Historysize[g], sizeof(historyentry));
		pthread_local_var->Nextfreehistoryindex[g] = 1;
	}

	pthread_local_var->Sortmultiples = (int *)calloc(pthread_local_var->maxdistance + 1, sizeof(int));

	int build_entire_code = 1;
	int badcharindex = 0;

	// 主循环优化
	for (int i = 0; i <= Length; i++)
	{
		// 使用查找表替代strchr
		if (i == 0 || !ACGT[(unsigned char)Sequence[i]])
		{
			badcharindex = i;
			build_entire_code = 1;

			// 查找下一个有效的序列段
			int g = 0;
			while (g < mintuplesize && i < Length)
			{
				i++;
				if (!ACGT[(unsigned char)Sequence[i]])
				{
					badcharindex = i;
					g = 0;
				}
				else
				{
					g++;
				}
			}
			if (g < mintuplesize)
				break;
		}

		// 构建tuple代码
		int code;
		if (build_entire_code)
		{
			code = 0;
			for (int g = badcharindex + 1; g <= i; g++)
			{
				code = (code << 2) + ro_vars->Index[(int)Sequence[g]];
			}
			build_entire_code = (i - badcharindex < maxtuplesize);
		}
		else
		{
			code = ((code % four_to_the[Tuplesize[NTS] - 1]) << 2) + ro_vars->Index[(int)Sequence[i]];
		}

		Tuplecode[NTS] = code;
		for (int h = NTS - 1; h >= 1; h--)
		{
			Tuplecode[h] = code % four_to_the[Tuplesize[h]];
		}

		// 处理每个tuple大小
		for (int g = 1; g <= NTS; g++)
		{
			const int tuple_size = Tuplesize[g];
			if (i - badcharindex < tuple_size)
				continue;

			const int tuple_code = Tuplecode[g];
			int *hash_table = Tuplehash[g];
			historyentry *history = History[g];
			const int historysize = pthread_local_var->Historysize[g];
			const int max_dist_g = pthread_local_var->Tuplemaxdistance[g];
			const int max_dist_g_minus = (g > 1) ? pthread_local_var->Tuplemaxdistance[g - 1] : 0;

			int y = hash_table[tuple_code];
			int h_idx = pthread_local_var->Nextfreehistoryindex[g];
			int j = h_idx + 1;

			if (j == historysize)
			{
				j = 1;
			}

			// 清除旧条目
			if (history[j].location != 0 && j == hash_table[history[j].code])
			{
				hash_table[history[j].code] = 0;
			}

			pthread_local_var->Nextfreehistoryindex[g] = j;
			hash_table[tuple_code] = h_idx;

			// 存储新条目
			history[h_idx].location = i;
			history[h_idx].previous = y;
			history[h_idx].code = tuple_code;

			int yy = h_idx;
			while (y != 0)
			{
				const int d = i - history[y].location;
				if (d > max_dist_g)
				{
					history[yy].previous = 0;
					y = 0;
				}
				else
				{
					yy = y;
					y = history[y].previous;

					if (d > max_dist_g_minus)
					{
						// 处理距离匹配
						add_tuple_match_to_Distance_entry(i, tuple_size, d, pthread_local_var->Distance, ro_vars);
						int found = search_for_distance_match_in_distanceseenarray(d, i, pthread_local_var);

						if (!found)
						{
							if (!pthread_local_var->Distance[d].linked)
							{
								link_Distance_window(d, pthread_local_var, ro_vars);
							}

							if (new_meet_criteria_3(d, i, tuple_size, pthread_local_var, ro_vars) &&
								(d <= 250 || search_for_range_in_bestperiodlist(i, d, pthread_local_var)))
							{

								pthread_local_var->Rows = 0;
								if (d <= SMALLDISTANCE)
								{
									newwrap(i, d, WITHOUTCONSENSUS, pthread_local_var, ro_vars);
									get_pair_alignment_with_copynumber(d, pthread_local_var, ro_vars);
								}
								else
								{
									const int dr = d_range(d, ro_vars);
									const int min_band = max(MINBANDRADIUS, dr);
									const int max_band = min(2 * min_band, d / 3);

									narrowbandwrap(i, d, min_band, max_band, WITHOUTCONSENSUS, RECENTERCRITERION, pthread_local_var, ro_vars);
									get_narrowband_pair_alignment_with_copynumber(d, max_band, LOCAL, pthread_local_var, ro_vars);
								}

								add_to_distanceseenarray(i, d, pthread_local_var->Maxrealrow, pthread_local_var->Maxscore, pthread_local_var);

								// 复制数检查
								const float cn = pthread_local_var->Copynumber;
								const int cn_check = (d <= 50 && cn < 1.9f) ||
													 (d > 50 && d <= 100 && cn < 1.9f - 0.002f * (d - 50)) ||
													 (d > 100 && cn < 1.8f);

								if (!cn_check)
								{
									int pass_multiples_test = multiples_criteria_4(d, pthread_local_var, ro_vars);
									add_to_bestperiodlist(d, pthread_local_var);

									if (pass_multiples_test)
									{
										// 获取共识序列
										pthread_local_var->Classlength = d;
										get_consensus(d, pthread_local_var);
										if (pthread_local_var->ConsClasslength != pthread_local_var->Classlength)
										{
											pthread_local_var->Classlength = pthread_local_var->ConsClasslength;
										}

										const int Classlength = pthread_local_var->Classlength;
										pthread_local_var->Rows = 0;

										if (Classlength <= SMALLDISTANCE)
										{
											newwrap(i, Classlength, WITHCONSENSUS, pthread_local_var, ro_vars);
											get_pair_alignment_with_copynumber(Classlength, pthread_local_var, ro_vars);
										}
										else
										{
											const int dr_cl = d_range(Classlength, ro_vars);
											const int min_band_cl = max(MINBANDRADIUS, dr_cl);
											const int max_band_cl = min(2 * min_band_cl, Classlength / 3);

											narrowbandwrap(i, Classlength, min_band_cl, max_band_cl,
														   WITHCONSENSUS, RECENTERCRITERION, pthread_local_var, ro_vars);
											get_narrowband_pair_alignment_with_copynumber(Classlength, max_band_cl, LOCAL, pthread_local_var, ro_vars);
										}

										add_to_distanceseenarray(i, d, pthread_local_var->Maxrealrow, pthread_local_var->Maxscore, pthread_local_var);
										adjust_bestperiod_entry(d, pthread_local_var);

										// 最终复制数检查
										const float cn_final = pthread_local_var->Copynumber;
										const int cn_final_check = (Classlength <= 50 && cn_final < 1.9f) ||
																   (Classlength > 50 && Classlength <= 100 && cn_final < 1.9f - 0.002f * (d - 50)) ||
																   (Classlength > 100 && cn_final < 1.8f);

										if (!cn_final_check)
										{
											if (Classlength >= pthread_local_var->Minsize &&
												pthread_local_var->AlignPair.score >= ro_vars->minscore)
											{
												trf_message("\nFound at i:%d original size:%d final size:%d",
															i, d, Classlength);
												pthread_local_var->OUTPUTcount++;
												myResult = get_statistics(pthread_local_var, myResult, ro_vars, i);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	free(pthread_local_var->Bandcenter);
	return myResult;
}

#endif
