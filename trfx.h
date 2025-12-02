#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include "bseq.h"

typedef struct
{
    unsigned int maxwraplength; //

    signed char match;
    signed char mismatch;
    unsigned int indel;
    unsigned int PM;
    unsigned int PI;
    unsigned int minscore;
    unsigned int maxperiod;
} trf_opt;

  
typedef struct
{
    unsigned int maxwraplength; //

    int alpha; /* match bonus */ // global read only											// global read only												// global read only
    int beta;                    /* mismatch penalty */
    int delta;                   /* indel penalty */
    int PM;                      // global read only
    int PI;                      // global read only											///(int)ceil(maxdistance/TAGSEP); /* last tag in list */  //gloabl read only constant
    int minscore;
    int maxperiod;
                             // global read only
    //int four_to_the[11];      // = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576}; // global read only
    int Min_Distance_Entries; // = 20;										 // minimum number of places to store a tuple match.  Usually this is the same as the  distance, but for small distances, we allow more tuple matches because we want see a significant number of matches */
    int Min_Distance_Window;  // = 20; /* minimum size of distance window. */ //
    signed char *SM;                  // = NULL;														 // global read only												 // global read only
    unsigned char *waitdata;
    int *sumdata; // first used in main, then pass to TRF, then read only ,  deternmined by user para
    int NTS;      // same as above
    int *Index;   // = NULL;														// glbal read only
    int Tuplesize[10];
} readonly_vars_struct; 


extern int mm_verbose;
extern double mm_realtime0;


#ifdef __cplusplus
extern "C"
{
#endif

    // void trf_opt_init(trf_opt *opt);
    int trf_search_file(const char *fn, const trf_opt *opt, int n_threads, long tbatch_size , readonly_vars_struct *ro_vars);
    void trf_init(trf_opt *opt);
    void init_sm(signed char match, signed char mismatch, readonly_vars_struct *ro_vars);
    void init_and_fill_coin_toss_stats2000_with_4tuplesizes(readonly_vars_struct *ro_vars);
    void init_index(readonly_vars_struct *ro_vars);
  //  void init_readonly(readonly_vars_struct *ro_vars, trf_opt *opt);
    void init_readonly(readonly_vars_struct *ro_vars, const trf_opt *opt);

    double cputime(void);
    double realtime(void);

#ifdef __cplusplus
}
#endif

#endif
