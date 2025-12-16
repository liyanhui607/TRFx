#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "trfx.h"
#include "misc.c"
#include "bseq.h"

#define min(a, b) (((a) <= (b)) ? (a) : (b))
extern  void init_d_index();

#define TRFx_VERSION "TRFx1.0"

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

size_t get_file_size(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
        return st.st_size;
    } else {
        fprintf(stderr, "Warning: Unable to get file size for %s\n", filename);
        return 0;
    }
}

size_t get_available_memory() {
    long pages = sysconf(_SC_PHYS_PAGES);    // 系统总物理页数
    long page_size = sysconf(_SC_PAGE_SIZE); // 系统页大小（字节）
    
    if (pages == -1 || page_size == -1) {
        fprintf(stderr, "Error: Unable to determine system memory\n");
        return 0;
    }
    
    return (size_t)(pages * page_size);  // 总物理内存
}

size_t get_optimal_batch_size(const char *filename) {
    size_t free_mem = get_available_memory();
    size_t max_mem = 1000 * 1024 * 1024 * 1024L; // 1000G
    size_t min_batch = 100 * 1024 *1024; // 最小100M
    size_t safe_batch = (free_mem * 0.5); // 使用50%内存

    // 添加错误检查
    if (free_mem == 0) {
        fprintf(stderr, "Warning: Unable to determine available memory, using default\n");
        safe_batch = max_mem;
    }

    // 获取输入文件大小
    size_t file_size = get_file_size(filename);//+1M
    if (file_size > 0) {
        //fprintf(stderr, "Input file size: %zu bytes (%.2f GB)\n", file_size, file_size / 1024.0 / 1024 / 1024);
    } else {
        fprintf(stderr, "Warning: Unable to determine file size, using max_mem\n");
        file_size = max_mem; // 如果无法获取文件大小，使用max_mem
    }

    // 在 safe_batch, file_size, max_mem 三者中取最小值
    size_t batch_size = safe_batch;
    if (file_size < batch_size) batch_size = file_size;
    if (max_mem < batch_size) batch_size = max_mem;
    
    // 添加下限保护
    if (batch_size < min_batch) {
        batch_size = min_batch;
    }
    
    return batch_size;
}


int main(int argc, char *argv[])
{
    // 初始化默认值和选项
    int n_threads = 3;
    liftrlimit();
    mm_realtime0 = realtime();

    // 使用复合字面量初始化选项
    trf_opt opt = {
        .maxwraplength = 2000000,
        .match = 2,
        .mismatch = -7,
        .indel = -7,
        .PM = 80,
        .PI = 10,
        .minscore = 50,
        .maxperiod = 2000
    };

    // 解析命令行选项
    int c;
    while ((c = getopt(argc, argv, "Vv:Aa:Bb:Dd:Mm:Ii:Ss:Pp:Tt:")) >= 0) {
        switch (c) {
            case 't': n_threads = atoi(optarg); break;
            case 'A': case 'a': opt.match = atoi(optarg); break;
            case 'B': case 'b': opt.mismatch = -atoi(optarg); break;
            case 'D': case 'd': opt.indel = -atoi(optarg); break;
            case 'M': case 'm': opt.PM = atoi(optarg); break;
            case 'I': case 'i': opt.PI = atoi(optarg); break;
            case 'S': case 's': opt.minscore = atoi(optarg); break;
            case 'P': case 'p': opt.maxperiod = atoi(optarg); break;
            case 'T': n_threads = atoi(optarg); break;
            case 'V': case 'v':
                puts(TRFx_VERSION);
                return 0;
        }
    }

    // 检查是否有输入文件
    if (optind >= argc) {
        fprintf(stderr, "Usage: trfx File -a Match -b Mismatch  -d Delta  -m  PM -i  PI  -s Minscore -p Maxperiod -t thread \n");
        fprintf(stderr, "Default: trfx inputFile -a 2 -b 7 -d 7 -m 80 -i 10 -s 50 -p 2000 -t 3\n");
        fprintf(stderr, "Default is OK in most time , So simply use: ./trfx ./inputFile -t 8 (number of threads)\n");
        fprintf(stderr, "Where: (all weights, penalties, and scores are positive)\n");
        fprintf(stderr, "  File = sequences input file\n");
        fprintf(stderr, "  Match  = matching weight [2]\n");
        fprintf(stderr, "  Mismatch  = mismatching penalty [7]\n");
        fprintf(stderr, "  Delta = indel penalty [7]\n");
        fprintf(stderr, "  PM = match probability (whole number) [80]\n");
        fprintf(stderr, "  PI = indel probability (whole number) [10]\n");
        fprintf(stderr, "  Minscore = minimum alignment score to report [50]\n");
        fprintf(stderr, "  MaxPeriod = maximum period size to report [2000]\n");
        fprintf(stderr, "  [options] = one or more of the following:\n");
        fprintf(stderr, "    -t INT     number of threads [%d]\n", n_threads);
        fprintf(stderr, "    -V         show version number\n");
        return 1;
    }

    // 获取第一个输入文件名并计算最优批处理大小
    const char *input_filename = argv[optind];
    size_t tbatch_size = get_optimal_batch_size(input_filename);
    fprintf(stderr, "batch_size:  (%.2f GB)\n", tbatch_size / 1024.0 / 1024 / 1024);

    // 初始化只读变量
    readonly_vars_struct ro_vars;
    init_readonly(&ro_vars, &opt);
    init_and_fill_coin_toss_stats2000_with_4tuplesizes(&ro_vars);
    init_index(&ro_vars);

    // 处理所有输入文件
    for (int i = optind; i < argc; ++i) {
        trf_search_file(argv[i], &opt, n_threads, tbatch_size, &ro_vars);
    }

    // 输出运行信息
    fprintf(stderr, "[M::%s] Version: %s\n", __func__, TRFx_VERSION);
    fprintf(stderr, "[M::%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, " %s", argv[i]);
    }
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", 
            __func__, realtime() - mm_realtime0, cputime());

    return 0;
}


