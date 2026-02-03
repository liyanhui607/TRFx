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
#define max(a, b) (((a) >= (b)) ? (a) : (b))
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


static void usage_mod(FILE *fp, trf_opt opt)
{
	fprintf(stderr, "Usage: trfx [options] <in.fa>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Parameters:\n");
	fprintf(stderr, "    -a INT     Match = matching weight [%d]\n", opt.match);
	fprintf(stderr, "    -b INT     Mismatch = mismatching penalty [%d]\n", -opt.mismatch);
	fprintf(stderr, "    -g INT     Delta = indel penalty [%d]\n", -opt.indel);
	fprintf(stderr, "    -k INT     PM = match probability (whole number; 75 or 80) [%d]\n", opt.PM);
	fprintf(stderr, "    -i INT     PI = indel probability (whole number; 10 or 20) [%d]\n", opt.PI);
	fprintf(stderr, "    -s INT     Minscore = minimum alignment score to report [%d]\n", opt.minscore);
	fprintf(stderr, "    -p INT     MaxPeriod = maximum period size to report, within [1,2000] [%d]\n", opt.maxperiod);
	fprintf(stderr, "    -l INT     maximum TR length expected (in millions) [%d]\n", (int)(opt.maxwraplength*1e-6+.499));
	fprintf(stderr, "    -t INT     number of threads [%d]\n", opt.n_threads);
    fprintf(stderr, "  more output formats:\n");
	fprintf(stderr, "    -n         output in the TRF NGS format\n");
	fprintf(stderr, "    -d         output in the TRF dat format\n");
	fprintf(stderr, "    -m         output masked sequence file \n");
	fprintf(stderr, "    -r         don't eliminate redundancy  \n");
	fprintf(stderr, "    -v         print versioning information\n");
	fprintf(stderr, "Notes:\n");
	fprintf(stderr, "  * BED output format (NB: length of pattern may differ from period):\n");
	fprintf(stderr, "      ctg start end period copyNum fracMatch fracGap score entroy pattern\n");
	fprintf(stderr, "  * TRF NGS output format:\n");
	fprintf(stderr, "      start end period copyNum patLen %%Match %%Gap score %%A %%C %%G %%T entroy pattern seq\n");
	fprintf(stderr, "  * Larger -l: faster but using more memory\n");
	fprintf(stderr, "  * Smaller -b: more sensitive but slower\n");
}

int main(int argc, char *argv[])
{
    // 初始化默认值和选项
    
    liftrlimit();
    mm_realtime0 = realtime();

    // 使用复合字面量初始化选项
    trf_opt opt = {
        .maxwraplength = 12000000,//12000000
        .match = 2,
        .mismatch = -7,
        .indel = -7,
        .PM = 80,
        .PI = 10,
        .minscore = 30, //30
        .maxperiod = 500,
        .masked = 0,
        .redundoff = 0, 
        .datafile = 0,
        .ngs = 0,
        .n_threads = 8,
        .ngsfilename = "",
        .datafilename = "",
        .maskfilename = ""
    };


    
    // 解析命令行选项
    int c;
    while ((c = getopt(argc, argv, "VvAa:Bb:Gg:Ii:Ss:Pp:Ll:Kk:Tt:MmNnDdRr")) >= 0) {
        switch (c) {
            
            case 'A': opt.match = atoi(optarg); break; 
            case 'a': opt.match = atoi(optarg); break;
            case 'B': opt.mismatch = -atoi(optarg); break;
            case 'b': opt.mismatch = -atoi(optarg); break;
            case 'G': opt.indel = -atoi(optarg); break;
            case 'g': opt.indel = -atoi(optarg); break;
            case 'K': opt.PM = atoi(optarg); break;
            case 'k': opt.PM = atoi(optarg); break;
            case 'I': opt.PI = atoi(optarg); break;
            case 'i': opt.PI = atoi(optarg); break;
            case 'S': opt.minscore = atoi(optarg); break;
            case 's': opt.minscore = atoi(optarg); break;
            case 'P': opt.maxperiod = atoi(optarg); break;
            case 'p': opt.maxperiod = atoi(optarg); break;
            case 'L': opt.maxwraplength = max(opt.maxwraplength, atoi(optarg) * 1000000); break;
            case 'l': opt.maxwraplength = max(opt.maxwraplength, atoi(optarg) * 1000000); break;
            case 'T': opt.n_threads = atoi(optarg); break;
            case 't': opt.n_threads = atoi(optarg); break;
            case 'M': opt.masked = 1; break;
            case 'm': opt.masked = 1; break;
            case 'N': opt.ngs = 1; break;
            case 'n': opt.ngs = 1; break;
            case 'R': opt.redundoff = 1; break;
            case 'r': opt.redundoff = 1; break;
            case 'D': opt.datafile = 1; break;
            case 'd': opt.datafile = 1; break;
            case 'V': puts(TRFx_VERSION); return 0;
            case 'v':
                puts(TRFx_VERSION);
                return 0;
        }
    }


    if (optind  >= argc) {
		usage_mod(stderr, opt);
		exit(1);
	}




    // 获取第一个输入文件名并计算最优批处理大小
    const char *input_filename = argv[optind];
    char  tmp_str[50];
    sprintf(tmp_str,"%d.%d.%d.%d.%d.%d.%d",
			opt.match,-opt.mismatch,-opt.indel,
			opt.PM,opt.PI,opt.minscore,opt.maxperiod);
    sprintf(opt.maskfilename,"%s.%s.mask",input_filename,tmp_str);
    sprintf(opt.datafilename,"%s.%s.dat",input_filename,tmp_str);
    sprintf(opt.ngsfilename,"%s.%s.ngs",input_filename,tmp_str);
    size_t tbatch_size = get_optimal_batch_size(input_filename);

    fprintf(stderr, "Current batch_size:  (%.2f GB)\n", tbatch_size / 1024.0 / 1024 / 1024);

    // 初始化只读变量
    readonly_vars_struct ro_vars;
    init_readonly(&ro_vars, &opt);
    init_and_fill_coin_toss_stats2000_with_4tuplesizes(&ro_vars);
    init_index(&ro_vars);

    // 处理所有输入文件
    for (int i = optind; i < argc; ++i) {
        trf_search_file(argv[i], &opt, tbatch_size, &ro_vars);
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


