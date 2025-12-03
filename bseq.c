#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bseq.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


extern unsigned char seq_nt4_table[256];

struct bseq_file_s
{
	int is_eof;
	gzFile fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const  char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0)
		return 0;
	fp = (bseq_file_t *)calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}


unsigned char *remove_spaces(const unsigned char *input, int *new_length)
{
    if (input == NULL)
    {
        return NULL;
    }

    // 计算去除空格后的字符串长度
    int length = 0;
    const unsigned char *p = input;
    while (*p != '\0')
    {
        length += (*p != ' ');
        p++;
    }

    // 分配对齐内存
    unsigned char *result = NULL;
    int err = posix_memalign((void **)&result, 64, (length + 1) * sizeof(unsigned char));
    if (err != 0 || result == NULL)
    {
        return NULL;
    }

    // 复制并处理字符
    unsigned char *dst = result;
    p = input;
    while (*p != '\0')
    {
        if (*p != ' ')
        {
            if (*p >= 'a' && *p <= 'z')
            {
                *dst++ = *p - 32;  // 小写转大写
            }
            else if (*p >= 'A' && *p <= 'Z')
            {
                *dst++ = *p;  // 直接复制大写
            }
            // 非字母字符被隐式跳过（不复制）
        }
        p++;
    }
    *dst = '\0';

    if (new_length)
    {
        *new_length = length;  // 注意：这里返回的是非空格字符数（包含非字母字符）
    }

    return result;
}


bseq1_t *bseq_read(bseq_file_t *fp, long chunk_size, int *n_)
{
	long size = 0, m, n;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0;
	seqs = 0;
	while (kseq_read(ks) >= 0)
	{
		bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (n >= m)
		{
			m = m ? m << 1 : 256;
			seqs = (bseq1_t *)realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = (unsigned char *) strdup((char *)ks->name.s);
  
		s->seq = remove_spaces(ks->seq.s, &s->l_seq); // get seq and length
		size += seqs[n++].l_seq;
		if (size >= chunk_size)
			break;
	}
	if (n == 0)
		fp->is_eof = 1;
	*n_ = n;
	return seqs;
}

int bseq_eof(bseq_file_t *fp)
{
	return fp->is_eof;
}
