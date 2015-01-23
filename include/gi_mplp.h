/*
 * =====================================================================================
 *
 *       Filename:  gi_mplp.h
 *
 *    Description:  The header file for pileup
 *
 *        Version:  1.0
 *        Created:  2015年01月14日 15时01分30秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  QING ZHANG (), zhqingaca@gmail.com
 *   Organization:  University of California, Los Angeles
 *
 * =====================================================================================
 */

#include "htslib/kstring.h"

#define MPLP_BCF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)
#define B2B_FMT_DP      (1<<0)
#define B2B_FMT_SP      (1<<1)
#define B2B_FMT_DV      (1<<2)
#define B2B_FMT_DP4     (1<<3)
#define B2B_FMT_DPR     (1<<4)
#define B2B_INFO_DPR    (1<<5)

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);


typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
    int rflag_require, rflag_filter;
    int openQ, extQ, tandemQ, min_support; // for indels
    double min_frac; // for indels
    char *reg, *pl_list, *fai_fname, *output_fname;
    faidx_t *fai;
    void *bed, *rghash;
    int argc;
    char **argv;
} mplp_conf_t;

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    bam_hdr_t *h;
    int ref_id;
    char *ref;
    const mplp_conf_t *conf;
} mplp_aux_t;


typedef struct {
    int n, m;
    char **smpl;
    void *rg2smid, *sm2id;
} bam_sample_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;
