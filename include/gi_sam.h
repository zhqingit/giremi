/*
 * =====================================================================================
 *
 *       Filename:  sam.h
 *
 *    Description:  Header file for the samtools
 *
 *        Version:  1.0
 *        Created:  2015年01月13日 10时18分04秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  QING ZHANG (), zhqingaca@gmail.com
 *   Organization:  University of California, Los Angeles
 *
 * =====================================================================================
 */

#include "htslib/sam.h"
//#include "bam.h"

/*! @typedef
  @abstract SAM/BAM file handler
  @field  type    type of the handler; bit 1 for BAM, 2 for reading and bit 3-4 for flag format
  @field  bam   BAM file handler; valid if (type&1) == 1
  @field  tamr  SAM file handler for reading; valid if type == 2
  @field  tamw  SAM file handler for writing; valid if type == 0
  @field  header  header struct
 */
typedef struct {
    samFile *file;
    struct { BGZF *bam; } x;  // Hack so that fp->x.bam still works
    bam_hdr_t *header;
} samfile_t;

