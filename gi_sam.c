/*
 * =====================================================================================
 *
 *       Filename:  gi_sam.c
 *
 *    Description:  The functions handling sam/bam files
 *
 *        Version:  1.0
 *        Created:  2015年01月07日 16时57分10秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  QING ZHANG (), zhqingaca@gmail.com
 *   Organization:  University of California, Los Angeles
 *
 * =====================================================================================
 */

/*!
  @abstract Verbose level between 0 and 3; 0 is supposed to disable all
  debugging information, though this may not have been implemented.
 */
#define bam_verbose hts_verbose

#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
//#include "htslib/faidx.h"
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "gi_sam.h"

samfile_t *samopen(const char *fn, const char *mode, const void *aux)
{
    // hts_open() is really sam_open(), except for #define games
    samFile *hts_fp = hts_open(fn, mode);
    if (hts_fp == NULL)  return NULL;

    samfile_t *fp = malloc(sizeof (samfile_t));
    fp->file = hts_fp;
    fp->x.bam = hts_fp->fp.bgzf;
    if (strchr(mode, 'r')) {
        if (aux) hts_set_fai_filename(fp->file, aux);
        fp->header = sam_hdr_read(fp->file);  // samclose() will free this
        if (fp->header->n_targets == 0 && bam_verbose >= 1)
            //fprintf(stderr, "[samopen] no @SQ lines in the header.\n");
            printf("%s","[samopen] no @SQ lines in the header.\n");
    }
    else {
        fp->header = (bam_hdr_t *)aux;  // For writing, we won't free it
        if (fp->file->is_bin || fp->file->is_cram || strchr(mode, 'h')) sam_hdr_write(fp->file, fp->header);
    }

    return fp;
}


int main(int argc, char *argv[]){

	printf("%s\t%s\t%s\n",argv[1],argv[2],argv[3]);
	//samfile_t fp = samopen(argv[1],argv[2],argv[3]);
	samfile_t *fp;
	fp = samopen(argv[1],argv[2],argv[3]);
	//printf("%d\n",argc);
}
