#!/usr/bin/env python2


import re
import sys
import os

SAMTOOLS = '/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Software/samtools'
GATK = '/mnt/dfc_data1/home/lixuefei/Pipeline/Exon/Software/gatk'

list_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',\
           'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',\
           'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',\
           'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
inbam = sys.argv[1]
lable = os.path.basename(sys.argv[1]).rstrip('.bam')

for chr_num in list_chr:
    s_sub = SAMTOOLS+' view -b -h {} {} > {}_{}.bam'.format(inbam, chr_num, lable, chr_num)
    s_sort = GATK+' SortSam --INPUT {}_{}.bam --OUTPUT {}_{}.sort.bam --SORT_ORDER coordinate'.format(lable, chr_num, lable, chr_num)
    s_dup = GATK+' MarkDuplicates --INPUT {}_{}.sort.bam --OUTPUT {}_{}.sort.dedup.bam --METRICS_FILE {}_{}.metrics --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000'.format(lable, chr_num, lable, chr_num, lable, chr_num)
    print s_sub+' && '+s_sort+' && '+s_dup
