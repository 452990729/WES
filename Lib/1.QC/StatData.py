#!/usr/bin/env python2


import sys
import re
import os
import argparse
import ConfigParser
import numpy as np
import gzip
import json
import pysam
from glob import glob
import pandas as pd

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../../Bin/config.ini')
BAMDSK = config.get('SOFTWARE', 'bamdst')

def GetTheLengthOfFile(file_in):
    count = 0
    if file_in.endswith('.gz'):
        f1 = gzip.open(file_in, 'rb')
    else:
        f1 = open(file_in, 'rb')
    while True:
        buffers = f1.read(1024 * 8192)
        if not buffers:
            break
        count += buffers.count('\n')
    f1.close()
    return count

def GetAlignmentOfBam(file_in):
    if not os.path.exists(file_in+'.bai'):
        pysam.index(file_in)
    samfile = pysam.AlignmentFile(file_in, 'rb')
    total = samfile.mapped
    return total

def GetSequenceQC(file_in):
    with open(file_in, 'r') as f:
        content = f.read()
    data = json.loads(content)
    raw = data['summary']['before_filtering']['total_reads']
    qc = data['summary']['after_filtering']['total_reads']
    q20 = data['summary']['after_filtering']['q20_rate']
    q30 = data['summary']['after_filtering']['q30_rate']
    gc = data['summary']['after_filtering']['gc_content']
    return raw,qc,q20,q30,gc

def GetCoverage(file_in, bedfile):
    depth = file_in.rstrip('.sort.dedup.bam')
    if not os.path.exists(depth):
        os.mkdir(depth)
    if not os.path.exists(depth+'/coverage.report'):
        os.system('{} -p {} -o {} {}'.format(BAMDSK, bedfile, depth, file_in))
    dict_tmp = {}
    with open(depth+'/coverage.report', 'r') as f:
        for line in f.readlines()[3:]:
            list_split = re.split('\t', line.strip())
            dict_tmp[list_split[0]] = list_split[1].rstrip('%')
    AlignedReads = int(dict_tmp['[Total] Mapped Reads'])
    UnAlignedReads = int(dict_tmp['[Total] Paired Reads'])-int(AlignedReads)
    AligendRate = dict_tmp['[Total] Fraction of Mapped Reads'].rstrip('%')
    DuplicateReads = int(dict_tmp['[Total] Paired Reads'])-int(dict_tmp['[Total] Read1(rmdup)'])-int(dict_tmp['[Total] Read2(rmdup)'])
    ReadsAfterDedup = int(dict_tmp['[Total] Read1(rmdup)'])+int(dict_tmp['[Total] Read2(rmdup)'])
    DuplicateRate = round(float(DuplicateReads*100)/int(dict_tmp['[Total] Paired Reads']), 2)
    TargetReads = int(dict_tmp['[Target] Target Reads'])
    UnTargetReads = int(dict_tmp['[Total] Paired Reads'])-int(dict_tmp['[Target] Target Reads'])
    TargetRate = dict_tmp['[Target] Fraction of Target Reads in all reads'].rstrip('%')
    AverageSequencingDepth = dict_tmp['[Target] Average depth']
    Coverage = dict_tmp['[Target] Coverage (>0x)'].rstrip('%')
    CoverageAtLeast4x = dict_tmp['[Target] Coverage (>=4x)'].rstrip('%')
    CoverageAtLeast10x = dict_tmp['[Target] Coverage (>=10x)'].rstrip('%')
    CoverageAtLeast30x = dict_tmp['[Target] Coverage (>=30x)'].rstrip('%')
    CoverageAtLeast100x = dict_tmp['[Target] Coverage (>=100x)'].rstrip('%')
    return AlignedReads,UnAlignedReads,AligendRate,DuplicateReads,ReadsAfterDedup,DuplicateRate,TargetReads,UnTargetReads,TargetRate,AverageSequencingDepth, Coverage, CoverageAtLeast4x, CoverageAtLeast10x, CoverageAtLeast30x, CoverageAtLeast100x

def GetDepth(file_in):
    list_dep = []
    with gzip.open(file_in, 'rb') as f:
        for line in f:
            if not line.startswith('#'):
                list_split = re.split('\t', line)
                list_dep.append(float(list_split[3]))
    [n5,n10,n25,n50,n75,n100] = list(np.percentile(list_dep, [5, 10, 25, 50, 75, 100]))
#    [n5,n10,n25,n50,n75,n100] = list(np.percentile(list_dep, [95, 90, 75, 50, 25, 0]))
    return n5,n10,n25,n50,n75,n100

def ExonStat(file_in):
    m = 0
    total = 0
    dep = 0
    cov = 0
    list_cov = []
    list_dep = []
    list_out = []
    with gzip.open(file_in, 'rb') as f:
        for line in f:
            if not line.startswith('#'):
                list_split = re.split('\t', line)
                m += 1
                total += int(list_split[2])-int(list_split[1])
                dep += float(list_split[3])
                cov += float(list_split[6])
                list_cov.append(float(list_split[6]))
                list_dep.append(float(list_split[3]))
    for i in [5,10,20,30,50,100,300]:
        list_filter = filter(lambda x:x>=i, list_dep)
        list_out.append(round(len(list_filter)*100/float(m), 2))
    list_out.append(round(dep/m, 2))
    for i in [5, 20,50, 80, 100]:
        list_filter = filter(lambda x:x>=i, list_cov)
        list_out.append(round(len(list_filter)*100/float(m), 2))
    list_out.append(round(cov/m, 2))
    list_out.append(total)
    list_out.append(round(total/m, 0))
    return list_out

def StatVCF(file_in):
    dict_tmp = {}
    with open(file_in, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                list_sample = re.split('\t', line.strip())[9:]
                for i in list_sample:
                    dict_tmp[i] = 0
            elif not line.startswith('#'):
                list_tmp = re.split('\t', line.strip())
                list_split = list_tmp[9:]
                if list_tmp[6] != 'Filter':
                    for i in range(len(list_split)):
                        if not list_split[i].startswith('0/0') and not list_split[i].startswith('./.'):
                            dict_tmp[list_sample[i]] += 1
                else:
                    print 1
    return dict_tmp

def StatTable(file_in):
    with open(file_in, 'r') as f:
        list_fl = f.readlines()
    list_sample = []
    list_tp = ['exonic', 'intergenic', 'intronic', 'splicing', 'UTR3', 'UTR5', 'downstream', 'upstream']
    for line in re.split('\t', list_fl[0].strip()):
        if line != 'level' and line != '#Chr':
            list_sample.append(line)
        else:
            break
    pd_feature = pd.DataFrame(index=list_sample, columns=['#Raw', 'exonic', 'intergenic', 'intronic', 'splicing', 'UTR3', 'UTR5', 'downstream', 'upstream', 'PositionOther'])
    pd_feature.loc[:,:] = 0
    for line in list_fl[1:]:
        list_split = re.split('\t', line)
        if list_split[-66] in list_tp:
            tp = list_split[-66]
        else:
            tp = 'PositionOther'
        for i in range(len(list_sample)):
            if list_split[i].startswith('0/0'):
                pd_feature.loc[list_sample[i], tp] += 0
                pd_feature.loc[list_sample[i], '#Raw'] += 0
            else:
                pd_feature.loc[list_sample[i], tp] += 1
                pd_feature.loc[list_sample[i], '#Raw'] += 1
    return pd_feature

def StatQC(path_in, bedfile):
    RawData = path_in+'/RawData'
    QC = path_in+'/1.QC'
    Align = path_in+'/2.Align'
    list_sample = []
    for fl in glob(RawData+'/*_1.fq.gz'):
        lb = re.split('_1\.fq\.gz', os.path.basename(fl))[0]
        list_sample.append(lb)
    pd_stat = pd.DataFrame(index=sorted(list_sample), columns=['#RawReads', '#ReadsAfterQC', 'QCRate(%)', 'Q20','Q30','GC','#AlignedReads', '#UnAlignedReads', 'AligendRate(%)', '#DuplicateReads', '#ReadsAfterDedup', 'DuplicateRate(%)', '#TargetReads', '#UnTargetReads', 'TargetRate(%)', 'AverageSequencingDepth(x)', 'Percentile5Depth', 'Percentile10Depth', 'Percentile25Depth', 'Percentile50Depth', 'Percentile75Depth', 'Percentile100Depth', 'Coverage(%)', 'CoverageAtLeast4x(%)', 'CoverageAtLeast10x(%)', 'CoverageAtLeast30x(%)', 'CoverageAtLeast100x(%)', 'ExonDepthAtLeast5x(%)', 'ExonDepthAtLeast10x(%)', 'ExonDepthAtLeast20x(%)', 'ExonDepthAtLeast30x(%)', 'ExonDepthAtLeast50x(%)', 'ExonDepthAtLeast100x(%)', 'ExonDepthAtLeast300x(%)', 'ExonAvergeDepth', 'ExonCoverageAtLeast5x', 'ExonCoverageAtLeast20x', 'ExonCoverageAtLeast50x', 'ExonCoverageAtLeast80x', 'ExonCoverageAtLeast100x', 'ExonAverageCoverage', 'ExonTotalLength', 'ExonAvarageLength'])
    for fl in list_sample:
        raw,qc,q20,q30,gc = GetSequenceQC(os.path.join(QC, fl+'_QC_report.json'))
        QCRate = round((float(qc)*100)/raw, 2)
        AlignedReads,UnAlignedReads,AligendRate,DuplicateReads,ReadsAfterDedup,DuplicateRate,TargetReads,UnTargetReads,TargetRate,AverageSequencingDepth,Coverage,CoverageAtLeast4x,CoverageAtLeast10x,CoverageAtLeast30x,CoverageAtLeast100x = GetCoverage(os.path.join(Align, fl+'.sort.dedup.bam'), bedfile)
        n5,n10,n25,n50,n75,n100 = GetDepth(Align+'/'+fl+'/'+'region.tsv.gz')
        [ExonDepthAtLeast5x,ExonDepthAtLeast10x,ExonDepthAtLeast20x,ExonDepthAtLeast30x,ExonDepthAtLeast50x,ExonDepthAtLeast100x,ExonDepthAtLeast300x,ExonAvergeDepth,ExonCoverageAtLeast5x,ExonCoverageAtLeast20x,ExonCoverageAtLeast50x,ExonCoverageAtLeast80x,ExonCoverageAtLeast100x,ExonAverageCoverage,ExonTotalLength,ExonAvarageLength] = ExonStat(Align+'/'+fl+'/'+'region.tsv.gz')
        pd_stat.loc[fl, :] = [raw,qc,QCRate,q20,q30,gc,AlignedReads,UnAlignedReads,AligendRate,DuplicateReads,ReadsAfterDedup,DuplicateRate,TargetReads,UnTargetReads,TargetRate,AverageSequencingDepth,n5,n10,n25,n50,n75,n100,Coverage,CoverageAtLeast4x,CoverageAtLeast10x,CoverageAtLeast30x,CoverageAtLeast100x,ExonDepthAtLeast5x,ExonDepthAtLeast10x,ExonDepthAtLeast20x,ExonDepthAtLeast30x,ExonDepthAtLeast50x,ExonDepthAtLeast100x,ExonDepthAtLeast300x,ExonAvergeDepth,ExonCoverageAtLeast5x,ExonCoverageAtLeast20x,ExonCoverageAtLeast50x,ExonCoverageAtLeast80x,ExonCoverageAtLeast100x,ExonAverageCoverage,ExonTotalLength,ExonAvarageLength]
    return pd_stat

def StatSNP(path_in):
    rawsnp = path_in+'/6.Final/AllSamplesnp.txt'
    roisnp = path_in+'/6.Final/FinalSamplesnpROI.txt'
    finalsnp = path_in+'/6.Final/FinalSamplesnp.txt'
    pd_feature = StatTable(rawsnp)
    pd_roi = StatTable(roisnp)
    pd_final = StatTable(finalsnp)
    pd_feature['#ROISNP'] = pd_roi['#Raw']
    pd_feature['#FinalSNP'] = pd_final['#Raw']
    return pd_feature

def StatINDEL(path_in):
    rawindel = path_in+'/6.Final/AllSampleindel.txt'
    roiindel = path_in+'/6.Final/FinalSampleindelROI.txt'
    finalindel = path_in+'/6.Final/FinalSampleindel.txt'
    pd_feature = StatTable(rawindel)
    pd_roi = StatTable(roiindel)
    pd_final = StatTable(finalindel)
    pd_feature['#ROIINDEL'] = pd_roi['#Raw']
    pd_feature['#FinalINDEL'] = pd_final['#Raw']
    return pd_feature

def main():
    parser = argparse.ArgumentParser(description="WES STAT pipeline")
    parser.add_argument('-i', help='the wes analysis main path', required=True)
    parser.add_argument('-b', help='the bedfile', default=BasePath+'/../../Database/hg19/Agilent_V5.ex_region.sort.bed')
    parser.add_argument('-o', help='the output file', default='./qcstat.xls')
    argv=vars(parser.parse_args())
    pd_qc = StatQC(argv['i'], argv['b'])
    qc_out = pd.ExcelWriter(argv['o'])
    pd_qc.loc[:, ['#RawReads', '#ReadsAfterQC', 'QCRate(%)', 'Q20', 'Q30', 'GC']].to_excel(qc_out, sheet_name='Sequencing QC', header=True, index=True)
    pd_qc.loc[:, ['#AlignedReads', '#UnAlignedReads', 'AligendRate(%)', '#DuplicateReads', '#ReadsAfterDedup', 'DuplicateRate(%)']].to_excel(qc_out, sheet_name='Alignment QC', header=True, index=True)
    pd_qc.loc[:, ['#TargetReads', '#UnTargetReads', 'TargetRate(%)']].to_excel(qc_out, sheet_name='Capture QC', header=True, index=True)
    pd_qc.loc[:,['ExonDepthAtLeast5x(%)', 'ExonDepthAtLeast10x(%)', 'ExonDepthAtLeast20x(%)', 'ExonDepthAtLeast30x(%)', 'ExonDepthAtLeast50x(%)', 'ExonDepthAtLeast100x(%)', 'ExonDepthAtLeast300x(%)', 'ExonAvergeDepth', 'ExonCoverageAtLeast5x', 'ExonCoverageAtLeast20x', 'ExonCoverageAtLeast50x', 'ExonCoverageAtLeast80x', 'ExonCoverageAtLeast100x', 'ExonAverageCoverage', 'ExonTotalLength', 'ExonAvarageLength']].to_excel(qc_out, sheet_name='Capture QC(exon level)', header=True, index=True)
    pd_qc.loc[:,['Percentile5Depth', 'Percentile10Depth', 'Percentile25Depth', 'Percentile50Depth', 'Percentile75Depth', 'Percentile100Depth']].to_excel(qc_out, sheet_name='Depth distribution (bp)', header=True, index=True)
    pd_qc.loc[:,['Coverage(%)', 'CoverageAtLeast4x(%)', 'CoverageAtLeast10x(%)', 'CoverageAtLeast30x(%)', 'CoverageAtLeast100x(%)']].to_excel(qc_out, sheet_name='Coverage distribution(bp)', header=True, index=True)
    pd_snp =StatSNP(argv['i'])
    pd_snp.to_excel(qc_out, sheet_name='SNP STAT', header=True, index=True)
    pd_indel = StatINDEL(argv['i'])
    pd_indel.to_excel(qc_out, sheet_name='INDEL STAT', header=True, index=True)
    qc_out.close()


if __name__ == '__main__':
    main()
