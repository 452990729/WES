#!/usr/bin/env python


import sys
import re
import argparse
import pandas as pd

def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None)
    return pd_data

def ReadVCF(file_in):
    dict_tmp = {}
    with open(file_in, 'r') as f:
        for line in f:
            if not line.startswith('##'):
                line = line.strip()
                if line.startswith('#CHROM'):
                    label = re.split('\t', line)[9:]
                else:
                    list_split = re.split('\t', line)
                    dict_tmp[list_split[0].lstrip('chr')+'_'+list_split[1]] = list_split[9:]
    return label, dict_tmp

def Merge(pd_data, label, dict_vcf, outfile):
    out = open(outfile, 'w')
    out.write('\t'.join(label+[str(i) for i in pd_data.columns])+'\n')
    for index in pd_data.index:
        lb = str(pd_data.loc[index, '#Chr'])+'_'+str(pd_data.loc[index, 'Start'])
        if lb in dict_vcf:
            out.write('\t'.join(dict_vcf[lb]+[str(i) for i in pd_data.loc[index,:]])+'\n')
        else:
            list_tmp = re.split('_', lb)
            for i in range(30):
                tlb = list_tmp[0]+'_'+str(int(list_tmp[1])-i)
                if tlb in dict_vcf:
                    out.write('\t'.join(dict_vcf[tlb]+[str(i) for i in pd_data.loc[index,:]])+'\n')
                    break
    out.close()

def main():
    parser = argparse.ArgumentParser(description="Add sample to result")
    parser.add_argument('-i', help='intervar result', required=True)
    parser.add_argument('-v', help='vcf file', required=True)
    parser.add_argument('-o', help='output file<<intervar_vcf.xls>>', default='intervar_vcf.xls')
    argv=vars(parser.parse_args())
    pd_data = ReadData(argv['i'])
    label, dict_vcf = ReadVCF(argv['v'])
    Merge(pd_data, label, dict_vcf, argv['o'])


if __name__ == '__main__':
    main()
