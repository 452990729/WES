#!/usr/bin/env python2


import sys
import re
import argparse
import pandas as pd


def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None)
    pd_data.index = ['-'.join([str(m).strip('chr') for m in pd_data.loc[i,:][:2]]) for i in pd_data.index]
    return pd_data

def MergeData(pd_annovar, pd_inter):
    cols_to_use  = filter(lambda x:x not in list(pd_inter.columns), list(pd_annovar.columns))
    pd_out = pd_inter.merge(pd_annovar[cols_to_use], left_index=True, right_index=True, how='inner')
    list_rm = ['avsnp147', 'Freq_ExAC_ALL', 'Freq_esp6500siv2_all', 'Freq_1000g2015aug_all', 'Otherinfo', 'Chr', 'Gene.refGene', 'GeneDetail.refGene']
    cols_remain = filter(lambda x:x not in list_rm, list(pd_out.columns))
    pd_out = pd_out[cols_remain]
    return pd_out

def main():
    parser = argparse.ArgumentParser(description="Merge annovar and intervar")
    parser.add_argument('-i', help='intervar file', required=True)
    parser.add_argument('-a', help='annovar file', required=True)
    parser.add_argument('-o', help='output file<<MergedMatrix.txt>>', default='MergedMatrix.txt')
    argv=vars(parser.parse_args())
    pd_annovar = ReadData(argv['a'])
    pd_inter = ReadData(argv['i'])
    pd_out = MergeData(pd_annovar, pd_inter)
    pd_out.to_csv(argv['o'], sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()



