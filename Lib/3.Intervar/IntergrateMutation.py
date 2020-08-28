#!/usr/bin/env python2


import sys
import re
import os
import argparse
import pandas as pd


BinPath = os.path.split(os.path.realpath(__file__))[0]
Database = BinPath+'/../../Database/'

def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None)
    return pd_data

def MergeData(pd_snp, pd_indel):
    list_column = list(pd_snp.columns)+['Mutation_type']
    list_label = ['SNP']*pd_snp.shape[0]+['INDEL']*pd_indel.shape[0]
    pd_out = pd_snp.append(pd_indel, ignore_index=True)
    pd_out['Mutation_type'] = list_label
    pd_out = pd_out.sort_values(by=['#Chr', 'Start'])
    return pd_out[list_column]

def AddKnowHpo(file_in, lb, pd_data):
    list_tmp = []
    list_out = []
    with open(file_in, 'r') as f:
        for line in f:
            list_tmp.append(line.strip())
    for gene in pd_data.loc[:,'Ref.Gene']:
        if gene in list_tmp:
            list_out.append('Yes')
        else:
            list_out.append('No')
    pd_data[lb] = list_out
    return pd_data

def main():
    parser = argparse.ArgumentParser(description="Merge indel and snp")
    parser.add_argument('-s', help='snp file', required=True)
    parser.add_argument('-i', help='indel file', required=True)
    parser.add_argument('-o', help='output file<<MergedMatrix.txt>>', default='MergedMatrix.txt')
    argv=vars(parser.parse_args())
    pd_snp = ReadData(argv['s'])
    pd_indel = ReadData(argv['i'])
    pd_out = MergeData(pd_snp, pd_indel)
    pd_out = AddKnowHpo(Database+'/Dermatology/hpo.assoicated.genes.lst', 'Dermatology_hpo', pd_out)
    pd_out = AddKnowHpo(Database+'/Dermatology/inherit.genes.lst', 'Dermatology_hgmd', pd_out)
    pd_out.to_csv(argv['o'], sep='\t', header=True, index=False)
    pd_out.to_excel(argv['o'].rstrip('txt')+'xlsx', header=True, index=False, engine='openpyxl')


if __name__ == '__main__':
    main()
