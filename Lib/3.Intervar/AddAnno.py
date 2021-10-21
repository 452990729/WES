#!/usr/bin/env python3
# -*- coding: utf-8 -*

import os
import sys
import re
import argparse
import pandas as pd

BinPath = os.path.split(os.path.realpath(__file__))[0]
Database = BinPath+'/../../Database/'

def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None, error_bad_lines=False)
    return pd_data

def AddHPO(pd_hpo, pd_anno):
    list_HPOTermID = []
    list_HPOTermName = []
    for index in pd_anno.index:
        gene = pd_anno.loc[index,'Ref.Gene']
        if gene in list(pd_hpo.loc[:,'GeneSymbol']):
            re_index = list(pd_hpo[pd_hpo['GeneSymbol']==gene].index)[0]
            list_HPOTermID.append(pd_hpo.loc[re_index, 'HPOTermID'])
            list_HPOTermName.append(pd_hpo.loc[re_index, 'HPOTermName'])
        else:
            list_HPOTermID.append('.')
            list_HPOTermName.append('.')
    pd_anno['HPOTermID'] = list_HPOTermID
    pd_anno['HPOTermName'] = list_HPOTermName
    return pd_anno

def AddCHPO(pd_chpo, pd_anno):
    list_chpo = []
    list_ref = list(pd_chpo['HPO编号'])
    pd_chpo.index = list_ref
    for i in pd_anno['HPOTermID']:
        list_split = re.split(',', i)
        list_tmp = []
        for m in list_split:
            if m in list_ref:
                list_tmp.append(str(pd_chpo.loc[m, '中文译名']))
        list_chpo.append(','.join(list_tmp))
    pd_anno['HPOTermChineseName'] = list_chpo
    return pd_anno

def AddHgmd(pd_hgmd, pd_anno):
    list_index = []
    for index in pd_hgmd.index:
        lb = pd_hgmd.loc[index, '#Chr']+'-'+str(pd_hgmd.loc[index, 'Start'])+'-'+str(pd_hgmd.loc[index, 'Stop'])+'-'+str(pd_hgmd.loc[index, 'Refer'])+'-'+str(pd_hgmd.loc[index, 'Call'])
        list_index.append(lb)
    pd_hgmd.index = list_index
    list_disease = []
    list_pmid = []
    for index in pd_anno.index:
        lb = 'chr'+str(pd_anno.loc[index, '#Chr'])+'-'+str(pd_anno.loc[index, 'Start'])+'-'+str(pd_anno.loc[index, 'End'])+'-'+str(pd_anno.loc[index, 'Ref'])+'-'+str(pd_anno.loc[index, 'Alt'])
        if lb in list(pd_hgmd.index):
            list_disease.append(pd_hgmd.loc[lb, 'disease'])
            list_pmid.append(pd_hgmd.loc[lb, 'pmid'])
        else:
            list_disease.append('.')
            list_pmid.append('.')
    pd_anno['HGMD_Disease'] = list_disease
    pd_anno['HGMD_Pmid'] = list_pmid
    return pd_anno

def main():
    parser = argparse.ArgumentParser(description="Anno HPO and HGMD")
    parser.add_argument('-i', help='the input annofile', required=True)
    parser.add_argument('-o', help='the output file', required=True)
    argv=vars(parser.parse_args())
    pd_anno = ReadData(argv['i'])
    pd_hpo = ReadData(os.path.join(Database, 'hg19/HPO.txt'))
    print(os.path.join(Database, 'hg19/chpo.2018.txt'))
    pd_chpo = ReadData(os.path.join(Database, 'hg19/chpo.2018.txt'))
    pd_hgmd = ReadData(os.path.join(Database, 'hg19/hgmd.2018.allmut.final.sort.0504.txt'))
    pd_anno = AddHPO(pd_hpo, pd_anno)
    pd_anno = AddCHPO(pd_chpo, pd_anno)
    pd_anno = AddHgmd(pd_hgmd, pd_anno)
    pd_anno.to_csv(argv['o'], sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()
