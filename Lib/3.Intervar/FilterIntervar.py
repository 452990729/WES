#!/usr/bin/env python2


import os
import sys
import re
import argparse
import pandas as pd


def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None)
    return pd_data

def GetROI(pd_data):
    Func_tp = ['exonic', 'splicing', 'UTR3', 'UTR5']
    pd_out = pd.DataFrame(columns=pd_data.columns)
    for index in pd_data.index:
        if pd_data.loc[index, 'Func.refGene'] in Func_tp:
            pd_out.loc[index, pd_data.columns] = pd_data.loc[index,:]
    return pd_out

def Classify(pd_data):
    l1 = pd.DataFrame(columns=['level']+list(pd_data.columns))
    for index in pd_data.index:
        if 'Pathogenic' in pd_data.loc[index, ' InterVar: InterVar and Evidence ']:
            l1.loc[index, pd_data.columns] = pd_data.loc[index,:]
            l1.loc[index, 'level'] = '1'
        elif 'Likely pathogenic' in  pd_data.loc[index, ' InterVar: InterVar and Evidence ']:
            l1.loc[index, pd_data.columns] = pd_data.loc[index,:]
            l1.loc[index, 'level'] = '2'
        elif 'Uncertain significance' in pd_data.loc[index, ' InterVar: InterVar and Evidence ']:
            l1.loc[index, pd_data.columns] = pd_data.loc[index,:]
            l1.loc[index, 'level'] = '3'
    return l1.sort_values(by='level', ascending=True)


def main():
    parser = argparse.ArgumentParser(description="Filter SNP")
    parser.add_argument('-i', help='input file', required=True)
    parser.add_argument('-o', help='the outputpath', default='./')
    argv=vars(parser.parse_args())
    pd_data = ReadData(argv['i'])
    lb = re.split('\.', os.path.basename(argv['i']))[1]
    pd_roi = GetROI(pd_data)
#    pd_roi.to_excel(os.path.join(argv['o'], 'Final{}ROI.xls'.format(lb)), header=True, index=False)
    pd_roi.to_csv(os.path.join(argv['o'], 'Final{}ROI.txt'.format(lb)), sep='\t', header=True, index=False)
#    lb = re.split('\.', os.path.basename(sys.argv[1]))[2]
    pd_out = Classify(pd_roi)
#    pd_out.to_excel(os.path.join(argv['o'], 'Final{}Result.xls'.format(lb)), header=True, index=False)
    pd_out.to_csv(os.path.join(argv['o'], 'Final{}Result.txt'.format(lb)), sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()
