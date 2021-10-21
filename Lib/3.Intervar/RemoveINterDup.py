#!/usr/bin/env python2


import sys
import re
import argparse
import pandas as pd

def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None)
    pd_out = pd_data.drop_duplicates()
    return pd_out

def MergeDup(pd_data, outfile):
    dict_tmp = {}
    dict_omim = {}
    for index in pd_data.index:
        label = str(pd_data.loc[index,'#Chr'])+'_'+str(pd_data.loc[index,'Start'])+\
                '_'+str(pd_data.loc[index,'End'])+'_'+pd_data.loc[index,'Ref']+'_'+\
                pd_data.loc[index,'Alt']
        if label not in dict_tmp and pd_data.loc[index,'Ref.Gene'] != 'NONE':
            dict_tmp[label] = [i if i != 'nan' else '' for i in [str(i) for i in pd_data.loc[index,:]]]
        if label not in dict_omim:
            dict_omim[label] = [pd_data.loc[index,'OMIM'],]
        else:
            dict_omim[label] += [pd_data.loc[index,'OMIM'],]
    pd_out = pd.DataFrame(columns=pd_data.columns)
    out = open(outfile, 'w')
    out.write('\t'.join(pd_data.columns)+'\n')
    for label in dict_tmp:
        list_out = dict_tmp[label]
        list_out[-4] = ','.join([i for i in set(dict_omim[label]) if i != '.'])
        out.write('\t'.join(list_out)+'\n')
#        pd_out.loc[i,'OMIM'] = ','.join(set(dict_omim[label]))
    out.close()

def main():
    parser = argparse.ArgumentParser(description="Filter Dup")
    parser.add_argument('-i', help='input file', required=True)
    parser.add_argument('-o', help='the output', required=True)
    argv=vars(parser.parse_args())
    pd_data = ReadData(argv['i'])
    MergeDup(pd_data, argv['o'])
#    pd_out.to_csv(argv['o'], sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()
