#!/usr/bin/env python3


import sys
import re
import argparse
import pandas as pd

def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None)
    return pd_data

def Jk(pd_data, index, string, fre):
    if pd_data.loc[index, string]=='.':
        return True
    else:
        if float(pd_data.loc[index, string])<fre:
            return True
        else:
            return False

def FilterMut(pd_data, fre):
    list_index = []
    for index in pd_data.index:
        if Jk(pd_data, index, 'ExAC_ALL', fre) and Jk(pd_data, index, 'esp6500siv2_all', fre) and Jk(pd_data, index, '1000g2015aug_all', fre):
            list_index.append(index)
    return pd_data.loc[list_index,:]

def main():
    parser = argparse.ArgumentParser(description="Filter mut frequency")
    parser.add_argument('-i', help='input intergrate mutation result', required=True)
    parser.add_argument('-f', help='frequency<<0.01>>', default='0.01')
    parser.add_argument('-o', help='output file<<FilterMutation.xlsx>>', default='FilterMutation.xlsx')
    argv=vars(parser.parse_args())
    pd_data = ReadData(argv['i'])
    pd_out = FilterMut(pd_data, float(argv['f']))
    pd_out.to_excel(argv['o'], header=True, index=False, engine='openpyxl')


if __name__ == '__main__':
    main()

