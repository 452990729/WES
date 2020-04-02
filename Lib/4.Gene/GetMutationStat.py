#!/usr/bin/env python2


import sys
import re
import pandas as pd


def ReadTable(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=None)
    return pd_data

def TransformData(pd_data):
    for index in pd_data.index:
        for column in pd_data.columns:
            if column != '#Chr':
                mut = re.split(':', pd_data.loc[index, column])[0]
                if mut  == './.':
                    pd_data.loc[index, column] = '-'
                elif mut == '0/0':
                    pd_data.loc[index, column] = '0'
                else:
                    pd_data.loc[index, column] = '1'
            else:
                break
    return pd_data

def main():
    pd_data = ReadTable(sys.argv[1])
    pd_out = TransformData(pd_data)
    pd_out.to_csv(sys.argv[2], sep='\t', header=True, index=True)


if __name__ == '__main__':
    main()
