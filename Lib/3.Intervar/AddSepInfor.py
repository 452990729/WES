#!/usr/bin/env python3


import sys
import re
import pandas as pd

def ReadData(file_in):
    pd_data = pd.read_excel(file_in, header=0, index_col=None)
    return pd_data

def GetHet(string_in):
    list_out = []
    list_split = re.split(':', string_in)
    a = list_split[0][0]
    b = list_split[0][2]
    list_s = re.split(',', list_split[1])
    r1 = list_s[0]
    r2 = list_s[1]
    if a == b:
        if a == '.':
            list_out.append('未覆盖')
        elif a == '0':
            list_out.append('未突变')
        else:
            list_out.append('纯合')
    else:
        if a == '0' or b == '0':
            list_out.append('杂合')
        else:
            list_out.append('多重杂合')
    fre = round(float(r2)/(float(r1)+float(r2)), 2)
    list_out += [r1, r2, str(fre)]
    return list_out

def ProcessPd(pd_data):
    list_columns = list(pd_data.columns)
    list_s = []
    for i in list_columns:
        if i != '#Chr':
            list_s.append(i)
        else:
            break
    list_i = []
    for m in list_s:
        list_i += [m+'_MutationType', m+'_Allele1', m+'_Allele2', m+'_AlleleFrequency']
    pd_add = pd.DataFrame(index = pd_data.index, columns=list_i)
    for index in pd_data.index:
        list_tmp = []
        for s in list_s:
            list_tmp += GetHet(pd_data.loc[index, s])
        pd_add.loc[index, :] = list_tmp
    for i in range(len(list_s), 5*len(list_s)):
        pd_data.insert(i, list_i[i-len(list_s)], pd_add.loc[:,list_i[i-len(list_s)]])
    return pd_data

def main():
    pd_data = ReadData(sys.argv[1])
    pd_out = ProcessPd(pd_data)
    lb = sys.argv[1].rstrip('.xlsx')
    pd_out.to_excel(lb+'.final.xlsx', header=True, index=None)


if __name__ == '__main__':
    main()
