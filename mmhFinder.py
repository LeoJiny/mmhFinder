'''
This program is used for 1st selection by SNP polymorphism
Input a csv file (all SNP)
filter those not satisfying condition
Output a csv file with fewer data
'''
import pandas as pd
import numpy as np
import os
import sys
from pandas import DataFrame
import vcf


def split(s):
    return s.split(' ')
'''
Turn VCF format file into CSV format file
'''
def vcf_to_csv(number):
    file_name = '/home/sdd/JYY/CHROM_' +number+'.vcf'
    vcf_reader = vcf.Reader(filename = file_name)
    all_data = list()
    for record in vcf_reader:
        former_data_index = record.CHROM + ' ' + record.ID + ' ' + str(record.POS) + ' ' \
                            + record.REF + ' ' + str(record.ALT[0]) + ' ' + str(record.QUAL)
        data_index = split(former_data_index)  # 将这个字符串变成一个列表 's'变成['s']
        for a in range(0, len(record.samples)):
            # 按下标访问Call，按.sample访问sample，按键访问FORMAT对应信息
            sample_data = record.samples[a]['GT']
            data_index.append(sample_data)
        # 变成列表以后才可以用来Dataframe
        all_data.append(data_index)
    return all_data



def creat_dt(name_list):
    dt = DataFrame(columns = name_list)
    return dt  #create a new dataframe with the same column name

'''
1st round filter only retains SNP with MAF 
'''
def select_SNP_poly(snplist,first_dt,deleted_dt,MAF):
    i = 0   #row instructor
    frequency = [1-MAF,MAF]
    while i < len(snplist):
        t = snplist.iloc[i][6:].value_counts('1/1')[0]
        if t >= min(frequency) and t <=max(frequency) :
            first_dt = first_dt.append(snplist.iloc[i])
            i += 1
        else:
            deleted_dt = deleted_dt.append(snplist.iloc[i])
            i += 1
    return [first_dt,deleted_dt]

def split_mSNPS(snplist,startrow,namelist):
    mSNP_blank_list = creat_dt(namelist)
    startpos = snplist.iloc[startrow][2]
    endpos = startpos +200
    while startrow < len(snplist) and snplist.iloc[startrow][2] <= endpos:
        mSNP_blank_list = mSNP_blank_list.append(snplist.iloc[startrow])
        startrow += 1
    return mSNP_blank_list

def return_first_index(condition,snplist):
    i = 0
    while i <= len(snplist):
        if snplist.iloc[i][0] == condition:
            return  i
        else:
            i = i+1
def return_last_index(condition,snplist):
    i = return_first_index(condition,snplist)
    while i <= len(snplist)-1:
        if snplist.iloc[i-1][0] == condition and snplist.iloc[i][0] != condition:
            return i-1
        elif i == len(snplist)-1:
            return i
        else:
            i += 1

def count_mSNPs(nchr,snplist,column_name,calculate_data):
    first_row =  return_first_index(nchr,snplist)   #return the beginning counter
    # eg firstrow = 100 then the whole search starts with snplist.iloc[100][2]
    #print(endpos)
    nlist = list(range(1,11))

    last_row = return_last_index(nchr,snplist)
    while first_row <= last_row :
        msnps = split_mSNPS(snplist,first_row)
        x = len(msnps)
        if x in nlist:
            rowname = 'n_' + str(x)
            calculate_data.loc[rowname][column_name] += 1
            first_row += 1
        else:
            calculate_data.loc['n_d10'][column_name] += 1
            first_row +=1
    return calculate_data


'''
用于记录微单倍型信息
'''
def return_mSNP_detail(msnpdt,detaildt):
    v = list()
    v.append(msnpdt.iloc[0,0]) #Chrome
    v.append(msnpdt.iloc[0,1]) #StartID
    v.append(msnpdt.iloc[len(msnpdt)-1,1])
    v.append(msnpdt.iloc[0,2])
    v.append(msnpdt.iloc[len(msnpdt)-1,2])
    v.append(len(msnpdt))
    n = msnpdt.iloc[len(msnpdt)-1,2] - msnpdt.iloc[0,2]
    v.append(n)
    x = return_NA(msnpdt)
    v.append(x)
    y = calculate_Ae(msnpdt)  #add
    v.append(y) #add
    detaildt.loc[detaildt.shape[0]] = v
    return detaildt



'''
Calculate Ae
'''
def calculate_Ae(msnpdt):
    j = 6  # 列计数器
    i = 0  # 行计数器
    singlemSNPlist = list()  # 用于存放单个品系的这一区段的微单倍型
    all_mSNP_list = list()  # 用于存放这一区段所有品系的微单倍型
    while j < int(msnpdt.shape[1]):
        while i < len(msnpdt):
            singlemSNPlist.append(msnpdt.iloc[i][j])
            i += 1
        j += 1
        i = 0
        a = ",".join(singlemSNPlist)
        all_mSNP_list.append(a)
        singlemSNPlist = list()
    y = 0 #(Pi)**2
    d = set(all_mSNP_list)
    for i in d:
        x = all_mSNP_list.count(i)   #Pi
        y += (x / 35) ** 2
    Ae = 1 / y
    return round(Ae,2)


'''
return number of alleles of each microhaplotype
'''

def return_NA(msnpdt):
    j = 6  # 列计数器
    i = 0  # 行计数器
    singlemSNPlist = list()  # 用于存放一个品系的这一区段的微单倍型
    all_mSNP_list = list()  # 用于存放这一区段所有品系的微单倍型
    while j < int(msnpdt.shape[1]):
        while i < len(msnpdt):
            singlemSNPlist.append(msnpdt.iloc[i][j])
            i += 1
        j += 1
        i = 0
        all_mSNP_list.append(singlemSNPlist)
        singlemSNPlist = list()
    stat_list = list()
    for element in all_mSNP_list:
        if element in stat_list:
            pass
        else:
            stat_list.append(element)
    return len(stat_list)

'''
Find Microhaplotypes
'''
def split_mSNPS(snplist,startrow,namelist):
    mSNP_blank_list = creat_dt(namelist)
    startpos = snplist.iloc[startrow][2]
    endpos = startpos +200
    while startrow < len(snplist) and snplist.iloc[startrow][2] <= endpos:
        mSNP_blank_list = mSNP_blank_list.append(snplist.iloc[startrow])
        startrow += 1
    return mSNP_blank_list
'''
Filter microhaplotypes with NA 
'''
def select_MEA(snplist,second_dt,namelist):
    startrow = 0   #根据染色体号给出第一行
    lastrow = len(snplist)   # 该染色体号的最后一行
    while startrow <lastrow:
        msnps = split_mSNPS(snplist, startrow,namelist)  # 从这个startrow开始找微单倍型
        #print(msnps)
        #temporary = pd.DataFrame(columns=['CHROM', 'StartID', 'EndID', 'Startpos',
                                           # 'EndPos', 'Distance','MEA','SNP amount'])
        AE = return_NA(msnps)
        if AE >= 3  and len(msnps) != 1:
            second_dt = second_dt.append(msnps)
            startrow = startrow +1
        else:
            startrow = startrow +1
    second_dt.drop_duplicates(inplace = True)
    return second_dt

