#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   2.2.CNVRRedundancy.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import numpy as np
import pandas as pd
import time
import argparse

def getPerLine(line):
    line = line.split('\t')
    line[7] = 0 if line[7] == "None" else line[7]
    chrom = line[0]
    start = int(line[1])
    end = int(line[2])
    effective_windows = int(line[3])
    length = end - start + 1
    weight_array = [int(line[3])] + [float(x) for x in line[4:]]
    gap = float(line[4])
    repeat = float(line[5])
    gc = float(line[6])
    kmer = float(line[7])
    sd = float(line[-1])
    average = float(line[-2])
    tmp = [float(x) for x in line[8:-2]]
    return tmp, chrom, start, end, effective_windows, length, weight_array, gap, repeat, gc, kmer, sd, average

def weight_average(weight_array_1, weight_array_2):
    weight = []
    effective_windows_1 = weight_array_1[0]
    effective_windows_2 = weight_array_2[0]
    sum_effective_windows = effective_windows_1 + effective_windows_2
    weight.append(sum_effective_windows)
    for i in range(1, len(weight_array_1)):
        tmp_value = (weight_array_1[i]*effective_windows_1 + weight_array_2[i]*effective_windows_2)/sum_effective_windows if sum_effective_windows else 0
        weight.append(tmp_value)
    return weight

def readPrimaryCNVR(primaryCNVRFile):
    hash_ = {}
    with open(primaryCNVRFile, 'r') as fr:
        header = fr.readline().strip()
        firstline = fr.readline().strip()
        tmp_start_array, tmp_chrom, tmp_start, tmp_end, tmp_effective_windows, tmp_length, tmp_weight_array1, tmp_gap, tmp_repeat, tmp_gc, tmp_kmer, tmp_sd, tmp_average = getPerLine(firstline)
        if tmp_chrom not in hash_:
            hash_[tmp_chrom] = {}
        hash_[tmp_chrom][tmp_start] = [tmp_end] + tmp_weight_array1
        for line in fr:
            line = line.strip()
            record = line.split("\t")
            if len(record) == 13:
                tmp, chrom, start, end, effective_windows, length, tmp_weight_array2, gap, repeat, gc, kmer, sd, average = getPerLine(line)
                tmp_cor = pd.Series(tmp).corr(pd.Series(tmp_start_array))
                if tmp_chrom == chrom and ((start-tmp_end+1) <= (tmp_length+length)*percent or (start-tmp_end+1) <= 300) and tmp_cor > correlation:
                    tmp_weight = weight_average(tmp_weight_array1, tmp_weight_array2)
                    hash_[tmp_chrom][tmp_start] = [end] + tmp_weight
                    tmp_start_array = tmp_weight[5:-2]
                    tmp_chrom = chrom
                    tmp_end = end
                    tmp_length = end - tmp_start + 1
                else:
                    tmp_weight_array1 = []
                    if chrom not in hash_:
                        hash_[chrom] = {}
                    hash_[chrom][start] = [end] + tmp_weight_array2
                    tmp_weight_array1 = tmp_weight_array2
                    tmp_start_array = tmp
                    tmp_chrom = chrom
                    tmp_end = end
                    tmp_start = start
                    tmp_length = start - end + 1
    return header, hash_

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='produce mergedCNVR file.')
    parser.add_argument('-i', type=str, help='Primary CNVR file.', required=True)
    parser.add_argument('-o', type=str, help='merged CNVR file generated.', required=True)
    parser.add_argument('-p', type=float, help='the percent for the distance of two adjacent CNV in their CNV length when merge (default: 0.1)', default=0.1)
    parser.add_argument('-cor', type=int, help='the pearson correlation of copy numbers between two adjacent CNV when merge (default: 0.05)', default=0.05)
    args = parser.parse_args()
    percent = args.p
    correlation = args.cor
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': Starting merging...')
    header, hash_ = readPrimaryCNVR(args.i)
    with open(args.o, 'w') as fw:
        fw.write(f'{header}\n')
        for chrom in sorted(hash_.keys()):
            for pos_start in hash_[chrom].keys():
                value = '\t'.join(map(str, hash_[chrom][pos_start]))
                fw.write(f'{chrom}\t{pos_start}\t{value}\n')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': ===== Done! ======')