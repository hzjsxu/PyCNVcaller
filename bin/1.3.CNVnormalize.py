#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   test.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import numpy as np
import argparse
import time

def getRawinfo(correctedCount, repeat_max):
    header = correctedCount.split('/')[-1]
    X_array = []
    global_array = []
    GC_region_average = {}
    genomeCopy= {}
    with open(correctedCount, 'r') as fr:
        for line in fr:
            line = line.strip().split('\t')
            if line[0] not in genomeCopy:
                genomeCopy[line[0]] = {}
            genomeCopy[line[0]][line[1]] = line[-1]
            if line[0] == sex:
                if int(line[4]) <= repeat_max and int(line[-1]) == 1:
                    X_array.append(int(line[2]))
            else:
                if int(line[4]) <= repeat_max and int(line[-1]) == 1:
                    global_array.append(int(line[2]))
                    if int(line[3]) not in GC_region_average:
                        GC_region_average[int(line[3])] = []
                    GC_region_average[int(line[3])].append(int(line[2]))
    return header, X_array, global_array, GC_region_average, genomeCopy

def getGlobalPercentile(global_array):
    global_array = sorted(global_array)
    len_global_array = len(global_array)
    median90 = global_array[int(len_global_array*0.9)]
    median50 = global_array[int(len_global_array*0.5)]
    median10 = global_array[int(len_global_array*0.1)]
    median8 = global_array[int(len_global_array*0.08)]
    median6 = global_array[int(len_global_array*0.06)]
    median4 = global_array[int(len_global_array*0.04)]
    median2 = global_array[int(len_global_array*0.02)]
    print(f'tail 2 percent window is lower than {median2}')
    print(f'tail 4 percent window is lower than {median4}')
    print(f'tail 6 percent window is lower than {median6}')
    print(f'tail 8 percent window is lower than {median8}')
    print(f'tail 10 percent window is lower than {median10}')
    print(f'tail 50 percent window is lower than {median50}')
    print(f'tail 90 percent window is lower than {median90}')
    return median50

def getXPercentile(X_array, sex='X'):
    Xmedian_50 = 0
    if len(X_array) >= 1:
        X_array = sorted(X_array)
        len_X_array = len(X_array)
        Xmedian_90 = X_array[int(len_X_array*0.9)]
        Xmedian_50 = X_array[int(len_X_array*0.5)]
        Xmedian_10 = X_array[int(len_X_array*0.1)]
        print(f'tail 10 percent {sex} window is lower than {Xmedian_10}')
        print(f'tail 50 percent {sex} window is lower than {Xmedian_50}')
        print(f'tail 90 percent {sex} window is lower than {Xmedian_90}')
    return Xmedian_50

def getSexCorrectFold(Xmedian_50, median50, sex='X'):
    sex_correct_fold = 0
    if Xmedian_50>median50*0.75 and Xmedian_50<median50*1.5:
        print(f'sex chromosome {sex} show similar coverage of autosome!')
        sex_correct_fold = 1
    elif Xmedian_50<median50*0.75 and Xmedian_50>median50*0.25:
        print(f'sex chromosome {sex} show half coverage of autosome!')
        sex_correct_fold = 2
    else:
        print(f'sex chromosome {sex} show unknown relationship with autosome?')
        sex_correct_fold = 1
    return sex_correct_fold

def getAveragePerGC(GC_region_average):
    sorted_tmp = []
    region_average_sd = {}
    for tmp_gc in sorted(GC_region_average.keys()):
        sum = 0
        sum_window = 0
        sorted_tmp = sorted(GC_region_average[tmp_gc])
        for i in range(int(l*len(sorted_tmp)), int((1-p)*len(sorted_tmp))+1):
            sum_window += 1
            sum += sorted_tmp[i]
        region_average_sd[tmp_gc] = round(sum/sum_window, 2) if sum_window>0 else 0
    print("calculate average value for each GC content region done!")
    return region_average_sd

def getCorrectPerGC(region_average_sd, correctedCountFile):
    standard_average = region_average_sd[windowSize*0.4]
    for tmp_gc in sorted(region_average_sd.keys()):
        region_average_sd[tmp_gc] = round(standard_average/region_average_sd[tmp_gc],2) if region_average_sd[tmp_gc]>0 else 0
    new_global_array = []
    clean_record = {}
    with open(correctedCountFile, 'r') as fr:
        for line in fr:
            line = line.split('\t')
            if int(line[3]) not in region_average_sd:
                region_average_sd[int(line[3])] = 0
            if line[0] == sex:
                line[2] = round(region_average_sd[int(line[3])]*int(line[2])*sex_correct_fold, 2)
            else:
                line[2] = round(region_average_sd[int(line[3])]*int(line[2]), 2)
                new_global_array.append(line[2])
            if line[0] not in clean_record:
                clean_record[line[0]] = {} 
            clean_record[line[0]][line[1]] = line[2]
    new_global_array = sorted(new_global_array)
    len_new_global_array = len(new_global_array)
    correct_median50 = new_global_array[int(len_new_global_array*0.5)]
    print("correct read count according to GC content done!")
    print(f'global median50 {correct_median50}')
    return clean_record, new_global_array, correct_median50

def calcGlobalMeanSd(new_global_array, clean_record, p=0.05, l=0.05):
    len_new_global_array = len(new_global_array)
    global_max = new_global_array[int((len_new_global_array-1)*(1-p))]
    global_min = new_global_array[int((len_new_global_array-1)*l)]
    global_array = []
    for tmp_chr in clean_record.keys():
        for tmp_pos in clean_record[tmp_chr]:
            if global_min < int(clean_record[tmp_chr][tmp_pos]) < global_max:
                global_array.append(int(clean_record[tmp_chr][tmp_pos]))
                
    global_average = round(np.mean(global_array), 2)
    global_sd = round(np.std(global_array), 2)
    print(f'global average: {global_average}\tglobal SD: {global_sd}')
    print(f'the {p}th percentile absolute reads count: {global_max}')
    print(f'the {l}th percentile absolute reads count: {global_min}')
    return global_average, global_sd

def writeout(outfile, clean_record, correct_median50):
    with open(outfile, 'w') as fw:
        for tmp_chr in sorted(clean_record.keys()):
            for tmp_pos in clean_record[tmp_chr].keys():
                copynumber = round(clean_record[tmp_chr][tmp_pos]/correct_median50, 2)
                fw.write(f'{tmp_chr}\t{tmp_pos}\t{copynumber}\t{genomeCopy[tmp_chr][tmp_pos]}\n') 

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description='GC correction and RD normalization')
    parse.add_argument('-w', type=int, help='window size consistent with reference database (default:1000)', default=1000)
    parse.add_argument('-l', type=float, help='the percentage of the windows with the highest absolute reads count will be excluded during calculating the average of reads count (default: 0.05)', default=0.05)
    parse.add_argument('-p', type=float, help='the percentage of the windows with the lowest absolute reads count will be excluded during calculating the average of reads count (default: 0.05)', default=0.05)
    parse.add_argument('-s', type=str, help='sex chromosome name, i.e. X,chrX,Z,chrZ', default='X')
    parse.add_argument('-i', type=str, help='absolute reads count file', required=True)
    args = parse.parse_args()
    l = args.l
    p = args.p
    sex = args.s
    windowSize = args.w
    CorrectRawCount = args.i
    repeat_max = 0.3*windowSize
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': ****** GC correction and RD normalization Start... ******')
    header, X_array, global_array, GC_region_average, genomeCopy = getRawinfo(correctedCount=CorrectRawCount, repeat_max=repeat_max)
    median50 = getGlobalPercentile(global_array=global_array)
    Xmedian_50 = getXPercentile(X_array=X_array)
    sex_correct_fold = getSexCorrectFold(Xmedian_50=Xmedian_50, median50=median50)
    region_average_sd = getAveragePerGC(GC_region_average=GC_region_average)
    GC_region_average = {}
    clean_record, new_global_array, correct_median50 = getCorrectPerGC(region_average_sd=region_average_sd, correctedCountFile=CorrectRawCount)
    global_average, global_sd = calcGlobalMeanSd(new_global_array=new_global_array, clean_record=clean_record)
    outfileName = header + '_mean_' + str(correct_median50) + '_SD_' + str(global_sd) + '_sex_' + str(sex_correct_fold)
    writeout(outfile=outfileName, clean_record=clean_record, correct_median50=correct_median50)
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': ****** Output normalization copy number done! ******')