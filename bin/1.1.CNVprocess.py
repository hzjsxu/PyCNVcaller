#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   1.1.CNVprocess.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import pysam
import os
import re
import time
import argparse

## eg. python 1.1.CNVprocess.py -b test.bam -r testDB_1000_500

# 根据headers中的SM标签提取样本名
def getSampleName(bamfile): 
    headers = pysam.view('-H', bamfile)
    sampleName = re.search(r'SM:(.*)\n', headers, re.M|re.I).group(1)
    return sampleName

# 对每个窗口的read进行计数
def getRawCount(bamfile, windowSize=1000, stepSize=500):
    windowSize=windowSize
    stepSize=stepSize
    content = pysam.view('-F', '0x504', bamfile, catch_stdout=True)
    content = content.strip().split('\n')
    uniq = {}
    multi = {}
    for fq_read in content:
        read = fq_read.split('\t')
        chrom = read[2]
        if chrom not in multi:
            multi[chrom] = {}
        if chrom not in uniq:
            uniq[chrom] = {}
        start = int(read[3])
        readLen = len(read[9])
        start += readLen/2
        start = int(start/windowSize)*windowSize+1
        end = int(start/windowSize+0.5)*windowSize-stepSize+1
        if "XA:Z" in fq_read:
            multi[chrom][start] = (multi[chrom][start]+1) if start in multi[chrom] else 1
            multi[chrom][end] = (multi[chrom][end]+1) if end in multi[chrom] else 1
        else:
            uniq[chrom][start] = (uniq[chrom][start]+1) if start in uniq[chrom] else 1
            uniq[chrom][end] = (uniq[chrom][end]+1) if end in uniq[chrom] else 1
    return uniq, multi

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='count the number of reads map to each window based on BWA mapped BAM files, need to install samtools software.')
    parser.add_argument('-b', type=str, help='BWA mapped bam file to be perocessed.', required=True)
    parser.add_argument('-r', type=str, help='reference DB file.', required=True)
    parser.add_argument('-w', type=int, help='window size (default: 1000)', default=1000)
    parser.add_argument('-s', type=int, help='step size (default: 500)', default=500)
    args = parser.parse_args()
    print(f'Window size: {args.w}bp')
    print(f'Step size: {args.s}bp')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': Extract Sample name...')
    sample = getSampleName(bamfile=args.b)
    print('***** Raw reads count Start! *****')
    uniq, multi = getRawCount(bamfile=args.b)
    with open(args.r, 'r') as fr:
        with open(sample+'_raw', 'w') as fw:
            fw.write(f'#Chr\tindex\tpos\tuniq\tmul\tGC\trepeat\n')
            for line in fr.readlines():
                genomeDB = line.strip().split('\t')
                multi_ = multi.get(genomeDB[0],{}).get(int(genomeDB[2]), 0)
                uniq_ = uniq.get(genomeDB[0],{}).get(int(genomeDB[2]), 0)
                fw.write(f'{genomeDB[0]}\t{genomeDB[1]}\t{genomeDB[2]}\t{uniq_}\t{multi_}\t{genomeDB[3]}\t{genomeDB[4]}\n')
    print('***** Raw reads count Done! *****')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()))
    
