#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   CNVReferenceDB.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import os
import re
import argparse
import time
from collections import OrderedDict, Counter

## eg. python CNVReferenceDB.py -i ~/apple/GDDH13_1-1_formatted.fasta -o testDB_1000_500 -s 500 -w 1000

## 读取fasta文件，以字典形式存储
def readFa(fafile):
    faSeqDict = OrderedDict()
    with open(fafile, 'r')  as f:
        chrom = f.readline().strip().split()[0][1:]
        seq = []
        for line in f:
            if line[0] != '>':
                seq.append(line.strip())
            else:
                faSeqDict[chrom] = ''.join(seq)
                chrom = line.strip().split()[0][1:]
                seq = []
        faSeqDict[chrom] = ''.join(seq)
    return faSeqDict


def main(fafile, outfile, windowsize=1000, stepsize=500, lower_GC=0.2, upper_GC=0.7, gap=0.5):
    faSeqDict = readFa(fafile)
    lower_GC_count = lower_GC * windowsize
    upper_GC_count = upper_GC * windowsize
    with open(outfile, 'w') as fw:
        for chrom in sorted(faSeqDict):
            winNumber = int(len(faSeqDict[chrom])/stepsize)
            seq = faSeqDict[chrom]
            i = 0
            for num in range(winNumber):
                position = num * stepsize
                subseq = seq[position:(position+windowsize)]
                GC_count = subseq.count('C') + subseq.count('G')
                if GC_count < lower_GC_count or GC_count > upper_GC_count:
                    continue
                repeat_count = subseq.count('a') + subseq.count('g') + subseq.count('c') + subseq.count('t')
                gap_content = round(subseq.count('N')/len(subseq), 2) if len(subseq)!=0 else 0
                if gap_content > 0.5:
                    continue
                i += 1
                fw.write(f'{chrom}\t{i}\t{position+1}\t{GC_count}\t{repeat_count}\t{gap_content}\n')

   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Format reference genome into state table of GC and repeat content per window, as removing gapped regions.')
    parser.add_argument('-i', type=str, help='Input an reference fasta file.', required=True)
    parser.add_argument('-w', type=int, help='window size (default: 1000)', default=1000)
    parser.add_argument('-s', type=int, help='step size (default: 500)', default=500)
    parser.add_argument('-l', type=float, help='GC contnet lower limit (default: 0.2)', default=0.2)
    parser.add_argument('-u', type=float, help='GC content upper limit (default: 0.7)', default=0.7)
    parser.add_argument('-g', type=float, help='gap content limit (default: 0.5)', default=0.5)
    args = parser.parse_args()
    print(f'Window size: {args.w}bp')
    print(f'Step size: {args.s}bp')
    print('***** Genome format Start! *****')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()))
    outfile = 'referenceDB.' + str(args.w) + '_' + str(args.s)
    main(fafile=args.i, outfile=outfile, windowsize=args.w, stepsize=args.s, lower_GC=args.l, upper_GC=args.u, gap=args.g)
    print('***** Genome format Done! *****')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()))
    
        