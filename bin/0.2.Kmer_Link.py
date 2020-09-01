#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   0.2.Kmer_Link.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

## eg. 
# 1. 建立索引文件 
#   sawriter ~/apple/GDDH13_1-1_formatted.fasta  得到sa索引文件
# 2. 比对
# blasr kmer.fa ~/apple/GDDH13_1-1_formatted.fasta --sa ~/apple/GDDH13_1-1_formatted.fasta.sa \
#               --out kmer.aln -m 5 --noSplitSubreads --minMatch 15 --maxMatch 20 --advanceHalf \
#               --advanceExactMatches 10 --fastMaxInterval --fastSDP --aggressiveIntervalCut --bestn 10 \
#               --nproc 10

# here put the import lib
import argparse
import time
import pandas as pd
import numpy as np

def rename(seqName):
    chrom= "_".join(seqName.split("_")[:-1])
    pos = int(seqName.split("_")[-1])
    return f'{chrom}:{pos}'

def main(blasrFile, outfile, windowSzie=1000, stepSize=500):
    col_names = ['qName', 'qLength', 'qStart', 'qEnd', 'tName', 'tStart', 'numMatch', 'numMismatch', 'numIns', 'numDel']
    col_index = [0, 1, 2, 3, 5, 7, 11, 12, 13, 14]
    df = pd.read_csv(blasrFile, sep='\s+', header=None, usecols=col_index, names=col_names, compression='infer')
    df['qName'] = df['qName'].map(rename)
    print(f'Parse blasr finished, records number is {df.shape[0]}.')
    df['coverage'] = (df['qEnd'] - df['qStart'])/df['qLength']
    df = df.loc[df['coverage']>=0.9, ['qName', 'qStart', 'tName', 'tStart', 'numMatch', 'numMismatch', 'numIns', 'numDel']]
    print(f'Drop coverage lower than 90%, {df.shape[0]} records remained.')
    df['identity'] = df['numMatch']/(df['numMatch'] + df['numMismatch'] + df['numIns'] + df['numDel'])
    df = df.loc[df['identity']>=0.97, ['qName', 'tName', 'tStart']]
    print(f'Drop identity lower than 97%, {df.shape[0]} records remained.')
    df['tNewstart'] = df['tStart'].map(lambda x: round(x/stepSize)*stepSize+1)
    df['tNewName'] = df['tName'].astype(str) + ':' + df['tNewstart'].astype(str)
    df = df.groupby('qName')[['qName', 'tNewName']].agg({'qName': np.unique,
                                                       'tNewName': lambda x: "\t".join(x)})
    
    with open(outfile, 'w') as fw:
        for line in df.values:
            if len(line[1].split()) > 1:
                outline = [str(x) for x in line]
                fw.write('\t'.join(outline) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='produce duplicated window record file from blasr outputfile.')
    parser.add_argument('-b', type=str, help='The blasr output file (m=5).', required=True)
    parser.add_argument('-o', type=str, help='output file generated.', required=True)
    parser.add_argument('-w', type=int, help='window size (default: 1000)', default=1000)
    parser.add_argument('-s', type=int, help='step size (default: 500)', default=500)
    args = parser.parse_args()
    print(f'Window size: {args.w}bp\tStep size: {args.s}bp')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': Generating duplicated window record file...')
    main(blasrFile=args.b, outfile=args.o, windowSzie=args.w, stepSize=args.s)
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': ===== Done! ======')
    