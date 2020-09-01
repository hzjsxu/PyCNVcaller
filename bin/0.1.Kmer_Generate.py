#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   0.1.Kmer_Generate.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import argparse
import time
from collections import OrderedDict

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

def main(fafile, outfile, windowSize=1000, stepSize=500):
    faSeqdict = readFa(fafile)
    with open(outfile, 'w') as f:
        for chrom, seq in faSeqdict.items():
            for start in range(0, len(seq), stepSize):
                f.write(f'>{chrom}_{start+1}\n{seq[start:start+windowSize]}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='produce kmer(<OUTFILE>) from refrence fasta(<FAFILE>) according to window size(<WINSIZE>).')
    parser.add_argument('-r', type=str, help='Input an reference fasta file.', required=True)
    parser.add_argument('-o', type=str, help='output file generated.', required=True)
    parser.add_argument('-w', type=int, help='window size (default: 1000)', default=1000)
    parser.add_argument('-s', type=int, help='step size (default: 500)', default=500)
    args = parser.parse_args()
    print(f'Window size: {args.w}bp')
    print(f'Step size: {args.s}bp')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': Generating kmer sequences for each window...')
    main(fafile=args.r, outfile=args.o)
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': ===== Done! ======')