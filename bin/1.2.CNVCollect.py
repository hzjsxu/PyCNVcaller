#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   1.2.CNVCollect.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import argparse
import time

def loadRaw(rawfile):
    rawDict = {}
    with open(rawfile, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line[0] not in rawDict:
                rawDict[line[0]] = {}
            rawDict[line[0]][line[2]] = int(line[3])+int(line[4])
    print(f'Raw reads count file loading finished!')
    return rawDict

def loadLink(linkfile):
    link = {}
    with open(linkfile, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chrom, pos = line[0].split(":")
            if chrom not in link:
                link[chrom] = {}
            link[chrom][pos] = line[1:]
    print(f'duplicated window record file loading finished!')
    return link

def main(rawfile, outfile):
    correct = 0
    with open(rawfile, 'r') as fr:
        with open(outfile, 'w') as fw:
            for line in fr:
                if line.startswith('#'):
                        continue
                line = line.strip().split('\t')
                absolute_cp = 0
                genome_copy = 0           
                if line[0] in link:
                    if line[2] in link[line[0]]:
                        correct += 1
                        for pos in link[line[0]][line[2]]:
                            genome_copy += 1
                            tmp_chr, tmp_pos = pos.split(':')
                            if tmp_pos in rawDict[tmp_chr]:
                                absolute_cp += rawDict[tmp_chr][tmp_pos]
                            else:
                                absolute_cp += 0
                        fw.write(f'{line[0]}\t{line[1]}\t{absolute_cp}\t{line[5]}\t{line[6]}\t{genome_copy}\n')
                    else:
                        fw.write(f'{line[0]}\t{line[1]}\t{int(line[3])+int(line[4])}\t{line[5]}\t{line[6]}\t1\n')
    return correct

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sum the raw reads count from the high similarity windows record in the duplicated window record file.')
    parser.add_argument('-raw', type=str, help='Raw reads count file.', required=True)
    parser.add_argument('-link', type=str, help='Duplicated window record file.', required=True)
    parser.add_argument('-o', type=str, help='output file.', required=True)
    args = parser.parse_args()
    rawDict = loadRaw(args.raw)
    link = loadLink(args.link)
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': Absolute reads counts correction...')
    correct = main(rawfile=args.raw, outfile=args.o)
    print(f'corrected window number: {correct}')
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), ': ===== Correction Done! ======')