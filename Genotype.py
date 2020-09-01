#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   Genotype.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import time
import argparse
import numpy as np
import pandas as pd
from itertools import chain
from multiprocessing import Pool
from collections import Counter
from sklearn.mixture import BayesianGaussianMixture
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score

def loadCNVRFile(cnvrFile):
    df = pd.read_table(cnvrFile, low_memory=False, header=0)
    info_df = df.iloc[:, :8]
    cnv_df = df.iloc[:, 8:-2]
    samples = cnv_df.columns
    avg_ay = df['average'].values
    sd_ay = df['sd'].values
    return info_df, cnv_df.values, list(samples), avg_ay, sd_ay

def genotype(cnvays):
    result = []
    n_com = 10 if cnvays.shape[1] > 10 else cnvays.shape[1]-1
    n_init = 3
    for cnvay in cnvays:
        cnv = [[x] for x in cnvay]
        gmm = GaussianMixture(n_components=n_com, n_init=n_init, max_iter=10000,
                              tol=1e-3, covariance_type='full').fit(cnv)
        labels = gmm.predict(cnv)
        normed_ay = np.arange(0, np.max(cnvay)+0.5, 0.5)
        swlables = {}
        for rawlabel in np.unique(labels):
            swlables[rawlabel] = normed_ay[np.argmin(np.abs(normed_ay - np.median(cnvay[labels == rawlabel])))]
        newlabels = [swlables[x] for x in labels]
        gtlabels = {0: 'dd', 0.5: 'Ad', 1: 'AA', 1.5: 'AB', 2: 'BB', 2.5: 'BC'}
        finalline = [gtlabels.get(x, 'M') for x in newlabels]
        if len(np.unique(finalline)) > 1:
            sc = silhouette_score(cnv, finalline, metric='euclidean')
            chs = calinski_harabasz_score(cnv, finalline)
        else:
            sc = np.nan
            chs = np.nan
        llh = gmm.score(cnv)
        finalline += [sc, chs, llh]
        result.append(finalline)
    return result

def stat_gt(ay):
    gts = ['dd', 'Ad', 'AA', 'AB', 'BB', 'BC', 'M']
    c = Counter(ay)
    counts = []
    for gt in gts:
        counts.append(c.get('gt', 0))
    return counts

def produce_vcf(df, cnvays, samples, outname):
    gt2alt = {
        'd': 'CN0',
        'A': 'CN1',
        'B': 'CN2',
        'C': 'CNH',
        'M': 'CNH'
    }
    with open(outname, 'w') as f:
        f.write('''##fileformat=VCFv4.2
##ALT=<ID=CN0,Description="Copy number allele: 0 copies">
##ALT=<ID=CN1,Description="Copy number allele: 1 copy">
##ALT=<ID=CN2,Description="Copy number allele: 2 copies">
##ALT=<ID=CNH,Description="Copy number allele: more than 2 copies">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CP,Number=1,Type=Float,Description="Normalized copy number">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SILHOUETTESCORE,Number=1,Type=Float,Description="silhouette score of genotype on this cnvr">
##INFO=<ID=CALINSKIHARABAZESCORE,Number=1,Type=Float,Description="calinski harabaz score of genotype on this cnvr">
##INFO=<ID=LOGLIKELIHOOD,Number=1,Type=Float,Description="log likelihood of genotype on this cnvr">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT''')
        for sm in samples:
            f.write(f'\t{sm}')
        f.write('\n')
        for row, cns in zip(df.values, cnvays):
            cns *= 2
            svtype = 'DEL' if np.mean(cns)<2 else 'DUP'
            chrom = row[0]
            start = row[1]
            end = row[2]
            gts = row[8: 8+len(samples)]
            allels = list(chain.from_iterable([x if x!= 'M' else 'MM' for x in gts]))
            allels = [gt2alt[x] for x in allels]
            alt = sorted(list(set(allels)))
            try:
                alt.remove('CN1')
            except Exception:
                pass
            allelsw = {'CN1': '0'}
            for index, gt in enumerate(alt, 1):
                allelsw[gt] = str(index)
            gts = []
            for a,b in zip(allels[::2], allels[1::2]):
                gts.append(f'{allelsw[a]}/{allelsw[b]}')
            fgts = []
            for gt, cn in zip(gts, cns):
                fgts.append(f'{gt}:{cn}')
            fgts = '\t'.join(fgts)
            alt = ','.join(alt) if alt else '.'
            silhouette_score = row[-12]
            calinski_harabasz_score = row[-11]
            log_likelihood = row[-10]
            f.write(f'{chrom}\t{start}\t{f"{chrom}:{start}-{end}"}\tA\t{alt}\t.\t.\tEND={end};SVTYPE={svtype};SILHOUETTESCORE={silhouette_score};CALINSKIHARABASZSCORE={calinski_harabasz_score};LOGLIKELIHOOD={log_likelihood}\tGT:CP\t{fgts}\n')

def produce_mergevcf(df, cnvays, samples, outname):
    gt2alt = {
        'd': 'CN0',
        'A': 'CN1',
        'B': 'CN2',
        'C': 'CN2',
        'M': 'CN2'
    }
    with open(outname, 'w') as f:
        f.write('''##fileformat=VCFv4.2
##ALT=<ID=CN0,Description="Copy number allele: 0 copies">
##ALT=<ID=CN1,Description="Copy number allele: 1 copy">
##ALT=<ID=CN2,Description="Copy number allele: 2 copies">
##ALT=<ID=CNH,Description="Copy number allele: more than 2 copies">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CP,Number=1,Type=Float,Description="Normalized copy number">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SILHOUETTESCORE,Number=1,Type=Float,Description="silhouette score of genotype on this cnvr">
##INFO=<ID=CALINSKIHARABAZESCORE,Number=1,Type=Float,Description="calinski harabaz score of genotype on this cnvr">
##INFO=<ID=LOGLIKELIHOOD,Number=1,Type=Float,Description="log likelihood of genotype on this cnvr">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT''')
        for sm in samples:
            f.write(f'\t{sm}')
        f.write('\n')
        for row, cns in zip(df.values, cnvays):
            svtype = 'DEL' if np.mean(cns)<2 else 'DUP'
            chrom = row[0]
            start = row[1]
            end = row[2]
            gts = row[8: 8+len(samples)]
            allels = list(chain.from_iterable([x if x!= 'M' else 'MM' for x in gts]))
            allels = [gt2alt[x] for x in allels]
            alt = sorted(list(set(allels)))
            try:
                alt.remove('CN1')
            except Exception:
                pass
            allelsw = {'CN1': '0'}
            for index, gt in enumerate(alt, 1):
                allelsw[gt] = str(index)
            gts = []
            for a,b in zip(allels[::2], allels[1::2]):
                gts.append(f'{allelsw[a]}/{allelsw[b]}')
            fgts = []
            for gt, cn in zip(gts, cns):
                fgts.append(f'{gt}:{cn}')
            fgts = '\t'.join(fgts)
            alt = ','.join(alt) if alt else '.'
            silhouette_score = row[-12]
            calinski_harabasz_score = row[-11]
            log_likelihood = row[-10]
            f.write(f'{chrom}\t{start}\t{f"{chrom}:{start}-{end}"}\tA\t{alt}\t.\t.\tEND={end};SVTYPE={svtype};SILHOUETTESCORE={silhouette_score};CALINSKIHARABASZSCORE={calinski_harabasz_score};LOGLIKELIHOOD={log_likelihood}\tGT:CP\t{fgts}\n')

def main(cnvfile, outprefix, merge, nproc):
    info_df, cnvays, samples, avg_ay, sd_ay = loadCNVRFile(cnvfile)
    sample_header = samples + ['silhouette_score', 'calinski_harabaz_score', 'Log_likelihood']
    if nproc > 1:
        pool = Pool(processes=nproc)
        gtdf = pd.DataFrame(np.vstack(pool.map(genotype, np.array_split(cnvays, nproc))))
    else:
        gtdf = pd.DataFrame(genotype(cnvays))
    gtdf.columns = sample_header
    df = pd.concat([info_df, gtdf], axis=1, join='inner')
    df['average'] = avg_ay
    df['sd'] = sd_ay
    df = pd.concat([df,
                    pd.DataFrame(df[samples].apply(stat_gt, axis=1).values.tolist(), columns=['dd', 'Ad', 'AA', 'AB', 'BB', 'BC', 'M'])],
                    axis=1)
    df.to_csv(f"{outprefix}.tsv", sep='\t', index=False, na_rep='NA')
    produce_vcf(df, cnvays, samples, f"{outprefix}.vcf")
    if merge:
        produce_mergevcf(df, cnvays, samples, f"{outprefix}_merge.vcf")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Clustering the input samples into genotypes uses Gaussian mixture modes.')
    parser.add_argument('-i', type=str, help='input cnvr file', required=True)
    parser.add_argument('-o', type=str, help='output prefix of genotyped file', required=True)
    parser.add_argument('--merge', type=bool, help='merge all duplication types to one type (default: False)' ,default=False)
    parser.add_argument('-t', type=int, help='number of process will be used (default: 1)', default=1)
    args = parser.parse_args()
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), 'Genotyping...')
    main(cnvfile=args.i, outprefix=args.o, merge=args.merge, nproc=args.t)
    print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()), 'Genotype Done!')