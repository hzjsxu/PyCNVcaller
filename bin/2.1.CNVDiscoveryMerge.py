#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   2.1.CNVDiscoveryMerge.py
@Author  :   jsxu 
@Contact :   hzaujsxu@163.com
'''

# here put the import lib
import pandas as pd
import numpy as np
import re
import argparse

def getFileterSample(excludelist):
	with open(excludelist, 'r') as fr:
		filter_sample = {}
		for line in fr:
			line = line.strip()
			if line:
				filter_sample[line] = 1
	return filter_sample

def getMedian75(normalizedFile, filter_sample, preHetCutoff):
	with open(normalizedFile, 'r') as fr:
		sv_array = []
		for line in fr:
			line = line.strip()
			if line and line not in filter_sample:
				average, sd = map(float, re.findall('(\d+\.?\d+)_SD_(\d+\.?\d+)', line)[0])
				sv = preHetCutoff*average/sd if sd>0 else 0
				sv_array.append(sv)
		sv_array = sorted(sv_array, reverse=True)
		median75 = sv_array[int(0.75*len(sv_array))]
	return median75

#####################################################
def getr(r, sample_num):
	if r:
		r = r
	elif sample_num <=30:
		r = 0.5
	elif 30 < sample_num <= 50:
		r = 0.4
	elif 50 < sample_num <= 100:
		r = 0.3
	elif 100 < sample_num <= 200:
		r = 0.2
	elif 200 < sample_num <= 500:
		r = 0.15
	else:
		r = 0.1
	print(f'pearson correlation: {r}')
	return r

#####################################################
def getInfo(normalizedFile, filter_sample, prehetCutoff, preHomoCutoff):
	hash_res_gain = {}
	hash_res_loss = {}
	hash_res_homo_gain = {}
	hash_res_homo_loss = {}
	genome_copy = {}
	effective_sample = 0
	sample_num = 0
	with open(normalizedFile, 'r') as fr:
		for sample in fr:
			sample = sample.strip()
			remove = 0
			if sample:
				sample_num += 1
				if sample in filter_sample:
					remove = 1
				if remove == 0:
					effective_sample += 1
				average, sd = map(float, re.findall('(\d+\.?\d+)_SD_(\d+\.?\d+)', sample)[0])
				sv = prehetCutoff*average/sd if sd>0 else 0
				hetcutoff = round(median75*sd/average, 2) if sv<=median75 else prehetCutoff
				print(f'{sample} hetcutoff: {hetcutoff}')
				with open(sample, 'r') as fp:
					for line in fp:
						line = line.strip().split('\t')
						if line[0] not in genome_copy:
							genome_copy[line[0]] = {}
						if line[1] not in genome_copy[line[0]]:
							genome_copy[line[0]][line[1]] = line[3]
						if remove == 1:
							if float(line[2])>=1+hetcutoff:
								if line[0] not in hash_res_gain:
									hash_res_gain[line[0]] = {}
								if line[1] not in hash_res_gain[line[0]]:
									hash_res_gain[line[0]][line[1]] = 0
								hash_res_gain[line[0]][line[1]] += 1
								if float(line[2])>=1+preHomoCutoff:
									if line[0] not in hash_res_homo_gain:
										hash_res_homo_gain[line[0]] = {}
									if line[1] not in hash_res_homo_gain[line[0]]:
										hash_res_homo_gain[line[0]][line[1]] = 0
									hash_res_homo_gain[line[0]][line[1]] += 1
							elif float(line[2])<=1-hetcutoff:
								if line[0] not in hash_res_loss:
									hash_res_loss[line[0]] = {}
								if line[1] not in hash_res_loss[line[0]]:
									hash_res_loss[line[0]][line[1]] = 0
								hash_res_loss[line[0]][line[1]] += 1
								if float(line[2])<=1-preHomoCutoff:
									if line[0] not in hash_res_homo_loss:
										hash_res_homo_loss[line[0]] = {}
									if line[1] not in hash_res_homo_loss[line[0]]:
										hash_res_homo_loss[line[0]][line[1]] = 0
									hash_res_homo_loss[line[0]][line[1]] += 1
							else:
								pass
	return hash_res_gain, hash_res_loss, hash_res_homo_gain, hash_res_homo_loss, genome_copy, effective_sample
##############################################################################################

def effective(refdbFile, f, m, hash_res_gain, hash_res_loss, hash_res_homo_gain, hash_res_homo_loss, effective_sample):
	refdb = {}
	gain_effective = {}
	loss_effective = {}
	with open(refdbFile, 'r') as fp:
		for record in fp:
			line = record.strip()
			if line:
				line = line.split('\t')
				if line[0] not in refdb:
					refdb[line[0]] = {}
				refdb[line[0]][int(line[1])] = record
				if line[0] not in gain_effective:
					gain_effective[line[0]] = {}
				if (line[1] in hash_res_gain[line[0]] and (hash_res_gain[line[0]][line[1]]/effective_sample) >= f) or (line[1] in hash_res_homo_gain[line[0]] and hash_res_homo_gain[line[0]][line[1]]>=m):
					hash_res_gain[line[0]][line[1]] = 1
					gain_effective[line[0]][line[1]] = 1
				else: 
					hash_res_gain[line[0]][line[1]] = 0
				if line[0] not in loss_effective:
					loss_effective[line[0]] = {}
				if (line[1] in hash_res_loss[line[0]] and (hash_res_loss[line[0]][line[1]]/effective_sample) >= f) or (line[1] in hash_res_homo_loss[line[0]] and hash_res_homo_loss[line[0]][line[1]]>=m):
					hash_res_loss[line[0]][line[1]] = 1
					loss_effective[line[0]][line[1]] = 1
				else: 
					hash_res_loss[line[0]][line[1]] = 0
	return refdb, gain_effective, loss_effective, hash_res_gain, hash_res_loss

#####################  merged cnv windows (as >=3 Yes windows in 5 windows) into Regions  #################
def getCNVregion(refdb, has_res_gain, has_res_loss):
	CNV_region_gain = {}
	CNV_region_loss = {}
	for tmp_chr in sorted(refdb.keys()):
		for i in sorted(refdb[tmp_chr].keys()):
			if i>2 and (str(i+1) in hash_res_gain[tmp_chr]) and (str(i+2) in hash_res_gain[tmp_chr]) and (str(i+1) in hash_res_loss[tmp_chr]) and (str(i+2) in hash_res_loss[tmp_chr]):
				if (hash_res_gain[tmp_chr][str(i-2)]+hash_res_gain[tmp_chr][str(i-1)]+hash_res_gain[tmp_chr][str(i)]+hash_res_gain[tmp_chr][str(i+1)]+hash_res_gain[tmp_chr][str(i+2)]) >= 3:
					if tmp_chr not in CNV_region_gain:
						CNV_region_gain[tmp_chr] = {}
					CNV_region_gain[tmp_chr][str(i)] = 'Y'
				if (hash_res_loss[tmp_chr][str(i-2)]+hash_res_loss[tmp_chr][str(i-1)]+hash_res_loss[tmp_chr][str(i)]+hash_res_loss[tmp_chr][str(i+1)]+hash_res_loss[tmp_chr][str(i+2)]) >= 3:
					if tmp_chr not in CNV_region_loss:
						CNV_region_loss[tmp_chr] = {}
					CNV_region_loss[tmp_chr][str(i)] = 'Y'
	return CNV_region_gain, CNV_region_loss

##############################################################################
def getNum(normalizedFile, CNV_region_gain, CNV_region_loss):
	hash_num = {}
	with open(normalizedFile, 'r') as fr:
		for sample in fr:
			sample = sample.strip()
			if sample:
				with open(sample, 'r') as fp:
					for line in fp:
						line = line.split('\t')
						if line[0] not in hash_num:
							hash_num[line[0]] = {}
						if (line[1] in CNV_region_gain[line[0]]) or (line[1] in CNV_region_loss[line[0]]):
							if line[1] not in hash_num[line[0]]:
								hash_num[line[0]][line[1]] = ''
							hash_num[line[0]][line[1]] += line[2]+'\t'   
	return hash_num

##############################################################################
def getTmpCNV(refdb, CNV_region_gain, gain_effective, CNV_region_loss, loss_effective):
	cnv_tmp_gain = {}
	cnv_tmp_loss = {}
	for tmp_chr in sorted(refdb.keys()):
		gain_index = gain_index2 = loss_index = loss_index2 = 0
		for i in sorted(refdb[tmp_chr].keys()):
			if str(i) in CNV_region_gain[tmp_chr]:
				gain_index += 1
				if gain_index == 1:
					tmp_gain_start = i
				if str(i) in gain_effective[tmp_chr]:
					gain_index2 += 1
				tmp_gain_end = i
			else:
				if gain_index2 >= 3:
					if tmp_chr not in cnv_tmp_gain:
						cnv_tmp_gain[tmp_chr] = {}
					cnv_tmp_gain[tmp_chr][tmp_gain_start] = tmp_gain_end
				gain_index = gain_index2 = 0
			
			if str(i) in CNV_region_loss[tmp_chr]:
				loss_index += 1
				if loss_index == 1:
					tmp_loss_start = i
				if str(i) in loss_effective[tmp_chr]:
					loss_index2 += 1
				tmp_loss_end = i
			else:
				if loss_index2 >= 3:
					if tmp_chr not in cnv_tmp_loss:
						cnv_tmp_loss[tmp_chr] = {}
					cnv_tmp_loss[tmp_chr][tmp_loss_start] = tmp_loss_end
				loss_index = loss_index2 = 0
		if gain_index2 >= 3:
			cnv_tmp_gain[tmp_chr][tmp_gain_start] = tmp_gain_end
		if loss_index2 >= 3:
			cnv_tmp_loss[tmp_chr][tmp_loss_start] = tmp_loss_end
	print("define effective window finished!")
	return cnv_tmp_gain, cnv_tmp_loss

##########################  CNV region refine  ##############################
def getGain(cnv_tmp_gain, gain_effective, hash_num,  r):
	cnv_gain = {}
	for tmp_chr in sorted(cnv_tmp_gain.keys()):
		tmp_effective = list(map(int, gain_effective[tmp_chr].keys()))
		for tmp_index in sorted(map(int, cnv_tmp_gain[tmp_chr].keys())):
			tmp_selected = [x for x in tmp_effective if x>=tmp_index and x<=cnv_tmp_gain[tmp_chr][tmp_index]]
			for tmp in range(len(tmp_selected)-2):
				tmp_cn_a = hash_num[tmp_chr][str(tmp_selected[tmp])].split('\t')
				tmp_cn_b = hash_num[tmp_chr][str(tmp_selected[tmp+2])].split('\t')
				tmp_cn_a = pd.Series([float(x) for x in tmp_cn_a if x])
				tmp_cn_b = pd.Series([float(x) for x in tmp_cn_b if x])
				tmp_correlation = tmp_cn_a.corr(tmp_cn_b)
				if tmp_chr not in cnv_gain:
					cnv_gain[tmp_chr] = {}
				if tmp_correlation >= r:
					cnv_gain[tmp_chr][tmp_selected[tmp]] = 1
					cnv_gain[tmp_chr][tmp_selected[tmp+2]] = 1
				else:
					if tmp_selected[tmp] in cnv_gain[tmp_chr]:
						del cnv_gain[tmp_chr][tmp_selected[tmp]]
	return cnv_gain

def getLoss(cnv_tmp_loss, loss_effective, hash_num, r):
	cnv_loss = {}
	for tmp_chr in sorted(cnv_tmp_loss.keys()):
		tmp_effective = list(map(int, loss_effective[tmp_chr].keys()))
		for tmp_index in sorted(map(int, cnv_tmp_loss[tmp_chr].keys())):
			tmp_selected = [x for x in tmp_effective if x>=tmp_index and x<=cnv_tmp_loss[tmp_chr][tmp_index]]
			for tmp in range(len(tmp_selected)-2):
				tmp_cn_a = hash_num[tmp_chr][str(tmp_selected[tmp])].split('\t')
				tmp_cn_b = hash_num[tmp_chr][str(tmp_selected[tmp+2])].split('\t')
				tmp_cn_a = [float(x) for x in tmp_cn_a if x]
				tmp_cn_b = [float(x) for x in tmp_cn_b if x]
				tmp_cn_a = pd.Series(tmp_cn_a)
				tmp_cn_b = pd.Series(tmp_cn_b)
				tmp_correlation = tmp_cn_a.corr(tmp_cn_b)
				if tmp_chr not in cnv_loss:
					cnv_loss[tmp_chr] = {}
				if tmp_correlation >= r:
					cnv_loss[tmp_chr][tmp_selected[tmp]] = 1
					cnv_loss[tmp_chr][tmp_selected[tmp+2]] = 1
				else:
					if tmp_selected[tmp] in cnv_loss[tmp_chr]:
						del cnv_loss[tmp_chr][tmp_selected[tmp]]
	return cnv_loss

##################################################################################
def getCNV(refdb, cnv_gain, cnv_loss, gain_effective, loss_effective):
	cnv = {}
	for tmp_chr in sorted(refdb.keys()):
		if tmp_chr not in cnv:
			cnv[tmp_chr] = {}
		gain_index = gain_index2 = loss_index = loss_index2 = 0
		for i in sorted(refdb[tmp_chr].keys()):
			if i in cnv_gain[tmp_chr]:
				gain_index += 1
				if gain_index == 1:
					tmp_gain_start = i
				if str(i) in gain_effective[tmp_chr]:
					gain_index2 += 1
				tmp_gain_end = i
			else:
				if gain_index2 >= 3:
					cnv[tmp_chr][tmp_gain_start] = tmp_gain_end
				gain_index = gain_index2 = 0
			if i in cnv_loss[tmp_chr]:
				loss_index += 1
				if loss_index == 1:
					tmp_loss_start = i
				if str(i) in loss_effective[tmp_chr]:
					loss_index2 += 1
				tmp_loss_end = i
			else:
				if loss_index >= 3:
					cnv[tmp_chr][tmp_loss_start] = tmp_loss_end
				loss_index = loss_index2 = 0
		if gain_index2 >= 3:
			cnv[tmp_chr][tmp_gain_start] = tmp_gain_end
		if loss_index2 >= 3:
			cnv[tmp_chr][tmp_loss_start] = tmp_loss_end
	print('Refine effective window finished!')
	return cnv
##################### compute the GC ratio, repeat ratio, Duplication per CNV region ###################
def info(cnv, gain_effective, loss_effective, genome_copy, windowsize):
	cnv_GC = {}
	cnv_repeat = {}
	cnv_gap = {}
	cnv_real_start = {}
	cnv_real_end = {}
	kmer_median = {}

	for tmp_chr in cnv.keys():
		if tmp_chr not in cnv_real_start:
			cnv_real_start[tmp_chr] = {}
		if tmp_chr not in cnv_real_end:
			cnv_real_end[tmp_chr] = {}
		if tmp_chr not in cnv_repeat:
			cnv_repeat[tmp_chr] = {}
		if tmp_chr not in cnv_gap:
			cnv_gap[tmp_chr] = {}
		if tmp_chr not in cnv_GC:
			cnv_GC[tmp_chr] = {}
		if tmp_chr not in kmer_median:
			kmer_median[tmp_chr] = {}
		for start_index in cnv[tmp_chr].keys():
			cnv_real_start[tmp_chr][start_index] = refdb[tmp_chr][start_index].split('\t')[2]
			tmp_median_kmer = []
			cnv_GC[tmp_chr][start_index] = 0
			cnv_repeat[tmp_chr][start_index] = 0
			cnv_gap[tmp_chr][start_index] = 0
			for tmp_index in range(start_index, cnv[tmp_chr][start_index]+1):
				if str(tmp_index) in gain_effective[tmp_chr] or str(tmp_index) in loss_effective[tmp_chr]:
					cnv_GC[tmp_chr][start_index] += int(refdb[tmp_chr][start_index].split('\t')[3])
					cnv_repeat[tmp_chr][start_index] += int(refdb[tmp_chr][start_index].split('\t')[4])
					cnv_gap[tmp_chr][start_index] += float(refdb[tmp_chr][start_index].split('\t')[5])
					if genome_copy[tmp_chr][str(tmp_index)]:
						tmp_median_kmer.append(genome_copy[tmp_chr][str(tmp_index)])
			if cnv[tmp_chr][start_index] not in refdb[tmp_chr]:
				continue
			cnv_real_end[tmp_chr][start_index] = int(refdb[tmp_chr][cnv[tmp_chr][start_index]].split('\t')[2])+windowsize-1
			tmp_median_kmer = sorted(map(int, tmp_median_kmer))
			# print(tmp_median_kmer)
			if tmp_median_kmer:
				kmer_median[tmp_chr][start_index] = tmp_median_kmer[int(0.5*len(tmp_median_kmer))]
	return cnv_GC, cnv_repeat, cnv_gap, cnv_real_start, cnv_real_end, kmer_median 

############## compute corrected coverage per sample per cnv region based on OUT file, using median50  ###########
def main(normalizedFile, outfile, cnv, gain_effective, loss_effective, cnv_repeat, cnv_GC, cnv_gap, cnv_real_start, cnv_real_end, kmer_median, windowsize):
	tmp_rd_sample = {}
	with open(outfile, 'w') as fw:
		fw.write(f'chr\tstart\tend\tnumber\tgap\trepeat\tgc\tkmer')
		with open(normalizedFile, 'r') as fr:
			tmp_file = [line.strip() for line in fr if line.strip()]
			sample = len(tmp_file)
			for num in range(sample):
				tmp_file_name = tmp_file[num]
				suffix = tmp_file_name.split('/')[-1]
				tmp = suffix.split('_')[0]
				fw.write(f'\t{tmp}')
		fw.write('\taverage\tsd\n')

		cnv_keys = sorted(cnv.keys())
		for tmp_chr in cnv_keys:
			for start_index in sorted(cnv[tmp_chr].keys()):
				tmp_effective_wind_num = 0
				for tmp_index in range(start_index, cnv[tmp_chr][start_index]+1):
					if str(tmp_index) in gain_effective[tmp_chr] or str(tmp_index) in loss_effective[tmp_chr]:
						tmp_effective_wind_num += 1

						for num in range(sample):
							tmp_rd_sample[num] = []
							if str(tmp_index) in hash_num[tmp_chr]:
								tmp_rd_sample[num].append(hash_num[tmp_chr][str(tmp_index)].split('\t')[num])
				repeat_ration = round(cnv_repeat[tmp_chr][start_index]/(tmp_effective_wind_num*windowsize),2) if tmp_effective_wind_num>0 else 0
				gc_ration = round(cnv_GC[tmp_chr][start_index]/(tmp_effective_wind_num*windowsize),2) if tmp_effective_wind_num>0 else 0
				gap_ration = round(cnv_gap[tmp_chr][start_index]/tmp_effective_wind_num, 2) if tmp_effective_wind_num>0 else 0
				if start_index in cnv_real_start[tmp_chr] and start_index in cnv_real_end[tmp_chr] and start_index in kmer_median[tmp_chr]:
					fw.write(f'{tmp_chr}\t{cnv_real_start[tmp_chr][start_index]}\t{cnv_real_end[tmp_chr][start_index]}\t{tmp_effective_wind_num}\t{gap_ration}\t{repeat_ration}\t{gc_ration}\t{kmer_median[tmp_chr][start_index]}')
				# print(hash_num[tmp_chr][str(tmp_index)])
				# print(tmp_rd_sample)
				cp_median_CNVR = []
				for num in range(sample):
					sample_CNVR_cp = sorted((float(x) for x in tmp_rd_sample[num] if x))
					cp_median = sample_CNVR_cp[int(0.5*len(sample_CNVR_cp))] if len(sample_CNVR_cp) else 0
					cp_median_CNVR.append(cp_median)
					fw.write(f'\t{cp_median}')
				average_cp = round(np.mean(cp_median_CNVR),2)
				sd_cp = round(np.std(cp_median_CNVR),2)
				fw.write(f'\t{average_cp}\t{sd_cp}\n')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='CNVR detected by scanning the population RD file with aberrantly RD, CNV allele frequency and significantly correlation with adjacent windows.')
	parser.add_argument('-ref', type=str, help='ReferenceDB file.', required=True)
	parser.add_argument('-n', type=str, help='normalized file list.', required=True)
	parser.add_argument('-e', type=str, help='excluded file list', required=True)
	parser.add_argument('-o', type=str, help='output primary CNVR file', required=True, default='PrimaryCNVR')
	parser.add_argument('-w', type=int, help='window size (bp) (default: 1000)', default=1000)
	parser.add_argument('-m', type=int, help='minimum number of homozygous gain/loss individuals when define a candidate CNV window (default: 3)', default=3)
	parser.add_argument('-f', type=float, help='minimum frequency of gain/loss individuals when define a candidate CNV window (default: 0.1)', default=0.1)
	parser.add_argument('-u', type=float, help='heterozygous gain/loss cutoff [ default 1+0.35 for heterozygous gain and 1-0.35 for heterozygous deletion ] (default: 0.35)', default=0.35)
	parser.add_argument('-v', type=float, help='homozygous gain/loss cutoff [ default 1+0.75 for homozygous gain and 1-0.75 for homozygous deletion ] (default: 0.75)', default=0.75)
	parser.add_argument('-r', type=float, help='minimum of pearson correlation coefficient between the two adjacent non-overlapping windows during CNVR discovery \
												0.5  for sample size (0, 30] \
				   								0.4  for sample size (30, 50] \
												0.3  for sample size (50, 100] \
												0.2  for sample size (100, 200] \
												0.15 for sample size (200, 500] \
												0.1  for sample size (500,+∞)')
	args = parser.parse_args()											
	windowsize, minFre, minNum, v, u, r = args.w, args.f, args.m, args.v, args.u ,args.r
	refDBFile, normalizedFile, excludeFile, outfile = args.ref, args.n, args.e, args.o
	r = getr(r)
	filter_sample = getFileterSample(excludelist = excludeFile)
	median75 = getMedian75(normalizedFile=normalizedFile, filter_sample=filter_sample, preHetCutoff=u)
	hash_res_gain, hash_res_loss, hash_res_homo_gain, hash_res_homo_loss, genome_copy, effective_sample = getInfo(normalizedFile=normalizedFile,
																												  filter_sample=filter_sample,
																												  prehetcutoff=u, prehomocutoff=v)
	refdb, gain_effective, loss_effective, hash_res_gain, hash_res_loss = effective(refdbFile=refDBFile, f=minFre, m=minNum, 
																					hash_res_gain=hash_res_gain,
																					hash_res_loss=hash_res_loss,
																					hash_res_homo_gain=hash_res_homo_gain,
																					hash_res_homo_loss=hash_res_homo_loss,
																					effective_sample=effective_sample)
	cnv_region_gain, cnv_region_loss = getCNVregion(refdb=refdb, has_res_gain=hash_res_gain, has_res_loss=hash_res_loss)
	hash_num = getNum(normalizeFile=normalizedFile, CNV_region_gain=cnv_region_gain, cnv_region_loss=cnv_region_loss)
	cnv_tmp_gain, cnv_tmp_loss = getTmpCNV(refdb=refdb, CNV_region_gain=cnv_region_gain, CNV_region_loss=cnv_region_loss,
														gain_effective=gain_effective, loss_effective=loss_effective)
	cnv_gain = getGain(cnv_tmp_gain=cnv_tmp_gain, gain_effective=gain_effective, hash_num=hash_num, r=r)
	cnv_loss = getLoss(cnv_tmp_loss=cnv_tmp_loss, loss_effective=loss_effective, hash_num=hash_num, r=r)
	cnv = getCNV(refdb=refdb, cnv_gain=cnv_gain, cnv_loss=cnv_loss, gain_effective=gain_effective, loss_effective=loss_effective)
	cnv_GC, cnv_repeat, cnv_gap, cnv_real_start, cnv_real_end, kmer_median = info(cnv=cnv, gain_effective=gain_effective, loss_effective=loss_effective,
																				  genome_copy=genome_copy, windowsize=windowsize)
	main(normalizedFile=normalizedFile, outfile=outfile, cnv=cnv, gain_effective=gain_effective, loss_effective=loss_effective,
										cnv_repeat=cnv_repeat, cnv_GC=cnv_GC, cnv_gap=cnv_gap, cnv_real_start=cnv_real_start, cnv_real_end=cnv_real_end,
										kmer_median=kmer_median, windowsize=windowsize)
	