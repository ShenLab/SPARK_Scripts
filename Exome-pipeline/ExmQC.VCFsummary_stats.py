#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# summarize vcf files
# Input: VCF file; Output: a table of vcf summary and several pdf visulization
#=========================================================================

from __future__ import division
from optparse import OptionParser
import matplotlib as mpl
mpl.use('Agg')
import gzip
import multiprocessing
import os
import re
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#=========================================================================
# Table data For Summary of the VCF.
# 0th Column is Totol summary by locus
# 1th Column is Total summary weighted by samples
# start from 2th Column is summary of each sample
#=========================================================================


class stats():
	def __init__(self, line, WhichCode):
		self.WhichCode = WhichCode
		self.Samples = line.split("\t")[7:]
		self.Samples[0] = 'Total(locus)'
		self.Samples[1] = 'SampleAVG'
		self.SamLength = len(self.Samples)
		self.SNVcount = [0] * self.SamLength
		self.InDelCount = [0] * self.SamLength
		self.KnownCount = [0] * self.SamLength
		self.NovelCount = [0] * self.SamLength
		self.knownTiCount = [0] * self.SamLength
		self.knownTvCount = [0] * self.SamLength
		self.novelTiCount = [0] * self.SamLength
		self.novelTvCount = [0] * self.SamLength
		self.silentCount = [0] * self.SamLength
		self.missenseCount = [0] * self.SamLength
		self.nonsenseCount = [0] * self.SamLength
		self.unknownCount = [0] * self.SamLength
		self.frameShiftCount = [0] * self.SamLength
		self.homozygousCount = [0] * self.SamLength
		self.heterozygousCount = [0] * self.SamLength
		self.referenceCount = [0] * self.SamLength
		self.KGRareCount = [0] * self.SamLength
		self.gnomADRareCount = [0] * self.SamLength
		self.RareCount = [0] * self.SamLength
		self.ExACRareCount = [0] * self.SamLength
		self.NoAAF = [0] * self.SamLength
		self.KnownTiTvRat = [0] * self.SamLength
		self.NovelTiTvRat = [0] * self.SamLength
		self.TotalTiTvRat = [0] * self.SamLength
		self.UnCalled = [0] * self.SamLength


def GetOptions():
	parser = OptionParser()
	parser.add_option("-v", "--vcf", dest="VCFfile",
			help="Input VCF file", metavar="VCFfile")
	parser.add_option("-q", "--qc_folder", dest="QC_folder",
			help="User specified name of folder for output to be written", metavar="QC_folder", default="QC_summary")
	parser.add_option("-p", "--num_procs", dest="Nprocs",
			help="Number of process running", metavar="Nprocs", default=1, type=int)
	(options, args) = parser.parse_args()
	QC_folder = os.path.abspath(options.QC_folder)
	if not os.path.exists(QC_folder):
		os.mkdir(QC_folder)
	return options.VCFfile, QC_folder, options.Nprocs


def GetInfoDict(linelist):
	INFOstring = linelist[7]
	INFOcolumnList = INFOstring.split(';')
	INFOdict = {}
	for item in INFOcolumnList:
		if '=' in item:
			Key, Value = item.split('=', 1)
			if Key not in INFOdict:
				INFOdict[Key] = Value
			else:
				INFOdict[Key] = INFOdict[Key] + ',' + Value
	return INFOdict


def InitStats(line):
	return stats(line, "Coding Variants"), stats(line, "Non-coding Variants"), stats(line, "All Variants")


def changevalue(value):
	if value != '.':
		return float(value)
	else:
		return 2


def GetFrequency(INFOdict):
	KGseq = [changevalue(rate) for rate in INFOdict.get(
		'1000g2015aug_all', 2).split(',')]
	gnomADseq = [changevalue(rate) for rate in INFOdict.get(
		'gnomAD_genome_ALL', 2).split(',')]
	ExACseq = [changevalue(rate)
			for rate in INFOdict.get('ExAC_ALL', 2).split(',')]
	KGseq.append(0)
	gnomADseq.append(0)
	ExACseq.append(0)
	KGscore = max(KGseq)
	gnomADscore = max(gnomADseq)
	ExACscore = max(ExACseq)
	return KGscore, gnomADscore, ExACscore

# Only get GT


def ParseFormat(Format, Sample):
	Dict = {}
	Formats = Format.split(":")
	Samples = Sample.split(":")
	for i in xrange(len(Formats)):
		try:
			Dict[Formats[i]] = Samples[i]
			if "GT" in Dict:
				Dict['GT'] = re.findall('[\d.]', Dict['GT'])
				return Dict
		except IndexError:
			# print Format,Sample
			Dict[Formats[i]] = '.'
	Dict['GT'] = re.findall('[\d.]', Dict['GT'])
	return Dict


def WriteSummary(STATS, outfname):
	Output = open(outfname, 'wb')
	for tmp in STATS:
		for i in range(0, tmp.SamLength):
			if int(tmp.knownTvCount[i]) > 0:
				tmp.KnownTiTvRat[i] = float(
						tmp.knownTiCount[i]) / float(tmp.knownTvCount[i])
				if int(tmp.novelTvCount[i]) > 0:
					tmp.NovelTiTvRat[i] = float(
							tmp.novelTiCount[i]) / float(tmp.novelTvCount[i])
					tmp.TotalTvCount = float(
							tmp.novelTvCount[i]) + float(tmp.knownTvCount[i])
					tmp.TotalTiCount = float(
							tmp.novelTiCount[i]) + float(tmp.knownTiCount[i])
					if tmp.TotalTvCount > 0:
						tmp.TotalTiTvRat[i] = float(
								tmp.TotalTiCount) / float(tmp.TotalTvCount)
						# write output
		Output.write("\t".join([str(tmp.WhichCode)]) + "\n")
		Output.write("\t".join(['Samples:'] + [str(i) for i in tmp.Samples]))
		Output.write("\t".join(['SNVs'] + [str(i)
			for i in tmp.SNVcount]) + "\n")
		Output.write("\t".join(['InDels'] + [str(i)
			for i in tmp.InDelCount]) + "\n")
		Output.write("\t".join(['Known'] + [str(i)
			for i in tmp.KnownCount]) + "\n")
		Output.write("\t".join(['Novel'] + [str(i)
			for i in tmp.NovelCount]) + "\n")
		Output.write("\t".join(['known Ti'] + [str(i)
			for i in tmp.knownTiCount]) + "\n")
		Output.write("\t".join(['known TV'] + [str(i)
			for i in tmp.knownTvCount]) + "\n")
		Output.write("\t".join(['known Ti/TV ratio'] + [str(i)
			for i in tmp.KnownTiTvRat]) + "\n")
		Output.write("\t".join(['novel Ti'] + [str(i)
			for i in tmp.novelTiCount]) + "\n")
		Output.write("\t".join(['novel TV'] + [str(i)
			for i in tmp.novelTvCount]) + "\n")
		Output.write("\t".join(['novel Ti/TV ratio'] + [str(i)
			for i in tmp.NovelTiTvRat]) + "\n")
		Output.write("\t".join(['total Ti/TV ratio'] + [str(i)
			for i in tmp.TotalTiTvRat]) + "\n")
		if str(tmp.WhichCode) != 'Non-coding Variants':
			Output.write("\t".join(['Silent'] + [str(i)
				for i in tmp.silentCount]) + "\n")
			Output.write("\t".join(['Missense'] + [str(i)
				for i in tmp.missenseCount]) + "\n")
			Output.write("\t".join(['Nonsense'] + [str(i)
				for i in tmp.nonsenseCount]) + "\n")
			Output.write("\t".join(['Unknown'] + [str(i)
				for i in tmp.unknownCount]) + "\n")
			Output.write("\t".join(['Frameshift'] + [str(i)
				for i in tmp.frameShiftCount]) + "\n")
			Output.write("\t".join(['Homozygous'] + [str(i)
				for i in tmp.homozygousCount]) + "\n")
			Output.write("\t".join(['Heterozygous'] + [str(i)
				for i in tmp.heterozygousCount]) + "\n")
			Output.write("\t".join(['Rare (AAF<0.01) - 1KG'] +
				[str(i) for i in tmp.KGRareCount]) + "\n")
			Output.write("\t".join(['Rare (AAF<0.01) - ExAC'] +
				[str(i) for i in tmp.ExACRareCount]) + "\n")
			Output.write("\t".join(['Rare (AAF<0.01) - gnomAD'] + [str(i)
				for i in tmp.gnomADRareCount]) + "\n")
			Output.write("\t".join(
				['Rare (AAF<0.01) - 1KG or gnomAD'] + [str(i) for i in tmp.RareCount]) + "\n")
			Output.write("\t".join(['No AAF'] + [str(i)
				for i in tmp.NoAAF]) + "\n")
			Output.write("\t".join(['not called'] + [str(i)
				for i in tmp.UnCalled]) + "\n")
			Output.write("\t\n")
	Output.close()


def PlotSummary(QC_folder, outname):
	print "Plotting"
	plot_caterogy = ['SNVs', 'InDels', 'known Ti/TV ratio',
			'novel Ti/TV ratio', 'Rare (AAF<0.01) - ExAC']
	with open(outname) as f:
		for line in f:
			line = line.strip()

			if line in ["Coding Variants", "Non-coding Variants", "All Variants"]:
				output_type = line
				continue
			if line == '':
				continue
			if line.startswith('Samples'):
				samples = line.split()
				print '# of samples in vcf:', len(samples) - 2
				continue

			info = line.split('\t')
			caterogy = info[0]
			if caterogy in plot_caterogy:
				print 'plot: {} - {}'.format(output_type, caterogy)
				hist = map(float, info[3:])
				print max(hist), hist.index(max(hist))
				pdfname = os.path.join(
						QC_folder, output_type + '-' + caterogy.replace('/', '') + '.pdf')
				pdf = PdfPages(pdfname)
				fig, ax = plt.subplots(dpi=100)
				rects1 = ax.hist(hist, bins=50)
				ax.spines['right'].set_visible(False)
				ax.spines['top'].set_visible(False)
				ax.xaxis.set_ticks_position('bottom')
				ax.yaxis.set_ticks_position('left')
				ax.set_title(output_type + ' ' + caterogy + '\n' + 'mean:' + str(np.round(np.mean(hist), decimals=3)) + ' ' +
						'median:' + str(np.round(np.median(hist), decimals=3)))
				plt.show()
				plt.tight_layout()
				pdf.savefig(bbox_inches='tight')
				pdf.close()
				plt.close()


def Summary(VCF, QCdir, Nprocs):
	Nucleotides = ['A', 'C', 'G', 'T']
	MutUnknown = ['unknown']
	MutSilent = ['synonymous_SNV']
	MutMissense = ['nonsynonymous_SNV']
	MutNonsense = ['stopgain', 'stoploss']
	MutNonframeshift = ['nonframeshift_insertion', 'nonframeshift_deletion']
	MutFrameshift = ['frameshift_insertion', 'frameshift_deletion']
	CodingCodes = ['splicing', 'exonic', 'exonic-splicing']

	print '=' * 50
	print 'QC summary start'

	if VCF.endswith('vcf.gz'):
		fin = gzip.open(VCF, 'rb')
	else:
		fin = open(VCF, 'rb')

	Count = 0
	for line in fin:
		if line.startswith('#CHROM'):
			coding, noncoding, All = InitStats(line)
		elif not line.startswith("#"):
			# Get Infos
			Count += 1
			linelist = line.split("\t")
			INFOdict = GetInfoDict(linelist)
			KGscore, gnomADscore, ExACscore = GetFrequency(INFOdict)

			# print KGscore,gnomADscore,ExACscore

			MutationFunct = str(INFOdict.get(
				'Func.refGene', 'none').split(',')[0])
			MutationClass = str(INFOdict.get(
				'ExonicFunc.refGene', 'none').split(',')[0])
			ID = str(linelist[2])
			REF = [str(item) for item in linelist[3].split(',')]
			ALT = [str(item) for item in linelist[4].split(',')]
			# Two Allele
			#if len(ALT) > 1:
			#	continue
			Known = (INFOdict['avsnp147'].split(',')[0] != '.' or ID != '.')
			Indel = not (all(i in Nucleotides for i in REF)
					and all(i in Nucleotides for i in ALT))
			Transition = (REF[0] == 'A' and ALT[0] == 'G') or (REF[0] == 'G' and ALT[0] == 'A') or (
					REF[0] == 'C' and ALT[0] == 'T') or (REF[0] == 'T' and ALT[0] == 'C')
			Transversion = ((REF[0] == 'C' and ALT[0] == 'A') or (REF[0] == 'A' and ALT[0] == 'C') or (REF[0] == 'G' and ALT[0] == 'C') or (REF[0] == 'C' and ALT[0] == 'G') or (
				REF[0] == 'T' and ALT[0] == 'A') or (REF[0] == 'A' and ALT[0] == 'T') or (REF[0] == 'T' and ALT[0] == 'G') or (REF[0] == 'G' and ALT[0] == 'T'))

			# Update stats
			for i in range(0, All.SamLength):
				ColNum = i + 7
				if i > 1:
					if './.' not in linelist[ColNum]:
						INFOSPL = ParseFormat(linelist[8], linelist[ColNum])
						if INFOSPL['GT'] != ['0', '0']:
							if '0' in INFOSPL['GT']:
								All.heterozygousCount[i] += 1
								if (MutationFunct in CodingCodes):
									coding.heterozygousCount[i] += 1
								else:
									noncoding.heterozygousCount[i] += 1
							else:
								All.homozygousCount[i] += 1
								if (MutationFunct in CodingCodes):
									coding.homozygousCount[i] += 1
								else:
									noncoding.homozygousCount[i] += 1
						GT = '/'.join(INFOSPL['GT'])
					else:
						All.UnCalled[i] += 1
						if (MutationFunct in CodingCodes):
							coding.UnCalled[i] += 1
						else:
							noncoding.UnCalled[i] += 1
						GT = './.'
				else:
					GT = '0/0'
				if (GT != '0/0' and GT != './.') or i == 0:
					if Indel:
						All.InDelCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.InDelCount[i] += 1
						else:
							noncoding.InDelCount[i] += 1
					else:
						All.SNVcount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.SNVcount[i] += 1
						else:
							noncoding.SNVcount[i] += 1
					if Known:
						All.KnownCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.KnownCount[i] += 1
						else:
							noncoding.KnownCount[i] += 1
						if Transition:
							All.knownTiCount[i] += 1
							if (MutationFunct in CodingCodes):
								coding.knownTiCount[i] += 1
							else:
								noncoding.knownTiCount[i] += 1
						if Transversion:
							All.knownTvCount[i] += 1
							if (MutationFunct in CodingCodes):
								coding.knownTvCount[i] += 1
							else:
								noncoding.knownTvCount[i] += 1
					else:
						All.NovelCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.NovelCount[i] += 1
						else:
							noncoding.NovelCount[i] += 1
						if Transition:
							All.novelTiCount[i] += 1
							if (MutationFunct in CodingCodes):
								coding.novelTiCount[i] += 1
							else:
								noncoding.novelTiCount[i] += 1
						if Transversion:
							All.novelTvCount[i] += 1
							if (MutationFunct in CodingCodes):
								coding.novelTvCount[i] += 1
							else:
								noncoding.novelTvCount[i] += 1
					if MutationClass in MutSilent:
						All.silentCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.silentCount[i] += 1
						else:
							noncoding.silentCount[i] += 1
					elif MutationClass in MutMissense:
						All.missenseCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.missenseCount[i] += 1
						else:
							noncoding.missenseCount[i] += 1
					elif MutationClass in MutNonsense:
						All.nonsenseCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.nonsenseCount[i] += 1
						else:
							noncoding.nonsenseCount[i] += 1
					elif MutationClass in MutUnknown:
						All.unknownCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.unknownCount[i] += 1
						else:
							noncoding.unknownCount[i] += 1
					elif MutationClass in MutFrameshift:
						All.frameShiftCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.frameShiftCount[i] += 1
						else:
							noncoding.frameShiftCount[i] += 1

					if KGscore <= 0.01:
						All.KGRareCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.KGRareCount[i] += 1
						else:
							noncoding.KGRareCount[i] += 1
					if gnomADscore <= 0.01:
						All.gnomADRareCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.gnomADRareCount[i] += 1
						else:
							noncoding.gnomADRareCount[i] += 1
					if ExACscore <= 0.01:
						All.ExACRareCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.ExACRareCount[i] += 1
						else:
							noncoding.ExACRareCount[i] += 1

					if KGscore <= 0.01 or gnomADscore <= 0.01:
						All.RareCount[i] += 1
						if (MutationFunct in CodingCodes):
							coding.RareCount[i] += 1
						else:
							noncoding.RareCount[i] += 1

					if ExACscore == 2 and KGscore == 2 and gnomADscore == 2:
						All.NoAAF[i] = All.NoAAF[i] + 1
						if (MutationFunct in CodingCodes):
							coding.NoAAF[i] += 1
						else:
							noncoding.NoAAF[i] += 1

			if Count > 0 and Count % 1000 == 0:
				print "Read %d Variants..." % Count
	print "\nRead All %d Variants" % Count
	# Loop over all variants Finished, update total summary weight by sample
	Len = len(All.Samples) - 2
	for i in xrange(2, len(All.Samples)):
		for group in [coding, noncoding, All]:
			group.SNVcount[1] = sum(group.SNVcount[2:])/Len
			group.InDelCount[1] = sum(group.InDelCount[2:])/Len
			group.KnownCount[1] = sum(group.KnownCount[2:])/Len
			group.NovelCount[1] = sum(group.NovelCount[2:])/Len
			group.knownTiCount[1] = sum(group.knownTiCount[2:])/Len
			group.knownTvCount[1] = sum(group.knownTvCount[2:])/Len
			group.novelTiCount[1] = sum(group.novelTiCount[2:])/Len
			group.novelTvCount[1] = sum(group.novelTvCount[2:])/Len
			group.silentCount[1] = sum(group.silentCount[2:])/Len
			group.missenseCount[1] = sum(group.missenseCount[2:])/Len
			group.nonsenseCount[1] = sum(group.nonsenseCount[2:])/Len
			group.unknownCount[1] = sum(group.unknownCount[2:])/Len
			group.frameShiftCount[1] = sum(group.frameShiftCount[2:])/Len
			group.homozygousCount[1] = sum(group.homozygousCount[2:])/Len
			group.heterozygousCount[1] = sum(group.heterozygousCount[2:])/Len
			group.referenceCount[1] = sum(group.referenceCount[2:])/Len
			group.KGRareCount[1] = sum(group.KGRareCount[2:])/Len
			group.gnomADRareCount[1] = sum(group.gnomADRareCount[2:])/Len
			group.RareCount[1] = sum(group.RareCount[2:])/Len
			group.ExACRareCount[1] = sum(group.ExACRareCount[2:])/Len
			group.NoAAF[1] = sum(group.NoAAF[2:])/Len
			group.UnCalled[1] = sum(group.UnCalled[2:])/Len

	STATS = [coding, noncoding, All]
	outfname = os.path.join(QCdir, 'sample.summary.txt')
	WriteSummary(STATS, outfname)

	print "QC Summary Finished"
	print '-' * 50
	return QCdir, outfname


def main():
	VCF, QCdir, Nprocs = GetOptions()
	QCdir, outfname = Summary(VCF, QCdir, Nprocs)
	PlotSummary(QCdir, os.path.join(QCdir, 'sample.summary.txt'))


if __name__ == '__main__':
	main()
