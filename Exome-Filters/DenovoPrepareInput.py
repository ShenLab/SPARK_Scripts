#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# PrepareInput.py
# Prepare bamout script and variant table for IGV shot
#=========================================================================

import argparse
import csv
import re

# SampleName = re.compile('')
SampleName = re.compile('([A-Za-z0-9-]+).bam')
EXM_REF = '/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh'
BAMOUT_CMD = '/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.15b.BamOut.sh'
RUN = 'nohup'


def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', '--csv', type=str, required=True,
			help='csv file contains variants to be IGV. Must have column "Sample", "Chrom", "Pos". ')
	parser.add_argument('-p', '--ped', type=str, required=True,
			help='Pedigree file contains family relationship of samples in csv file.')
	parser.add_argument('-b', '--bam', type=str, required=True,
			help='File contains bam locations. ')
	parser.add_argument('-o', '--output', type=str,
			help='Output 1 of this script, Input variant list for IGV. ')
	parser.add_argument('-s', '--script', type=str,
			help='Output 2 of this script, Script to run GATK Bamout. ')
	args = parser.parse_args()
	if args.output == None:
		args.output = args.csv.rstrip('.csv') + '.varlist'
	if args.script == None:
		args.script = 'run_bamout.sh'
	return args.csv, args.ped, args.bam, args.output, args.script


class Proband:
	def __init__(self, Sample, Father, Mother):
		self.Sample = Sample
		self.Father = Father
		self.Mother = Mother


class Pedigree:
	def __init__(self, pedfile):
		fin = open(pedfile, 'rb')
		self.Probands = {}
		for l in fin:
			if l.startswith('#'):
				continue
			llist = l.strip().split('\t')
			fam, sample, father, mother, sex, phenotype = llist[:6]
			if father != "0" and mother != "0" and phenotype == "2":  # Proband
				# print sample, father, mother
				self.Probands[sample] = Proband(sample, father, mother)


class BAM:
	def __init__(self, fpath):
		self.FullPath = fpath.strip()
		self.BamName = self.FullPath.split('/')[-1]
		self.Sample = SampleName.search(self.BamName).group(1)
		print self.BamName, self.Sample


class BamLocation:
	def __init__(self, bamlocationfile):
		fin = open(bamlocationfile, 'rb')
		self.Bams = {}
		for l in fin:
			if 'bamout' in l:
				continue
			bam = BAM(l)
			self.Bams[bam.BamName] = bam


class Variant:
	def __init__(self, header, row):
		self.Chrom = row[header.index('Chrom')]
		self.Pos = int(row[header.index('Pos')])
		self.Sample = row[header.index('Sample')]
		self.start = max(0, self.Pos - 150)
		self.end = self.Pos + 150

	def addPed(self, pedigree):
		self.Father = pedigree.Probands[self.Sample].Father
		self.Mother = pedigree.Probands[self.Sample].Mother

	def addBam(self, bamlocation):
		print self.Sample
		for bam in bamlocation.Bams.values():
			if self.Sample == bam.Sample:
				self.SampleBam = bam
				self.SampleBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
			if self.Father == bam.Sample:
				self.FatherBam = bam
				self.FatherBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
			if self.Mother == bam.Sample:
				self.MotherBam = bam
				self.MotherBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
		try:
			for var in [self.SampleBam, self.SampleBamout, self.FatherBam, self.FatherBamout, self.MotherBam, self.MotherBamout ]:
				if var not in locals():
					print 'Cant find {} in locals.'.format(var)
					exit()
		except:
			print self.Sample
	def OutVar(self):
		return '{}\t{}\t{}\n'.format(self.Chrom, self.Pos, ','.join([self.SampleBam.BamName, self.SampleBamout, self.FatherBam.BamName, self.FatherBamout, self.MotherBam.BamName, self.MotherBamout]))
	# return '{}\t{}\t{}\n'.format(self.Chrom, self.Pos, ','.join([self.SampleBam.BamName, self.FatherBam.BamName, self.MotherBam.BamName]))
		# return '{}\t{}\t{}\n'.format(self.Chrom, self.Pos, ','.join([self.SampleBam.BamName, self.SampleBamout, self.FatherBam.BamName, self.FatherBamout, self.MotherBam.BamName, self.MotherBamout]))

class Sample:
	def __init__(self, variant):
		self.Sample = variant.Sample
		self.Father = variant.Father
		self.Mother = variant.Mother
		self.SampleBam = variant.SampleBam
		self.FatherBam = variant.FatherBam
		self.MotherBam = variant.MotherBam
		self.SampleBamout = variant.SampleBamout 
		self.FatherBamout = variant.FatherBamout
		self.MotherBamout = variant.MotherBamout
		# self.variants = []
		chrom, start, end = variant.Chrom, variant.start, variant.end 
		self.Intervals = ['{}:{}-{}'.format(str(chrom), str(start), str(end))]
	def AddVar(self, variant):
		# self.vairants.append(variant)
		chrom, start, end = variant.Chrom, variant.start, variant.end
		self.Intervals.append('{}:{}-{}'.format(str(chrom), str(start), str(end))) 
	def FlushBamout(self):
		# SampleCMD =  ''.format(self.SampleBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
		if RUN == 'qsub':
			SampleCMD =  'qsub $CMD -r {} -i {} -t \"{}\"\n'.format('$REF', self.SampleBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
			FatherCMD =  'qsub $CMD -r {} -i {} -t \"{}\"\n'.format('$REF', self.FatherBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
			MotherCMD =  'qsub $CMD -r {} -i {} -t \"{}\"\n'.format('$REF', self.MotherBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
		else:
			SampleCMD =  'nohup $CMD -r {} -i {} -t \"{}\" &\n'.format('$REF', self.SampleBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
			FatherCMD =  'nohup $CMD -r {} -i {} -t \"{}\" &\n'.format('$REF', self.FatherBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
			MotherCMD =  'nohup $CMD -r {} -i {} -t \"{}\" &\n'.format('$REF', self.MotherBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
		return SampleCMD, FatherCMD, MotherCMD

def LoadCSV(_csv, ped, bam, output):
	pedigree = Pedigree(ped)
	Bams = BamLocation(bam)
	res = []
	fout = open(output,'wb')
	with open(_csv, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		header = reader.next()
		for row in reader:
			var = Variant(header, row)
			var.addPed(pedigree)
			var.addBam(Bams)
			res.append(var)
			fout.write(var.OutVar())
	fout.close()
	return res

def ReduceVarian2Sample(variants):
	res = {}
	for variant in variants:
		if variant.Sample not in res:
			res[variant.Sample] = Sample(variant)
		else:
			res[variant.Sample].AddVar(variant)
	return res.values()

def FlushBamout(variants, script):
	fout = open(script,'wb')
	fout.write('REF={}\n'.format(EXM_REF))
	fout.write('CMD={}\n\n'.format(BAMOUT_CMD)) 
	Samples = ReduceVarian2Sample(variants)
	for sample in Samples:
		cmd1, cmd2, cmd3 = sample.FlushBamout()
		fout.write('{}\n{}\n{}\n\n'.format(cmd1, cmd2, cmd3))
	fout.close()

def main():
	csv, ped, bam, output, script = GetOptions()
	variants = LoadCSV(csv, ped, bam, output)
	FlushBamout(variants, script)
	return


if __name__ == '__main__':
	main()
