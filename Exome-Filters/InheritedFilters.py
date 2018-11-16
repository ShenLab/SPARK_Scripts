#!/home/local/users/jw/anaconda2/bin/python
import time
import argparse
import yaml
import csv
import re
import os
import pprint

YML_FILTER_FILE = '/home/local/users/jw/CUMC/Exome-Filters-Jiayao/ALL_FILTER.yml'
VQSR = re.compile('([\d.]+)')
pp = pprint.PrettyPrinter(indent=4)


class YML_Filter():
	def __init__(self):
		with open(YML_FILTER_FILE, 'rb') as ymlfile:
			self.cfg = yaml.load(ymlfile)
		self.INFO = self.cfg['INFO']
		self.READS = self.cfg['READS']
		self.Filter = self.cfg['FILTER']
		self.SNP = self.cfg['SNP']
		self.INDEL = self.cfg['INDEL']

	def show(self):
		pp.pprint(self.INFO)
		pp.pprint(self.READS)
		pp.pprint(self.FILTER)
		pp.pprint(self.SNP)
		pp.pprint(self.INDEL)


class Sample():
	# GT:AD:DP:GQ:PL    0/0:7,0:7:18:0,18,270
	def __init__(self, Format, tmp):
		self.Format = tmp
		self.formats = Format.split(':')
		self.info = tmp.split(':')
		if self.info[0] == './.':
			self.GT = [0,0]
			return
		else:
			self.GT = map(int, re.findall('[\d.]', self.info[0]))
		self.fmt = {}
		for i, item in enumerate(self.formats):
			self.fmt[item] = self.info[i]
			if 'AD' in self.fmt:
				tmp = self.fmt['AD'].split(',')
				self.AD = [0] * 2
				self.AD[0] = float(tmp[self.GT[0]])
				self.AD[1] = float(tmp[self.GT[1]])

	def show(self):
		return '{}: {}:{}:{}\t'.format(self.name, self.fmt['GT'], self.fmt['AD'], self.fmt['GQ'])


class Variant():
	def __init__(self, record, headers):
		record = record.strip().split('\t')
		self.headers = headers
		self.List = record
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = record[0:9]
		self.Chrom = CHROM
		self.Pos = POS
		self.Id = ID
		self.Ref = REF
		self.Alt = ALT
		self.Alts = ALT.split(',')
		self.Alleles = [self.Ref] + self.Alts
		self.Qual = QUAL
		self.Filter = FILTER
		self.Info_str = INFO
		self.GetInfo(INFO)
		self.Format = FORMAT

	def GetInfo(self, INFO):
		res = {}
		tmp = INFO.split(';')
		for kv in tmp:
			try:
				k, v = kv.split('=')
				if k not in res:
					res[k] = v.split(',')
				else:
					res[k].extend(v.split(','))
			except:
				continue
		self.Info = res

	def CheckGT(self, idxAllele, Samples, Filters):
		AVGAD = [0, 0]
		AD1, AD2 = 0,0
		AVGPL = 0
		Total = 0
		for sample in Samples:
			if sample.info[0] == './.':
				continue
			if sample.AD[0] != 0:
				AVGAD[0] += float(sample.AD[0])
				AD1 += 1
			if sample.AD[1] != 0:
				AVGAD[1] += float(sample.AD[1])
				AD2 += 1
			AVGPL += float(sample.fmt['GQ'])
			Total += 1
		try:
			AVGAD[0] = AVGAD[0] / AD1 
		except:
			AVGAD[0] = 10
		try:
			AVGAD[1] = AVGAD[1] / AD2 
		except:
			AVGAD[1] = 10
		AVGPL = AVGPL / Total
		if Filters.READS['min_proband_AD'] != None:
			if (AVGAD[0] < Filters.READS['min_proband_AD'] and AVGAD[1] < Filters.READS['min_proband_AD']):
				return False
		if Filters.READS['min_proband_PL'] != None:
			if AVGPL < Filters.READS['min_proband_PL']:
				return False
		return True

	def CheckInfo(self, idxAllele, Filters):
		idx = idxAllele
		if Filters.INFO['exon'] != None and Filters.INFO['exon'] == True:
			if Filters.INFO['exon_flag'] != None:
				VarFunc = self.Info['Func.refGene'][idx]
				# VarFunc = self.Info.get('Func.refGene',[None]*len(self.Alts))[idx]
				# print VarFunc
				if VarFunc not in Filters.INFO['exon_flag']:
					return False
		#	print VarFunc, Proband.show(), Father.show(), Mother.show()
		# if Filters.INFO['max_AC'] != None:
		#	print self.Info['AC'][idx]
		#	if int(self.Info['AC'][idx]) > Filters.INFO['max_AC']:
		#		return False
		# print VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(),
		# Mother.show()
		if Filters.INFO['max_ExAC'] != None:
			if AF(self.Info['ExAC_ALL'][idx]) > Filters.INFO['max_ExAC']:
				return False
		if Filters.INFO['max_gnomAD'] != None:
			if AF(self.Info['gnomAD_genome_ALL'][idx]) > Filters.INFO['max_gnomAD']:
				return False
		if Filters.INFO['max_1KG'] != None:
			if AF(self.Info['1000g2015aug_all'][idx]) > Filters.INFO['max_1KG']:
				return False
		#print AF(self.Info['ExAC_ALL'][idx]), VarFunc,
		   # self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show()
		if Filters.INFO['excluded_gene'] != None:
			for gene in Filters.INFO['excluded_gene']:
				if gene in self.Info['Gene.refGene'][idx]:
					return False
		if Filters.INFO['excluded_chrom'] != None:
			if self.Chrom in Filters.INFO['excluded_chrom']:
				return False
		# print self.Info['Gene.refGene'][idx], AF(self.Info['ExAC_ALL'][idx]),
		# VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(),
		# Mother.show()
		try:
			if Filters.INFO['max_seqdup'] != None:
				segdupScore = self.Info['genomicSuperDups'][idx]
				if segdupScore != '.':
					# print self.Info['genomicSuperDups'][idx]
					segdupScore = re.search(
							'Score:(\d+?\.?\d+)', segdupScore).group(1)
					# print segdupScore
					if float(segdupScore) >= Filters.INFO['max_seqdup']:
						return False
		except:
			pass
		# print self.Info['Gene.refGene'][idx], AF(self.Info['ExAC_ALL'][idx]),
		# VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(),
		# Mother.show(), segdupScore
		return True

	def isSNP(self, idxAllele):
		if len(self.Ref) == 1 and len(self.Alleles[idxAllele + 1]) == 1:
			return True
		else:
			return False

	def isINDEL(self, idxAllele):
		return not self.isSNP(idxAllele)

	def CheckSNP(self, idxAllele, Filters):
		idx = idxAllele
		if Filters.SNP['max_FS'] != None:
			# print self.Info['FS']
			if float(self.Info['FS'][0]) > Filters.SNP['max_FS']:
				return False
		if Filters.SNP['min_QD'] != None:
			if float(self.Info['QD'][0]) < Filters.SNP['min_QD']:
				return False
		if Filters.SNP.get('min_ReadPosRankSum', None) != None:
			if float(self.Info['ReadPosRankSum']) > Filters.SNP['min_ReadPosRankSum']:
				return False
		return True

	def CheckINDEL(self, idxAllele, Filters):
		idx = idxAllele
		# print self.Alleles, self.Info['FS'], self.Info['QD'],
		# self.Info['ReadPosRankSum']
		if Filters.INDEL['max_FS'] != None:
			if float(self.Info['FS'][0]) > Filters.INDEL['max_FS']:
				return False
		if Filters.INDEL['min_QD'] != None:
			if float(self.Info['QD'][0]) < Filters.INDEL['min_QD']:
				return False
		if Filters.INDEL['min_ReadPosRankSum'] != None:
			try:
				if float(self.Info['ReadPosRankSum'][0]) > Filters.INDEL['min_ReadPosRankSum']:
					return False
			except KeyError:
				pass
		return True

	def CheckFilter(self, Filters):
		if self.Filter == 'PASS' or self.Filter == '.':
			return True
		v1, v2 = VQSR.findall(self.Filter)
		if float(v2) > Filters.Filter['VQSRSNP']:
			return False
		elif float(v2) > Filters.Filter['VQSRINDEL']:
			return False
		else:
			return True

	def FindAllele(self):
		Samples = []
		AllelePool = []
		for item in self.List[9:]:
			indi = Sample(self.Format, item)
			Samples.append(indi)
			AllelePool.extend(indi.GT)
		AllelePool = set(AllelePool)
		if len(AllelePool) != 2:
			return None, Samples
		return sorted(list(AllelePool))[1], Samples

	def CheckInherited(self, Filters):
		if not self.CheckFilter(Filters):
			return False
		else:
			pass
		idxAllele, Samples = self.FindAllele()
		if idxAllele == None:
			return False
		idxAllele -= 1
		if not self.CheckGT(idxAllele, Samples, Filters):
			return False
		else:
			# print 'pass GT'
			pass
		if not self.CheckInfo(idxAllele, Filters):
			return False
		else:
			# print 'pass info'
			pass
		if not ((self.isINDEL(idxAllele) and self.CheckINDEL(idxAllele, Filters)) or (self.isSNP(idxAllele) and self.CheckSNP(idxAllele, Filters))):
			return False
		else:
			#print 'pass snp/indel'
			pass
		return True


def AF(CHR):
	if CHR == '.':
		return 0
	else:
		return float(CHR)

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', '--vcf', type=str, help='VCF record')
	parser.add_argument('-o', '--out', type=str, help='Output Prifix')
	args = parser.parse_args()
	return args.vcf, args.out

def FilterInherited(VCF, OUT):
	Filters = YML_Filter()	
	fin = open(VCF, 'rb')
	fout = open(OUT, 'wb')
	for l in fin:
		if l.startswith('##'):
			fout.write(l)
		elif l.startswith('#'):
			fout.write(l)
			Header = l.strip().split('\t')
		else:
			var = Variant(l, Header)
			if var.CheckInherited(Filters):
				fout.write(l)
			elif len(var.Ref) != len(var.Alt):
				print l.strip()
	fin.close()
	fout.close()

def main():
	VCF, OUT = GetOptions()
	FilterInherited(VCF, OUT)
	print 'Inherited Filtering Done'
	return


if __name__ == '__main__':
	main()
