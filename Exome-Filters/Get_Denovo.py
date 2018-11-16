#!/home/yufengshen/anaconda2/bin/python
# Given set of vcf files, get the avg of N denovo and N inherited variants
import argparse
import re
import os
import utils
import time
import yaml
import csv
import pprint

CSV_HEADER = ['Sample', 'Phenotype', 'Chrom', 'Pos', 'Ref', 'Alt', 'AC', 'Mappability', 'Gene', 'GeneName', 'GeneFunc', 'ExonicFunc', 'AAchange', 'ExAC', 'gnomAD', 'VarType', 'REVEL', 'MPC','MCAP', 'MetaSVM', 'CADD', 'PP2', '1KG', 'mis_z', 'lof_z','pLI', 'pRec','HeartRank', 'LungRank', 'BrainRank','Filter', 'QUAL', 'ProbandGT', 'FatherGT', 'MotherGT', 'Relateness']
VQSR = re.compile('([\d.]+)')
pp = pprint.PrettyPrinter(indent=4)

GENENAME = '/home/local/users/jw/resources/Annotations/protein-coding_gene.txt'
EXACSCORE = '/home/local/users/jw/resources/Annotations/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'
HEARTRANK = '/home/local/users/jw/resources/Annotations/diaphragm_rank.csv'
LUNGRANK = '/home/local/users/jw/resources/Annotations/Lung_rank_RNAseq_Asselin-Labat-GPL13112_human.csv'
MOUSEBRIANRANK = '/home/local/users/jw/resources/Annotations/mousebrain.csv'

ALL_FILTER='/home/local/users/jw/CUMC/Exome-Filters-Jiayao/ALL_FILTER.yml'
DENOVO_FILTER='/home/local/users/jw/CUMC/Exome-Filters-Jiayao/DENOVO_FILTER.yml'



class GENE_ANNOTATION:
	def __init__(self, GeneName=True, ExAC=True, Diaphragm=True, Lung=True, MouseBrain=True):
		self.Genes = {}
		if GeneName == True:
			self.Load_GeneName()
		if ExAC == True:
			self.Load_ExAC()
		if Diaphragm == True:
			self.Load_Diaphragm()
		if Lung == True:
			self.Load_Lung()
		if MouseBrain == True:
			self.Load_MouseBrain()

	def Load_GeneName(self):
		stime = time.time()
		fin = open(GENENAME, 'rb')
		Header = fin.readline().strip().split('\t')
		idx_symbol = Header.index('symbol')
		idx_name = Header.index('name')
		for row in fin:
			row = row.strip().split('\t')
			gene, geneName = row[idx_symbol], row[idx_name]
			if gene not in self.Genes:
				self.Genes[gene] = GENE(gene)
			self.Genes[gene].Name = geneName
		print 'Finished Reading Gene Score MouseBrain %.3f' % (time.time() - stime)

	def Load_ExAC(self):
		stime = time.time()
		fin = open(EXACSCORE, 'rb')
		Header = fin.readline().strip().split('\t')
		idx_gene = Header.index('gene')
		idx_Mis_z, idx_Lof_z, idx_pLI, idx_pRec = Header.index(
				'mis_z'), Header.index('lof_z'), Header.index('pLI'), Header.index('pRec')
		for l in fin:
			llist = l.strip().split('\t')
			Gene, Mis_z, Lof_z, pLI, pRec = llist[idx_gene], llist[
					idx_Mis_z], llist[idx_Lof_z], llist[idx_pLI], llist[idx_pRec]
			if Gene not in self.Genes:
				self.Genes[Gene] = GENE(Gene)
			self.Genes[Gene].Mis_z = Mis_z
			self.Genes[Gene].Lof_z = Lof_z
			self.Genes[Gene].pLI = pLI
			self.Genes[Gene].pRec = pRec
		print 'Finished Reading Gene Score ExAC %.3f' % -(stime - time.time())

	def Load_Diaphragm(self):
		stime = time.time()
		Heart = open(HEARTRANK, 'rb')
		Reader = csv.reader(Heart)
		Header = Reader.next()
		idx_heart = Header.index('rank')
		for row in Reader:
			Gene, Rank = row[0], row[idx_heart]
			if Gene not in self.Genes:
				self.Genes[Gene] = GENE(Gene)
			self.Genes[Gene].DiaphragmRank = str(100 - float(Rank))
		print 'Finished Reading Gene Score Diaphragm %.3f' % (time.time() - stime)

	def Load_Lung(self):
		stime = time.time()
		fin = open(LUNGRANK, 'rb')
		Reader = csv.reader(fin)
		Header = Reader.next()
		idx_rank = Header.index('Control-Stroma rank')
		self.Lung = {}
		for row in Reader:
			Gene = row[1]
			Rank = row[idx_rank]
			if Gene not in self.Genes:
				self.Genes[Gene] = GENE(Gene)
			try:
				self.Genes[Gene].LungRank = str(100 - float(Rank))
			except:
				self.Genes[Gene].LungRank = "0"
		print 'Finished Reading Gene Score Lung %.3f' % (time.time() - stime)

	def Load_MouseBrain(self):
		stime = time.time()
		Brain = open(MOUSEBRIANRANK, 'rb')
		Reader = csv.reader(Brain)
		Header = Reader.next()
		idx_Brain = Header.index('brain_rank')
		self.MouseBrain = {}
		for row in Reader:
			Gene, Rank = row[0], row[idx_Brain]
			if Gene not in self.Genes:
				self.Genes[Gene] = GENE(Gene)
			self.Genes[Gene].MouseBrainRank = Rank
		print 'Finished Reading Gene Score MouseBrain %.3f' % (time.time() - stime)

class GENE:
	def __init__(self, symbol):
		self.Symbol = symbol
		self.Name = '.'
		self.pLI = '.'
		self.Mis_z = '.'
		self.Lof_z = '.'
		self.pRec = '.'
		self.DiaphragmRank = '.'
		self.LungRank = '.'
		self.MouseBrainRank = '.'

class YML_Filter():
	def __init__(self, yml_file):
		print yml_file
		yml_dict = yaml.load(open(yml_file,'rb'))
		print yml_dict
		self.INFO = yml_dict['INFO']
		self.READS = yml_dict['READS']
		self.FILTER = yml_dict['FILTER']
		if self.FILTER == None:
			self.FILTER = {}
		self.SNP = yml_dict['SNP']
		self.INDEL = yml_dict['INDEL']

	def show(self):
		pp.pprint(self.INFO)
		pp.pprint(self.READS)
		pp.pprint(self.FILTER)
		pp.pprint(self.SNP)
		pp.pprint(self.INDEL)

class Sample():
	# GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
	def __init__(self, name, Format, tmp):
		self.name = name
		self.Format = tmp
		self.formats = Format.split(':')
		self.info = tmp.split(':')
		self.GT = map(int, re.findall('[\d.]', self.info[0]))
		self.fmt = {}
		for i, item in enumerate(self.formats):
			self.fmt[item] = self.info[i]
		if 'AD' in self.fmt:
			tmp = self.fmt['AD'].split(',')
			self.AD = [0]*2
			self.AD[0] = float(tmp[self.GT[0]])
			self.AD[1] = float(tmp[self.GT[1]])
	def show(self):
		return '{}: {}:{}:{}\t'.format(self.name,self.fmt['GT'],self.fmt['AD'],self.fmt['GQ'])

class Variant():
	def __init__(self, record, headers):
		record = record.strip().split('\t')
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
				k,v = kv.split('=')
				if k not in res:
					res[k] = v.split(',')
				else:
					res[k].extend(v.split(','))
			except:
				continue
		self.Info = res

	def CheckGT(self, Proband, Father, Mother, Filters):
		#print Filters.READS
		#exit()
		Pool = set([Father.GT[0], Father.GT[1], Mother.GT[0], Mother.GT[1]])
		if Proband.GT[0] not in Pool or Proband.GT[1] not in Pool:
			#print Proband.GT, Father.GT, Mother.GT
			if Filters.READS['min_proband_AD'] != None: 
				if (Proband.AD[0] < Filters.READS['min_proband_AD'] or Proband.AD[1] < Filters.READS['min_proband_AD']) :
					return False
				#print Proband.AD
			if Filters.READS['min_proband_PL'] != None :
				if (Proband.fmt['GQ']) < Filters.READS['min_proband_PL']:
					return False
			if Filters.READS['min_proband_alt_freq'] != None :
				if  (Proband.AD[0]/float(Proband.fmt['DP']) < Filters.READS['min_proband_alt_freq'] or Proband.AD[1]/float(Proband.fmt['DP']) < Filters.READS['min_proband_alt_freq']):
					return False
			if Filters.READS['min_parents_DP'] != None :
				try:
					if  (int(Father.fmt['DP']) < Filters.READS['min_parents_DP'] or int(Mother.fmt['DP'])< Filters.READS['min_parents_DP']):
						return False
				except ValueError:
					return False
			if Filters.READS['min_parents_GQ'] != None :
				if  (Proband.fmt['GQ']) < Filters.READS['min_parents_GQ']:
					return False
			if Filters.READS['min_parents_ref_freq'] != None :
				if (float(Father.AD[0])/float(Father.fmt['DP']) < Filters.READS['min_parents_ref_freq'] and Mother.AD[0] / float(Mother.fmt['DP']) < Filters.READS['min_parents_ref_freq'] ):
					return False
			#print Proband.show(), Father.show(), Mother.show()
			return True
		else:
			return False
	def CheckInfo(self, Proband, Father, Mother, Filters):
		idx = Proband.GT[1]-1
		if Filters.INFO['exon'] != None and Filters.INFO['exon'] == True:
			if Filters.INFO['exon_flag'] != None:
				VarFunc = self.Info['Func.refGene'][idx]
				#VarFunc = self.Info.get('Func.refGene',[None]*len(self.Alts))[idx]
				#print VarFunc
				if VarFunc not in Filters.INFO['exon_flag']:
					return False
			#print VarFunc, Proband.show(), Father.show(), Mother.show()
		if Filters.INFO['max_AC'] != None:
			if  int(self.Info['AC'][idx]) > Filters.INFO['max_AC']:
				return False
		#print VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show()
		if Filters.INFO['max_ExAC'] != None:
			if AF(self.Info['ExAC_ALL'][idx]) > Filters.INFO['max_ExAC']:
				return False
		if Filters.INFO['max_gnomAD'] != None:
			if AF(self.Info['gnomAD_genome_ALL'][idx]) > Filters.INFO['max_gnomAD']:
				return False
		if Filters.INFO['max_1KG'] != None:
			if AF(self.Info['1000g2015aug_all'][idx]) > Filters.INFO['max_1KG']:
				return False
	   # print AF(self.Info['ExAC_ALL'][idx]), VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show()
		if Filters.INFO['excluded_gene'] != None:
			for gene in Filters.INFO['excluded_gene']:
				if gene in self.Info['Gene.refGene'][idx]:
					return False
		if Filters.INFO['excluded_chrom'] != None:
			if self.Chrom in Filters.INFO['excluded_chrom']:
				return False
		#print self.Info['Gene.refGene'][idx], AF(self.Info['ExAC_ALL'][idx]), VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show()
		if Filters.INFO.get('max_seqdup',None) != None:
			segdupScore = self.Info['genomicSuperDups'][idx]
			if segdupScore != '.':
				#print self.Info['genomicSuperDups'][idx]
				segdupScore = re.search('Score:(\d+?\.?\d+)',segdupScore).group(1)
				#print segdupScore
				if float(segdupScore) >= Filters.INFO['max_seqdup']:
					return False
		if Filters.INFO.get('min_Mappability',None) != None:
			Mappability = self.Info['Mappability'][0]
			if float(Mappability) < float(Filters.INFO['min_Mappability']):
				return False
		#print self.Info['Gene.refGene'][idx], AF(self.Info['ExAC_ALL'][idx]), VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show(), segdupScore
		return True

	def isSNP(self, Proband):
		if len(self.Alleles[Proband.GT[0]]) == 1 and len(self.Alleles[Proband.GT[1]]) == 1:
			return True
		else:
			return False

	def CheckSNP(self, Proband, Filters):
		idx = Proband.GT[1] - 1
		if Filters.SNP['max_FS'] != None:
			#print self.Info['FS']
			if float(self.Info['FS'][0]) > Filters.SNP['max_FS']:
				return False
		if Filters.SNP['min_QD'] != None:
			if float(self.Info['QD'][0]) < Filters.SNP['min_QD']:
				return False
		if Filters.SNP.get('min_ReadPosRankSum',None) != None:
			if float(self.Info['ReadPosRankSum']) > Filters.SNP['min_ReadPosRankSum']:
				return False
		return True

	def CheckINDEL(self, Proband, Filters):
		idx = Proband.GT[1] - 1
		#print self.Alleles, self.Info['FS'], self.Info['QD'], self.Info['ReadPosRankSum']
		if Filters.INDEL['max_FS'] != None:
			if float(self.Info['FS'][0]) > Filters.INDEL['max_FS']:
				return False
		if Filters.INDEL['min_QD'] != None:
			if float(self.Info['QD'][0]) < Filters.INDEL['min_QD']:
				return False
		if Filters.INDEL['min_ReadPosRankSum'] != None:
			try:
				if float(self.Info['ReadPosRankSum'][0]) < Filters.INDEL['min_ReadPosRankSum']:
					return False
			except KeyError:
				pass
		return True

	def isINDEL(self, Proband):
		return not self.isSNP(Proband)

	def CheckFILTER(self, Proband, Filters):
		if self.Filter == 'PASS' or self.Filter == '.':
			return True
		try:
			v1, v2 = VQSR.findall(self.Filter)
			if Filters.FILTER.get('VQSRSNP', None) != None:
				if float(v2) > Filters.FILTER.get('VQSRSNP', None):
					return False
			elif Filters.FILTER.get('VQSRINDEL', None) != None:
				if float(v2) > Filters.FILTER.get('VQSRINDEL', None):
					return False
			return True
		except Exception, e:
			print str(e)
			print self.Filter

	def CheckDeNovo(self, headers, Pedigree, Filters):
		Proband = self.List[headers.index(Pedigree.Proband.Sample)]
		Father = self.List[headers.index(Pedigree.Father.Sample)]
		Mother = self.List[headers.index(Pedigree.Mother.Sample)]
		try:    
			Proband = Sample(Pedigree.Proband.Sample, self.Format, Proband)
			Father = Sample(Pedigree.Father.Sample, self.Format, Father)
			Mother = Sample(Pedigree.Mother.Sample, self.Format, Mother)
		except ValueError:
			#print self.List[headers.index(Pedigree.Proband.Sample)], self.List[headers.index(Pedigree.Father.Sample)], self.List[headers.index(Pedigree.Mother.Sample)] 
			return False, None, None, None 
		if not self.CheckGT(Proband, Father, Mother, Filters):
			return False, None, None, None
		#if not self.CheckFILTER(Proband, Filters):
		#	return False, None, None, None
		else:
			#print 'pass GT'
			pass 
		if not self.CheckInfo(Proband, Father, Mother, Filters):
			return False, None, None, None
		else:
			#print 'pass info'
			pass
		if not ( (self.isINDEL(Proband) and self.CheckINDEL(Proband, Filters)) or (self.isSNP(Proband) and self.CheckSNP(Proband, Filters)) ):
			return False, None, None, None
		else:
			#print 'pass snp/indel'
			pass
		return True, Proband, Father, Mother

	def show(self):
		print '\t'.join(self.List)

	def OutAsCSV(self, Proband, Father, Mother, Pedigree, genescore):
		idx = Proband.GT[1]-1
		Gene = self.Info['Gene.refGene'][Proband.GT[1]-1]
		try:
			GeneName = genescore[Gene].Name
		except:
			GeneName = '.'
			tmp = GENE(Gene)
			genescore[tmp.Symbol] = tmp
		VarType = self.GetVarType(self.Info['Func.refGene'][idx],self.Info['ExonicFunc.refGene'][idx],self.Info['REVEL'][idx])
		mis_z = genescore[Gene].Mis_z
		lof_z = genescore[Gene].Lof_z
		pLI = genescore[Gene].pLI
		pRec = genescore[Gene].pRec
		LungRank = genescore[Gene].LungRank
		HeartRank = genescore[Gene].DiaphragmRank
		BrainRank = genescore[Gene].MouseBrainRank

		#print Proband
		Indi = Pedigree.GetIndi(Proband.name)
		Phenotype = Indi.PhenotypeDetail
		Relateness = Indi.Relateness

		#'Sample', 'Phenotype', 'Chrom', 'Pos', 'Ref', 'Alt', 'AC', 'Mappability', Gene', 'GeneName', 'GeneFunc', 'ExonicFunc', 'AAchange', 'ExAC', 'gnomAD', 'VarType', 'REVEL', 'MPC', 'MCAP', 'MetaSVM', 'CADD', 'PP2', '1KG', 'mis_z', 'lof_z','pLI', 'pRec','HeartRank', 'LungRank', 'BrainRank','Filter', 'QUAL', 'ProbandGT', 'FatherGT', 'MotherGT', 'Relateness'
		return [Proband.name, Phenotype, self.Chrom, self.Pos, self.Ref, self.Alt, ','.join(str(x) for x in self.Info['AC']), ','.join(self.Info.get("Mappability",".")), Gene, GeneName,','.join( self.Info['Func.refGene']), ','.join(self.Info['ExonicFunc.refGene']), ','.join(self.Info['AAChange.refGene']),','.join(str(AF(x)) for x in self.Info['ExAC_ALL']), ','.join(str(AF(x)) for x in self.Info['gnomAD_genome_ALL']), VarType, ','.join(self.Info['REVEL']), ','.join(self.Info.get("MPC",".")), ','.join(self.Info['MCAP']), ','.join(self.Info['MetaSVM_pred']),','.join(self.Info['CADD_phred']),','.join(self.Info['Polyphen2_HDIV_pred']) , ','.join(str(AF(x)) for x in self.Info['1000g2015aug_all']), mis_z, lof_z, pLI, pRec, HeartRank, LungRank, BrainRank, self.Filter, self.Qual, Proband.Format, Father.Format, Mother.Format ,Relateness   ]

	def GetVarType(self, GeneFunc, ExonicFunc, REVEL):
		if GeneFunc in ['splicing', 'exonic_splicing']:
			return "LGD"
		elif ExonicFunc in ['stoploss', 'stopgain', 'frameshift_insertion', 'frameshift_deletion', 'frameshift_block_substitution']:
			return "LGD"
		elif GeneFunc == "exonic" and ExonicFunc in ['nonsynonymous_SNV','unknown']:
			try:
				if float(REVEL) >= 0.5:
					return "D-mis"
			except:
				return "mis"
			else:
				return "mis"
		elif GeneFunc == "exonic" and ExonicFunc == 'synonymous_SNV':
			return 'silent'
		elif GeneFunc == "exonic" and ExonicFunc in ['nonframeshift_insertion', 'nonframeshift_deletion']:
			return 'inframe'
		else:
			print GeneFunc, ExonicFunc
			return '.'


class Individual():
	def __init__(self, List, Header):
		self.Fam, self.Sample, self.Father, self.Mother, self.Gender, self.Pheno = List[:6]
		#self.Relationship = List[Header.index('Relationship')]
		#self.PhenotypeDetail = List[Header.index('Disease')] + '/' + List[Header.index('DistinguishingFeatures')]
		try:
			self.PhenotypeDetail = List[Header.index('PhenotypeDetail')] 
		except:
			self.PhenotypeDetail = "."
		try:
			self.Relateness = List[Header.index('Relateness')]
		except:
			self.Relateness = "."
		try:
			self.SexCheck = List[Header.index('SexCheck')]
		except:
			self.SexCheck = "."
	def show(self):
		print self.Fam, self.Sample, self.Father, self.Mother

class Pedigree():
	def __init__(self, PedFil):
		fin = open(PedFil, 'r')
		self.individuals = []
		for l in fin:
			if l.startswith('#'):
				Header = l.strip().split('\t') 
			indi = l.strip().split('\t')
			indi = Individual(indi, Header)
			indi.show()
			self.individuals.append(indi)
		self.Proband, self.Father, self.Mother = None, None, None
		for ind in self.individuals:
			#if ind.Relationship == 'Proband':
			#	self.Proband = ind
			if ind.Father != "0" and ind.Mother != "0":
				self.Proband = ind
		for ind in self.individuals:
			if self.Proband.Father == ind.Sample:
				self.Father = ind
			if self.Proband.Mother == ind.Sample:
				self.Mother = ind
	def show(self):
		print 'Proband:%s\tFather:%s\tMother:%s' % (self.Proband.Sample, self.Father.Sample, self.Mother.Sample)
	def GetIndi(self, Sample):
		for Indi in self.individuals:
			if Indi.Sample == Sample:
				return Indi
		return None

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--dir", type=str, required=True, help="<Required> Dir contains family Ped and VCF files")
	parser.add_argument("-f", "--filter", type=str, help="ymal file contains filter criteria")
	parser.add_argument("-o", "--output", type=str, help="Out put file")
	args = parser.parse_args()
	if args.filter == None:
		args.filter = DENOVO_FILTER
	if args.output == None:
		args.output = "Denovo.csv"
	return args.dir, args.filter, args.output

def MatchVcfPed(vcfs, peds):
	res = []
	for vcf in vcfs:
		for ped in peds:
			if utils.GetBaseName(vcf) == utils.GetBaseName(ped):
				print utils.GetBaseName(vcf), utils.GetBaseName(ped)
				res.append((vcf, ped))
				break
	print "Finded vcf-ped pairs:"
	for v, p in res:
		print v, p
		print '\n'
	return res

def AF(CHR):
	if CHR == '.':
		return 0
	else:
		return float(CHR)

def For_one_vcf(VCF, Ped, Writer, Filters, CSV_HEADER, gene_score):
	Ped = Pedigree(Ped)
	fin = open(VCF, 'rb')
	for line in fin:
		if line.startswith("##"):
			# Meta info
			continue 
		elif line.startswith('#'):
			# header
			header = line.strip().split('\t')
		else:
			# variants line
			var = Variant(line, header)
			isDenovo, Proband, Father, Mother = var.CheckDeNovo(header, Ped, Filters)
			if isDenovo:
				#var.show()
				Writer.writerow(var.OutAsCSV(Proband, Father, Mother, Ped, gene_score))

def Call_DeNovo(InpDir, yaml_fname, outname):
	Filters = YML_Filter(yaml_fname)
	genescore = GENE_ANNOTATION()
	genescore = genescore.Genes
	vcfs = utils.get_files(InpDir, '.vcf')
	peds = utils.get_files(InpDir, '.ped')
	vcf_peds = MatchVcfPed(vcfs, peds)
	OutFil = open(outname,'wb')
	Writer = csv.writer(OutFil, delimiter=',')
	#CSV_HEADER[0] = '#'+CSV_HEADER[0]
	Writer.writerow(CSV_HEADER)
	for vcf, ped in vcf_peds:
		print "Processing vcf: %s\tpedigree: %s" % (vcf, ped)
		s_time = time.time()
		#if not vcf.endswith('CARE19.vcf'):
		#    continue
		For_one_vcf(vcf, ped, Writer, Filters, CSV_HEADER, genescore)
		print "%s Finished, used %.3fs"%(vcf, time.time()-s_time)


def main():
	dirc, yaml_fname, outname = GetOptions()
	dirc = os.path.abspath(dirc)
	if not dirc.endswith('/'):
		dirc += '/'
	print "Vcf list is", dirc
	Call_DeNovo(dirc, yaml_fname, outname)


if __name__ == '__main__':
	main()
