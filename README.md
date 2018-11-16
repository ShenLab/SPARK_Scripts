# SPARK_Scripts
Scripts used in SPARK project
For variant calling/de novo filtering/annotation
# Variant Calling
Variant Calling pipeline based on GATK Best Practices v3.6 (https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145). 
Exome-pipeline
Dependencies are stored in WES_Pipeline_References.b37.biocluster.sh
Functions are stored in exome.lib.sh 
## 1. BWA&MarkDuplicates .fastq -> .bam
	ExmAln.1a.Align_Fastq_to_Bam_with_BWAmem.sh (BWA V0.7.15; picard V2.7.1)
##2. HaplotypeCaller .bam -> .g.vcf
	ExmAln.2.HaplotypeCaller_GVCFmode.sh (GATK V3.6)
##3. JointCall with GenotypeGVCFs .g.vcf -> vcf
	ExmVC.1hc.GenotypeGVCFs.sh  (GATK V3.6)
	ExmVC.2.MergeVCF.sh
##4. Variant Quailty Recalibration with VQSR (this step seems already been removed from current version of best practices)
	ExmVC.3.RecalibrateVariantQuality.sh (GATK V3.6)
##5. Annotation with ANNOVAR
	ExmVC.4.AnnotateVCF.sh
##6. Some QC
	ExmVC.5.MakeKinTestFilesFromVCF.sh (Kinship test, vcftools v0.1.15, PLINK1.9)
	ExmAln.8.DepthofCoverage.sh (D15 vs Mean)
	ExmQC.VCFsummary_stats.py (some variants statistics for vcf)

# De novo Filtering
Exome-Filters
##1. Seperate cohort VCF into trio
	adhoc.1.vcf2trio.py

##2. Filter by Genotype and other heuristic filtering
	Get_Denovo.py #this script take trio vcfs and produce de novo variants according to filters stored in .yml files
	DENOVO_FILTER.yml # More stringent filters
	DENOVO_FILTER.Tier2.yml # Less stringent filters

##3. IGV inspection
	https://github.com/ShenLab/igv-classifier.git

