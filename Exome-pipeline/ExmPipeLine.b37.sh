#!/bin/bash
RefFil=$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
FastQC=$HOME/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.11.FastQC.sh
BWA=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAln.1a.Align_Fastq_to_Bam_with_BWAmem.sh
REMAP=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAln.1b.ReAlign_Bam_with_BWAmem.sh
MergeBAM=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
DoC=$HOME/CUMC/Exome-pipeline-Jiayao/ExmAln.8a.DepthofCoverage.sh
BAMQC=$HOME/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.17.CollectInsertSizeMetrics.sh
HapCaller=$HOME/CUMC/Exome-pipeline-Jiayao/ExmAln.2.HaplotypeCaller_GVCFmode.sh
JointGT=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.1hc.GenotypeGVCFs.sh
VCFMerge=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
VQSR=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.4.RecalibrateVariantQuality.sh
Annovar=$HOME/CUMC/Exome-pipeline-Jiayao/Test.AnnotateVCF_direct.sh

ProjectHome=
ProjectName=
BedFil=
FastQList=
RawBamList=
BamList=${ProjectHome}/src/${ProjectName}.bam.list 
GVCFList=${ProjectHome}/src/${ProjectName}.gvcf.list
SplitedDir=${ProjectHome}/JointGenotyping/${ProjectName}.splitfiles
RawVCF=${ProjectHome}/JointGenotyping/${ProjectName}.rawvariants.vcf.gz
VQSRVCF=${ProjectHome}/JointGenotyping/${ProjectName}.rawvariants.recalibrated.vcf
AnnovarVCF=${ProjectHome}/JointGenotyping/${ProjectName}.rawvariants.recalibrated.vcf.hg19_multianno.vcf.gz

#==========================================================================================================
#FastQC
FASTQCDIR=${ProjectHome}/FASTQC
mkdir -p $FASTQCDIR; cd $FASTQCDIR
Input=$BamList
NJob=`wc -l $Input|cut -f 1 -d ' '`
echo $NJob
seq $NJob | parallel -j 10 --eta $FastQC -i $BamList -r $RefFil -a {} -d `pwd`
#==========================================================================================================

#==========================================================================================================
#HaploytypeCaller
GVCF=${ProjectHome}/GVCF
mkdir -p $GVCF; cd $GVCF
NJob=`wc -l $BamList|cut -f 1 -d ' '`
echo $NJob
seq $NJob | parallel -j 15 --eta $HapCaller -i $BamList -r $RefFil -t $BedFil -a {}
find `pwd` -name '*.g.vcf.gz' > $GVCFList
#==========================================================================================================

#mkdir -p ${ProjectHome}/JointGenotyping
#==========================================================================================================
#Joint Genotyping
cd ${ProjectHome}/JointGenotyping
seq 10 | parallel -j 10 --eta sh $JointGT -i $GVCFList -r $RefFil -a {} -j 20 -t $BedFil -n $ProjectName
#==========================================================================================================

#==========================================================================================================
#BAM QC and STATS
BAMQC_DIR=${ProjectHome}/BamQC
mkdir -p $BAMQC_DIR; cd $BAMQC_DIR
NJob=`wc -l $BamList|cut -f 1 -d ' '`
echo $NJob
nohup seq $NJob | parallel -j 5 --eta $BAMQC -i $BamList -r $RefFil -a {} &
#==========================================================================================================


#==========================================================================================================
#DoC
DoC_DIR=${ProjectHome}/DoC
mkdir -p $DoC_DIR; cd $DoC_DIR
NJob=`wc -l $BamList|cut -f 1 -d ' '`
echo $NJob
nohup seq $NJob | parallel -j 5 --eta $DoC -i $BamList -r $RefFil -t $BedFil -a {} &
#==========================================================================================================

#==========================================================================================================
#MergeVCF
echo $SplitedDir
$VCFMerge -i $SplitedDir -r $RefFil 
#==========================================================================================================

#==========================================================================================================
#VQSR
echo $RawVcf
$VQSR -i $RawVCF -r $RefFil 
#==========================================================================================================

#==========================================================================================================
#Annovar
if [ -e $VQSRVCF ]
then
   echo $VQSRVCF
   echo $Annovar -i $VQSRVCF -r $RefFil 
else
	$Annovar -i $RawVCF -r $RefFil 
fi
#==========================================================================================================
