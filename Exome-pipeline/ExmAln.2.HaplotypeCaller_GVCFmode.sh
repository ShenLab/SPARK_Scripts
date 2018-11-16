#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N HaplotypeCaller_GVCFmode
#$ -l h_rt=108:00:00
#$ -l h_vmem=20G
#$ -cwd

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#    InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    LogFil - (optional) - File for logging progress
#    Flag - B - BadET - prevent GATK from phoning home
#    Flag - F - Fix mis-encoded base quality scores - see GATK manual. GATK will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $DBSNP - dbSNP vcf from GATK
# $HAPMAP - hapmap vcf from GATKf
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>
# bgzip 
# tabix

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
(-t <X>-<Y> [if providing a list]) ExmVC.1.HaplotypeCaller_GVCFmode.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -FBH

     -i (required) - Path to Bam file for variant calling or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -t (required) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability)
     -l (optional) - Log file
     -B (flag) - Prevent GATK from phoning home
     -F (flag) - Fix mis-encoded base quality scores - see GATK manual
     -H (flag) - echo this message and exit
"

BadET="false"
FixMisencoded="false"

while getopts i:r:a:n:l:t:BFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        a) ArrNum="$OPTARG";; 
        l) LogFil="$OPTARG";;
        t) TgtBed="$OPTARG";;
        B) BadET="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]] || [[ -z "$TgtBed" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi


if [[ -z "${ArrNum}" ]]
then
    ArrNum=$SGE_TASK_ID
fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set local Variables
funcGetTargetFile
InpFil=`readlink -f $InpFil` #resolve absolute path to bam
BamFil=$(tail -n+$ArrNum $InpFil | head -n 1) 
BamNam=`basename $BamFil | sed s/.bam//`
BamNam=${BamNam/.bam/} # a name for the output files
if [[ -z $LogFil ]]; then LogFil=$BamNam.HCgVCF.log; fi # a name for the log file
VcfFil=$BamNam.g.vcf #Output File
GatkLog=$BamNam.HCgVCF.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.HCgVCF.temp.log #temporary log file
TmpDir=$BamNam.HCgVCF.tempdir; mkdir -p $TmpDir #temporary directory
infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A HaplotypeScore -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions -A FisherStrand -A InbreedingCoeff -A ClippingRankSumTest -A DepthPerSampleHC" #Annotation fields to output into vcf files

echo "Reference Genome File is $REF"

#Start Log File
ProcessName="Genomic VCF generatation with GATK HaplotypeCaller" # Description of the script - used in log
funcWriteStartLog

##Run genomic VCF generation
StepName="gVCF generation with GATK HaplotypeCaller"
StepCmd="java -Xmx16G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T HaplotypeCaller
 -R $REF
 -L $TgtBed
 -I $BamFil
 --genotyping_mode DISCOVERY
 -stand_emit_conf 10
 -stand_call_conf 30
 --emitRefConfidence GVCF
 --variant_index_type LINEAR
 --variant_index_parameter 128000
 -o $VcfFil
 -D $DBSNP
 --comp:HapMapV3 $HAPMAP 
 -pairHMM VECTOR_LOGLESS_CACHING
 -rf BadCigar
 $infofields
 --filter_mismatching_base_and_quals
 --interval_padding 100
 --dontUseSoftClippedBases
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#--dontUseSoftClippedBases 


##gzip and index the gVCF
StepName="gzip and index the gVCF"
StepCmd="bgzip -f $VcfFil; tabix -f -p vcf $VcfFil.gz"
funcRunStep
rm $VcfFil.idx

#End Log
funcWriteEndLog

#End
