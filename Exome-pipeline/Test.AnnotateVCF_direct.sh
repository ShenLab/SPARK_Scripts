#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N ANNOVAR 
#$ -l h_rt=256:00:00
#$ -l h_vmem=100G
#$ -cwd



#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#    InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    Flag - C - FullCadd - Annotate with full CADD database. The default is to use the caddgt10 database, which contains only variants with CADD scores that are within top 10% percentile. Using the full CADD database significantly increases the amount of time required for annotation, especially for larger vcfs (can mean the difference between 30 mins and several hours)
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $ANNHDB - directory containing databases for annovar

#list of required tools:
# annovar <http://www.openbioinformatics.org/annovar/> <http://www.openbioinformatics.org/annovar/annovar_download_form.php>
# N.B. : The perl script "table_annovar_cadd.pl", which is used below, is a modified version of the table_annovar.pl script that was released independent of the main bundle on 24th February 2014 (see annovar homepage).  The "_cadd" version has added lines to allow for the inclusion of the phred-scaled cadd score from the cadd or caddgt10 annovar databases. In the normal perl script only the raw cadd scores are added to the annotation.

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmVC.5.AnnotatewithANNOVAR.sh -i <InputFile> -r <reference_file> -l <logfile> -PH

     -i (required) - Path to VCF file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -C (flag) - Annotate with full CADD database
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -X (flag) - Do not run Variant Quality Score Recalibration
     -B (flag) - Prevent GATK from phoning home - only if calling pipeline
     -H (flag) - echo this message and exit
"

PipeLine="false"
FullCadd="false"
NoRecal="false"

while getopts i:r:t:l:m:CXFPBH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        t) Threads="$OPTARG";;
        l) LogFil="$OPTARG";;
		m) Mappability="$OPTARG";;
        C) FullCadd="true";;
        P) PipeLine="true";;
        X) NoRecal="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
#if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi
if [[ ! -e "$Threads" ]]; then
	Threads=1
fi
if [[ ! -e "$Mappability" ]]; then
	Mappability=$MAPPABILITY_FIL
fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set local Variables
##Set local parameters
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
VcfFil=$InpFil #input vcf file
VcfNam=`basename $VcfFil |  sed s/.gz$//| sed s/.vcf$// | sed s/.rawvariants$//` #basename for outputs
if [[ -z $LogFil ]]; then LogFil=$VcfNam.AnnVCF.log; fi # a name for the log file
TmpLog=$VcfNam.AnnVCF.temp.log #temporary log file
TmpVar=$VcfNam.tempvar
AnnFil=$VcfNam.annovar
AnnModRscript=$VcfNam.ModifyAnnovarTable.rscript.r
SnpEffFil=$VcfNam.SnpEff.vcf #SnEff annotations files
VcfFilAnn=$VcfNam.Ann.vcf # annovar annotated VCF output file
VcfFilSnF=$VcfNam.SnF.vcf # SnpEff annotated VCF output file
VcfFilOut=$VcfNam.annotated.vcf # final annotated output file

TmpDir=$VcfNam.AnnVCF.tempdir; mkdir -p $TmpDir #temporary directory

#check the vcf file to see if it is zipped 
FilTyp=${VcfFil##*.}

#Start Log File
ProcessName="Annotate VCF" # Description of the script - used in log
funcWriteStartLog


##Run Annovar to Annotate VCF file
StepName="Build Annotation table using ANNOVAR"
StepCmd="table_annovar.pl $VcfFil $ANNHDB --buildver $BUILD  --tempdir $TmpDir --remove -protocol refGene,gnomad_exome,gnomad_genome,1000g2015aug_all,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_afr,1000g2015aug_sas,exac03,dbnsfp33a,cadd13gt10,cosmic70,genomicSuperDups,mcap,revel,avsnp147 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,r,f,f,f -otherinfo  -nastring .  -vcfinput "
if [[ "$FullCadd" == "true" ]]; then 
    StepCmd=${StepCmd/cadd13gt10/cadd13}
    echo "  Using full CADD database..." >> $TmpLog
fi
funcRunStep
#mv $VcfFil.hg19_multianno.vcf $basename.hg19_multianno.vcf
VcfFilOut=$VcfFil.hg19_multianno.vcf

sed -i -e 's/\\x3d/:/g' $VcfFilOut
sed -i -e 's/\\x3b/-/g' $VcfFilOut
bgzip -f $VcfFilOut
tabix -f -p vcf $VcfFilOut.gz
rm $VcfFil.avinput
#
StepName="Chenge invalid char"
StepCmd="sed -i -e 's/\x3d/:/g' $VcfFilOut;
	     sed -i -e 's/\x3b/-/g' $VcfFilOut;
		 bgzip -f $VcfFilOut;
		 tabix -f -p vcf $VcfFilOut.gz;"

#End Log
#funcRunStep
funcWriteEndLog

#Cleanup

#rm -f $VcfNam*invalid_input $VcfNam*bylocus*
