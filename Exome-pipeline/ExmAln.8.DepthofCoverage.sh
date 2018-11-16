#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N DoC 
#$ -l h_rt=80:00:00
#$ -l h_vmem=8G
#$ -cwd


#This script takes a bam file and generates depth of coverage statistics using GATK
#    InpFil - (required) - Path to Bam file to be aligned or a file containing a list of bam files one per line (file names must end ".list")
#            if it is a list then call the job as an array job with -t 1:n where n is the number of bams
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
#    LogFil - (optional) - File for logging progress
#    Flag - B - BadET - prevent GATK from phoning home
#    Flag - F - Fix mis-encoded base quality scores - see GATK manual. GATK will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $TARGET - exome capture intervals bed file or other target file (must end ".bed")
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
 (-t <X>-<Y> [if providing a list]) ExmAln.8a.DepthofCoverage.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -DCBH

     -i (required) - Path to Bam file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -t (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
     -l (optional) - Log file
     -D (flag) - keep full Depth of Coverage file
     -C (flag) - Allow BadCigar - see GATK documentation - allows reads that GATK interprets as indicating a malformed file, e.g. reads starting with a deletion
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

BadCigar="false"
BadEt="false"
FullDoC="false"

#get arguments
while getopts i:r:a:t:l:DCBFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        a) ArrNum="$OPTARG";; 
        t) TgtBed="$OPTARG";; 
        l) LogFil="$OPTARG";;
        D) FullDoC="true";;
        F) FixMisencoded="true";;
        C) BadCigar="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]] ; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil


if [[ -z "${ArrNum}" ]]
then
    ArrNum=$SGE_TASK_ID
fi

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"
funcGetTargetFile
#Set Local Variables
InpFil=`readlink -f $InpFil` #resolve absolute path to bam
BamFil=$(tail -n+$ArrNum $InpFil | head -n 1) 
BamNam=`basename $BamFil | sed s/.bam//` #a name to use for the various files
if [[ -z $LogFil ]];then LogFil=$BamNam.DoC.log; fi # a name for the log file
OutFil=$BamNam.DoC #prefix used in names of output files
GatkLog=$BamNam.DoC.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.DoC.temp.log #temporary log file 
TmpDir=$BamNam.DoC.tempdir; mkdir -p $TmpDir #temporary directory

#Start Log
ProcessName="Depth of Coverage with GATK" # Description of the script - used in log
funcWriteStartLog

#Calculate depth of coverage statistics
StepName="Calculate depth of coverage statistics using GATK DepthOfCoverage" # Description of this step - used in log
StepCmd="java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T DepthOfCoverage
 -R $REF
 -I $BamFil
 -o $OutFil
 -L $TgtBed
 -ct 1 
 -ct 5
 -ct 10
 -ct 15
 -ct 20
 --countType COUNT_FRAGMENTS 
 --omitDepthOutputAtEachBase
 --omitIntervalStatistics
 --omitLocusTable
 --minMappingQuality 20
 --filter_mismatching_base_and_quals
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#End Log
funcWriteEndLog
