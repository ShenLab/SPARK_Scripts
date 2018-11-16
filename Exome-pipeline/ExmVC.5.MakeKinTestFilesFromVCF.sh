#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N KINSHIP 
#$ -l h_rt=256:00:00
#$ -l h_vmem=50G
#$ -cwd

#This script takes a VCF file and generates files to check the familial relationships and sex of the samples
#    InpFil - (required) - Path to VCF file or a list of VCF Files to be recalibrated
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    OutNam - (optional) - Name for output files
#    LogFil - (optional) - File for logging progress
#    Help - H - (flag) - get usage info#rmation

#list of required vairables in reference file:
# None

#list of required tools:
# samtools
# vcftools
# R

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="ExmVC.3.RecalibrateVariantQuality.sh -i <InputFile> -r <reference_file> -o <output_file> -l <logfile> -PABH

     -i (required) - Path to VCF file
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -o (optional) - Base for output file names <default is VCF file name>
     -l (optional) - Log file
     -H (flag) - echo this message and exit
"

while getopts i:r:o:l:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        o) OutNam="$OPTARG";; 
        l) LogFil="$OPTARG";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
echo $InpFil
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]];  then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil


#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#set local variables
FILTER_COMMON=$HOME/CUMC/Exome-Filters-Jiayao/Filter_WGS_common.py
InpFil=`readlink -f $InpFil` #resolve absolute path to vcf
VcfFil=$(basename $InpFil | sed s/.gz//g| sed s/.vcf//g)
VcfFil="${VcfFil}.common.vcf"
FilterStep="python $FILTER_COMMON -v $InpFil -o $VcfFil"
eval $FilterStep
if [[ ! $OutNam ]]; then OutNam=`basename $InpFil | sed s/.gz$// | sed s/.vcf$// | sed s/.hardfiltered$// `; fi
if [[ ! $LogFil ]];then LogFil=$OutNam.kinshipfiles.log; fi
TmpLog=Temp.$LogFil
TmpRscript=Temp.rPlotHist.$OutNam.R # temporary R script to generate sex histogram
TmpRelatScript=Temp.Relatedeness.$OutNam.sh # temporary bash script to run vcftools relatedness2
TmpMissingScript=Temp.Missingness.$OutNam.sh # temporary bash script to run vcftools missingness
mkdir -p stdostde # output directory for std error/output from vcftools bash scripts
mkdir -p KinshipFiles # output directory

SampList=Temp.$OutNam.`date "+%j%y%H%M%N"`.samplist #file to hold list of samples in vcf
MrgList=Temp.$OutNam.splitplink.merge.list #a file to hold the merge list for plink if it is necessary to convert the vcf in steps

#check the vcf file to see if it is zipped 
FilTyp=${VcfFil##*.}

#Start Log File
ProcessName="Generate pedigree analysis files" # Description of the script - used in log
funcWriteStartLog

#get relatedness via KING algorithm
StepName="Get relatedness via KING algorithm in vcftools" # Description of this step - used in log
echo "#!/bin/bash" > $TmpRelatScript
echo "vcftools --vcf $VcfFil --relatedness2 --out KinshipFiles/$OutNam" >> $TmpRelatScript
if [[ $FilTyp == "gz" ]]; then 
    sed s/--vcf/--gzvcf/g $TmpRelatScript > $TmpRelatScript.2
    mv -f $TmpRelatScript.2 $TmpRelatScript
fi
chmod 775 $TmpRelatScript
StepCmd="./$TmpRelatScript > $TmpRelatScript.o 2>&1 "
funcRunStep

#get missingness for each individual
StepName="Get missingness with vcftools" # Description of this step - used in log
echo "#!/bin/bash" > $TmpMissingScript
echo "vcftools --vcf $VcfFil --missing-indv --out KinshipFiles/$OutNam" >> $TmpMissingScript
if [[ $FilTyp == "gz" ]]; then 
    sed s/--vcf/--gzvcf/g $TmpMissingScript > $TmpMissingScript.2
    mv -f $TmpMissingScript.2 $TmpMissingScript
fi
chmod 775 $TmpMissingScript
StepCmd="./$TmpMissingScript  > $TmpMissingScript.o 2>&1 "
funcRunStep


#convert to plink - the vcftools vcf --> plink tool cannot run with > 1000 samples on the cluster due to limitations on temporary files. Therefore, for larger cohorts, it is necessary to split and remerge the plink files.
less $VcfFil  | grep -m 1 "^#CHROM" | cut -f 10- | tr "\t" "\n" > $SampList
LEN=`cat $SampList | wc -l`
echo "There are $LEN samples in the vcf" >> $TmpLog
if [[ $LEN -lt 1000 ]]; then
    StepName="Convert vcf to plink ped/map using vcftools"
    StepCmd="vcftools --vcf $VcfFil --plink --out Temp.$OutNam"
    if [[ $FilTyp == "gz" ]]; then StepCmd=`echo $StepCmd | sed s/--vcf/--gzvcf/g`; fi
    funcRunStep
    StepName="Convert ped/map to bed/bim/fam using plink"
    StepCmd="plink --file Temp.$OutNam --make-bed --out KinshipFiles/$OutNam"
    funcRunStep
else
    split -l 900 $SampList $SampList.split.
    echo "Convert the vcf in "`ls | grep $SampList.split | wc -l`" steps then merge:"
    for i in $SampList.split*; do
        NAM=${i##*split.}
        StepName="Convert vcf to plink ped/map using vcftools - step $NAM"
        StepCmd="vcftools --vcf $VcfFil --keep $i --plink --out Temp.$OutNam.splitplink.$NAM"
        if [[ $FilTyp == "gz" ]]; then StepCmd=`echo $StepCmd | sed s/--vcf/--gzvcf/g`; fi
        funcRunStep
    done
    ls -r | grep splitplink | grep -E "ped$|map$" | paste - - > $MrgList
    echo "Merge list:" >> $TmpLog
    cat $MrgList >> $TmpLog
    StepName="Merge ped/map split files and convert to bed/bim/fam using plink"
    StepCmd="plink --merge-list $MrgList --make-bed --out KinshipFiles/$OutNam"
    funcRunStep
fi

#split the pseudo-autosomal X into "XY" for sex check
StepName="Split the pseudo-autosomal X into XY for sex check with plink"
StepCmd="plink --bfile KinshipFiles/$OutNam --split-x hg19 no-fail --make-bed --out KinshipFiles/$OutNam"
funcRunStep

#run sex check
StepName="Run imputation of sex using plink"
StepCmd="plink --bfile KinshipFiles/$OutNam --impute-sex --make-bed --out KinshipFiles/$OutNam.sexcheck"
funcRunStep
mv KinshipFiles/$OutNam.sexcheck.sexcheck KinshipFiles/$OutNam.sexcheck


StepName="Plot histogram of F statistics"
echo "png(\"KinshipFiles/$OutNam.plink_sexcheck_histogram.png\")
 hist(read.table(\"KinshipFiles/$OutNam.sexcheck\", header=T)[,\"F\"], main=\"$OutNam plink sexcheck - F statistic\", xlab=\"F\")
 dev.off()" > $TmpRscript
StepCmd="Rscript $TmpRscript"
funcRunStep



#Run IBD with rare variants
StepName="Run IBD with rare variants using plink"
StepCmd="plink --bfile KinshipFiles/$OutNam --maf 0.01 --genome --out KinshipFiles/$OutNam"
funcRunStep

#get missingness
StepName="Get missingness with plink"
StepCmd="plink --bfile KinshipFiles/$OutNam --missing --out KinshipFiles/$OutNam"
funcRunStep

#get missingness
StepName="Run relationship Inference by KING"
StepCmd="king -b KinshipFiles/$OutNam.bed --kinship --prefix $OutNam.kin"
funcRunStep

#End Log
funcWriteEndLog

#Clean up
rm -f Temp.$OutNam.* $OutNam*.sexcheck.* $OutNam*nosex $OutNam.ped $OutNam.map  $OutNam.log $OutNam.lmiss $OutNam*~ $TmpRscript KinshipFiles/*~ $TmpRelatScript $TmpMissingScript
