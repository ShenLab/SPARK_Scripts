#Library of functions used throughout Exome analysis scriptd

#-------------------------------------------------------------------------------------------------------
#Function to set the target file location when given a code present as a variable in the reference file
funcGetTargetFile (){
    if [[ "$TGTCODES" == *"$TgtBed"* ]];then
        eval TgtBed=\$$TgtBed
    fi
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#Function to get input file name from a list of files in an array job
funcFilfromList() {
ChecList=${InpFil##*.}
if [[ "$ChecList" == "list" ]];then
    echo $ChecList
    InpFil=$(head -n $ArrNum $InpFil | tail -n 1)
fi
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#Function to enter information about the script initiation into the log
funcWriteStartLog () {
uname -a >> $TmpLog
echo "Start "$ProcessName" - $0:`date`" >> $TmpLog
echo " Job name: "$JOB_NAME >> $TmpLog
echo " Job ID: "$JOB_ID >> $TmpLog
if [[ -n "$SGE_TASK_ID" ]]; then echo " Array task ID: "$SGE_TASK_ID >> $TmpLog; fi
echo " Input File: "$InpFil >> $TmpLog
if [[ -n "$BamFil" ]]; then echo " Bam File: "$BamFil >> $TmpLog; fi
if [[ -n "$BamNam" ]]; then echo " Base name for outputs: $BamNam" >> $TmpLog; fi
if [[ -n "$VcfFil" ]]; then echo " Vcf File: "$VcfFil >> $TmpLog; fi
if [[ -n "$VcfNam" ]]; then echo " Base name for outputs: $VcfNam" >> $TmpLog; fi
if [[ -n "$Chr" ]]; then echo " Chromosome: "$Chr >> $TmpLog; fi
if [[ -n "$TgtBed" ]]; then echo " Target Intervals File: "$TgtBed >> $TmpLog; fi
if [[ -n "$RefFil" ]]; then echo " Pipeline Reference File: "$RefFil >> $TmpLog; fi
echo "----------------------------------------------------------------" >> $TmpLog
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function to log the start of each step wtihin a script
funcLogStepStart () { echo "- Start $StepName `date`...">> $TmpLog ; } 
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#func to trim GATK output log and write it to the temp log
funcTrimGATKlog (){
    echo "  --- GATK output log for $StepName ----------------" >> $TmpLog
    grep -vE "ProgressMeter - *[dc0-9XY]|Copyright|INITIALIZATION COMPLETE|----|For support and documentation|Done preparing for traversal|^WARN[[:print:]]*SnpEff|ReadShardBalancer|this tool is currently set to genotype at most 6 alternate alleles in a given context" $GatkLog | awk '{ print "\t\t"$0 }' >> $TmpLog
    echo "  --- --- --- --- --- --- ---" >> $TmpLog
    rm $GatkLog
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function checks that the step has completed successfully, if not it writes and error message to the log and exits the script, otherwise it logs the completion of the step
funcLogStepFinit () { 
if [[ $? -ne 0 ]]; then #check exit status and if error then...
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
    echo "     $StepName failed `date`" >> $TmpLog
    #qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID:" >> $TmpLog #get cluster usage stats
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
    echo "=================================================================" >> $TmpLog
    funcTrimGATKlog
    grep "ERROR MESSAGE" $SGE_STDERR_PATH | awk '{ print "\t\t"$0 }' >> $TmpLog
    cat $TmpLog >> $LogFil
    rm $TmpLog
    exit 1
fi
if [[ "$StepCmd" == *GenomeAnalysisTK.jar* ]]; then funcTrimGATKlog; fi
echo "- End $StepName `date`...">> $TmpLog # if no error log the completion of the step
echo "-----------------------------------------------------------------------" >> $TmpLog
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function to run and log the initiation/completion/failure of each step in the script
funcRunStep (){
funcLogStepStart
echo $StepCmd >> $TmpLog
eval $StepCmd
funcLogStepFinit
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function to log the end of each script and transfer the contents of temporary log file to the main log file
funcWriteEndLog () {
echo "End "$ProcessName" $0:`date`" >> $TmpLog
SrchTrm="usage"
if [[ "$SGE_TASK_ID" -gt 0 ]]; then SrchTrm=$SrchTrm" *$SGE_TASK_ID:"; fi
echo "===========================================================================================" >> $TmpLog
echo "" >> $TmpLog
cat $TmpLog >> $LogFil
rm -r $TmpLog $TmpDir
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# function for common additional arguments to GATK
funcGatkAddArguments (){
if [[ "$AllowMisencoded" == "true" ]]; then StepCmd=$StepCmd" -allowPotentiallyMisencodedQuals"; fi
if [[ "$FixMisencoded" == "true" ]]; then StepCmd=$StepCmd" -fixMisencodedQuals"; fi
if [[ "$BadCigar" == "true" ]]; then StepCmd=$StepCmd" -rf BadCigar"; fi
if [[ "$BadET" == "true" ]]; then StepCmd=$StepCmd" -et NO_ET -K $ETKEY"; fi
}

#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# function for calling next step in pipeline
funcPipeLine (){
if [[ "$PipeLine" == "true" ]]; then
    mkdir -p stdostde
    echo "- Call $NextJob `date`:" >> $TmpLog
    echo "    "$QsubCmd  >> $TmpLog
    eval $QsubCmd >> $TmpLog
    echo "----------------------------------------------------------------" >> $TmpLog
fi
}
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
