#!/bin/bash

REF=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
InpFil=/home/local/users/jw/CDH_MIPS/CDH_UW_control/CDH_UW.bam.list
BED=/home/local/users/jw/CDH_MIPS/Modified_bed.bed
CMD=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAln.8a.DepthofCoverage.sh

NJob=`wc -l $InpFil|cut -f 1 -d ' '`
echo "$NJob jobs.."
seq 1 $NJob | parallel -j 30 --eta sh $CMD -i $InpFil -r $REF -t $BED -a {} &

