#!/bin/bash

REF=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
INPUT=$1
BED=/home/local/users/jw/resources/GIAB/NA12878/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed
CMD=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAln.1a.Align_Fastq_to_Bam_with_BWAmem.sh

num=`wc -l $INPUT |cut -f 1 -d ' '`
echo "$num jobs"
#nohup $CMD -r $REF -i $INPUT &
seq 1 $num | parallel -j 20 --eta sh $CMD -i $INPUT -r $REF -a {} &
