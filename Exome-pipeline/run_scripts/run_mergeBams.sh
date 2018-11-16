#!/bin/bash

REF=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
BAM=$1
BED=/home/local/users/jw/resources/GIAB/NA12878/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed
CMD=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.2.MergeBams.sh


nohup $CMD -r $REF -i $BAM &
