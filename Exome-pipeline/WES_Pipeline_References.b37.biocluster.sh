## Resource Directories
export EXOMPPLN="$HOME/CUMC/Exome-pipeline-Jiayao"
export EXOMRES="/share/data/resoures" # Directory containing resources/references for pipeline

#jar files  and directories for software
#GATKJAR="$HOME/bin/GenomeAnalysisTK.jar" #Current GATK jar file
GATKJAR="$HOME/bin/GenomeAnalysisTK-3.6.jar" #Current GATK jar file
#GATKJAR="$HOME/bin/gatk-package-4.0.5.0-local.jar" #Current GATK jar file

PICARD="$HOME/bin/picard.jar" #directory containing Picard jar files
SNPEFF="$HOME/bin/snpEff.jar" # Current snpEff jar file
FREEBAYES="$HOME/software_pkg/freebayes/bin/freebayes"
SAMTOOLS="$HOME/bin/samtools"
BCFTOOLS="$HOME/bin/bcftools"
PLATYPUS="$HOME/software_pkg/Platypus/bin/Platypus.py"

## References
export BUILD="hg19" # shorthand for build
export REF="$EXOMRES/reference_genomes/hg19/hg19.fasta" # human 1000 genome assembly from GATK
export HAPMAP="$EXOMRES/references/b37/hapmap_3.3.b37.vcf" # hapmap vcf from GATK
export INDEL="$EXOMRES/references/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf" # Gold standard INDEL reference from GATK
export TGVCF="$EXOMRES/references/b37/1000G_omni2.5.b37.sites.vcf" 
export INDEL1KG="$EXOMRES/references/b37/1000G_phase1.indels.b37.vcf" # INDEL reference from 1000 genomes
export DBSNP="$EXOMRES/references/b37/dbsnp_138.b37.vcf" # dbSNP vcf from GATK
export ONEKG="$EXOMRES/references/b37/1000G_phase1.snps.high_confidence.b37.vcf" # 1000 genome SNPs vcf
export ANNOVAR="$HOME/software_pkg/annovar"
export ANNHDB="/share/shenlab/ANNOVAR_DATA/humandb"

export STHSH="$EXOMRES/references/b37/stampy_b37" # hash file for Stampy - omit ".sthash" extension for compatibility with Stampy
export STIDX="$EXOMRES/references/b37/stampy_b37" # genome index file for Stampy - omit ".stidx" extension for compatibility with Stampy
#GATK no-phone-home key

#Capture Kit Target Files
