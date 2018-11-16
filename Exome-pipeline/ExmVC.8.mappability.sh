while getopts i:l:t:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        l) Length="$OPTARG";;
        t) BED="$OPTARG";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
#if [[ ! -z "$InpFil" ]] || [[ ! -z "$Length" ]] || [[ -z "$BED" ]]; then
# echo "Missing/Incorrect required arguments"
# echo "provided arguments: -i $InpFil -r $Length -t $BED"
# echo "usage: $usage"
# exit
#fi

PipePath="$HOME/CUMC/Exome-pipeline-Jiayao"
source $PipePath/mappability_generate.sh
function annotatedwithMappability {
   len=$1
   fMAP=$2
   fvcf=$3
   fbname=$(basename $fvcf)
   fbname=${fbname%.vcf.gz}
   fbname=${fbname%.vcf}
   fout="${fbname}.mappability.vcf.gz"
   echo $fMAP

  if [ ! -f $fMAP ]; then
         Redo_mappability $len
   fi
   keywd="Mappability"
   des="Mappability is obtained in file $fMAP"
   fhead="header_add.txt"

   echo "##INFO=<ID=$keywd,Number=1,Type=String,Description=\"$des\">" > $fhead
   bcftools annotate -a $fMAP  -c CHROM,FROM,TO,Mappability  -h  $fhead  $fvcf |bgzip -c > $fout
}


echo "input format annotatedwithMappability len *bed.gz  *.vcf.gz"
echo " bash $PipePath/mappability.annotation.sh 152 /home/local/ARCS/nz2274/Resources/mappability/hg19_152bp_mappability.bed.gz  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
annotatedwithMappability   $InpFil $Length $BED
if [ ! -f $2 ];then
   echo "$2 cannot find!\n"
fi

if [ ! -f $3 ]; then
  echo "$3 cannot find!\n"
fi
