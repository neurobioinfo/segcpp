#! /bin/bash

module load StdEnv/2020  gcc/9.3.0; module load bcftools/1.16

read INPUT_VCF OUTPUT_VCF <<< $@
if [[ ! -s $INPUT_VCF ]]; then >&2 echo ""; failure=1; fi
if [[ $failure -eq 1 ]]; then >&2 echo ""; exit 42; fi
OUTPUT_VCF=${OUTPUT_VCF:-${INPUT_VCF/vcf/segcpp_formatted.vcf}}
OUTPUT_VCF=${OUTPUT_VCF%.gz}

n=$(bcftools +split-vep -l $INPUT_VCF |tail -1|cut -f1)
bcftools +split-vep --columns 1-$n --select worst $INPUT_VCF |bcftools annotate -x INFO/CSQ --output $OUTPUT_VCF -
