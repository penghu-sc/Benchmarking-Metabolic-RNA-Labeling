#!/bin/bash
STAR_path=star_index
gtf=Danio_rerio.GRCz11.108.gtf
fq_dir=$1
outdir=$2
cd ${outdir}
mkdir -p 1_align 2_consensus 3_count_filter 4_estimate  5_count_unfilt estimate_ctrl count_ctrl
cd ${fq_dir}
fq1=$3
fq2=$4
cellnum=$5
name=$6

mkdir -p ${outdir}/1_align/${name}
# align
dynast align -i ${STAR_path} -t 40 -o ${outdir}/1_align/${name} -x dropseq ${fq2} ${fq1} --STAR-overrides " --soloCellFilter CellRanger2.2 ${cellnum} 0.99 10 --limitBAMsortRAM 18082363863" --tmp ${outdir}/1_align/${name}/tmp

if [ $? -eq 0 ]; then
  echo $name 'align done'
else
  echo $name "Error: align failed"
  exit 1  
fi
# count 
mkdir -p ${outdir}/5_count_unfilt/${name}
dynast count -t 40 -g ${gtf} --barcode-tag CB --umi-tag UB \
${outdir}/1_align/${name}/Aligned.sortedByCoord.out.bam \
-o ${outdir}/5_count_unfilt/${name} \
--conversion TC --gene-names --verbose \
--barcodes ${outdir}/1_align/${name}/Solo.out/Gene/filtered/barcodes.tsv --tmp ${outdir}/5_count_unfilt/${name}/tmp

if [ $? -eq 0 ]; then
  echo $name 'count done'
else
  echo $name "Error: count failed"
  exit 1  
fi

#consensus
mkdir -p ${outdir}/2_consensus/${name}
dynast consensus -t 40 -g ${gtf} \
--barcodes ${outdir}/1_align/${name}/Solo.out/Gene/filtered/barcodes.tsv \
--barcode-tag CB --umi-tag UB ${outdir}/1_align/${name}/Aligned.sortedByCoord.out.bam \
-o ${outdir}/2_consensus/${name} --tmp ${outdir}/2_consensus/${name}/tmp

if [ $? -eq 0 ]; then
  echo $name 'consensus done'
else
  echo $name "Error: consensus failed"
  exit 1  
fi

# 
if [[ $name == *"Ctrl"* ]]; then
  
  echo "Name contains 'Ctrl', proceeding with the steps."
  
    align_bam=${outdir}/2_consensus/${name}/consensus.bam
    
    #to get snp.csv
    mkdir ${outdir}/count_ctrl/${name}_for_snp
    dynast count --control -t 10 --snp-threshold 0.5 -o ${outdir}/count_ctrl/${name}_for_snp --conversion TC -g ${gtf} ${align_bam}
    if [ $? -eq 0 ]; then
      echo $name 'count snp.csv done'
    else
      echo $name "Error: count snp.csv failed"
      exit 1  
    fi
    #to get p_e.csv
    dynast estimate --control -o ${outdir}/estimate_ctrl/${name}_for_pe ${outdir}/count_ctrl/${name}_for_snp -t 10
    if [ $? -eq 0 ]; then
      echo $name 'estimate p_e.csv done'
    else
      echo $name "Error: estimate p_e.csv failed"
      exit 1  
    fi
else
  
  echo "Name does not contain 'Ctrl', skipping the steps."
fi

mkdir -p ${outdir}/3_count_filter/${name}
dynast count --snp-csv ${outdir}/count_ctrl/snps.csv \
-t 40 -g ${gtf} --barcode-tag CB --umi-tag UB \
${outdir}/2_consensus/${name}/consensus.bam \
-o ${outdir}/3_count_filter/${name} \
--conversion TC --gene-names --verbose \
--barcodes ${outdir}/1_align/${name}/Solo.out/Gene/filtered/barcodes.tsv --tmp ${outdir}/3_count_filter/${name}/tmp
 
if [ $? -eq 0 ]; then
  echo $name 'count done'
else
  echo $name "Error: count failed"
  exit 1  
fi

#estimate
dynast estimate ${outdir}/3_count_filter/${name} -t 35 \
-o ${outdir}/4_estimate/${name} \
--p-e ${outdir}/estimate_ctrl/p_e.csv --method alpha  --gene-names
