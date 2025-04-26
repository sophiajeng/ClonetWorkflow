#!/usr/bin/bash
while IFS=$'\t' read -r -a myArray
do
	echo "${myArray[0]/.vcf/.vep_for_aseq.vcf}"
	echo "${myArray[1]}"
	echo "${myArray[2]}"
	srun --input none -A mcweeney_lab --time=1440 -p exacloud /home/exacloud/gscratch/mcweeney_lab/jengs/software/ASEQ/aseq-v1.1.11-source/ASEQ vcf=${myArray[0]/.vcf/.vep_for_aseq.vcf} bam=/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/bam_files/${myArray[1]}.b37.bam out=/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/ASEQ
        srun --input none -A mcweeney_lab --time=1440 -p exacloud /home/exacloud/gscratch/mcweeney_lab/jengs/software/ASEQ/aseq-v1.1.11-source/ASEQ vcf=${myArray[0]/.vcf/.vep_for_aseq.vcf} bam=/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/bam_files/${myArray[2]}.b37.bam out=/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/ASEQ
done < /home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/file_for_aseq.txt

#/home/exacloud/gscratch/mcweeney_lab/jengs/software/ASEQ/aseq-v1.1.11-source/ASEQ vcflist=vcf.list bamlist=bam.list
#bam_files='/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/bam_files/*bam'
#for eachfile in $bam_files
#do
#        echo ${eachfile}
#	srun --mem=64G -A mcweeney_lab -p exacloud --time=1440 /home/exacloud/gscratch/mcweeney_lab/jengs/software/ASEQ/aseq-v1.1.11-source/ASEQ mode=PILEUP vcf=dbsnp.137.UCSC.Coding.vcf bam=${eachfile} out=ASEQ
#done
