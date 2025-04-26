#!/usr/bin/bash
source /home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/activate vep-99.1
vcf_files='/home/groups/mcweeney_lab/jengs/HNSCC/mutect2_vcf_files/mutect2_vcf_files/*vep.vcf'
for each_file in $vcf_files
do
  echo ${each_file}
  new_file=${each_file/.vcf/_for_aseq.vcf}
  echo ${new_file}
  srun -A mcweeney_lab --time=1440 -p exacloud filter_vep -i ${each_file}  --filter "gnomAD_AF > 0.05" | /home/exacloud/gscratch/mcweeney_lab/jengs/software/bcftools-1.20/bcftools view -i'FILTER~"germline"' |  /home/exacloud/gscratch/mcweeney_lab/jengs/software/bcftools-1.20/bcftools view -i'GT=="0/1"' > ${new_file}
done

