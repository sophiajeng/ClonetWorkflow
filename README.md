# Clonet Workflow

Clonet takes in as input pileup files generate from ASEQ.

Pileups should be generated relative to heterozygous SNPs in the normal sample.  

To create these vcf files, we needed to filter each Mutcet2 vcf file for 'germline' and normal genotype (GT) is '0/1' (use filter_mutect_files.sh)

We further restricted by annotating with the gnomad frequency and selecting those with frequency > 5%.  This was done using VEP.

Once we had vcfs for ASEQ, run ASEQ using ASEQ.sh.  This takes in as input the normal and tumour sample file names (files_for_aseq.txt)

After the mpileups from ASEQ have been created, use clonet.wdl to run Clonet.  Be sure to modify parameters like file locations.
