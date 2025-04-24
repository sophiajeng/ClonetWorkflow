workflow implement_clonet {


     #should look like vcf_file\ttumour_name\tnormal_name

     Array[Array[File]] input_bams_vcfs = read_tsv("/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/file_for_wdl.txt")
     File rscript = "wdl_sequencing/vizome_formatting.R"
     String build = "GRCh37"
     String? addtl_flag
     #you will need to modify the clonet scripts for specific files and parameter
     File clonet_rscript="/home/exacloud/gscratch/mcweeney_lab/jengs/software/CLONET/CLONET_HNSCC.R"
     File clonet_config_rscript="/home/exacloud/gscratch/mcweeney_lab/jengs/software/CLONET/Examples/CLONETv2/ConfigurationFile/CLONET.Config.HNSCC.R"
     #aseq_output should be the same directory as in the clonet files
     String aseq_output = "/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/ASEQ_wdl/"
 
     scatter (input_bams_vcf in input_bams_vcfs){

        #call VEPAnnotate{
        #    input:
        #        input_vcf=input_bams_vcf[0],
        #        build = build,
        #        addtl_flag=addtl_flag
        #}

        #call ParseVEP{
        #    input:
        #        annot_vcf=VEPAnnotate.out_vcf,
        #        rscript=rscript
        #}
        #call filterVCF{
	#    input:
                #annot_vcf=VEPAnnotate.out_vcf
       #         annot_vcf=input_bams_vcf[0]
       #}
       call filterGnomAD{
	    input:
		annot_vcf=input_bams_vcf[0]

       }
       call filterGermline{
	    input:
                 annot_vcf=filterGnomAD.out_vcf
	}
	call filterSomatic {
		input:
			annot_vcf=filterGermline.out_vcf
	}
       call runASEQ{
            input:
                filt_vcf=filterSomatic.out_vcf,
                tumor_name=input_bams_vcf[1],
		normal_name=input_bams_vcf[2],
		aseq_output=aseq_output
       }
     }

     #call SummarizeRData{
     #   input:
     #       inp_rdata=ParseVEP.out_rdata,
     #       rscript=rscript
     #}
     call runClonet{
        input:
            runClonetScript = clonet_rscript,
            configClonetScript = clonet_config_rscript
     }
    output{
        #Array[File] out_vcfs = VEPAnnotate.out_vcf
        #Array[File] out_rdatas = ParseVEP.out_rdata
        #File sum_rdata= SummarizeRData.out_rdata
    }

}


task runClonet{
    File runClonetScript
    File configClonetScript

    command{
        set -e

        /home/exacloud/gscratch/mcweeney_lab/jengs/software/R-4.4.1/bin/Rscript ${runClonetScript} ${configClonetScript}
    }

    runtime {
        runtime_minutes: 10
        requested_memory_per_core: "5G"
        cpus: 1
        maxRetries: 3
    }

    #output{
    #    File out_rdata = "${out_rdata_str}"
    #}
}


task SummarizeRData{

    Array[File] inp_rdata
    File rscript
    String out_rdata_str = "summarized_genos.RData"

    command{
        set -e

        R CMD BATCH --vanilla '--args cl_rbind ${out_rdata_str} ${sep=" " inp_rdata}' ${rscript}
    }

    runtime {
        runtime_minutes: 10
        requested_memory_per_core: "5G"
        cpus: 1
        maxRetries: 3
    }

    output{
        File out_rdata = "${out_rdata_str}"
    }
}

task ParseVEP{

    File annot_vcf
    File rscript

    String out_rdata_str = basename(annot_vcf, ".vcf") + ".RData"

    command {
        set -e

        R CMD BATCH --vanilla '--args cl_parse_vcf --vcf=${annot_vcf} --output=${out_rdata_str}' ${rscript}
    }

    runtime {
        runtime_minutes: 100
        requested_memory_per_core: "3G"
        cpus: 1
        maxRetries: 3
    }

    output{
        File out_rdata="${out_rdata_str}"
    }

}

task VEPAnnotate{

    File input_vcf
    String build
    String? addtl_flag

    String out_vcf_str = basename(input_vcf, ".vcf") + ".vep.vcf"

    command {
        set -e

        vep -i ${input_vcf} \
        --fork 4 \
        --species homo_sapiens --assembly ${build} ${addtl_flag} --verbose \
        --no_stats --no_progress \
        --cache --offline --format vcf --vcf --force_overwrite --dir_cache /home/exacloud/gscratch/mcweeney_lab/resources/gatk_v4/vep  \
        --check_existing --clin_sig_allele 1 --hgvs \
        --sift b --polyphen b --ccds --symbol --numbers --canonical --protein --biotype --uniprot --tsl --af --af_1kg --af_esp --af_gnomad --variant_class \
        --total_length --allele_number --hgvsg --flag_pick_allele --xref_refseq --shift_hgvs 1 --no_escape --failed 1  \
        --pick_order canonical,tsl,biotype,rank,ccds,length \
        --output_file ${out_vcf_str}
    }
    runtime {
        runtime_minutes: 100
        requested_memory_per_core: "2G"
        cpus: 5
        maxRetries: 3
    }

    output {
        File out_vcf = "${out_vcf_str}"
    }

}

task filterGnomAD{

    File annot_vcf

    String out_vcf_str = basename(annot_vcf, ".vcf") + ".gnomAD_filt.vcf"

    command {
	filter_vep -i ${annot_vcf} --filter "gnomAD_AF > 0.05" -o ${out_vcf_str}
    }
    runtime {
        runtime_minutes: 100
        requested_memory_per_core: "2G"
        cpus: 5
        maxRetries: 3
    }

    output {
        File out_vcf = "${out_vcf_str}"
    }

}
task filterGermline{

    File annot_vcf

    String out_vcf_str = basename(annot_vcf, ".vcf") + ".germline.vcf"

    command {
        /home/exacloud/gscratch/mcweeney_lab/jengs/software/bcftools-1.20/bcftools view -i'FILTER~"germline"' ${annot_vcf} -o ${out_vcf_str}
    }
    runtime {
        runtime_minutes: 100
        requested_memory_per_core: "2G"
        cpus: 5
        maxRetries: 3
    }

    output {
        File out_vcf = "${out_vcf_str}"
    }

}

task filterSomatic{

    File annot_vcf
    
    String out_vcf_str = basename(annot_vcf, ".vcf") + ".somatic.vcf"

    command {
        /home/exacloud/gscratch/mcweeney_lab/jengs/software/bcftools-1.20/bcftools view -i'GT=="0/1"' ${annot_vcf} -o ${out_vcf_str}
    }
    runtime {
        runtime_minutes: 100
        requested_memory_per_core: "2G"
        cpus: 5
        maxRetries: 3 
    }  
    
    output {
        File out_vcf = "${out_vcf_str}"
    }   
        
}       

task runASEQ{

    File filt_vcf
    String tumor_name
    String normal_name
    String aseq_output
    
#    String aseq_out = basename(tumor_bam_file, ".bam") + ".for_aseq.vcf"

    command {
    	/home/exacloud/gscratch/mcweeney_lab/jengs/software/ASEQ/aseq-v1.1.11-source/ASEQ vcf=${filt_vcf} bam=/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/bam_files/${tumor_name}.b37.bam out=${aseq_output}
        /home/exacloud/gscratch/mcweeney_lab/jengs/software/ASEQ/aseq-v1.1.11-source/ASEQ vcf=${filt_vcf}  bam=/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/bam_files/${normal_name}.b37.bam out=${aseq_output}
    }
    runtime {
        runtime_minutes: 100
        requested_memory_per_core: "2G"
        cpus: 5
        maxRetries: 3
    }

    #output {
    #    File out_vcf = "${out_vcf_str}"
    #}

}

