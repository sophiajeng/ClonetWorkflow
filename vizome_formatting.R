library(data.table)
library(stringr)

args <- commandArgs(TRUE)

.my.short.prots <- c(Ala="A", Arg="R", Asn="N", Asp="D", Asx ="B", Cys="C", Glu="E", Gln ="Q", Glx="Z", Gly="G", His= "H", Ile ="I", Leu= "L",
                     Lys ="K", Met= "M", Phe ="F", Pro= "P", Ser= "S", Thr= "T", Trp= "W", Tyr= "Y", Val= "V", Xxx ="X", Ter ="*" )


#taken from http://dec2015.archive.ensembl.org/info/genome/variation/predicted_data.html#consequences
.ens.cons.tab <- "transcript_ablation 	A feature ablation whereby the deleted region includes a transcript feature 	SO:0001893 	Transcript ablation
splice_acceptor_variant 	A splice variant that changes the 2 base region at the 3' end of an intron 	SO:0001574 	Essential splice site
splice_donor_variant 	A splice variant that changes the 2 base region at the 5' end of an intron 	SO:0001575 	Essential splice site
stop_gained 	A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript 	SO:0001587 	Stop gained
frameshift_variant 	A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three 	SO:0001589 	Frameshift coding
stop_lost 	A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript 	SO:0001578 	Stop lost
start_lost  A codon variant that changes at least one base of the canonical start codo	SO:0002012	Start lost
transcript_amplification 	A feature amplification of a region containing a transcript 	SO:0001889 	Transcript amplification
inframe_insertion 	An inframe non synonymous variant that inserts bases into in the coding sequence 	SO:0001821 	Non synonymous coding
inframe_deletion 	An inframe non synonymous variant that deletes bases from the coding sequence 	SO:0001822 	Non synonymous coding
missense_variant 	A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved 	SO:0001583 	Non synonymous coding
protein_altering_variant	A sequence_variant which is predicted to change the protein encoded in the coding sequence	SO:0001818	Protein altering variant
splice_region_variant 	A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron 	SO:0001630 	Splice site
incomplete_terminal_codon_variant 	A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed 	SO:0001626 	Partial codon
stop_retained_variant 	A sequence variant where at least one base in the terminator codon is changed, but the terminator remains 	SO:0001567 	Synonymous coding
synonymous_variant 	A sequence variant where there is no resulting change to the encoded amino acid 	SO:0001819 	Synonymous coding
coding_sequence_variant 	A sequence variant that changes the coding sequence 	SO:0001580 	Coding unknown
mature_miRNA_variant 	A transcript variant located with the sequence of the mature miRNA 	SO:0001620 	Within mature miRNA
5_prime_UTR_variant 	A UTR variant of the 5' UTR 	SO:0001623 	5prime UTR
3_prime_UTR_variant 	A UTR variant of the 3' UTR 	SO:0001624 	3prime UTR
non_coding_transcript_exon_variant 	A sequence variant that changes non-coding exon sequence in a non-coding transcript 	SO:0001792 	Within non coding gene
intron_variant 	A transcript variant occurring within an intron 	SO:0001627 	Intronic
NMD_transcript_variant 	A variant in a transcript that is the target of NMD 	SO:0001621 	NMD transcript
non_coding_transcript_variant 	A transcript variant of a non coding RNA gene 	SO:0001619 	Within non coding gene
upstream_gene_variant 	A sequence variant located 5' of a gene 	SO:0001631 	Upstream
downstream_gene_variant 	A sequence variant located 3' of a gene 	SO:0001632 	Downstream
TFBS_ablation 	A feature ablation whereby the deleted region includes a transcription factor binding site 	SO:0001895 	Tfbs ablation
TFBS_amplification 	A feature amplification of a region containing a transcription factor binding site 	SO:0001892 	Tfbs amplification
TF_binding_site_variant 	A sequence variant located within a transcription factor binding site 	SO:0001782 	Regulatory region
regulatory_region_ablation 	A feature ablation whereby the deleted region includes a regulatory region 	SO:0001894 	Regulatory region ablation
regulatory_region_amplification 	A feature amplification of a region containing a regulatory region 	SO:0001891 	Regulatory region amplification
feature_elongation 	A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence 	SO:0001907 	Feature elongation
regulatory_region_variant 	A sequence variant located within a regulatory region 	SO:0001566 	Regulatory region
feature_truncation 	A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence 	SO:0001906 	Feature truncation
intergenic_variant 	A sequence variant located in the intergenic region, between genes 	SO:0001628 	Intergenic"

.my.consequence.order <- function(){
  
  cons.lines <- strsplit(.ens.cons.tab, "\\n")
  
  return(sapply(strsplit(cons.lines[[1]], "\\s+"), "[", 1))
}

.fix.protein.ids <- function(csq.df){
  
  #approach adapted from vcf2maf
  
  proteins <- as.character(csq.df$HGVSp)
  
  # Remove transcript ID from HGVS codon/protein changes, to make it easier on the eye
  proteins <- sub("^.*:", "", proteins, perl=T)
  
  # Remove the prefixed HGVSc code in HGVSp, if found
  proteins <- sapply(regmatches(proteins, regexec("\\(*p\\.\\S+\\)*",proteins)), "[", 1)
  
  proteins <- gsub("[\\(\\)]", "", proteins)
  
  # Create a shorter HGVS protein format using 1-letter codes
  
  for (i in names(.my.short.prots)){
    proteins <- gsub(i, .my.short.prots[i], proteins)
  }
  
  # Fix HGVSp_Short,for splice acceptor/donor variants
  
  splice.pos <- csq.df$Variant_Classification %in% c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")
  
  if (sum(splice.pos) > 0){
    
    c.pos <- as.numeric(sapply(regmatches(as.character(csq.df$HGVSc), regexec("c\\.(\\d+)", as.character(csq.df$HGVSc))), "[", 2))
    
    c.pos <- ifelse(is.na(c.pos) ==F & c.pos < 0, 1, c.pos)
    
    proteins <- ifelse((splice.pos == T) & (is.na(c.pos) ==F) , paste0("p.X", sprintf("%.0f", (c.pos + c.pos %% 3)/3), "_splice"), proteins)
  }
  
  # Fix HGVSp_Short for Silent mutations, so it mentions the amino-acid and position
  
  p.pos <- as.numeric(sapply(regmatches(as.character(csq.df$Protein_position), regexec("^(\\d+)\\/\\d+$", as.character(csq.df$Protein_position))), "[", 2))
  
  fin.prots <- ifelse(csq.df$Variant_Classification == "synonymous_variant", paste0("p.", csq.df$Amino_acids, p.pos,csq.df$Amino_acids) , proteins)
  
  return(fin.prots)
  
}


.extract.csq <- function(vcf.dt, header, should.pick=T){
  
  info.split <- str_split(vcf.dt$INFO, "[;,]", n = Inf, simplify = FALSE)
  
  csq.lines <- lapply(seq_along(info.split), function(x){
    
    use.x <- info.split[[x]]
    
    which.csq <- grep("CSQ=", use.x)
    stopifnot(length(which.csq)==1)
    
    use.x[which.csq] <- sub("CSQ=", "",use.x[which.csq] )
    
    csq.end <- max(grep("\\|", use.x))
    
    list(csq=paste(paste0(x, "|", use.x[which.csq:csq.end]), collapse="\n"))
    
  })
  
  csq.dt <- fread(paste(sapply(csq.lines, "[[", "csq"), collapse="\n"), sep="|", header=F, colClasses = c(V1="character"))
  
  csq.info <- gsub("\"", "", header[grep("CSQ", header)])
  
  stopifnot(length(csq.info)==1)
  
  csq.header <- strsplit(regmatches(csq.info, regexec("Format:\\s+(.+)>", csq.info))[[1]][2], "\\|")[[1]]
  csq.header <- append("vcfRow", csq.header)
  
  names(csq.dt) <- csq.header
  
  if ("ALLELE_NUM" %in% names(csq.dt)==F){
    warning("ALLELE_NUM column not found, assuming no multi-allelic variants")
    csq.dt[,allele:="1"]
  }else{
    csq.dt[,allele:=as.character(ALLELE_NUM)]
  }
  
  csq.dt[,Allele:=NULL]
  
  if ("PICK" %in% names(csq.dt) && should.pick==T){
    
    csq.dt <- csq.dt[is.na(PICK)==F]
    
    #pick a most deletarious consequence
    tmp <- cbind(csq.dt[,.(vcfRow,allele)],csq.dt[,tstrsplit(Consequence, "&")])
    melt.csq <- melt(id.vars=c("vcfRow","allele"),data=tmp)[is.na(value)==F]
    melt.csq[, val_fac:=factor(value, levels=.my.consequence.order(), ordered=T)]
    
    var.class <- melt.csq[,.(Variant_Classification=value[order(val_fac)][1]),by=.(vcfRow,allele)]
    
    stopifnot(var.class[,.N] == csq.dt[,.N])
    
    csq.dt <- merge(csq.dt, var.class[,.(vcfRow, allele, Variant_Classification)], by=c("vcfRow","allele"))
    
  }else{
    
    csq.dt[,Variant_Classification:=sapply(strsplit(Consequence, "&"), function(x) x[order(factor(x, levels=.my.consequence.order(), ordered=T))][1])]
    
  }
  
  csq.dt[,HGVSp_short:=.fix.protein.ids(.SD)]
  
  csq.dt
}

.extract.gt <- function(vcf.dt, samples.to.include){
  
  if (is.null(names(samples.to.include))){
    sample.list <- setNames(samples.to.include, rep("tumor", length(samples.to.include)))
  }else{
    if (all(names(samples.to.include) %in% c("tumor","normal"))==F){
      stop("samples.to.include needs to be of the form: c(tumor='sample_1', normal='sample_2')")
    }else{
      sample.list <- samples.to.include
    }
  }
  
  split.format <- split(vcf.dt, by="FORMAT")
  
  vcf.gts <- do.call(rbind, lapply(split.format, function(x){
    
    use.x <- copy(x)
    
    rbindlist(lapply(sample.list, function(i){
      geno.names <- strsplit(use.x[["FORMAT"]][1], split=":")[[1]]
      #Some unified genotyper results with only a ./. or an incomplete record
      inp.str <- paste(ifelse(grepl("./.", use.x[[i]], fixed=T), paste("./.", paste(rep(NA, length(geno.names)-1), collapse=":"), sep=":") , use.x[[i]]), collapse="\n")
      #cause fread needs to see a '\n'
      if (grepl("\n", inp.str)==F){
        inp.str <- paste0(inp.str, "\n")
      }
      geno.info <- fread(inp.str, sep=":", header=F)
      names(geno.info) <- geno.names
      
      geno.counts <- melt(id.vars="vcfRow",data=cbind(use.x[,.(vcfRow)],geno.info[,tstrsplit(AD, ",")]), value.name = "parsed_alt_counts")
      geno.counts[,alt_reads:=as.integer(parsed_alt_counts)]
      geno.counts[,allele:=as.character(as.integer(sub("V","", variable))-1L)]
      geno.counts[,`:=`(variable=NULL, parsed_alt_counts=NULL)]
      geno.counts <- geno.counts[is.na(alt_reads)==F]
      
      cbind(samples=i, geno.counts)
      
    }), idcol="type")
    
  }))
  
  vcf.gts[,total_reads:=sum(alt_reads), by=.(samples, vcfRow)]
  
  vcf.list <- split(vcf.gts, by="type")
  
  if (length(vcf.list)==2){
    
    comb.geno <- merge(vcf.list$tumor, vcf.list$normal[,.(vcfRow, allele, normal_alt_reads=alt_reads, normal_total_reads=total_reads)], by=c("vcfRow","allele"))
    
  }else{
    comb.geno <- vcf.list[[1]][,.(vcfRow, allele,  type, samples, alt_reads, total_reads, normal_alt_reads=0L, normal_total_reads=0L)]
  }
  
  comb.geno[,type:=NULL]
  
  comb.geno
  
}


vcf.samples <- function(vcf.file){
  
  header <- system(paste0("head -n 1000 ",vcf.file," | grep '#' "), intern=T)
  
  header.line <- strsplit(header[length(header)], "\t")[[1]]
  
  header.line[(which(header.line=="FORMAT")+1):length(header.line)]
}


#need to have --output --vcf
#optional --somatic=[true,false] --exclude=[,]
cl_parse_vcf <- function(...){
  
  args <- list(...)
  
  base.args <- sapply(args, function(x) strsplit(sub("--", "", x), "=")[[1]])
  
  clean.args <- setNames(base.args[2,], base.args[1,])
  
  samps <- vcf.samples(as.character(clean.args["vcf"]))
  
  names(samps) <- rep("tumor", length(samps))
  
  #figure out normal sample from the vcf header
  #looks like normal_sample=7255537e-e793-42c4-b556-f1f4b3f5e342
  
  vcf.head <- system(paste0("head -n 1000 ",as.character(clean.args["vcf"])," | grep '#' "), intern=T)
  norm.samp <- vcf.head[grepl("normal_sample=", vcf.head)]
  
  if (length(norm.samp) > 0){
     names(samps)[samps == sub("##normal_sample=", "", norm.samp, fixed=T)] <- "normal"
  }
  
  if ("exclude" %in% names(clean.args)){
    exclude <- strsplit(clean.args["exclude"], ",")[[1]]
  }else{
    exclude <- character()
  }
  
  res.dt <- parse.vcf(vcf.file=as.character(clean.args["vcf"]), samples.to.include = samps, exclude.filter = exclude)
  
  save(res.dt, file=as.character(clean.args["output"]))
  
}

#assume the first is the outputname
cl_rbind <- function(...){
  
  args <- list(...)
  
  out.name <- args[[1]]
  
  args <- args[-1]
  
  sum.dt <- do.call(rbind, lapply(args, function(x){
    
    local({get(load(x))})
    
  }))
  
  save(sum.dt, file=out.name)
  
}


parse.vcf <- function(vcf.file, samples.to.include, exclude.filter=character()){
  
  header <- system(paste0("head -n 1000 ",vcf.file," | grep '#' "), intern=T)
  
  vcf.dt <- fread(vcf.file, sep="\t", skip=length(header)-1, header=T)
  
  if (missing(samples.to.include) || is.null(samples.to.include)){
    samples.to.include <- names(vcf.dt)[(which(names(vcf.dt)=="FORMAT")+1):ncol(vcf.dt)]
  }
  
  vcf.dt <- vcf.dt[sapply(strsplit(FILTER,";"), function(x) any(x %in% exclude.filter))==F]
  
  vcf.dt[,vcfRow:=as.character(seq_len(.N))]
  
  geno.dt <- .extract.gt(vcf.dt, samples.to.include)
  
  csq.dt <- .extract.csq(vcf.dt, header, should.pick=T)
  
  base.vcf <- vcf.dt[,.(vcfRow, seqnames=`#CHROM`, start_position=POS, end_position=POS+nchar(REF)-1, ref=REF, alt=ALT, filter=FILTER)]
  
  alt.dt <- melt(id.vars="vcfRow",data=cbind(base.vcf[,.(vcfRow)],base.vcf[,tstrsplit(alt, ",")]), value.name = "alt")
  alt.dt[,allele:=as.character(as.integer(sub("V","", variable)))]
  alt.dt[,variable:=NULL]
  alt.dt <- alt.dt[is.na(alt)==F]
  
  base.vcf[,alt:=NULL]
  
  vcf.dt.m <- merge(base.vcf, alt.dt, by=c("vcfRow"))
  
  vcf.dt.m <- merge(vcf.dt.m, geno.dt, by=c("vcfRow", "allele"))
  
  stopifnot(all(vcf.dt.m[,.N,samples][,N] == csq.dt[,.N]))
  
  vcf.dt.m <- merge(vcf.dt.m, csq.dt, by=c("vcfRow", "allele"))
  
  vcf.dt.m[,`:=`(t_vaf=ifelse(total_reads==0, 0, alt_reads/total_reads), n_vaf=ifelse(normal_total_reads==0, 0, normal_alt_reads/normal_total_reads))]
  
  ret.dt <- vcf.dt.m[,c("samples", "seqnames", "start_position", "end_position", "ref", "alt", "filter", "alt_reads", "total_reads", "normal_alt_reads", "normal_total_reads",
                        "t_vaf", "n_vaf", setdiff(names(csq.dt), c("vcfRow","allele"))), with=F]
  
  names(ret.dt) <- tolower(names(ret.dt))
  
  ret.dt
  
}

#can be invoked as a script for use with WDL files
if (length(args) > 0){
  if (exists(args[1])){
    do.call(args[1], as.list(args[-1]))
  }
}

