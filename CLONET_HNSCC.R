#!/usr/bin/env Rscript

######################################
## Load configuration file
cmd_args <- commandArgs(trailingOnly = T)
if (length(cmd_args) != 1){
  cat(paste0("[",Sys.time() ,"] ERROR - Usage: CLONET.R <config_file>\n"))
  quit()
}

######################################
## Plot intro message
ClonET_version<-"2.0.0"
cat("\n")
cat(paste0("[",Sys.time() ,"] CLONET v",ClonET_version,"\n\n",sep=""))
#cat(paste("[",Sys.time() ,"] Working directory: ",getwd(),"\n",sep=""))

######################################
## Load configuration file
configFile <- cmd_args[1]
cat(paste0("[",Sys.time() ,"] Load configuration file: ",configFile,"\n",sep=""))
source(configFile)
## check if defined bamChr parameters
if (!exists("bamChr")){
  bamChr <- F
}

######################################
## Check R version
minMajorVersion <- 2
minMinorVersion <- 7 
currect_R_version <- R.version
if ( ( as.numeric(currect_R_version$major) < minMajorVersion)  || 
     ( as.numeric(currect_R_version$major) == minMajorVersion && as.numeric(currect_R_version$minor) < minMinorVersion) ){
  cat(paste0("[",Sys.time() ,"] ERROR - CLONET require at least version ",minMajorVersion,".",minMinorVersion," but you have version ",currect_R_version$major,".",currect_R_version$minor, "\n"))
  quit()
}

######################################
## Check installed packages
cat(paste0("[",Sys.time() ,"] Check installed packages\n",sep=""))
pkg.list = c("parallel","pso","dgof","sets")
pkg.list.match=match(pkg.list,installed.packages()[,"Package"],nomatch=FALSE)
pkg.list.match = !as.logical(pkg.list.match)
if(sum(pkg.list.match)>0) { install.packages(pkg.list[pkg.list.match],repos="http://cran.r-project.org") } 
library(parallel,warn.conflicts=F)
library(pso,warn.conflicts=F)
library(dgof,warn.conflicts=F)
library(sets,warn.conflicts=F)
options(scipen=99999)

######################################
## Load CLONET basic functions
cat(paste0("[",Sys.time() ,"] Load CLONET basic functions in ", path_to_CLONET_functions,"\n",sep=""))
source(paste0(path_to_CLONET_functions,"Functions/CLONET.BasicFunctions.R"))
source(paste0(path_to_CLONET_functions,"Functions/CLONET.extendedBetaTable.Functions.R"))


######################################
## Check output directory
cat(paste("[",Sys.time() ,"] Output directory ",output_DIR,"\n",sep=""))
if (!file.exists(output_DIR)){
  cat(paste("[",Sys.time() ,"] WARNING - Output directory does not exit and will be created\n",sep=""))
  dir.create(path=output_DIR)
}
cat(paste("[",Sys.time() ,"] Create subfolders\n",sep=""))
analysis_dir <- paste0(output_DIR,"RData/")
if (!file.exists(analysis_dir)){
  dir.create(path=analysis_dir)
}
results_dir <- paste0(output_DIR,"Results/")
if (!file.exists(results_dir)){
  dir.create(path=results_dir)
}

#####################################
## plot parameters
cat("\n\n######################################################################\n")
#cat(paste("[",Sys.time() ,"] Sample info file ",sampleInfoFile,"\n",sep=""))
#cat(paste("[",Sys.time() ,"] Segments file ",segmentListFile ,"\n",sep=""))
#cat(paste("[",Sys.time() ,"] Pileup folder ",pileup_dir,"\n",sep=""))
cat(paste("[",Sys.time() ,"] Minimum per base tumor coverage ",minCoverage,"\n",sep=""))
cat(paste("[",Sys.time() ,"] Number of samples processed in parallel ", NsamplesToProcessInParallel,"\n",sep=""))
cat(paste("[",Sys.time() ,"] Number of cores per sample ", perSampleCores,"\n",sep=""))

#####################################
## Computing per sample beta tables
cat("\n\n######################################################################\n")
if (sampleInfoFile == ""){
  cat(paste("[",Sys.time() ,"] Sample info file not provided\n",sep=""))
}else{
  cat(paste("[",Sys.time() ,"] Sample info file ",sampleInfoFile,"\n",sep=""))
  samplesInfo <- read.csv(sampleInfoFile,sep="\t",header=T,check.names=F,as.is=T)  
  samplesToAnalyse <- samplesInfo$Tumor.Bam.Name
}

if (segmentListFile == ""){
  cat(paste("[",Sys.time() ,"] Segments file not provided\n",sep=""))
}else{
  
  ## determine if bed or seg
  FileExt <- substr(segmentListFile, nchar(segmentListFile)-3+1, nchar(segmentListFile))
  if (! FileExt %in% c("seg","bed")){
    stop("File segmentListFile has to be extension seg or bed but has extension ",FileExt)
  }
  mutationList <- read.csv(segmentListFile ,sep="\t",header=F,check.names=F,as.is=T)  
  
  ## if the file has the header reload with the header
  if (mutationList[1,1] != mutationList[2,1]){
    mutationList <- read.csv(segmentListFile ,sep="\t",header=T,check.names=F,as.is=T)  
  }
  
  ## rename the columns
  if (FileExt == "seg"){
    colnames(mutationList)  <- c("sample","chr","start","end","NotUsed","log2","NotUsed2")
  }else{
    colnames(mutationList)  <- c("chr","start","end","log2","sample")
  }
  cat(paste("[",Sys.time() ,"] Segments file ",segmentListFile ,"\n",sep=""))
}

if (pileup_dir == ""){
  cat(paste("[",Sys.time() ,"] Pileup folder not provided\n",sep=""))
}else{
  cat(paste("[",Sys.time() ,"] Pileup folder ",pileup_dir,"\n",sep=""))
}



if (1 %in% stages){
  
  if (!exists("betaCompute.method") || betaCompute.method == "GB"){
    cat(paste("[",Sys.time() ,"] Stage 1 - Compute beta table with GB method\n",sep=""))
    source(paste0(path_to_CLONET_functions,"Functions/CLONET.SampleAnalysis.R") )
  }else if (betaCompute.method == "STM"){
    cat(paste("[",Sys.time() ,"] Stage 1 - Compute beta table with STM method\n",sep=""))
    source(paste0(path_to_CLONET_functions,"Functions/CLONET.SampleAnalysis.STMmethod.HNSCC.R") )
  }else{
    stop("Beta computing method ",betaCompute.method," not defined.")
  }
  
}else{
  cat(paste("[",Sys.time() ,"] Skip stage 1 - Compute beta table\n",sep=""))  
}



#####################################
## Aggregate beta tables
cat("\n\n######################################################################\n")

AnalysisFile_suffix <- ".contaminationAnalysis.RData"
ListAnalysisFiles <- paste0(analysis_dir,list.files(analysis_dir,pattern=AnalysisFile_suffix))
betaTable_file <- paste0(results_dir,"betaTable.txt")

if (2 %in% stages){
  
  if (!exists("betaCompute.method") || betaCompute.method == "GB"){
    cat(paste("[",Sys.time() ,"] Stage 2 - Aggregate beta tables\n",sep=""))
    source(paste0(path_to_CLONET_functions,"Functions/CLONET.AggregateBetaValuesFromAnalysisFiles.R"))
  }else if (betaCompute.method == "STM"){
    cat(paste("[",Sys.time() ,"] Skip stage 2 - STM method create betaTable directly... Saving it...\n",sep=""))
    write.table(betaTable,betaTable_file ,quote=F,sep="\t",row.names=F,col.names=T)
  }
  
}else{
  cat(paste("[",Sys.time() ,"] Skip stage 2 - Aggregate beta tables\n",sep=""))
  cat(paste("[",Sys.time() ,"]      Use beta table ",betaTable_file,"\n",sep=""))
  betaTable <- read.delim(betaTable_file, as.is=T)
}

#####################################
## Ploidy correction
cat("\n\n######################################################################\n")

ploidyTable_file <- paste0(results_dir,"ploidyTable.txt")

if (3 %in% stages){
  cat(paste("[",Sys.time() ,"] Stage 3 - Compute ploidy correction\n",sep=""))
  source(paste0(path_to_CLONET_functions,"Functions/CLONET.ComputePloidyShift.R"))
}else{
  cat(paste("[",Sys.time() ,"] Skip stage 3 - Compute ploidy correction\n",sep=""))
  cat(paste("[",Sys.time() ,"]      Use ploidy table ",ploidyTable_file,"\n",sep=""))
  ploidyTable <- read.delim(ploidyTable_file,as.is=T)
}


#####################################
## Global admixture
cat("\n\n######################################################################\n")

errorTable <- read.csv(errorTable_file,sep="\t",check.names=F,as.is=T,header=T)
admTable_file <- paste0(results_dir,"globalAdmTable.txt")


if (4 %in% stages){
  cat(paste("[",Sys.time() ,"] Stage 4 - Compute global DNA admixture\n",sep=""))  
  if ( !exists("adm.global.method") ){ # for compatibility with old config file
    source(paste0(path_to_CLONET_functions,"Functions/CLONET.CreateGlobalAdmTable.R"))  
  }else if ( adm.global.method == "1D" ){
    source(paste0(path_to_CLONET_functions,"Functions/CLONET.CreateGlobalAdmTable.R"))  
  }else if (adm.global.method == "2D" ){
    source(paste0(path_to_CLONET_functions,"Functions/CLONET.CreateGlobalAdmTable.2D.R"))  
  }else if (adm.global.method == "SNV" ){
    ## check if variables are defined (for back compatibility)
#     if ( exists("PMreadCounts_file") && exists("minCov_SNV") && 
#          exists("minAltReads_SNV") && exists("minAltReads_SNV" &&
#          exists("bwRange"))){
      PMreadCounts <- read.delim(PMreadCounts_file, as.is=T)
      source(paste0(path_to_CLONET_functions,"Functions/CLONET.CreateGlobalAdmTable.SNV.R"))  
  }else{
      stop("Adm.global computing modality ",adm.global.method," not defined.")
  }
}else{
  cat(paste("[",Sys.time() ,"] Skip stage 4 - Compute global DNA admixture\n",sep=""))
  cat(paste("[",Sys.time() ,"]      Use global DNA admixture table ",admTable_file,"\n",sep=""))
  globalAdm <- read.delim(admTable_file,as.is=T)
}

#####################################
## Clonality table
cat("\n\n######################################################################\n")

clonalityTable_file <- paste0(results_dir,"clonalityTable.txt")
reportFile <- paste0(results_dir,"clonality.report.pdf")

if (5 %in% stages){
  cat(paste("[",Sys.time() ,"] Stage 5 - Compute clonality of each valid segment\n",sep=""))
  source(paste0(path_to_CLONET_functions,"Functions/CLONET.ComputeClonalityTable.R"))
}else{
  cat(paste("[",Sys.time() ,"] Skip stage 5 - Compute clonality of each valid segment\n",sep=""))
}
  

#####################################
## Allele specific
cat("\n\n######################################################################\n")

AITable_file <- paste0(results_dir,"allelicImbalanceTable.txt")
AIreportFile <- paste0(results_dir,"allelicImbalance.report.pdf")

if (6 %in% stages){
  cat(paste("[",Sys.time() ,"] Stage 6 - Allele specific analysis\n",sep=""))
  source(paste0(path_to_CLONET_functions,"Functions/CLONET.ComputeAlleleSpecificTable.R"))
}else{
  cat(paste("[",Sys.time() ,"] Skip stage 6 - Allele specific analysis\n",sep=""))
}


cat("\n")
cat(paste0("[",Sys.time() ,"] Computation ended\n"))
