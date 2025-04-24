## Path to CLONET
path_to_CLONET_functions <- "/home/exacloud/gscratch/mcweeney_lab/jengs/software/CLONET/"

## Output directory
output_DIR <- "/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/CLONETv2/"

# information about sample names
sampleInfoFile <- "/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/sample_info_file.csv"

# list of segment 
segmentListFile <- "/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/all_purecn_10000_maxclonal03.seg"

# folder with informative SNPs 
##NOT correct but RData are already present
pileup_dir <- "/home/exacloud/gscratch/mcweeney_lab/jengs/HNSCC/ASEQ_wdl/"

# Suffix of the informative SNPs pileup
PaPI_Suffix <- ".PILEUP.ASEQ"

## Path to error table
errorTable_file = "/home/exacloud/gscratch/mcweeney_lab/jengs/software/CLONET/Examples/CLONETv2/BasicData/errorTable.csv"

## Method to compete beta
# "GB" : classic model defined in http://www.genomebiology.com/2014/15/8/439 (only for backward compatibility)
# "STM" : new method based on error estimation on matched normal published in http://stm.sciencemag.org/content/6/254/254ra125.short (RECOMMENDED)
betaCompute.method <- "STM"

## Method to compete Adm.global
# "1D" : classic model that use only deletions (only for backward compatibility)
# "2D" : model based on cnA vs cnB transformation (RECOMMENDED)
adm.global.method <- "2D"

# minimum tumor coverage to consider an informative SNP as valid
minCoverage <- 10

# number of samples to process in parallel
NsamplesToProcessInParallel <- 2

# number of cores assigned to the analysis of each process
perSampleCores <- 2

# min number of informative SNPs for a genomic segment
min_nsnps <- 10
#min_nsnps <- 2
# min mean coverage of a genomic segment
min_cov <- 20
#min_cov <- 1

# minimum value of beta above which the two alleles are present in the same number (cn equal to 1+1 or 2+2 or 3+3 ...) to account ref map bias
equalCN.betaThr <- 0.90

# Homozygous deletions threshold (change only if you know what you are doing)
maxHomoDels <- 0.01

## Parameters of a valid deletion used to compute Adm.global
# log2 thresholds
deletionsLog2Levels <- c(-1,-0.25)

## Percentage of used deletions to compute Adm.global varibility interval (change only if you know what you are doing)
alphaPar <- 0.90

## Threshold on clonality value (change only if you know what you are doing)
clonalityThreshold <- 0.85

## Threshold on beta value (change only if you know what you are doing)
betaThreshold <- 0.85

## Stages to perform
# 1. analyse single sample and produe RData
# 2. create beta table aggragating single samples analysis
# 3. compute ploidy shift
# 4. compute global admixture
# 5. compute clonality table
# 6. compute allele specific copy number table
stages <- c(3,4)

## CLONET can try to compute reference mapping bias and to adjust beta estimation
computeRefMapBias <- F

## number of significant digits for the beta table (NA for reporting all digits)
n.digits <- 3

## maximum distance from (cnA, cnB) integer copy number value
## Dafault 0.5 corresponds to round cnA and cnB
AllelicImbalanceTh <- 0.5

