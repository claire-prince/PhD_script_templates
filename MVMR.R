cran_packages <- c("readr", "remotes", "rlang", "tibble", "vroom", "data.table", "tidyverse")
cran_installed_packages <- cran_packages %in% rownames(installed.packages())
if (any(cran_installed_packages == FALSE)) {
  install.packages(cran_packages, repos = "https://www.stats.bris.ac.uk/R/")
}
lapply(cran_packages, library, character.only = TRUE)
github_packages <- c("MRCIEU/TwoSampleMR", "mrcieu/ieugwasr")
remotes::install_github(github_packages)
devtools::install_github("explodecomputer/genetics.binaRies", "WSpiller/MVMR")
library(TwoSampleMR)
library(ieugwasr)
library(MVMR)
library(genetics.binaRies)

# Source function script from https://marinalearning.netlify.app/2021/03/22/setting-up-multivariable-mendelian-randomization-analysis/
source("~/function.R")

# Read in exposure, adjustment and outcome variable intended for analysis
ref=as.matrix(read.csv("~/analysis.csv"))

# Set where data is 
data_path="~/data/"

# Set where phenotypic correlation data is 
pheno_path="~/pheno_correlations/"

# Set where to save output
output_path="~/outputs/"

for (row in 1:nrow(ref)) {
    print(row)
    # Pull out variables from row of reference for analysis
    exp  <- ref[row, "exp"]
    adj <- ref[row, "adj"]
    out <- ref[row, "out"]
    exp_lab <- ref[row,"Exposure"]
    adj_lab <- ref[row,"Adjustment"]
    out_lab <- ref[row, "Outcome"]
    file_name <- ref[row, "file_name"]

    print(exp_lab)
    print(adj_lab)
    print(out_lab)

    # Read in top hits for exposure and adjustment
    exp_tophits <- read_delim(paste0(data_path, exp, "_tophits.tsv"), delim=",")
    adj_tophits <- read_delim(paste0(data_path, adj, "_tophits.tsv"), delim=",")

    # Read in full gwas for exposure and adjustment
    file=paste0(data_path, exp, "_gwas.txt.gz")
    exp_gwas <- as_tibble(fread(file, sep=","))
    file=paste0(data_path, adj,"_gwas.txt.gz")
    adj_gwas <- as_tibble(fread(file, sep=","))

    # Combine tophits in list
    tophits_list <- list(exp_tophits, adj_tophits)

    # Combine full gwas' in list
    full_gwas_list <- list(exp_gwas, adj_gwas)

    # Create exposure_dat, i.e. obtain instruments for each exposure
    exposure_dat <- get_mv_exposures(tophits_list, full_gwas_list)

    # Match up id.exposure and exposurein exposure_dat
    exposure_dat$id.exposure=exposure_dat$exposure

    # Read in outcome gwas based on SNPs in exposure_dat
    outcome_dat <- read_outcome_data(
        filename=paste0(data_path,file_name, ".txt.gz"),
        snps = exposure_dat$SNP)
    

    # Once the data has been obtained, harmonise so that exposure_dat and outcome_dat are on the same reference allele
    rawdat_mvmr <- mv_harmonise_data(exposure_dat, outcome_dat)

    #ensure the exposure of interest is first in the data
    rawdat_mvmr$exposure_beta <- rawdat_mvmr$exposure_beta[, c(which(colnames(rawdat_mvmr$exposure_beta)==exp_lab), 
                                                                which(colnames(rawdat_mvmr$exposure_beta)==adj_lab))] 
    rawdat_mvmr$exposure_se <- rawdat_mvmr$exposure_se[, c(which(colnames(rawdat_mvmr$exposure_se)==exp_lab), 
                                                                which(colnames(rawdat_mvmr$exposure_se)==adj_lab))] 
    rawdat_mvmr$exposure_pval <- rawdat_mvmr$exposure_pval[, c(which(colnames(rawdat_mvmr$exposure_pval)==exp_lab), 
                                                                which(colnames(rawdat_mvmr$exposure_pval)==adj_lab))] 


    # Format data for MVMR
    ## rawdat_mvmr should have the format exp1_betas, exp2_betas, exp1_se, exp2_se, out_betas, , out_se, SNP
    F.data <- format_mvmr(BXGs = rawdat_mvmr$exposure_beta,
                        BYG = rawdat_mvmr$outcome_beta,
                        seBXGs = rawdat_mvmr$exposure_se,
                        seBYG = rawdat_mvmr$outcome_se,
                        RSID = rownames(rawdat_mvmr$exposure_beta))
    rownames(F.data)=NULL

    # Read in phenotypic correlation
    mvmrcovmatrix=as.matrix(read.table(paste0(pheno_path, exp, "_", adj, "_corr.txt"), sep=",", header=F))

    # Create correlation dataframe 
    Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,6:7])

    # Test instrument strength
    sres <- strength_mvmr(r_input = F.data, gencov = Xcovmat)  

    # test for heterogeneity
    pres <- data.frame(pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat))


    # Estimate causal effect i.e. run MVMR using the IVW method
    mvmr <- ivw_mvmr(r_input = F.data)

    # Run qhet if instrument strength is less than 10, or evidence for pleiotropy is less than pval=0.05
    ## Qhet does not performed if instrument strength is less than 5, in these cases qhet will not be run
    if (sres[1]>10 & pres["Qpval"]>0.05) {
        print("Strong instruments and no evidence of heterogeneity")
        next
        } 
    if (sres[1]<4) {
        print("Instruments too weak to run Qhet")
        next 
        }   
    if ((sres[1]<10 | sres[1]>=4) | pres["Qpval"]<0.05) {
        print("Run qhet")
        qhet <- qhet_mvmr(F.data, mvmrcovmatrix, CI = T, iterations = 1000)
        qhet_vector <- c(exp_lab, adj_lab, out_lab, qhet[1,])
        qhet_df=rbind(qhet_df, qhet_vector)
        colnames(qhet_df) <- c("Exposure", "Adjustment", "Outcome", c(colnames(qhet)))
        } 

    # format data for MR Egger analysis
    MRMVInputObject <- mr_mvinput(bx = rawdat_mvmr$exposure_beta,
                                bxse = rawdat_mvmr$exposure_se,
                                by = rawdat_mvmr$outcome_beta, 
                                byse = rawdat_mvmr$outcome_se,
	 	        snps = rownames(rawdat_mvmr$exposure_beta)) 
  
  
    MRMVObject <- mr_mvegger(MRMVInputObject, orientate = 1)
    egger_results <- t(rbind(MRMVObject@Exposure, MRMVObject@Estimate,  MRMVObject@StdError.Est, MRMVObject@CILower.Est,MRMVObject@CIUpper.Est, MRMVObject@Pvalue.Est))
    colnames(egger_results) <- c("Exposure", "Estimate", "SE", "CIL", "CIU","Pvalue")

}