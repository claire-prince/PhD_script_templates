cran_packages <- c("remotes", "rlang", "tidyverse", "devtools")
cran_installed_packages <- cran_packages %in% rownames(installed.packages())
if (any(cran_installed_packages == FALSE)) {
  install.packages(cran_packages)
}
lapply(cran_packages, library, character.only = TRUE)

devtools::install_github("rondolab/MR-PRESSO")
remotes::install_github("MRCIEU/TwoSampleMR")
lapply(c("TwoSampleMR", "MRPRESSO"), library, character.only = TRUE)

# Read in exposure and outcome variable intended for analysis
ref=as.matrix(read.csv("~/analysis.csv"))

# Set where data is 
data_path="~/data"

# Set where phenotypic correlation data is 
pheno_path="~/pheno_correlations"

# Set where to save output
output_path="~/outputs"

for (row in 1:nrow(df)) { # nolint: seq_linter.
    # Pull out variables from row for analysis
    exp  <- df[row, "exp"]
    out <- df[row, "out"]
    exp_lab <- df[row,"Exposure"]
    out_lab <- df[row, "Outcome"]

    file_name <- df[row, "file_name"]

    print(exp_lab)
    print(out_lab)

    # Read in tophits for exposure
    exp_tophits <- read_delim(paste0(data_path, exp, "_tophits.tsv"), delim = ",") # nolint: line_length_linter.

    # Read in outcome gwas based on SNPs in exp_tophits
    outcome_dat <- read_outcome_data(
        filename=paste0(data_path, file_name, "_gwas.txt.gz"),
        snps = exp_tophits$SNP
    )

    dat <- harmonise_data(exp_tophits, outcome_dat)
    res <- mr(dat)

    results <- generate_odds_ratios(res)

    #generating F statistic
    f_stat=mean((dat$beta.exposure^2)/(dat$se.exposure^2))
    if (nrow(dat)<=1){next}

    # perform and output heterogeneity test from TwoSampleMR
    het=mr_heterogeneity(dat)

    # perform and output egger intercept test pleiotropy from TwoSampleMR
    p_test=mr_pleiotropy_test(dat)

    #perform and output mr presso from MR-PRESSO
    presso=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)
    presso_res=presso$`Main MR results`

    presso_gtest=presso$`MR-PRESSO results`$`Global Test`
    presso_dtest=presso$`MR-PRESSO results`$`Distortion Test`
    ##create data frame of outlier and distorition information
    ### if no outliers and therefore not distortion results
    if (is.null(presso_dtest)) {
        presso_outliers <- cbind(presso_gtest$`RSSobs`,presso_gtest$`Pvalue`)
        colnames(presso_outliers) <- c("RSSobs", "Pvalue")
    ### else include distortion results
    } else {
        presso_outliers <- cbind(presso_gtest$`RSSobs`,presso_gtest$`Pvalue`, presso_dtest$`Distortion Coefficient`, presso_dtest$`Pvalue`, length(presso_dtest$`Outliers Indices`))
        colnames(presso_outliers) <- c("RSSobs", "Pvalue", "Distortion Coefficient", "Pvalue", "Number of outliers")
    }

}   