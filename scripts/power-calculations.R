# ------------------------------------------------------------------
# Power calcualtions
# ------------------------------------------------------------------

## Aim: to calculate power to run one-sample and two-sample MR of vitb12 levels on DNAm

## outputs:
## Power to do 1SMR of MATERNAL vitb12 on NEONATAL dnam
## Power to do 2SMR of NEONATAL vitb12 on NEONATAL dnam

## Date: 2021-06-10

## pkgs
library(tidyverse) # tidy code and data
library(TwoSampleMR) # Helpful for functions?
library(readxl) # reading in excel spreadsheets - gwas data!

## data
vitb12_gwas <- read_xlsx("data/vitb12-gwas-data.xlsx")
vitb12_ewas_mat <- read_csv("data/MA.b12.maternal.main.model.GSM.no.MARBLES.csv")
vitb12_ewas_neo <- read_csv("data/MA.b12.newborn.main.model.GSM.csv")
vitb12_dist_mat <- read_xlsx("data/b12-distribution-edit.xlsx", sheet = "maternal")
load("data/aries-antenatal-variances.RData")
cpg_var <- data.frame(cpg = names(vars), varcpg = vars)

ewas_res_mat <- vitb12_ewas_mat %>% 
	dplyr::filter(P.value.FDR < 0.05) %>% 
	left_join(cpg_var, by = c("MarkerName" = "cpg"))

ewas_res_neo <- vitb12_ewas_neo %>%
	dplyr::filter(P.value.FDR < 0.05)

# ------------------------------------------------------------------
# Formatting exposure data
# ------------------------------------------------------------------

get_se <- function(p, beta) abs(beta / qnorm(p/2))

vitb12_gwas <- vitb12_gwas %>%
	mutate(original_se = get_se(pval, original_beta), 
		   se = get_se(pval, beta))

exposure <- format_data(vitb12_gwas)

# ------------------------------------------------------------------
# Functions for 2SMR power calculations
# ------------------------------------------------------------------

## Calculate r2
calc_r2 <- function(b, eaf) 2 * b^2 * eaf * (1 - eaf)

## Calculate F-statistic
calc_f <- function(r2, n, k) r2 * (n - 1 - k) / ((1 - r2) * k) # k = number of snps

## Calculate power
calc_power_2smr <- function(n, r2, beta, alpha) pnorm(sqrt(n * r2) * beta - qnorm(1 - alpha / 2))

# ------------------------------------------------------------------
# 2SMR
# ------------------------------------------------------------------

exposure_clumped <- clump_data(exposure, clump_r2 = 0.001, clump_kb = 10000)
nrow(exposure_clumped) # 11
nrow(exposure) # 11 -- no SNPs removed

## getting exposure info for easy use
eaf <- exposure_clumped[, "eaf.exposure"]
b <- exposure_clumped[, "beta.exposure"]
se <- exposure_clumped[, "se.exposure"]
p <- exposure_clumped[, "pval.exposure"]
snp <- exposure_clumped[, "SNP"]
vitb12_gwas_n <- 45576 # n comes from vitb12 gwas paper (PMID = 23754956)

## individual r2 and F-stats
r2 <- calc_r2(b, eaf)
f_stat <- calc_f(r2, n = vitb12_gwas_n, k = 1) 

## power in the 2smr when using each exposure SNP individually
n_out <- 25000 # roughly (from godmc website)
power <- lapply(ewas_res_neo$Effect, function(x) {
	tibble(snp = snp, power = calc_power_2smr(n_out, r2, beta = x, alpha = 0.05))
})
names(power) <- ewas_res_neo$MarkerName

power_ind <- bind_rows(power, .id = "cpg")

arrange(power_ind, desc(power))

## comb r2 and F-stats
all_r2 <- sum(r2)
all_f_stat <- calc_f(all_r2, n = vitb12_gwas_n, k = length(snp)) 

## power in the 2smr when combining the exposure SNPs (e.g. in meta-analysis or as a GRS)
all_power <- map_dfr(1:nrow(ewas_res_neo), function(x) {
	tibble(cpg = ewas_res_neo$MarkerName[x], 
		   power = calc_power_2smr(n_out, all_r2, beta = ewas_res_neo$Effect[x], alpha = 0.05))
})

## output results
out_2smr <- list(power_comb = all_power, power_ind = power_ind, 
				 r2 = all_r2, f_stat = all_f_stat)

save(out_2smr, file = "results/power-2smr.RData")

# ------------------------------------------------------------------
# Functions for 1SMR power calculations
# ------------------------------------------------------------------

## Function taken from https://github.com/kn3in/mRnd/blob/master/functions.R (which feeds into https://shiny.cnsgenomics.com/mRnd/)
calc_power_1smr <- function(N, alpha, byx, bOLS, R2xz, varx, vary, epower) {
    
    threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
    f.value <- 1 + N * R2xz / (1 - R2xz)
    con <- (bOLS - byx) * varx # covariance due to YX confounding
    vey <- vary - byx * varx * (2 * bOLS - byx)
    
    if (vey < 0) {
    
        data.frame(Error = "Error: Invalid input. The provided parameters result in a negative estimate for variance of the error term in the two-stage least squares model.")
    
    } else {

        if (is.na(epower)) {
        
            b2sls <- byx + con / (N * R2xz)
            v2sls <- vey / (N * R2xz * varx)
            NCP <- b2sls^2 / v2sls
            # 2-sided test
            power <- 1 - pchisq(threschi, 1, NCP)
            data.frame(power = power, f_stat = f.value)
        
        } else {
        
            # Calculation of sample size given power
            z1 <- qnorm(1 - alpha / 2)
            z2 <- qnorm(epower)
            Z  <- (z1 + z2)^2
            # Solve quadratic equation in N
            a <- (byx * R2xz)^2
            b <- R2xz * (2 * byx * con - Z * vey / varx)
            c <- con^2
            N1 <- ceiling((-b + sqrt(b^2 - 4 * a * c)) / (2 * a))
            data.frame(n_required = N1)
        
        }
    }
}

# ------------------------------------------------------------------
# 1SMR power calculations
# ------------------------------------------------------------------

n <- seq(1000, 3000, 500) # N (taken from analysis plan)
alpha <- 0.05 # alpha level
byx <- ewas_res_mat$Effect # unknown true causal effect (taken as the EWAS estimate)
bOLS <- ewas_res_mat$Effect * 1.5 # effect including that of confounders (taken as EWAS estimate)
R2xz <- sum(calc_r2(b, eaf)) # variation of exposure explained by SNPs
varx <- 1 # variance of exposure - this is one as the vitb12 levels were scaled
vary <- ewas_res_mat$varcpg # variance of outcome. This is taken from ALSPAC samples
epower <- NA # power level (estimating this so set to NA)

power_res <- lapply(1:nrow(ewas_res_mat), function(x) {
	x <- lapply(n, function(N) calc_power_1smr(N, alpha, byx[x], bOLS[x], R2xz, varx, vary[x], epower))
	names(x) <- n
	out <- bind_rows(x, .id="n")
	return(out)
})
names(power_res) <- ewas_res_mat$MarkerName
power_1smr <- bind_rows(power_res, .id="cpg")

power_1smr %>% 
	arrange(desc(power)) %>% 
	head

ori <- power_1smr

## sample sizes needed to run analyses
epower <- 0.8 # 80% power
sample_res <- map_dfr(1:nrow(ewas_res_mat), function(x) {
	tibble(cpg = ewas_res_mat$MarkerName[x], 
		   calc_power_1smr(N = NULL, alpha, byx[x], bOLS[x], R2xz, varx, vary[x], epower))
})

## output results 
out_1smr <- list(power_res = power_1smr, sample_res = sample_res)

save(out_1smr, file = "results/power-1smr.RData")

