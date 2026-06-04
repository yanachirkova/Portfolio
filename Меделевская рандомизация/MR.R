library(TwoSampleMR)
library(data.table)
setwd('/mnt/data/chirkovayd/MR')

exp_raw <- fread('intraocular.tsv')
View(exp_raw)

# select only significant variants
exp_raw <- subset(exp_raw,exp_raw$p_value<5e-8)
exp_dat <- format_data(as.data.frame(exp_raw),
                        type = "exposure",
                        snp_col = "rsid",
                        beta_col = "beta",
                        se_col = "standard_error",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        eaf_col = "effect_allele_frequency",
                        pval_col = "p_value"
)

# clumping
clumped_exp <- clump_data(exp_dat,clump_r2=0.01,pop="EAS")

# outcome
out_raw <- fread("glaucoma.tsv")
names(out_raw)
head(out_raw, 5)

out_raw$beta <- log(out_raw$odds_ratio)
out_raw$SE <- (log(out_raw$ci_upper) - log(out_raw$ci_lower)) / (2 * 1.96)

out_dat <- format_data(as.data.frame(out_raw),
                       type = "outcome",
                       snp_col = "rsid",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele",
                       eaf_col = "effect_allele_frequency",
                       pval_col = "p_value",
)

harmonized_data <- harmonise_data(clumped_exp,out_dat,action=1)

# выполняем MR
res <- mr(harmonized_data)
p1 <- mr_scatter_plot(res, harmonized_data)
p1[[1]]

# тест на гетерогенность
heterogenity <- mr_heterogeneity(harmonized_data)
heterogenity

# Intercept in MR-Egger
mr_pleiotropy_test(harmonized_data)

# single SNP MR
res_single <- mr_singlesnp(harmonized_data)
res_single
p2 <- mr_forest_plot(res_single)
p2[[1]]
p4 <- mr_funnel_plot(res_single)
p4[[1]]

# leave-one-out MR
res_loo <- mr_leaveoneout(harmonized_data)
res_loo
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

