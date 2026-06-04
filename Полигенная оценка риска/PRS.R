setwd("/mnt/data/chirkovayd/PRSise/")
library(data.table)
library(fmsb)
library(pROC)
library(ggplot2)

data <- fread('penncath.csv')
cat("Распределение CAD:\n"); print(table(data$CAD))
cat("\nСоотношение случаев/контролей:", 
    round(table(data$CAD)[2]/table(data$CAD)[1], 2))
cat("\n\nПропущенные значения:\n")
print(colSums(is.na(data)))

# проверка дубликатов
bim <- fread('imputed.QC.bim')
colnames(bim) <- c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2')
snp_counts <- bim[, .N, by = SNP]
duplicated_snps <- snp_counts[N > 1, unique(SNP)]
unique_snps <- snp_counts[N == 1, .(SNP)]
write.table(unique_snps,'imputed.QC.unique.snps',
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE, sep = '\t')

length(duplicated_snps)
nrow(unique_snps)

# проверка гетерозиготности
het <- fread('imputed.QC.het')
het$F <- (het$'O(HOM)'- het$'E(HOM)') / (het$'N(NM)'- het$'E(HOM)')
mean_F <- mean(het$F, na.rm = TRUE)
sd_F <- sd(het$F, na.rm = TRUE)
valid_het <- het[F <= (mean_F + 3*sd_F) & F >= (mean_F- 3*sd_F)]
write.table(valid_het[, .(FID, IID)],'imputed.QC.valid.het',
              col.names = FALSE, row.names = FALSE, quote = FALSE,
              sep = '\t')
valid_het_list <- fread('imputed.QC.valid.het', header = FALSE)
colnames(valid_het_list) <- c('FID', 'IID')

cat("Среднее F:", mean_F, "\n")
cat("SD F:", sd_F, "\n")
cat("Аномальных образцов:", sum(het$F < (mean(het$F) - 3*sd(het$F)) | 
                                  het$F > (mean(het$F) + 3*sd(het$F))))

# удаление родственников
rel_samples <- fread('imputed.QC.rel.id', header = FALSE)
excluded_rel_samples <- fread('imputed.QC.rel.excluded', header = FALSE)
colnames(rel_samples) <- c('FID', 'IID')
final_valid <- merge(valid_het_list, rel_samples, by = c('FID', 'IID'))
write.table(final_valid,'imputed.QC.final',
              col.names = FALSE, row.names = FALSE, quote = FALSE,
              sep = '\t')

nrow(final_valid)
bim <- fread('imputed.QC.filtered.bim', header = FALSE)

# Фильтрация base данных
colnames(bim) <- c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2')
bim[, A1 := toupper(A1)]
bim[, A2 := toupper(A2)]

base <- fread('CAD_META')
dup_counts <- base[, .N, by = oldID]
cat("Дубликатов rsID:", sum(dup_counts$N > 1), "\n")

base <- base[!duplicated(oldID)]
base[, Allele1 := toupper(Allele1)]
base[, Allele2 := toupper(Allele2)]
base <- base[!is.na(oldID) & oldID != ""]
base <- base[!is.na(Freq1) & Freq1 >= 0.01 & Freq1 <= 0.99]
cat("SNP после фильтрации:", nrow(base), "\n")
fwrite(base, 'CAD_META.filtered.txt', sep = '\t', quote = FALSE)

# Согласование кодировки аллелей
complement <- function(x) {
  switch(x, "A" = "T", "C" = "G", "T" = "A", "G" = "C", return(NA))
}
View(base)
View(bim)
info <- merge(bim, base, by.x = c("SNP", "CHR", "BP"),
                by.y = c("oldID", "CHR", "BP"), all.x = FALSE)
View(info)
info_match <- info[A1 == Allele1 & A2 == Allele2]
info$C_A1 <- sapply(info$A1, complement)
info$C_A2 <- sapply(info$A2, complement)
info_complement <- info[C_A1 == Allele1 & C_A2 == Allele2]
info_recode <- info[A1 == Allele2 & A2 == Allele1]
info_crecode <- info[C_A1 == Allele2 & C_A2 == Allele1]

compatible_snps <- unique(c(info_match$SNP, info_complement$SNP,
                               info_recode$SNP, info_crecode$SNP))
mismatch_snps <- bim[!SNP %in% compatible_snps, SNP]

if (length(mismatch_snps) > 0) {
  write.table(data.table(SNP = mismatch_snps),
              "imputed.mismatch",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}
cat("Совпадающие аллели имеют ", nrow(info_match), "SNP\n")
cat("Strand flipping требуется для ", nrow(info_complement), "SNP\n")
cat("Recoding требуется для ", nrow(info_recode), "SNP\n")
cat("Совсем не совпадающие аллели имеют ", length(mismatch_snps), "SNP\n")

# разделение на выборки
pheno <- fread('penncath.csv')
cases <- pheno[CAD == 1]
controls <- pheno[CAD == 0]
set.seed(123)
train_cases_idx <- sample(1:nrow(cases), size = round(0.7 * nrow(cases)))
train_cases <- cases[train_cases_idx]
val_cases <- cases[-train_cases_idx]

train_controls_idx <- sample(1:nrow(controls), size = round(0.7 * nrow(controls)))
train_controls <- controls[train_controls_idx]
val_controls <- controls[-train_controls_idx]

train_samples <- rbind(train_cases, train_controls)
val_samples <- rbind(val_cases, val_controls)
write.table(train_samples[, .(FamID, FamID)], 'train_samples.txt',
            col.names = F, row.names = F, quote = F, sep = '\t')
write.table(val_samples[, .(FamID, FamID)], 'val_samples.txt',
              col.names = F, row.names = F, quote = F, sep = '\t')

cat("В тренировочной выборке", nrow(fread('imputed.train.bim', header = FALSE)), "SNP и",
    nrow(fread('imputed.train.fam', header = FALSE)), "образцов.\n")
cat("В валидационной выборке", nrow(fread('imputed.val.bim', header = FALSE)), "SNP и",
    nrow(fread('imputed.val.fam', header = FALSE)), "образцов.\n")

# создание файла ковариат
pca <- fread("imputed.QC.eigenvec")
View(pca)
eigenvals <- fread("imputed.QC.eigenval")

plot(1:20, eigenvals$V1[1:20], 
     type = "b", 
     xlab = "Главная компонента", 
     ylab = "Собственное значение",
     main = "Scree plot: определение числа PC")
abline(v = which(diff(eigenvals$V1[1:20]) < 0.1*max(eigenvals$V1))[1] + 1, 
       col = "red", lty = 2)

covariate_file <- pheno[, .(FamID, FamID, sex, age)]
colnames(covariate_file) <- c("FID", "IID", "Sex", "Age")
pcs <- fread("imputed.QC.eigenvec", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:(ncol(pcs)-2)))
pcs <- pcs[, 1:9]
covariate_file <- merge(covariate_file, pcs, by = c("FID", "IID"),
                          all.x = TRUE)
write.table(covariate_file, "covariate_file.txt",
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep
              = "\t")
View(covariate_file)


# Анализ резултатов PRS
prs_results <- fread("train.PRS.prsice")
View(prs_results)
best_model <- prs_results[which.max(R2)]
cat("Best-fit threshold:", best_model$Threshold, "\n")

cat("Best-fit R²:", best_model$R2, "\n")

cat("SNP used in best model:", best_model$Num_SNP, "\n")

print(prs_results[, .(Threshold, R2, Num_SNP)])
best_model

# создание файла с весами SNP
best_file <- fread('train.PRS.summary', header = TRUE)
p_threshold <- best_file[1, Threshold]

prs_snps <- fread('train.PRS.snp', header = TRUE)
snps_below_threshold <- prs_snps[prs_snps$P < p_threshold, ]
rsids_to_extract <- snps_below_threshold$SNP

base_filtered <- fread('CAD_META.filtered.txt')
output_df <- base_filtered[base_filtered$oldID %in% rsids_to_extract,
                           c("oldID", "Allele1", "Effect")]
colnames(output_df) <- c("rsid", "effect_allele", "weight")
output_df <- output_df[!duplicated(rsid)]
write.table(output_df,'prs_CAD.txt', row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = '\t')

# оценка предсказательной способности на обучающей выборке
train_prs <- fread('train.PRS.best')
metadata <- fread('penncath.csv')
covariates <- fread('covariate_file.txt')

train_data <- merge(train_prs, covariates, by = "IID")
train_data <- merge(train_data, metadata, by.x = "IID", by.y = "FamID")

#train_data <- merge(train_prs, metadata, by.x = 'IID', by.y = 'FamID', all.x = TRUE)
model_train <- glm(CAD ~ PRS + sex + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, 
                   data = train_data, family = binomial)
summary_model <- summary(model_train)
summary(model_train)
prs_pvalue <- coef(summary_model)["PRS", "Pr(>|z|)"]
cat("P-value для PRS:", format(prs_pvalue, scientific = TRUE, digits = 3), "\n")
if(prs_pvalue < 0.05) {
  cat("PRS является значимым предиктором CAD (p < 0.05)\n")
} else {
  cat("PRS НЕ является значимым предиктором CAD (p >= 0.05)\n")
}

prs_coef <- coef(model_train)["PRS"]
prs_or <- exp(prs_coef)
prs_ci <- exp(confint(model_train)["PRS", ])
cat("Коэффициент β (log OR):", round(prs_coef, 4), "\n")
cat("Odds Ratio (OR):", round(prs_or, 4), "\n")
cat("95% доверительный интервал OR: [", round(prs_ci[1], 4), ",", round(prs_ci[2], 4), "]\n")

NagelkerkeR2(model_train)$R2
detach("package:fmsb", unload = TRUE)
roc_train <- roc(train_data$CAD, train_data$PRS, quiet = TRUE)
auc_train <- auc(roc_train)
auc_train
plot(roc_train, main = 'ROC Curve- Training Set')

# оценка на валидационной выборке
val_prs <- fread('val.PRS.profile')
colnames(val_prs)[colnames(val_prs) == "SCORE"] <- "PRS"

val_data <- merge(val_prs, covariates, by = "IID")
val_data <- merge(val_data, metadata, by.x = "IID", by.y = "FamID")

#val_data <- merge(val_prs, metadata, by.x = 'IID', by.y = 'FamID', all.x = TRUE)

model_val <- glm(CAD ~ PRS + sex + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, 
                 data = val_data, family = binomial)
summery_val_model <- summary(model_val)
summary(model_val)

prs_pvalue_val <- coef(summery_val_model)["PRS", "Pr(>|z|)"]
cat("P-value для PRS:", format(prs_pvalue_val, scientific = TRUE, digits = 3), "\n")
if(prs_pvalue_val < 0.05) {
  cat("PRS является значимым предиктором CAD (p < 0.05)\n")
} else {
  cat("PRS НЕ является значимым предиктором CAD (p >= 0.05)\n")
}

prs_coef_val <- coef(model_val)["PRS"]
prs_or_val <- exp(prs_coef_val)
prs_ci_val <- exp(confint(model_val)["PRS", ])
cat("Коэффициент β (log OR):", round(prs_coef_val, 4), "\n")
cat("Odds Ratio (OR):", round(prs_or_val, 4), "\n")
cat("95% доверительный интервал OR: [", round(prs_ci_val[1], 4), ",", round(prs_ci_val[2], 4), "]\n")

library(fmsb)
NagelkerkeR2(model_val)$R2
detach("package:fmsb", unload = TRUE)
roc_val <- roc(val_data$CAD, val_data$PRS, quiet = TRUE)
auc_val <- auc(roc_val)
auc_val
plot(roc_val, main = 'ROC Curve- Validation Set')

# распределения PRS
val_data$CAD_factor <- factor(val_data$CAD, levels = c(0, 1),
                                labels = c("Controls", "Cases"))
ggplot(val_data, aes(x = PRS, fill = CAD_factor)) + geom_density(alpha = 0.6) +
  labs(x = 'Polygenic Risk Score (PRS)', y = 'Density',
     title = 'PRS Distribution- Validation Set', fill = 'Group') +
  theme_minimal()
