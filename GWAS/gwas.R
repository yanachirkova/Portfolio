library(data.table)
library(ggplot2)
library(dplyr)
library(bigsnpr)
library(bigstatsr)


head(fam <- fread('gwas/penncath.fam'))

head(bim <- fread('gwas/penncath.bim'))

df <- fread('gwas/penncath.csv')
sum(is.na(df))
df %>% 
  count(CAD)

snp_readBed("gwas/penncath.bed")
penncath <- snp_attach("gwas/penncath.rds")
str(penncath)
penncath$map$chromosome |> table()
unique(penncath$map$chromosome)

# Доля пропущенных значений по SNP
counts <- big_counts(penncath$genotypes)
missing_rate <- counts[4, ] / colSums(counts)
sum(missing_rate > 0.50)
sum(missing_rate > 0.10)

print(t(counts[, 1:10]))
View(t(counts))

# Доля гетерозиготности
sample_counts <- big_counts(penncath$genotypes, byrow = TRUE)
geno_matrix <- sample_counts[1:3, ]
het_matrix <- sweep(geno_matrix, MARGIN = 1, STATS = c(0, 1, 0), FUN = "*")
het_counts <- colSums(het_matrix)
total_snps_per_sample <- colSums(geno_matrix)
het_rate <- het_counts / total_snps_per_sample
sum(het_rate < 0.25)
sum(het_rate > 0.30)

# Подсчет MAF
calculate_maf <- function(G) {
  # 1. Получаем счетчики генотипов
  cnt <- big_counts(G)  # матрица 4 x n_SNP
  
  # 2. Всего успешных генотипов (без NA)
  total <- colSums(cnt[1:3, ])
  
  # 3. Частота альтернативного аллеля
  #    cnt[1,] = гомозиготы ref (0)
  #    cnt[2,] = гетерозиготы (1)  
  #    cnt[3,] = гомозиготы alt (2)
  freq_alt <- (cnt[2, ] + 2 * cnt[3, ]) / (2 * total)
  
  # 4. MAF = min(freq_alt, 1 - freq_alt)
  pmin(freq_alt, 1 - freq_alt)
}
maf_values <- calculate_maf(penncath$genotypes)

summary(maf_values)
cat("Всего SNP:", length(maf_values), "\n")
cat("MAF < 0.01:", sum(maf_values < 0.01, na.rm = TRUE), "\n")
cat("MAF < 0.05:", sum(maf_values < 0.05, na.rm = TRUE), "\n")
maf_001 <- sum(maf_values < 0.01, na.rm = TRUE)
maf_005 <- sum(maf_values < 0.05, na.rm = TRUE)
cat("Доля SNP с MAF < 0.01:", round(maf_001/length(maf_values)*100, 2), "%\n")
cat("Доля SNP с MAF < 0.05:", round(maf_005/length(maf_values)*100, 2), "%\n")

# Проверка равновесия Харди-Вайнберга
hwe_simple <- function(G) {
  cnt <- big_counts(G)
  n <- colSums(cnt[1:3, ])
  p <- (cnt[2, ] + 2 * cnt[3, ]) / (2 * n)
  obs_het <- cnt[2, ] / n
  exp_het <- 2 * p * (1 - p)
  
  # Хи-квадрат тест
  chi2 <- n * (obs_het - exp_het)^2 / (exp_het * (1 - exp_het))
  
  # p-value
  pchisq(chi2, df = 1, lower.tail = FALSE)
}

hwe_pvals <- hwe_simple(penncath$genotypes)

hwe_deviations <- sum(hwe_pvals < 1e-6, na.rm = TRUE)
cat("SNP с отклонением от HWE (p < 1e-6):", hwe_deviations, "\n")
cat("Всего SNP:", length(hwe_pvals), "\n")
cat("Доля отклоняющихся SNP:", 
    round(hwe_deviations / length(hwe_pvals) * 100, 3), "%\n")

# Фильтрация
system("plink --bfile gwas/penncath --make-bed --out gwas/penncath_fixed")

qc_result <- snp_plinkQC(
  plink.path = "plink",
  prefix.in = "gwas/penncath_fixed",  # исправленный файл
  maf = 0.01,
  geno = 0.05,
  hwe = 1e-6,
  mind = 0.05,
  autosome.only = TRUE,
  prefix.out = "gwas/penncath_qc",
  verbose = TRUE
)

snp_readBed("gwas/penncath_qc.bed")
penncath_qc <- snp_attach("gwas/penncath_qc.rds")

cat("=== СРАВНЕНИЕ ДО И ПОСЛЕ QC ===\n")
cat("ДО фильтрации:\n")
cat("  Образцов:", nrow(penncath$genotypes), "\n")
cat("  SNP:", ncol(penncath$genotypes), "\n\n")

cat("ПОСЛЕ фильтрации:\n")
cat("  Образцов:", nrow(penncath_qc$genotypes), "\n")
cat("  SNP:", ncol(penncath_qc$genotypes), "\n\n")

# Импутация пропусков
qc_counts <- big_counts(penncath_qc$genotypes)
snp_with_na <- sum(qc_counts[4, ] > 0)
cat("SNP с пропущенными значениями:", snp_with_na, "\n")
cat("Всего SNP:", ncol(qc_counts), "\n")

prop_snp_with_na <- snp_with_na / ncol(qc_counts)
cat("Доля SNP с пропусками:", round(prop_snp_with_na * 100, 2), "%\n")

# Call rate по SNP
call_rate_snp <- 1 - (qc_counts[4, ] / nrow(penncath_qc$genotypes))

low_call_95 <- sum(call_rate_snp < 0.95)
low_call_99 <- sum(call_rate_snp < 0.99)

cat("SNP с call rate < 95%:", low_call_95, "\n")
cat("SNP с call rate < 99%:", low_call_99, "\n")

qc_sample_counts <- big_counts(penncath_qc$genotypes, byrow = TRUE)
samples_with_na <- sum(qc_sample_counts[4, ] > 0)
cat("Индивидуумов с пропущенными значениями:", samples_with_na, "\n")
cat("Всего индивидуумов:", ncol(qc_sample_counts), "\n")

# Call rate по индивидуумам
call_rate_sample <- 1 - (qc_sample_counts[4, ] / ncol(penncath_qc$genotypes))

cat("Образцы с call rate < 95%:", sum(call_rate_sample < 0.95), "\n")
cat("Образцы с call rate < 99%:", sum(call_rate_sample < 0.99), "\n")

imputed <- snp_fastImputeSimple(penncath_qc$genotypes, method = "mode")
#CODE_IMPUTE_PRED <- 3
cnt_before <- big_counts(penncath_qc$genotypes)
cnt_after <- big_counts(imputed)  # или imputed$genotypes

cat("NA до импутации:", sum(cnt_before[4, ]), "\n")
cat("NA после импутации:", sum(cnt_after[4, ]), "\n")

penncath_imputed <- list(
  genotypes = imputed,          # импутированная матрица
  fam = penncath_qc$fam,           # информация об образцах
  map = penncath_qc$map            # информация о SNP
)
class(penncath_imputed) <- "bigSNP"

snp_writeBed(penncath_imputed, 
             bedfile = "gwas/penncath_imputed.bed")
snp_readBed("gwas/penncath_imputed.bed")
penncath_imp <- snp_attach("gwas/penncath_imputed.rds")

# РСA
start_time <- Sys.time()

svd_result <- big_SVD(
  X = penncath_imp$genotypes,  # импутированные генотипы
  fun.scaling = big_scale(),   # центрирование и масштабирование
  k = 10                      # 10 главных компонент            
)

end_time <- Sys.time()
computation_time <- end_time - start_time
cat("SVD занимает ", computation_time, "\n")

str(svd_result)

saveRDS(svd_result, file = "gwas/svd_result.rds")

# scree plot
plot(svd_result, type = "screeplot")

# График PCA scores
plot(svd_result, type = "scores")

pca_df <- data.frame(
  FamID = penncath_imp$fam$family.ID,
  PC1 = svd_result$u[,1],
  PC2 = svd_result$u[,2]
)

pca_df <- left_join(
  pca_df,                     # левая таблица (все образцы)
  df,                   # правая таблица (клинические данные)
  by = "FamID"                   # ключ для соединения
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = sex, shape = as.factor(CAD))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA: пол и CAD статус",
       x = "PC1", y = "PC2") +
  scale_shape_discrete(name = "CAD", labels = c("Control", "Case")) +
  theme_minimal()

# Главные компоненты = u × diag(d)
PC <- sweep(svd_result$u, MARGIN = 2, STATS = svd_result$d, FUN = "*")

colnames(PC) <- paste0("PC", 1:ncol(PC))

for_gwas <- penncath_imp$fam %>%
  left_join(df, by = c("family.ID" = "FamID"))

# GWAS
gwas_result <- big_univLogReg(
  X = penncath_imp$genotypes,  # импутированные генотипы
  y01.train = for_gwas$CAD,              # фенотип (0/1)
  covar.train = PC[, 1:5],          # коварианты (первые 5 PC)
  ncores = nb_cores()           # использует все доступные ядра
)
saveRDS(gwas_result, file = "gwas/gwas_results.rds")
str(gwas_result)

# QQ-plot
snp_qq(gwas_result)

# Manhattan plot
manh_plot <- snp_manhattan(
  gwas = gwas_result,                     # результаты GWAS
  infos.chr = penncath_imp$map$chromosome, # хромосомы
  infos.pos = penncath_imp$map$physical.pos, # позиции
  colors = c("#1f78b4", "#a6cee3"),       # два цвета для четных/нечетных хромосом
  npoints = 20000,                        # ограничиваем точки для скорости
  coeff = 0.8                             # размер текста
)

manh_plot <- manh_plot +
  geom_hline(
    yintercept = -log10(5e-8),           # порог genome-wide значимости
    color = "red",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = -log10(5e-6),           # порог для предположительных результатов
    color = "orange",
    linetype = "dotted",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  labs(
    title = "Manhattan Plot GWAS",
    subtitle = paste("Всего SNP:", length(gwas_result$score)),
    x = "Хромосома",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50")
  )

print(manh_plot)

pvals <- predict(gwas_result, log10 = FALSE)

cat("=== РЕЗУЛЬТАТЫ ПО ПОРОГАМ ===\n")
cat("1. SNP с p < 5×10⁻⁸ (genome-wide значимость):", 
    sum(pvals < 5e-8, na.rm = TRUE), "\n")

cat("2. SNP с p < 5×10⁻⁶ (предположительные ассоциации):", 
    sum(pvals < 5e-6, na.rm = TRUE), "\n")

cat("3. Всего SNP в анализе:", length(pvals), "\n")
cat("   - Доля genome-wide значимых:", 
    round(sum(pvals < 5e-8)/length(pvals)*100, 4), "%\n")
cat("   - Доля предположительных:", 
    round(sum(pvals < 5e-6)/length(pvals)*100, 4), "%\n")

# Создание таблицы результатов GWAS
gwas_table <- data.frame(
  # 1. Идентификаторы SNP
  SNP_ID = penncath_imp$map$marker.ID,
  Chromosome = penncath_imp$map$chromosome,
  Position = penncath_imp$map$physical.pos,
  
  # 2. Статистики ассоциации
  Beta = gwas_result$estim,            # коэффициент β (log OR)
  SE = gwas_result$std.err,            # стандартная ошибка
  Score = gwas_result$score,           # χ² статистика
  P_value = predict(gwas_result, log10 = FALSE),  # p-значение
  
  # 3. Производные статистики
  logP = -log10(predict(gwas_result, log10 = FALSE)),  # -log10(p)
  OR = exp(gwas_result$estim),         # отношение шансов
  
  # 4. Дополнительная информация
  Allele1 = penncath_imp$map$allele1,   # референсный аллель
  Allele2 = penncath_imp$map$allele2,   # альтернативный аллель
  
  stringsAsFactors = FALSE
)

gwas_table <- gwas_table[order(gwas_table$P_value), ]

cat("Топ-10 наиболее значимых SNP:\n")
print(head(gwas_table, 10))

write.csv(gwas_table, file = "gwas/gwas_results_table.csv", row.names = FALSE)


cat("Всего SNP в таблице:", nrow(gwas_table), "\n")
cat("Genome-wide значимых (p < 5e-8):", sum(gwas_table$P_value < 5e-8), "\n")
cat("Предположительных (p < 5e-6):", sum(gwas_table$P_value < 5e-6), "\n")
cat("Самый значимый SNP:", gwas_table$SNP_ID[1], 
    "p =", format(gwas_table$P_value[1], scientific = TRUE), "\n")
cat("Максимальный OR:", round(max(gwas_table$OR, na.rm = TRUE), 2), "\n")
