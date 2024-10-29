df_edu = read.table(snakemake@input[["GWAS_edu"]], header = T)
snp_df = read.csv(snakemake@input[["snp_map"]], sep = '\t')
df_edu$Predictor <- sapply(df_edu$variant, function(x) {paste(strsplit(as.character(x), ":")[[1]][1:2], collapse = ":")})
df_edu = merge(df_edu, snp_df, by.x = "Predictor", by.y = "id")
df_edu["Direction"] = sign(df_edu$beta)
df_edu["Stat"] = (df_edu$beta/df_edu$se)**2
df_edu$Z = df_edu$beta/df_edu$se

df_edu = df_edu[, c("rsid.x", "alt", "ref", "Z", "n_complete_samples", "pval")]
names(df_edu) = c("SNP", "A1", "A2", "Z", "N", "p")

write.table(df_edu, snakemake@output[[1]], sep = '\t', row.names = F, quote = F)
