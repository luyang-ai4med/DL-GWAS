df = read.table(snakemake@input[["GWAS"]], header = T)
snp_df = read.csv(snakemake@input[["snp_map"]], sep = '\t')

df = merge(df, snp_df, by.x = "Predictor", by.y = "id")

df$Z = sqrt(df$Stat)*df$Direction
df$p = 2*pnorm(abs(df$Z), lower.tail = F)

df = df[, c("rsid", "A1", "A2", "Z", "n", "p")]
names(df) = c("SNP", "A1", "A2", "Z", "N", "p")

write.table(df, snakemake@output[[1]], sep = '\t', row.names = F, quote = F)
