library(data.table)

# formulas from here: http://csg.sph.umich.edu/abecasis/publications/pdf/Bioinformatics.vol.26-pp.2190.pdf

df_1 = fread(snakemake@input[["quant"]])
df_2 = fread(snakemake@input[["binary"]])

df = merge(df_1, df_2, by = "Predictor", suffixes = c("_1", "_2"))
# harmonize to df_1
df = df[((df$A1_1 == df$A1_2) & (df$A2_1 == df$A2_2)) | ((df$A2_1 == df$A1_2) & (df$A1_1 == df$A2_2)), ]
df[((df$A2_1 == df$A1_2) & (df$A1_1 == df$A2_2)) == T, "Direction_2"] = -df[((df$A2_1 == df$A1_2) & (df$A1_1 == df$A2_2)) == T, "Direction_2"]

df$Z_1 = sqrt(df$Stat_1)*df$Direction_1
df$Z_2 = sqrt(df$Stat_2)*df$Direction_2

df$w_1 = sqrt(df$n_1)
df$w_2 = sqrt(df$n_2)

# meta-analysis
df$n = df$n_1 + df$n_2
df$Z = (df$Z_1*df$w_1 + df$Z_2*df$w_2)/sqrt(df$w_1**2 + df$w_2**2)
df$Stat = df$Z**2
df$Direction = sign(df$Z)
df$P = 2*pnorm(abs(df$Z), lower.tail = F)

# save meta-analysis in LDAK summaries format
df_sum = df[, c("Predictor", "A1_1", "A2_2", "Direction", "Stat", "n")]
names(df_sum) = c("Predictor", "A1", "A2", "Direction", "Stat", "n")

write.table(df_sum, snakemake@output[[1]], row.names = FALSE, quote = FALSE)
write.table(df[, c("Predictor", "P")], snakemake@output[[2]], row.names = FALSE, quote = FALSE)