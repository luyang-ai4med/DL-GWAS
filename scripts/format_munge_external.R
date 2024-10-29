df_external = read.table(snakemake@input[["GWAS_external"]], header = T)
snp_df = read.csv(snakemake@input[["snp_map"]], sep = '\t')
df_external = merge(df_external, snp_df, by.x = "SNP", by.y = "rsid")
df_external["Direction"] = sign(df_external$b)
df_external["Stat"] = (df_external$b/df_external$se)**2
df_external$Z = df_external$b/df_external$se
df_external =df_external[, c("SNP", "A1", "A2", "Z", "N", "p")]
write.table(df_external, snakemake@output[[1]], sep = '\t', row.names = F, quote = F)
