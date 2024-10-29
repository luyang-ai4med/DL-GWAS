pheno = snakemake@wildcards[["pheno"]]
cov_df = read.table(snakemake@input[["cov"]], header = T)
sample_df = read.table(snakemake@input[["samples"]])
names(sample_df) = c("FID", "IID")
pheno_df = read.table(snakemake@input[["pheno_file"]], sep=',', header=T)
pheno_sex = pheno_df[pheno_df$index == as.integer(pheno), "sex"]
if (pheno_sex == "F") {
    cov_df = cov_df[cov_df$sex==0, ]
    sample_df = sample_df[sample_df$FID %in% cov_df$FID, ]
} else if (pheno_sex == "M")  {
    cov_df = cov_df[cov_df$sex==1, ]
    sample_df = sample_df[sample_df$FID %in% cov_df$FID, ]
}

TRAIN_PATH = paste0(snakemake@params[["PATH"]], pheno, "/", pheno, "_train_popdx_all.csv")
VAL_PATH = paste0(snakemake@params[["PATH"]], pheno, "/", pheno, "_val_popdx_all.csv")
TRAIN_PATH_SAMPLE = paste0(snakemake@params[["PATH"]], pheno, "/", pheno, "_train_patients.txt")
TEST_PATH_SAMPLE = paste0(snakemake@params[["PATH"]], pheno, "/", pheno, "_test_patients.txt")
VAL_PATH_SAMPLE = paste0(snakemake@params[["PATH"]], pheno, "/", pheno, "_val_patients.txt")

train_df = read.csv(TRAIN_PATH)
train_s = read.table(TRAIN_PATH_SAMPLE)
test_s = read.table(TEST_PATH_SAMPLE)
names(train_s) = c("FID")
names(test_s) = c("FID")
train_s = rbind(train_s, test_s)
train_df = cbind(train_s, train_df)
train_df = merge(train_df, cov_df, by = "FID")
train_df = train_df[train_df$FID %in% sample_df$FID, ]

val_df = read.csv(VAL_PATH)
val_s = read.table(VAL_PATH_SAMPLE)
names(val_s) = c("FID")
val_df = cbind(val_s, val_df)
val_df = merge(val_df, cov_df, by = "FID")
val_df = val_df[val_df$FID %in% sample_df$FID, ]

#### Training set ####

# prepare quantitative phenotypes for controls
train_q_df = train_df[train_df$true_labels_2019 == 0, ]
train_q_df$pheno = scale(train_q_df$POPDx)
train_q_df$pheno_res = residuals(lm(paste0("pheno~sex+age+", paste0("PC", 1:20, collapse="+")),
                                    data=train_q_df, na.action=na.exclude))
write.table(train_q_df[, c("FID", "IID", "pheno_res")], snakemake@output[["train_q"]], quote = F, row.names = F)

# prepare binary phenotype & covariates file
train_df$binary = 1
train_df[train_df$true_labels_2019 == 1, "binary"] = 2
write.table(train_df[, c("FID", "IID", "binary")], snakemake@output[["train_b"]], quote = F, row.names = F)

if (pheno_sex != "B") {
    write.table(train_df[, c("FID", "IID", "age", paste0("PC", 1:20))], snakemake@output[["train_cov"]], quote = F, row.names = F)} else { 
    write.table(train_df[, c("FID", "IID", "sex", "age", paste0("PC", 1:20))], snakemake@output[["train_cov"]], quote = F, row.names = F)}



#### Training & validation set combined ####
val_q_df = val_df[val_df$true_labels_2019 == 0, ]

val_df$binary = 1
val_df[val_df$true_labels_2019 == 1, "binary"] = 2

# prepare quantitative phenotypes for controls
cols = c("FID", "IID", "POPDx", "sex", "age", paste0("PC", 1:20))
all_q_df = rbind(train_q_df[, cols], val_q_df[, cols])
all_q_df$pheno = scale(all_q_df$POPDx)
all_q_df$pheno_res = residuals(lm(paste0("pheno~sex+age+", paste0("PC", 1:20, collapse="+")),
                                    data=all_q_df, na.action=na.exclude))
write.table(all_q_df[, c("FID", "IID", "pheno_res")], snakemake@output[["all_q"]], quote = F, row.names = F)

# prepare quantitative phenotypes for all
cols = c("FID", "IID", "POPDx", "sex", "age", paste0("PC", 1:20))
all_df = rbind(train_df[, cols], val_df[, cols])
all_df$pheno = scale(all_df$POPDx)
all_df$pheno_res = residuals(lm(paste0("pheno~sex+age+", paste0("PC", 1:20, collapse="+")),
                                    data=all_df, na.action=na.exclude))
write.table(all_df[, c("FID", "IID", "pheno_res")], snakemake@output[["all_all_q"]], quote = F, row.names = F)


# prepare binary phenotype & covariates file
cols = c("FID", "IID", "binary", "sex", "age", paste0("PC", 1:20))
all_df = rbind(train_df[, cols], val_df[, cols])
write.table(all_df[, c("FID", "IID", "binary")], snakemake@output[["all_b"]], quote = F, row.names = F)
if (pheno_sex != "B") {
    write.table(all_df[, c("FID", "IID", "age", paste0("PC", 1:20))], snakemake@output[["all_cov"]], quote = F, row.names = F)} else {write.table(all_df[, c("FID", "IID", "sex", "age", paste0("PC", 1:20))], snakemake@output[["all_cov"]], quote = F, row.names = F)}


