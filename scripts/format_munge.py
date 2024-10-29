import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input["GWAS"], delim_whitespace=True)

df["Z"] = np.sqrt(df.Stat)*df.Direction

df = df.rename(columns = {"Predictor": "SNP", "n": "N"})

df[["SNP", "A1", "A2", "Z", "N"]].to_csv(snakemake.output[0], sep = '\t', index = False)
