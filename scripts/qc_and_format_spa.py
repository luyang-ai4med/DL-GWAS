import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input["gwas"], delim_whitespace=True) # A1 == ALT == Effect allele
df["Stat"] = df.Z**2
df["Direction"] = np.sign(df.Z)

df[["Predictor", "A1", "A2", "Direction", "Stat", "n"]].to_csv(snakemake.output[0], sep = ' ', index = False)
