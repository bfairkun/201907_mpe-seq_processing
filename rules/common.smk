import pandas as pd
import os

###### Config file and sample sheets #####
configfile: "config.yaml"

output_dir = config["output_dir"]
samples = pd.read_csv(config["samples"],sep='\t', index_col=0)

# # How to access values in samples.tsv

# print(samples)
# print( expand("Hello {sample}", sample=samples.index) )
# print( samples.at["A", "R1"] )

def GetRNASeqFastq(wildcards):
    R1 = samples.at[wildcards.sample, "R1"]
    R2 = samples.at[wildcards.sample, "R2"]
    return([R1, R2])

