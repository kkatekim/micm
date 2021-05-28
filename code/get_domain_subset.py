import pandas as pd
from pathlib import Path

import utils

file_p = Path("../data/datasets/prot2hg_1938_112019.csv")
file_v = Path("../data/datasets/2020_07_22_ALS_rare_variants_0.01_with_gnomAD_NFE_MPCvalues.txt")

protein = pd.read_csv(file_p, delimiter=";")
variants = pd.read_csv(file_v, delimiter="\t")

print(protein.tail())

chr_v = variants["chr"].tolist()
pos_v = variants["pos"].tolist()

chr_p = protein["chrom"].tolist()
# get chr and pos
# compare to chr and domain range in als
# if true, add row and gene/domain name?

df = pd.DataFrame()

for i in range(len(chr_v)):
    for j in range(len(chr_p)):

        pass
    
