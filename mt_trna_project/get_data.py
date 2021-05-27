import pandas as pd
from pathlib import Path

import subset


file = Path("/home/kate/2020_07_22_ALS_rare_variants_0.01_with_gnomAD_NFE_MPCvalues.txt")

mt = subset.merge_genes(file, "../data/Human.MitoCarta2.0.xls", "mt", "Symbol", "A Human MitoCarta2.0")
trna = subset.merge_genes(file, "../data/tRNAs_ALS.xlsx", "trna", "Gene_name")
mt = subset.final_table("../data/merged_mt.csv", "mt")
trna = subset.final_table("../data/merged_trna.csv", "trna")