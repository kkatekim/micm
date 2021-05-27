import pandas as pd
from pathlib import Path

import utils


file = Path("/home/kate/research_data/2020_07_22_ALS_rare_variants_0.01_with_gnomAD_NFE_MPCvalues.txt")

mt = utils.merge_genes(file, "../data/Human.MitoCarta2.0.xls", "mt", "Symbol", "A Human MitoCarta2.0")
trna = utils.merge_genes(file, "../data/tRNAs_ALS.xlsx", "trna", "Gene_name")
mt = utils.final_table("../data/merged_mt.csv", "mt")
trna = utils.final_table("../data/merged_trna.csv", "trna")