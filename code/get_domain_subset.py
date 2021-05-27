import pandas as pd
from pathlib import Path

import utils

file = Path("/home/kate/research_data/prot2hg_1938_112019.csv")

df = pd.read_csv(file, delimiter=";")

