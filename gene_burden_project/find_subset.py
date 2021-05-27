import pandas as pd
from pathlib import Path

def find_subset(df, column_name, condition):
    '''finds the rows in the column with the specified value'''
    return df[df[column_name] == condition]

