import pandas as pd
from pathlib import Path

def get_geneIDs(file, gene_type, column_name, sheet=None):
    '''Extracts genes names as series, writes them into csv, and returns it'''
    
    filename = Path(file)
    out_file = Path("../data/extracted_{}.csv".format(gene_type))

    if sheet is None:
        df = pd.read_excel(filename)
    else:
        df = pd.read_excel(filename, sheet_name=sheet)

    gene_IDs = df[column_name].dropna().reset_index(drop=True)
    gene_IDs.to_csv(out_file, index=False)

    return gene_IDs

def merge_genes(file, gene_file, gene_type, column_name, sheet=None):
    '''Finds genes of interest, writes them into csv, and returns df'''

    out_file = Path("../data/merged_{}.csv".format(gene_type))

    df = pd.read_csv(file, delimiter="\t")

    gene = get_geneIDs(gene_file, gene_type, column_name, sheet)
    gene.rename("gene", inplace=True)

    gene = df.merge(gene, how="inner", on="gene")
    gene.to_csv(out_file, index=False)

    return gene

def find_subset(df, column_name, condition):
    '''finds the rows in the column with the specified value'''
    return df[df[column_name] == condition]

def final_table(file, gene_type):
    ''' generates df with the count for case/control + variant type for each gene'''

    out_file = Path("../data/{}_case_control_variants_per_gene.csv".format(gene_type))

    df = pd.read_csv(file)
    genes = df["gene"].drop_duplicates().reset_index(drop=True)

    rows = []

    for val in genes.tolist():
        df1 = find_subset(df, "gene", val)
        
        syn = find_subset(df1, "consequence", "synonymous")
        mis = find_subset(df1, "consequence", "missense")
        ptv = find_subset(df1, "consequence", "lof")

        syn_cases = find_subset(syn, "case_control", "Case")
        mis_cases = find_subset(mis, "case_control", "Case")
        ptv_cases = find_subset(ptv, "case_control", "Case")

        syn_control = find_subset(syn, "case_control", "Control")
        mis_control = find_subset(mis, "case_control", "Control")
        ptv_control = find_subset(ptv, "case_control", "Control")

        sample_syn = syn.drop_duplicates(subset="sample")
        sample_mis = mis.drop_duplicates(subset="sample")
        sample_ptv = ptv.drop_duplicates(subset="sample")
        
        new_row = {}
        new_row["gene"] = val
        new_row["cases_synonymous"] = syn_cases.shape[0] 
        new_row["cases_missense"] = mis_cases.shape[0] 
        new_row["cases_PTVs"] = ptv_cases.shape[0]
        new_row["control_synonymous"] = syn_control.shape[0] 
        new_row["control_missense"] = mis_control.shape[0] 
        new_row["control_PTVs"] = ptv_control.shape[0]
        new_row["sample_synonymous"] = sample_syn.shape[0]
        new_row["sample_missense"] = sample_mis.shape[0]
        new_row["sample_PTVs"] = sample_ptv.shape[0]

        rows.append(new_row)

    df = pd.DataFrame(rows)
    df.to_csv(out_file, index=False)

    return df