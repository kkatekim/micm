import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.stats as stats
import hail as hl

# helper functions. read docstring for usage. TODO: finish docstring
# TODO: change functions to have parameter df and return df. save files automatically in a newly created folder

def get_geneIDs(file, gene_type, column_name, sheet=None):
    '''Extracts genes names as series, writes them into csv, and returns it'''
    
    out_file = Path("../data/{}_gene_IDs.csv".format(gene_type))

    if sheet is None:
        df = pd.read_excel(file)
    else:
        df = pd.read_excel(file, sheet_name=sheet)

    gene_IDs = df[column_name].dropna().reset_index(drop=True)
    gene_IDs.to_csv(out_file, index=False)

    return gene_IDs


def merge_genes(file, gene_file, gene_type, column_name, sheet=None):
    '''Finds genes of interest, writes them into csv, and returns df'''

    gene_file = Path(gene_file)
    out_file = Path("../data/merged_{}.csv".format(gene_type))

    df = pd.read_csv(file, delimiter="\t")

    gene = get_geneIDs(gene_file, gene_type, column_name, sheet)
    gene.rename("gene", inplace=True)

    gene = df.merge(gene, how="inner", on="gene")
    gene.to_csv(out_file, index=False)

    return gene


def find_subset(df, column_name, factor, condition):
    '''finds the rows in the column with the specified value'''

    if condition == "=":
        return df[df[column_name] == factor]
    
    elif condition == "!=":
        return df[df[column_name] != factor]
    
    elif condition == "<":
        return df[df[column_name] < factor]

    elif condition == ">":
        return df[df[column_name] > factor]

    elif condition == "<=":
        return df[df[column_name] <= factor]

    elif condition == ">=":
        return df[df[column_name] >= factor]


def final_table(file, gene_type):
    ''' generates df with the count for case/control + variant type for each gene'''

    out_file = Path("../data/{}_case_control_variants_per_gene.csv".format(gene_type))

    df = pd.read_csv(file)

    # list of all genes to use in for loop
    genes = df["gene"].drop_duplicates().reset_index(drop=True)

    rows = []

    # TODO: rewrite this using groupby and size()

    for val in genes.tolist():
        df1 = find_subset(df, "gene", val, "=")
        
        syn = find_subset(df1, "consequence", "synonymous", "=")
        mis = find_subset(df1, "consequence", "missense", "=")
        ptv = find_subset(df1, "consequence", "lof", "=")

        syn_cases = find_subset(syn, "case_control", "Case", "=")
        mis_cases = find_subset(mis, "case_control", "Case", "=")
        ptv_cases = find_subset(ptv, "case_control", "Case", "=")

        syn_control = find_subset(syn, "case_control", "Control", "=")
        mis_control = find_subset(mis, "case_control", "Control", "=")
        ptv_control = find_subset(ptv, "case_control", "Control", "=")

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


def get_chr_num(chrom):
    '''subsets hg19 chrom string to find just chrom number''' 
    if "_" in chrom:
        return chrom[3:chrom.index("_")]
    return chrom[3:]


def fishers_test(df, variant):
    ''' runs fisher's test on one column from dataframe (ie syn, mis, PTV)'''

    p_list = []
    oddsratio_list = []
    lowci_list = []
    highci_list = []

    case = "cases_" + variant
    control = "control_" + variant

    for i in range(df.shape[0]):

        case_carrier = df[case].iloc[i].astype(np.int32)
        control_carrier = df[control].iloc[i].astype(np.int32)
        case_noncarrier = int(3864 - df[case].iloc[i].astype(np.int32))
        control_noncarrier = int(7839 - df[control].iloc[i].astype(np.int32))

        result = hl.eval(hl.fisher_exact_test(case_carrier, case_noncarrier, control_carrier, control_noncarrier))
        p_list.append(hl.eval(result["p_value"]))
        oddsratio_list.append(hl.eval(result["odds_ratio"]))
        lowci_list.append(hl.eval(result["ci_95_lower"]))
        highci_list.append(hl.eval(result["ci_95_upper"]))

    return (p_list, oddsratio_list, lowci_list, highci_list)


def all_stats(in_file, out_file):
    ''' runs fisher's on whole dataframe and writes it to file'''

    df = pd.read_csv(in_file)

    variants = ["synonymous", "missense", "PTVs"]

    p_list = []
    oddsratio_list = []
    lowci_list = []
    highci_list = []

    for variant in variants:
        col_p = "pval_" + variant
        col_or = "OR_" + variant
        col_lowci = "lowci_" + variant
        col_highci = "highci_" + variant

        p_variant, or_variant, low_ci, high_ci = fishers_test(df, variant)

        df[col_p] = p_variant
        df[col_or] = or_variant
        df[col_lowci] = low_ci
        df[col_highci] = high_ci

    df.sort_values(by="gene", inplace=True)
    df.to_csv(out_file, index=False)

    return df


def plot_qq(filename, name):
    '''plots QQ graph of expected vs observed using plt scatter plot. doesn't include CI rn'''

    df = pd.read_csv(filename)

    variants = ["synonymous", "missense", "PTVs"]

    for variant in variants:

        # finds relevant column names
        col_p = "pval_" + variant
        col_low = "lowci_" + variant
        col_high = "highci_" + variant

        # removes all NaN
        df.dropna(subset=[col_p, col_low, col_high], inplace=True)
        df.sort_values(col_p, inplace=True)

        # -log to get observed values
        p_list = df[col_p].to_numpy()
        obs = -1 * np.log10(p_list)
        exp = -1 * np.log10(np.arange(1, p_list.shape[0]+1) / p_list.shape[0])
        
        fig, ax = plt.subplots()
        ax.scatter(exp, obs)
        xpoints = ypoints = ax.get_xlim()
        ax.plot(xpoints, ypoints, color='black', scalex=False, scaley=False)
        plt.title("{} {} Q-Q Plot".format(name, variant))
        plt.xlabel("Expected value")
        plt.ylabel("Observed value")
        plt.savefig("../data/datasets/figures/{}_{}_QQ.png".format(name, variant))
        plt.show()


def find_significant_genes(df, alpha, variant, name):
    '''find genes less than p value based on mutation variant'''
    
    # multiple test correction, divide alpha by # of genes
    adjusted_p = alpha/df.shape[0]

    col = "pval_" + variant
    df.sort_values(col, ignore_index=True, inplace=True)

    significant_genes = find_subset(df, col, adjusted_p, "<=")
    significant_genes = find_subset(significant_genes, "pval_synonymous", adjusted_p, ">")

    significant_genes.to_csv("../data/datasets/{}_{}_significant_genes.csv".format(name, variant), index=False)
    
    return significant_genes