import pandas as pd
import utils 

def remove_significant_syn(df):

    adjusted_p = 0.05/df.shape[0]
    return utils.find_subset(df, "pval_synonymous", adjusted_p, ">=")


def plot_qq_for_all_variants(df, title, remove_syn=False):

    df1 = df
    variants = ["synonymous", "missense","PTVs", "missense_MPC>=2"]

    if remove_syn:
        df1 = remove_significant_syn(df1)
    
    for variant in variants:
        #print(df.shape[0], df["pval_{}".format(variant)].dropna().shape[0], no_syn["pval_{}".format(variant)].dropna().shape[0])
        utils.plot_qq(df1, title, variant)


if __name__ == "__main__":

    df = pd.read_csv("../data/proteinDomain/variants_missense_mpc_greaterequal_2_fishers.csv")
    no_syn = remove_significant_syn(df)
    print(df.shape[0], no_syn.shape[0])
    print(df.dropna(subset=["pval_missense_MPC>=2"]).shape[0], no_syn.dropna(subset=["pval_missense_MPC>=2"]).shape[0])
    plot_qq_for_all_variants(df, "Variants")
    plot_qq_for_all_variants(df, "Variants (significant synonymous genes removed)", True)

    df1 = pd.read_csv("../data/proteinDomain/variants_af_nfe_mpc_combined_fishers.csv")
    no_syn1 = remove_significant_syn(df1)
    print(df1.shape[0], no_syn1.shape[0])
    print(df1.dropna(subset=["pval_missense_MPC>=2"]).shape[0], no_syn1.dropna(subset=["pval_missense_MPC>=2"]).shape[0])
    plot_qq_for_all_variants(df, "Variants AF_NFE=0")
    plot_qq_for_all_variants(df, "Variants AF_NFE=0 (significant synonymous genes removed)", True)       