import pandas as pd

import utils

def find_af_nfe_equal_0(df, out_file):
    df = utils.find_subset(df, "AF_NFE", 0, "equal")
    df.to_csv(out_file, index=False)

def get_cases_controls(filename, gene_type):
    df = utils.final_table(filename, gene_type)

def run_fishers(filename, out_file):
    utils.all_stats(filename, out_file)

def get_af_nfe_equals_0():
    df = pd.read_csv("../data/datasets/variants_in_domain.csv")
    find_af_nfe_equal_0(df, "../data/datasets/variants_with_af_nfe_0.csv")
    get_cases_controls("../data/datasets/variants_with_af_nfe_0.csv", "variants_with_0_af_nfe")
    run_fishers("../data/variants_with_0_af_nfe_case_control_variants_per_gene.csv", "../data/datasets/variants_with_0_af_nfe_fishers.csv")

def get_mpc_greaterequal_2():
    df = pd.read_csv("../data/datasets/variants_in_domain.csv")
    df = utils.find_subset(df, "MPC", 2, "greater equal")
    df.to_csv("../data/datasets/variants_mpc_greaterequal2.csv", index=False)
    df = utils.final_table("../data/datasets/variants_mpc_greaterequal2.csv","variants_mpc")
    utils.all_stats("../data/variants_mpc_case_control_variants_per_gene.csv", "../data/datasets/variants_mpc_greaterequal_2_fishers.csv")

def merge_missense_mpc(mpc_file, all_file, out_file):
    df = pd.read_csv(mpc_file)

    cols = ["gene", "cases_missense", "control_missense", "sample_missense", 
            "pval_missense", "OR_missense", "lowci_missense", "highci_missense"]
    df = df[cols]
    df.rename(columns={"cases_missense": "cases_missense_MPC>=2", "control_missense": "control_missense_MPC>=2",
                            "sample_missense": "sample_missense_MPC>=2", "pval_missense": "pval_missense_MPC>=2",
                            "OR_missense": "OR_missense_MPC>=2", "lowci_missense": "lowci_missense_MPC>=2",
                            "highci_missense": "highci_missense_MPC>=2"}, inplace=True)

    df1 = pd.read_csv(all_file)
    df = pd.merge(df1, df, how="outer", on="gene")
    df.to_csv(out_file, index=False)

def get_mpc_greaterequal_2_for_af_nfe_0():
    af_nfe = pd.read_csv("../data/datasets/variants_with_af_nfe_0.csv")
    mpc = pd.read_csv("../data/datasets/variants_mpc_greaterequal2.csv")

    # merge the af_nfe df with the mpc since i want to find the intersection of these
    af_nfe_mpc = pd.merge(af_nfe, mpc, how="inner")
    af_nfe_mpc.to_csv("../data/datasets/variants_af_nfe_mpc.csv", index=False)
    df = utils.final_table("../data/datasets/variants_af_nfe_mpc.csv","variants_af_nfe_mpc")
    utils.all_stats("../data/variants_af_nfe_mpc_case_control_variants_per_gene.csv", "../data/datasets/variants_af_nfe_mpc_fishers.csv")


if __name__ == "__main__":
    merge_missense_mpc("../data/datasets/variants_mpc_greaterequal_2_fishers.csv",
                            "../data/datasets/protein_variants_fishers.csv",
                            "../data/datasets/variants_missense_mpc_greaterequal_2_fishers.csv")

    merge_missense_mpc("../data/datasets/variants_af_nfe_mpc_fishers.csv",
                             "../data/datasets/variants_with_0_af_nfe_fishers.csv",
                             "../data/datasets/variants_af_nfe_mpc_combined_fishers.csv")