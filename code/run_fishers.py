import pandas as pd
import argparse
from pathlib import Path
import hail as hl

import utils

'''Runs fisher's test on variants found in protein domain (includes MPC>=2 and AF_NFE=0).'''

def calculate_fishers(case_carrier, control_carrier, case_noncarrier, control_noncarrier):
    '''Runs fisher's test one one gene and returns a list with pvalue, OR, and CI.'''
    result = hl.eval(hl.fisher_exact_test(case_carrier, case_noncarrier, control_carrier, control_noncarrier))
    return [result["p_value"], result["odds_ratio"], result["ci_95_lower"], result["ci_95_upper"]]


def fishers_test(df, variant):
    '''Runs fisher's test on specified variant and returns df with pval, OR, and CI added.'''
    case = "case_" + variant
    control = "control_" + variant

    fishers_df = pd.DataFrame(index=df.index)
    fishers_df["case_carrier"] = df[case].values
    fishers_df["control_carrier"] = df[control].values
    fishers_df["case_noncarrier"] = 3864 - df[case].values
    fishers_df["control_noncarrier"] = 7839 - df[control].values 
    fishers_df = fishers_df.astype(np.int32)

    # col names
    col_p = "pval_" + variant
    col_or = "OR_" + variant
    col_lowci = "lowci_" + variant
    col_highci = "highci_" + variant
    
    fishers_df[[col_p, col_or, col_lowci, col_highci]] = fishers_df.apply(lambda x: calculate_fishers(
                            x.case_carrier, x.control_carrier, x.case_noncarrier, x.control_noncarrier),
                            axis=1, result_type="expand")

    return fishers_df.drop(labels=["case_carrier", "control_carrier", "case_noncarrier", "control_noncarrier"], axis=1)


def run_fishers_on_variants(df, mpc=False):
    '''Runs fishers test on all variants.'''
    if not mpc:
        variants = ["synonymous", "missense", "PTVs"]
        new_df = utils.get_case_control_per_variant(df)
    else:
        variants = ["missense_mpc>=2"]
        new_df = utils.get_case_control_per_variant(df, mpc)

    fishers_list = []

    for variant in variants:
        tmp = utils.fishers_test(new_df, variant)
        fishers_list.append(tmp)

    return new_df.join(pd.concat(fishers_list, axis=1))


def merge_variants_mpc_fishers(all_df, mpc_df, file_name=None):
    '''Merges df of all variants and MPC>=2 variants and writes to file (if needed)..'''
    combined_df = pd.merge(all_df, mpc_df, on="gene", how="outer")

    if file_name is not None:
        out_file = Path("../data/summaryData/{}_fishers.csv".format(file_name))
        combined_df.to_csv(out_file)

    return combined_df


if __name__ == "__main__":

    # first create case control count
    # find mpc >= 0
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default=None)
    parser.add_argument("-o", "--output", type=str, default=None)
    args = parser.parse_args()

    if args.input is not None:
        docs_file = Path(args.input)
    else:
        docs_file = (
            Path(__file__)
            .resolve()
            .parents[1]
            .joinpath("data", "proteinDomain", "variants_in_protein_domain.csv")
        )
        
    if args.output is not None:
        filename = args.outout
    else:
        filename = "protein_variants"

    df = pd.read_csv(docs_file)
    
    all_df = run_fishers_on_variants(df)
    mpc_df = run_fishers_on_variants(utils.find_subset(df, "MPC", 2, ">="), True)
    combined = merge_variants_mpc_fishers(all_df, mpc_df, filename)
    
    afe_df = utils.find_subset(df, "AF_NFE", 0, "=")
    afe_all_df =  run_fishers_on_variants(afe_df)
    afe_mpc_df = run_fishers_on_variants(utils.find_subset(afe_df, "MPC", 2, ">="), True)
    afe_combined = merge_variants_mpc_fishers(afe_all_df, afe_mpc_df, "afe_{}".format(filename))

    print(df.shape[0])
    print(combined.shape[0])
    print(afe_combined.shape[0])