#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd


def parse_args(args):
    """Argument Parser.
    Args:
        args (list): Arguments passed to script.
    Returns:
        Namespace: parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--max_csslevel_file",
        help="MaxCSS output file from GUNC.",
        required=True,
        metavar="\b",
    )
    parser.add_argument(
        "-d",
        "--gunc_detailed_output_dir",
        help="GUNC detailed output dir.",
        required=True,
        metavar="\b",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        help="File in which to write outputfile with added score.",
        required=True,
        metavar="\b",
    )
    return parser.parse_args()


def get_pass_GUNC5_score(detail_file):
    df = pd.read_csv(detail_file, sep="\t", header=0)
    df = df.drop(["pass.GUNC"], axis=1)
    max_CSS = df.iloc[[0]].to_dict("records")[0]
    df = df[df["taxonomic_level"] != "species"]
    df = df[df["contamination_portion"] > 0.05]
    if len(df) > 0:
        max_CSSidx = df["clade_separation_score"].idxmax()
        if not pd.isna(max_CSSidx):
            max_CSS = df.loc[[max_CSSidx]].to_dict("records")[0]
            if max_CSS["clade_separation_score"] > 0.45:
                max_CSS["pass.GUNC5"] = False
                return max_CSS
    max_CSS["pass.GUNC5"] = True
    return max_CSS


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    max_csslevel_file = pd.read_csv(args.max_csslevel_file, sep="\t", header=0)
    max_csslevel_file = max_csslevel_file.to_dict("index")

    for row in max_csslevel_file:
        pass_gunc = max_csslevel_file[row]["pass.GUNC"]
        if pass_gunc or pd.isna(pass_gunc):
            max_csslevel_file[row]["pass.GUNC5"] = True
        elif not pass_gunc:
            bin_name = max_csslevel_file[row]["genome"]
            detail_file = os.path.join(
                args.gunc_detailed_output_dir,
                bin_name + ".progenomes_2.1.all_levels.tsv",
            )
            new_row = get_pass_GUNC5_score(detail_file)
            max_csslevel_file[row] = new_row
    out_df = pd.DataFrame.from_dict(max_csslevel_file, orient="index")
    out_df.to_csv(args.output_file, index=False, sep="\t", na_rep="nan")
