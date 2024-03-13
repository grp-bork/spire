#!/usr/bin/env python3
import sys
import argparse


def parse_args(args):
    """Argument Parser.
    Args:
        args (list): Arguments passed to script.
    Returns:
        Namespace: parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--reference",
        help="Filter input fastqs using this reference.",
        required=True,
        metavar="\b",
    )
    return parser.parse_args()


raw_data_filter_qc_smoothtrim = """ngless "1.0"
import "parallel" version "0.6"
import "mocat" version "0.0"

sample = ARGV[2]
input = load_mocat_sample(ARGV[1] + '/' + sample)

input = preprocess(input) using |read|:
    read = smoothtrim(read, min_quality=20, window=4)
    if len(read) < 45:
        discard

write(input, ofile=sample+'/'+sample+'.filtered.fq.gz')
write(qcstats({fastq}), ofile=sample+'.qc_stats')
"""


def ngless_filter_template(reference):
    return f"""ngless "1.0"
import "parallel" version "0.6"
import "mocat" version "0.0"

sample = ARGV[2]
input = load_mocat_sample(ARGV[1] + '/' + sample)

input = preprocess(input, keep_singles=True) using |read|:
    read = smoothtrim(read, min_quality=20, window=4)
    if len(read) < 45:
        discard
mapped = map(input, fafile='{reference}')

mapped = select(mapped) using |mr|:
    mr = mr.filter(min_match_size=45, min_identity_pc=90, action={{unmatch}})
    if mr.flag({{mapped}}):
        discard

write(as_reads(mapped), ofile=sample+'/'+sample+'.filtered.fq.gz')
"""


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    with open("raw_data_filter.ngl", "w") as f:
        reference_fa = args.reference
        if reference_fa:
            f.write(ngless_filter_template(reference_fa))
        else:
            f.write(raw_data_filter_qc_smoothtrim)
