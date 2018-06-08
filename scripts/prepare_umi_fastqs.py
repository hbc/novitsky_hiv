#!/bin/bash
"""Extract UMIs and prepare fastqs from multiple samples.

Top level script to manage preparation of ready to run bcbio fastqs
organized by UMIs.

This takes an Excel spreadsheet with Sample IDs, read 1 and read 2 tags, and
a directory of files, and processes them into UMI tagged fastqs and bcbio
template CSV files.

Usage:
  bcbio_python prepare_umi_fastqs.py <sample_summary.xlxs> <data directory>
"""
import glob
import os
import subprocess
import sys

import pandas as pd

from bcbio import utils

def main(sample_file, data_dir):
    if sample_file.endswith("csv"):
        df = pd.read_csv(sample_file)
    else:
        xl = pd.ExcelFile(sample_file)
        df = xl.parse(xl.sheet_names[0])
    umi_files = []
    do_rc = 1
    for i, row in df.iterrows():
        base_fastq = os.path.join(data_dir, "%s*" % row["Sample_Project"], "BaseSpace*",
                                  "%s*" % row["Sample_ID"], row["Sample_ID"])
        r1_tag = row["N_index1_setup"][0]
        r2_tag = row["N_index2_setup"][0]
        try:
            fq1 = glob.glob(base_fastq + "*_R1_*.fastq.gz")[0]
            fq2 = glob.glob(base_fastq + "*_R2_*.fastq.gz")[0]
        except IndexError:
            print(base_fastq + "*_R1_*.fastq.gz")
            raise
        sample = "%s_%s" % (row["Sample_Project"], row["Sample_ID"])
        out_dir = utils.safe_makedir(os.path.join(os.getcwd(), "umis", sample))
        out_fq1 = os.path.join(out_dir, "%s_R1.fq.gz" % sample)
        out_fq2 = os.path.join(out_dir, "%s_R2.fq.gz" % sample)
        out_umi = os.path.join(out_dir, "%s-umicounts.csv" % sample)
        if not utils.file_exists(out_umi):
            cmd = [sys.executable, os.path.join(os.path.dirname(__file__, ), "prep_umi_from_adapters.py"),
                   out_dir, sample, fq1, fq2, r1_tag, r2_tag, str(do_rc)]
            print(" ".join(cmd))
            subprocess.check_call(cmd)
        umi_files.append((sample, out_fq1, out_fq2, out_umi))

if __name__ == "__main__":
    main(*sys.argv[1:])
