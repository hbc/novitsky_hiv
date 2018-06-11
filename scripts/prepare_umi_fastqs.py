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
import csv
import glob
import gzip
import os
import subprocess
import sys

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import editdistance

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
    for sample, fq1, fq2, umi_file in umi_files:
        for umi_group in prepare_umi_groups(umi_file):
            out_fq1 = extract_umi_group(fq1, umi_group)
            out_fq2 = extract_umi_group(fq2, umi_group)
            print(sample, umi_group[0], out_fq1, out_fq2)

def extract_umi_group(in_fq, umi_group):
    base_dir, base_name = os.path.split(in_fq)
    out_dir = utils.safe_makedir(os.path.join(base_dir, "split"))
    base_name = utils.splitext_plus(base_name)[0]
    out_fq = os.path.join(out_dir, "%s_%s.fq" % (base_name, umi_group[0]))
    if not utils.file_exists(out_fq + ".gz.gbi"):
        with gzip.open(in_fq) as in_fqh:
            with open(out_fq, "w") as out_fqh:
                for name, seq, qual in FastqGeneralIterator(in_fqh):
                    if any([name.find("UMI_%s" % u) >= 0 for u in umi_group]):
                     out_fqh.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        subprocess.check_call(["bgzip", "-f", out_fq])
        subprocess.check_call(["grabix", "index", "%s.gz" % out_fq])
    return out_fq

def prepare_umi_groups(umi_file):
    """Group together UMIs with less than the allowed number of edits.
    """
    min_umi_count = 50
    min_add_multiplier = 5
    allowed_edits = 2
    max_groups = 10
    umis = []
    counts = []
    with open(umi_file) as in_handle:
        reader = csv.reader(in_handle)
        for umi, count in reader:
            count = int(count)
            if count > min_umi_count:
                added = False
                for i, gs in enumerate(umis):
                    if not added:
                        for g in gs:
                            if editdistance.eval(umi, g) <= allowed_edits:
                                umis[i].append(umi)
                                counts[i] += count
                                added = True
                                break
                if not added:
                    umis.append([umi])
                    counts.append(count)
    out = []
    for count, umis in sorted(zip(counts, umis), reverse=True):
        if count > min_umi_count * min_add_multiplier:
            out.append(umis)
        if len(out) >= max_groups:
            break
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
