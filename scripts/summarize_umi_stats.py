#!/usr/bin/env python
"""Summarize UMI statistics and common counts from aligned BAM files.

Simple version of fgbio's consensus creation to count tags.
"""
import csv
import collections
import os
import sys

import pysam

def main(*in_files):
    with open("top_umi_counts.csv", "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["sample", "umi", "position", "count"])
        for in_file in in_files:
            _count_umis(in_file, writer)

def _count_umis(in_file, writer):
    umi_counts = collections.defaultdict(int)
    with pysam.AlignmentFile(in_file, "rb", check_sq=False) as bam_iter:
        for rec in bam_iter:
            umi = _get_umi_tag(rec)
            if umi and not rec.is_unmapped:
                chrom = bam_iter.getrname(rec.reference_id)
                pos = rec.reference_start
                key = (chrom, pos, umi)
                umi_counts[key] += 1
    umis = [(v, c, p, u) for ((c, p, u), v) in umi_counts.items()]
    umis.sort(reverse=True)
    print(os.path.basename(in_file))
    for count, chrom, pos, umi in umis[:50:]:
        writer.writerow([os.path.basename(in_file).replace("-sort.bam", ""), umi, pos, count])

def _get_umi_tag(rec):
    """Handle UMI and duplex tag retrieval.
    """
    for tag in ["RX", "XC"]:
        try:
            return rec.get_tag(tag)
        except KeyError:
            pass

if __name__ == "__main__":
    main(*sys.argv[1:])
