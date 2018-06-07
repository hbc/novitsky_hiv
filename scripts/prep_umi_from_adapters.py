#!/bin/env python
"""Prepare fastq pairs with name embedded UMIs and adapters stripped.

Output names with barcodes look like:

@M03802:306:000000000-BCD3W:1:1101:15879:1332:UMI_NNNNNNNNN:SAMPLE_NNNN-NNNN 1:N:0:1

Adapter sequences are:

- Read 1 -- index (4 or 5N) - {6bp index} TCTGAAAATCCATATAACACTCCAGTATTTGC - insert
- Read 2 -- index (4 or 5N) - {6bp index} TGGTATCGAAGTCATCCTGCTAG - UMI (10N) - TGGAGTTCATACCCCATCCAAAG - insert

Usage:
    prep_umi_from_adapters.py <out directory> <flowcell name> <fastq read1> <fastq read2>
"""
import bz2
import collections
import gzip
import os
import sys
import subprocess

from Bio.SeqIO.QualityIO import FastqGeneralIterator

def safe_open(f):
    if f.endswith(".gz"):
        return gzip.open(f)
    elif f.endswith(".bz2"):
        return bz2.BZ2File(f)
    else:
        return open(f)

def main(out_dir, fcname, fq1, fq2, r1_tag=4, r2_tag=4):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_fq1 = os.path.join(out_dir, "%s_R1.fq" % fcname)
    out_fq2 = os.path.join(out_dir, "%s_R2.fq" % fcname)
    umis = collections.defaultdict(list)
    with safe_open(fq1) as fq1h:
        with safe_open(fq2) as fq2h:
            with open(out_fq1, "w") as out_fq1h:
                with open(out_fq2, "w") as out_fq2h:
                    it1 = FastqGeneralIterator(fq1h)
                    it2 = FastqGeneralIterator(fq2h)
                    for title1, seq1, qual1 in it1:
                        title2, seq2, qual2 = it2.next()
                        out = extract_umi(title1, seq1, qual1, title2, seq2, qual2, r1_tag, r2_tag)
                        if out:
                            o_name1, o_seq1, o_qual1, o_name2, o_seq2, o_qual2, umi = out
                            out_fq1h.write("@%s\n%s\n+\n%s\n" % (o_name1, o_seq1, o_qual1))
                            out_fq2h.write("@%s\n%s\n+\n%s\n" % (o_name2, o_seq2, o_qual2))
                            umis[umi] += 1
    subprocess.check_call(["bgzip", out_fq1])
    subprocess.check_call(["grabix", "index", "%s.gz" % out_fq1])
    subprocess.check_call(["bgzip", out_fq2])
    subprocess.check_call(["grabix", "index", "%s.gz" % out_fq2])
    with open(os.path.join(out_dir, "%s-umicounts.csv" % fcname), "w") as out_handle:
        for i, (count, umi) in enumerate(sorted([(c, u) for (u, c) in umis.items()], reverse=True)):
            if i < 6:
                print(umi, count)
            out_handle.write("%s,%s\n" % (umi, count))

def extract_umi(n1, s1, q1, n2, s2, q2, r1_tag, r2_tag):
    r1_adapter = "TCTGAAAATCCATATAACACTCCAGTATTTGC"
    r2a_adapter = "TGGTATCGAAGTCATCCTGCTAG"
    r2b_adapter = "TGGAGTTCATACCCCATCCAAAG"
    umi_length = 10
    index_length = 6

    r1_pos = s1.find(r1_adapter)
    r2a_pos = s2.find(r2a_adapter)
    r2b_pos = s2.find(r2b_adapter)
    if r1_pos >= 0 and r2a_pos >= 0 and r2b_pos > r2a_pos:
        tag1 = s1[max([0, r1_pos - r1_tag - index_length]):r1_pos]
        out_s1 = s1[r1_pos + len(r1_adapter):]
        out_q1 = q1[r1_pos + len(r1_adapter):]
        tag2 = s2[max([0, r2a_pos - r2_tag - index_length]):r2a_pos]
        umi = s2[r2a_pos + len(r2a_adapter):r2b_pos]
        if len(umi) == umi_length:
            out_s2 = s2[r2b_pos + len(r2b_adapter):]
            out_q2 = q2[r2b_pos + len(r2b_adapter):]
            out_parts1 = n1.split()
            out_parts1[0] = "%s:UMI_%s:SAMPLE_%s-%s" % (out_parts1[0], umi, tag1, tag2)
            out_n1 = " ".join(out_parts1)
            out_parts2 = n2.split()
            out_parts2[0] = "%s:UMI_%s:SAMPLE_%s-%s" % (out_parts2[0], umi, tag1, tag2)
            out_n2 = " ".join(out_parts2)
            return out_n1, out_s1, out_q1, out_n2, out_s2, out_q2
    #else:
    #    print(r1_pos, r2a_pos, r2b_pos)

if __name__ == "__main__":
    main(*sys.argv[1:])
