#!/usr/bin/env python3
"""Normalize biallelic SNP VCF REF/ALT against a FASTA reference.

For biallelic SNPs, records are handled as follows:
- REF equals FASTA base: keep unchanged
- ALT equals FASTA base: swap REF/ALT and flip diploid/haploid GT alleles
- Neither allele equals FASTA base: keep unchanged and count as unfixable

This is intended for PLINK-like VCFs whose REF may be a provisional allele.
Existing functional annotations such as ANN are retained but may describe the
pre-normalized allele orientation; downstream comparisons should account for
that limitation or re-run the external annotator after normalization.
"""

from __future__ import annotations

import argparse
import gzip
import json
from typing import Dict, List, TextIO

from pyfaidx import Fasta


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", compresslevel=1)
    return open(path, mode)


def flip_gt(gt: str) -> str:
    if gt in {".", "./.", ".|."}:
        return gt
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)
    flipped = []
    for part in parts:
        if part == "0":
            flipped.append("1")
        elif part == "1":
            flipped.append("0")
        else:
            flipped.append(part)
    return sep.join(flipped)


def flip_sample(sample: str, format_keys: List[str]) -> str:
    values = sample.split(":")
    try:
        gt_idx = format_keys.index("GT")
    except ValueError:
        return sample
    if gt_idx >= len(values):
        return sample
    values[gt_idx] = flip_gt(values[gt_idx])
    return ":".join(values)


def normalize(args: argparse.Namespace) -> Dict[str, object]:
    fasta = Fasta(args.reference)
    stats: Dict[str, object] = {
        "records_seen": 0,
        "records_written": 0,
        "ref_match": 0,
        "ref_alt_swapped": 0,
        "unfixable": 0,
        "non_snp_or_multiallelic": 0,
        "missing_contig": 0,
        "examples": [],
    }
    with open_text(args.vcf, "r") as inp, open_text(args.output, "w") as out:
        for line in inp:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                out.write('##INFO=<ID=REFNORM_STATUS,Number=1,Type=String,Description="REF normalization status before DeNovoPath scoring: match, swap, unfixable, non_snp_or_multiallelic, missing_contig">\n')
                out.write('##DeNovoPathRefNormalization="biallelic SNP REF/ALT normalized against FASTA; ANN retained from input"\n')
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                out.write(line)
                continue
            if args.limit and int(stats["records_seen"]) >= args.limit:
                break
            stats["records_seen"] = int(stats["records_seen"]) + 1
            chrom, pos_text, ref, alt = fields[0], fields[1], fields[3], fields[4]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                stats["non_snp_or_multiallelic"] = int(stats["non_snp_or_multiallelic"]) + 1
                fields[7] = append_info(fields[7], "REFNORM_STATUS", "non_snp_or_multiallelic")
                out.write("\t".join(fields) + "\n")
                stats["records_written"] = int(stats["records_written"]) + 1
                continue
            try:
                fasta_ref = str(fasta[chrom][int(pos_text) - 1 : int(pos_text)]).upper()
            except KeyError:
                stats["missing_contig"] = int(stats["missing_contig"]) + 1
                fields[7] = append_info(fields[7], "REFNORM_STATUS", "missing_contig")
                out.write("\t".join(fields) + "\n")
                stats["records_written"] = int(stats["records_written"]) + 1
                continue
            ref_u = ref.upper()
            alt_u = alt.upper()
            if ref_u == fasta_ref:
                stats["ref_match"] = int(stats["ref_match"]) + 1
                fields[7] = append_info(fields[7], "REFNORM_STATUS", "match")
            elif alt_u == fasta_ref:
                fields[3], fields[4] = alt, ref
                fields[7] = append_info(fields[7], "REFNORM_STATUS", "swap")
                if len(fields) > 8:
                    format_keys = fields[8].split(":")
                    fields[9:] = [flip_sample(sample, format_keys) for sample in fields[9:]]
                stats["ref_alt_swapped"] = int(stats["ref_alt_swapped"]) + 1
            else:
                stats["unfixable"] = int(stats["unfixable"]) + 1
                fields[7] = append_info(fields[7], "REFNORM_STATUS", "unfixable")
                examples = stats["examples"]
                if isinstance(examples, list) and len(examples) < 10:
                    examples.append(f"{chrom}:{pos_text} VCF_REF={ref} VCF_ALT={alt} FASTA_REF={fasta_ref}")
            out.write("\t".join(fields) + "\n")
            stats["records_written"] = int(stats["records_written"]) + 1
    return stats


def append_info(info: str, key: str, value: str) -> str:
    parts = [] if info in {"", "."} else [part for part in info.split(";") if not part.startswith(f"{key}=")]
    parts.append(f"{key}={value}")
    return ";".join(parts) if parts else "."


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    stats = normalize(args)
    with open(args.summary_out, "w") as out:
        json.dump(stats, out, indent=2, sort_keys=True)
        out.write("\n")


if __name__ == "__main__":
    main()
