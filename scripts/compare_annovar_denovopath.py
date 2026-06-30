#!/usr/bin/env python3
"""Compare ANNOVAR refGene fields with DeNovoPath DNP consequence fields."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, TextIO


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_text: str) -> Dict[str, str]:
    info: Dict[str, str] = {}
    if not info_text or info_text == ".":
        return info
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "1"
    return info


def split_values(info: Dict[str, str], key: str, n_alts: int, default: str = ".") -> List[str]:
    raw = info.get(key, default)
    if raw in ("", "."):
        return [default] * n_alts
    values = raw.split(",")
    if len(values) == n_alts:
        return values
    if len(values) == 1:
        return values * n_alts
    if len(values) < n_alts:
        return values + [default] * (n_alts - len(values))
    return values[:n_alts]


def decode_annovar_value(value: str) -> str:
    return value.replace("\\x3b", ";")


def annovar_category(func: str, exonic_func: str) -> str:
    func = decode_annovar_value(func or ".").lower()
    exonic_func = decode_annovar_value(exonic_func or ".").lower()
    funcs = set(part for part in func.split(";") if part and part != ".")
    exonic = set(part for part in exonic_func.split(";") if part and part != ".")
    if "splicing" in funcs:
        return "splice"
    if "exonic" in funcs:
        if any("frameshift" in item for item in exonic):
            return "frameshift"
        if any("stopgain" in item or "stoploss" in item for item in exonic):
            return "stop_altering"
        if any("nonsynonymous" in item for item in exonic):
            return "missense"
        if any("synonymous" in item for item in exonic):
            return "synonymous"
        return "coding"
    if "utr5" in funcs or "utr3" in funcs:
        return "utr"
    if "intronic" in funcs:
        return "intron"
    if "upstream" in funcs:
        return "promoter"
    if "downstream" in funcs:
        return "downstream"
    if "intergenic" in funcs:
        return "intergenic"
    return "unknown"


def denovopath_category(consequence: str) -> str:
    consequence = (consequence or ".").lower()
    if consequence in {"frameshift"}:
        return "frameshift"
    if consequence in {
        "stop_gained",
        "stop_gained_early",
        "stop_gained_terminal",
        "stop_lost",
        "stop_lost_readthrough",
        "start_lost",
    }:
        return "stop_altering"
    if consequence == "missense":
        return "missense"
    if consequence in {"synonymous", "stop_retained"}:
        return "synonymous"
    if consequence in {"exon_boundary_disruption", "inframe_deletion", "inframe_insertion", "cds_complex"}:
        return "coding"
    if "splice" in consequence:
        return "splice"
    if consequence in {"utr5", "utr3"}:
        return "utr"
    if consequence == "intron":
        return "intron"
    if consequence == "promoter":
        return "promoter"
    if consequence == "intergenic":
        return "intergenic"
    return "unknown"


def category_match(annovar_cat: str, dnp_cat: str) -> bool:
    if annovar_cat == dnp_cat:
        return True
    if annovar_cat == "downstream" and dnp_cat in {"intergenic", "promoter"}:
        return True
    if annovar_cat == "coding" and dnp_cat in {"missense", "synonymous", "frameshift", "stop_altering", "coding"}:
        return True
    return False


def normalize_gene(value: str) -> str:
    value = decode_annovar_value(value or ".")
    genes = [item for item in value.split(";") if item and item not in {".", "NONE"}]
    return "|".join(sorted(set(genes))) if genes else "."


def dnp_gene_match(annovar_gene: str, dnp_gene: str) -> bool:
    ann_genes = set(normalize_gene(annovar_gene).split("|")) - {"."}
    dnp_genes = set((dnp_gene or ".").split("&")) - {"."}
    return bool(ann_genes and dnp_genes and ann_genes.intersection(dnp_genes))


def compare(vcf: str, max_records: int = 0) -> tuple[dict, list[dict]]:
    records = 0
    allele_rows: list[dict] = []
    ann_counter: Counter[str] = Counter()
    dnp_counter: Counter[str] = Counter()
    pair_counter: Counter[str] = Counter()
    matches = 0
    gene_comparable = 0
    gene_matches = 0
    coding_comparable = 0
    coding_matches = 0
    mismatches: list[dict] = []

    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            records += 1
            chrom, pos, ref, alts = fields[0], fields[1], fields[3], fields[4].split(",")
            info = parse_info(fields[7])
            n_alts = len(alts)
            dnp_conseq = split_values(info, "DNP_CONSEQ", n_alts)
            dnp_gene = split_values(info, "DNP_GENE", n_alts)
            dnp_score = split_values(info, "DNP_SCORE", n_alts)
            dnp_level = split_values(info, "DNP_LEVEL", n_alts)
            func = decode_annovar_value(info.get("Func.refGene", "."))
            gene = decode_annovar_value(info.get("Gene.refGene", "."))
            exonic_func = decode_annovar_value(info.get("ExonicFunc.refGene", "."))
            aa_change = decode_annovar_value(info.get("AAChange.refGene", "."))
            ann_cat = annovar_category(func, exonic_func)
            for idx, alt in enumerate(alts):
                dnp_cat = denovopath_category(dnp_conseq[idx])
                matched = category_match(ann_cat, dnp_cat)
                matches += int(matched)
                ann_counter[ann_cat] += 1
                dnp_counter[dnp_cat] += 1
                pair_counter[f"{ann_cat}->{dnp_cat}"] += 1
                gene_match = dnp_gene_match(gene, dnp_gene[idx])
                if normalize_gene(gene) != "." and dnp_gene[idx] != ".":
                    gene_comparable += 1
                    gene_matches += int(gene_match)
                if ann_cat in {"missense", "synonymous", "frameshift", "stop_altering", "coding"}:
                    coding_comparable += 1
                    coding_matches += int(matched)
                row = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "annovar_func": func,
                    "annovar_exonic_func": exonic_func,
                    "annovar_gene": normalize_gene(gene),
                    "annovar_aa": aa_change,
                    "annovar_category": ann_cat,
                    "dnp_consequence": dnp_conseq[idx],
                    "dnp_category": dnp_cat,
                    "dnp_gene": dnp_gene[idx],
                    "dnp_score": dnp_score[idx],
                    "dnp_level": dnp_level[idx],
                    "category_match": "1" if matched else "0",
                    "gene_match": "1" if gene_match else "0",
                }
                allele_rows.append(row)
                if not matched and len(mismatches) < 30:
                    mismatches.append(row)
            if max_records and records >= max_records:
                break

    total_alleles = len(allele_rows)
    summary = {
        "records": records,
        "alt_alleles": total_alleles,
        "category_matches": matches,
        "category_match_rate": round(matches / total_alleles, 4) if total_alleles else 0.0,
        "coding_comparable_alleles": coding_comparable,
        "coding_match_rate": round(coding_matches / coding_comparable, 4) if coding_comparable else None,
        "gene_comparable_alleles": gene_comparable,
        "gene_match_rate": round(gene_matches / gene_comparable, 4) if gene_comparable else None,
        "annovar_category_counts": dict(ann_counter.most_common()),
        "denovopath_category_counts": dict(dnp_counter.most_common()),
        "category_pairs": dict(pair_counter.most_common()),
        "example_mismatches": mismatches,
    }
    return summary, allele_rows


def write_tsv(rows: Sequence[dict], path: str) -> None:
    if not rows:
        Path(path).write_text("")
        return
    fields = list(rows[0].keys())
    with open(path, "w") as out:
        out.write("\t".join(fields) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="DeNovoPath-scored VCF retaining ANNOVAR INFO fields")
    parser.add_argument("--summary-out", required=True, help="JSON comparison summary")
    parser.add_argument("--details-out", required=True, help="Per-ALT TSV comparison details")
    parser.add_argument("--max-records", type=int, default=0, help="Optional record limit")
    args = parser.parse_args()
    summary, rows = compare(args.vcf, args.max_records)
    with open(args.summary_out, "w") as out:
        json.dump(summary, out, indent=2, sort_keys=True)
        out.write("\n")
    write_tsv(rows, args.details_out)


if __name__ == "__main__":
    main()
