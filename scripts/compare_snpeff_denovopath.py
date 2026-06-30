#!/usr/bin/env python3
"""Compare snpEff ANN annotations with DeNovoPath DNP consequence fields."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, TextIO


def open_text(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


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
    if raw in {"", "."}:
        return [default] * n_alts
    values = raw.split(",")
    if len(values) == n_alts:
        return values
    if len(values) == 1:
        return values * n_alts
    if len(values) < n_alts:
        return values + [default] * (n_alts - len(values))
    return values[:n_alts]


def snpeff_category(annotation: str) -> str:
    terms = {item.strip().lower() for item in (annotation or "").split("&") if item.strip()}
    if not terms:
        return "unknown"
    if "frameshift_variant" in terms:
        return "frameshift"
    if terms.intersection({"stop_gained", "stop_lost", "start_lost", "initiator_codon_variant"}):
        return "stop_altering"
    if "missense_variant" in terms:
        return "missense"
    if terms.intersection({"synonymous_variant", "stop_retained_variant"}):
        return "synonymous"
    if any("splice" in item for item in terms):
        return "splice"
    if terms.intersection(
        {
            "conservative_inframe_insertion",
            "disruptive_inframe_insertion",
            "conservative_inframe_deletion",
            "disruptive_inframe_deletion",
            "coding_sequence_variant",
            "protein_protein_contact",
            "structural_interaction_variant",
        }
    ):
        return "coding"
    if any("utr" in item for item in terms):
        return "utr"
    if terms.intersection({"intron_variant"}):
        return "intron"
    if terms.intersection({"upstream_gene_variant"}):
        return "promoter"
    if terms.intersection({"downstream_gene_variant"}):
        return "downstream"
    if terms.intersection({"intergenic_region"}):
        return "intergenic"
    return "unknown"


def denovopath_category(consequence: str) -> str:
    consequence = (consequence or ".").lower()
    if consequence == "frameshift":
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


def category_match(snpeff_cat: str, dnp_cat: str) -> bool:
    if snpeff_cat == dnp_cat:
        return True
    if snpeff_cat == "downstream" and dnp_cat in {"intergenic", "promoter", "utr"}:
        return True
    if snpeff_cat == "coding" and dnp_cat in {"missense", "synonymous", "frameshift", "stop_altering", "coding"}:
        return True
    return False


def parse_ann_entries(info: Dict[str, str], alts: Sequence[str]) -> Dict[str, List[Dict[str, str]]]:
    by_alt: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    raw = info.get("ANN", "")
    if not raw or raw == ".":
        return by_alt
    for entry in raw.split(","):
        fields = entry.split("|")
        if len(fields) < 16:
            fields += [""] * (16 - len(fields))
        row = {
            "allele": fields[0],
            "annotation": fields[1],
            "impact": fields[2],
            "gene_name": fields[3],
            "gene_id": fields[4],
            "feature_type": fields[5],
            "feature_id": fields[6],
            "biotype": fields[7],
            "rank": fields[8],
            "hgvs_c": fields[9],
            "hgvs_p": fields[10],
            "distance": fields[14],
        }
        by_alt[row["allele"]].append(row)
    for alt in alts:
        by_alt.setdefault(alt, [])
    return by_alt


def severity_rank(category: str) -> int:
    ranks = {
        "stop_altering": 8,
        "frameshift": 7,
        "splice": 6,
        "missense": 5,
        "coding": 4,
        "synonymous": 3,
        "utr": 2,
        "intron": 1,
        "promoter": 1,
        "downstream": 0,
        "intergenic": 0,
        "unknown": -1,
    }
    return ranks.get(category, -1)


def choose_best_ann(entries: Sequence[Dict[str, str]]) -> Dict[str, str]:
    if not entries:
        return {
            "annotation": ".",
            "category": "unknown",
            "gene_id": ".",
            "gene_name": ".",
            "feature_id": ".",
            "hgvs_c": ".",
            "hgvs_p": ".",
        }
    return max(entries, key=lambda row: severity_rank(snpeff_category(row.get("annotation", ""))))


def normalize_gene(value: str) -> str:
    if not value or value == ".":
        return "."
    return value


def gene_match(snpeff_gene: str, dnp_gene: str) -> bool:
    snp = normalize_gene(snpeff_gene)
    dnp = normalize_gene(dnp_gene)
    if snp == "." or dnp == ".":
        return False
    return snp in set(dnp.split("&"))


def update_group(summary: Dict[str, object], row: Dict[str, str]) -> None:
    summary["alt_alleles"] += 1
    summary["category_matches"] += int(row["category_match"] == "1")
    summary["snpeff_category_counts"][row["snpeff_category"]] += 1
    summary["denovopath_category_counts"][row["dnp_category"]] += 1
    summary["category_pairs"][f"{row['snpeff_category']}->{row['dnp_category']}"] += 1
    if row["snpeff_gene"] != "." and row["dnp_gene"] != ".":
        summary["gene_comparable_alleles"] += 1
        summary["gene_matches"] += int(row["gene_match"] == "1")
    if row["snpeff_category"] in {"missense", "synonymous", "frameshift", "stop_altering", "coding"}:
        summary["coding_comparable_alleles"] += 1
        summary["coding_matches"] += int(row["category_match"] == "1")


def empty_group() -> Dict[str, object]:
    return {
        "alt_alleles": 0,
        "category_matches": 0,
        "gene_comparable_alleles": 0,
        "gene_matches": 0,
        "coding_comparable_alleles": 0,
        "coding_matches": 0,
        "snpeff_category_counts": Counter(),
        "denovopath_category_counts": Counter(),
        "category_pairs": Counter(),
    }


def finalize_group(group: Dict[str, object]) -> Dict[str, object]:
    total = int(group["alt_alleles"])
    gene_total = int(group["gene_comparable_alleles"])
    coding_total = int(group["coding_comparable_alleles"])
    return {
        "alt_alleles": total,
        "category_matches": int(group["category_matches"]),
        "category_match_rate": round(int(group["category_matches"]) / total, 4) if total else 0.0,
        "gene_comparable_alleles": gene_total,
        "gene_match_rate": round(int(group["gene_matches"]) / gene_total, 4) if gene_total else None,
        "coding_comparable_alleles": coding_total,
        "coding_match_rate": round(int(group["coding_matches"]) / coding_total, 4) if coding_total else None,
        "snpeff_category_counts": dict(group["snpeff_category_counts"].most_common()),
        "denovopath_category_counts": dict(group["denovopath_category_counts"].most_common()),
        "category_pairs": dict(group["category_pairs"].most_common()),
    }


def compare(vcf: str, max_records: int = 0) -> tuple[dict, list[dict]]:
    records = 0
    rows: list[dict] = []
    groups = {"all": empty_group(), "ref_match": empty_group(), "ref_swap": empty_group(), "ref_other": empty_group()}
    mismatches: list[dict] = []
    refnorm_counter: Counter[str] = Counter()
    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            records += 1
            if max_records and records > max_records:
                break
            chrom, pos, ref, alts = fields[0], fields[1], fields[3], fields[4].split(",")
            info = parse_info(fields[7])
            refnorm = info.get("REFNORM_STATUS", "unmarked")
            refnorm_counter[refnorm] += 1
            n_alts = len(alts)
            dnp_conseq = split_values(info, "DNP_CONSEQ", n_alts)
            dnp_gene = split_values(info, "DNP_GENE", n_alts)
            dnp_score = split_values(info, "DNP_SCORE", n_alts)
            dnp_level = split_values(info, "DNP_LEVEL", n_alts)
            ann_by_alt = parse_ann_entries(info, alts)
            for idx, alt in enumerate(alts):
                best = choose_best_ann(ann_by_alt.get(alt, []))
                snp_cat = snpeff_category(best.get("annotation", ""))
                dnp_cat = denovopath_category(dnp_conseq[idx])
                matched = category_match(snp_cat, dnp_cat)
                genes_match = gene_match(best.get("gene_id", "."), dnp_gene[idx])
                row = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "refnorm_status": refnorm,
                    "snpeff_annotation": best.get("annotation", "."),
                    "snpeff_category": snp_cat,
                    "snpeff_impact": best.get("impact", "."),
                    "snpeff_gene": best.get("gene_id", "."),
                    "snpeff_gene_name": best.get("gene_name", "."),
                    "snpeff_feature": best.get("feature_id", "."),
                    "snpeff_hgvs_c": best.get("hgvs_c", "."),
                    "snpeff_hgvs_p": best.get("hgvs_p", "."),
                    "dnp_consequence": dnp_conseq[idx],
                    "dnp_category": dnp_cat,
                    "dnp_gene": dnp_gene[idx],
                    "dnp_score": dnp_score[idx],
                    "dnp_level": dnp_level[idx],
                    "category_match": "1" if matched else "0",
                    "gene_match": "1" if genes_match else "0",
                }
                rows.append(row)
                update_group(groups["all"], row)
                if refnorm == "match":
                    update_group(groups["ref_match"], row)
                elif refnorm == "swap":
                    update_group(groups["ref_swap"], row)
                else:
                    update_group(groups["ref_other"], row)
                if not matched and len(mismatches) < 30:
                    mismatches.append(row)
    summary = {
        "records": records,
        "refnorm_status_counts": dict(refnorm_counter.most_common()),
        "all": finalize_group(groups["all"]),
        "ref_match": finalize_group(groups["ref_match"]),
        "ref_swap": finalize_group(groups["ref_swap"]),
        "ref_other": finalize_group(groups["ref_other"]),
        "example_mismatches": mismatches,
    }
    return summary, rows


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
    parser.add_argument("--vcf", required=True, help="DeNovoPath-scored VCF retaining snpEff ANN")
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--details-out", required=True)
    parser.add_argument("--max-records", type=int, default=0)
    args = parser.parse_args()
    summary, rows = compare(args.vcf, args.max_records)
    with open(args.summary_out, "w") as out:
        json.dump(summary, out, indent=2, sort_keys=True)
        out.write("\n")
    write_tsv(rows, args.details_out)


if __name__ == "__main__":
    main()
