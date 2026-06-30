#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import heapq
import itertools
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple


def open_text(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def parse_info(info_text: str) -> Dict[str, str]:
    info = {}
    if info_text in {"", "."}:
        return info
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "1"
    return info


def split_info_values(info: Dict[str, str], key: str, n: int, default: str = ".") -> List[str]:
    value = info.get(key, default)
    values = value.split(",") if value != "." else [default]
    if len(values) < n:
        values.extend([default] * (n - len(values)))
    return values[:n]


def to_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except ValueError:
        return default


def to_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except ValueError:
        return default


def row_rank_key(row: Dict[str, str]) -> Tuple[float, float, float]:
    return (to_float(row["score"]), to_float(row["confidence"]), to_float(row["qc"]))


def export_ranked(args: argparse.Namespace) -> None:
    rows = []
    top_heap = []
    counter = itertools.count()
    gene_stats = defaultdict(
        lambda: {
            "n_variants": 0,
            "n_high": 0,
            "n_moderate": 0,
            "max_score": 0.0,
            "sum_score": 0.0,
            "frameshift": 0,
            "stop_gained": 0,
            "missense": 0,
            "splice": 0,
        }
    )

    with open_text(args.vcf, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, var_id, ref, alt_text, qual, filt, info_text = fields[:8]
            alts = alt_text.split(",")
            info = parse_info(info_text)
            scores = split_info_values(info, "DNP_SCORE", len(alts), "0")
            impacts = split_info_values(info, "DNP_IMPACT", len(alts), "0")
            prots = split_info_values(info, "DNP_PROT", len(alts), "0")
            granthams = split_info_values(info, "DNP_GRANTHAM", len(alts), "0")
            blosums = split_info_values(info, "DNP_BLOSUM", len(alts), "0")
            codonuses = split_info_values(info, "DNP_CODONUSE", len(alts), "0")
            protctxs = split_info_values(info, "DNP_PROTCTX", len(alts), "0")
            structs = split_info_values(info, "DNP_STRUCT", len(alts), "0")
            af_structs = split_info_values(info, "DNP_AFSTRUCT", len(alts), "0")
            esms = split_info_values(info, "DNP_ESM", len(alts), "0")
            protein_lms = split_info_values(info, "DNP_PROTLM", len(alts), "0")
            domains = split_info_values(info, "DNP_DOMAIN", len(alts), "0")
            splices = split_info_values(info, "DNP_SPLICE", len(alts), "0")
            splice_motifs = split_info_values(info, "DNP_SPLICE_MOTIF", len(alts), "0")
            splice_pwms = split_info_values(info, "DNP_SPLICE_PWM", len(alts), "0")
            splice_maxents = split_info_values(info, "DNP_SPLICE_MAXENT", len(alts), "0")
            splice_auxs = split_info_values(info, "DNP_SPLICE_AUX", len(alts), "0")
            splice_eses = split_info_values(info, "DNP_SPLICE_ESE", len(alts), "0")
            utrs = split_info_values(info, "DNP_UTR", len(alts), "0")
            rnafolds = split_info_values(info, "DNP_RNAFOLD", len(alts), "0")
            mirnas = split_info_values(info, "DNP_MIRNA", len(alts), "0")
            promoters = split_info_values(info, "DNP_PROM", len(alts), "0")
            seqs = split_info_values(info, "DNP_SEQ", len(alts), "0")
            kmers = split_info_values(info, "DNP_KMER", len(alts), "0")
            repeats = split_info_values(info, "DNP_REPEAT", len(alts), "0")
            mutctxs = split_info_values(info, "DNP_MUTCTX", len(alts), "0")
            dna_lms = split_info_values(info, "DNP_DNALM", len(alts), "0")
            cohorts = split_info_values(info, "DNP_COHORT", len(alts), "0")
            hwes = split_info_values(info, "DNP_HWE", len(alts), "0")
            het_obs = split_info_values(info, "DNP_HETOBS", len(alts), "0")
            het_exp = split_info_values(info, "DNP_HETEXP", len(alts), "0")
            het_dev = split_info_values(info, "DNP_HETDEV", len(alts), "0")
            fis_values = split_info_values(info, "DNP_FIS", len(alts), "0")
            fsts = split_info_values(info, "DNP_FST", len(alts), "0")
            case_ctrls = split_info_values(info, "DNP_CASECTRL", len(alts), "0")
            pis = split_info_values(info, "DNP_PI", len(alts), "0")
            thetas = split_info_values(info, "DNP_THETA", len(alts), "0")
            tajds = split_info_values(info, "DNP_TAJD", len(alts), "0")
            lds = split_info_values(info, "DNP_LD", len(alts), "0")
            haps = split_info_values(info, "DNP_HAP", len(alts), "0")
            gene_lofs = split_info_values(info, "DNP_GENELOF", len(alts), "0")
            gene_missenses = split_info_values(info, "DNP_GENEMIS", len(alts), "0")
            gene_constraints = split_info_values(info, "DNP_GENECON", len(alts), "0")
            qcs = split_info_values(info, "DNP_QC", len(alts), "0")
            confs = split_info_values(info, "DNP_CONF", len(alts), "0")
            mls = split_info_values(info, "DNP_ML", len(alts), "0")
            calibrateds = split_info_values(info, "DNP_CAL", len(alts), "0")
            uncertainties = split_info_values(info, "DNP_UNCERT", len(alts), "0")
            oods = split_info_values(info, "DNP_OOD", len(alts), "0")
            levels = split_info_values(info, "DNP_LEVEL", len(alts), ".")
            consequences = split_info_values(info, "DNP_CONSEQ", len(alts), ".")
            genes = split_info_values(info, "DNP_GENE", len(alts), ".")
            txs = split_info_values(info, "DNP_TX", len(alts), ".")
            all_txs = split_info_values(info, "DNP_ALLTX", len(alts), ".")
            aas = split_info_values(info, "DNP_AA", len(alts), ".")
            codons = split_info_values(info, "DNP_CODON", len(alts), ".")
            hgvs = split_info_values(info, "DNP_HGVS", len(alts), ".")
            norms = split_info_values(info, "DNP_NORM", len(alts), ".")
            domain_ids = split_info_values(info, "DNP_DOMID", len(alts), ".")
            structure_ids = split_info_values(info, "DNP_AFID", len(alts), ".")
            esm_ids = split_info_values(info, "DNP_ESMID", len(alts), ".")
            mirna_ids = split_info_values(info, "DNP_MIRID", len(alts), ".")
            feature_importances = split_info_values(info, "DNP_FEATIMP", len(alts), ".")
            ctx96s = split_info_values(info, "DNP_96CTX", len(alts), ".")
            subafs = split_info_values(info, "DNP_SUBAF", len(alts), ".")
            privates = split_info_values(info, "DNP_PRIVATE", len(alts), ".")
            caseafs = split_info_values(info, "DNP_CASEAF", len(alts), ".")
            acs = split_info_values(info, "DNP_AC", len(alts), "0")
            ans = split_info_values(info, "DNP_AN", len(alts), "0")
            afs = split_info_values(info, "DNP_AF", len(alts), "0")
            maf_bins = split_info_values(info, "DNP_MAFBIN", len(alts), ".")
            carrs = split_info_values(info, "DNP_CARR", len(alts), "0")

            for i, alt in enumerate(alts):
                score = to_float(scores[i])
                qc = to_float(qcs[i])
                ac = to_int(acs[i])
                if score < args.min_score:
                    continue
                if qc < args.min_qc:
                    continue
                if ac < args.min_ac:
                    continue
                gene = genes[i]
                consequence = consequences[i]
                row = {
                    "chrom": chrom,
                    "pos": pos,
                    "id": var_id,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "filter": filt,
                    "score": f"{score:.4f}",
                    "impact": impacts[i],
                    "protein": prots[i],
                    "grantham": granthams[i],
                    "blosum": blosums[i],
                    "codon_usage": codonuses[i],
                    "protein_context": protctxs[i],
                    "protein_structure": structs[i],
                    "protein_structure_model": af_structs[i],
                    "protein_esm": esms[i],
                    "protein_lm": protein_lms[i],
                    "protein_domain": domains[i],
                    "splice": splices[i],
                    "splice_motif": splice_motifs[i],
                    "splice_pwm": splice_pwms[i],
                    "splice_maxent": splice_maxents[i],
                    "splice_aux": splice_auxs[i],
                    "splice_ese": splice_eses[i],
                    "utr": utrs[i],
                    "rnafold": rnafolds[i],
                    "mirna": mirnas[i],
                    "promoter": promoters[i],
                    "sequence": seqs[i],
                    "kmer": kmers[i],
                    "repeat": repeats[i],
                    "mutation_context_score": mutctxs[i],
                    "dna_lm": dna_lms[i],
                    "cohort": cohorts[i],
                    "hwe": hwes[i],
                    "heterozygosity_observed": het_obs[i],
                    "heterozygosity_expected": het_exp[i],
                    "heterozygosity_deviation": het_dev[i],
                    "fis": fis_values[i],
                    "fst": fsts[i],
                    "case_control": case_ctrls[i],
                    "pi": pis[i],
                    "theta": thetas[i],
                    "tajima_d": tajds[i],
                    "ld": lds[i],
                    "haplotype": haps[i],
                    "gene_lof_oe": gene_lofs[i],
                    "gene_missense_oe": gene_missenses[i],
                    "gene_constraint": gene_constraints[i],
                    "qc": qcs[i],
                    "confidence": confs[i],
                    "ml": mls[i],
                    "calibrated": calibrateds[i],
                    "uncertainty": uncertainties[i],
                    "ood": oods[i],
                    "level": levels[i],
                    "consequence": consequence,
                    "gene": gene,
                    "transcript": txs[i],
                    "all_transcripts": all_txs[i],
                    "aa_change": aas[i],
                    "codon_change": codons[i],
                    "hgvs_like": hgvs[i],
                    "normalized_variant": norms[i],
                    "protein_domain_id": domain_ids[i],
                    "protein_structure_id": structure_ids[i],
                    "protein_esm_id": esm_ids[i],
                    "mirna_id": mirna_ids[i],
                    "feature_importance": feature_importances[i],
                    "context_96": ctx96s[i],
                    "subpopulation_af": subafs[i],
                    "private_shared": privates[i],
                    "case_control_af": caseafs[i],
                    "ac": str(ac),
                    "an": ans[i],
                    "af": afs[i],
                    "maf_bin": maf_bins[i],
                    "carriers": carrs[i],
                }
                if args.top > 0:
                    entry = (row_rank_key(row), next(counter), row)
                    if len(top_heap) < args.top:
                        heapq.heappush(top_heap, entry)
                    elif entry[0] > top_heap[0][0]:
                        heapq.heapreplace(top_heap, entry)
                else:
                    rows.append(row)
                if gene != ".":
                    stat = gene_stats[gene]
                    stat["n_variants"] += 1
                    stat["sum_score"] += score
                    stat["max_score"] = max(stat["max_score"], score)
                    if levels[i] == "HIGH":
                        stat["n_high"] += 1
                    if levels[i] == "MODERATE":
                        stat["n_moderate"] += 1
                    if consequence in stat:
                        stat[consequence] += 1
                    if consequence.startswith("stop_gained_"):
                        stat["stop_gained"] += 1
                    if consequence.startswith("splice") or consequence == "exon_boundary_disruption":
                        stat["splice"] += 1

    if args.top > 0:
        rows = [entry[2] for entry in sorted(top_heap, key=lambda entry: (entry[0], entry[1]), reverse=True)]
    else:
        rows.sort(key=row_rank_key, reverse=True)

    variant_columns = [
        "chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "qual",
        "filter",
        "score",
        "impact",
        "protein",
        "grantham",
        "blosum",
        "codon_usage",
        "protein_context",
        "protein_structure",
        "protein_structure_model",
        "protein_esm",
        "protein_lm",
        "protein_domain",
        "splice",
        "splice_motif",
        "splice_pwm",
        "splice_maxent",
        "splice_aux",
        "splice_ese",
        "utr",
        "rnafold",
        "mirna",
        "promoter",
        "sequence",
        "kmer",
        "repeat",
        "mutation_context_score",
        "dna_lm",
        "cohort",
        "hwe",
        "heterozygosity_observed",
        "heterozygosity_expected",
        "heterozygosity_deviation",
        "fis",
        "fst",
        "case_control",
        "pi",
        "theta",
        "tajima_d",
        "ld",
        "haplotype",
        "gene_lof_oe",
        "gene_missense_oe",
        "gene_constraint",
        "qc",
        "confidence",
        "ml",
        "calibrated",
        "uncertainty",
        "ood",
        "level",
        "consequence",
        "gene",
        "transcript",
        "all_transcripts",
        "aa_change",
        "codon_change",
        "hgvs_like",
        "normalized_variant",
        "protein_domain_id",
        "protein_structure_id",
        "protein_esm_id",
        "mirna_id",
        "feature_importance",
        "context_96",
        "subpopulation_af",
        "private_shared",
        "case_control_af",
        "ac",
        "an",
        "af",
        "maf_bin",
        "carriers",
    ]
    with open(args.variants_out, "w") as out:
        out.write("\t".join(variant_columns) + "\n")
        for row in rows:
            out.write("\t".join(row[column] for column in variant_columns) + "\n")

    if args.genes_out:
        gene_rows = []
        for gene, stat in gene_stats.items():
            n = stat["n_variants"]
            if n == 0:
                continue
            gene_rows.append(
                {
                    "gene": gene,
                    "n_variants": str(n),
                    "n_high": str(stat["n_high"]),
                    "n_moderate": str(stat["n_moderate"]),
                    "max_score": f"{stat['max_score']:.4f}",
                    "mean_score": f"{stat['sum_score'] / n:.4f}",
                    "frameshift": str(stat["frameshift"]),
                    "stop_gained": str(stat["stop_gained"]),
                    "missense": str(stat["missense"]),
                    "splice": str(stat["splice"]),
                }
            )
        gene_rows.sort(
            key=lambda row: (
                int(row["n_high"]),
                to_float(row["max_score"]),
                int(row["n_moderate"]),
                int(row["n_variants"]),
            ),
            reverse=True,
        )
        gene_columns = [
            "gene",
            "n_variants",
            "n_high",
            "n_moderate",
            "max_score",
            "mean_score",
            "frameshift",
            "stop_gained",
            "missense",
            "splice",
        ]
        with open(args.genes_out, "w") as out:
            out.write("\t".join(gene_columns) + "\n")
            for row in gene_rows:
                out.write("\t".join(row[column] for column in gene_columns) + "\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export ranked DeNovoPath variants and gene burden TSVs.")
    parser.add_argument("--vcf", required=True, help="Scored VCF/VCF.GZ containing DNP_* INFO fields")
    parser.add_argument("--variants-out", required=True, help="Ranked variant TSV output")
    parser.add_argument("--genes-out", help="Gene burden TSV output")
    parser.add_argument("--min-score", type=float, default=0.0, help="Minimum DNP_SCORE to export")
    parser.add_argument("--min-qc", type=float, default=0.0, help="Minimum DNP_QC to export")
    parser.add_argument("--min-ac", type=int, default=1, help="Minimum DNP_AC to export; default keeps carried ALT alleles only")
    parser.add_argument("--top", type=int, default=0, help="Keep only top N variants after filtering")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.top < 0:
        raise SystemExit("--top must be >= 0")
    export_ranked(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
