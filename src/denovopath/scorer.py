from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pyfaidx import Fasta
import yaml


DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

AA_HYDROPATHY = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}

AA_CHARGE = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1,
    "H": 0.5,
    "A": 0,
    "C": 0,
    "F": 0,
    "G": 0,
    "I": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "S": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0,
}

AA_VOLUME = {
    "A": 88.6,
    "R": 173.4,
    "N": 114.1,
    "D": 111.1,
    "C": 108.5,
    "Q": 143.8,
    "E": 138.4,
    "G": 60.1,
    "H": 153.2,
    "I": 166.7,
    "L": 166.7,
    "K": 168.6,
    "M": 162.9,
    "F": 189.9,
    "P": 112.7,
    "S": 89.0,
    "T": 116.1,
    "W": 227.8,
    "Y": 193.6,
    "V": 140.0,
}

CONSEQUENCE_BASE_SCORE = {
    "stop_gained": 1.00,
    "stop_gained_early": 1.00,
    "stop_gained_terminal": 0.82,
    "frameshift": 0.98,
    "stop_lost": 0.95,
    "stop_lost_readthrough": 0.95,
    "start_lost": 0.92,
    "splice_acceptor_donor": 0.92,
    "exon_boundary_disruption": 0.90,
    "stop_retained": 0.04,
    "inframe_indel": 0.66,
    "inframe_deletion": 0.68,
    "inframe_insertion": 0.62,
    "missense": 0.58,
    "splice_region": 0.56,
    "cds_complex": 0.50,
    "utr5": 0.26,
    "utr3": 0.22,
    "promoter": 0.18,
    "intron": 0.10,
    "synonymous": 0.05,
    "intergenic": 0.03,
    "unknown": 0.08,
}

GRANTHAM_DISTANCE = {
    ("A", "R"): 112, ("A", "N"): 111, ("A", "D"): 126, ("A", "C"): 195, ("A", "Q"): 91,
    ("A", "E"): 107, ("A", "G"): 60, ("A", "H"): 86, ("A", "I"): 94, ("A", "L"): 96,
    ("A", "K"): 106, ("A", "M"): 84, ("A", "F"): 113, ("A", "P"): 27, ("A", "S"): 99,
    ("A", "T"): 58, ("A", "W"): 148, ("A", "Y"): 112, ("A", "V"): 64,
    ("R", "N"): 86, ("R", "D"): 96, ("R", "C"): 180, ("R", "Q"): 43, ("R", "E"): 54,
    ("R", "G"): 125, ("R", "H"): 29, ("R", "I"): 97, ("R", "L"): 102, ("R", "K"): 26,
    ("R", "M"): 91, ("R", "F"): 97, ("R", "P"): 103, ("R", "S"): 110, ("R", "T"): 71,
    ("R", "W"): 101, ("R", "Y"): 77, ("R", "V"): 96,
    ("N", "D"): 23, ("N", "C"): 139, ("N", "Q"): 46, ("N", "E"): 42, ("N", "G"): 80,
    ("N", "H"): 68, ("N", "I"): 149, ("N", "L"): 153, ("N", "K"): 94, ("N", "M"): 142,
    ("N", "F"): 158, ("N", "P"): 91, ("N", "S"): 46, ("N", "T"): 65, ("N", "W"): 174,
    ("N", "Y"): 143, ("N", "V"): 133,
    ("D", "C"): 154, ("D", "Q"): 61, ("D", "E"): 45, ("D", "G"): 94, ("D", "H"): 81,
    ("D", "I"): 168, ("D", "L"): 172, ("D", "K"): 101, ("D", "M"): 160, ("D", "F"): 177,
    ("D", "P"): 108, ("D", "S"): 65, ("D", "T"): 85, ("D", "W"): 181, ("D", "Y"): 160,
    ("D", "V"): 152,
    ("C", "Q"): 154, ("C", "E"): 170, ("C", "G"): 159, ("C", "H"): 174, ("C", "I"): 198,
    ("C", "L"): 198, ("C", "K"): 202, ("C", "M"): 196, ("C", "F"): 205, ("C", "P"): 169,
    ("C", "S"): 112, ("C", "T"): 149, ("C", "W"): 215, ("C", "Y"): 194, ("C", "V"): 192,
    ("Q", "E"): 29, ("Q", "G"): 87, ("Q", "H"): 24, ("Q", "I"): 109, ("Q", "L"): 113,
    ("Q", "K"): 53, ("Q", "M"): 101, ("Q", "F"): 116, ("Q", "P"): 76, ("Q", "S"): 68,
    ("Q", "T"): 42, ("Q", "W"): 130, ("Q", "Y"): 99, ("Q", "V"): 96,
    ("E", "G"): 98, ("E", "H"): 40, ("E", "I"): 134, ("E", "L"): 138, ("E", "K"): 56,
    ("E", "M"): 126, ("E", "F"): 140, ("E", "P"): 93, ("E", "S"): 80, ("E", "T"): 65,
    ("E", "W"): 152, ("E", "Y"): 122, ("E", "V"): 121,
    ("G", "H"): 98, ("G", "I"): 135, ("G", "L"): 138, ("G", "K"): 127, ("G", "M"): 127,
    ("G", "F"): 153, ("G", "P"): 42, ("G", "S"): 56, ("G", "T"): 59, ("G", "W"): 184,
    ("G", "Y"): 147, ("G", "V"): 109,
    ("H", "I"): 94, ("H", "L"): 99, ("H", "K"): 32, ("H", "M"): 87, ("H", "F"): 100,
    ("H", "P"): 77, ("H", "S"): 89, ("H", "T"): 47, ("H", "W"): 115, ("H", "Y"): 83,
    ("H", "V"): 84,
    ("I", "L"): 5, ("I", "K"): 102, ("I", "M"): 10, ("I", "F"): 21, ("I", "P"): 95,
    ("I", "S"): 142, ("I", "T"): 89, ("I", "W"): 61, ("I", "Y"): 33, ("I", "V"): 29,
    ("L", "K"): 107, ("L", "M"): 15, ("L", "F"): 22, ("L", "P"): 98, ("L", "S"): 145,
    ("L", "T"): 92, ("L", "W"): 61, ("L", "Y"): 36, ("L", "V"): 32,
    ("K", "M"): 95, ("K", "F"): 102, ("K", "P"): 103, ("K", "S"): 121, ("K", "T"): 78,
    ("K", "W"): 110, ("K", "Y"): 85, ("K", "V"): 97,
    ("M", "F"): 28, ("M", "P"): 87, ("M", "S"): 135, ("M", "T"): 81, ("M", "W"): 67,
    ("M", "Y"): 36, ("M", "V"): 21,
    ("F", "P"): 114, ("F", "S"): 155, ("F", "T"): 103, ("F", "W"): 40, ("F", "Y"): 22,
    ("F", "V"): 50,
    ("P", "S"): 74, ("P", "T"): 38, ("P", "W"): 147, ("P", "Y"): 110, ("P", "V"): 68,
    ("S", "T"): 58, ("S", "W"): 177, ("S", "Y"): 144, ("S", "V"): 124,
    ("T", "W"): 128, ("T", "Y"): 92, ("T", "V"): 69,
    ("W", "Y"): 37, ("W", "V"): 88,
    ("Y", "V"): 55,
}

BLOSUM62_ROWS = {
    "A": " 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0",
    "R": "-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3",
    "N": "-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3",
    "D": "-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3",
    "C": " 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
    "Q": "-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2",
    "E": "-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2",
    "G": " 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3",
    "H": "-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3",
    "I": "-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3",
    "L": "-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1",
    "K": "-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2",
    "M": "-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1",
    "F": "-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1",
    "P": "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2",
    "S": " 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2",
    "T": " 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0",
    "W": "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3",
    "Y": "-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1",
    "V": " 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4",
}
AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")
DEFAULT_MAX_CODON_USAGE_TRAINING_CODONS = 1_000_000
DEFAULT_MAX_PROTEIN_LM_TRAINING_RESIDUES = 500_000
DEFAULT_MAX_SPLICE_PWM_TRAINING_INTRONS = 20_000
FAST_GZIP_COMPRESSLEVEL = 1
TANDEM_REPEAT_PATTERNS = {
    unit_len: re.compile(r"([ACGT]{%d})(?:\1){2,}" % unit_len)
    for unit_len in range(1, 7)
}
BLOSUM62 = {
    (row_aa, col_aa): int(value)
    for row_aa, row in BLOSUM62_ROWS.items()
    for col_aa, value in zip(AA_ORDER, row.split())
}

ML_FEATURES = [
    "impact_score",
    "protein_score",
    "grantham_score",
    "blosum_score",
    "codon_usage_score",
    "protein_context_score",
    "protein_structure_score",
    "protein_structure_model_score",
    "protein_esm_score",
    "protein_lm_score",
    "protein_domain_score",
    "splice_score",
    "splice_motif_score",
    "splice_pwm_score",
    "splice_maxent_score",
    "splice_aux_score",
    "splice_ese_score",
    "utr_score",
    "rnafold_score",
    "mirna_score",
    "promoter_score",
    "sequence_score",
    "kmer_score",
    "repeat_score",
    "mutation_context_score",
    "dna_lm_score",
    "cohort_score",
    "hwe_score",
    "heterozygosity_deviation_score",
    "fst_score",
    "case_control_score",
    "window_pi",
    "window_theta",
    "window_tajima_d",
    "window_ld",
    "window_haplotype",
    "gene_lof_oe",
    "gene_missense_oe",
    "gene_constraint_score",
    "qc_score",
    "confidence_score",
]


@dataclass
class ScoreConfig:
    score_weights: Dict[str, float] = field(
        default_factory=lambda: {
            "impact": 0.42,
            "protein": 0.28,
            "splice": 0.15,
            "sequence": 0.10,
        }
    )
    cohort_weights: Dict[str, float] = field(
        default_factory=lambda: {
            "single_sample": 0.02,
            "small_cohort": 0.05,
            "large_cohort": 0.12,
        }
    )
    level_thresholds: Dict[str, float] = field(
        default_factory=lambda: {
            "high": 0.80,
            "moderate": 0.50,
            "low": 0.20,
        }
    )
    confidence_weights: Dict[str, float] = field(
        default_factory=lambda: {
            "qc": 0.70,
            "score_separation": 0.30,
        }
    )
    protein_weights: Dict[str, float] = field(
        default_factory=lambda: {
            "physchem": 0.50,
            "grantham": 0.25,
            "blosum": 0.25,
            "codon_usage": 0.10,
            "protein_context": 0.10,
            "protein_structure": 0.12,
            "protein_structure_model": 0.15,
            "protein_esm": 0.15,
            "protein_lm": 0.10,
            "protein_domain": 0.12,
        }
    )
    sequence_weights: Dict[str, float] = field(
        default_factory=lambda: {
            "context": 0.50,
            "kmer": 0.20,
            "repeat": 0.15,
            "mutctx": 0.15,
            "dna_lm": 0.15,
        }
    )
    max_signal_weight: float = 0.55
    weighted_signal_weight: float = 0.45
    transcript_priority: str = "score"


@dataclass
class Region:
    chrom: str
    start: int = 1
    end: Optional[int] = None

    def contains(self, chrom: str, pos: int) -> bool:
        if chrom != self.chrom:
            return False
        if pos < self.start:
            return False
        if self.end is not None and pos > self.end:
            return False
        return True


def _merge_float_dict(base: Dict[str, float], updates: object, section: str) -> Dict[str, float]:
    merged = dict(base)
    if updates is None:
        return merged
    if not isinstance(updates, dict):
        raise ValueError(f"Config section '{section}' must be a mapping")
    for key, value in updates.items():
        if key not in merged:
            raise ValueError(f"Unknown config key '{section}.{key}'")
        try:
            merged[key] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Config value '{section}.{key}' must be numeric") from exc
    return merged


def load_score_config(path: Optional[str]) -> ScoreConfig:
    config = ScoreConfig()
    if not path:
        return config
    with open(path) as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("Config file must contain a YAML mapping")
    config.score_weights = _merge_float_dict(config.score_weights, data.get("score_weights"), "score_weights")
    config.cohort_weights = _merge_float_dict(config.cohort_weights, data.get("cohort_weights"), "cohort_weights")
    config.level_thresholds = _merge_float_dict(config.level_thresholds, data.get("level_thresholds"), "level_thresholds")
    config.confidence_weights = _merge_float_dict(
        config.confidence_weights, data.get("confidence_weights"), "confidence_weights"
    )
    config.protein_weights = _merge_float_dict(config.protein_weights, data.get("protein_weights"), "protein_weights")
    config.sequence_weights = _merge_float_dict(config.sequence_weights, data.get("sequence_weights"), "sequence_weights")
    if "max_signal_weight" in data:
        config.max_signal_weight = float(data["max_signal_weight"])
    if "weighted_signal_weight" in data:
        config.weighted_signal_weight = float(data["weighted_signal_weight"])
    if "transcript_priority" in data:
        config.transcript_priority = str(data["transcript_priority"])
    valid_priorities = {"score", "longest_cds", "longest_transcript", "first"}
    if config.transcript_priority not in valid_priorities:
        allowed = ", ".join(sorted(valid_priorities))
        raise ValueError(f"transcript_priority must be one of: {allowed}")
    if config.max_signal_weight < 0 or config.weighted_signal_weight < 0:
        raise ValueError("Signal blending weights must be non-negative")
    if config.max_signal_weight + config.weighted_signal_weight <= 0:
        raise ValueError("At least one signal blending weight must be positive")
    for section_name, values in (
        ("score_weights", config.score_weights),
        ("cohort_weights", config.cohort_weights),
        ("confidence_weights", config.confidence_weights),
        ("protein_weights", config.protein_weights),
        ("sequence_weights", config.sequence_weights),
    ):
        if any(value < 0 for value in values.values()):
            raise ValueError(f"{section_name} values must be non-negative")
        if sum(values.values()) <= 0:
            raise ValueError(f"At least one {section_name} value must be positive")
    high = config.level_thresholds["high"]
    moderate = config.level_thresholds["moderate"]
    low = config.level_thresholds["low"]
    if not (0 <= low <= moderate <= high <= 1):
        raise ValueError("level_thresholds must satisfy 0 <= low <= moderate <= high <= 1")
    return config


def parse_region(text: str) -> Region:
    text = text.replace(",", "").strip()
    if not text:
        raise ValueError("Empty region string")
    if ":" not in text:
        return Region(chrom=text)
    chrom, span = text.split(":", 1)
    if not chrom or not span:
        raise ValueError(f"Invalid region '{text}'")
    if "-" in span:
        start_text, end_text = span.split("-", 1)
        start = int(start_text) if start_text else 1
        end = int(end_text) if end_text else None
    else:
        start = int(span)
        end = start
    if start < 1:
        raise ValueError(f"Region start must be >= 1 in '{text}'")
    if end is not None and end < start:
        raise ValueError(f"Region end must be >= start in '{text}'")
    return Region(chrom=chrom, start=start, end=end)


def in_selected_regions(chrom: str, pos: int, regions: Sequence[Region]) -> bool:
    if not regions:
        return True
    return any(region.contains(chrom, pos) for region in regions)


def active_methods(
    n_samples: int,
    has_sample_info: bool = False,
    has_phenotype: bool = False,
    has_protein_domains: bool = False,
    has_protein_structures: bool = False,
    has_protein_esm_scores: bool = False,
    has_mirna_sites: bool = False,
    has_ml_model: bool = False,
    has_population_windows: bool = True,
    has_gene_constraint: bool = True,
) -> List[str]:
    methods = [
        "impact_consequence",
        "protein_physchem",
        "protein_grantham",
        "protein_blosum62",
        "protein_codon_usage",
        "protein_low_complexity_context",
        "protein_structure_heuristic",
        "protein_lm_kmer_delta_proxy",
        "splice_boundary",
        "splice_motif_delta",
        "splice_species_pwm",
        "splice_maxent_like",
        "splice_branch_polypyrimidine",
        "splice_exonic_enhancer_silencer",
        "utr_start_polyadenylation",
        "promoter_core_motifs",
        "annotation_minimal_allele_normalization",
        "annotation_stop_altering_subclasses",
        "annotation_exon_boundary_spanning",
        "sequence_context",
        "sequence_kmer_trinucleotide",
        "sequence_repeat_low_mappability_proxy",
        "sequence_96_context",
        "sequence_dna_lm_kmer_delta",
        "utr_rnafold_delta_g_heuristic",
        "qc_support",
        "confidence",
        "summary_score",
    ]
    if has_protein_domains:
        methods.append("protein_domain_annotation")
    if has_protein_structures:
        methods.append("protein_structure_alphafold_esmfold_annotation")
    if has_protein_esm_scores:
        methods.append("protein_esm2_precomputed_delta")
    if has_mirna_sites:
        methods.append("regulatory_mirna_seed_disruption")
    if n_samples <= 0:
        methods.append("cohort_disabled_no_samples")
    elif n_samples == 1:
        methods.append("cohort_single_sample_gt")
    elif n_samples < 10:
        methods.append("cohort_small_sample_carrier_pattern")
    else:
        methods.append("cohort_large_sample_frequency")
        methods.append("cohort_large_sample_maf_binning")
        methods.append("cohort_large_sample_hwe_deviation")
        methods.append("cohort_large_sample_heterozygosity")
        methods.append("cohort_large_sample_inbreeding_coefficient")
        if has_population_windows:
            methods.append("cohort_window_pi_theta")
            methods.append("cohort_window_tajima_d")
            methods.append("cohort_window_ld_haplotype_proxy")
        if has_gene_constraint:
            methods.append("cohort_gene_constraint_proxy")
    if has_sample_info and n_samples >= 2:
        methods.append("cohort_group_subpopulation_af")
        methods.append("cohort_group_private_shared")
        methods.append("cohort_group_fst")
    if has_phenotype and n_samples >= 2:
        methods.append("cohort_case_control_enrichment")
    if has_ml_model:
        methods.append("ml_json_model_inference")
        methods.append("ml_calibration")
        methods.append("ml_uncertainty_ood")
        methods.append("ml_feature_importance")
    return methods


METHOD_SOURCE_TYPES = [
    {
        "field": "DNP_SCORE",
        "source_type": "deterministic_summary",
        "description": "Weighted summary of deterministic annotation, protein, splice, sequence, cohort, and optional ML evidence.",
    },
    {
        "field": "DNP_STRUCT",
        "source_type": "heuristic_proxy",
        "description": "Protein structural-context heuristic from amino-acid properties and local context, not a folded structure model.",
    },
    {
        "field": "DNP_AFSTRUCT",
        "source_type": "external_import",
        "description": "Optional AlphaFold/ESMFold-like per-residue annotation imported from --protein-structures.",
    },
    {
        "field": "DNP_ESM",
        "source_type": "precomputed_import",
        "description": "Optional precomputed ESM-2-like substitution score imported from --protein-lm-scores.",
    },
    {
        "field": "DNP_PROTLM",
        "source_type": "heuristic_proxy",
        "description": "Species-local protein k-mer language-model delta proxy trained from supplied protein/CDS FASTA.",
    },
    {
        "field": "DNP_RNAFOLD",
        "source_type": "heuristic_proxy",
        "description": "Local UTR base-pairing delta-G proxy, not a ViennaRNA minimum-free-energy calculation.",
    },
    {
        "field": "DNP_DNALM",
        "source_type": "reference_trained_proxy",
        "description": "Reference-trained k-mer DNA language-model likelihood delta proxy, not a transformer model.",
    },
    {
        "field": "DNP_ML",
        "source_type": "optional_json_model",
        "description": "Optional portable JSON ML model trained from deterministic DeNovoPath INFO features.",
    },
]


def cohort_gating_note(
    n_samples: int,
    has_population_windows: bool = True,
    has_gene_constraint: bool = True,
) -> str:
    if n_samples <= 0:
        return "No VCF samples were detected; cohort statistics are disabled."
    if n_samples == 1:
        return "Single-sample VCF: only low-weight genotype-state evidence is enabled; population statistics require larger cohorts."
    if n_samples < 10:
        return "Small cohort VCF: carrier-pattern and AC/AN/AF summaries are enabled; HWE, diversity, LD, and gene-constraint proxies require at least 10 samples."
    extras = []
    if has_population_windows:
        extras.append("window diversity/LD/haplotype")
    if has_gene_constraint:
        extras.append("gene-constraint")
    if extras:
        return (
            "Large cohort mode: sample-count-gated population frequency, HWE, heterozygosity, "
            + ", ".join(extras)
            + " proxies are enabled."
        )
    return "Large cohort mode: sample-count-gated population frequency, HWE, and heterozygosity are enabled; optional window and gene-constraint prescans are disabled."


@dataclass
class Transcript:
    tx_id: str
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    cds_segments: List[Tuple[int, int]] = field(default_factory=list)
    utr5_segments: List[Tuple[int, int]] = field(default_factory=list)
    utr3_segments: List[Tuple[int, int]] = field(default_factory=list)

    def transcript_order_cds(self) -> List[Tuple[int, int]]:
        segs = sorted(self.cds_segments)
        if self.strand == "-":
            segs.reverse()
        return segs

    def cds_length(self) -> int:
        return sum(end - start + 1 for start, end in self.cds_segments)

    def transcript_length(self) -> int:
        return self.end - self.start + 1

    def exon_segments(self) -> List[Tuple[int, int]]:
        return sorted(self.cds_segments + self.utr5_segments + self.utr3_segments)

    def contains(self, pos: int) -> bool:
        return self.start <= pos <= self.end

    def overlaps(self, start: int, end: int) -> bool:
        return self.start <= end and self.end >= start

    def region_at(self, pos: int) -> str:
        if any(start <= pos <= end for start, end in self.cds_segments):
            return "CDS"
        if any(start <= pos <= end for start, end in self.utr5_segments):
            return "utr5"
        if any(start <= pos <= end for start, end in self.utr3_segments):
            return "utr3"
        if self.contains(pos):
            return "intron"
        return "intergenic"

    def cds_offset(self, pos: int) -> Optional[int]:
        offset = 0
        for start, end in self.transcript_order_cds():
            seg_len = end - start + 1
            if start <= pos <= end:
                if self.strand == "+":
                    return offset + (pos - start)
                return offset + (end - pos)
            offset += seg_len
        return None

    def min_exon_boundary_distance(self, pos: int) -> Optional[int]:
        distances = []
        for start, end in self.exon_segments():
            distances.append(abs(pos - start))
            distances.append(abs(pos - end))
        return min(distances) if distances else None

    def overlaps_exon_boundary(self, start: int, end: int) -> bool:
        if end < start:
            start, end = end, start
        transcript_start = max(start, self.start)
        transcript_end = min(end, self.end)
        if transcript_end < transcript_start:
            return False
        exon_overlap = 0
        for exon_start, exon_end in self.exon_segments():
            ov_start = max(transcript_start, exon_start)
            ov_end = min(transcript_end, exon_end)
            if ov_end >= ov_start:
                exon_overlap += ov_end - ov_start + 1
        if exon_overlap <= 0:
            return False
        transcript_overlap = transcript_end - transcript_start + 1
        return exon_overlap < transcript_overlap

    def boundary_overlaps_cds(self, start: int, end: int) -> bool:
        if end < start:
            start, end = end, start
        for cds_start, cds_end in self.cds_segments:
            if max(start, cds_start) <= min(end, cds_end):
                return True
        return False


@dataclass
class AltScore:
    score: float
    impact_score: float
    protein_score: float
    grantham_score: float
    blosum_score: float
    codon_usage_score: float
    protein_context_score: float
    protein_structure_score: float
    protein_structure_model_score: float
    protein_esm_score: float
    protein_lm_score: float
    protein_domain_score: float
    splice_score: float
    splice_motif_score: float
    splice_pwm_score: float
    splice_maxent_score: float
    splice_aux_score: float
    splice_ese_score: float
    utr_score: float
    rnafold_score: float
    mirna_score: float
    promoter_score: float
    sequence_score: float
    kmer_score: float
    repeat_score: float
    mutation_context_score: float
    dna_lm_score: float
    cohort_score: float
    hwe_score: float
    heterozygosity_observed: float
    heterozygosity_expected: float
    heterozygosity_deviation_score: float
    inbreeding_coefficient: float
    fst_score: float
    case_control_score: float
    window_pi: float
    window_theta: float
    window_tajima_d: float
    window_ld: float
    window_haplotype: float
    gene_lof_oe: float
    gene_missense_oe: float
    gene_constraint_score: float
    qc_score: float
    confidence_score: float
    consequence: str
    level: str
    ml_score: float = 0.0
    calibrated_score: float = 0.0
    uncertainty_score: float = 0.0
    ood_score: float = 0.0
    feature_importance: str = "."
    gene_id: str = "."
    tx_id: str = "."
    all_transcripts: str = "."
    aa_change: str = "."
    codon_change: str = "."
    hgvs_change: str = "."
    normalized_variant: str = "."
    protein_domain_label: str = "."
    protein_structure_label: str = "."
    protein_esm_label: str = "."
    mirna_label: str = "."
    mutation_context: str = "."
    maf_bin: str = "."
    group_af: str = "."
    private_shared: str = "."
    case_control_af: str = "."
    ac: int = 0
    an: int = 0
    af: float = 0.0
    carriers: int = 0
    n_het: int = 0
    n_hom_alt: int = 0
    n_missing: int = 0


@dataclass
class SampleInfo:
    groups: Dict[str, str]
    phenotypes: Dict[str, str] = field(default_factory=dict)

    def group_indices(self, sample_names: Sequence[str]) -> Dict[str, List[int]]:
        indices: Dict[str, List[int]] = defaultdict(list)
        for idx, sample in enumerate(sample_names):
            group = self.groups.get(sample)
            if group:
                indices[group].append(idx)
        return dict(sorted(indices.items()))

    def phenotype_indices(self, sample_names: Sequence[str]) -> Dict[str, List[int]]:
        indices: Dict[str, List[int]] = {"case": [], "control": []}
        for idx, sample in enumerate(sample_names):
            phenotype = self.phenotypes.get(sample)
            if phenotype in indices:
                indices[phenotype].append(idx)
        return {key: value for key, value in indices.items() if value}

    def has_phenotype(self) -> bool:
        values = set(self.phenotypes.values())
        return "case" in values and "control" in values


@dataclass
class ProteinDomain:
    tx_id: str
    start: int
    end: int
    domain_id: str
    score: float = 0.65

    def contains(self, aa_pos: int) -> bool:
        return self.start <= aa_pos <= self.end

    def distance(self, aa_pos: int) -> int:
        if self.contains(aa_pos):
            return 0
        if aa_pos < self.start:
            return self.start - aa_pos
        return aa_pos - self.end


class ProteinDomainIndex:
    def __init__(self, path: Optional[str] = None) -> None:
        self.by_tx: Dict[str, List[ProteinDomain]] = defaultdict(list)
        if path:
            self._load(path)

    def _load(self, path: str) -> None:
        opener = gzip.open if path.endswith(".gz") else open
        delimiter = "," if path.endswith(".csv") else "\t"
        with opener(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                lowered = {
                    str(key).strip().lower(): str(value).strip()
                    for key, value in row.items()
                    if key is not None and value is not None
                }
                tx_id = first_present(lowered, ("transcript", "tx", "tx_id", "protein", "protein_id", "id"))
                start_text = first_present(lowered, ("start", "aa_start", "begin", "from"))
                end_text = first_present(lowered, ("end", "aa_end", "stop", "to"))
                domain_id = first_present(lowered, ("domain", "domain_id", "name", "label", "description")) or "domain"
                if not tx_id or not start_text or not end_text:
                    continue
                try:
                    start = int(float(start_text))
                    end = int(float(end_text))
                    score_text = first_present(lowered, ("score", "weight", "confidence"))
                    score = float(score_text) if score_text else 0.65
                except ValueError:
                    continue
                if start <= 0 or end <= 0:
                    continue
                if end < start:
                    start, end = end, start
                self.by_tx[tx_id].append(
                    ProteinDomain(
                        tx_id=tx_id,
                        start=start,
                        end=end,
                        domain_id=domain_id.replace(" ", "_"),
                        score=clamp(score),
                    )
                )
        for domains in self.by_tx.values():
            domains.sort(key=lambda item: (item.start, item.end, item.domain_id))

    def __len__(self) -> int:
        return sum(len(items) for items in self.by_tx.values())

    def score(self, tx_id: str, aa_pos: int) -> Tuple[float, str]:
        if aa_pos <= 0:
            return 0.0, "."
        domains = self.by_tx.get(tx_id, [])
        if not domains:
            return 0.0, "."
        best_score = 0.0
        best_label = "."
        for domain in domains:
            distance = domain.distance(aa_pos)
            if distance == 0:
                score = max(0.45, domain.score)
            elif distance <= 5:
                score = max(0.10, domain.score * (1.0 - distance / 6.0) * 0.45)
            else:
                score = 0.0
            if score > best_score:
                best_score = score
                best_label = f"{domain.domain_id}:{domain.start}-{domain.end}"
        return clamp(best_score), best_label


@dataclass
class ProteinStructureFeature:
    tx_id: str
    start: int
    end: int
    plddt: float = 70.0
    rsa: Optional[float] = None
    ss: str = "."
    source: str = "structure"

    def contains(self, aa_pos: int) -> bool:
        return self.start <= aa_pos <= self.end

    def confidence(self) -> float:
        return clamp(self.plddt / 100.0)

    def buriedness(self) -> float:
        if self.rsa is None:
            return 0.50
        return clamp(1.0 - self.rsa)

    def label(self, aa_pos: int) -> str:
        parts = [
            self.source,
            str(aa_pos),
            f"plddt{self.plddt:.0f}",
        ]
        if self.rsa is not None:
            parts.append(f"rsa{self.rsa:.2f}")
        if self.ss and self.ss != ".":
            parts.append(f"ss{self.ss}")
        return ":".join(parts)


class ProteinStructureIndex:
    def __init__(self, path: Optional[str] = None) -> None:
        self.by_tx: Dict[str, List[ProteinStructureFeature]] = defaultdict(list)
        if path:
            self._load(path)

    def _load(self, path: str) -> None:
        opener = gzip.open if path.endswith(".gz") else open
        delimiter = "," if path.endswith(".csv") else "\t"
        with opener(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                lowered = {
                    str(key).strip().lower(): str(value).strip()
                    for key, value in row.items()
                    if key is not None and value is not None
                }
                tx_id = first_present(lowered, ("transcript", "tx", "tx_id", "protein", "protein_id", "id"))
                pos_text = first_present(lowered, ("aa_pos", "position", "pos", "residue", "residue_index"))
                start_text = first_present(lowered, ("start", "aa_start", "begin", "from"))
                end_text = first_present(lowered, ("end", "aa_end", "stop", "to"))
                if pos_text and not start_text:
                    start_text = pos_text
                if pos_text and not end_text:
                    end_text = pos_text
                if not tx_id or not start_text or not end_text:
                    continue
                try:
                    start = int(float(start_text))
                    end = int(float(end_text))
                    plddt_text = first_present(lowered, ("plddt", "confidence", "score", "quality"))
                    plddt = float(plddt_text) if plddt_text else 70.0
                    if plddt <= 1.0:
                        plddt *= 100.0
                    rsa_text = first_present(lowered, ("rsa", "asa_norm", "relative_asa", "relative_sasa"))
                    exposure_text = first_present(lowered, ("exposure", "sasa_norm", "surface_exposure"))
                    rsa = float(rsa_text) if rsa_text else (float(exposure_text) if exposure_text else None)
                except ValueError:
                    continue
                if start <= 0 or end <= 0:
                    continue
                if end < start:
                    start, end = end, start
                ss = first_present(lowered, ("ss", "secstruct", "secondary_structure", "secondary")) or "."
                source = first_present(lowered, ("source", "model", "method")) or "structure"
                self.by_tx[tx_id].append(
                    ProteinStructureFeature(
                        tx_id=tx_id,
                        start=start,
                        end=end,
                        plddt=max(0.0, min(100.0, plddt)),
                        rsa=clamp(rsa) if rsa is not None else None,
                        ss=ss.replace(" ", "_"),
                        source=source.replace(" ", "_"),
                    )
                )
        for records in self.by_tx.values():
            records.sort(key=lambda item: (item.start, item.end, item.source))

    def __len__(self) -> int:
        return sum(len(items) for items in self.by_tx.values())

    def score(self, tx_id: str, aa_pos: int, ref_aa: str, alt_aa: str) -> Tuple[float, str]:
        records = [item for item in self.by_tx.get(tx_id, []) if item.contains(aa_pos)]
        if aa_pos <= 0 or not records:
            return 0.0, "."
        best_score = 0.0
        best_label = "."
        for record in records:
            score = self._score_record(record, ref_aa, alt_aa)
            if score > best_score:
                best_score = score
                best_label = record.label(aa_pos)
        return clamp(best_score), best_label

    @staticmethod
    def _score_record(record: ProteinStructureFeature, ref_aa: str, alt_aa: str) -> float:
        confidence = record.confidence()
        buried = record.buriedness()
        if ref_aa in AA_ORDER and alt_aa in AA_ORDER and ref_aa != alt_aa:
            charge_change = min(1.0, abs(AA_CHARGE[ref_aa] - AA_CHARGE[alt_aa]) / 2.0)
            hydropathy_change = min(1.0, abs(AA_HYDROPATHY[ref_aa] - AA_HYDROPATHY[alt_aa]) / 9.0)
            ref_polar = ref_aa in set("RNDQEHKSTY")
            alt_polar = alt_aa in set("RNDQEHKSTY")
            polarity_change = 1.0 if ref_polar != alt_polar else 0.0
            change = max(charge_change, hydropathy_change, polarity_change)
            special = 0.0
            if "P" in {ref_aa, alt_aa}:
                special = max(special, 0.55)
            if "G" in {ref_aa, alt_aa}:
                special = max(special, 0.35)
            if "C" in {ref_aa, alt_aa}:
                special = max(special, 0.30)
        else:
            change = 0.55
            special = 0.25
        ss = record.ss.upper()
        secondary = 0.0
        if ss.startswith(("H", "E", "B")):
            if alt_aa == "P":
                secondary = 0.75
            elif alt_aa == "G":
                secondary = 0.45
            elif change > 0.6:
                secondary = 0.35
        buried_change = buried * change
        return clamp(confidence * (0.55 * buried_change + 0.25 * secondary + 0.20 * special))


@dataclass
class ProteinLmScoreRecord:
    tx_id: str
    start: int
    end: int
    ref_aa: str = "."
    alt_aa: str = "."
    score: float = 0.0
    label: str = "ESM2"

    def matches(self, aa_pos: int, ref_aa: str, alt_aa: str) -> bool:
        if not (self.start <= aa_pos <= self.end):
            return False
        ref_match = self.ref_aa in {"", "."} or self.ref_aa == ref_aa
        alt_match = self.alt_aa in {"", "."} or self.alt_aa == alt_aa
        return ref_match and alt_match


class ProteinLmScoreIndex:
    def __init__(self, path: Optional[str] = None) -> None:
        self.by_tx: Dict[str, List[ProteinLmScoreRecord]] = defaultdict(list)
        if path:
            self._load(path)

    def _load(self, path: str) -> None:
        opener = gzip.open if path.endswith(".gz") else open
        delimiter = "," if path.endswith(".csv") else "\t"
        with opener(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                lowered = {
                    str(key).strip().lower(): str(value).strip()
                    for key, value in row.items()
                    if key is not None and value is not None
                }
                tx_id = first_present(lowered, ("transcript", "tx", "tx_id", "protein", "protein_id", "id"))
                pos_text = first_present(lowered, ("aa_pos", "position", "pos", "residue", "residue_index"))
                start_text = first_present(lowered, ("start", "aa_start", "begin", "from"))
                end_text = first_present(lowered, ("end", "aa_end", "stop", "to"))
                if pos_text and not start_text:
                    start_text = pos_text
                if pos_text and not end_text:
                    end_text = pos_text
                score_text = first_present(
                    lowered,
                    ("score", "esm", "esm_score", "esm2", "esm2_score", "delta", "delta_score", "pathogenicity"),
                )
                if not tx_id or not start_text or not end_text or not score_text:
                    continue
                try:
                    start = int(float(start_text))
                    end = int(float(end_text))
                    score = float(score_text)
                except ValueError:
                    continue
                if start <= 0 or end <= 0:
                    continue
                if end < start:
                    start, end = end, start
                ref_aa = (first_present(lowered, ("ref_aa", "ref", "reference", "wildtype", "wt")) or ".").upper()
                alt_aa = (first_present(lowered, ("alt_aa", "alt", "alternate", "mutant", "mt")) or ".").upper()
                label = first_present(lowered, ("label", "source", "model", "method", "id")) or "ESM2"
                self.by_tx[tx_id].append(
                    ProteinLmScoreRecord(
                        tx_id=tx_id,
                        start=start,
                        end=end,
                        ref_aa=ref_aa,
                        alt_aa=alt_aa,
                        score=clamp(score),
                        label=label.replace(" ", "_"),
                    )
                )
        for records in self.by_tx.values():
            records.sort(key=lambda item: (item.start, item.end, item.ref_aa, item.alt_aa, item.label))

    def __len__(self) -> int:
        return sum(len(items) for items in self.by_tx.values())

    def score(self, tx_id: str, aa_pos: int, ref_aa: str, alt_aa: str) -> Tuple[float, str]:
        records = self.by_tx.get(tx_id, [])
        if aa_pos <= 0 or not records:
            return 0.0, "."
        best_score = 0.0
        best_label = "."
        for record in records:
            if record.matches(aa_pos, ref_aa, alt_aa) and record.score > best_score:
                best_score = record.score
                allele = f"{record.ref_aa}{aa_pos}{record.alt_aa}" if record.ref_aa not in {"", "."} else str(aa_pos)
                best_label = f"{record.label}:{allele}"
        return clamp(best_score), best_label


def first_present(row: Dict[str, str], aliases: Sequence[str]) -> str:
    for alias in aliases:
        value = row.get(alias, "")
        if value:
            return value
    return ""


@dataclass
class MirnaSite:
    chrom: str
    start: int
    end: int
    mirna_id: str
    tx_id: str = "."
    site_id: str = "."
    score: float = 0.65

    def overlaps(self, start: int, end: int) -> bool:
        return self.start <= end and self.end >= start

    def distance(self, start: int, end: int) -> int:
        if self.overlaps(start, end):
            return 0
        if end < self.start:
            return self.start - end
        return start - self.end


class MirnaSiteIndex:
    def __init__(self, path: Optional[str] = None) -> None:
        self.by_chrom: Dict[str, List[MirnaSite]] = defaultdict(list)
        if path:
            self._load(path)

    def _load(self, path: str) -> None:
        opener = gzip.open if path.endswith(".gz") else open
        delimiter = "," if path.endswith(".csv") else "\t"
        with opener(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                lowered = {
                    str(key).strip().lower(): str(value).strip()
                    for key, value in row.items()
                    if key is not None and value is not None
                }
                chrom = first_present(lowered, ("chrom", "chr", "seqid", "sequence", "contig"))
                start_text = first_present(lowered, ("start", "target_start", "site_start", "begin", "from"))
                end_text = first_present(lowered, ("end", "target_end", "site_end", "stop", "to"))
                if not chrom or not start_text or not end_text:
                    continue
                try:
                    start = int(float(start_text))
                    end = int(float(end_text))
                    score_text = first_present(lowered, ("score", "weight", "confidence", "penalty"))
                    raw_score = float(score_text) if score_text else 0.65
                except ValueError:
                    continue
                if start <= 0 or end <= 0:
                    continue
                if end < start:
                    start, end = end, start
                mirna_id = first_present(lowered, ("mirna", "mirna_id", "mir", "name", "id")) or "miRNA"
                tx_id = first_present(lowered, ("transcript", "tx", "tx_id", "target", "target_id")) or "."
                site_id = first_present(lowered, ("site", "site_id", "label", "seed", "seed_id")) or "."
                score = clamp(raw_score)
                self.by_chrom[chrom].append(
                    MirnaSite(
                        chrom=chrom,
                        start=start,
                        end=end,
                        mirna_id=mirna_id.replace(" ", "_"),
                        tx_id=tx_id,
                        site_id=site_id.replace(" ", "_"),
                        score=score,
                    )
                )
        for sites in self.by_chrom.values():
            sites.sort(key=lambda item: (item.start, item.end, item.mirna_id, item.tx_id))

    def __len__(self) -> int:
        return sum(len(items) for items in self.by_chrom.values())

    def score(self, tx: Transcript, start: int, end: int, region: str) -> Tuple[float, str]:
        if region not in {"utr5", "utr3"}:
            return 0.0, "."
        if end < start:
            start, end = end, start
        candidates = self.by_chrom.get(tx.chrom, [])
        best_score = 0.0
        best_label = "."
        for site in candidates:
            if site.end < start - 3:
                continue
            if site.start > end + 3:
                break
            if site.tx_id not in {".", tx.tx_id, tx.gene_id}:
                continue
            distance = site.distance(start, end)
            if distance == 0:
                score = max(0.45, site.score)
            elif distance <= 3:
                score = max(0.12, site.score * (1.0 - distance / 4.0) * 0.50)
            else:
                continue
            if score > best_score:
                best_score = score
                site_tag = site.site_id if site.site_id != "." else f"{site.start}-{site.end}"
                best_label = f"{site.mirna_id}:{site_tag}"
        return clamp(best_score), best_label


@dataclass
class WindowPopulationStat:
    pi: float = 0.0
    theta: float = 0.0
    tajima_d: float = 0.0
    ld: float = 0.0
    haplotype: float = 0.0
    segregating_sites: int = 0


@dataclass
class WindowAccumulator:
    pi_sum: float = 0.0
    segregating_sites: int = 0
    an_sum: int = 0
    variant_sites: int = 0
    last_dosages: Optional[List[Optional[int]]] = None
    ld_sum: float = 0.0
    ld_max: float = 0.0
    ld_pairs: int = 0

    def add_site(self, counts: List[Dict[str, int]], dosages: Optional[List[Optional[int]]] = None) -> None:
        if not counts:
            return
        an_values = [item["an"] for item in counts if item["an"] > 0]
        if not an_values:
            return
        an = min(an_values)
        if an < 2:
            return
        alt_counts = [max(0, min(item["ac"], an)) for item in counts]
        total_alt = min(sum(alt_counts), an)
        ref_count = max(0, an - total_alt)
        allele_counts = [ref_count] + alt_counts
        pi_contrib = 1.0 - sum((count / an) ** 2 for count in allele_counts if count > 0)
        if pi_contrib <= 0.0:
            return
        self.pi_sum += pi_contrib
        self.segregating_sites += 1
        self.an_sum += an
        self.variant_sites += 1
        if dosages is not None:
            if self.last_dosages is not None:
                r2_value = dosage_r2(self.last_dosages, dosages)
                if r2_value is not None:
                    self.ld_sum += r2_value
                    self.ld_max = max(self.ld_max, r2_value)
                    self.ld_pairs += 1
            self.last_dosages = dosages

    def finalize(self) -> WindowPopulationStat:
        if self.segregating_sites <= 0 or self.variant_sites <= 0:
            return WindowPopulationStat()
        n = max(2, round(self.an_sum / self.variant_sites))
        a1 = sum(1.0 / i for i in range(1, n))
        if a1 <= 0:
            return WindowPopulationStat(pi=self.pi_sum, segregating_sites=self.segregating_sites)
        theta = self.segregating_sites / a1
        tajima_d = 0.0
        if self.segregating_sites >= 2:
            a2 = sum(1.0 / (i * i) for i in range(1, n))
            b1 = (n + 1) / (3 * (n - 1))
            b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
            c1 = b1 - (1 / a1)
            c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
            e1 = c1 / a1
            e2 = c2 / ((a1 * a1) + a2)
            variance = e1 * self.segregating_sites + e2 * self.segregating_sites * (self.segregating_sites - 1)
            if variance > 0:
                tajima_d = (self.pi_sum - theta) / math.sqrt(variance)
        return WindowPopulationStat(
            pi=self.pi_sum,
            theta=theta,
            tajima_d=tajima_d,
            ld=self.ld_max,
            haplotype=clamp(1.0 - (self.ld_sum / self.ld_pairs)) if self.ld_pairs else 0.0,
            segregating_sites=self.segregating_sites,
        )


@dataclass
class GeneConstraintStat:
    lof_oe: float = 0.0
    missense_oe: float = 0.0
    constraint_score: float = 0.0
    observed_lof: int = 0
    observed_missense: int = 0


def clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1.0 / (1.0 + z)
    z = math.exp(value)
    return z / (1.0 + z)


def safe_logit(probability: float) -> float:
    probability = clamp(probability, 1e-6, 1.0 - 1e-6)
    return math.log(probability / (1.0 - probability))


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1].upper()


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs = {}
    for item in attr_text.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attrs[key] = value
    return attrs


def load_fasta_records(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    records: Dict[str, List[str]] = {}
    current = None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line.upper())
    return {key: "".join(value) for key, value in records.items()}


class AnnotationIndex:
    def __init__(self, gff_path: str, bin_size: int = 100_000, promoter_upstream: int = 2000) -> None:
        self.bin_size = bin_size
        self.promoter_upstream = promoter_upstream
        self.transcripts: Dict[str, Transcript] = {}
        self.gene_by_id: Dict[str, Tuple[str, int, int, str]] = {}
        self.tx_bins: Dict[str, Dict[int, List[str]]] = defaultdict(lambda: defaultdict(list))
        self.promoter_bins: Dict[str, Dict[int, List[str]]] = defaultdict(lambda: defaultdict(list))
        self._load(gff_path)
        self._build_bins()

    def _load(self, gff_path: str) -> None:
        pending_children: List[Tuple[str, int, int, str, str]] = []
        opener = gzip.open if gff_path.endswith(".gz") else open
        with opener(gff_path, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _source, feature, start, end, _score, strand, _phase, attrs_text = fields
                start_i, end_i = int(start), int(end)
                attrs = parse_attrs(attrs_text)
                if feature == "gene":
                    gene_id = attrs.get("ID", ".")
                    self.gene_by_id[gene_id] = (chrom, start_i, end_i, strand)
                elif feature in {"mRNA", "transcript"}:
                    tx_id = attrs.get("ID")
                    if not tx_id:
                        continue
                    gene_id = attrs.get("Parent", tx_id.rsplit(".", 1)[0])
                    self.transcripts[tx_id] = Transcript(
                        tx_id=tx_id,
                        gene_id=gene_id,
                        chrom=chrom,
                        start=start_i,
                        end=end_i,
                        strand=strand,
                    )
                elif feature in {"CDS", "five_prime_UTR", "three_prime_UTR"}:
                    parent = attrs.get("Parent")
                    if not parent:
                        continue
                    for tx_id in parent.split(","):
                        pending_children.append((feature, start_i, end_i, chrom, tx_id))

        for _feature, _start, _end, _chrom, tx_id in pending_children:
            if tx_id in self.transcripts or tx_id not in self.gene_by_id:
                continue
            chrom, start, end, strand = self.gene_by_id[tx_id]
            self.transcripts[tx_id] = Transcript(
                tx_id=tx_id,
                gene_id=tx_id,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
            )

        for feature, start, end, _chrom, tx_id in pending_children:
            tx = self.transcripts.get(tx_id)
            if not tx:
                continue
            if feature == "CDS":
                tx.cds_segments.append((start, end))
            elif feature == "five_prime_UTR":
                tx.utr5_segments.append((start, end))
            elif feature == "three_prime_UTR":
                tx.utr3_segments.append((start, end))

    def _build_bins(self) -> None:
        for tx_id, tx in self.transcripts.items():
            first_bin = tx.start // self.bin_size
            last_bin = tx.end // self.bin_size
            for bin_id in range(first_bin, last_bin + 1):
                self.tx_bins[tx.chrom][bin_id].append(tx_id)
            promoter_start, promoter_end = self.promoter_interval(tx)
            first_promoter_bin = promoter_start // self.bin_size
            last_promoter_bin = promoter_end // self.bin_size
            for bin_id in range(first_promoter_bin, last_promoter_bin + 1):
                self.promoter_bins[tx.chrom][bin_id].append(tx_id)

    def query(self, chrom: str, pos: int) -> List[Transcript]:
        bin_id = pos // self.bin_size
        tx_ids = self.tx_bins.get(chrom, {}).get(bin_id, [])
        return [self.transcripts[tx_id] for tx_id in tx_ids if self.transcripts[tx_id].contains(pos)]

    def query_span(self, chrom: str, start: int, end: int) -> List[Transcript]:
        if end < start:
            start, end = end, start
        first_bin = start // self.bin_size
        last_bin = end // self.bin_size
        seen = set()
        hits = []
        for bin_id in range(first_bin, last_bin + 1):
            for tx_id in self.tx_bins.get(chrom, {}).get(bin_id, []):
                if tx_id in seen:
                    continue
                seen.add(tx_id)
                tx = self.transcripts[tx_id]
                if tx.overlaps(start, end):
                    hits.append(tx)
        return hits

    def promoter_interval(self, tx: Transcript) -> Tuple[int, int]:
        if tx.strand == "-":
            return tx.end + 1, tx.end + self.promoter_upstream
        return max(1, tx.start - self.promoter_upstream), tx.start - 1

    def query_promoter(self, chrom: str, pos: int) -> List[Transcript]:
        bin_id = pos // self.bin_size
        tx_ids = self.promoter_bins.get(chrom, {}).get(bin_id, [])
        hits = []
        for tx_id in tx_ids:
            tx = self.transcripts[tx_id]
            start, end = self.promoter_interval(tx)
            if start <= pos <= end:
                hits.append(tx)
        return hits


def translate_codon(codon: str) -> str:
    if len(codon) != 3 or any(base not in "ACGT" for base in codon.upper()):
        return "X"
    return CODON_TABLE.get(codon.upper(), "X")


def aa_property_score(ref_aa: str, alt_aa: str) -> float:
    if ref_aa == alt_aa:
        return 0.02
    if "*" in {ref_aa, alt_aa}:
        return 1.0 if alt_aa == "*" else 0.92
    if ref_aa not in AA_HYDROPATHY or alt_aa not in AA_HYDROPATHY:
        return 0.35
    hydro = abs(AA_HYDROPATHY[ref_aa] - AA_HYDROPATHY[alt_aa]) / 9.0
    charge = min(1.0, abs(AA_CHARGE[ref_aa] - AA_CHARGE[alt_aa]) / 2.0)
    volume = min(1.0, abs(AA_VOLUME[ref_aa] - AA_VOLUME[alt_aa]) / 170.0)
    special = 0.0
    if "P" in {ref_aa, alt_aa}:
        special += 0.12
    if "C" in {ref_aa, alt_aa} and ref_aa != alt_aa:
        special += 0.10
    if "G" in {ref_aa, alt_aa} and ref_aa != alt_aa:
        special += 0.06
    return clamp(0.35 * hydro + 0.30 * charge + 0.25 * volume + special)


def grantham_score(ref_aa: str, alt_aa: str) -> float:
    if ref_aa == alt_aa:
        return 0.0
    key = (ref_aa, alt_aa)
    reverse_key = (alt_aa, ref_aa)
    distance = GRANTHAM_DISTANCE.get(key, GRANTHAM_DISTANCE.get(reverse_key))
    if distance is None:
        return 0.0
    return clamp(distance / 215.0)


def blosum62_score(ref_aa: str, alt_aa: str) -> float:
    if ref_aa == alt_aa:
        return 0.0
    raw = BLOSUM62.get((ref_aa, alt_aa))
    if raw is None:
        return 0.0
    return clamp((4.0 - raw) / 8.0)


def weighted_mean(scores: Dict[str, float], weights: Dict[str, float]) -> float:
    weight_sum = sum(weights.get(key, 0.0) for key in scores)
    if weight_sum <= 0:
        return 0.0
    return clamp(sum(scores[key] * weights.get(key, 0.0) for key in scores) / weight_sum)


def build_codon_usage(
    cds_records: Dict[str, str],
    max_training_codons: int = DEFAULT_MAX_CODON_USAGE_TRAINING_CODONS,
) -> Dict[str, float]:
    aa_codon_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    remaining = max(0, int(max_training_codons))
    if not cds_records or remaining <= 0:
        return {}
    for seq in cds_records.values():
        seq = seq.upper()
        for idx in range(0, len(seq) - 2, 3):
            if remaining <= 0:
                break
            codon = seq[idx : idx + 3]
            if len(codon) != 3 or any(base not in "ACGT" for base in codon):
                continue
            aa = translate_codon(codon)
            if aa == "X":
                continue
            aa_codon_counts[aa][codon] += 1
            remaining -= 1
        if remaining <= 0:
            break
    usage = {}
    for aa, counts in aa_codon_counts.items():
        total = sum(counts.values())
        if total <= 0:
            continue
        for codon, count in counts.items():
            usage[codon] = count / total
    return usage


def codon_usage_bias_score(ref_codon: str, alt_codon: str, usage: Dict[str, float]) -> float:
    ref_codon = ref_codon.upper()
    alt_codon = alt_codon.upper()
    if ref_codon == alt_codon:
        return 0.0
    ref_aa = translate_codon(ref_codon)
    alt_aa = translate_codon(alt_codon)
    if "X" in {ref_aa, alt_aa} or "*" in {ref_aa, alt_aa}:
        return 0.0
    ref_rel = usage.get(ref_codon, 0.0)
    alt_rel = usage.get(alt_codon, 0.0)
    if alt_rel <= 0:
        alt_rel = 0.05
    if ref_aa == alt_aa:
        rarity = 1.0 - alt_rel
        directional_loss = max(0.0, ref_rel - alt_rel)
        return clamp(0.70 * rarity + 0.30 * directional_loss)
    return clamp(0.30 * (1.0 - alt_rel))


def trinucleotide_kmer_score(seq: str, variant_offset: int, ref: str, alt: str) -> float:
    seq = seq.upper()
    ref = ref.upper()
    alt = alt.upper()
    if not seq or variant_offset < 0 or variant_offset >= len(seq):
        return 0.0
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        return clamp(0.08 + 0.08 * min(1.0, abs(len(alt) - len(ref)) / 10.0))
    if variant_offset == 0 or variant_offset >= len(seq) - 1:
        return 0.0
    ref_tri = seq[variant_offset - 1 : variant_offset + 2]
    if len(ref_tri) != 3 or any(base not in "ACGT" for base in ref_tri):
        return 0.0
    mutated = seq[:variant_offset] + alt + seq[variant_offset + 1 :]
    alt_tri = mutated[variant_offset - 1 : variant_offset + 2]
    kmers = [seq[idx : idx + 3] for idx in range(0, len(seq) - 2)]
    valid_kmers = [kmer for kmer in kmers if all(base in "ACGT" for base in kmer)]
    if not valid_kmers:
        return 0.0
    total = len(valid_kmers)
    ref_freq = valid_kmers.count(ref_tri) / total
    alt_freq = valid_kmers.count(alt_tri) / total
    alt_rarity = 1.0 - alt_freq
    novelty = 1.0 if alt_freq == 0 else 0.0
    directional_change = clamp(abs(ref_freq - alt_freq) * 8.0)
    cpg_change = 1.0 if ("CG" in ref_tri) != ("CG" in alt_tri) else 0.0
    return clamp(0.35 * alt_rarity + 0.30 * novelty + 0.20 * directional_change + 0.15 * cpg_change)


def repeat_low_mappability_proxy_score(seq: str) -> float:
    seq = "".join(base for base in seq.upper() if base in "ACGT")
    if not seq:
        return 0.0
    entropy_score = 1.0 - shannon_entropy(seq)
    homopolymer_score = min(1.0, max(0, max_homopolymer(seq) - 5) / 10)
    tandem_score = tandem_repeat_score(seq)
    return clamp(0.35 * tandem_score + 0.35 * entropy_score + 0.30 * homopolymer_score)


def mutation_96_context(seq: str, variant_offset: int, ref: str, alt: str) -> Tuple[str, float]:
    seq = seq.upper()
    ref = ref.upper()
    alt = alt.upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        return ".", 0.05
    if variant_offset <= 0 or variant_offset >= len(seq) - 1:
        return ".", 0.0
    tri = seq[variant_offset - 1 : variant_offset + 2]
    if len(tri) != 3 or any(base not in "ACGT" for base in tri):
        return ".", 0.0
    if tri[1] != ref:
        tri = tri[0] + ref + tri[2]
    context = tri
    ref_base = ref
    alt_base = alt
    if ref_base in {"A", "G"}:
        context = revcomp(tri)
        ref_base = ref_base.translate(DNA_COMP).upper()
        alt_base = alt_base.translate(DNA_COMP).upper()
    label = f"{context[0]}[{ref_base}>{alt_base}]{context[2]}"
    score = mutation_context_score(context, ref_base, alt_base)
    return label, score


def mutation_context_score(context: str, ref: str, alt: str) -> float:
    if ref == alt or ref not in "CT" or alt not in "ACGT":
        return 0.0
    substitution = f"{ref}>{alt}"
    transition = substitution in {"C>T", "T>C"}
    cpg = context[1:3] == "CG" if ref == "C" else context[0:2] == "CA"
    if cpg and substitution == "C>T":
        return 0.75
    if cpg:
        return 0.50
    if substitution in {"C>A", "C>G"} and context[2] == "G":
        return 0.45
    if transition:
        return 0.25
    return 0.35


def rna_pairing_energy_proxy(seq: str) -> float:
    rna = seq.upper().replace("T", "U")
    rna = "".join(base for base in rna if base in "ACGU")
    n = len(rna)
    if n < 12:
        return 0.0
    pair_energy = {
        ("G", "C"): -3.0,
        ("C", "G"): -3.0,
        ("A", "U"): -2.0,
        ("U", "A"): -2.0,
        ("G", "U"): -1.0,
        ("U", "G"): -1.0,
    }
    weighted_pairs = 0.0
    for left in range(n - 4):
        left_base = rna[left]
        for right in range(left + 4, n):
            energy = pair_energy.get((left_base, rna[right]))
            if energy is None:
                continue
            span = right - left
            if span < 6:
                distance_weight = 0.45
            elif span <= 35:
                distance_weight = 1.0
            else:
                distance_weight = max(0.20, 1.0 - (span - 35) / 80.0)
            weighted_pairs += abs(energy) * distance_weight
    gc_fraction = (rna.count("G") + rna.count("C")) / n
    gc_adjustment = 0.85 + 0.30 * gc_fraction
    return -(weighted_pairs / max(1.0, n / 2.0)) * gc_adjustment


class DnaKmerLanguageModel:
    def __init__(self, reference: Fasta, k: int = 5, max_training_bases: int = 250_000) -> None:
        self.k = max(2, int(k))
        self.pseudo = 0.25
        self.context_counts: Dict[str, Counter[str]] = defaultdict(Counter)
        self._train(reference, max_training_bases)

    def _train(self, reference: Fasta, max_training_bases: int) -> None:
        remaining = max(0, int(max_training_bases))
        if remaining <= 0:
            return
        try:
            names = list(reference.keys())
        except Exception:
            names = []
        if not names:
            return
        per_contig = max(2000, remaining // max(1, min(len(names), 80)))
        for name in names:
            if remaining <= 0:
                break
            try:
                chrom_len = len(reference[name])
                take = min(chrom_len, per_contig, remaining)
                if take <= 0:
                    continue
                seq = reference[name][0:take].seq.upper()
            except Exception:
                continue
            self.add_sequence(seq)
            remaining -= take

    def add_sequence(self, seq: str) -> None:
        seq = "".join(base for base in seq.upper() if base in "ACGT")
        if len(seq) < self.k:
            return
        context_len = self.k - 1
        for idx in range(0, len(seq) - context_len):
            context = seq[idx : idx + context_len]
            base = seq[idx + context_len]
            self.context_counts[context][base] += 1

    def avg_neg_log_prob(self, seq: str) -> float:
        seq = "".join(base for base in seq.upper() if base in "ACGT")
        if len(seq) < self.k:
            return 0.0
        context_len = self.k - 1
        total_score = 0.0
        n_terms = 0
        for idx in range(0, len(seq) - context_len):
            context = seq[idx : idx + context_len]
            base = seq[idx + context_len]
            total_score += self.neg_log_prob(context, base)
            n_terms += 1
        return total_score / n_terms if n_terms else 0.0

    def neg_log_prob(self, context: str, base: str) -> float:
        counts = self.context_counts.get(context)
        if counts:
            denom = sum(counts.values()) + 4.0 * self.pseudo
            prob = (counts.get(base, 0.0) + self.pseudo) / denom
        else:
            prob = 0.25
        return -math.log(prob)

    def delta_score(self, seq: str, variant_offset: int, ref: str, alt: str) -> float:
        if not seq or variant_offset < 0 or variant_offset > len(seq):
            return 0.0
        ref = ref.upper()
        alt = alt.upper()
        if any(base not in "ACGT" for base in ref + alt):
            return 0.0
        ref_end = min(len(seq), variant_offset + len(ref))
        mutated = seq[:variant_offset] + alt + seq[ref_end:]
        if len(mutated) < self.k:
            return 0.0
        if len(ref) == len(alt) and all(base in "ACGT" for base in seq):
            context_len = self.k - 1
            n_terms = max(1, len(seq) - context_len)
            var_end = variant_offset + len(ref) - 1
            first_term = max(0, variant_offset - context_len)
            last_term = min(n_terms - 1, var_end)
            delta_sum = 0.0
            for idx in range(first_term, last_term + 1):
                ref_context = seq[idx : idx + context_len]
                ref_base = seq[idx + context_len]
                alt_context = mutated[idx : idx + context_len]
                alt_base = mutated[idx + context_len]
                delta_sum += self.neg_log_prob(alt_context, alt_base) - self.neg_log_prob(ref_context, ref_base)
            delta = delta_sum / n_terms
        else:
            ref_nll = self.avg_neg_log_prob(seq)
            alt_nll = self.avg_neg_log_prob(mutated)
            delta = alt_nll - ref_nll
        increase = max(0.0, delta)
        absolute = abs(delta)
        length_boost = min(0.15, abs(len(alt) - len(ref)) / 30.0)
        return clamp(0.75 * min(1.0, increase / 0.45) + 0.25 * min(1.0, absolute / 0.60) + length_boost)


class ProteinKmerLanguageModel:
    def __init__(
        self,
        records: Dict[str, str],
        k: int = 3,
        max_training_residues: int = DEFAULT_MAX_PROTEIN_LM_TRAINING_RESIDUES,
    ) -> None:
        self.k = k
        self.alphabet = set(AA_ORDER)
        self.prefix_counts: Dict[str, Counter[str]] = defaultdict(Counter)
        self.prefix_totals: Dict[str, int] = defaultdict(int)
        self.background = Counter()
        self.background_total = 0
        remaining = max(0, int(max_training_residues))
        for seq in records.values():
            if remaining <= 0:
                break
            used = self._add_sequence(seq, remaining)
            remaining -= used

    def _clean(self, seq: str) -> str:
        return "".join(aa for aa in seq.upper() if aa in self.alphabet)

    def _add_sequence(self, seq: str, max_residues: Optional[int] = None) -> int:
        seq = self._clean(seq)
        if max_residues is not None:
            seq = seq[: max(0, int(max_residues))]
        for aa in seq:
            self.background[aa] += 1
            self.background_total += 1
        if len(seq) <= self.k:
            return len(seq)
        for idx in range(self.k - 1, len(seq)):
            prefix = seq[idx - self.k + 1 : idx]
            aa = seq[idx]
            self.prefix_counts[prefix][aa] += 1
            self.prefix_totals[prefix] += 1
        return len(seq)

    def __len__(self) -> int:
        return len(self.prefix_totals)

    def avg_neg_log_prob(self, seq: str) -> float:
        seq = self._clean(seq)
        if not seq:
            return 0.0
        vocab = len(AA_ORDER)
        values = []
        if len(seq) < self.k or not self.prefix_totals:
            denom = self.background_total + vocab
            for aa in seq:
                prob = (self.background[aa] + 1.0) / denom if denom else 1.0 / vocab
                values.append(-math.log(prob))
            return sum(values) / len(values) if values else 0.0
        for idx in range(self.k - 1, len(seq)):
            prefix = seq[idx - self.k + 1 : idx]
            aa = seq[idx]
            denom = self.prefix_totals.get(prefix, 0) + vocab
            count = self.prefix_counts.get(prefix, Counter()).get(aa, 0)
            prob = (count + 1.0) / denom
            values.append(-math.log(prob))
        return sum(values) / len(values) if values else 0.0

    def delta_score(self, protein: str, aa_pos: int, ref_aa: str, alt_aa: str, window: int = 24) -> float:
        if not protein or aa_pos < 1 or ref_aa == alt_aa:
            return 0.0
        if ref_aa not in self.alphabet or alt_aa not in self.alphabet:
            return 0.0
        protein = protein.upper()
        center = aa_pos - 1
        if center >= len(protein):
            return 0.0
        start = max(0, center - window)
        end = min(len(protein), center + window + 1)
        local_offset = center - start
        seq = protein[start:end]
        if local_offset >= len(seq) or seq[local_offset] != ref_aa:
            seq = seq[:local_offset] + ref_aa + seq[local_offset + 1 :]
        mutated = seq[:local_offset] + alt_aa + seq[local_offset + 1 :]
        ref_nll = self.avg_neg_log_prob(seq)
        alt_nll = self.avg_neg_log_prob(mutated)
        increase = max(0.0, alt_nll - ref_nll)
        absolute = abs(alt_nll - ref_nll)
        return clamp(0.70 * min(1.0, increase / 0.75) + 0.30 * min(1.0, absolute / 1.10))


class MlModel:
    def __init__(self, path: Optional[str] = None) -> None:
        self.path = path
        self.loaded = False
        self.model_type = "logistic_json"
        self.features = list(ML_FEATURES)
        self.coef: List[float] = []
        self.intercept = 0.0
        self.learning_rate = 1.0
        self.trees: List[Dict[str, List[float]]] = []
        self.mean: List[float] = []
        self.scale: List[float] = []
        self.calibration_scale = 1.0
        self.calibration_shift = 0.0
        self.feature_importance: List[Tuple[str, float]] = []
        if path:
            self._load(path)

    def _load(self, path: str) -> None:
        with open(path) as handle:
            data = json.load(handle)
        self.model_type = str(data.get("model_type", "logistic_json"))
        features = data.get("features", ML_FEATURES)
        if not isinstance(features, list):
            raise ValueError("ML model JSON must contain list field 'features'")
        self.features = [str(item) for item in features]
        self.intercept = float(data.get("intercept", 0.0))
        if self.model_type == "sklearn_gradient_boosting_json":
            trees = data.get("trees")
            if not isinstance(trees, list) or not trees:
                raise ValueError("GBM model JSON must contain non-empty list field 'trees'")
            self.learning_rate = float(data.get("learning_rate", 1.0))
            self.trees = []
            for tree in trees:
                if not isinstance(tree, dict):
                    raise ValueError("Each GBM tree must be a JSON object")
                required = ("children_left", "children_right", "feature", "threshold", "value")
                if any(key not in tree for key in required):
                    raise ValueError("Each GBM tree must contain children_left, children_right, feature, threshold, and value")
                self.trees.append(
                    {
                        "children_left": [int(item) for item in tree["children_left"]],
                        "children_right": [int(item) for item in tree["children_right"]],
                        "feature": [int(item) for item in tree["feature"]],
                        "threshold": [float(item) for item in tree["threshold"]],
                        "value": [float(item) for item in tree["value"]],
                    }
                )
            self.coef = []
        else:
            coef = data.get("coef")
            if not isinstance(coef, list):
                raise ValueError("Logistic ML model JSON must contain list field 'coef'")
            if len(self.features) != len(coef):
                raise ValueError("ML model feature and coefficient lengths differ")
            self.coef = [float(item) for item in coef]
        self.mean = [float(item) for item in data.get("mean", [0.0] * len(self.features))]
        self.scale = [max(1e-6, float(item)) for item in data.get("scale", [1.0] * len(self.features))]
        if len(self.mean) != len(self.features) or len(self.scale) != len(self.features):
            raise ValueError("ML model mean/scale lengths must match features")
        calibration = data.get("calibration", {}) if isinstance(data.get("calibration", {}), dict) else {}
        self.calibration_scale = float(calibration.get("scale", 1.0))
        self.calibration_shift = float(calibration.get("shift", 0.0))
        importance = data.get("feature_importance")
        if isinstance(importance, dict):
            self.feature_importance = sorted(
                ((str(key), float(value)) for key, value in importance.items()),
                key=lambda item: item[1],
                reverse=True,
            )
        else:
            default_values = [abs(value) for value in self.coef] if self.coef else [1.0] * len(self.features)
            self.feature_importance = sorted(
                zip(self.features, default_values),
                key=lambda item: item[1],
                reverse=True,
            )
        self.loaded = True

    def __len__(self) -> int:
        return 1 if self.loaded else 0

    def predict(self, feature_values: Dict[str, float]) -> Tuple[float, float, float, float]:
        if not self.loaded:
            return 0.0, 0.0, 0.0, 0.0
        z = self.intercept
        abs_zscores = []
        zscores = []
        for idx, feature in enumerate(self.features):
            raw = float(feature_values.get(feature, 0.0))
            zscore = (raw - self.mean[idx]) / self.scale[idx]
            zscores.append(zscore)
            if self.model_type != "sklearn_gradient_boosting_json":
                z += self.coef[idx] * zscore
            abs_zscores.append(abs(zscore))
        if self.model_type == "sklearn_gradient_boosting_json":
            z = self.intercept + self.learning_rate * sum(self._tree_value(tree, zscores) for tree in self.trees)
        ml_score = sigmoid(z)
        logit_value = safe_logit(ml_score)
        calibrated = sigmoid(self.calibration_scale * logit_value + self.calibration_shift)
        mean_abs = sum(abs_zscores) / len(abs_zscores) if abs_zscores else 0.0
        max_abs = max(abs_zscores) if abs_zscores else 0.0
        ood = clamp(0.65 * min(1.0, mean_abs / 4.0) + 0.35 * min(1.0, max_abs / 8.0))
        uncertainty = clamp((1.0 - abs(calibrated - 0.5) * 2.0) * 0.75 + ood * 0.25)
        return clamp(ml_score), clamp(calibrated), uncertainty, ood

    @staticmethod
    def _tree_value(tree: Dict[str, List[float]], values: Sequence[float]) -> float:
        node = 0
        children_left = tree["children_left"]
        children_right = tree["children_right"]
        features = tree["feature"]
        thresholds = tree["threshold"]
        outputs = tree["value"]
        while 0 <= node < len(outputs):
            left = int(children_left[node])
            right = int(children_right[node])
            if left < 0 and right < 0:
                return float(outputs[node])
            feature_idx = int(features[node])
            observed = values[feature_idx] if 0 <= feature_idx < len(values) else 0.0
            node = left if observed <= float(thresholds[node]) else right
        return 0.0

    def format_feature_importance(self, feature_values: Dict[str, float], top_n: int = 5) -> str:
        if not self.loaded:
            return "."
        if self.model_type == "sklearn_gradient_boosting_json":
            ranked = self.feature_importance[:top_n]
            parts = [f"{feature}:{importance:.3f}" for feature, importance in ranked]
            return "|".join(parts) if parts else "."
        contributions = []
        for idx, feature in enumerate(self.features):
            raw = float(feature_values.get(feature, 0.0))
            zscore = (raw - self.mean[idx]) / self.scale[idx]
            contribution = self.coef[idx] * zscore
            contributions.append((feature, contribution, abs(contribution)))
        if not contributions:
            return "."
        ranked = sorted(contributions, key=lambda item: item[2], reverse=True)[:top_n]
        parts = [f"{feature}:{contribution:+.3f}" for feature, contribution, _abs_value in ranked]
        return "|".join(parts) if parts else "."


def normalize_allele(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    ref = ref.upper()
    alt = alt.upper()
    while ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
        ref = ref[:-1]
        alt = alt[:-1]
    while ref and alt and ref[0] == alt[0] and (len(ref) > 1 or len(alt) > 1):
        pos += 1
        ref = ref[1:]
        alt = alt[1:]
    return pos, ref, alt


def normalized_variant_label(pos: int, ref: str, alt: str) -> str:
    norm_pos, norm_ref, norm_alt = normalize_allele(pos, ref, alt)
    return f"{norm_pos}:{norm_ref or '-'}>{norm_alt or '-'}"


def changed_ref_span(pos: int, ref: str, alt: str) -> Tuple[int, int]:
    norm_pos, norm_ref, _norm_alt = normalize_allele(pos, ref, alt)
    if norm_ref:
        return norm_pos, norm_pos + len(norm_ref) - 1
    return norm_pos, norm_pos


def tandem_repeat_score(seq: str, max_unit: int = 6) -> float:
    seq = "".join(base for base in seq.upper() if base in "ACGT")
    if len(seq) < 6:
        return 0.0
    best_span = 0
    for unit_len in range(1, max_unit + 1):
        if len(seq) < unit_len * 3:
            continue
        pattern = TANDEM_REPEAT_PATTERNS.get(unit_len)
        if not pattern:
            continue
        for match in pattern.finditer(seq):
            best_span = max(best_span, match.end() - match.start())
    return clamp(max(0, best_span - 6) / 30)


def shannon_entropy(seq: str) -> float:
    seq = [base for base in seq.upper() if base in "ACGT"]
    if not seq:
        return 0.0
    counts = Counter(seq)
    total = len(seq)
    entropy = -sum((count / total) * math.log2(count / total) for count in counts.values())
    return entropy / 2.0


def max_homopolymer(seq: str) -> int:
    best = 0
    current = 0
    last = None
    for base in seq.upper():
        if base == last and base in "ACGT":
            current += 1
        else:
            current = 1
            last = base
        best = max(best, current)
    return best


def normalized_aa_entropy(seq: str) -> float:
    seq = [aa for aa in seq.upper() if aa in AA_ORDER]
    if not seq:
        return 0.0
    counts = Counter(seq)
    total = len(seq)
    entropy = -sum((count / total) * math.log2(count / total) for count in counts.values())
    return entropy / math.log2(20)


def max_residue_run(seq: str) -> int:
    best = 0
    current = 0
    last = None
    for aa in seq.upper():
        if aa == last and aa in AA_ORDER:
            current += 1
        else:
            current = 1
            last = aa
        best = max(best, current)
    return best


def protein_context_score(protein: str, aa_pos: int, window: int = 12) -> float:
    if not protein or aa_pos < 1:
        return 0.0
    center = aa_pos - 1
    if center >= len(protein):
        return 0.0
    start = max(0, center - window)
    end = min(len(protein), center + window + 1)
    seq = protein[start:end].upper()
    residues = [aa for aa in seq if aa in AA_ORDER]
    if not residues:
        return 0.0
    low_complexity = 1.0 - normalized_aa_entropy("".join(residues))
    run_score = clamp(max(0, max_residue_run("".join(residues)) - 3) / 8.0)
    disorder_residues = set("PGSQEKR")
    order_residues = set("WFYILVMC")
    disorder_fraction = sum(1 for aa in residues if aa in disorder_residues) / len(residues)
    order_fraction = sum(1 for aa in residues if aa in order_residues) / len(residues)
    disorder_bias = clamp(disorder_fraction - 0.35 + max(0.0, 0.30 - order_fraction))
    return clamp(0.40 * low_complexity + 0.25 * run_score + 0.35 * disorder_bias)


def protein_structure_score(protein: str, aa_pos: int, ref_aa: str, alt_aa: str, window: int = 7) -> float:
    if not protein or aa_pos < 1 or ref_aa == alt_aa:
        return 0.0
    center = aa_pos - 1
    if center >= len(protein) or ref_aa not in AA_ORDER or alt_aa not in AA_ORDER:
        return 0.0
    start = max(0, center - window)
    end = min(len(protein), center + window + 1)
    residues = [aa for aa in protein[start:end].upper() if aa in AA_ORDER]
    if not residues:
        return 0.0

    mean_hydropathy = sum(AA_HYDROPATHY[aa] for aa in residues) / len(residues)
    buried_score = clamp((mean_hydropathy + 1.5) / 5.5)
    surface_score = 1.0 - buried_score

    helix_pref = set("ALMEQKRH")
    sheet_pref = set("VIFYWT")
    turn_breakers = set("PGNSD")
    ref_helix = 1.0 if ref_aa in helix_pref else 0.0
    alt_helix = 1.0 if alt_aa in helix_pref else 0.0
    ref_sheet = 1.0 if ref_aa in sheet_pref else 0.0
    alt_sheet = 1.0 if alt_aa in sheet_pref else 0.0
    secondary_shift = clamp(0.5 * abs(ref_helix - alt_helix) + 0.5 * abs(ref_sheet - alt_sheet))
    breaker_gain = 1.0 if alt_aa in turn_breakers and ref_aa not in turn_breakers else 0.0
    breaker_loss = 0.5 if ref_aa in turn_breakers and alt_aa not in turn_breakers else 0.0

    charge_change = min(1.0, abs(AA_CHARGE[ref_aa] - AA_CHARGE[alt_aa]) / 2.0)
    ref_polar = ref_aa in set("RNDQEHKSTY")
    alt_polar = alt_aa in set("RNDQEHKSTY")
    polarity_change = 1.0 if ref_polar != alt_polar else 0.0
    buried_charged_risk = buried_score * max(charge_change, polarity_change)

    gly_pro_cys_risk = 0.0
    if "P" in {ref_aa, alt_aa} and ref_aa != alt_aa:
        gly_pro_cys_risk += 0.30
    if "G" in {ref_aa, alt_aa} and ref_aa != alt_aa:
        gly_pro_cys_risk += 0.18
    if "C" in {ref_aa, alt_aa} and ref_aa != alt_aa:
        gly_pro_cys_risk += 0.22

    exposure_penalty = 0.20 * surface_score * charge_change
    return clamp(
        0.42 * buried_charged_risk
        + 0.24 * secondary_shift
        + 0.18 * breaker_gain
        + 0.08 * breaker_loss
        + 0.18 * gly_pro_cys_risk
        + exposure_penalty
    )


def genotype_counts(format_keys: List[str], sample_values: List[str], alt_count: int) -> List[Dict[str, int]]:
    counts = [
        {
            "ac": 0,
            "an": 0,
            "carriers": 0,
            "n_het": 0,
            "n_hom_alt": 0,
            "n_hom_ref": 0,
            "n_called_diploid": 0,
            "n_missing": 0,
        }
        for _ in range(alt_count)
    ]
    try:
        gt_idx = format_keys.index("GT")
    except ValueError:
        return counts

    for sample_text in sample_values:
        values = sample_text.split(":")
        gt = values[gt_idx] if gt_idx < len(values) else "."
        if gt in {".", "./.", ".|."}:
            for item in counts:
                item["n_missing"] += 1
            continue
        alleles = re.split(r"[/|]", gt)
        called = [int(a) for a in alleles if a.isdigit()]
        if not called:
            for item in counts:
                item["n_missing"] += 1
            continue
        for item in counts:
            item["an"] += len(called)
        for alt_idx in range(1, alt_count + 1):
            n_alt = sum(1 for allele in called if allele == alt_idx)
            item = counts[alt_idx - 1]
            item["ac"] += n_alt
            is_ref_current_alt_diploid = len(called) == 2 and all(allele in {0, alt_idx} for allele in called)
            if is_ref_current_alt_diploid:
                item["n_called_diploid"] += 1
            if n_alt:
                item["carriers"] += 1
            if is_ref_current_alt_diploid and n_alt == 1:
                item["n_het"] += 1
            if is_ref_current_alt_diploid and n_alt == len(called):
                item["n_hom_alt"] += 1
            if is_ref_current_alt_diploid and n_alt == 0:
                item["n_hom_ref"] += 1
    return counts


def genotype_alt_dosages(format_keys: List[str], sample_values: List[str], alt_idx: int = 1) -> List[Optional[int]]:
    try:
        gt_idx = format_keys.index("GT")
    except ValueError:
        return [None for _ in sample_values]
    dosages: List[Optional[int]] = []
    for sample_text in sample_values:
        values = sample_text.split(":")
        gt = values[gt_idx] if gt_idx < len(values) else "."
        if gt in {".", "./.", ".|."}:
            dosages.append(None)
            continue
        alleles = re.split(r"[/|]", gt)
        called = [int(a) for a in alleles if a.isdigit()]
        if len(called) != 2:
            dosages.append(None)
            continue
        dosages.append(sum(1 for allele in called if allele == alt_idx))
    return dosages


def dosage_r2(left: Sequence[Optional[int]], right: Sequence[Optional[int]]) -> Optional[float]:
    pairs = [(float(a), float(b)) for a, b in zip(left, right) if a is not None and b is not None]
    if len(pairs) < 4:
        return None
    xs = [item[0] for item in pairs]
    ys = [item[1] for item in pairs]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    if var_x <= 0.0 or var_y <= 0.0:
        return None
    cov = sum((x - mean_x) * (y - mean_y) for x, y in pairs)
    return clamp((cov * cov) / (var_x * var_y))


def normalize_sample_info_key(row: Dict[str, str], aliases: Sequence[str]) -> str:
    lowered = {key.strip().lower(): value.strip() for key, value in row.items() if key is not None and value is not None}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    return ""


def normalize_phenotype(value: str) -> str:
    text = value.strip().lower()
    if text in {"case", "affected", "disease", "trait", "1", "yes", "true"}:
        return "case"
    if text in {"control", "unaffected", "normal", "0", "no", "false"}:
        return "control"
    return ""


def load_sample_info(path: Optional[str]) -> Optional[SampleInfo]:
    if not path:
        return None
    delimiter = "," if path.endswith(".csv") else "\t"
    groups: Dict[str, str] = {}
    phenotypes: Dict[str, str] = {}
    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            sample = normalize_sample_info_key(row, ("sample", "sample_id", "sampleid", "id"))
            group = normalize_sample_info_key(row, ("group", "population", "pop", "subpopulation", "subpop"))
            if not sample or not group:
                continue
            groups[sample] = group.replace(" ", "_")
            phenotype_text = normalize_sample_info_key(row, ("phenotype", "status", "case_control", "casecontrol", "trait"))
            phenotype = normalize_phenotype(phenotype_text)
            if phenotype:
                phenotypes[sample] = phenotype
    if not groups:
        raise ValueError("sample_info must contain at least sample and group columns with one matching row")
    return SampleInfo(groups=groups, phenotypes=phenotypes)


def group_population_scores(
    format_keys: List[str],
    sample_values: List[str],
    alt_count: int,
    sample_names: Sequence[str],
    sample_info: Optional[SampleInfo],
) -> List[Dict[str, object]]:
    empty = [
        {
            "group_af": ".",
            "private_shared": ".",
            "fst_score": 0.0,
            "case_control_af": ".",
            "case_control_score": 0.0,
        }
        for _ in range(alt_count)
    ]
    if not sample_info or not sample_names or not sample_values:
        return empty

    group_indices = sample_info.group_indices(sample_names)
    group_counts = {
        group: genotype_counts(format_keys, [sample_values[idx] for idx in indices], alt_count)
        for group, indices in group_indices.items()
        if indices
    }
    phenotype_indices = sample_info.phenotype_indices(sample_names)
    phenotype_counts = {
        phenotype: genotype_counts(format_keys, [sample_values[idx] for idx in indices], alt_count)
        for phenotype, indices in phenotype_indices.items()
        if indices
    }

    results = []
    for alt_idx in range(alt_count):
        group_af_parts = []
        present_groups = []
        group_an_total = 0
        group_ac_total = 0
        heterozygosity_sum = 0.0
        private_shared = "."
        if len(group_counts) >= 2:
            for group, counts_by_alt in group_counts.items():
                count = counts_by_alt[alt_idx]
                an = count["an"]
                ac = count["ac"]
                af = ac / an if an else 0.0
                group_af_parts.append(f"{group}:{af:.4f}")
                if ac > 0:
                    present_groups.append(group)
                if an:
                    group_an_total += an
                    group_ac_total += ac
                    heterozygosity_sum += an * 2.0 * af * (1.0 - af)

        if len(group_counts) >= 2 and not present_groups:
            private_shared = "absent"
        elif len(group_counts) >= 2 and len(present_groups) == 1:
            private_shared = f"private:{present_groups[0]}"
        elif len(group_counts) >= 2 and len(present_groups) == len(group_counts):
            private_shared = "shared_all"
        elif len(group_counts) >= 2:
            private_shared = "shared:" + "|".join(present_groups)
        else:
            private_shared = "."

        fst_score = 0.0
        if len(group_counts) >= 2 and group_an_total and 0 < group_ac_total < group_an_total:
            pbar = group_ac_total / group_an_total
            ht = 2.0 * pbar * (1.0 - pbar)
            hs = heterozygosity_sum / group_an_total if group_an_total else 0.0
            fst_score = clamp((ht - hs) / ht) if ht > 0 else 0.0

        case_control_af = "."
        case_control_score = 0.0
        if "case" in phenotype_counts and "control" in phenotype_counts:
            case_count = phenotype_counts["case"][alt_idx]
            control_count = phenotype_counts["control"][alt_idx]
            case_af = case_count["ac"] / case_count["an"] if case_count["an"] else 0.0
            control_af = control_count["ac"] / control_count["an"] if control_count["an"] else 0.0
            case_control_af = f"case:{case_af:.4f}|control:{control_af:.4f}"
            case_control_score = abs(case_af - control_af)

        results.append(
            {
                "group_af": "|".join(group_af_parts) if group_af_parts else ".",
                "private_shared": private_shared,
                "fst_score": fst_score,
                "case_control_af": case_control_af,
                "case_control_score": case_control_score,
            }
        )
    return results


def window_index(pos: int, window_size: int) -> int:
    return max(0, (pos - 1) // window_size)


def precompute_population_windows(
    vcf_path: str,
    regions: Sequence[Region],
    window_size: int,
    limit: int = 0,
) -> Tuple[Dict[Tuple[str, int], WindowPopulationStat], int]:
    if window_size <= 0:
        return {}, 0
    accumulators: Dict[Tuple[str, int], WindowAccumulator] = defaultdict(WindowAccumulator)
    n_samples = 0
    records = 0
    with open_text(vcf_path, "r") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                n_samples = max(0, len(header_fields) - 9)
                if n_samples < 10:
                    return {}, n_samples
                continue
            if line.startswith("#"):
                continue
            if n_samples < 10:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            if not in_selected_regions(chrom, pos, regions):
                continue
            alts = fields[4].split(",")
            format_keys = fields[8].split(":")
            sample_values = fields[9:]
            counts = genotype_counts(format_keys, sample_values, len(alts))
            dosages = genotype_alt_dosages(format_keys, sample_values, 1)
            accumulators[(chrom, window_index(pos, window_size))].add_site(counts, dosages)
            records += 1
            if limit and records >= limit:
                break
    return {key: accumulator.finalize() for key, accumulator in accumulators.items()}, n_samples


LOF_CONSEQUENCES = {
    "stop_gained",
    "stop_gained_early",
    "stop_gained_terminal",
    "frameshift",
    "stop_lost",
    "stop_lost_readthrough",
    "start_lost",
    "splice_acceptor_donor",
    "exon_boundary_disruption",
}


def gene_cds_lengths(scorer: "DeNovoPathScorer") -> Dict[str, int]:
    lengths: Dict[str, int] = defaultdict(int)
    for tx in scorer.annotation.transcripts.values():
        lengths[tx.gene_id] = max(lengths[tx.gene_id], tx.cds_length())
    return dict(lengths)


def finalize_gene_constraint(
    observed: Dict[str, Counter[str]],
    lengths: Dict[str, int],
) -> Dict[str, GeneConstraintStat]:
    stats: Dict[str, GeneConstraintStat] = {}
    for gene, counts in observed.items():
        length = max(1, lengths.get(gene, 1000))
        expected_lof = max(0.25, length / 3000.0)
        expected_missense = max(0.50, length / 1000.0)
        lof_obs = counts.get("lof", 0)
        missense_obs = counts.get("missense", 0)
        lof_oe = min(9.9999, lof_obs / expected_lof)
        missense_oe = min(9.9999, missense_obs / expected_missense)
        combined_oe = (2.0 * lof_oe + missense_oe) / 3.0
        stats[gene] = GeneConstraintStat(
            lof_oe=lof_oe,
            missense_oe=missense_oe,
            constraint_score=clamp(1.0 - combined_oe),
            observed_lof=lof_obs,
            observed_missense=missense_obs,
        )
    return stats


def precompute_gene_constraint(
    vcf_path: str,
    regions: Sequence[Region],
    scorer: "DeNovoPathScorer",
    limit: int = 0,
) -> Tuple[Dict[str, GeneConstraintStat], int]:
    n_samples = 0
    records = 0
    observed: Dict[str, Counter[str]] = defaultdict(Counter)
    with open_text(vcf_path, "r") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                n_samples = max(0, len(header_fields) - 9)
                if n_samples < 10:
                    return {}, n_samples
                continue
            if line.startswith("#"):
                continue
            if n_samples < 10:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            if not in_selected_regions(chrom, pos, regions):
                continue
            alts = fields[4].split(",")
            format_keys = fields[8].split(":")
            sample_values = fields[9:]
            counts = genotype_counts(format_keys, sample_values, len(alts))
            for alt_idx, alt in enumerate(alts):
                if counts[alt_idx]["ac"] <= 0:
                    continue
                ann = scorer.annotation_score(chrom, pos, fields[3], alt)
                gene = str(ann.get("gene_id", "."))
                if gene == ".":
                    continue
                consequence = str(ann.get("consequence", "."))
                if consequence in LOF_CONSEQUENCES:
                    observed[gene]["lof"] += 1
                elif consequence == "missense":
                    observed[gene]["missense"] += 1
            records += 1
            if limit and records >= limit:
                break
    return finalize_gene_constraint(observed, gene_cds_lengths(scorer)), n_samples


class DeNovoPathScorer:
    def __init__(
        self,
        reference_fasta: str,
        gff: str,
        cds_fasta: Optional[str] = None,
        protein_fasta: Optional[str] = None,
        protein_domains: Optional[str] = None,
        protein_structures: Optional[str] = None,
        protein_lm_scores: Optional[str] = None,
        mirna_sites: Optional[str] = None,
        ml_model: Optional[str] = None,
        window: int = 50,
        config: Optional[ScoreConfig] = None,
    ) -> None:
        self.reference = Fasta(reference_fasta, sequence_always_upper=True)
        self.annotation = AnnotationIndex(gff)
        self.cds = load_fasta_records(cds_fasta)
        self.proteins = load_fasta_records(protein_fasta)
        self.protein_domains = ProteinDomainIndex(protein_domains)
        self.protein_structures = ProteinStructureIndex(protein_structures)
        self.protein_lm_scores = ProteinLmScoreIndex(protein_lm_scores)
        self.mirna_sites = MirnaSiteIndex(mirna_sites)
        self.ml_model = MlModel(ml_model)
        protein_lm_records = self.proteins or {
            tx_id: "".join(translate_codon(seq[idx : idx + 3]) for idx in range(0, len(seq) - 2, 3))
            for tx_id, seq in self.cds.items()
        }
        self.protein_lm = ProteinKmerLanguageModel(protein_lm_records)
        self.codon_usage = build_codon_usage(self.cds)
        self.dna_lm = DnaKmerLanguageModel(self.reference)
        self.window = window
        self.config = config or ScoreConfig()
        self.splice_pwms = self._build_splice_pwms()

    def validate_reference_allele(self, chrom: str, pos: int, ref: str) -> Tuple[str, str]:
        ref = (ref or "").upper()
        if not ref or any(base not in "ACGTN" for base in ref):
            return "unchecked", "."
        try:
            chrom_len = len(self.reference[chrom])
            if pos < 1 or pos + len(ref) - 1 > chrom_len:
                return "unchecked", "."
            observed = self.reference[chrom][pos - 1 : pos - 1 + len(ref)].seq.upper()
        except Exception:
            return "unchecked", "."
        if observed == ref:
            return "match", observed
        return "mismatch", observed or "."

    def score_record(
        self,
        fields: List[str],
        n_samples: int,
        sample_names: Optional[Sequence[str]] = None,
        sample_info: Optional[SampleInfo] = None,
    ) -> List[AltScore]:
        chrom, pos_text, _id, ref, alt_text = fields[:5]
        pos = int(pos_text)
        alts = alt_text.split(",")
        format_keys = fields[8].split(":") if len(fields) > 8 else []
        sample_values = fields[9:] if len(fields) > 9 else []
        gt_counts = genotype_counts(format_keys, sample_values, len(alts))
        population_scores = group_population_scores(
            format_keys,
            sample_values,
            len(alts),
            sample_names or [],
            sample_info,
        )

        scores = []
        for alt_idx, alt in enumerate(alts):
            cohort = gt_counts[alt_idx]
            population = population_scores[alt_idx]
            seq_score, kmer_score, repeat_score, mutctx_score, dna_lm_score, mutctx_label = self.sequence_scores(
                chrom, pos, ref, alt
            )
            ann = self.annotation_score(chrom, pos, ref, alt)
            cohort_score = self.cohort_score(cohort, max(1, n_samples))
            qc_score = self.qc_score(fields, format_keys, sample_values, alt_idx + 1)
            final_score = self.combine_scores(
                ann["impact_score"],
                ann["protein_score"],
                ann["splice_score"],
                seq_score,
                cohort_score,
                n_samples,
            )
            an = cohort["an"]
            ac = cohort["ac"]
            af = ac / an if an else 0.0
            het_obs, het_exp, het_dev, fis = self.heterozygosity_stats(cohort, n_samples)
            scores.append(
                AltScore(
                    score=final_score,
                    impact_score=ann["impact_score"],
                    protein_score=ann["protein_score"],
                    grantham_score=ann["grantham_score"],
                    blosum_score=ann["blosum_score"],
                    codon_usage_score=ann["codon_usage_score"],
                    protein_context_score=ann["protein_context_score"],
                    protein_structure_score=ann["protein_structure_score"],
                    protein_structure_model_score=ann["protein_structure_model_score"],
                    protein_esm_score=ann["protein_esm_score"],
                    protein_lm_score=ann["protein_lm_score"],
                    protein_domain_score=ann["protein_domain_score"],
                    splice_score=ann["splice_score"],
                    splice_motif_score=ann["splice_motif_score"],
                    splice_pwm_score=ann["splice_pwm_score"],
                    splice_maxent_score=ann["splice_maxent_score"],
                    splice_aux_score=ann["splice_aux_score"],
                    splice_ese_score=ann["splice_ese_score"],
                    utr_score=ann["utr_score"],
                    rnafold_score=ann["rnafold_score"],
                    mirna_score=ann["mirna_score"],
                    promoter_score=ann["promoter_score"],
                    sequence_score=seq_score,
                    kmer_score=kmer_score,
                    repeat_score=repeat_score,
                    mutation_context_score=mutctx_score,
                    dna_lm_score=dna_lm_score,
                    cohort_score=cohort_score,
                    hwe_score=self.hwe_deviation_score(cohort, n_samples),
                    heterozygosity_observed=het_obs,
                    heterozygosity_expected=het_exp,
                    heterozygosity_deviation_score=het_dev,
                    inbreeding_coefficient=fis,
                    fst_score=float(population["fst_score"]),
                    case_control_score=float(population["case_control_score"]),
                    window_pi=0.0,
                    window_theta=0.0,
                    window_tajima_d=0.0,
                    window_ld=0.0,
                    window_haplotype=0.0,
                    gene_lof_oe=0.0,
                    gene_missense_oe=0.0,
                    gene_constraint_score=0.0,
                    qc_score=qc_score,
                    confidence_score=self.confidence_score(final_score, qc_score),
                    consequence=ann["consequence"],
                    level=self.level(final_score),
                    gene_id=ann["gene_id"],
                    tx_id=ann["tx_id"],
                    all_transcripts=ann["all_transcripts"],
                    aa_change=ann["aa_change"],
                    codon_change=ann["codon_change"],
                    hgvs_change=ann["hgvs_change"],
                    normalized_variant=normalized_variant_label(pos, ref, alt),
                    protein_domain_label=ann["protein_domain_label"],
                    protein_structure_label=ann["protein_structure_label"],
                    protein_esm_label=ann["protein_esm_label"],
                    mirna_label=ann["mirna_label"],
                    mutation_context=mutctx_label,
                    maf_bin=self.maf_bin_label(ac, an, n_samples),
                    group_af=str(population["group_af"]),
                    private_shared=str(population["private_shared"]),
                    case_control_af=str(population["case_control_af"]),
                    ac=ac,
                    an=an,
                    af=af,
                    carriers=cohort["carriers"],
                    n_het=cohort["n_het"],
                    n_hom_alt=cohort["n_hom_alt"],
                    n_missing=cohort["n_missing"],
                )
            )
        return scores

    @staticmethod
    def ml_feature_values(score: AltScore) -> Dict[str, float]:
        return {name: float(getattr(score, name, 0.0)) for name in ML_FEATURES}

    def apply_ml_scores(self, score: AltScore) -> None:
        feature_values = self.ml_feature_values(score)
        (
            score.ml_score,
            score.calibrated_score,
            score.uncertainty_score,
            score.ood_score,
        ) = self.ml_model.predict(feature_values)
        score.feature_importance = self.ml_model.format_feature_importance(feature_values)

    def sequence_score(self, chrom: str, pos: int, ref: str, alt: str) -> float:
        return self.sequence_scores(chrom, pos, ref, alt)[0]

    def sequence_scores(self, chrom: str, pos: int, ref: str, alt: str) -> Tuple[float, float, float, float, float, str]:
        try:
            chrom_len = len(self.reference[chrom])
            start = max(1, pos - self.window)
            end = min(chrom_len, pos + len(ref) + self.window - 1)
            seq = self.reference[chrom][start - 1 : end].seq.upper()
        except Exception:
            return 0.20, 0.0, 0.0, 0.0, 0.0, "."
        valid = [base for base in seq if base in "ACGT"]
        if not valid:
            return 0.20, 0.0, 0.0, 0.0, 0.0, "."
        gc = (valid.count("G") + valid.count("C")) / len(valid)
        entropy = shannon_entropy(seq)
        homopolymer = max_homopolymer(seq)
        gc_extreme = min(1.0, abs(gc - 0.45) / 0.45)
        low_complexity = 1.0 - entropy
        homopolymer_score = min(1.0, max(0, homopolymer - 5) / 10)
        length_score = min(1.0, abs(len(alt) - len(ref)) / 12)
        if len(ref) == 1 and len(alt) == 1:
            mut_score = 0.10 if {ref.upper(), alt.upper()} in [{"A", "G"}, {"C", "T"}] else 0.18
        else:
            mut_score = 0.24 + 0.35 * length_score
        context_score = clamp(
            0.28 * low_complexity
            + 0.22 * gc_extreme
            + 0.25 * homopolymer_score
            + 0.25 * mut_score
        )
        kmer_score = trinucleotide_kmer_score(seq, pos - start, ref, alt)
        repeat_score = repeat_low_mappability_proxy_score(seq)
        mutctx_label, mutctx_score = mutation_96_context(seq, pos - start, ref, alt)
        dna_lm_score = self.dna_lm.delta_score(seq, pos - start, ref, alt)
        sequence_score = weighted_mean(
            {
                "context": context_score,
                "kmer": kmer_score,
                "repeat": repeat_score,
                "mutctx": mutctx_score,
                "dna_lm": dna_lm_score,
            },
            self.config.sequence_weights,
        )
        return sequence_score, kmer_score, repeat_score, mutctx_score, dna_lm_score, mutctx_label

    def cohort_score(self, cohort: Dict[str, int], n_samples: int) -> float:
        an = cohort["an"]
        ac = cohort["ac"]
        if not an or not ac:
            return 0.0
        af = ac / an
        if n_samples <= 1:
            return 0.18 if cohort["n_hom_alt"] else 0.12
        if n_samples < 10:
            if cohort["carriers"] == 1:
                return 0.35
            if af >= 0.95:
                return 0.04
            return 0.18
        if af >= 0.95:
            return 0.03
        if af >= 0.50:
            return 0.05
        return clamp(1.0 - (af / 0.10), 0.05, 0.75)

    @staticmethod
    def hwe_deviation_score(cohort: Dict[str, int], n_samples: int) -> float:
        called = cohort.get("n_called_diploid", 0)
        if n_samples < 10 or called < 10:
            return 0.0
        an = called * 2
        ac = cohort["n_het"] + 2 * cohort["n_hom_alt"]
        if an <= 0 or ac <= 0 or ac >= an:
            return 0.0
        p = ac / an
        q = 1.0 - p
        expected = {
            "n_hom_ref": called * q * q,
            "n_het": called * 2.0 * p * q,
            "n_hom_alt": called * p * p,
        }
        if expected["n_het"] < 1.0:
            return 0.0
        observed = {
            "n_hom_ref": cohort.get("n_hom_ref", 0),
            "n_het": cohort["n_het"],
            "n_hom_alt": cohort["n_hom_alt"],
        }
        chi_square = 0.0
        for key, exp_value in expected.items():
            if exp_value > 0:
                chi_square += (observed[key] - exp_value) ** 2 / exp_value
        return clamp(chi_square / (chi_square + 10.0))

    @staticmethod
    def heterozygosity_stats(cohort: Dict[str, int], n_samples: int) -> Tuple[float, float, float, float]:
        called = cohort.get("n_called_diploid", 0)
        if n_samples < 10 or called < 10:
            return 0.0, 0.0, 0.0, 0.0
        an = called * 2
        ac = cohort["n_het"] + 2 * cohort["n_hom_alt"]
        if an <= 0 or ac <= 0 or ac >= an:
            return 0.0, 0.0, 0.0, 0.0
        p = ac / an
        expected = 2.0 * p * (1.0 - p)
        observed = cohort["n_het"] / called
        deviation = abs(observed - expected)
        fis = 0.0 if expected <= 0 else clamp(1.0 - (observed / expected), -1.0, 1.0)
        return observed, expected, deviation, fis

    @staticmethod
    def maf_bin_label(ac: int, an: int, n_samples: int) -> str:
        if n_samples <= 0 or an <= 0:
            return "no_genotypes"
        if n_samples == 1:
            return "single_sample"
        if n_samples < 10:
            return "small_cohort"
        if ac <= 0:
            return "absent"
        if ac >= an or (an - ac) <= 1 or (ac / an) >= 0.95:
            return "fixed_or_near_fixed"
        af = ac / an
        if ac == 1:
            return "singleton"
        if af < 0.001:
            return "ultra_rare"
        if af < 0.01:
            return "rare"
        if af < 0.05:
            return "low_frequency"
        return "common"

    def qc_score(
        self,
        fields: List[str],
        format_keys: List[str],
        sample_values: List[str],
        alt_index: int,
    ) -> float:
        filter_text = fields[6]
        filter_score = 1.0 if filter_text in {".", "PASS"} else 0.55
        try:
            qual = float(fields[5])
            qual_score = clamp(math.log10(max(qual, 0.0) + 1.0) / 3.0)
        except ValueError:
            qual_score = 0.50

        if not sample_values or "GT" not in format_keys:
            return clamp(0.65 * qual_score + 0.35 * filter_score)

        sample_scores = []
        key_to_idx = {key: idx for idx, key in enumerate(format_keys)}
        for sample_text in sample_values:
            values = sample_text.split(":")
            gt = self._fmt_value(values, key_to_idx, "GT")
            if not self._gt_carries_alt(gt, alt_index):
                continue
            dp = self._to_float(self._fmt_value(values, key_to_idx, "DP"))
            gq = self._to_float(self._fmt_value(values, key_to_idx, "GQ"))
            ad_text = self._fmt_value(values, key_to_idx, "AD")
            sb_text = self._fmt_value(values, key_to_idx, "SB")

            dp_score = clamp((dp or 0.0) / 20.0)
            gq_score = clamp((gq or 0.0) / 60.0)
            ad_score, ab_score = self._allele_depth_scores(ad_text, gt, alt_index)
            strand_score = self._strand_score(sb_text)
            sample_scores.append(
                0.25 * dp_score
                + 0.25 * gq_score
                + 0.20 * ad_score
                + 0.20 * ab_score
                + 0.10 * strand_score
            )

        if not sample_scores:
            carrier_score = 0.35
        else:
            carrier_score = sum(sample_scores) / len(sample_scores)
        return clamp((0.15 * qual_score + 0.20 * filter_score + 0.65 * carrier_score) * filter_score)

    @staticmethod
    def _fmt_value(values: List[str], key_to_idx: Dict[str, int], key: str) -> str:
        idx = key_to_idx.get(key)
        if idx is None or idx >= len(values):
            return "."
        return values[idx]

    @staticmethod
    def _to_float(value: str) -> Optional[float]:
        if value in {"", "."}:
            return None
        try:
            return float(value)
        except ValueError:
            return None

    @staticmethod
    def _gt_carries_alt(gt: str, alt_index: int) -> bool:
        if gt in {"", ".", "./.", ".|."}:
            return False
        return str(alt_index) in re.split(r"[/|]", gt)

    @staticmethod
    def _allele_depth_scores(ad_text: str, gt: str, alt_index: int) -> Tuple[float, float]:
        if ad_text in {"", "."}:
            return 0.50, 0.50
        try:
            depths = [float(item) if item != "." else 0.0 for item in ad_text.split(",")]
        except ValueError:
            return 0.50, 0.50
        if alt_index >= len(depths):
            return 0.35, 0.35
        total = sum(depths)
        alt_depth = depths[alt_index]
        if total <= 0:
            return 0.20, 0.20
        ad_score = clamp(alt_depth / 8.0)
        called = [allele for allele in re.split(r"[/|]", gt) if allele.isdigit()]
        if not called:
            return ad_score, 0.50
        expected = sum(1 for allele in called if int(allele) == alt_index) / len(called)
        observed = alt_depth / total
        tolerance = 0.35 if 0.25 <= expected <= 0.75 else 0.25
        ab_score = 1.0 - clamp(abs(observed - expected) / tolerance)
        return ad_score, ab_score

    @staticmethod
    def _strand_score(sb_text: str) -> float:
        if sb_text in {"", "."}:
            return 0.50
        try:
            values = [float(item) if item != "." else 0.0 for item in sb_text.split(",")]
        except ValueError:
            return 0.50
        if len(values) < 4:
            return 0.50
        alt_fwd, alt_rev = values[2], values[3]
        total = alt_fwd + alt_rev
        if total <= 0:
            return 0.50
        balance = min(alt_fwd, alt_rev) / total
        return clamp(balance / 0.20)

    def confidence_score(self, final_score: float, qc_score: float) -> float:
        separation = abs(final_score - 0.5) * 2.0
        weights = self.config.confidence_weights
        weight_sum = sum(weights.values()) or 1.0
        return clamp((weights["qc"] * qc_score + weights["score_separation"] * separation) / weight_sum)

    def level(self, final_score: float) -> str:
        thresholds = self.config.level_thresholds
        if final_score >= thresholds["high"]:
            return "HIGH"
        if final_score >= thresholds["moderate"]:
            return "MODERATE"
        if final_score >= thresholds["low"]:
            return "LOW"
        return "MINIMAL"

    def annotation_score(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, object]:
        candidates = []
        span_start, span_end = changed_ref_span(pos, ref, alt)
        for tx in self.annotation.query_span(chrom, span_start, span_end):
            candidates.append(self._score_transcript(tx, pos, ref, alt))
        if not candidates:
            for tx in self.annotation.query_promoter(chrom, pos):
                candidates.append(self._score_promoter(tx, pos, ref, alt))
        if not candidates:
            return {
                "consequence": "intergenic",
                "impact_score": CONSEQUENCE_BASE_SCORE["intergenic"],
                "protein_score": 0.0,
                "grantham_score": 0.0,
                "blosum_score": 0.0,
                "codon_usage_score": 0.0,
                "protein_context_score": 0.0,
                "protein_structure_score": 0.0,
                "protein_structure_model_score": 0.0,
                "protein_esm_score": 0.0,
                "protein_lm_score": 0.0,
                "protein_domain_score": 0.0,
                "protein_domain_label": ".",
                "protein_structure_label": ".",
                "protein_esm_label": ".",
                "splice_score": 0.0,
                "splice_motif_score": 0.0,
                "splice_pwm_score": 0.0,
                "splice_maxent_score": 0.0,
                "splice_aux_score": 0.0,
                "splice_ese_score": 0.0,
                "utr_score": 0.0,
                "rnafold_score": 0.0,
                "mirna_score": 0.0,
                "mirna_label": ".",
                "promoter_score": 0.0,
                "gene_id": ".",
                "tx_id": ".",
                "all_transcripts": ".",
                "aa_change": ".",
                "codon_change": ".",
                "hgvs_change": ".",
            }
        selected = self.select_annotation_candidate(candidates)
        selected.setdefault("protein_lm_score", 0.0)
        selected.setdefault("protein_structure_model_score", 0.0)
        selected.setdefault("protein_esm_score", 0.0)
        selected.setdefault("protein_domain_score", 0.0)
        selected.setdefault("protein_domain_label", ".")
        selected.setdefault("protein_structure_label", ".")
        selected.setdefault("protein_esm_label", ".")
        selected.setdefault("rnafold_score", 0.0)
        selected.setdefault("mirna_score", 0.0)
        selected.setdefault("mirna_label", ".")
        selected["all_transcripts"] = self.format_all_transcripts(candidates)
        return selected

    @staticmethod
    def format_all_transcripts(candidates: List[Dict[str, object]]) -> str:
        entries = []
        seen = set()
        for item in sorted(
            candidates,
            key=lambda value: (
                str(value.get("gene_id", ".")),
                str(value.get("tx_id", ".")),
                str(value.get("consequence", ".")),
            ),
        ):
            gene = sanitize_info_value(str(item.get("gene_id", ".")))
            tx = sanitize_info_value(str(item.get("tx_id", ".")))
            consequence = sanitize_info_value(str(item.get("consequence", ".")))
            entry = f"{gene}:{tx}:{consequence}"
            if entry not in seen:
                seen.add(entry)
                entries.append(entry)
        return "|".join(entries) if entries else "."

    def select_annotation_candidate(self, candidates: List[Dict[str, object]]) -> Dict[str, object]:
        def score_key(item: Dict[str, object]) -> Tuple[float, float, float, int, int]:
            return (
                float(item["impact_score"]),
                float(item["protein_score"]),
                float(item["splice_score"]),
                int(item.get("cds_length", 0)),
                int(item.get("transcript_length", 0)),
            )

        priority = self.config.transcript_priority
        if priority == "longest_cds":
            return max(candidates, key=lambda item: (int(item.get("cds_length", 0)), score_key(item)))
        if priority == "longest_transcript":
            return max(candidates, key=lambda item: (int(item.get("transcript_length", 0)), score_key(item)))
        if priority == "first":
            return candidates[0]
        return max(candidates, key=score_key)

    def _score_transcript(self, tx: Transcript, pos: int, ref: str, alt: str) -> Dict[str, object]:
        region = tx.region_at(pos)
        span_start, span_end = changed_ref_span(pos, ref, alt)
        boundary_splice_score = self.splice_score(tx, pos, region)
        motif_splice_score = self.splice_motif_score(tx, pos, ref, alt)
        pwm_splice_score = self.splice_pwm_score(tx, pos, ref, alt)
        maxent_splice_score = self.splice_maxent_score(tx, pos, ref, alt)
        aux_splice_score = self.splice_aux_score(tx, pos, ref, alt)
        ese_splice_score = self.splice_ese_score(tx, pos, ref, alt, region)
        splice_score = max(
            boundary_splice_score,
            motif_splice_score,
            pwm_splice_score,
            maxent_splice_score,
            aux_splice_score,
            ese_splice_score,
        )
        utr_score = self.utr_score(tx, pos, ref, alt, region)
        rnafold_score = self.rnafold_delta_g_score(tx, pos, ref, alt, region)
        mirna_score, mirna_label = self.mirna_sites.score(tx, span_start, span_end, region)
        utr_score = max(utr_score, rnafold_score, mirna_score)
        if tx.overlaps_exon_boundary(span_start, span_end):
            boundary_score = 0.95 if tx.boundary_overlaps_cds(span_start, span_end) else 0.78
            splice_score = max(splice_score, boundary_score)
            return {
                "consequence": "exon_boundary_disruption",
                "impact_score": max(CONSEQUENCE_BASE_SCORE["exon_boundary_disruption"], splice_score),
                "protein_score": 0.35 if tx.boundary_overlaps_cds(span_start, span_end) else 0.0,
                "grantham_score": 0.0,
                "blosum_score": 0.0,
                "codon_usage_score": 0.0,
                "protein_context_score": 0.0,
                "protein_structure_score": 0.0,
                "protein_structure_model_score": 0.0,
                "protein_esm_score": 0.0,
                "protein_lm_score": 0.0,
                "protein_domain_score": 0.0,
                "protein_domain_label": ".",
                "protein_structure_label": ".",
                "protein_esm_label": ".",
                "splice_score": splice_score,
                "splice_motif_score": motif_splice_score,
                "splice_pwm_score": pwm_splice_score,
                "splice_maxent_score": maxent_splice_score,
                "splice_aux_score": aux_splice_score,
                "splice_ese_score": ese_splice_score,
                "utr_score": utr_score,
                "rnafold_score": rnafold_score,
                "mirna_score": mirna_score,
                "mirna_label": mirna_label,
                "promoter_score": 0.0,
                "gene_id": tx.gene_id,
                "tx_id": tx.tx_id,
                "aa_change": ".",
                "codon_change": ".",
                "hgvs_change": f"{tx.tx_id}:exon_boundary:{span_start}_{span_end}",
                "cds_length": tx.cds_length(),
                "transcript_length": tx.transcript_length(),
            }
        if region == "CDS":
            coding = self.coding_score(tx, pos, ref, alt)
            coding["splice_score"] = max(splice_score, float(coding["splice_score"]))
            coding["splice_motif_score"] = motif_splice_score
            coding["splice_pwm_score"] = pwm_splice_score
            coding["splice_maxent_score"] = maxent_splice_score
            coding["splice_aux_score"] = aux_splice_score
            coding["splice_ese_score"] = ese_splice_score
            coding["utr_score"] = utr_score
            coding["rnafold_score"] = rnafold_score
            coding["mirna_score"] = mirna_score
            coding["mirna_label"] = mirna_label
            coding["promoter_score"] = 0.0
            return coding
        consequence = region
        if region == "intron" and splice_score >= 0.90:
            consequence = "splice_acceptor_donor"
        elif region == "intron" and splice_score >= 0.35:
            consequence = "splice_region"
        return {
            "consequence": consequence,
            "impact_score": max(CONSEQUENCE_BASE_SCORE.get(consequence, 0.08), splice_score, utr_score),
            "protein_score": 0.0,
            "grantham_score": 0.0,
            "blosum_score": 0.0,
            "codon_usage_score": 0.0,
            "protein_context_score": 0.0,
            "protein_structure_score": 0.0,
            "protein_structure_model_score": 0.0,
            "protein_esm_score": 0.0,
            "protein_lm_score": 0.0,
            "protein_domain_score": 0.0,
            "protein_domain_label": ".",
            "protein_structure_label": ".",
            "protein_esm_label": ".",
            "splice_score": splice_score,
            "splice_motif_score": motif_splice_score,
            "splice_pwm_score": pwm_splice_score,
            "splice_maxent_score": maxent_splice_score,
            "splice_aux_score": aux_splice_score,
            "splice_ese_score": ese_splice_score,
            "utr_score": utr_score,
            "rnafold_score": rnafold_score,
            "mirna_score": mirna_score,
            "mirna_label": mirna_label,
            "promoter_score": 0.0,
            "gene_id": tx.gene_id,
            "tx_id": tx.tx_id,
            "aa_change": ".",
            "codon_change": ".",
            "hgvs_change": ".",
            "cds_length": tx.cds_length(),
            "transcript_length": tx.transcript_length(),
        }

    def protein_sequence(self, tx_id: str) -> str:
        protein = self.proteins.get(tx_id)
        if protein:
            return protein.rstrip("*")
        cds_seq = self.cds.get(tx_id, "")
        aas = []
        for idx in range(0, len(cds_seq) - 2, 3):
            aa = translate_codon(cds_seq[idx : idx + 3])
            if aa == "X":
                continue
            if aa == "*":
                break
            aas.append(aa)
        return "".join(aas)

    def _build_splice_pwms(
        self,
        motif_len: int = 6,
        max_training_introns: int = DEFAULT_MAX_SPLICE_PWM_TRAINING_INTRONS,
    ) -> Dict[str, List[Dict[str, float]]]:
        counts = {
            "donor": [Counter({base: 1.0 for base in "ACGT"}) for _ in range(motif_len)],
            "acceptor": [Counter({base: 1.0 for base in "ACGT"}) for _ in range(motif_len)],
        }
        introns_used = 0
        for tx in self.annotation.transcripts.values():
            if introns_used >= max_training_introns:
                break
            exons = tx.exon_segments()
            if len(exons) < 2:
                continue
            for left_exon, right_exon in zip(exons, exons[1:]):
                if introns_used >= max_training_introns:
                    break
                intron_start = left_exon[1] + 1
                intron_end = right_exon[0] - 1
                if intron_end - intron_start + 1 < motif_len:
                    continue
                if tx.strand == "-":
                    donor = self._oriented_reference(tx.chrom, intron_end - motif_len + 1, intron_end, "-")
                    acceptor = self._oriented_reference(tx.chrom, intron_start, intron_start + motif_len - 1, "-")
                else:
                    donor = self._oriented_reference(tx.chrom, intron_start, intron_start + motif_len - 1, "+")
                    acceptor = self._oriented_reference(tx.chrom, intron_end - motif_len + 1, intron_end, "+")
                for label, seq in (("donor", donor), ("acceptor", acceptor)):
                    if len(seq) != motif_len or any(base not in "ACGT" for base in seq):
                        continue
                    for idx, base in enumerate(seq):
                        counts[label][idx][base] += 1.0
                introns_used += 1
        pwms = {}
        for label, label_counts in counts.items():
            pwm = []
            for pos_counts in label_counts:
                total = sum(pos_counts.values())
                pwm.append({base: pos_counts[base] / total for base in "ACGT"})
            pwms[label] = pwm
        return pwms

    def _oriented_reference(self, chrom: str, start: int, end: int, strand: str) -> str:
        if start < 1 or end < start:
            return ""
        try:
            chrom_len = len(self.reference[chrom])
            if end > chrom_len:
                return ""
            seq = self.reference[chrom][start - 1 : end].seq.upper()
        except Exception:
            return ""
        return revcomp(seq) if strand == "-" else seq

    def splice_pwm_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> float:
        motif_len = len(self.splice_pwms.get("donor", []))
        if motif_len == 0:
            return 0.0
        exons = tx.exon_segments()
        if len(exons) < 2:
            return 0.0
        best = 0.0
        for left_exon, right_exon in zip(exons, exons[1:]):
            intron_start = left_exon[1] + 1
            intron_end = right_exon[0] - 1
            if intron_end - intron_start + 1 < motif_len:
                continue
            if tx.strand == "-":
                donor = self._splice_pwm_delta(tx.chrom, intron_end - motif_len + 1, intron_end, "-", pos, ref, alt, "donor")
                acceptor = self._splice_pwm_delta(
                    tx.chrom, intron_start, intron_start + motif_len - 1, "-", pos, ref, alt, "acceptor"
                )
            else:
                donor = self._splice_pwm_delta(tx.chrom, intron_start, intron_start + motif_len - 1, "+", pos, ref, alt, "donor")
                acceptor = self._splice_pwm_delta(
                    tx.chrom, intron_end - motif_len + 1, intron_end, "+", pos, ref, alt, "acceptor"
                )
            best = max(best, donor, acceptor)
        return best

    def _splice_pwm_delta(
        self,
        chrom: str,
        window_start: int,
        window_end: int,
        strand: str,
        pos: int,
        ref: str,
        alt: str,
        label: str,
    ) -> float:
        variant_start = pos
        variant_end = pos + len(ref) - 1
        if variant_end < window_start or variant_start > window_end:
            return 0.0
        try:
            window_seq = self.reference[chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return 0.0
        oriented_before = revcomp(window_seq) if strand == "-" else window_seq
        before = self._score_splice_pwm(oriented_before, label)
        if len(ref) != len(alt):
            return clamp(before)
        mutated = list(window_seq)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(mutated):
                mutated[local_idx] = base
        oriented_after = revcomp("".join(mutated)) if strand == "-" else "".join(mutated)
        after = self._score_splice_pwm(oriented_after, label)
        return clamp((before - after) / 0.25)

    def _score_splice_pwm(self, seq: str, label: str) -> float:
        pwm = self.splice_pwms.get(label, [])
        if len(seq) != len(pwm) or not pwm:
            return 0.0
        probs = []
        for idx, base in enumerate(seq.upper()):
            if base not in "ACGT":
                return 0.0
            probs.append(pwm[idx].get(base, 0.0))
        return sum(probs) / len(probs)

    def splice_maxent_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> float:
        exons = tx.exon_segments()
        if len(exons) < 2:
            return 0.0
        best = 0.0
        for left_exon, right_exon in zip(exons, exons[1:]):
            intron_start = left_exon[1] + 1
            intron_end = right_exon[0] - 1
            if intron_end - intron_start + 1 < 6:
                continue
            if tx.strand == "-":
                donor_start, donor_end = intron_end - 5, intron_end + 3
                acceptor_start, acceptor_end = intron_start - 3, intron_start + 19
            else:
                donor_start, donor_end = intron_start - 3, intron_start + 5
                acceptor_start, acceptor_end = intron_end - 19, intron_end + 3
            best = max(
                best,
                self._splice_maxent_delta(tx.chrom, donor_start, donor_end, tx.strand, pos, ref, alt, "donor"),
                self._splice_maxent_delta(tx.chrom, acceptor_start, acceptor_end, tx.strand, pos, ref, alt, "acceptor"),
            )
        return best

    def _splice_maxent_delta(
        self,
        chrom: str,
        window_start: int,
        window_end: int,
        strand: str,
        pos: int,
        ref: str,
        alt: str,
        label: str,
    ) -> float:
        variant_end = pos + len(ref) - 1
        if window_start < 1 or window_end < window_start or variant_end < window_start or pos > window_end:
            return 0.0
        try:
            chrom_len = len(self.reference[chrom])
            if window_end > chrom_len:
                return 0.0
            raw_before = self.reference[chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return 0.0
        oriented_before = revcomp(raw_before) if strand == "-" else raw_before
        before = self._score_maxent_like_window(oriented_before, label)
        if len(ref) != len(alt):
            return clamp(0.20 + 0.65 * before)
        raw_after = list(raw_before)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(raw_after):
                raw_after[local_idx] = base
        oriented_after = revcomp("".join(raw_after)) if strand == "-" else "".join(raw_after)
        after = self._score_maxent_like_window(oriented_after, label)
        disruption = max(0.0, before - after)
        creation = max(0.0, after - before)
        if disruption > 0:
            return clamp(0.12 + 0.88 * disruption)
        if creation > 0:
            return clamp(0.04 + 0.35 * creation)
        return clamp(0.03 * before)

    @staticmethod
    def _score_maxent_like_window(seq: str, label: str) -> float:
        seq = seq.upper()
        if any(base not in "ACGT" for base in seq):
            return 0.0
        if label == "donor":
            if len(seq) != 9:
                return 0.0
            weights = [0.06, 0.06, 0.08, 0.24, 0.24, 0.10, 0.08, 0.08, 0.06]
            preferences = [
                {"A", "G"},
                {"A", "C", "G", "T"},
                {"G"},
                {"G"},
                {"T"},
                {"A", "G"},
                {"A", "G"},
                {"G"},
                {"T"},
            ]
            score = sum(weight for base, weight, pref in zip(seq, weights, preferences) if base in pref)
            dinuc_bonus = 0.18 if seq[3:5] == "GT" else 0.0
            return clamp(score + dinuc_bonus)
        if label == "acceptor":
            if len(seq) != 23:
                return 0.0
            polypy = seq[:18]
            py_frac = sum(1 for base in polypy if base in {"C", "T"}) / len(polypy)
            py_run = 0
            current = 0
            for base in polypy:
                if base in {"C", "T"}:
                    current += 1
                    py_run = max(py_run, current)
                else:
                    current = 0
            ag = 1.0 if seq[18:20] == "AG" else DeNovoPathScorer._motif_match_fraction(seq[18:20], "AG")
            exon_bonus = 0.0
            if seq[20] in {"G", "A"}:
                exon_bonus += 0.08
            if seq[22] == "G":
                exon_bonus += 0.06
            return clamp(0.45 * ag + 0.30 * py_frac + 0.15 * (py_run / len(polypy)) + exon_bonus)
        return 0.0

    def splice_score(self, tx: Transcript, pos: int, region: str) -> float:
        dist = tx.min_exon_boundary_distance(pos)
        if dist is None:
            return 0.0
        if region == "intron":
            if dist <= 2:
                return 0.95
            if dist <= 8:
                return 0.65
            if dist <= 20:
                return 0.35
        elif region in {"CDS", "utr5", "utr3"}:
            if dist <= 2:
                return 0.45
            if dist <= 8:
                return 0.25
        return 0.0

    def splice_motif_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> float:
        exons = tx.exon_segments()
        if len(exons) < 2:
            return 0.0
        best = 0.0
        for left_exon, right_exon in zip(exons, exons[1:]):
            intron_start = left_exon[1] + 1
            intron_end = right_exon[0] - 1
            if intron_end - intron_start + 1 < 2:
                continue
            if tx.strand == "-":
                donor = self._splice_motif_delta(tx.chrom, intron_end - 1, intron_end, "-", pos, ref, alt, "GT")
                acceptor = self._splice_motif_delta(tx.chrom, intron_start, intron_start + 1, "-", pos, ref, alt, "AG")
            else:
                donor = self._splice_motif_delta(tx.chrom, intron_start, intron_start + 1, "+", pos, ref, alt, "GT")
                acceptor = self._splice_motif_delta(tx.chrom, intron_end - 1, intron_end, "+", pos, ref, alt, "AG")
            best = max(best, donor, acceptor)
        return best

    def _splice_motif_delta(
        self,
        chrom: str,
        motif_start: int,
        motif_end: int,
        strand: str,
        pos: int,
        ref: str,
        alt: str,
        expected: str,
    ) -> float:
        variant_start = pos
        variant_end = pos + len(ref) - 1
        if variant_end < motif_start or variant_start > motif_end:
            return 0.0
        try:
            motif_seq = self.reference[chrom][motif_start - 1 : motif_end].seq.upper()
        except Exception:
            return 0.0
        if strand == "-":
            motif_before = revcomp(motif_seq)
        else:
            motif_before = motif_seq
        before = self._motif_match_fraction(motif_before, expected)
        if len(ref) != len(alt):
            return 1.0 if before >= 1.0 else 0.65 * before

        mutated = list(motif_seq)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            if motif_start <= genome_pos <= motif_end:
                mutated[genome_pos - motif_start] = base
        motif_after_seq = "".join(mutated)
        motif_after = revcomp(motif_after_seq) if strand == "-" else motif_after_seq
        after = self._motif_match_fraction(motif_after, expected)
        disruption = max(0.0, before - after)
        if before >= 1.0 and after < 1.0:
            return clamp(0.85 + 0.15 * disruption)
        return clamp(disruption)

    @staticmethod
    def _motif_match_fraction(seq: str, expected: str) -> float:
        if len(seq) != len(expected):
            return 0.0
        matches = sum(1 for observed, wanted in zip(seq.upper(), expected.upper()) if observed == wanted)
        return matches / len(expected)

    def splice_aux_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> float:
        exons = tx.exon_segments()
        if len(exons) < 2:
            return 0.0
        best = 0.0
        for left_exon, right_exon in zip(exons, exons[1:]):
            intron_start = left_exon[1] + 1
            intron_end = right_exon[0] - 1
            if intron_end - intron_start + 1 < 25:
                continue
            if tx.strand == "-":
                branch_start = max(intron_start, intron_start + 18)
                branch_end = min(intron_end, intron_start + 40)
                polypy_start = max(intron_start, intron_start + 5)
                polypy_end = min(intron_end, intron_start + 20)
            else:
                branch_start = max(intron_start, intron_end - 40)
                branch_end = min(intron_end, intron_end - 18)
                polypy_start = max(intron_start, intron_end - 20)
                polypy_end = min(intron_end, intron_end - 5)
            best = max(
                best,
                self._splice_branch_delta(tx.chrom, branch_start, branch_end, tx.strand, pos, ref, alt),
                self._splice_polypy_delta(tx.chrom, polypy_start, polypy_end, tx.strand, pos, ref, alt),
            )
        return best

    def _splice_branch_delta(
        self,
        chrom: str,
        window_start: int,
        window_end: int,
        strand: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> float:
        if window_end < window_start or pos + len(ref) - 1 < window_start or pos > window_end:
            return 0.0
        try:
            raw_before = self.reference[chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return 0.0
        before = self._branchpoint_strength(raw_before, strand)
        if len(ref) != len(alt):
            return clamp(0.15 + 0.55 * before)
        raw_after = list(raw_before)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(raw_after):
                raw_after[local_idx] = base
        after = self._branchpoint_strength("".join(raw_after), strand)
        disruption = max(0.0, before - after)
        creation = max(0.0, after - before)
        if disruption > 0:
            return clamp(0.18 + 0.62 * disruption)
        if creation > 0:
            return clamp(0.08 + 0.30 * creation)
        return clamp(0.05 * before)

    def _splice_polypy_delta(
        self,
        chrom: str,
        window_start: int,
        window_end: int,
        strand: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> float:
        if window_end < window_start or pos + len(ref) - 1 < window_start or pos > window_end:
            return 0.0
        try:
            raw_before = self.reference[chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return 0.0
        before = self._polypyrimidine_strength(raw_before, strand)
        if len(ref) != len(alt):
            return clamp(0.12 + 0.50 * before)
        raw_after = list(raw_before)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(raw_after):
                raw_after[local_idx] = base
        after = self._polypyrimidine_strength("".join(raw_after), strand)
        disruption = max(0.0, before - after)
        creation = max(0.0, after - before)
        if disruption > 0:
            return clamp(0.12 + 0.55 * disruption)
        if creation > 0:
            return clamp(0.05 + 0.20 * creation)
        return clamp(0.04 * before)

    @staticmethod
    def _branchpoint_strength(raw_seq: str, strand: str) -> float:
        seq = revcomp(raw_seq) if strand == "-" else raw_seq.upper()
        motifs = {
            "YTNAY": 1.00,
            "CTRAY": 0.92,
            "CTAAC": 0.88,
            "CTAAT": 0.78,
        }
        best = 0.0
        for motif, weight in motifs.items():
            best = max(best, weight * DeNovoPathScorer._best_iupac_motif_fraction(seq, motif))
        return clamp(best)

    @staticmethod
    def _polypyrimidine_strength(raw_seq: str, strand: str) -> float:
        seq = revcomp(raw_seq) if strand == "-" else raw_seq.upper()
        valid = [base for base in seq if base in "ACGT"]
        if not valid:
            return 0.0
        pyrimidines = sum(1 for base in valid if base in {"C", "T"})
        max_run = 0
        current = 0
        for base in valid:
            if base in {"C", "T"}:
                current += 1
                max_run = max(max_run, current)
            else:
                current = 0
        return clamp(0.65 * (pyrimidines / len(valid)) + 0.35 * (max_run / len(valid)))

    def splice_ese_score(self, tx: Transcript, pos: int, ref: str, alt: str, region: str) -> float:
        if region not in {"CDS", "utr5", "utr3"}:
            return 0.0
        variant_end = pos + len(ref) - 1
        exon_hits = [segment for segment in tx.exon_segments() if not (variant_end < segment[0] or pos > segment[1])]
        if not exon_hits:
            return 0.0
        exon_start = min(segment[0] for segment in exon_hits)
        exon_end = max(segment[1] for segment in exon_hits)
        window_start = max(exon_start, pos - 12)
        window_end = min(exon_end, variant_end + 12)
        try:
            raw_before = self.reference[tx.chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return 0.0
        before_ese, before_ess = self._ese_ess_strength(raw_before, tx.strand)
        if len(ref) != len(alt):
            return clamp(0.12 + 0.35 * max(before_ese, before_ess))
        raw_after = list(raw_before)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(raw_after):
                raw_after[local_idx] = base
        after_ese, after_ess = self._ese_ess_strength("".join(raw_after), tx.strand)
        enhancer_loss = max(0.0, before_ese - after_ese)
        silencer_gain = max(0.0, after_ess - before_ess)
        enhancer_gain = max(0.0, after_ese - before_ese)
        silencer_loss = max(0.0, before_ess - after_ess)
        harmful = max(enhancer_loss, silencer_gain)
        protective = max(enhancer_gain, silencer_loss)
        if harmful > 0:
            return clamp(0.10 + 0.70 * harmful)
        if protective > 0:
            return clamp(0.04 + 0.12 * protective)
        return clamp(0.03 * max(before_ese, before_ess))

    @staticmethod
    def _ese_ess_strength(raw_seq: str, strand: str) -> Tuple[float, float]:
        seq = revcomp(raw_seq) if strand == "-" else raw_seq.upper()
        ese_motifs = {
            "GAAGAA": 1.00,
            "AAGAAG": 0.95,
            "GAAAGA": 0.90,
            "CAGGAA": 0.82,
            "ACGCGG": 0.72,
            "GAR": 0.45,
        }
        ess_motifs = {
            "TTTGGG": 1.00,
            "TTCCTT": 0.90,
            "CTAGTT": 0.82,
            "TAGGGA": 0.78,
            "GGGG": 0.58,
            "TTTT": 0.55,
        }
        ese = 0.0
        ess = 0.0
        for motif, weight in ese_motifs.items():
            ese = max(ese, weight * DeNovoPathScorer._has_iupac_motif(seq, motif))
        for motif, weight in ess_motifs.items():
            ess = max(ess, weight * DeNovoPathScorer._has_iupac_motif(seq, motif))
        return clamp(ese), clamp(ess)

    @staticmethod
    def _has_iupac_motif(seq: str, motif: str) -> float:
        allowed = {
            "A": {"A"},
            "C": {"C"},
            "G": {"G"},
            "T": {"T"},
            "W": {"A", "T"},
            "S": {"C", "G"},
            "R": {"A", "G"},
            "Y": {"C", "T"},
            "N": {"A", "C", "G", "T"},
        }
        seq = seq.upper()
        motif = motif.upper()
        if len(seq) < len(motif):
            return 0.0
        for start in range(0, len(seq) - len(motif) + 1):
            if all(seq[start + idx] in allowed.get(code, {code}) for idx, code in enumerate(motif)):
                return 1.0
        return 0.0

    def utr_score(self, tx: Transcript, pos: int, ref: str, alt: str, region: str) -> float:
        if region == "utr5":
            return self.utr5_start_context_score(tx, pos, ref, alt)
        if region == "utr3":
            return self.utr3_polya_score(tx, pos, ref, alt)
        return 0.0

    def rnafold_delta_g_score(self, tx: Transcript, pos: int, ref: str, alt: str, region: str) -> float:
        if region not in {"utr5", "utr3"}:
            return 0.0
        before, after = self._oriented_variant_window(tx, pos, ref, alt, flank=45)
        if len(before) < 12 or len(after) < 12:
            return 0.0
        before_energy = rna_pairing_energy_proxy(before)
        after_energy = rna_pairing_energy_proxy(after)
        delta = abs(after_energy - before_energy)
        disruption = max(0.0, after_energy - before_energy)
        creation = max(0.0, before_energy - after_energy)
        baseline = max(1.0, abs(before_energy))
        relative = delta / baseline
        indel_boost = min(0.20, abs(len(after) - len(before)) / 25.0)
        return clamp(0.55 * min(1.0, delta / 4.0) + 0.25 * min(1.0, relative) + 0.15 * min(1.0, disruption / 3.0) + 0.05 * min(1.0, creation / 3.0) + indel_boost)

    def _oriented_variant_window(
        self,
        tx: Transcript,
        pos: int,
        ref: str,
        alt: str,
        flank: int = 45,
    ) -> Tuple[str, str]:
        try:
            chrom_len = len(self.reference[tx.chrom])
            window_start = max(1, pos - flank)
            window_end = min(chrom_len, pos + len(ref) + flank - 1)
            before = self.reference[tx.chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return "", ""
        rel_start = max(0, pos - window_start)
        rel_end = min(len(before), rel_start + len(ref))
        after = before[:rel_start] + alt.upper() + before[rel_end:]
        if tx.strand == "-":
            return revcomp(before), revcomp(after)
        return before, after

    def utr5_start_context_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> float:
        cds_order = tx.transcript_order_cds()
        if not cds_order:
            return 0.0
        first_cds = cds_order[0]
        if tx.strand == "-":
            start_codon_anchor = first_cds[1]
            window_start = start_codon_anchor - 2
            window_end = start_codon_anchor + 6
        else:
            start_codon_anchor = first_cds[0]
            window_start = start_codon_anchor - 6
            window_end = start_codon_anchor + 2
        if pos + len(ref) - 1 < window_start or pos > window_end:
            distance = abs(pos - start_codon_anchor)
            if distance <= 10:
                return 0.25
            if distance <= 50:
                return 0.12
            return 0.04
        before = self._oriented_reference(tx.chrom, window_start, window_end, tx.strand)
        if len(before) != 9:
            return 0.20
        before_score = self._start_context_strength(before)
        if len(ref) != len(alt):
            return clamp(0.45 + 0.35 * before_score)
        try:
            raw = list(self.reference[tx.chrom][window_start - 1 : window_end].seq.upper())
        except Exception:
            return 0.20
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(raw):
                raw[local_idx] = base
        after_seq = "".join(raw)
        after = revcomp(after_seq) if tx.strand == "-" else after_seq
        after_score = self._start_context_strength(after)
        disruption = max(0.0, before_score - after_score)
        return clamp(max(0.18, 0.30 * before_score + 0.70 * disruption))

    @staticmethod
    def _start_context_strength(seq: str) -> float:
        seq = seq.upper()
        if len(seq) != 9:
            return 0.0
        start_match = sum(1 for observed, expected in zip(seq[6:9], "ATG") if observed == expected) / 3.0
        minus3 = 1.0 if seq[3] in {"A", "G"} else 0.0
        minus1 = 1.0 if seq[5] in {"A", "C"} else 0.0
        return clamp(0.55 * start_match + 0.30 * minus3 + 0.15 * minus1)

    def utr3_polya_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> float:
        try:
            chrom_len = len(self.reference[tx.chrom])
            window_start = max(1, pos - 25)
            window_end = min(chrom_len, pos + len(ref) + 25)
            raw_before = self.reference[tx.chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return 0.08
        before = revcomp(raw_before) if tx.strand == "-" else raw_before
        before_strength = self._polya_strength(before)
        if len(ref) != len(alt):
            return clamp(0.25 + 0.45 * before_strength)
        raw_after = list(raw_before)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(raw_after):
                raw_after[local_idx] = base
        after_seq = "".join(raw_after)
        after = revcomp(after_seq) if tx.strand == "-" else after_seq
        after_strength = self._polya_strength(after)
        disruption = max(0.0, before_strength - after_strength)
        creation = max(0.0, after_strength - before_strength)
        return clamp(max(0.08, 0.55 * disruption + 0.25 * creation))

    @staticmethod
    def _polya_strength(seq: str) -> float:
        motifs = {"AATAAA": 1.0, "ATTAAA": 0.85, "TATAAA": 0.45, "AGTAAA": 0.40}
        best = 0.0
        seq = seq.upper()
        for motif, weight in motifs.items():
            if motif in seq:
                best = max(best, weight)
        return best

    def _score_promoter(self, tx: Transcript, pos: int, ref: str, alt: str) -> Dict[str, object]:
        promoter_score = self.promoter_score(tx, pos, ref, alt)
        return {
            "consequence": "promoter",
            "impact_score": max(CONSEQUENCE_BASE_SCORE["promoter"], promoter_score),
            "protein_score": 0.0,
            "grantham_score": 0.0,
            "blosum_score": 0.0,
            "codon_usage_score": 0.0,
            "protein_context_score": 0.0,
            "protein_structure_score": 0.0,
            "protein_structure_model_score": 0.0,
            "protein_esm_score": 0.0,
            "protein_lm_score": 0.0,
            "protein_domain_score": 0.0,
            "protein_domain_label": ".",
            "protein_structure_label": ".",
            "protein_esm_label": ".",
            "splice_score": 0.0,
            "splice_motif_score": 0.0,
            "splice_pwm_score": 0.0,
            "splice_maxent_score": 0.0,
            "splice_aux_score": 0.0,
            "splice_ese_score": 0.0,
            "utr_score": 0.0,
            "rnafold_score": 0.0,
            "mirna_score": 0.0,
            "mirna_label": ".",
            "promoter_score": promoter_score,
            "gene_id": tx.gene_id,
            "tx_id": tx.tx_id,
            "aa_change": ".",
            "codon_change": ".",
            "hgvs_change": ".",
            "cds_length": tx.cds_length(),
            "transcript_length": tx.transcript_length(),
        }

    def promoter_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> float:
        promoter_start, promoter_end = self.annotation.promoter_interval(tx)
        variant_end = pos + len(ref) - 1
        if variant_end < promoter_start or pos > promoter_end:
            return 0.0
        motif_flank = 8
        window_start = max(promoter_start, pos - motif_flank)
        window_end = min(promoter_end, variant_end + motif_flank)
        try:
            raw_before = self.reference[tx.chrom][window_start - 1 : window_end].seq.upper()
        except Exception:
            return 0.08
        if len(ref) != len(alt):
            distance = self._promoter_distance_from_tss(tx, pos)
            proximity = 1.0 - clamp((distance - 1) / max(1, self.annotation.promoter_upstream))
            return clamp(0.18 + 0.35 * proximity + 0.30 * self._promoter_motif_strength(raw_before, tx.strand))
        raw_after = list(raw_before)
        for idx, base in enumerate(alt.upper()):
            genome_pos = pos + idx
            local_idx = genome_pos - window_start
            if 0 <= local_idx < len(raw_after):
                raw_after[local_idx] = base
        before_strength = self._promoter_motif_strength(raw_before, tx.strand)
        after_strength = self._promoter_motif_strength("".join(raw_after), tx.strand)
        disruption = max(0.0, before_strength - after_strength)
        creation = max(0.0, after_strength - before_strength)
        distance = self._promoter_distance_from_tss(tx, pos)
        proximity = 1.0 - clamp((distance - 1) / max(1, self.annotation.promoter_upstream))
        if disruption > 0:
            return clamp(0.18 + 0.62 * disruption + 0.20 * proximity)
        if creation > 0:
            return clamp(0.10 + 0.35 * creation + 0.10 * proximity)
        if before_strength > 0:
            return clamp(0.08 + 0.10 * before_strength + 0.08 * proximity)
        return clamp(0.04 + 0.06 * proximity)

    @staticmethod
    def _promoter_motif_strength(raw_seq: str, strand: str) -> float:
        seq = revcomp(raw_seq) if strand == "-" else raw_seq.upper()
        motif_groups = [
            (("TATAAA", "TATATA", "TATAWA"), 1.0),
            (("CCAAT", "CAAT"), 0.72),
            (("GGGCGG", "CCGCCC", "GCGGGG"), 0.58),
        ]
        best = 0.0
        for motifs, weight in motif_groups:
            for motif in motifs:
                best = max(best, weight * DeNovoPathScorer._best_iupac_motif_fraction(seq, motif))
        return clamp(best)

    @staticmethod
    def _best_iupac_motif_fraction(seq: str, motif: str) -> float:
        allowed = {
            "A": {"A"},
            "C": {"C"},
            "G": {"G"},
            "T": {"T"},
            "W": {"A", "T"},
            "S": {"C", "G"},
            "R": {"A", "G"},
            "Y": {"C", "T"},
            "N": {"A", "C", "G", "T"},
        }
        seq = seq.upper()
        motif = motif.upper()
        if len(seq) < len(motif):
            return 0.0
        best = 0.0
        for start in range(0, len(seq) - len(motif) + 1):
            matches = 0
            for idx, code in enumerate(motif):
                if seq[start + idx] in allowed.get(code, {code}):
                    matches += 1
            best = max(best, matches / len(motif))
        return best

    @staticmethod
    def _promoter_distance_from_tss(tx: Transcript, pos: int) -> int:
        tss = tx.end if tx.strand == "-" else tx.start
        return abs(tss - pos)

    def coding_score(self, tx: Transcript, pos: int, ref: str, alt: str) -> Dict[str, object]:
        norm_pos, norm_ref, norm_alt = normalize_allele(pos, ref, alt)
        offset = self._coding_event_offset(tx, norm_pos, norm_ref)
        cds_seq = self.cds.get(tx.tx_id)
        if offset is None or not cds_seq:
            return self._coding_result(tx, "cds_complex", 0.50, 0.35)
        approx_aa_pos = offset // 3 + 1
        approx_context = protein_context_score(self.protein_sequence(tx.tx_id), approx_aa_pos)
        approx_domain_score, approx_domain_label = self.protein_domains.score(tx.tx_id, approx_aa_pos)
        approx_struct_model_score, approx_struct_model_label = self.protein_structures.score(
            tx.tx_id, approx_aa_pos, "X", "X"
        )
        approx_esm_score, approx_esm_label = self.protein_lm_scores.score(tx.tx_id, approx_aa_pos, ".", ".")
        if len(norm_ref) != len(norm_alt):
            delta = len(norm_alt) - len(norm_ref)
            hgvs_change = self._hgvs_like_indel(tx, offset, norm_ref, norm_alt)
            if delta % 3:
                return self._coding_result(
                    tx,
                    "frameshift",
                    0.98,
                    1.00,
                    protein_context_value=approx_context,
                    protein_structure_model_value=approx_struct_model_score,
                    protein_structure_label=approx_struct_model_label,
                    protein_esm_value=approx_esm_score,
                    protein_esm_label=approx_esm_label,
                    protein_domain_value=approx_domain_score,
                    protein_domain_label=approx_domain_label,
                    hgvs_change=hgvs_change,
                )
            if delta > 0:
                return self._coding_result(
                    tx,
                    "inframe_insertion",
                    0.62,
                    0.62,
                    protein_context_value=approx_context,
                    protein_structure_model_value=approx_struct_model_score,
                    protein_structure_label=approx_struct_model_label,
                    protein_esm_value=approx_esm_score,
                    protein_esm_label=approx_esm_label,
                    protein_domain_value=approx_domain_score,
                    protein_domain_label=approx_domain_label,
                    hgvs_change=hgvs_change,
                )
            return self._coding_result(
                tx,
                "inframe_deletion",
                0.68,
                0.68,
                protein_context_value=approx_context,
                protein_structure_model_value=approx_struct_model_score,
                protein_structure_label=approx_struct_model_label,
                protein_esm_value=approx_esm_score,
                protein_esm_label=approx_esm_label,
                protein_domain_value=approx_domain_score,
                protein_domain_label=approx_domain_label,
                hgvs_change=hgvs_change,
            )

        mutated = list(cds_seq)
        changed_offsets = []
        for i, base in enumerate(norm_alt.upper()):
            genomic_pos = norm_pos + i
            idx = tx.cds_offset(genomic_pos)
            if idx is None or idx >= len(mutated):
                return self._coding_result(tx, "cds_complex", 0.50, 0.35)
            mutated[idx] = base if tx.strand == "+" else base.translate(DNA_COMP).upper()
            changed_offsets.append(idx)
        codon_starts = {(idx // 3) * 3 for idx in changed_offsets}
        if len(codon_starts) != 1:
            return self.multi_codon_substitution_score(tx, cds_seq, mutated, codon_starts, changed_offsets)
        codon_start = codon_starts.pop()
        ref_codon = cds_seq[codon_start : codon_start + 3].upper()
        alt_codon = "".join(mutated[codon_start : codon_start + 3]).upper()
        ref_aa = translate_codon(ref_codon)
        alt_aa = translate_codon(alt_codon)
        codon_usage = codon_usage_bias_score(ref_codon, alt_codon, self.codon_usage)
        aa_pos = codon_start // 3 + 1
        protein_seq = self.protein_sequence(tx.tx_id)
        prot_context = protein_context_score(protein_seq, aa_pos)
        prot_structure = protein_structure_score(protein_seq, aa_pos, ref_aa, alt_aa)
        prot_structure_model, prot_structure_label = self.protein_structures.score(tx.tx_id, aa_pos, ref_aa, alt_aa)
        prot_esm, prot_esm_label = self.protein_lm_scores.score(tx.tx_id, aa_pos, ref_aa, alt_aa)
        prot_lm = self.protein_lm.delta_score(protein_seq, aa_pos, ref_aa, alt_aa)
        domain_score, domain_label = self.protein_domains.score(tx.tx_id, aa_pos)
        codon_change = f"{ref_codon}>{alt_codon}"
        aa_change = f"p.{ref_aa}{aa_pos}{alt_aa}"
        hgvs_change = self._hgvs_like_substitution(tx, cds_seq, mutated, changed_offsets, aa_change)
        if aa_pos == 1 and ref_aa == "M" and alt_aa != "M":
            g_score = grantham_score(ref_aa, alt_aa)
            b_score = blosum62_score(ref_aa, alt_aa)
            prot_score = max(
                0.88,
                self.protein_substitution_score(
                    ref_aa,
                    alt_aa,
                    g_score,
                    b_score,
                    codon_usage,
                    prot_context,
                    prot_structure,
                    prot_structure_model,
                    prot_esm,
                    prot_lm,
                    domain_score,
                ),
            )
            return self._coding_result(
                tx,
                "start_lost",
                0.92,
                prot_score,
                aa_change,
                codon_change,
                g_score,
                b_score,
                codon_usage,
                prot_context,
                prot_structure,
                prot_structure_model,
                prot_esm,
                prot_lm,
                domain_score,
                domain_label,
                prot_structure_label,
                prot_esm_label,
                hgvs_change,
            )
        if ref_aa == "*" and alt_aa == "*":
            return self._coding_result(
                tx,
                "stop_retained",
                0.04,
                0.02,
                aa_change,
                codon_change,
                protein_context_value=prot_context,
                protein_structure_model_value=prot_structure_model,
                protein_structure_label=prot_structure_label,
                protein_esm_value=prot_esm,
                protein_esm_label=prot_esm_label,
                protein_domain_value=domain_score,
                protein_domain_label=domain_label,
                hgvs_change=hgvs_change,
            )
        if ref_aa == alt_aa:
            prot_score = max(0.02, 0.18 * codon_usage)
            return self._coding_result(
                tx,
                "synonymous",
                0.05,
                prot_score,
                aa_change,
                codon_change,
                codon_usage_value=codon_usage,
                protein_context_value=prot_context,
                protein_structure_value=prot_structure,
                protein_structure_model_value=prot_structure_model,
                protein_esm_value=prot_esm,
                protein_lm_value=prot_lm,
                protein_domain_value=domain_score,
                protein_domain_label=domain_label,
                protein_structure_label=prot_structure_label,
                protein_esm_label=prot_esm_label,
                hgvs_change=hgvs_change,
            )
        if alt_aa == "*":
            consequence, impact_score, protein_score = self.stop_altering_subclass(
                "stop_gained", aa_pos, len(cds_seq) // 3
            )
            return self._coding_result(
                tx,
                consequence,
                impact_score,
                protein_score,
                aa_change,
                codon_change,
                protein_context_value=prot_context,
                protein_structure_model_value=prot_structure_model,
                protein_esm_value=prot_esm,
                protein_lm_value=prot_lm,
                protein_domain_value=domain_score,
                protein_domain_label=domain_label,
                protein_structure_label=prot_structure_label,
                protein_esm_label=prot_esm_label,
                hgvs_change=hgvs_change,
            )
        if ref_aa == "*":
            return self._coding_result(
                tx,
                "stop_lost_readthrough",
                0.95,
                0.92,
                aa_change,
                codon_change,
                protein_context_value=prot_context,
                protein_structure_model_value=prot_structure_model,
                protein_esm_value=prot_esm,
                protein_lm_value=prot_lm,
                protein_domain_value=domain_score,
                protein_domain_label=domain_label,
                protein_structure_label=prot_structure_label,
                protein_esm_label=prot_esm_label,
                hgvs_change=hgvs_change,
            )
        g_score = grantham_score(ref_aa, alt_aa)
        b_score = blosum62_score(ref_aa, alt_aa)
        prot_score = self.protein_substitution_score(
            ref_aa,
            alt_aa,
            g_score,
            b_score,
            codon_usage,
            prot_context,
            prot_structure,
            prot_structure_model,
            prot_esm,
            prot_lm,
            domain_score,
        )
        return self._coding_result(
            tx,
            "missense",
            0.58,
            prot_score,
            aa_change,
            codon_change,
            g_score,
            b_score,
            codon_usage,
            prot_context,
            prot_structure,
            prot_structure_model,
            prot_esm,
            prot_lm,
            domain_score,
            domain_label,
            prot_structure_label,
            prot_esm_label,
            hgvs_change,
        )

    @staticmethod
    def _coding_event_offset(tx: Transcript, pos: int, ref: str) -> Optional[int]:
        offset = tx.cds_offset(pos)
        if offset is not None:
            return offset
        if ref:
            return None
        prev_offset = tx.cds_offset(pos - 1)
        if prev_offset is not None:
            return prev_offset + 1
        next_offset = tx.cds_offset(pos + 1)
        if next_offset is not None:
            return max(0, next_offset - 1)
        return None

    def multi_codon_substitution_score(
        self,
        tx: Transcript,
        cds_seq: str,
        mutated: List[str],
        codon_starts: Iterable[int],
        changed_offsets: List[int],
    ) -> Dict[str, object]:
        events = []
        for codon_start in sorted(codon_starts):
            ref_codon = cds_seq[codon_start : codon_start + 3].upper()
            alt_codon = "".join(mutated[codon_start : codon_start + 3]).upper()
            if len(ref_codon) != 3 or len(alt_codon) != 3:
                return self._coding_result(tx, "cds_complex", 0.50, 0.35)
            ref_aa = translate_codon(ref_codon)
            alt_aa = translate_codon(alt_codon)
            aa_pos = codon_start // 3 + 1
            codon_usage = codon_usage_bias_score(ref_codon, alt_codon, self.codon_usage)
            protein_seq = self.protein_sequence(tx.tx_id)
            prot_context = protein_context_score(protein_seq, aa_pos)
            prot_structure = protein_structure_score(protein_seq, aa_pos, ref_aa, alt_aa)
            prot_structure_model, prot_structure_label = self.protein_structures.score(tx.tx_id, aa_pos, ref_aa, alt_aa)
            prot_esm, prot_esm_label = self.protein_lm_scores.score(tx.tx_id, aa_pos, ref_aa, alt_aa)
            prot_lm = self.protein_lm.delta_score(protein_seq, aa_pos, ref_aa, alt_aa)
            domain_score, domain_label = self.protein_domains.score(tx.tx_id, aa_pos)
            event = {
                "consequence": "synonymous",
                "impact_score": 0.05,
                "protein_score": max(0.02, 0.18 * codon_usage),
                "grantham_score": 0.0,
                "blosum_score": 0.0,
                "codon_usage_score": codon_usage,
                "protein_context_score": prot_context,
                "protein_structure_score": prot_structure,
                "protein_structure_model_score": prot_structure_model,
                "protein_esm_score": prot_esm,
                "protein_lm_score": prot_lm,
                "protein_domain_score": domain_score,
                "protein_domain_label": domain_label,
                "protein_structure_label": prot_structure_label,
                "protein_esm_label": prot_esm_label,
                "aa_change": f"p.{ref_aa}{aa_pos}{alt_aa}",
                "codon_change": f"{ref_codon}>{alt_codon}",
            }
            if aa_pos == 1 and ref_aa == "M" and alt_aa != "M":
                g_score = grantham_score(ref_aa, alt_aa)
                b_score = blosum62_score(ref_aa, alt_aa)
                prot_score = max(
                    0.88,
                    self.protein_substitution_score(
                        ref_aa,
                        alt_aa,
                        g_score,
                        b_score,
                        codon_usage,
                        prot_context,
                        prot_structure,
                        prot_structure_model,
                        prot_esm,
                        prot_lm,
                        domain_score,
                    ),
                )
                event.update(
                    {
                        "consequence": "start_lost",
                        "impact_score": 0.92,
                        "protein_score": prot_score,
                        "grantham_score": g_score,
                        "blosum_score": b_score,
                    }
                )
            elif ref_aa == "*" and alt_aa == "*":
                event.update({"consequence": "stop_retained", "impact_score": 0.04, "protein_score": 0.02})
            elif ref_aa == alt_aa:
                pass
            elif alt_aa == "*":
                consequence, impact_score, protein_score = self.stop_altering_subclass(
                    "stop_gained", aa_pos, len(cds_seq) // 3
                )
                event.update({"consequence": consequence, "impact_score": impact_score, "protein_score": protein_score})
            elif ref_aa == "*":
                event.update({"consequence": "stop_lost_readthrough", "impact_score": 0.95, "protein_score": 0.92})
            else:
                g_score = grantham_score(ref_aa, alt_aa)
                b_score = blosum62_score(ref_aa, alt_aa)
                prot_score = self.protein_substitution_score(
                    ref_aa,
                    alt_aa,
                    g_score,
                    b_score,
                    codon_usage,
                    prot_context,
                    prot_structure,
                    prot_structure_model,
                    prot_esm,
                    prot_lm,
                    domain_score,
                )
                event.update(
                    {
                        "consequence": "missense",
                        "impact_score": 0.58,
                        "protein_score": prot_score,
                        "grantham_score": g_score,
                        "blosum_score": b_score,
                    }
                )
            events.append(event)

        if not events:
            return self._coding_result(tx, "cds_complex", 0.50, 0.35)

        best = max(events, key=lambda item: (item["impact_score"], item["protein_score"]))
        aa_change = "|".join(str(item["aa_change"]) for item in events)
        codon_change = "|".join(str(item["codon_change"]) for item in events)
        return self._coding_result(
            tx,
            str(best["consequence"]),
            float(best["impact_score"]),
            max(float(item["protein_score"]) for item in events),
            aa_change,
            codon_change,
            max(float(item["grantham_score"]) for item in events),
            max(float(item["blosum_score"]) for item in events),
            max(float(item["codon_usage_score"]) for item in events),
            max(float(item["protein_context_score"]) for item in events),
            max(float(item["protein_structure_score"]) for item in events),
            max(float(item["protein_structure_model_score"]) for item in events),
            max(float(item["protein_esm_score"]) for item in events),
            max(float(item["protein_lm_score"]) for item in events),
            max(float(item["protein_domain_score"]) for item in events),
            self._best_domain_label(events),
            self._best_structure_label(events),
            self._best_esm_label(events),
            self._hgvs_like_substitution(tx, cds_seq, mutated, changed_offsets, aa_change),
        )

    @staticmethod
    def _best_domain_label(events: Sequence[Dict[str, object]]) -> str:
        best_score = 0.0
        best_label = "."
        for event in events:
            score = float(event.get("protein_domain_score", 0.0))
            label = str(event.get("protein_domain_label", "."))
            if score > best_score and label != ".":
                best_score = score
                best_label = label
        return best_label

    @staticmethod
    def _best_structure_label(events: Sequence[Dict[str, object]]) -> str:
        best_score = 0.0
        best_label = "."
        for event in events:
            score = float(event.get("protein_structure_model_score", 0.0))
            label = str(event.get("protein_structure_label", "."))
            if score > best_score and label != ".":
                best_score = score
                best_label = label
        return best_label

    @staticmethod
    def _best_esm_label(events: Sequence[Dict[str, object]]) -> str:
        best_score = 0.0
        best_label = "."
        for event in events:
            score = float(event.get("protein_esm_score", 0.0))
            label = str(event.get("protein_esm_label", "."))
            if score > best_score and label != ".":
                best_score = score
                best_label = label
        return best_label

    @staticmethod
    def stop_altering_subclass(event: str, aa_pos: int, coding_codons: int) -> Tuple[str, float, float]:
        if event == "stop_gained":
            if coding_codons < 10:
                return "stop_gained", 1.00, 1.00
            remaining_codons = max(0, coding_codons - aa_pos)
            terminal_window = max(2, math.ceil(coding_codons * 0.10))
            if remaining_codons <= terminal_window:
                return "stop_gained_terminal", 0.82, 0.78
            return "stop_gained_early", 1.00, 1.00
        if event == "stop_lost":
            return "stop_lost_readthrough", 0.95, 0.92
        return event, CONSEQUENCE_BASE_SCORE.get(event, 0.08), 0.0

    @staticmethod
    def _hgvs_like_substitution(
        tx: Transcript,
        cds_seq: str,
        mutated: List[str],
        changed_offsets: List[int],
        aa_change: str,
    ) -> str:
        cdna_changes = []
        for idx in sorted(set(changed_offsets)):
            if idx >= len(cds_seq) or idx >= len(mutated):
                continue
            ref_base = cds_seq[idx].upper()
            alt_base = mutated[idx].upper()
            if ref_base == alt_base:
                continue
            cdna_changes.append(f"c.{idx + 1}{ref_base}>{alt_base}")
        if not cdna_changes:
            return "."
        return f"{tx.tx_id}:{'|'.join(cdna_changes)}:{aa_change}"

    @staticmethod
    def _hgvs_like_indel(tx: Transcript, offset: int, ref: str, alt: str) -> str:
        start = offset + 1
        end = offset + max(1, len(ref))
        if tx.strand == "-":
            ref_tx = ref.translate(DNA_COMP).upper()[::-1]
            alt_tx = alt.translate(DNA_COMP).upper()[::-1]
        else:
            ref_tx = ref.upper()
            alt_tx = alt.upper()
        if not ref_tx and alt_tx:
            left = max(1, start - 1)
            event = f"c.{left}_{start}ins{alt_tx}"
        elif ref_tx and not alt_tx:
            event = f"c.{start}_{end}del"
        elif len(ref_tx) > len(alt_tx):
            event = f"c.{start}_{end}delins{alt_tx}"
        elif len(alt_tx) > len(ref_tx):
            event = f"c.{start}_{end}delins{alt_tx}"
        else:
            event = f"c.{start}_{end}delins{alt_tx}"
        return f"{tx.tx_id}:{event}:p.?"

    def protein_substitution_score(
        self,
        ref_aa: str,
        alt_aa: str,
        g_score: float,
        b_score: float,
        codon_usage: float,
        prot_context: float,
        prot_structure: float,
        prot_structure_model: float,
        prot_esm: float,
        prot_lm: float,
        protein_domain: float,
    ) -> float:
        return weighted_mean(
            {
                "physchem": aa_property_score(ref_aa, alt_aa),
                "grantham": g_score,
                "blosum": b_score,
                "codon_usage": codon_usage,
                "protein_context": prot_context,
                "protein_structure": prot_structure,
                "protein_structure_model": prot_structure_model,
                "protein_esm": prot_esm,
                "protein_lm": prot_lm,
                "protein_domain": protein_domain,
            },
            self.config.protein_weights,
        )

    def _coding_result(
        self,
        tx: Transcript,
        consequence: str,
        impact_score: float,
        protein_score: float,
        aa_change: str = ".",
        codon_change: str = ".",
        grantham_value: float = 0.0,
        blosum_value: float = 0.0,
        codon_usage_value: float = 0.0,
        protein_context_value: float = 0.0,
        protein_structure_value: float = 0.0,
        protein_structure_model_value: float = 0.0,
        protein_esm_value: float = 0.0,
        protein_lm_value: float = 0.0,
        protein_domain_value: float = 0.0,
        protein_domain_label: str = ".",
        protein_structure_label: str = ".",
        protein_esm_label: str = ".",
        hgvs_change: str = ".",
    ) -> Dict[str, object]:
        return {
            "consequence": consequence,
            "impact_score": impact_score,
            "protein_score": protein_score,
            "grantham_score": grantham_value,
            "blosum_score": blosum_value,
            "codon_usage_score": codon_usage_value,
            "protein_context_score": protein_context_value,
            "protein_structure_score": protein_structure_value,
            "protein_structure_model_score": protein_structure_model_value,
            "protein_esm_score": protein_esm_value,
            "protein_lm_score": protein_lm_value,
            "protein_domain_score": protein_domain_value,
            "protein_domain_label": protein_domain_label,
            "protein_structure_label": protein_structure_label,
            "protein_esm_label": protein_esm_label,
            "splice_score": 0.0,
            "splice_motif_score": 0.0,
            "splice_pwm_score": 0.0,
            "splice_maxent_score": 0.0,
            "splice_aux_score": 0.0,
            "splice_ese_score": 0.0,
            "utr_score": 0.0,
            "rnafold_score": 0.0,
            "mirna_score": 0.0,
            "mirna_label": ".",
            "promoter_score": 0.0,
            "gene_id": tx.gene_id,
            "tx_id": tx.tx_id,
            "aa_change": aa_change,
            "codon_change": codon_change,
            "hgvs_change": hgvs_change,
            "cds_length": tx.cds_length(),
            "transcript_length": tx.transcript_length(),
        }

    def combine_scores(
        self,
        impact: float,
        protein: float,
        splice: float,
        sequence: float,
        cohort: float,
        n_samples: int,
    ) -> float:
        cohort_weight = (
            self.config.cohort_weights["single_sample"]
            if n_samples <= 1
            else (
                self.config.cohort_weights["small_cohort"]
                if n_samples < 10
                else self.config.cohort_weights["large_cohort"]
            )
        )
        raw_weights = dict(self.config.score_weights)
        raw_weights["cohort"] = cohort_weight
        weight_sum = sum(raw_weights.values())
        weighted = (
            raw_weights["impact"] * impact
            + raw_weights["protein"] * protein
            + raw_weights["splice"] * splice
            + raw_weights["sequence"] * sequence
            + raw_weights["cohort"] * cohort
        ) / weight_sum
        max_signal = max(impact, protein, splice)
        blend_sum = self.config.max_signal_weight + self.config.weighted_signal_weight
        return clamp(
            (
                self.config.max_signal_weight * max_signal
                + self.config.weighted_signal_weight * weighted
            )
            / blend_sum
        )


def fmt_float(value: float) -> str:
    return f"{value:.4f}"


def sanitize_info_value(value: str) -> str:
    if value in {"", None}:  # type: ignore[comparison-overlap]
        return "."
    return re.sub(r"[\s,;=]", "_", str(value))


INFO_HEADERS = [
    '##INFO=<ID=DNP_SCORE,Number=A,Type=Float,Description="DeNovoPath weighted summary deleteriousness score, 0-1">',
    '##INFO=<ID=DNP_IMPACT,Number=A,Type=Float,Description="DeNovoPath annotation/consequence rule score, 0-1">',
    '##INFO=<ID=DNP_PROT,Number=A,Type=Float,Description="DeNovoPath coding protein physicochemical score, 0-1">',
    '##INFO=<ID=DNP_GRANTHAM,Number=A,Type=Float,Description="DeNovoPath normalized Grantham amino-acid substitution distance score, 0-1">',
    '##INFO=<ID=DNP_BLOSUM,Number=A,Type=Float,Description="DeNovoPath normalized BLOSUM62 amino-acid substitution severity score, 0-1">',
    '##INFO=<ID=DNP_CODONUSE,Number=A,Type=Float,Description="DeNovoPath CDS-derived codon usage bias change score, 0-1">',
    '##INFO=<ID=DNP_PROTCTX,Number=A,Type=Float,Description="DeNovoPath protein low-complexity/disorder-context score around the affected amino acid, 0-1">',
    '##INFO=<ID=DNP_STRUCT,Number=A,Type=Float,Description="DeNovoPath heuristic protein structural-context score for secondary-structure shift, solvent exposure, and buried charged-change risk, 0-1">',
    '##INFO=<ID=DNP_AFSTRUCT,Number=A,Type=Float,Description="DeNovoPath optional AlphaFold/ESMFold-derived per-residue structural-context score from --protein-structures, 0-1">',
    '##INFO=<ID=DNP_ESM,Number=A,Type=Float,Description="DeNovoPath optional precomputed ESM-2 protein language-model substitution score from --protein-lm-scores, 0-1">',
    '##INFO=<ID=DNP_PROTLM,Number=A,Type=Float,Description="DeNovoPath protein k-mer language-model ref-vs-alt likelihood delta proxy score, 0-1">',
    '##INFO=<ID=DNP_DOMAIN,Number=A,Type=Float,Description="DeNovoPath protein domain overlap/proximity score from optional --protein-domains annotation, 0-1">',
    '##INFO=<ID=DNP_SPLICE,Number=A,Type=Float,Description="DeNovoPath splice-boundary proximity score, 0-1">',
    '##INFO=<ID=DNP_SPLICE_MOTIF,Number=A,Type=Float,Description="DeNovoPath canonical splice donor/acceptor motif disruption score, 0-1">',
    '##INFO=<ID=DNP_SPLICE_PWM,Number=A,Type=Float,Description="DeNovoPath species-specific splice donor/acceptor PWM disruption score trained from GFF/reference, 0-1">',
    '##INFO=<ID=DNP_SPLICE_MAXENT,Number=A,Type=Float,Description="DeNovoPath MaxEntScan-like heuristic donor/acceptor splice-window log-odds disruption score, 0-1">',
    '##INFO=<ID=DNP_SPLICE_AUX,Number=A,Type=Float,Description="DeNovoPath auxiliary splice score for branch-point and polypyrimidine-tract disruption or creation, 0-1">',
    '##INFO=<ID=DNP_SPLICE_ESE,Number=A,Type=Float,Description="DeNovoPath exonic splicing enhancer/silencer motif disruption or creation score, 0-1">',
    '##INFO=<ID=DNP_UTR,Number=A,Type=Float,Description="DeNovoPath UTR regulatory-context score including start/Kozak proximity and polyadenylation motif changes, 0-1">',
    '##INFO=<ID=DNP_RNAFOLD,Number=A,Type=Float,Description="DeNovoPath heuristic local RNA folding delta-G proxy score for UTR variants, 0-1">',
    '##INFO=<ID=DNP_MIRNA,Number=A,Type=Float,Description="DeNovoPath optional miRNA target/seed overlap disruption score from --mirna-sites, 0-1">',
    '##INFO=<ID=DNP_PROM,Number=A,Type=Float,Description="DeNovoPath promoter core motif score including TATA/CAAT/GC-box disruption or creation, 0-1">',
    '##INFO=<ID=DNP_SEQ,Number=A,Type=Float,Description="DeNovoPath local sequence-context disruption score, 0-1">',
    '##INFO=<ID=DNP_KMER,Number=A,Type=Float,Description="DeNovoPath local trinucleotide/k-mer context change score, 0-1">',
    '##INFO=<ID=DNP_REPEAT,Number=A,Type=Float,Description="DeNovoPath reference-only repeat/low-complexity proxy score for low-mappability context, 0-1">',
    '##INFO=<ID=DNP_MUTCTX,Number=A,Type=Float,Description="DeNovoPath mutation 96-context heuristic score, 0-1">',
    '##INFO=<ID=DNP_DNALM,Number=A,Type=Float,Description="DeNovoPath reference-trained k-mer DNA language-model likelihood delta score, 0-1">',
    '##INFO=<ID=DNP_COHORT,Number=A,Type=Float,Description="DeNovoPath genotype cohort rarity/carrier-pattern score, 0-1; weakly weighted for tiny cohorts">',
    '##INFO=<ID=DNP_HWE,Number=A,Type=Float,Description="DeNovoPath Hardy-Weinberg genotype deviation score from REF/current-ALT diploid GT classes for cohorts with >=10 samples, 0-1">',
    '##INFO=<ID=DNP_HETOBS,Number=A,Type=Float,Description="DeNovoPath observed heterozygosity from REF/current-ALT diploid GT classes for cohorts with >=10 samples, 0-1">',
    '##INFO=<ID=DNP_HETEXP,Number=A,Type=Float,Description="DeNovoPath expected heterozygosity 2pq from GT-derived ALT frequency for cohorts with >=10 samples, 0-1">',
    '##INFO=<ID=DNP_HETDEV,Number=A,Type=Float,Description="DeNovoPath observed-vs-expected heterozygosity deviation score for cohorts with >=10 samples, 0-1">',
    '##INFO=<ID=DNP_FIS,Number=A,Type=Float,Description="DeNovoPath bounded inbreeding coefficient estimate 1-Hobs/Hexp for cohorts with >=10 samples, -1 to 1">',
    '##INFO=<ID=DNP_FST,Number=A,Type=Float,Description="DeNovoPath group differentiation Fst-style score from optional sample_info groups, 0-1">',
    '##INFO=<ID=DNP_CASECTRL,Number=A,Type=Float,Description="DeNovoPath case-control ALT frequency difference score from optional sample_info phenotypes, 0-1">',
    '##INFO=<ID=DNP_PI,Number=A,Type=Float,Description="DeNovoPath fixed-window nucleotide diversity pi from GT-derived allele frequencies for cohorts with >=10 samples">',
    '##INFO=<ID=DNP_THETA,Number=A,Type=Float,Description="DeNovoPath fixed-window Watterson theta from segregating sites for cohorts with >=10 samples">',
    '##INFO=<ID=DNP_TAJD,Number=A,Type=Float,Description="DeNovoPath fixed-window Tajima D from pi and Watterson theta for cohorts with >=10 samples">',
    '##INFO=<ID=DNP_LD,Number=A,Type=Float,Description="DeNovoPath fixed-window adjacent-variant genotype dosage LD proxy as max r2 for cohorts with >=10 samples, 0-1">',
    '##INFO=<ID=DNP_HAP,Number=A,Type=Float,Description="DeNovoPath fixed-window haplotype diversity proxy as 1-mean adjacent LD for cohorts with >=10 samples, 0-1">',
    '##INFO=<ID=DNP_GENELOF,Number=A,Type=Float,Description="DeNovoPath gene-level observed/expected LoF proxy from current cohort VCF for cohorts with >=10 samples">',
    '##INFO=<ID=DNP_GENEMIS,Number=A,Type=Float,Description="DeNovoPath gene-level observed/expected missense proxy from current cohort VCF for cohorts with >=10 samples">',
    '##INFO=<ID=DNP_GENECON,Number=A,Type=Float,Description="DeNovoPath gene-level LoF/missense constraint proxy from current cohort VCF for cohorts with >=10 samples, 0-1">',
    '##INFO=<ID=DNP_QC,Number=A,Type=Float,Description="DeNovoPath genotype/call quality support score for each ALT, 0-1">',
    '##INFO=<ID=DNP_CONF,Number=A,Type=Float,Description="DeNovoPath prediction confidence combining DNP_QC and score separation, 0-1">',
    '##INFO=<ID=DNP_ML,Number=A,Type=Float,Description="DeNovoPath optional JSON ML-model deleteriousness score from deterministic features, 0-1">',
    '##INFO=<ID=DNP_CAL,Number=A,Type=Float,Description="DeNovoPath optional calibrated ML score after model-supplied logit scaling, 0-1">',
    '##INFO=<ID=DNP_UNCERT,Number=A,Type=Float,Description="DeNovoPath optional ML uncertainty proxy combining calibrated-score ambiguity and OOD score, 0-1">',
    '##INFO=<ID=DNP_OOD,Number=A,Type=Float,Description="DeNovoPath optional ML out-of-distribution proxy from standardized feature distance, 0-1">',
    '##INFO=<ID=DNP_LEVEL,Number=A,Type=String,Description="DeNovoPath summary impact level: HIGH, MODERATE, LOW, or MINIMAL">',
    '##INFO=<ID=DNP_CONSEQ,Number=A,Type=String,Description="DeNovoPath strongest predicted consequence per ALT">',
    '##INFO=<ID=DNP_GENE,Number=A,Type=String,Description="Gene ID supporting DNP_CONSEQ per ALT">',
    '##INFO=<ID=DNP_TX,Number=A,Type=String,Description="Transcript ID supporting DNP_CONSEQ per ALT">',
    '##INFO=<ID=DNP_ALLTX,Number=A,Type=String,Description="All transcript consequences considered per ALT as gene:transcript:consequence entries separated by |">',
    '##INFO=<ID=DNP_AA,Number=A,Type=String,Description="Predicted amino-acid change per ALT when available">',
    '##INFO=<ID=DNP_CODON,Number=A,Type=String,Description="Predicted codon change per ALT when available">',
    '##INFO=<ID=DNP_HGVS,Number=A,Type=String,Description="HGVS-like CDS/protein change per ALT when available; heuristic, not fully normalized">',
    '##INFO=<ID=DNP_NORM,Number=A,Type=String,Description="Minimal allele representation per ALT as POS:REF>ALT after trimming shared prefix/suffix; empty allele is -">',
    '##INFO=<ID=DNP_DOMID,Number=A,Type=String,Description="Protein domain annotation supporting DNP_DOMAIN per ALT when optional --protein-domains is supplied">',
    '##INFO=<ID=DNP_AFID,Number=A,Type=String,Description="Protein structure annotation supporting DNP_AFSTRUCT per ALT when optional --protein-structures is supplied">',
    '##INFO=<ID=DNP_ESMID,Number=A,Type=String,Description="Precomputed ESM-2 label supporting DNP_ESM per ALT when optional --protein-lm-scores is supplied">',
    '##INFO=<ID=DNP_MIRID,Number=A,Type=String,Description="miRNA target/seed annotation supporting DNP_MIRNA per ALT when optional --mirna-sites is supplied">',
    '##INFO=<ID=DNP_FEATIMP,Number=A,Type=String,Description="Top per-ALT DeNovoPath ML feature contributions as feature:signed_contribution entries separated by |">',
    '##INFO=<ID=DNP_96CTX,Number=A,Type=String,Description="Pyrimidine-oriented 96-context mutation label per ALT, for example A[C>T]G">',
    '##INFO=<ID=DNP_SUBAF,Number=A,Type=String,Description="Subpopulation ALT AF from optional sample_info groups as group:AF entries separated by |">',
    '##INFO=<ID=DNP_PRIVATE,Number=A,Type=String,Description="Group private/shared ALT status from optional sample_info groups">',
    '##INFO=<ID=DNP_CASEAF,Number=A,Type=String,Description="Case/control ALT AF from optional sample_info phenotypes as case:AF|control:AF">',
    '##INFO=<ID=DNP_AC,Number=A,Type=Integer,Description="ALT allele count recomputed from GT fields">',
    '##INFO=<ID=DNP_AN,Number=A,Type=Integer,Description="Allele number recomputed from GT fields">',
    '##INFO=<ID=DNP_AF,Number=A,Type=Float,Description="ALT allele frequency recomputed from GT fields">',
    '##INFO=<ID=DNP_MAFBIN,Number=A,Type=String,Description="Sample-count-aware minor allele frequency bin per ALT; true MAF bins require >=10 samples">',
    '##INFO=<ID=DNP_CARR,Number=A,Type=Integer,Description="Number of samples carrying each ALT allele">',
    '##INFO=<ID=DNP_HET,Number=A,Type=Integer,Description="Number of REF/current-ALT heterozygous samples for each ALT allele; genotypes containing another ALT are excluded from this class">',
    '##INFO=<ID=DNP_HOMALT,Number=A,Type=Integer,Description="Number of current-ALT homozygous samples for each ALT allele; genotypes containing another ALT are excluded from HWE genotype classes">',
    '##INFO=<ID=DNP_MISS,Number=A,Type=Integer,Description="Number of samples missing GT for each ALT allele">',
]


def append_scores_to_info(info: str, scores: List[AltScore]) -> str:
    additions = {
        "DNP_SCORE": ",".join(fmt_float(item.score) for item in scores),
        "DNP_IMPACT": ",".join(fmt_float(item.impact_score) for item in scores),
        "DNP_PROT": ",".join(fmt_float(item.protein_score) for item in scores),
        "DNP_GRANTHAM": ",".join(fmt_float(item.grantham_score) for item in scores),
        "DNP_BLOSUM": ",".join(fmt_float(item.blosum_score) for item in scores),
        "DNP_CODONUSE": ",".join(fmt_float(item.codon_usage_score) for item in scores),
        "DNP_PROTCTX": ",".join(fmt_float(item.protein_context_score) for item in scores),
        "DNP_STRUCT": ",".join(fmt_float(item.protein_structure_score) for item in scores),
        "DNP_AFSTRUCT": ",".join(fmt_float(item.protein_structure_model_score) for item in scores),
        "DNP_ESM": ",".join(fmt_float(item.protein_esm_score) for item in scores),
        "DNP_PROTLM": ",".join(fmt_float(item.protein_lm_score) for item in scores),
        "DNP_DOMAIN": ",".join(fmt_float(item.protein_domain_score) for item in scores),
        "DNP_SPLICE": ",".join(fmt_float(item.splice_score) for item in scores),
        "DNP_SPLICE_MOTIF": ",".join(fmt_float(item.splice_motif_score) for item in scores),
        "DNP_SPLICE_PWM": ",".join(fmt_float(item.splice_pwm_score) for item in scores),
        "DNP_SPLICE_MAXENT": ",".join(fmt_float(item.splice_maxent_score) for item in scores),
        "DNP_SPLICE_AUX": ",".join(fmt_float(item.splice_aux_score) for item in scores),
        "DNP_SPLICE_ESE": ",".join(fmt_float(item.splice_ese_score) for item in scores),
        "DNP_UTR": ",".join(fmt_float(item.utr_score) for item in scores),
        "DNP_RNAFOLD": ",".join(fmt_float(item.rnafold_score) for item in scores),
        "DNP_MIRNA": ",".join(fmt_float(item.mirna_score) for item in scores),
        "DNP_PROM": ",".join(fmt_float(item.promoter_score) for item in scores),
        "DNP_SEQ": ",".join(fmt_float(item.sequence_score) for item in scores),
        "DNP_KMER": ",".join(fmt_float(item.kmer_score) for item in scores),
        "DNP_REPEAT": ",".join(fmt_float(item.repeat_score) for item in scores),
        "DNP_MUTCTX": ",".join(fmt_float(item.mutation_context_score) for item in scores),
        "DNP_DNALM": ",".join(fmt_float(item.dna_lm_score) for item in scores),
        "DNP_COHORT": ",".join(fmt_float(item.cohort_score) for item in scores),
        "DNP_HWE": ",".join(fmt_float(item.hwe_score) for item in scores),
        "DNP_HETOBS": ",".join(fmt_float(item.heterozygosity_observed) for item in scores),
        "DNP_HETEXP": ",".join(fmt_float(item.heterozygosity_expected) for item in scores),
        "DNP_HETDEV": ",".join(fmt_float(item.heterozygosity_deviation_score) for item in scores),
        "DNP_FIS": ",".join(fmt_float(item.inbreeding_coefficient) for item in scores),
        "DNP_FST": ",".join(fmt_float(item.fst_score) for item in scores),
        "DNP_CASECTRL": ",".join(fmt_float(item.case_control_score) for item in scores),
        "DNP_PI": ",".join(fmt_float(item.window_pi) for item in scores),
        "DNP_THETA": ",".join(fmt_float(item.window_theta) for item in scores),
        "DNP_TAJD": ",".join(fmt_float(item.window_tajima_d) for item in scores),
        "DNP_LD": ",".join(fmt_float(item.window_ld) for item in scores),
        "DNP_HAP": ",".join(fmt_float(item.window_haplotype) for item in scores),
        "DNP_GENELOF": ",".join(fmt_float(item.gene_lof_oe) for item in scores),
        "DNP_GENEMIS": ",".join(fmt_float(item.gene_missense_oe) for item in scores),
        "DNP_GENECON": ",".join(fmt_float(item.gene_constraint_score) for item in scores),
        "DNP_QC": ",".join(fmt_float(item.qc_score) for item in scores),
        "DNP_CONF": ",".join(fmt_float(item.confidence_score) for item in scores),
        "DNP_ML": ",".join(fmt_float(item.ml_score) for item in scores),
        "DNP_CAL": ",".join(fmt_float(item.calibrated_score) for item in scores),
        "DNP_UNCERT": ",".join(fmt_float(item.uncertainty_score) for item in scores),
        "DNP_OOD": ",".join(fmt_float(item.ood_score) for item in scores),
        "DNP_LEVEL": ",".join(sanitize_info_value(item.level) for item in scores),
        "DNP_CONSEQ": ",".join(sanitize_info_value(item.consequence) for item in scores),
        "DNP_GENE": ",".join(sanitize_info_value(item.gene_id) for item in scores),
        "DNP_TX": ",".join(sanitize_info_value(item.tx_id) for item in scores),
        "DNP_ALLTX": ",".join(sanitize_info_value(item.all_transcripts) for item in scores),
        "DNP_AA": ",".join(sanitize_info_value(item.aa_change) for item in scores),
        "DNP_CODON": ",".join(sanitize_info_value(item.codon_change) for item in scores),
        "DNP_HGVS": ",".join(sanitize_info_value(item.hgvs_change) for item in scores),
        "DNP_NORM": ",".join(sanitize_info_value(item.normalized_variant) for item in scores),
        "DNP_DOMID": ",".join(sanitize_info_value(item.protein_domain_label) for item in scores),
        "DNP_AFID": ",".join(sanitize_info_value(item.protein_structure_label) for item in scores),
        "DNP_ESMID": ",".join(sanitize_info_value(item.protein_esm_label) for item in scores),
        "DNP_MIRID": ",".join(sanitize_info_value(item.mirna_label) for item in scores),
        "DNP_FEATIMP": ",".join(sanitize_info_value(item.feature_importance) for item in scores),
        "DNP_96CTX": ",".join(sanitize_info_value(item.mutation_context) for item in scores),
        "DNP_SUBAF": ",".join(sanitize_info_value(item.group_af) for item in scores),
        "DNP_PRIVATE": ",".join(sanitize_info_value(item.private_shared) for item in scores),
        "DNP_CASEAF": ",".join(sanitize_info_value(item.case_control_af) for item in scores),
        "DNP_AC": ",".join(str(item.ac) for item in scores),
        "DNP_AN": ",".join(str(item.an) for item in scores),
        "DNP_AF": ",".join(fmt_float(item.af) for item in scores),
        "DNP_MAFBIN": ",".join(sanitize_info_value(item.maf_bin) for item in scores),
        "DNP_CARR": ",".join(str(item.carriers) for item in scores),
        "DNP_HET": ",".join(str(item.n_het) for item in scores),
        "DNP_HOMALT": ",".join(str(item.n_hom_alt) for item in scores),
        "DNP_MISS": ",".join(str(item.n_missing) for item in scores),
    }
    base = [] if info in {".", ""} else [part for part in info.split(";") if not part.startswith("DNP_")]
    base.extend(f"{key}={value}" for key, value in additions.items())
    return ";".join(base) if base else "."


def open_text(path: str, mode: str):
    if path.endswith(".gz"):
        if any(flag in mode for flag in ("w", "a", "x")):
            return gzip.open(path, mode + "t", compresslevel=FAST_GZIP_COMPRESSLEVEL)
        return gzip.open(path, mode + "t")
    return open(path, mode)


def peak_rss_mb() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    rss = float(getattr(usage, "ru_maxrss", 0.0))
    if sys.platform == "darwin":
        return rss / (1024.0 * 1024.0)
    return rss / 1024.0


def maybe_index_vcf(path: str, mode: str = "auto") -> Dict[str, object]:
    mode = mode or "auto"
    result: Dict[str, object] = {
        "requested": mode,
        "status": "skipped",
        "index_path": None,
        "message": "",
    }
    if mode == "never":
        result["status"] = "disabled"
        result["message"] = "indexing disabled"
        return result
    if not path.endswith(".vcf.gz"):
        result["message"] = "output is not .vcf.gz"
        if mode == "always":
            raise RuntimeError(str(result["message"]))
        return result
    bgzip = shutil.which("bgzip")
    tabix = shutil.which("tabix")
    if not bgzip or not tabix:
        result["message"] = "bgzip and/or tabix not found"
        if mode == "always":
            raise RuntimeError(str(result["message"]))
        return result

    tmp_plain = None
    tmp_bgz = f"{path}.bgzip.tmp"
    try:
        with tempfile.NamedTemporaryFile("wt", delete=False, suffix=".vcf") as plain:
            tmp_plain = plain.name
            with gzip.open(path, "rt") as source:
                shutil.copyfileobj(source, plain)
        with open(tmp_bgz, "wb") as bgz_out:
            subprocess.run([bgzip, "-c", tmp_plain], stdout=bgz_out, check=True)
        os.replace(tmp_bgz, path)
        subprocess.run([tabix, "-f", "-p", "vcf", path], check=True)
        result["status"] = "indexed"
        result["index_path"] = f"{path}.tbi"
        result["message"] = "bgzip recompressed and tabix indexed"
        return result
    except Exception as exc:
        result["status"] = "failed"
        result["message"] = str(exc)
        if mode == "always":
            raise
        return result
    finally:
        if tmp_plain and os.path.exists(tmp_plain):
            os.unlink(tmp_plain)
        if os.path.exists(tmp_bgz):
            os.unlink(tmp_bgz)


def build_run_summary(
    args: argparse.Namespace,
    n_samples: int,
    sample_info: Optional[SampleInfo],
    pop_window_size: int,
    gene_constraint_stats: Dict[str, "GeneConstraintStat"],
    input_records: int,
    records: int,
    skipped_region: int,
    alt_records: int,
    index_result: Dict[str, object],
    prescan_seconds: float,
    scoring_seconds: float,
    index_seconds: float,
    elapsed_seconds: float,
    level_counts: Counter[str],
    consequence_counts: Counter[str],
    protein_domain_records: int = 0,
    protein_structure_records: int = 0,
    protein_lm_score_records: int = 0,
    mirna_site_records: int = 0,
    ml_model_loaded: bool = False,
    ml_model_type: str = "",
    ml_model_features: Sequence[str] = (),
    ml_feature_importance: Sequence[Tuple[str, float]] = (),
    ref_match_records: int = 0,
    ref_mismatch_records: int = 0,
    ref_unchecked_records: int = 0,
    ref_mismatch_examples: Sequence[str] = (),
    population_windows_enabled: bool = True,
    gene_constraint_enabled: bool = True,
) -> Dict[str, object]:
    records_per_second = records / scoring_seconds if scoring_seconds > 0 else 0.0
    alt_alleles_per_second = alt_records / scoring_seconds if scoring_seconds > 0 else 0.0
    has_sample_info = sample_info is not None
    has_phenotype = bool(sample_info and sample_info.has_phenotype())
    return {
        "input_vcf": args.vcf,
        "output_vcf": args.output,
        "config": args.config,
        "sample_info": getattr(args, "sample_info", None),
        "protein_domains": getattr(args, "protein_domains", None),
        "protein_domain_records": protein_domain_records,
        "protein_structures": getattr(args, "protein_structures", None),
        "protein_structure_records": protein_structure_records,
        "protein_lm_scores": getattr(args, "protein_lm_scores", None),
        "protein_lm_score_records": protein_lm_score_records,
        "mirna_sites": getattr(args, "mirna_sites", None),
        "mirna_site_records": mirna_site_records,
        "ml_model": getattr(args, "ml_model", None),
        "ml_model_loaded": ml_model_loaded,
        "ml_model_type": ml_model_type,
        "ml_model_features": list(ml_model_features),
        "ml_model_top_features": [
            {"feature": feature, "importance": round(float(importance), 6)}
            for feature, importance in list(ml_feature_importance)[:10]
        ],
        "population_window_size": pop_window_size,
        "population_windows_enabled": population_windows_enabled,
        "gene_constraint_enabled": gene_constraint_enabled,
        "gene_constraint_genes": len(gene_constraint_stats),
        "regions": args.region or [],
        "active_methods": active_methods(
            n_samples,
            has_sample_info,
            has_phenotype,
            has_protein_domains=protein_domain_records > 0,
            has_protein_structures=protein_structure_records > 0,
            has_protein_esm_scores=protein_lm_score_records > 0,
            has_mirna_sites=mirna_site_records > 0,
            has_ml_model=ml_model_loaded,
            has_population_windows=population_windows_enabled,
            has_gene_constraint=gene_constraint_enabled,
        ),
        "method_source_types": METHOD_SOURCE_TYPES,
        "sample_count_gating": cohort_gating_note(
            n_samples,
            population_windows_enabled,
            gene_constraint_enabled,
        ),
        "n_samples": n_samples,
        "sample_info_groups": sorted(set(sample_info.groups.values())) if sample_info else [],
        "sample_info_has_case_control": has_phenotype,
        "validation": {
            "reference": {
                "checked_records": ref_match_records + ref_mismatch_records,
                "matching_records": ref_match_records,
                "mismatch_records": ref_mismatch_records,
                "unchecked_records": ref_unchecked_records,
                "mismatch_examples": list(ref_mismatch_examples),
            }
        },
        "input_records_seen": input_records,
        "records_scored": records,
        "records_skipped_by_region": skipped_region,
        "alt_alleles_scored": alt_records,
        "index": index_result,
        "benchmark": {
            "elapsed_seconds": round(elapsed_seconds, 6),
            "prescan_seconds": round(prescan_seconds, 6),
            "scoring_seconds": round(scoring_seconds, 6),
            "index_seconds": round(index_seconds, 6),
            "records_per_second": round(records_per_second, 3),
            "alt_alleles_per_second": round(alt_alleles_per_second, 3),
            "peak_rss_mb": round(peak_rss_mb(), 3),
        },
        "performance_settings": {
            "gzip_compresslevel": FAST_GZIP_COMPRESSLEVEL,
            "max_codon_usage_training_codons": DEFAULT_MAX_CODON_USAGE_TRAINING_CODONS,
            "max_protein_lm_training_residues": DEFAULT_MAX_PROTEIN_LM_TRAINING_RESIDUES,
            "max_splice_pwm_training_introns": DEFAULT_MAX_SPLICE_PWM_TRAINING_INTRONS,
        },
        "level_counts": dict(sorted(level_counts.items())),
        "consequence_counts": dict(consequence_counts.most_common()),
    }


def html_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ", ".join(html_value(item) for item in value)
    if isinstance(value, dict):
        return ", ".join(f"{html_value(key)}: {html_value(val)}" for key, val in value.items())
    return html.escape(str(value), quote=True)


def html_table(title: str, rows: Sequence[Tuple[object, object]], empty_message: str = "None") -> str:
    body = []
    for key, value in rows:
        body.append(f"<tr><th>{html_value(key)}</th><td>{html_value(value)}</td></tr>")
    if not body:
        body.append(f"<tr><td colspan=\"2\">{html_value(empty_message)}</td></tr>")
    return (
        f"<section><h2>{html_value(title)}</h2>"
        "<table><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


CHART_COLORS = [
    "#0f766e",
    "#2563eb",
    "#b45309",
    "#be123c",
    "#6d28d9",
    "#047857",
    "#7c2d12",
    "#475569",
    "#0891b2",
    "#a16207",
]

LEVEL_COLORS = {
    "HIGH": "#be123c",
    "MODERATE": "#b45309",
    "LOW": "#2563eb",
    "MINIMAL": "#64748b",
}


def html_int(value: object) -> int:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0


def html_float(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0.0


def fmt_number(value: object) -> str:
    number = html_float(value)
    if abs(number - int(number)) < 0.000001:
        return f"{int(number):,}"
    return f"{number:,.3f}".rstrip("0").rstrip(".")


def fmt_duration(seconds: object) -> str:
    seconds_f = html_float(seconds)
    if seconds_f >= 3600:
        return f"{seconds_f / 3600:.2f} h"
    if seconds_f >= 60:
        return f"{seconds_f / 60:.1f} min"
    return f"{seconds_f:.2f} s"


def html_metric_card(label: str, value: object, detail: object = "") -> str:
    return (
        "<div class=\"metric-card\">"
        f"<div class=\"metric-label\">{html_value(label)}</div>"
        f"<div class=\"metric-value\">{html_value(value)}</div>"
        f"<div class=\"metric-detail\">{html_value(detail)}</div>"
        "</div>"
    )


def normalized_count_items(counts: object, preferred_order: Sequence[str] = ()) -> List[Tuple[str, int]]:
    if not isinstance(counts, dict):
        return []
    items: List[Tuple[str, int]] = []
    used = set()
    for key in preferred_order:
        if key in counts:
            items.append((key, html_int(counts.get(key))))
            used.add(key)
    for key, value in sorted(counts.items(), key=lambda item: html_int(item[1]), reverse=True):
        key_s = str(key)
        if key_s not in used:
            items.append((key_s, html_int(value)))
    return [(key, value) for key, value in items if value > 0]


def html_donut_chart(title: str, items: Sequence[Tuple[str, int]], color_map: Optional[Dict[str, str]] = None) -> str:
    total = sum(value for _key, value in items)
    if total <= 0:
        return (
            f"<section class=\"chart-panel\"><h2>{html_value(title)}</h2>"
            "<p class=\"empty-note\">No chart data available.</p></section>"
        )
    offset = 25.0
    circles = []
    legend = []
    for idx, (label, value) in enumerate(items):
        pct = (value / total) * 100.0
        color = (color_map or {}).get(label, CHART_COLORS[idx % len(CHART_COLORS)])
        circles.append(
            "<circle class=\"donut-segment\" cx=\"90\" cy=\"90\" r=\"62\" "
            f"stroke=\"{color}\" stroke-dasharray=\"{pct:.4f} {100 - pct:.4f}\" "
            f"stroke-dashoffset=\"{offset:.4f}\" pathLength=\"100\"></circle>"
        )
        legend.append(
            "<li>"
            f"<span class=\"swatch\" style=\"background:{color}\"></span>"
            f"<span>{html_value(label)}</span>"
            f"<strong>{fmt_number(value)}</strong>"
            f"<em>{pct:.1f}%</em>"
            "</li>"
        )
        offset -= pct
    return (
        f"<section class=\"chart-panel\"><h2>{html_value(title)}</h2>"
        "<div class=\"donut-layout\">"
        "<svg class=\"donut\" viewBox=\"0 0 180 180\" role=\"img\" "
        f"aria-label=\"{html_value(title)} donut chart\">"
        "<circle class=\"donut-track\" cx=\"90\" cy=\"90\" r=\"62\"></circle>"
        + "".join(circles)
        + f"<text x=\"90\" y=\"84\" text-anchor=\"middle\" class=\"donut-total\">{fmt_number(total)}</text>"
        + "<text x=\"90\" y=\"105\" text-anchor=\"middle\" class=\"donut-caption\">total</text>"
        + "</svg>"
        + "<ul class=\"chart-legend\">"
        + "".join(legend)
        + "</ul></div></section>"
    )


def html_bar_chart(title: str, items: Sequence[Tuple[str, int]], max_items: int = 12) -> str:
    chart_items = list(items)[:max_items]
    max_value = max((value for _key, value in chart_items), default=0)
    if max_value <= 0:
        return (
            f"<section class=\"chart-panel wide\"><h2>{html_value(title)}</h2>"
            "<p class=\"empty-note\">No chart data available.</p></section>"
        )
    rows = []
    for idx, (label, value) in enumerate(chart_items):
        pct = (value / max_value) * 100.0
        color = CHART_COLORS[idx % len(CHART_COLORS)]
        rows.append(
            "<div class=\"bar-row\">"
            f"<div class=\"bar-label\" title=\"{html_value(label)}\">{html_value(label)}</div>"
            "<div class=\"bar-track\">"
            f"<div class=\"bar-fill\" style=\"width:{pct:.3f}%;background:{color}\"></div>"
            "</div>"
            f"<div class=\"bar-value\">{fmt_number(value)}</div>"
            "</div>"
        )
    return (
        f"<section class=\"chart-panel wide\"><h2>{html_value(title)}</h2>"
        "<div class=\"bar-chart\" role=\"img\" "
        f"aria-label=\"{html_value(title)} horizontal bar chart\">"
        + "".join(rows)
        + "</div></section>"
    )


def html_stacked_bar(title: str, items: Sequence[Tuple[str, int]], color_map: Optional[Dict[str, str]] = None) -> str:
    total = sum(value for _key, value in items)
    if total <= 0:
        return (
            f"<section class=\"chart-panel\"><h2>{html_value(title)}</h2>"
            "<p class=\"empty-note\">No chart data available.</p></section>"
        )
    segments = []
    legend = []
    for idx, (label, value) in enumerate(items):
        pct = (value / total) * 100.0
        color = (color_map or {}).get(label, CHART_COLORS[idx % len(CHART_COLORS)])
        segments.append(
            f"<div class=\"stack-segment\" style=\"width:{pct:.4f}%;background:{color}\" "
            f"title=\"{html_value(label)}: {fmt_number(value)} ({pct:.1f}%)\"></div>"
        )
        legend.append(
            "<li>"
            f"<span class=\"swatch\" style=\"background:{color}\"></span>"
            f"<span>{html_value(label)}</span>"
            f"<strong>{fmt_number(value)}</strong>"
            f"<em>{pct:.1f}%</em>"
            "</li>"
        )
    return (
        f"<section class=\"chart-panel\"><h2>{html_value(title)}</h2>"
        "<div class=\"stacked-bar\" role=\"img\" "
        f"aria-label=\"{html_value(title)} stacked bar chart\">"
        + "".join(segments)
        + "</div><ul class=\"chart-legend compact\">"
        + "".join(legend)
        + "</ul></section>"
    )


def consequence_group_counts(consequence_counts: object) -> List[Tuple[str, int]]:
    if not isinstance(consequence_counts, dict):
        return []
    coding = {
        "missense",
        "synonymous",
        "frameshift",
        "inframe_deletion",
        "inframe_insertion",
        "stop_gained",
        "stop_gained_early",
        "stop_gained_terminal",
        "stop_lost_readthrough",
        "stop_retained",
        "start_lost",
        "cds_complex",
        "exon_boundary_disruption",
    }
    splice = {"splice_acceptor_donor", "splice_region"}
    regulatory = {"promoter", "utr5", "utr3"}
    grouped = Counter()
    for key, value in consequence_counts.items():
        key_s = str(key)
        count = html_int(value)
        if key_s in coding:
            grouped["coding"] += count
        elif key_s in splice:
            grouped["splice"] += count
        elif key_s in regulatory:
            grouped["regulatory"] += count
        elif key_s == "intron":
            grouped["intron"] += count
        elif key_s == "intergenic":
            grouped["intergenic"] += count
        else:
            grouped["other"] += count
    order = ["coding", "splice", "regulatory", "intron", "intergenic", "other"]
    return [(key, grouped[key]) for key in order if grouped[key] > 0]


def html_filter_strategy_section(
    summary: Dict[str, object],
    level_counts: Sequence[Tuple[str, int]],
    consequence_counts: Sequence[Tuple[str, int]],
) -> str:
    levels = dict(level_counts)
    consequences = dict(consequence_counts)
    high = levels.get("HIGH", 0)
    moderate = levels.get("MODERATE", 0)
    low = levels.get("LOW", 0)
    coding_terms = [
        "missense",
        "frameshift",
        "stop_gained_early",
        "stop_gained_terminal",
        "start_lost",
        "splice_acceptor_donor",
        "exon_boundary_disruption",
    ]
    coding_signal = sum(consequences.get(term, 0) for term in coding_terms)
    regulatory_signal = consequences.get("promoter", 0) + consequences.get("utr5", 0) + consequences.get("utr3", 0)
    splice_signal = consequences.get("splice_acceptor_donor", 0) + consequences.get("splice_region", 0)
    sample_gating = summary.get("sample_count_gating") or "Sample-count gating was not recorded."
    strategies = [
        (
            "Strict high-confidence triage",
            f"{fmt_number(high)} HIGH alleles",
            "Start with DNP_LEVEL=HIGH, then prioritize frameshift, stop/start altering, splice donor/acceptor, and exon-boundary consequences. Use this for the smallest manual-review queue.",
        ),
        (
            "Broad candidate discovery",
            f"{fmt_number(high + moderate)} HIGH+MODERATE alleles",
            "Use DNP_LEVEL in HIGH or MODERATE when recall matters. Follow with consequence class, confidence, and sample-specific genotype evidence before biological interpretation.",
        ),
        (
            "Coding-impact focus",
            f"{fmt_number(coding_signal)} coding/splice-like alleles",
            "Filter for missense, frameshift, stop/start altering, splice_acceptor_donor, and exon_boundary_disruption. This is the most direct protein-effect review strategy.",
        ),
        (
            "Regulatory and splicing focus",
            f"{fmt_number(regulatory_signal + splice_signal)} regulatory/splice alleles",
            "Filter promoter, UTR, splice_region, and splice_acceptor_donor records. Inspect DNP_PROM, DNP_UTR, DNP_SPLICE, DNP_SPLICE_PWM, and DNP_SPLICE_AUX before ranking.",
        ),
        (
            "Low-priority background",
            f"{fmt_number(low)} LOW alleles",
            "LOW and MINIMAL variants are usually background for broad scans. Keep them for completeness, but avoid manual review unless a gene, interval, or sample hypothesis requires it.",
        ),
    ]
    cards = []
    for title, estimate, description in strategies:
        cards.append(
            "<article class=\"strategy-card\">"
            f"<h3>{html_value(title)}</h3>"
            f"<div class=\"strategy-count\">{html_value(estimate)}</div>"
            f"<p>{html_value(description)}</p>"
            "</article>"
        )
    return (
        "<section class=\"strategy-section\">"
        "<div class=\"section-heading\"><h2>Filtering Strategy Guide</h2>"
        "<p>These strategies turn the summary counts into practical review queues. DNP_SCORE is a ranking score, not a calibrated pathogenicity probability.</p></div>"
        "<div class=\"strategy-grid\">"
        + "".join(cards)
        + "</div>"
        "<div class=\"callout\"><strong>Sample-count gating:</strong> "
        + html_value(sample_gating)
        + "</div>"
        "</section>"
    )


def render_html_report(summary: Dict[str, object], path: str) -> None:
    active_methods = summary.get("active_methods", [])
    level_counts = summary.get("level_counts", {})
    consequence_counts = summary.get("consequence_counts", {})
    benchmark = summary.get("benchmark", {})
    index_info = summary.get("index", {})
    validation = summary.get("validation", {})
    reference_validation = validation.get("reference", {}) if isinstance(validation, dict) else {}
    method_source_types = summary.get("method_source_types", [])

    active_rows = [(idx + 1, method) for idx, method in enumerate(active_methods if isinstance(active_methods, list) else [])]
    level_items = normalized_count_items(level_counts, ("HIGH", "MODERATE", "LOW", "MINIMAL"))
    consequence_items = normalized_count_items(consequence_counts)
    level_rows = level_items
    consequence_rows = consequence_items[:25]
    benchmark_rows = list(benchmark.items()) if isinstance(benchmark, dict) else []
    index_rows = list(index_info.items()) if isinstance(index_info, dict) else []
    validation_rows = list(reference_validation.items()) if isinstance(reference_validation, dict) else []
    source_rows = []
    if isinstance(method_source_types, list):
        for item in method_source_types:
            if isinstance(item, dict):
                source_rows.append(
                    (
                        item.get("field", "."),
                        f"{item.get('source_type', '.')}: {item.get('description', '.')}",
                    )
                )

    overview_rows = [
        ("Input VCF", summary.get("input_vcf")),
        ("Output VCF", summary.get("output_vcf")),
        ("Config", summary.get("config")),
        ("Sample info", summary.get("sample_info")),
        ("Protein domains", summary.get("protein_domains")),
        ("Protein domain records", summary.get("protein_domain_records")),
        ("Protein structures", summary.get("protein_structures")),
        ("Protein structure records", summary.get("protein_structure_records")),
        ("Protein LM scores", summary.get("protein_lm_scores")),
        ("Protein LM score records", summary.get("protein_lm_score_records")),
        ("miRNA sites", summary.get("mirna_sites")),
        ("miRNA site records", summary.get("mirna_site_records")),
        ("ML model", summary.get("ml_model")),
        ("ML model loaded", summary.get("ml_model_loaded")),
        ("ML model type", summary.get("ml_model_type")),
        ("ML model top features", summary.get("ml_model_top_features")),
        ("Samples", summary.get("n_samples")),
        ("Sample-count gating", summary.get("sample_count_gating")),
        ("Regions", summary.get("regions")),
        ("Input records seen", summary.get("input_records_seen")),
        ("Records scored", summary.get("records_scored")),
        ("Records skipped by region", summary.get("records_skipped_by_region")),
        ("ALT alleles scored", summary.get("alt_alleles_scored")),
        ("REF mismatch records", reference_validation.get("mismatch_records") if isinstance(reference_validation, dict) else None),
        ("Population window size", summary.get("population_window_size")),
        ("Gene constraint genes", summary.get("gene_constraint_genes")),
        ("Sample info groups", summary.get("sample_info_groups")),
        ("Case-control enabled", summary.get("sample_info_has_case_control")),
    ]

    input_vcf = summary.get("input_vcf") or "VCF"
    records_scored = html_int(summary.get("records_scored"))
    alt_alleles = html_int(summary.get("alt_alleles_scored"))
    n_samples = html_int(summary.get("n_samples"))
    mismatch_records = html_int(reference_validation.get("mismatch_records") if isinstance(reference_validation, dict) else 0)
    unchecked_records = html_int(reference_validation.get("unchecked_records") if isinstance(reference_validation, dict) else 0)
    elapsed_seconds = html_float(benchmark.get("elapsed_seconds") if isinstance(benchmark, dict) else 0.0)
    records_per_second = html_float(benchmark.get("records_per_second") if isinstance(benchmark, dict) else 0.0)
    peak_rss = html_float(benchmark.get("peak_rss_mb") if isinstance(benchmark, dict) else 0.0)
    validation_chart_items = [
        ("matching", html_int(reference_validation.get("matching_records") if isinstance(reference_validation, dict) else 0)),
        ("mismatch", mismatch_records),
        ("unchecked", unchecked_records),
    ]
    validation_colors = {"matching": "#0f766e", "mismatch": "#be123c", "unchecked": "#b45309"}
    benchmark_chart_items = []
    if isinstance(benchmark, dict):
        benchmark_chart_items = [
            (
                "scoring_seconds",
                max(0, int(round(html_float(benchmark.get("scoring_seconds"))))),
            ),
            (
                "prescan_seconds",
                max(0, int(round(html_float(benchmark.get("prescan_seconds"))))),
            ),
            (
                "index_seconds",
                max(0, int(round(html_float(benchmark.get("index_seconds"))))),
            ),
        ]
    consequence_group_items = consequence_group_counts(consequence_counts)

    metric_cards = [
        html_metric_card("Records scored", fmt_number(records_scored), "input records processed"),
        html_metric_card("ALT alleles", fmt_number(alt_alleles), "allele-level scores"),
        html_metric_card("Samples", fmt_number(n_samples), summary.get("sample_count_gating", "")),
        html_metric_card("REF mismatches", fmt_number(mismatch_records), f"{fmt_number(unchecked_records)} unchecked"),
        html_metric_card("Elapsed time", fmt_duration(elapsed_seconds), f"{fmt_number(records_per_second)} records/s"),
        html_metric_card("Peak RSS", f"{fmt_number(peak_rss)} MB", "maximum resident memory"),
    ]

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DeNovoPath Run Report</title>
<style>
:root {{
  color-scheme: light;
  --fg: #17202a;
  --muted: #5b6876;
  --soft: #eef3f7;
  --line: #d9e2ea;
  --bg: #f5f7fa;
  --panel: #ffffff;
  --ink: #0f172a;
  --accent: #0f766e;
  --accent-2: #2563eb;
  --warn: #b45309;
  --danger: #be123c;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: Arial, Helvetica, sans-serif;
  color: var(--fg);
  background: var(--bg);
  line-height: 1.5;
}}
main {{ max-width: 1180px; margin: 0 auto; padding: 28px 18px 48px; }}
header {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 6px;
  padding: 22px;
  margin-bottom: 18px;
}}
.eyebrow {{ margin: 0 0 6px; color: var(--accent); font-weight: 700; font-size: 13px; text-transform: uppercase; }}
h1 {{ margin: 0; font-size: 30px; line-height: 1.2; color: var(--ink); }}
h2 {{ margin: 0 0 12px; font-size: 18px; line-height: 1.3; color: var(--ink); }}
h3 {{ margin: 0 0 8px; font-size: 15px; line-height: 1.35; color: var(--ink); }}
p {{ margin: 0; color: var(--muted); }}
.subtitle {{ margin-top: 8px; max-width: 840px; }}
.run-meta {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 16px; }}
.chip {{
  display: inline-flex;
  align-items: center;
  min-height: 28px;
  padding: 4px 10px;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: #f8fafc;
  color: #334155;
  font-size: 13px;
}}
section {{
  margin-top: 16px;
  padding: 16px;
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 6px;
}}
.section-heading {{ margin-bottom: 14px; }}
.metric-grid {{
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 12px;
  margin: 16px 0;
}}
.metric-card {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 6px;
  padding: 14px;
  min-height: 106px;
}}
.metric-label {{ color: var(--muted); font-size: 13px; font-weight: 700; }}
.metric-value {{ margin-top: 6px; color: var(--ink); font-size: 26px; font-weight: 700; line-height: 1.15; }}
.metric-detail {{ margin-top: 6px; color: var(--muted); font-size: 13px; overflow-wrap: anywhere; }}
.chart-grid {{
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 16px;
  align-items: stretch;
}}
.chart-panel {{ margin-top: 0; min-height: 260px; }}
.chart-panel.wide {{ grid-column: 1 / -1; }}
.donut-layout {{ display: grid; grid-template-columns: 190px 1fr; gap: 16px; align-items: center; }}
.donut {{ width: 180px; height: 180px; }}
.donut-track {{ fill: none; stroke: #e5edf3; stroke-width: 22; }}
.donut-segment {{ fill: none; stroke-width: 22; transform: rotate(-90deg); transform-origin: 90px 90px; }}
.donut-total {{ font-size: 18px; font-weight: 700; fill: var(--ink); }}
.donut-caption {{ font-size: 11px; fill: var(--muted); }}
.chart-legend {{ list-style: none; margin: 0; padding: 0; display: grid; gap: 8px; }}
.chart-legend li {{ display: grid; grid-template-columns: 12px 1fr auto auto; gap: 8px; align-items: center; color: #334155; font-size: 13px; }}
.chart-legend.compact {{ margin-top: 12px; }}
.swatch {{ width: 12px; height: 12px; border-radius: 3px; display: inline-block; }}
.chart-legend strong {{ color: var(--ink); font-weight: 700; font-variant-numeric: tabular-nums; }}
.chart-legend em {{ color: var(--muted); font-style: normal; font-variant-numeric: tabular-nums; }}
.bar-chart {{ display: grid; gap: 10px; }}
.bar-row {{ display: grid; grid-template-columns: minmax(140px, 220px) 1fr minmax(76px, auto); gap: 10px; align-items: center; }}
.bar-label {{ color: #334155; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
.bar-track {{ height: 12px; border-radius: 999px; background: #e8eef4; overflow: hidden; }}
.bar-fill {{ height: 100%; border-radius: 999px; min-width: 2px; }}
.bar-value {{ text-align: right; color: var(--ink); font-weight: 700; font-variant-numeric: tabular-nums; }}
.stacked-bar {{ display: flex; height: 24px; overflow: hidden; border-radius: 6px; background: #e8eef4; }}
.stack-segment {{ min-width: 0; }}
.strategy-grid {{ display: grid; grid-template-columns: repeat(5, minmax(0, 1fr)); gap: 12px; }}
.strategy-card {{
  border: 1px solid var(--line);
  border-radius: 6px;
  padding: 14px;
  background: #fbfdff;
  min-height: 190px;
}}
.strategy-count {{ margin: 8px 0; color: var(--accent-2); font-weight: 700; font-size: 18px; }}
.strategy-card p {{ font-size: 13px; }}
.callout {{
  margin-top: 12px;
  padding: 12px 14px;
  border-left: 4px solid var(--accent);
  background: #ecfdf5;
  color: #164e40;
  border-radius: 6px;
}}
table {{ width: 100%; border-collapse: collapse; table-layout: fixed; }}
th, td {{ padding: 8px 10px; border-top: 1px solid var(--line); vertical-align: top; overflow-wrap: anywhere; text-align: left; }}
tbody tr:first-child th, tbody tr:first-child td {{ border-top: 0; }}
th {{ width: 30%; color: var(--muted); font-weight: 700; }}
.method-list {{ columns: 2 320px; margin: 0; padding-left: 22px; }}
.method-list li {{ margin: 0 0 6px; break-inside: avoid; }}
.empty-note {{ color: var(--muted); }}
@media (max-width: 900px) {{
  .metric-grid, .chart-grid, .strategy-grid {{ grid-template-columns: 1fr; }}
  .chart-panel.wide {{ grid-column: auto; }}
  .donut-layout {{ grid-template-columns: 1fr; }}
  .bar-row {{ grid-template-columns: 1fr; gap: 6px; }}
  .bar-value {{ text-align: left; }}
  th {{ width: 42%; }}
}}
</style>
</head>
<body>
<main>
<header>
<p class="eyebrow">DeNovoPath Report</p>
<h1>DeNovoPath Run Report</h1>
<p class="subtitle">Sample-count-aware VCF deleteriousness scoring summary for {html_value(input_vcf)}.</p>
<div class="run-meta">
<span class="chip">Samples: {fmt_number(n_samples)}</span>
<span class="chip">Records: {fmt_number(records_scored)}</span>
<span class="chip">REF mismatches: {fmt_number(mismatch_records)}</span>
<span class="chip">Elapsed: {fmt_duration(elapsed_seconds)}</span>
</div>
</header>
"""
    document += "<section><div class=\"section-heading\"><h2>Run At A Glance</h2><p>Core output, validation, and performance indicators.</p></div><div class=\"metric-grid\">"
    document += "".join(metric_cards)
    document += "</div></section>"
    document += "<div class=\"chart-grid\">"
    document += html_donut_chart("Impact Level Distribution", level_items, LEVEL_COLORS)
    document += html_stacked_bar("Reference Validation", validation_chart_items, validation_colors)
    document += html_bar_chart("Top Consequence Classes", consequence_items, max_items=14)
    document += html_stacked_bar("Annotation Category Mix", consequence_group_items)
    document += html_stacked_bar("benchmark Runtime Split", benchmark_chart_items)
    document += "</div>"
    document += html_filter_strategy_section(summary, level_items, consequence_items)
    document += html_table("Overview", overview_rows)
    document += html_table("Reference Validation", validation_rows)
    document += html_table("Method Source Types", source_rows)
    document += "<section><h2>Active Methods</h2>"
    if active_rows:
        document += "<ol class=\"method-list\">" + "".join(f"<li>{html_value(method)}</li>" for _, method in active_rows) + "</ol>"
    else:
        document += "<p>No active methods recorded.</p>"
    document += "</section>"
    document += html_table("Level Counts", level_rows)
    document += html_table("Consequence Counts", consequence_rows)
    document += html_table("benchmark", benchmark_rows)
    document += html_table("Index", index_rows)
    document += """
</main>
</body>
</html>
"""
    with open(path, "w") as out:
        out.write(document)


def score_vcf(args: argparse.Namespace) -> None:
    run_start = time.perf_counter()
    config = load_score_config(args.config)
    regions = [parse_region(item) for item in (args.region or [])]
    sample_info = load_sample_info(getattr(args, "sample_info", None))
    pop_window_size = max(0, int(getattr(args, "pop_window", 100000)))
    population_windows_enabled = pop_window_size > 0
    gene_constraint_enabled = not bool(getattr(args, "skip_gene_constraint", False))
    prescan_start = time.perf_counter()
    population_windows, precomputed_n_samples = precompute_population_windows(
        args.vcf,
        regions,
        pop_window_size,
        limit=args.limit,
    )
    scorer = DeNovoPathScorer(
        reference_fasta=args.reference,
        gff=args.gff,
        cds_fasta=args.cds,
        protein_fasta=args.pep,
        protein_domains=getattr(args, "protein_domains", None),
        protein_structures=getattr(args, "protein_structures", None),
        protein_lm_scores=getattr(args, "protein_lm_scores", None),
        mirna_sites=getattr(args, "mirna_sites", None),
        ml_model=getattr(args, "ml_model", None),
        window=args.window,
        config=config,
    )
    if gene_constraint_enabled:
        gene_constraint_stats, constraint_n_samples = precompute_gene_constraint(
            args.vcf,
            regions,
            scorer,
            limit=args.limit,
        )
    else:
        gene_constraint_stats, constraint_n_samples = {}, 0
    prescan_seconds = time.perf_counter() - prescan_start
    n_samples = 0
    sample_names: List[str] = []
    records = 0
    input_records = 0
    skipped_region = 0
    alt_records = 0
    ref_match_records = 0
    ref_mismatch_records = 0
    ref_unchecked_records = 0
    ref_mismatch_examples: List[str] = []
    level_counts: Counter[str] = Counter()
    consequence_counts: Counter[str] = Counter()
    wrote_info_headers = False
    scoring_start = time.perf_counter()
    with open_text(args.vcf, "r") as inp, open_text(args.output, "w") as out:
        for line in inp:
            if line.startswith("##INFO=<ID=DNP_"):
                continue
            if line.startswith("##DeNovoPathActiveMethods="):
                continue
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                sample_names = header_fields[9:]
                n_samples = max(0, len(header_fields) - 9)
                if precomputed_n_samples and precomputed_n_samples != n_samples:
                    n_samples = precomputed_n_samples
                if constraint_n_samples and constraint_n_samples != n_samples:
                    n_samples = constraint_n_samples
                if not wrote_info_headers:
                    for header in INFO_HEADERS:
                        out.write(header + "\n")
                    out.write(
                        "##DeNovoPathActiveMethods="
                        + ",".join(
                            active_methods(
                                n_samples,
                                sample_info is not None,
                                bool(sample_info and sample_info.has_phenotype()),
                                has_protein_domains=bool(len(scorer.protein_domains)),
                                has_protein_structures=bool(len(scorer.protein_structures)),
                                has_protein_esm_scores=bool(len(scorer.protein_lm_scores)),
                                has_mirna_sites=bool(len(scorer.mirna_sites)),
                                has_ml_model=bool(len(scorer.ml_model)),
                                has_population_windows=population_windows_enabled,
                                has_gene_constraint=gene_constraint_enabled,
                            )
                        )
                        + "\n"
                    )
                    wrote_info_headers = True
                out.write(line)
                continue
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                out.write(line)
                continue
            input_records += 1
            chrom = fields[0]
            pos = int(fields[1])
            if not in_selected_regions(chrom, pos, regions):
                skipped_region += 1
                continue
            ref_status, observed_ref = scorer.validate_reference_allele(chrom, pos, fields[3])
            if ref_status == "match":
                ref_match_records += 1
            elif ref_status == "mismatch":
                ref_mismatch_records += 1
                example = f"{chrom}:{pos} VCF_REF={fields[3]} FASTA_REF={observed_ref}"
                if len(ref_mismatch_examples) < 10:
                    ref_mismatch_examples.append(example)
                if getattr(args, "strict_ref", False):
                    raise ValueError(f"VCF REF does not match reference FASTA: {example}")
            else:
                ref_unchecked_records += 1
            scores = scorer.score_record(fields, n_samples, sample_names=sample_names, sample_info=sample_info)
            window_stat = (
                population_windows.get((chrom, window_index(pos, pop_window_size)), WindowPopulationStat())
                if pop_window_size > 0
                else WindowPopulationStat()
            )
            for score in scores:
                score.window_pi = window_stat.pi
                score.window_theta = window_stat.theta
                score.window_tajima_d = window_stat.tajima_d
                score.window_ld = window_stat.ld
                score.window_haplotype = window_stat.haplotype
                gene_stat = gene_constraint_stats.get(score.gene_id, GeneConstraintStat())
                score.gene_lof_oe = gene_stat.lof_oe
                score.gene_missense_oe = gene_stat.missense_oe
                score.gene_constraint_score = gene_stat.constraint_score
                scorer.apply_ml_scores(score)
            alt_records += len(scores)
            for score in scores:
                level_counts[score.level] += 1
                consequence_counts[score.consequence] += 1
            fields[7] = append_scores_to_info(fields[7], scores)
            out.write("\t".join(fields) + "\n")
            records += 1
            if args.limit and records >= args.limit:
                break
            if args.progress and records % args.progress == 0:
                print(f"processed {records} records", file=sys.stderr)
    scoring_seconds = time.perf_counter() - scoring_start
    index_start = time.perf_counter()
    index_result = maybe_index_vcf(args.output, getattr(args, "index_output", "auto"))
    index_seconds = time.perf_counter() - index_start
    elapsed_seconds = time.perf_counter() - run_start
    summary_path = getattr(args, "summary", None)
    html_report_path = getattr(args, "html_report", None)
    if summary_path or html_report_path:
        summary = build_run_summary(
            args,
            n_samples,
            sample_info,
            pop_window_size,
            gene_constraint_stats,
            input_records,
            records,
            skipped_region,
            alt_records,
            index_result,
            prescan_seconds,
            scoring_seconds,
            index_seconds,
            elapsed_seconds,
            level_counts,
            consequence_counts,
            protein_domain_records=len(scorer.protein_domains),
            protein_structure_records=len(scorer.protein_structures),
            protein_lm_score_records=len(scorer.protein_lm_scores),
            mirna_site_records=len(scorer.mirna_sites),
            ml_model_loaded=bool(len(scorer.ml_model)),
            ml_model_type=scorer.ml_model.model_type if len(scorer.ml_model) else "",
            ml_model_features=scorer.ml_model.features,
            ml_feature_importance=scorer.ml_model.feature_importance,
            ref_match_records=ref_match_records,
            ref_mismatch_records=ref_mismatch_records,
            ref_unchecked_records=ref_unchecked_records,
            ref_mismatch_examples=ref_mismatch_examples,
            population_windows_enabled=population_windows_enabled,
            gene_constraint_enabled=gene_constraint_enabled,
        )
    else:
        summary = None
    if summary_path and summary is not None:
        with open(summary_path, "w") as out:
            json.dump(summary, out, indent=2, sort_keys=True)
            out.write("\n")
    if html_report_path and summary is not None:
        render_html_report(summary, html_report_path)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Annotate VCF INFO with DeNovoPath scores.")
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("--gff", required=True, help="GFF3/GTF-like gene annotation")
    parser.add_argument("--cds", help="Transcript CDS FASTA; headers should match mRNA IDs")
    parser.add_argument("--pep", help="Protein FASTA; currently loaded for ID validation/future extensions")
    parser.add_argument(
        "--protein-domains",
        help="Optional TSV/CSV protein domain annotation with transcript,start,end,domain[,score] columns",
    )
    parser.add_argument(
        "--protein-structures",
        help="Optional TSV/CSV AlphaFold/ESMFold-like residue structure annotation with transcript,aa_pos/start/end,plddt[,rsa,ss] columns",
    )
    parser.add_argument(
        "--protein-lm-scores",
        help="Optional TSV/CSV precomputed ESM-2-like protein substitution scores with transcript,aa_pos,ref,alt,score columns",
    )
    parser.add_argument(
        "--mirna-sites",
        help="Optional TSV/CSV miRNA target/seed sites with chrom,start,end,mirna[,transcript,site,score] columns",
    )
    parser.add_argument(
        "--ml-model",
        help="Optional DeNovoPath JSON ML model for DNP_ML/DNP_CAL/DNP_UNCERT/DNP_OOD/DNP_FEATIMP",
    )
    parser.add_argument(
        "--sample-info",
        help="Optional TSV/CSV with sample and group columns, plus optional phenotype/status/case_control",
    )
    parser.add_argument("--output", required=True, help="Output VCF/VCF.GZ")
    parser.add_argument("--window", type=int, default=50, help="Sequence-context half-window size")
    parser.add_argument(
        "--pop-window",
        type=int,
        default=100000,
        help="Fixed window size for pi/theta/Tajima D cohort statistics; 0 disables",
    )
    parser.add_argument(
        "--skip-gene-constraint",
        action="store_true",
        help="Skip cohort gene-constraint prescan and leave DNP_GENELOF/DNP_GENEMIS/DNP_GENECON at 0",
    )
    parser.add_argument("--config", help="Optional YAML config for score weights and level thresholds")
    parser.add_argument(
        "--region",
        action="append",
        help="Optional region filter, e.g. chr1 or chr1:1000-2000; may be supplied multiple times",
    )
    parser.add_argument("--summary", help="Optional JSON summary output path")
    parser.add_argument("--html-report", help="Optional self-contained HTML run report path")
    parser.add_argument(
        "--index-output",
        choices=["auto", "never", "always"],
        default="auto",
        help="Optionally bgzip-recompress and tabix-index .vcf.gz output when tools are available",
    )
    parser.add_argument("--limit", type=int, default=0, help="Stop after N records, for smoke tests")
    parser.add_argument("--progress", type=int, default=100000, help="Progress interval; 0 disables")
    parser.add_argument(
        "--strict-ref",
        action="store_true",
        help="Fail when VCF REF does not match the reference FASTA; by default mismatches are counted in the summary",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    score_vcf(args)
    return 0
