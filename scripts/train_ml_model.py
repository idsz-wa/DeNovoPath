#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from denovopath.scorer import ML_FEATURES  # noqa: E402


FEATURE_INFO_KEYS = {
    "impact_score": "DNP_IMPACT",
    "protein_score": "DNP_PROT",
    "grantham_score": "DNP_GRANTHAM",
    "blosum_score": "DNP_BLOSUM",
    "codon_usage_score": "DNP_CODONUSE",
    "protein_context_score": "DNP_PROTCTX",
    "protein_structure_score": "DNP_STRUCT",
    "protein_structure_model_score": "DNP_AFSTRUCT",
    "protein_esm_score": "DNP_ESM",
    "protein_lm_score": "DNP_PROTLM",
    "protein_domain_score": "DNP_DOMAIN",
    "splice_score": "DNP_SPLICE",
    "splice_motif_score": "DNP_SPLICE_MOTIF",
    "splice_pwm_score": "DNP_SPLICE_PWM",
    "splice_maxent_score": "DNP_SPLICE_MAXENT",
    "splice_aux_score": "DNP_SPLICE_AUX",
    "splice_ese_score": "DNP_SPLICE_ESE",
    "utr_score": "DNP_UTR",
    "rnafold_score": "DNP_RNAFOLD",
    "mirna_score": "DNP_MIRNA",
    "promoter_score": "DNP_PROM",
    "sequence_score": "DNP_SEQ",
    "kmer_score": "DNP_KMER",
    "repeat_score": "DNP_REPEAT",
    "mutation_context_score": "DNP_MUTCTX",
    "dna_lm_score": "DNP_DNALM",
    "cohort_score": "DNP_COHORT",
    "hwe_score": "DNP_HWE",
    "heterozygosity_deviation_score": "DNP_HETDEV",
    "fst_score": "DNP_FST",
    "case_control_score": "DNP_CASECTRL",
    "window_pi": "DNP_PI",
    "window_theta": "DNP_THETA",
    "window_tajima_d": "DNP_TAJD",
    "window_ld": "DNP_LD",
    "window_haplotype": "DNP_HAP",
    "gene_lof_oe": "DNP_GENELOF",
    "gene_missense_oe": "DNP_GENEMIS",
    "gene_constraint_score": "DNP_GENECON",
    "qc_score": "DNP_QC",
    "confidence_score": "DNP_CONF",
}

HIGH_CONSEQUENCES = {
    "frameshift",
    "stop_gained",
    "stop_gained_early",
    "stop_lost",
    "stop_lost_readthrough",
    "start_lost",
    "splice_acceptor_donor",
    "exon_boundary_disruption",
}


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


def split_info_values(info: Dict[str, str], key: str, n: int, default: str = "0") -> List[str]:
    value = info.get(key, default)
    values = value.split(",") if value != "." else [default]
    if len(values) < n:
        values.extend([default] * (n - len(values)))
    return values[:n]


def to_float(value: str, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(result) or math.isinf(result):
        return default
    return result


def read_labels(path: Optional[str]) -> Dict[Tuple[str, str, str, str], int]:
    if not path:
        return {}
    labels: Dict[Tuple[str, str, str, str], int] = {}
    delimiter = "," if path.endswith(".csv") else "\t"
    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            lowered = {str(key).strip().lower(): str(value).strip() for key, value in row.items() if key}
            chrom = lowered.get("chrom") or lowered.get("#chrom")
            pos = lowered.get("pos") or lowered.get("position")
            ref = lowered.get("ref")
            alt = lowered.get("alt")
            label = lowered.get("label") or lowered.get("y") or lowered.get("target")
            if not chrom or not pos or not ref or not alt or label is None:
                continue
            label_value = str(label).lower()
            labels[(chrom, pos, ref, alt)] = 1 if label_value in {"1", "true", "yes", "pathogenic", "deleterious"} else 0
    return labels


def vcf_samples(path: str) -> int:
    with open_text(path, "r") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return max(0, len(line.rstrip("\n").split("\t")) - 9)
    return 0


def pseudo_label(
    score: float,
    qc: float,
    consequence: str,
    high_threshold: float,
    low_threshold: float,
    min_qc: float,
) -> Optional[int]:
    if qc < min_qc:
        return None
    if score >= high_threshold or consequence in HIGH_CONSEQUENCES:
        return 1
    if score <= low_threshold and consequence in {"intergenic", "intron", "synonymous", "utr3", "utr5", "promoter"}:
        return 0
    return None


def collect_examples(args: argparse.Namespace) -> Tuple[List[List[float]], List[int], Dict[str, object]]:
    labels = read_labels(args.labels)
    n_samples = vcf_samples(args.vcf)
    allow_pseudo = not labels and n_samples >= args.min_samples_for_pseudo
    if not labels and not allow_pseudo:
        raise SystemExit(
            "Pseudo-label training is disabled for this VCF because sample count "
            f"({n_samples}) is below --min-samples-for-pseudo ({args.min_samples_for_pseudo}). "
            "Provide --labels or lower the threshold intentionally."
        )

    x_rows: List[List[float]] = []
    y_rows: List[int] = []
    input_records = 0
    alt_records = 0
    labeled_records = 0
    skipped_unlabeled = 0
    explicit_labeled_records = 0
    pseudo_labeled_records = 0

    with open_text(args.vcf, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            input_records += 1
            chrom, pos, _var_id, ref, alt_text = fields[:5]
            alts = alt_text.split(",")
            info = parse_info(fields[7])
            feature_values = {
                feature: split_info_values(info, info_key, len(alts), "0")
                for feature, info_key in FEATURE_INFO_KEYS.items()
            }
            scores = split_info_values(info, "DNP_SCORE", len(alts), "0")
            qcs = split_info_values(info, "DNP_QC", len(alts), "0")
            consequences = split_info_values(info, "DNP_CONSEQ", len(alts), ".")
            for idx, alt in enumerate(alts):
                alt_records += 1
                row = [to_float(feature_values[feature][idx]) for feature in args.features]
                explicit_label = labels.get((chrom, pos, ref, alt))
                label_source = "explicit"
                if explicit_label is None:
                    if not allow_pseudo:
                        skipped_unlabeled += 1
                        continue
                    label_source = "pseudo"
                    explicit_label = pseudo_label(
                        to_float(scores[idx]),
                        to_float(qcs[idx]),
                        consequences[idx],
                        args.high_threshold,
                        args.low_threshold,
                        args.min_qc,
                    )
                if explicit_label is None:
                    skipped_unlabeled += 1
                    continue
                x_rows.append(row)
                y_rows.append(int(explicit_label))
                labeled_records += 1
                if label_source == "explicit":
                    explicit_labeled_records += 1
                else:
                    pseudo_labeled_records += 1
                if args.max_examples and labeled_records >= args.max_examples:
                    break
            if args.max_examples and labeled_records >= args.max_examples:
                break

    summary = {
        "input_vcf": args.vcf,
        "label_source": "labels" if labels else "pseudo_labels",
        "n_samples": n_samples,
        "input_records": input_records,
        "alt_records": alt_records,
        "labeled_alt_records": labeled_records,
        "label_source_counts": {
            "explicit": explicit_labeled_records,
            "pseudo": pseudo_labeled_records,
        },
        "skipped_unlabeled_alt_records": skipped_unlabeled,
        "positive_labels": int(sum(y_rows)),
        "negative_labels": int(len(y_rows) - sum(y_rows)),
    }
    return x_rows, y_rows, summary


def standardize_matrix(x_rows: List[List[float]]) -> Tuple[List[List[float]], List[float], List[float]]:
    if not x_rows:
        return [], [], []
    n_features = len(x_rows[0])
    means = []
    scales = []
    for col in range(n_features):
        values = [row[col] for row in x_rows]
        mean = sum(values) / len(values)
        variance = sum((value - mean) ** 2 for value in values) / max(1, len(values) - 1)
        scale = max(1e-6, math.sqrt(variance))
        means.append(mean)
        scales.append(scale)
    standardized = [[(row[col] - means[col]) / scales[col] for col in range(n_features)] for row in x_rows]
    return standardized, means, scales


def fit_calibration(logits: Sequence[float], labels: Sequence[int]) -> Dict[str, float]:
    if len(set(labels)) < 2 or len(labels) < 20:
        return {"scale": 1.0, "shift": 0.0}
    try:
        from sklearn.linear_model import LogisticRegression
    except ImportError:
        return {"scale": 1.0, "shift": 0.0}
    model = LogisticRegression(solver="lbfgs")
    model.fit([[value] for value in logits], labels)
    return {"scale": float(model.coef_[0][0]), "shift": float(model.intercept_[0])}


def balanced_sample_weight(labels: Sequence[int]) -> List[float]:
    n = len(labels)
    positives = sum(1 for item in labels if item == 1)
    negatives = n - positives
    if positives == 0 or negatives == 0:
        return [1.0] * n
    return [n / (2.0 * positives) if item == 1 else n / (2.0 * negatives) for item in labels]


def safe_logit_local(probability: float) -> float:
    probability = max(1e-6, min(1.0 - 1e-6, probability))
    return math.log(probability / (1.0 - probability))


def export_sklearn_tree(tree) -> Dict[str, List[float]]:
    raw_value = tree.tree_.value
    return {
        "children_left": [int(item) for item in tree.tree_.children_left.tolist()],
        "children_right": [int(item) for item in tree.tree_.children_right.tolist()],
        "feature": [int(item) for item in tree.tree_.feature.tolist()],
        "threshold": [float(item) for item in tree.tree_.threshold.tolist()],
        "value": [float(item[0][0]) for item in raw_value.tolist()],
    }


def initial_log_odds(model, n_features: int) -> float:
    try:
        probabilities = model.init_.predict_proba([[0.0] * n_features])[0]
        classes = list(model.classes_)
        idx = classes.index(1)
        return safe_logit_local(float(probabilities[idx]))
    except Exception:
        return 0.0


def train_model(args: argparse.Namespace) -> Dict[str, object]:
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.model_selection import train_test_split
        from sklearn.metrics import average_precision_score, brier_score_loss, matthews_corrcoef, roc_auc_score
    except ImportError as exc:
        raise SystemExit("scikit-learn is required in the conda py3 environment for model training") from exc

    x_rows, y_rows, summary = collect_examples(args)
    positives = int(sum(y_rows))
    negatives = int(len(y_rows) - positives)
    if positives < args.min_positive or negatives < args.min_negative:
        raise SystemExit(
            "Not enough labeled examples to train: "
            f"positive={positives}, negative={negatives}, "
            f"required positive>={args.min_positive}, negative>={args.min_negative}"
        )

    x_scaled, means, scales = standardize_matrix(x_rows)
    calibration = {"scale": 1.0, "shift": 0.0}
    metrics: Dict[str, object] = {}
    train_x = x_scaled
    train_y = y_rows
    if 0.0 < args.validation_fraction < 0.5 and len(y_rows) >= 40 and positives >= 5 and negatives >= 5:
        train_x, valid_x, train_y, valid_y = train_test_split(
            x_scaled,
            y_rows,
            test_size=args.validation_fraction,
            stratify=y_rows,
            random_state=args.random_state,
        )
        if args.backend == "sklearn_gbm":
            model_for_valid = GradientBoostingClassifier(
                n_estimators=args.n_estimators,
                learning_rate=args.learning_rate,
                max_depth=args.max_depth,
                random_state=args.random_state,
            )
            model_for_valid.fit(train_x, train_y, sample_weight=balanced_sample_weight(train_y))
        else:
            model_for_valid = LogisticRegression(max_iter=args.max_iter, class_weight=args.class_weight, random_state=args.random_state)
            model_for_valid.fit(train_x, train_y)
        valid_logits = model_for_valid.decision_function(valid_x)
        calibration = fit_calibration(list(valid_logits), valid_y)
        valid_probs = model_for_valid.predict_proba(valid_x)[:, 1]
        metrics["validation_records"] = len(valid_y)
        if len(set(valid_y)) >= 2:
            metrics["validation_auc"] = float(roc_auc_score(valid_y, valid_probs))
            metrics["validation_pr_auc"] = float(average_precision_score(valid_y, valid_probs))
            metrics["validation_mcc"] = float(matthews_corrcoef(valid_y, [1 if value >= 0.5 else 0 for value in valid_probs]))
            metrics["validation_brier"] = float(brier_score_loss(valid_y, valid_probs))

    metrics["training_records"] = len(y_rows)
    metrics["training_positive"] = positives
    metrics["training_negative"] = negatives
    metrics["backend"] = args.backend

    if args.backend == "sklearn_gbm":
        model = GradientBoostingClassifier(
            n_estimators=args.n_estimators,
            learning_rate=args.learning_rate,
            max_depth=args.max_depth,
            random_state=args.random_state,
        )
        model.fit(x_scaled, y_rows, sample_weight=balanced_sample_weight(y_rows))
        feature_importance = {
            feature: float(model.feature_importances_[idx])
            for idx, feature in enumerate(args.features)
        }
        return {
            "model_type": "sklearn_gradient_boosting_json",
            "features": args.features,
            "trees": [export_sklearn_tree(estimator[0]) for estimator in model.estimators_],
            "intercept": initial_log_odds(model, len(args.features)),
            "learning_rate": float(model.learning_rate),
            "mean": means,
            "scale": scales,
            "calibration": calibration,
            "feature_importance": feature_importance,
            "training_summary": summary,
            "metrics": metrics,
        }

    model = LogisticRegression(max_iter=args.max_iter, class_weight=args.class_weight, random_state=args.random_state)
    model.fit(x_scaled, y_rows)
    coef = [float(value) for value in model.coef_[0]]
    intercept = float(model.intercept_[0])
    feature_importance = {feature: abs(coef[idx]) for idx, feature in enumerate(args.features)}

    return {
        "model_type": "sklearn_logistic_regression_json",
        "features": args.features,
        "coef": coef,
        "intercept": intercept,
        "mean": means,
        "scale": scales,
        "calibration": calibration,
        "feature_importance": feature_importance,
        "training_summary": summary,
        "metrics": metrics,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train a DeNovoPath JSON ML model from a scored VCF.")
    parser.add_argument("--vcf", required=True, help="Scored VCF/VCF.GZ containing DNP_* INFO fields")
    parser.add_argument("--output", required=True, help="Output JSON model path")
    parser.add_argument("--labels", help="Optional TSV/CSV labels with chrom,pos,ref,alt,label columns")
    parser.add_argument(
        "--backend",
        choices=["logistic", "sklearn_gbm"],
        default="logistic",
        help="Training backend; sklearn_gbm exports a portable JSON gradient-boosted tree ensemble",
    )
    parser.add_argument(
        "--features",
        nargs="+",
        default=list(ML_FEATURES),
        choices=list(ML_FEATURES),
        help="Feature names to train on; defaults to all DeNovoPath deterministic ML features",
    )
    parser.add_argument("--high-threshold", type=float, default=0.80, help="Pseudo-positive DNP_SCORE threshold")
    parser.add_argument("--low-threshold", type=float, default=0.20, help="Pseudo-negative DNP_SCORE threshold")
    parser.add_argument("--min-qc", type=float, default=0.20, help="Minimum DNP_QC for pseudo-label examples")
    parser.add_argument(
        "--min-samples-for-pseudo",
        type=int,
        default=10,
        help="Minimum VCF sample count required before pseudo-label training is allowed",
    )
    parser.add_argument("--min-positive", type=int, default=25, help="Minimum positive examples required")
    parser.add_argument("--min-negative", type=int, default=25, help="Minimum negative examples required")
    parser.add_argument("--max-examples", type=int, default=0, help="Maximum labeled ALT examples to use; 0 means all")
    parser.add_argument("--validation-fraction", type=float, default=0.20, help="Held-out fraction for calibration when enough labels exist")
    parser.add_argument("--class-weight", default="balanced", choices=["balanced"], help="LogisticRegression class_weight")
    parser.add_argument("--max-iter", type=int, default=1000, help="LogisticRegression max_iter")
    parser.add_argument("--n-estimators", type=int, default=60, help="Number of trees for --backend sklearn_gbm")
    parser.add_argument("--learning-rate", type=float, default=0.05, help="Learning rate for --backend sklearn_gbm")
    parser.add_argument("--max-depth", type=int, default=2, help="Tree max_depth for --backend sklearn_gbm")
    parser.add_argument("--random-state", type=int, default=13, help="Random seed")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    model = train_model(args)
    with open(args.output, "w") as out:
        json.dump(model, out, indent=2, sort_keys=True)
        out.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
