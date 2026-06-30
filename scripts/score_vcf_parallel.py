#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from denovopath.scorer import build_arg_parser as build_score_arg_parser
from denovopath.scorer import maybe_index_vcf, open_text, render_html_report


@dataclass(frozen=True)
class Shard:
    index: int
    input_path: str
    output_path: str
    summary_path: str
    records: int


def _open_input(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _open_output(path: str):
    return gzip.open(path, "wt", compresslevel=1) if path.endswith(".gz") else open(path, "w")


def split_vcf_by_records(
    input_vcf: str,
    temp_dir: str,
    records_per_shard: int,
    limit: int = 0,
) -> List[Shard]:
    if records_per_shard <= 0:
        raise ValueError("--records-per-shard must be positive")

    shards: List[Shard] = []
    header_lines: List[str] = []
    writer = None
    shard_records = 0
    total_records = 0

    def close_writer() -> None:
        nonlocal writer
        if writer is not None:
            writer.close()
            writer = None

    def open_next_shard() -> None:
        nonlocal writer, shard_records
        close_writer()
        idx = len(shards)
        input_path = os.path.join(temp_dir, f"shard_{idx:05d}.vcf.gz")
        output_path = os.path.join(temp_dir, f"shard_{idx:05d}.scored.vcf.gz")
        summary_path = os.path.join(temp_dir, f"shard_{idx:05d}.summary.json")
        shards.append(Shard(idx, input_path, output_path, summary_path, 0))
        writer = _open_output(input_path)
        for header in header_lines:
            writer.write(header)
        shard_records = 0

    try:
        with _open_input(input_vcf) as handle:
            for line in handle:
                if line.startswith("#"):
                    header_lines.append(line)
                    continue
                if limit and total_records >= limit:
                    break
                if writer is None or shard_records >= records_per_shard:
                    open_next_shard()
                assert writer is not None
                writer.write(line)
                shard_records += 1
                total_records += 1
                current = shards[-1]
                shards[-1] = Shard(
                    current.index,
                    current.input_path,
                    current.output_path,
                    current.summary_path,
                    current.records + 1,
                )
    finally:
        close_writer()

    return [shard for shard in shards if shard.records > 0]


def build_worker_args(args: argparse.Namespace, shard: Shard) -> Dict[str, object]:
    return {
        "vcf": shard.input_path,
        "reference": args.reference,
        "gff": args.gff,
        "cds": args.cds,
        "pep": args.pep,
        "protein_domains": getattr(args, "protein_domains", None),
        "protein_structures": getattr(args, "protein_structures", None),
        "protein_lm_scores": getattr(args, "protein_lm_scores", None),
        "mirna_sites": getattr(args, "mirna_sites", None),
        "ml_model": getattr(args, "ml_model", None),
        "sample_info": getattr(args, "sample_info", None),
        "output": shard.output_path,
        "window": args.window,
        "pop_window": args.pop_window,
        "skip_gene_constraint": args.skip_gene_constraint,
        "config": args.config,
        "region": args.region,
        "summary": shard.summary_path,
        "html_report": None,
        "index_output": "never",
        "limit": 0,
        "progress": 0,
        "strict_ref": getattr(args, "strict_ref", False),
    }


def run_worker(worker_args: Dict[str, object]) -> None:
    from denovopath.scorer import score_vcf

    score_vcf(argparse.Namespace(**worker_args))


def merge_vcf_outputs(shards: Sequence[Shard], output_path: str) -> None:
    wrote_header = False
    with open_text(output_path, "w") as out:
        for shard in shards:
            with _open_input(shard.output_path) as handle:
                for line in handle:
                    if line.startswith("#"):
                        if not wrote_header:
                            out.write(line)
                            if line.startswith("#CHROM"):
                                wrote_header = True
                        continue
                    out.write(line)


def _sum_counter(summaries: Iterable[Dict[str, object]], key: str) -> Dict[str, int]:
    counter: Counter[str] = Counter()
    for summary in summaries:
        data = summary.get(key, {})
        if isinstance(data, dict):
            counter.update({str(k): int(v) for k, v in data.items()})
    return dict(counter.most_common())


def _reference_counts(summaries: Sequence[Dict[str, object]]) -> Dict[str, object]:
    checked = matching = mismatch = unchecked = 0
    examples: List[str] = []
    for summary in summaries:
        ref = summary.get("validation", {}).get("reference", {}) if isinstance(summary.get("validation"), dict) else {}
        if not isinstance(ref, dict):
            continue
        checked += int(ref.get("checked_records", 0))
        matching += int(ref.get("matching_records", 0))
        mismatch += int(ref.get("mismatch_records", 0))
        unchecked += int(ref.get("unchecked_records", 0))
        for item in ref.get("mismatch_examples", []) or []:
            if len(examples) < 10:
                examples.append(str(item))
    return {
        "checked_records": checked,
        "matching_records": matching,
        "mismatch_records": mismatch,
        "unchecked_records": unchecked,
        "mismatch_examples": examples,
    }


def aggregate_summaries(
    args: argparse.Namespace,
    shard_summaries: Sequence[Dict[str, object]],
    shards: Sequence[Shard],
    index_result: Dict[str, object],
    elapsed_seconds: float,
) -> Dict[str, object]:
    if not shard_summaries:
        raise ValueError("No shard summaries to aggregate")
    first = dict(shard_summaries[0])
    input_records = sum(int(item.get("input_records_seen", 0)) for item in shard_summaries)
    records = sum(int(item.get("records_scored", 0)) for item in shard_summaries)
    skipped = sum(int(item.get("records_skipped_by_region", 0)) for item in shard_summaries)
    alt_records = sum(int(item.get("alt_alleles_scored", 0)) for item in shard_summaries)
    prescan_seconds = sum(float(item.get("benchmark", {}).get("prescan_seconds", 0.0)) for item in shard_summaries)
    scoring_seconds = sum(float(item.get("benchmark", {}).get("scoring_seconds", 0.0)) for item in shard_summaries)
    worker_elapsed = max(float(item.get("benchmark", {}).get("elapsed_seconds", 0.0)) for item in shard_summaries)
    worker_peak_values = [
        float(item.get("benchmark", {}).get("peak_rss_mb", 0.0))
        for item in shard_summaries
        if isinstance(item.get("benchmark"), dict)
    ]
    max_worker_rss = max(worker_peak_values) if worker_peak_values else 0.0
    mean_worker_rss = sum(worker_peak_values) / len(worker_peak_values) if worker_peak_values else 0.0
    concurrent_workers = min(args.jobs, len(shards))
    estimated_parallel_peak_rss = max_worker_rss * concurrent_workers
    records_per_second = records / elapsed_seconds if elapsed_seconds > 0 else 0.0
    alt_per_second = alt_records / elapsed_seconds if elapsed_seconds > 0 else 0.0

    benchmark = dict(first.get("benchmark", {})) if isinstance(first.get("benchmark"), dict) else {}
    benchmark.update(
        {
            "elapsed_seconds": round(elapsed_seconds, 6),
            "prescan_seconds": round(prescan_seconds, 6),
            "scoring_seconds": round(scoring_seconds, 6),
            "index_seconds": float(index_result.get("index_seconds", 0.0) or 0.0),
            "records_per_second": round(records_per_second, 3),
            "alt_alleles_per_second": round(alt_per_second, 3),
            "worker_max_elapsed_seconds": round(worker_elapsed, 6),
            "worker_scoring_seconds_sum": round(scoring_seconds, 6),
            "peak_rss_mb": round(estimated_parallel_peak_rss, 3),
            "worker_peak_rss_mb_max": round(max_worker_rss, 3),
            "worker_peak_rss_mb_mean": round(mean_worker_rss, 3),
            "parallel_peak_rss_note": "Estimated as max worker RSS multiplied by concurrent worker count.",
        }
    )

    first.update(
        {
            "input_vcf": args.vcf,
            "output_vcf": args.output,
            "config": args.config,
            "regions": args.region or [],
            "input_records_seen": input_records,
            "records_scored": records,
            "records_skipped_by_region": skipped,
            "alt_alleles_scored": alt_records,
            "index": {key: value for key, value in index_result.items() if key != "index_seconds"},
            "benchmark": benchmark,
            "level_counts": _sum_counter(shard_summaries, "level_counts"),
            "consequence_counts": _sum_counter(shard_summaries, "consequence_counts"),
            "validation": {"reference": _reference_counts(shard_summaries)},
            "parallel": {
                "enabled": True,
                "jobs": args.jobs,
                "records_per_shard": args.records_per_shard,
                "shards": len(shards),
                "shard_records": [shard.records for shard in shards],
                "mode": "record_chunk_fast_mode",
                "requires_pop_window_zero": True,
                "requires_skip_gene_constraint": True,
            },
        }
    )
    return first


def load_json(path: str) -> Dict[str, object]:
    with open(path) as handle:
        return json.load(handle)


def score_vcf_parallel(args: argparse.Namespace) -> None:
    if args.jobs <= 1:
        from denovopath.scorer import score_vcf

        score_vcf(args)
        return
    if args.pop_window != 0 or not args.skip_gene_constraint:
        raise ValueError(
            "Parallel sharding currently requires --pop-window 0 --skip-gene-constraint "
            "so per-shard results are identical to single-process fast mode."
        )

    start = time.perf_counter()
    base_temp_dir = args.temp_dir or os.path.dirname(os.path.abspath(args.output)) or "."
    os.makedirs(base_temp_dir, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix="denovopath_shards_", dir=base_temp_dir)
    try:
        shards = split_vcf_by_records(args.vcf, work_dir, args.records_per_shard, limit=args.limit)
        if not shards:
            raise ValueError("Input VCF contains no records to score")
        worker_args = [build_worker_args(args, shard) for shard in shards]
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_shard: Dict[Future[None], Shard] = {
                executor.submit(run_worker, item): shard for item, shard in zip(worker_args, shards)
            }
            for future in as_completed(future_to_shard):
                shard = future_to_shard[future]
                try:
                    future.result()
                except Exception as exc:
                    raise RuntimeError(f"Shard {shard.index} failed") from exc

        merge_vcf_outputs(shards, args.output)
        index_start = time.perf_counter()
        index_result = maybe_index_vcf(args.output, args.index_output)
        index_result["index_seconds"] = round(time.perf_counter() - index_start, 6)
        elapsed = time.perf_counter() - start
        shard_summaries = [load_json(shard.summary_path) for shard in shards]
        summary = aggregate_summaries(args, shard_summaries, shards, index_result, elapsed)
        if args.summary:
            with open(args.summary, "w") as out:
                json.dump(summary, out, indent=2, sort_keys=True)
                out.write("\n")
        if args.html_report:
            render_html_report(summary, args.html_report)
    finally:
        if args.keep_shards:
            print(f"Kept shard work directory: {work_dir}", file=sys.stderr)
        else:
            shutil.rmtree(work_dir, ignore_errors=True)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = build_score_arg_parser()
    parser.description = "Annotate VCF INFO with DeNovoPath scores using record-sharded parallel fast mode."
    parser.add_argument("--jobs", type=int, default=2, help="Number of shard workers to run in parallel")
    parser.add_argument(
        "--records-per-shard",
        type=int,
        default=50000,
        help="Number of input records per temporary shard",
    )
    parser.add_argument("--temp-dir", help="Directory for temporary shard input/output files")
    parser.add_argument("--keep-shards", action="store_true", help="Keep temporary shard files for debugging")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    score_vcf_parallel(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
