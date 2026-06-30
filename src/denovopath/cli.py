from __future__ import annotations

import argparse
import sys
from typing import Callable, Optional, Sequence

from denovopath import __version__
from denovopath import scorer
from scripts import export_ranked, score_vcf_parallel, train_ml_model


CommandMain = Callable[[Optional[Sequence[str]]], int]


def _dispatch(command_main: CommandMain, argv: Sequence[str]) -> int:
    return command_main(list(argv))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="denovopath",
        description="DeNovoPath command-line interface.",
    )
    parser.add_argument("--version", action="version", version=f"denovopath {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    score_parser = subparsers.add_parser(
        "score-vcf",
        help="Annotate a VCF with DeNovoPath DNP_* INFO fields.",
        add_help=False,
    )
    score_parser.add_argument("args", nargs=argparse.REMAINDER)

    predict_parser = subparsers.add_parser(
        "predict",
        help="Alias for score-vcf, matching the target PRD CLI.",
        add_help=False,
    )
    predict_parser.add_argument("args", nargs=argparse.REMAINDER)

    parallel_parser = subparsers.add_parser(
        "parallel-score-vcf",
        help="Annotate a VCF with record-sharded parallel fast-mode scoring.",
        add_help=False,
    )
    parallel_parser.add_argument("args", nargs=argparse.REMAINDER)

    train_parser = subparsers.add_parser(
        "train",
        help="Train a portable JSON ML model from a scored VCF.",
        add_help=False,
    )
    train_parser.add_argument("args", nargs=argparse.REMAINDER)

    export_parser = subparsers.add_parser(
        "export-ranked",
        help="Export ranked variant and optional gene TSVs from a scored VCF.",
        add_help=False,
    )
    export_parser.add_argument("args", nargs=argparse.REMAINDER)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in {"score-vcf", "predict"}:
        return _dispatch(scorer.main, argv[1:])
    if argv and argv[0] == "parallel-score-vcf":
        return _dispatch(score_vcf_parallel.main, argv[1:])
    if argv and argv[0] == "train":
        return _dispatch(train_ml_model.main, argv[1:])
    if argv and argv[0] == "export-ranked":
        return _dispatch(export_ranked.main, argv[1:])
    parser = build_arg_parser()
    parser.parse_args(argv)
    parser.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
