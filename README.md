# DeNovoPath

DeNovoPath is a reference-local VCF annotation and prioritization framework for population-scale variants, with a focus on non-model species where curated variant-effect databases are often unavailable.

The software annotates VCF records using local genome resources and writes interpretable `DNP_*` INFO fields. It combines transcript consequence annotation, protein-effect proxies, splice and regulatory context, sequence-context evidence, cohort statistics, quality-control signals and confidence estimates into variant-level scores and gene-level ranking outputs.

## Key Features

- Reference-local annotation from FASTA, GFF, CDS FASTA and protein FASTA.
- Consequence annotation for coding, splice-adjacent, UTR, promoter, intronic and intergenic variants.
- Interpretable `DNP_SCORE`, `DNP_LEVEL`, confidence and evidence component fields.
- Cohort-aware features including AC, AN, AF, carrier counts, heterozygosity and HWE-related summaries.
- REF validation against the supplied reference genome.
- REF/ALT normalization utility for biallelic SNP VCFs with provisional reference alleles.
- Optional record-sharded parallel fast mode for large VCFs.
- Ranked variant and gene TSV export.
- Optional portable ML model training and inference hooks.

## Repository Layout

```text
.
├── config/                  # Default scoring configuration
├── scripts/                 # Command-line helper scripts
├── src/denovopath/          # Python package source
├── tests/                   # Unit and smoke tests with synthetic data
├── pyproject.toml           # Python package metadata
├── requirements.txt         # Runtime dependencies
├── requirements-dev.txt     # Optional development dependencies
└── README.md
```

Large VCF, FASTA, GFF and result files are intentionally not included in this upload package. Put local datasets under `data/` and generated outputs under `results/`; both directories are ignored by Git.

## Installation

Create and activate a Python environment. Python 3.10 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Runtime dependencies are:

- `pyfaidx`
- `PyYAML`

For optional ML training support:

```bash
python -m pip install -r requirements-dev.txt
```

If you use conda, the equivalent workflow is:

```bash
conda create -n denovopath python=3.10 -y
conda activate denovopath
python -m pip install -e .
```

## Quick Start

Prepare the following input files:

- VCF or VCF.GZ
- reference genome FASTA
- GFF/GFF3/GTF-like annotation
- CDS FASTA
- protein FASTA
- optional YAML config, or use `config/default.yaml`

Run single-process scoring:

```bash
denovopath score-vcf \
  --vcf data/input.vcf.gz \
  --reference data/reference.fa \
  --gff data/genes.gff3 \
  --cds data/cds.fa \
  --pep data/proteins.fa \
  --config config/default.yaml \
  --summary results/input.denovopath.summary.json \
  --html-report results/input.denovopath.report.html \
  --output results/input.denovopath.vcf.gz
```

The same command can be run through the legacy script entry point:

```bash
python scripts/score_vcf.py \
  --vcf data/input.vcf.gz \
  --reference data/reference.fa \
  --gff data/genes.gff3 \
  --cds data/cds.fa \
  --pep data/proteins.fa \
  --config config/default.yaml \
  --summary results/input.denovopath.summary.json \
  --output results/input.denovopath.vcf.gz
```

## Parallel Fast Mode

For large VCFs, DeNovoPath can split records into shards and score shards in parallel. Parallel sharding currently requires disabling cross-record window statistics and gene-constraint prescans:

```bash
denovopath parallel-score-vcf \
  --vcf data/input.vcf.gz \
  --reference data/reference.fa \
  --gff data/genes.gff3 \
  --cds data/cds.fa \
  --pep data/proteins.fa \
  --config config/default.yaml \
  --pop-window 0 \
  --skip-gene-constraint \
  --jobs 4 \
  --records-per-shard 50000 \
  --temp-dir /tmp \
  --summary results/input.parallel.summary.json \
  --html-report results/input.parallel.report.html \
  --output results/input.parallel.denovopath.vcf.gz
```

This mode preserves the primary annotation, scoring, allele-frequency, HWE, heterozygosity, QC and confidence outputs. It disables optional fixed-window population-diversity statistics and gene-constraint prescans.

## REF/ALT Normalization

Some population VCFs contain provisional REF alleles. For biallelic SNPs, DeNovoPath provides a normalization helper:

```bash
python scripts/normalize_biallelic_snp_ref.py \
  --vcf data/input.vcf.gz \
  --reference data/reference.fa \
  --output results/input.refnorm.vcf.gz \
  --summary-out results/input.refnorm.summary.json
```

Records are marked with `REFNORM_STATUS`:

- `match`: VCF REF already matches the FASTA base.
- `swap`: VCF ALT matches the FASTA base, so REF and ALT were swapped and genotype allele codes were flipped.
- `unfixable`: neither REF nor ALT matches the FASTA base.
- `non_snp_or_multiallelic`: record was not a biallelic SNP.
- `missing_contig`: contig was not found in the FASTA.

Existing external annotation fields such as `ANN` are retained, but they may describe the pre-normalized ALT allele after a swap. Re-run the external annotator on the normalized VCF when strict allele-level comparison is required.

## Ranked Output Export

After scoring, export ranked variants and optional gene-level summaries:

```bash
denovopath export-ranked \
  --vcf results/input.denovopath.vcf.gz \
  --variants-out results/input.denovopath.top.tsv \
  --genes-out results/input.denovopath.genes.tsv \
  --top 1000
```

## Command Overview

After `python -m pip install -e .`, these entry points are available:

```bash
denovopath score-vcf --help
denovopath parallel-score-vcf --help
denovopath predict --help
denovopath export-ranked --help
denovopath train --help
denovopath-parallel --help
```

`denovopath predict` is an alias for `score-vcf`.

## Testing

The test suite uses synthetic in-memory data and does not require external datasets:

```bash
python -m unittest tests.test_denovopath
```

## Data Policy

Do not commit large genome, VCF or result files to this repository. Use the following local layout:

```text
data/
  input.vcf.gz
  reference.fa
  genes.gff3
  cds.fa
  proteins.fa

results/
  input.denovopath.vcf.gz
  input.denovopath.summary.json
  input.denovopath.top.tsv
  input.denovopath.genes.tsv
```

For public releases, deposit large datasets in an appropriate data repository and cite the accession or DOI in project documentation.

## License

MIT license 

## Citation

If you use DeNovoPath in a publication, cite the corresponding paper, repository release or archived software record once available.
