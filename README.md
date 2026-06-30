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

## Understanding DeNovoPath Scores

`DNP_SCORE` is a `0-1` variant-prioritization score. A higher value means that the variant has stronger reference-local evidence for functional impact and should be considered earlier during candidate review. It is a ranking score, not a calibrated pathogenicity probability and not a substitute for experimental validation.

The default score is a weighted summary of several evidence layers defined in `config/default.yaml`:

| Evidence layer | Main field | Default weight | Meaning |
|---|---:|---:|---|
| Functional consequence | `DNP_IMPACT` | 0.42 | Severity implied by the annotated consequence, such as intergenic, promoter, missense, stop-gained, frameshift or splice-site classes |
| Protein impact | `DNP_PROT` | 0.28 | Coding and protein-level evidence, including amino-acid substitution severity, codon usage, protein context, structure, domain and optional protein-model inputs |
| Splice impact | `DNP_SPLICE` | 0.15 | Splice-boundary, donor/acceptor motif, PWM, MaxEnt-like, branch-point, polypyrimidine-tract and ESE/ESS evidence |
| Sequence context | `DNP_SEQ` | 0.10 | Local sequence-context evidence, including k-mer changes, repeat/low-complexity context, mutation context and DNA k-mer language-model proxy signals |
| Cohort evidence | `DNP_COHORT` | sample-count dependent | GT-derived rarity, carrier pattern, allele frequency and population-statistic evidence |

The cohort evidence weight depends on sample count:

| VCF context | Cohort weight |
|---|---:|
| Single-sample VCF | 0.02 |
| Small-cohort VCF | 0.05 |
| Large-cohort VCF | 0.12 |

This sample-count gating is intentional. Genotype and population statistics are weak evidence in a single sample, but become more useful when a cohort is large enough to support frequency, heterozygosity, HWE, group differentiation or case/control summaries.

### Score Levels

`DNP_LEVEL` is a coarse interpretation of `DNP_SCORE`:

| Level | Default score range | Suggested interpretation |
|---|---:|---|
| `HIGH` | `>= 0.80` | High-priority candidate for manual review or validation |
| `MODERATE` | `0.50-0.80` | Candidate worth retaining in broader screens, especially with relevant gene or phenotype context |
| `LOW` | `0.20-0.50` | Lower-priority candidate; may still be useful in targeted gene analyses |
| `MINIMAL` | `< 0.20` | Lowest-priority under the available evidence |

### Main Evidence Fields

Use the component fields to understand why a variant scored highly:

| Field | What it indicates |
|---|---|
| `DNP_CONSEQ` | Strongest predicted consequence per ALT allele |
| `DNP_GENE`, `DNP_TX` | Gene and transcript supporting the strongest consequence |
| `DNP_IMPACT` | Rule-based impact from consequence and region context |
| `DNP_PROT` | Combined protein-level impact score |
| `DNP_SPLICE` | Combined splice-impact score |
| `DNP_SEQ` | Combined local sequence-context score |
| `DNP_COHORT` | Genotype cohort rarity/carrier-pattern score |
| `DNP_QC` | Genotype/call-quality support score |
| `DNP_CONF` | Confidence score combining `DNP_QC` and score separation |
| `DNP_AC`, `DNP_AN`, `DNP_AF` | GT-derived ALT allele count, allele number and allele frequency |
| `DNP_CARR`, `DNP_HET`, `DNP_HOMALT`, `DNP_MISS` | Carrier, heterozygous, homozygous-ALT and missing-genotype counts |

Protein subfields include `DNP_GRANTHAM`, `DNP_BLOSUM`, `DNP_CODONUSE`, `DNP_PROTCTX`, `DNP_STRUCT`, `DNP_DOMAIN`, `DNP_AFSTRUCT`, `DNP_ESM` and `DNP_PROTLM`. Splice subfields include `DNP_SPLICE_MOTIF`, `DNP_SPLICE_PWM`, `DNP_SPLICE_MAXENT`, `DNP_SPLICE_AUX` and `DNP_SPLICE_ESE`. Regulatory and sequence-context subfields include `DNP_PROM`, `DNP_UTR`, `DNP_RNAFOLD`, `DNP_MIRNA`, `DNP_KMER`, `DNP_REPEAT`, `DNP_MUTCTX`, `DNP_DNALM` and `DNP_96CTX`.

Optional cohort/group fields include `DNP_HWE`, `DNP_HETOBS`, `DNP_HETEXP`, `DNP_HETDEV`, `DNP_FIS`, `DNP_FST`, `DNP_CASECTRL`, `DNP_SUBAF`, `DNP_PRIVATE`, `DNP_PI`, `DNP_THETA`, `DNP_TAJD`, `DNP_LD`, `DNP_HAP`, `DNP_GENELOF`, `DNP_GENEMIS` and `DNP_GENECON`. These fields are only meaningful when the VCF and optional metadata support them.

### Practical Interpretation

A recommended review order is:

1. Sort or filter by `DNP_LEVEL` and `DNP_SCORE`.
2. Inspect `DNP_CONSEQ` to determine whether the signal is coding, splice, regulatory or non-coding.
3. Check the main component fields: `DNP_IMPACT`, `DNP_PROT`, `DNP_SPLICE`, `DNP_SEQ` and `DNP_COHORT`.
4. Use `DNP_QC` and `DNP_CONF` to remove low-quality or low-confidence candidates.
5. Interpret candidates together with `DNP_GENE`, `DNP_TX`, `DNP_AF`, carrier counts, gene annotation and downstream biological context.

Examples:

- A `HIGH` variant with `DNP_CONSEQ=stop_gained_early`, high `DNP_IMPACT`, high `DNP_PROT`, low `DNP_AF` and high `DNP_CONF` is a strong loss-of-function candidate.
- A `MODERATE` missense variant with high `DNP_PROT`, high `DNP_GRANTHAM` or high `DNP_DOMAIN` may be worth retaining if the gene is biologically relevant.
- A high `DNP_SPLICE` variant near a donor or acceptor motif can be important even if it is not a missense or stop-gained variant.
- A promoter or UTR variant may have a lower total score than a coding loss-of-function event, but can still be useful in targeted regulatory analyses if `DNP_PROM`, `DNP_UTR` or related regulatory evidence is high.
- A high `DNP_SCORE` with low `DNP_QC` or low `DNP_CONF` should be treated cautiously until the original genotype and read-level evidence are reviewed.

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
