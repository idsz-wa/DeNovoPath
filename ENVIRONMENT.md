# Environment and Reproducibility Notes

## Python Environment

Recommended:

```bash
conda create -n denovopath python=3.10 -y
conda activate denovopath
python -m pip install --upgrade pip
python -m pip install -e .
```

Minimal runtime dependencies:

```bash
python -m pip install -r requirements.txt
```

Optional development and ML-training dependencies:

```bash
python -m pip install -r requirements-dev.txt
```

## Required Input Types

DeNovoPath scoring requires:

- VCF or compressed VCF (`.vcf`, `.vcf.gz`)
- reference genome FASTA
- GFF/GFF3/GTF-like gene annotation
- output VCF path

Recommended for richer annotation:

- CDS FASTA
- protein FASTA
- YAML scoring configuration

Optional:

- sample metadata
- protein domain table
- protein structure feature table
- protein language-model score table
- miRNA target-site interval table
- portable JSON ML model

## Local Directory Convention

The GitHub repository should not include large datasets. Use:

```text
data/       # local inputs, ignored by Git
results/    # local outputs, ignored by Git
tmp/        # temporary shard files, ignored by Git
```

## Verification

Run:

```bash
python -m unittest tests.test_denovopath
```

The tests build synthetic FASTA, GFF, VCF and optional annotation resources in a temporary directory. They are intended to verify core annotation, scoring, CLI dispatch, ranked export, REF validation, HTML reporting and parallel fast-mode behavior without external data.
