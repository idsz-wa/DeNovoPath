# Data Placement

This repository intentionally excludes large input and result files.

Place local input files under `data/`:

```text
data/
  input.vcf.gz
  reference.fa
  genes.gff3
  cds.fa
  proteins.fa
```

Place generated outputs under `results/`:

```text
results/
  input.denovopath.vcf.gz
  input.denovopath.summary.json
  input.denovopath.report.html
  input.denovopath.top.tsv
  input.denovopath.genes.tsv
```

Both directories are ignored by `.gitignore`.

For reproducible publications, deposit large datasets in a public archive such as Zenodo, Figshare, NCBI SRA, ENA or a field-specific repository, then document the accession or DOI in the repository README or release notes.
