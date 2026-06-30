# GitHub Release Checklist

Before pushing this directory to GitHub:

1. Confirm the repository contains only source code, tests, configuration and documentation.
2. Confirm large datasets and generated outputs are not present:

   ```bash
   find . -type f \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.fa" -o -name "*.fasta" -o -name "*.gff" -o -name "*.gff3" -o -name "*.bam" -o -name "*.cram" \)
   ```

3. Run tests:

   ```bash
   python -m unittest tests.test_denovopath
   ```

4. Build a wheel locally:

   ```bash
   python -m pip wheel --no-build-isolation --no-deps --wheel-dir /tmp/denovopath_wheel_check .
   ```

5. Review `README.md`, `ENVIRONMENT.md`, `DATA.md` and `LICENSE`.
6. Create the GitHub repository and push:

   ```bash
   git init
   git add .
   git commit -m "Initial DeNovoPath release"
   git branch -M main
   git remote add origin <repository-url>
   git push -u origin main
   ```

7. After pushing, create a tagged release when the public version is ready.
