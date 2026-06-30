from __future__ import annotations

import gzip
import json
import os
import sys
import tempfile
import unittest


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from denovopath.scorer import DeNovoPathScorer, MlModel, SampleInfo, ScoreConfig, active_methods, genotype_counts, parse_region, score_vcf  # noqa: E402
from denovopath.cli import main as denovopath_cli_main  # noqa: E402
from scripts.export_ranked import export_ranked, parse_info  # noqa: E402
from scripts.score_vcf_parallel import score_vcf_parallel  # noqa: E402
from scripts.train_ml_model import main as train_ml_model_main  # noqa: E402


class Args:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class DeNovoPathTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.work = self.tmp.name
        self.reference = os.path.join(self.work, "ref.fa")
        self.gff = os.path.join(self.work, "genes.gff")
        self.cds = os.path.join(self.work, "cds.fa")
        self.pep = os.path.join(self.work, "pep.fa")
        self.domains = os.path.join(self.work, "domains.tsv")
        self.structures = os.path.join(self.work, "structures.tsv")
        self.esm_scores = os.path.join(self.work, "esm_scores.tsv")
        self.mirna_sites = os.path.join(self.work, "mirna_sites.tsv")
        self.ml_model = os.path.join(self.work, "ml_model.json")
        self.vcf = os.path.join(self.work, "input.vcf")
        self.scored_vcf = os.path.join(self.work, "scored.vcf.gz")
        self.summary = os.path.join(self.work, "summary.json")
        self.html_report = os.path.join(self.work, "summary.html")
        self.html_only_vcf = os.path.join(self.work, "html_only.vcf.gz")
        self.html_only_report = os.path.join(self.work, "html_only.html")
        self.variants = os.path.join(self.work, "variants.tsv")
        self.genes = os.path.join(self.work, "genes.tsv")
        chr5_seq = list("A" * 52)
        chr5_seq[0:3] = list("ATG")
        chr5_seq[19:24] = list("CTAAC")
        chr5_seq[33:45] = list("T" * 12)
        chr5_seq[47:49] = list("AG")
        chr5_seq[49:52] = list("AAA")
        chr6_seq = "AAAGAAGAATTT"
        chr9_seq = "ATG" + "GAA" * 10 + "TAA"

        # Positions 10-18 are the CDS ATGTCAAAG for transcript tx1.
        with open(self.reference, "w") as out:
            out.write(">chr1\n")
            out.write("ACGTACGTAATGTCAAAGACGTACGTACGTACGTACGT\n")
            out.write(">chr2\n")
            out.write("ATGGTCCAGAAA\n")
            out.write(">chr3\n")
            out.write("AAAGAAATGCCC\n")
            out.write(">chr4\n")
            out.write("AAATATAAACCCCCCCCCCATGTCAAAG\n")
            out.write(">chr5\n")
            out.write("".join(chr5_seq) + "\n")
            out.write(">chr6\n")
            out.write(chr6_seq + "\n")
            out.write(">chr7\n")
            out.write("ATGAAACCCGGG\n")
            out.write(">chr8\n")
            out.write("ATATATATATATATATATATATAT\n")
            out.write(">chr9\n")
            out.write(chr9_seq + "\n")
            out.write(">chr10\n")
            out.write("ATGAAACCCGGGTAAATTTT\n")
            out.write(">chr11\n")
            out.write("ATGTCAAAG\n")
        with open(self.gff, "w") as out:
            out.write("chr1\tDNP\tgene\t10\t18\t.\t+\t.\tID=gene1\n")
            out.write("chr1\tDNP\tmRNA\t10\t18\t.\t+\t.\tID=tx1;Parent=gene1\n")
            out.write("chr1\tDNP\tCDS\t10\t18\t.\t+\t0\tID=cds1;Parent=tx1\n")
            out.write("chr1\tDNP\tgene\t20\t22\t.\t+\t.\tID=gene2\n")
            out.write("chr1\tDNP\tmRNA\t20\t22\t.\t+\t.\tID=tx2;Parent=gene2\n")
            out.write("chr1\tDNP\tCDS\t20\t22\t.\t+\t0\tID=cds2;Parent=tx2\n")
            out.write("chr2\tDNP\tgene\t1\t12\t.\t+\t.\tID=gene3\n")
            out.write("chr2\tDNP\tmRNA\t1\t12\t.\t+\t.\tID=tx3;Parent=gene3\n")
            out.write("chr2\tDNP\tCDS\t1\t3\t.\t+\t0\tID=cds3a;Parent=tx3\n")
            out.write("chr2\tDNP\tCDS\t10\t12\t.\t+\t0\tID=cds3b;Parent=tx3\n")
            out.write("chr3\tDNP\tgene\t1\t12\t.\t+\t.\tID=gene4\n")
            out.write("chr3\tDNP\tmRNA\t1\t12\t.\t+\t.\tID=tx4;Parent=gene4\n")
            out.write("chr3\tDNP\tfive_prime_UTR\t1\t6\t.\t+\t.\tID=utr4;Parent=tx4\n")
            out.write("chr3\tDNP\tCDS\t7\t12\t.\t+\t0\tID=cds4;Parent=tx4\n")
            out.write("chr4\tDNP\tgene\t20\t28\t.\t+\t.\tID=gene5\n")
            out.write("chr4\tDNP\tmRNA\t20\t28\t.\t+\t.\tID=tx5;Parent=gene5\n")
            out.write("chr4\tDNP\tCDS\t20\t28\t.\t+\t0\tID=cds5;Parent=tx5\n")
            out.write("chr5\tDNP\tgene\t1\t52\t.\t+\t.\tID=gene6\n")
            out.write("chr5\tDNP\tmRNA\t1\t52\t.\t+\t.\tID=tx6;Parent=gene6\n")
            out.write("chr5\tDNP\tCDS\t1\t3\t.\t+\t0\tID=cds6a;Parent=tx6\n")
            out.write("chr5\tDNP\tCDS\t50\t52\t.\t+\t0\tID=cds6b;Parent=tx6\n")
            out.write("chr6\tDNP\tgene\t1\t12\t.\t+\t.\tID=gene7\n")
            out.write("chr6\tDNP\tmRNA\t1\t12\t.\t+\t.\tID=tx7;Parent=gene7\n")
            out.write("chr6\tDNP\tCDS\t1\t12\t.\t+\t0\tID=cds7;Parent=tx7\n")
            out.write("chr7\tDNP\tgene\t1\t12\t.\t+\t.\tID=gene8\n")
            out.write("chr7\tDNP\tmRNA\t1\t6\t.\t+\t.\tID=tx8;Parent=gene8\n")
            out.write("chr7\tDNP\tCDS\t1\t6\t.\t+\t0\tID=cds8;Parent=tx8\n")
            out.write("chr7\tDNP\tmRNA\t1\t12\t.\t+\t.\tID=tx9;Parent=gene8\n")
            out.write("chr7\tDNP\tCDS\t1\t12\t.\t+\t0\tID=cds9;Parent=tx9\n")
            out.write("chr9\tDNP\tgene\t1\t36\t.\t+\t.\tID=gene9\n")
            out.write("chr9\tDNP\tmRNA\t1\t36\t.\t+\t.\tID=tx10;Parent=gene9\n")
            out.write("chr9\tDNP\tCDS\t1\t36\t.\t+\t0\tID=cds10;Parent=tx10\n")
            out.write("chr10\tDNP\tgene\t1\t15\t.\t-\t.\tID=gene10\n")
            out.write("chr10\tDNP\tmRNA\t1\t15\t.\t-\t.\tID=tx11;Parent=gene10\n")
            out.write("chr10\tDNP\tCDS\t1\t15\t.\t-\t0\tID=cds11;Parent=tx11\n")
            out.write("chr11\tDNP\tgene\t1\t9\t.\t+\t.\tID=gene11\n")
            out.write("chr11\tDNP\tCDS\t1\t9\t.\t+\t0\tID=cds12;Parent=gene11\n")
        with open(self.cds, "w") as out:
            out.write(">tx1\nATGTCAAAG\n")
            out.write(">tx2\nTGA\n")
            out.write(">tx3\nATGAAA\n")
            out.write(">tx4\nATGCCC\n")
            out.write(">tx5\nATGTCAAAG\n")
            out.write(">tx6\nATGAAA\n")
            out.write(">tx7\nAAAGAAGAATTT\n")
            out.write(">tx8\nATGAAA\n")
            out.write(">tx9\nATGAAACCCGGG\n")
            out.write(">tx10\n" + chr9_seq + "\n")
            out.write(">tx11\nATGAAACCCGGGTAA\n")
            out.write(">gene11\nATGTCAAAG\n")
        with open(self.pep, "w") as out:
            out.write(">tx1\nMSKPPPPPPPPPPSSSSS\n")
            out.write(">tx2\n*\n")
            out.write(">tx3\nMK\n")
            out.write(">tx5\nMSK\n")
            out.write(">tx6\nMK\n")
            out.write(">tx7\nKEEF\n")
            out.write(">tx8\nMK\n")
            out.write(">tx9\nMKPG\n")
            out.write(">tx10\nMEEEEEEEEEE*\n")
            out.write(">tx11\nMKPG*\n")
            out.write(">gene11\nMSK\n")
        with open(self.domains, "w") as out:
            out.write("transcript\tstart\tend\tdomain\tscore\n")
            out.write("tx1\t2\t3\tkinase_like\t0.9\n")
        with open(self.structures, "w") as out:
            out.write("transcript\taa_pos\tplddt\trsa\tss\tsource\n")
            out.write("tx10\t2\t95\t0.05\tH\tAlphaFold_like\n")
        with open(self.esm_scores, "w") as out:
            out.write("transcript\taa_pos\tref\talt\tscore\tlabel\n")
            out.write("tx10\t2\tE\tV\t0.87\tESM2_test\n")
            out.write("tx1\t2\tS\t*\t0.91\tESM2_stop\n")
        with open(self.mirna_sites, "w") as out:
            out.write("chrom\tstart\tend\ttranscript\tmirna\tsite\tscore\n")
            out.write("chr3\t2\t5\ttx4\tmiR156\tseed_2_5\t0.88\n")
        with open(self.ml_model, "w") as out:
            json.dump(
                {
                    "features": ["impact_score", "protein_score", "splice_score", "sequence_score", "qc_score"],
                    "coef": [2.0, 1.0, 0.8, 0.4, 0.2],
                    "intercept": -0.3,
                    "mean": [0.0, 0.0, 0.0, 0.0, 0.0],
                    "scale": [1.0, 1.0, 1.0, 1.0, 1.0],
                    "calibration": {"scale": 1.0, "shift": 0.0},
                },
                out,
            )
        with open(self.vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
            out.write("chr1\t14\tstop\tC\tG\t100\tPASS\t.\tGT:DP:AD:GQ\t0/1:20:10,10:60\n")
            out.write("chr1\t15\tsyn\tA\tG\t100\tPASS\t.\tGT:DP:AD:GQ\t0/1:20:10,10:60\n")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_coding_consequences_are_ranked(self) -> None:
        scorer = DeNovoPathScorer(self.reference, self.gff, cds_fasta=self.cds, protein_fasta=self.pep)
        stop = scorer.score_record(
            ["chr1", "14", "stop", "C", "G", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        start_lost = scorer.score_record(
            ["chr1", "10", "start", "A", "T", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        synonymous = scorer.score_record(
            ["chr1", "15", "syn", "A", "G", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        multi_codon_mnv = scorer.score_record(
            ["chr1", "15", "mnv2codon", "AA", "GC", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        inframe_insertion = scorer.score_record(
            ["chr1", "15", "ins", "A", "AGGG", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        inframe_deletion = scorer.score_record(
            ["chr1", "15", "del", "ATCA", "A", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        stop_retained = scorer.score_record(
            ["chr1", "22", "stop_retained", "G", "A", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        stop_gained_early = scorer.score_record(
            ["chr9", "4", "stop_early", "G", "T", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        stop_gained_terminal = scorer.score_record(
            ["chr9", "31", "stop_terminal", "G", "T", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        stop_lost_readthrough = scorer.score_record(
            ["chr9", "34", "stop_lost", "T", "C", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        structural_missense = scorer.score_record(
            ["chr9", "5", "struct", "A", "T", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        exon_boundary = scorer.score_record(
            ["chr2", "2", "exon_boundary", "TGG", "T", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        splice_motif = scorer.score_record(
            ["chr2", "4", "splice", "G", "A", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        utr5 = scorer.score_record(
            ["chr3", "4", "utr5", "G", "C", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        promoter = scorer.score_record(
            ["chr4", "4", "prom", "T", "C", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        branchpoint = scorer.score_record(
            ["chr5", "20", "branch", "C", "A", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        exonic_splice = scorer.score_record(
            ["chr6", "5", "ese", "A", "C", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        repeat_context = scorer.score_record(
            ["chr8", "7", "repeat", "A", "G", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        self.assertEqual(stop.consequence, "stop_gained")
        self.assertEqual(start_lost.consequence, "start_lost")
        self.assertEqual(synonymous.consequence, "synonymous")
        self.assertEqual(multi_codon_mnv.consequence, "missense")
        self.assertEqual(multi_codon_mnv.aa_change, "p.S2S|p.K3Q")
        self.assertEqual(multi_codon_mnv.codon_change, "TCA>TCG|AAG>CAG")
        self.assertEqual(multi_codon_mnv.hgvs_change, "tx1:c.6A>G|c.7A>C:p.S2S|p.K3Q")
        self.assertEqual(inframe_insertion.consequence, "inframe_insertion")
        self.assertEqual(inframe_deletion.consequence, "inframe_deletion")
        self.assertEqual(inframe_insertion.normalized_variant, "16:->GGG")
        self.assertEqual(inframe_deletion.normalized_variant, "15:ATC>-")
        self.assertEqual(inframe_deletion.hgvs_change, "tx1:c.6_8del:p.?")
        self.assertEqual(stop_retained.consequence, "stop_retained")
        self.assertEqual(stop_gained_early.consequence, "stop_gained_early")
        self.assertEqual(stop_gained_terminal.consequence, "stop_gained_terminal")
        self.assertEqual(stop_lost_readthrough.consequence, "stop_lost_readthrough")
        self.assertGreater(stop_gained_early.score, stop_gained_terminal.score)
        self.assertEqual(structural_missense.consequence, "missense")
        self.assertGreater(structural_missense.protein_structure_score, 0.0)
        self.assertEqual(exon_boundary.consequence, "exon_boundary_disruption")
        self.assertEqual(exon_boundary.normalized_variant, "3:GG>-")
        self.assertGreaterEqual(exon_boundary.splice_score, 0.90)
        self.assertIn("exon_boundary", exon_boundary.hgvs_change)
        self.assertEqual(splice_motif.consequence, "splice_acceptor_donor")
        self.assertEqual(utr5.consequence, "utr5")
        self.assertEqual(promoter.consequence, "promoter")
        self.assertEqual(branchpoint.consequence, "splice_region")
        self.assertIn(exonic_splice.consequence, {"missense", "synonymous"})
        self.assertEqual(stop.level, "HIGH")
        self.assertGreater(stop.score, synonymous.score)
        self.assertGreater(start_lost.grantham_score, 0.0)
        self.assertGreater(start_lost.blosum_score, 0.0)
        self.assertGreater(synonymous.codon_usage_score, 0.0)
        self.assertGreater(multi_codon_mnv.grantham_score, 0.0)
        self.assertGreater(synonymous.protein_context_score, 0.2)
        self.assertGreater(structural_missense.protein_lm_score, 0.0)
        self.assertGreater(stop.kmer_score, 0.0)
        self.assertEqual(stop.mutation_context, "T[C>G]A")
        self.assertGreater(stop.mutation_context_score, 0.0)
        self.assertGreaterEqual(stop.dna_lm_score, 0.0)
        self.assertGreaterEqual(splice_motif.splice_motif_score, 0.9)
        self.assertGreater(splice_motif.splice_pwm_score, 0.0)
        self.assertGreater(splice_motif.splice_maxent_score, 0.2)
        self.assertGreater(branchpoint.splice_aux_score, 0.25)
        self.assertGreater(exonic_splice.splice_ese_score, 0.25)
        self.assertGreater(repeat_context.repeat_score, 0.25)
        self.assertGreater(utr5.utr_score, 0.2)
        self.assertGreater(utr5.rnafold_score, 0.0)
        self.assertGreater(promoter.promoter_score, 0.25)
        self.assertEqual(promoter.gene_id, "gene5")
        self.assertEqual(promoter.tx_id, "tx5")
        self.assertEqual(stop.gene_id, "gene1")
        self.assertEqual(stop.tx_id, "tx1")
        self.assertEqual(stop.hgvs_change, "tx1:c.5C>G:p.S2*")

        domain_scorer = DeNovoPathScorer(
            self.reference,
            self.gff,
            cds_fasta=self.cds,
            protein_fasta=self.pep,
            protein_domains=self.domains,
            protein_structures=self.structures,
            protein_lm_scores=self.esm_scores,
            mirna_sites=self.mirna_sites,
        )
        structured_missense = domain_scorer.score_record(
            ["chr9", "5", "struct", "A", "T", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        self.assertGreater(structured_missense.protein_structure_model_score, 0.0)
        self.assertIn("AlphaFold_like", structured_missense.protein_structure_label)
        self.assertGreater(structured_missense.protein_esm_score, 0.8)
        self.assertIn("ESM2_test", structured_missense.protein_esm_label)
        domain_stop = domain_scorer.score_record(
            ["chr1", "14", "stop", "C", "G", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        self.assertGreater(domain_stop.protein_domain_score, 0.8)
        self.assertEqual(domain_stop.protein_domain_label, "kinase_like:2-3")
        self.assertGreater(domain_stop.protein_esm_score, 0.8)
        self.assertIn("ESM2_stop", domain_stop.protein_esm_label)
        mirna_utr = domain_scorer.score_record(
            ["chr3", "4", "utr5_mirna", "G", "C", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        self.assertGreater(mirna_utr.mirna_score, 0.8)
        self.assertEqual(mirna_utr.mirna_label, "miR156:seed_2_5")

    def test_promoter_scoring_handles_short_scaffold_out_of_bounds_variant(self) -> None:
        scorer = DeNovoPathScorer(self.reference, self.gff, cds_fasta=self.cds, protein_fasta=self.pep)
        score = scorer.score_record(
            ["chr10", "30", "short_scaffold_promoter", "A", "C", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        self.assertEqual(score.consequence, "promoter")
        self.assertEqual(score.gene_id, "gene10")
        self.assertEqual(score.tx_id, "tx11")
        self.assertGreaterEqual(score.promoter_score, 0.0)

    def test_gff_with_cds_parent_gene_without_mrna_is_annotated(self) -> None:
        scorer = DeNovoPathScorer(self.reference, self.gff, cds_fasta=self.cds, protein_fasta=self.pep)
        score = scorer.score_record(
            ["chr11", "5", "gene_parent_cds", "C", "G", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        self.assertEqual(score.gene_id, "gene11")
        self.assertEqual(score.tx_id, "gene11")
        self.assertIn(score.consequence, {"missense", "synonymous", "stop_gained"})
        self.assertNotEqual(score.consequence, "intergenic")

    def test_configurable_transcript_priority(self) -> None:
        scorer = DeNovoPathScorer(
            self.reference,
            self.gff,
            cds_fasta=self.cds,
            protein_fasta=self.pep,
            config=ScoreConfig(transcript_priority="longest_cds"),
        )
        selected = scorer.score_record(
            ["chr7", "6", "shared", "A", "G", "100", "PASS", ".", "GT:DP:AD:GQ", "0/1:20:10,10:60"],
            n_samples=1,
        )[0]
        self.assertEqual(selected.consequence, "synonymous")
        self.assertEqual(selected.tx_id, "tx9")
        self.assertIn("gene8:tx8:synonymous", selected.all_transcripts)
        self.assertIn("gene8:tx9:synonymous", selected.all_transcripts)

    def test_region_scoring_summary_and_rank_export(self) -> None:
        score_vcf(
            Args(
                vcf=self.vcf,
                reference=self.reference,
                gff=self.gff,
                cds=self.cds,
                pep=self.pep,
                protein_domains=self.domains,
                protein_structures=self.structures,
                protein_lm_scores=self.esm_scores,
                mirna_sites=self.mirna_sites,
                ml_model=self.ml_model,
                output=self.scored_vcf,
                window=10,
                config=None,
                region=["chr1:10-14"],
                summary=self.summary,
                html_report=self.html_report,
                index_output="never",
                limit=0,
                progress=0,
            )
        )
        with gzip.open(self.scored_vcf, "rt") as handle:
            lines = handle.readlines()
        records = [line for line in lines if not line.startswith("#")]
        self.assertEqual(len(records), 1)
        self.assertIn("DNP_CONSEQ=stop_gained", records[0])
        self.assertIn("DNP_SCORE=", records[0])
        self.assertIn("DNP_ALLTX=", records[0])
        self.assertIn("DNP_HGVS=", records[0])
        self.assertIn("DNP_NORM=", records[0])
        self.assertIn("DNP_GRANTHAM=", records[0])
        self.assertIn("DNP_BLOSUM=", records[0])
        self.assertIn("DNP_CODONUSE=", records[0])
        self.assertIn("DNP_PROTCTX=", records[0])
        self.assertIn("DNP_STRUCT=", records[0])
        self.assertIn("DNP_AFSTRUCT=", records[0])
        self.assertIn("DNP_ESM=", records[0])
        self.assertIn("DNP_PROTLM=", records[0])
        self.assertIn("DNP_DOMAIN=", records[0])
        self.assertIn("DNP_DOMID=kinase_like:2-3", records[0])
        self.assertIn("DNP_AFID=", records[0])
        self.assertIn("DNP_ESMID=ESM2_stop:S2*", records[0])
        self.assertIn("DNP_SPLICE_MOTIF=", records[0])
        self.assertIn("DNP_SPLICE_PWM=", records[0])
        self.assertIn("DNP_SPLICE_MAXENT=", records[0])
        self.assertIn("DNP_SPLICE_AUX=", records[0])
        self.assertIn("DNP_SPLICE_ESE=", records[0])
        self.assertIn("DNP_UTR=", records[0])
        self.assertIn("DNP_RNAFOLD=", records[0])
        self.assertIn("DNP_MIRNA=", records[0])
        self.assertIn("DNP_MIRID=", records[0])
        self.assertIn("DNP_PROM=", records[0])
        self.assertIn("DNP_KMER=", records[0])
        self.assertIn("DNP_REPEAT=", records[0])
        self.assertIn("DNP_MUTCTX=", records[0])
        self.assertIn("DNP_DNALM=", records[0])
        self.assertIn("DNP_96CTX=", records[0])
        self.assertIn("DNP_MAFBIN=", records[0])
        self.assertIn("DNP_HWE=", records[0])
        self.assertIn("DNP_HETOBS=", records[0])
        self.assertIn("DNP_HETEXP=", records[0])
        self.assertIn("DNP_HETDEV=", records[0])
        self.assertIn("DNP_FIS=", records[0])
        self.assertIn("DNP_SUBAF=", records[0])
        self.assertIn("DNP_PRIVATE=", records[0])
        self.assertIn("DNP_FST=", records[0])
        self.assertIn("DNP_CASEAF=", records[0])
        self.assertIn("DNP_CASECTRL=", records[0])
        self.assertIn("DNP_PI=", records[0])
        self.assertIn("DNP_THETA=", records[0])
        self.assertIn("DNP_TAJD=", records[0])
        self.assertIn("DNP_LD=", records[0])
        self.assertIn("DNP_HAP=", records[0])
        self.assertIn("DNP_GENELOF=", records[0])
        self.assertIn("DNP_GENEMIS=", records[0])
        self.assertIn("DNP_GENECON=", records[0])
        self.assertIn("DNP_ML=", records[0])
        self.assertIn("DNP_CAL=", records[0])
        self.assertIn("DNP_UNCERT=", records[0])
        self.assertIn("DNP_OOD=", records[0])
        self.assertIn("DNP_FEATIMP=", records[0])
        self.assertTrue(any(line.startswith("##DeNovoPathActiveMethods=") for line in lines))
        with open(self.summary) as handle:
            summary_text = handle.read()
        self.assertIn('"records_scored": 1', summary_text)
        self.assertIn('"records_skipped_by_region": 1', summary_text)
        self.assertIn('"validation":', summary_text)
        self.assertIn('"mismatch_records": 0', summary_text)
        self.assertIn('"method_source_types":', summary_text)
        self.assertIn('"sample_count_gating":', summary_text)
        self.assertIn('"benchmark":', summary_text)
        self.assertIn('"elapsed_seconds":', summary_text)
        self.assertIn('"peak_rss_mb":', summary_text)
        self.assertIn('"index":', summary_text)
        self.assertIn('"status": "disabled"', summary_text)
        self.assertIn("cohort_single_sample_gt", summary_text)
        self.assertIn("protein_domain_annotation", summary_text)
        self.assertIn('"protein_domain_records": 1', summary_text)
        self.assertIn("protein_structure_alphafold_esmfold_annotation", summary_text)
        self.assertIn('"protein_structure_records": 1', summary_text)
        self.assertIn("protein_esm2_precomputed_delta", summary_text)
        self.assertIn('"protein_lm_score_records": 2', summary_text)
        self.assertIn("regulatory_mirna_seed_disruption", summary_text)
        self.assertIn('"mirna_site_records": 1', summary_text)
        self.assertIn("ml_json_model_inference", summary_text)
        self.assertIn('"ml_model_loaded": true', summary_text)
        with open(self.html_report) as handle:
            html_report = handle.read()
        self.assertIn("DeNovoPath Run Report", html_report)
        self.assertIn("Reference Validation", html_report)
        self.assertIn("Method Source Types", html_report)
        self.assertIn("Impact Level Distribution", html_report)
        self.assertIn("Top Consequence Classes", html_report)
        self.assertIn("Filtering Strategy Guide", html_report)
        self.assertIn("Strict high-confidence triage", html_report)
        self.assertIn("Sample-count gating", html_report)
        self.assertIn("Active Methods", html_report)
        self.assertIn("benchmark", html_report)
        self.assertIn("cohort_single_sample_gt", html_report)
        self.assertIn("protein_domain_annotation", html_report)
        self.assertIn("protein_structure_alphafold_esmfold_annotation", html_report)
        self.assertIn("protein_esm2_precomputed_delta", html_report)
        self.assertIn("regulatory_mirna_seed_disruption", html_report)
        self.assertIn("ml_json_model_inference", html_report)
        score_vcf(
            Args(
                vcf=self.vcf,
                reference=self.reference,
                gff=self.gff,
                cds=self.cds,
                pep=self.pep,
                output=self.html_only_vcf,
                window=10,
                config=None,
                region=["chr1:10-14"],
                summary=None,
                html_report=self.html_only_report,
                index_output="never",
                limit=0,
                progress=0,
            )
        )
        self.assertTrue(os.path.exists(self.html_only_report))

        export_ranked(
            Args(
                vcf=self.scored_vcf,
                variants_out=self.variants,
                genes_out=self.genes,
                min_score=0.0,
                min_qc=0.0,
                min_ac=1,
                top=1,
            )
        )
        with open(self.variants) as handle:
            variant_lines = handle.readlines()
        self.assertEqual(len(variant_lines), 2)
        self.assertIn("\tgrantham\tblosum\tcodon_usage\tprotein_context\t", variant_lines[0])
        self.assertIn(
            "\tprotein_context\tprotein_structure\tprotein_structure_model\tprotein_esm\tprotein_lm\tprotein_domain\tsplice\t",
            variant_lines[0],
        )
        self.assertIn(
            "\tsplice\tsplice_motif\tsplice_pwm\tsplice_maxent\tsplice_aux\tsplice_ese\tutr\trnafold\tmirna\tpromoter\tsequence\tkmer\t",
            variant_lines[0],
        )
        self.assertIn(
            "\tkmer\trepeat\tmutation_context_score\tdna_lm\tcohort\thwe\theterozygosity_observed\theterozygosity_expected\theterozygosity_deviation\tfis\tfst\tcase_control\tpi\ttheta\ttajima_d\tld\thaplotype\tgene_lof_oe\tgene_missense_oe\tgene_constraint\tqc\tconfidence\tml\tcalibrated\tuncertainty\tood\t",
            variant_lines[0],
        )
        self.assertIn(
            "\thgvs_like\tnormalized_variant\tprotein_domain_id\tprotein_structure_id\tprotein_esm_id\tmirna_id\tfeature_importance\tcontext_96\tsubpopulation_af\tprivate_shared\tcase_control_af\t",
            variant_lines[0],
        )
        self.assertIn("\taf\tmaf_bin\tcarriers", variant_lines[0])
        self.assertIn("\tgene\ttranscript\tall_transcripts\taa_change\t", variant_lines[0])
        self.assertIn("\tcodon_change\thgvs_like\tnormalized_variant\t", variant_lines[0])
        self.assertIn("\tstop_gained\tgene1\t", variant_lines[1])

    def test_cli_score_vcf_wrapper(self) -> None:
        cli_output = os.path.join(self.work, "cli.scored.vcf.gz")
        cli_summary = os.path.join(self.work, "cli.summary.json")
        rc = denovopath_cli_main(
            [
                "score-vcf",
                "--vcf",
                self.vcf,
                "--reference",
                self.reference,
                "--gff",
                self.gff,
                "--cds",
                self.cds,
                "--pep",
                self.pep,
                "--output",
                cli_output,
                "--region",
                "chr1:10-14",
                "--summary",
                cli_summary,
                "--index-output",
                "never",
                "--progress",
                "0",
            ]
        )
        self.assertEqual(rc, 0)
        with gzip.open(cli_output, "rt") as handle:
            self.assertIn("##DeNovoPathActiveMethods=", handle.read())
        with open(cli_summary) as handle:
            summary = json.load(handle)
        self.assertEqual(summary["records_scored"], 1)
        self.assertIn("cohort_single_sample_gt", summary["active_methods"])

    def test_reference_mismatch_is_reported_and_can_be_strict(self) -> None:
        bad_vcf = os.path.join(self.work, "bad_ref.vcf")
        bad_output = os.path.join(self.work, "bad_ref.scored.vcf.gz")
        bad_summary = os.path.join(self.work, "bad_ref.summary.json")
        with open(bad_vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
            out.write("chr1\t14\tbad\tA\tG\t100\tPASS\t.\tGT\t0/1\n")
        args = Args(
            vcf=bad_vcf,
            reference=self.reference,
            gff=self.gff,
            cds=self.cds,
            pep=self.pep,
            output=bad_output,
            window=10,
            config=None,
            region=None,
            summary=bad_summary,
            html_report=None,
            index_output="never",
            limit=0,
            progress=0,
        )
        score_vcf(args)
        with open(bad_summary) as handle:
            summary = json.load(handle)
        self.assertEqual(summary["validation"]["reference"]["mismatch_records"], 1)
        self.assertIn("chr1:14", summary["validation"]["reference"]["mismatch_examples"][0])
        args.strict_ref = True
        with self.assertRaises(ValueError):
            score_vcf(args)

    def test_train_ml_model_from_explicit_labels(self) -> None:
        output = os.path.join(self.work, "train_input.scored.vcf.gz")
        labels = os.path.join(self.work, "labels.tsv")
        model_json = os.path.join(self.work, "trained_model.json")
        gbm_model_json = os.path.join(self.work, "trained_gbm_model.json")
        score_vcf(
            Args(
                vcf=self.vcf,
                reference=self.reference,
                gff=self.gff,
                cds=self.cds,
                pep=self.pep,
                output=output,
                window=10,
                config=None,
                region=None,
                summary=None,
                html_report=None,
                index_output="never",
                limit=0,
                progress=0,
            )
        )
        with open(labels, "w") as out:
            out.write("chrom\tpos\tref\talt\tlabel\n")
            out.write("chr1\t14\tC\tG\t1\n")
            out.write("chr1\t15\tA\tG\t0\n")
        rc = train_ml_model_main(
            [
                "--vcf",
                output,
                "--labels",
                labels,
                "--output",
                model_json,
                "--features",
                "impact_score",
                "protein_score",
                "qc_score",
                "confidence_score",
                "--min-positive",
                "1",
                "--min-negative",
                "1",
                "--validation-fraction",
                "0",
            ]
        )
        self.assertEqual(rc, 0)
        with open(model_json) as handle:
            model = json.load(handle)
        self.assertEqual(model["model_type"], "sklearn_logistic_regression_json")
        self.assertEqual(model["features"], ["impact_score", "protein_score", "qc_score", "confidence_score"])
        self.assertEqual(len(model["coef"]), 4)
        self.assertIn("feature_importance", model)
        self.assertEqual(model["training_summary"]["label_source_counts"]["explicit"], 2)
        self.assertEqual(model["training_summary"]["label_source_counts"]["pseudo"], 0)

        rc = train_ml_model_main(
            [
                "--vcf",
                output,
                "--labels",
                labels,
                "--output",
                gbm_model_json,
                "--backend",
                "sklearn_gbm",
                "--features",
                "impact_score",
                "protein_score",
                "qc_score",
                "confidence_score",
                "--min-positive",
                "1",
                "--min-negative",
                "1",
                "--validation-fraction",
                "0",
                "--n-estimators",
                "3",
                "--max-depth",
                "1",
            ]
        )
        self.assertEqual(rc, 0)
        with open(gbm_model_json) as handle:
            gbm_model = json.load(handle)
        self.assertEqual(gbm_model["model_type"], "sklearn_gradient_boosting_json")
        self.assertEqual(gbm_model["features"], ["impact_score", "protein_score", "qc_score", "confidence_score"])
        self.assertGreater(len(gbm_model["trees"]), 0)
        self.assertIn("metrics", gbm_model)
        self.assertIn("backend", gbm_model["metrics"])
        ml_model = MlModel(gbm_model_json)
        ml_score, calibrated, uncertainty, ood = ml_model.predict(
            {"impact_score": 0.9, "protein_score": 0.8, "qc_score": 0.9, "confidence_score": 0.8}
        )
        self.assertGreaterEqual(ml_score, 0.0)
        self.assertLessEqual(calibrated, 1.0)
        self.assertGreaterEqual(uncertainty, 0.0)
        self.assertGreaterEqual(ood, 0.0)

    def test_train_ml_model_validation_metrics(self) -> None:
        scored_vcf = os.path.join(self.work, "metric_input.scored.vcf")
        labels = os.path.join(self.work, "metric_labels.tsv")
        model_json = os.path.join(self.work, "metric_model.json")
        with open(scored_vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
            for idx in range(50):
                if idx % 2 == 0:
                    info = "DNP_SCORE=0.9000;DNP_IMPACT=0.9000;DNP_PROT=0.8000;DNP_QC=0.9000;DNP_CONF=0.8500;DNP_CONSEQ=stop_gained"
                    out.write(f"chr1\t14\tpos{idx}\tC\tG\t100\tPASS\t{info}\tGT\t0/1\n")
                else:
                    info = "DNP_SCORE=0.0500;DNP_IMPACT=0.0500;DNP_PROT=0.0500;DNP_QC=0.9000;DNP_CONF=0.8500;DNP_CONSEQ=synonymous"
                    out.write(f"chr1\t15\tneg{idx}\tA\tG\t100\tPASS\t{info}\tGT\t0/1\n")
        with open(labels, "w") as out:
            out.write("chrom\tpos\tref\talt\tlabel\n")
            out.write("chr1\t14\tC\tG\t1\n")
            out.write("chr1\t15\tA\tG\t0\n")
        rc = train_ml_model_main(
            [
                "--vcf",
                scored_vcf,
                "--labels",
                labels,
                "--output",
                model_json,
                "--features",
                "impact_score",
                "protein_score",
                "qc_score",
                "confidence_score",
                "--min-positive",
                "1",
                "--min-negative",
                "1",
                "--validation-fraction",
                "0.25",
            ]
        )
        self.assertEqual(rc, 0)
        with open(model_json) as handle:
            model = json.load(handle)
        self.assertEqual(model["metrics"]["validation_records"], 13)
        self.assertIn("validation_auc", model["metrics"])
        self.assertIn("validation_pr_auc", model["metrics"])
        self.assertIn("validation_mcc", model["metrics"])
        self.assertIn("validation_brier", model["metrics"])
        self.assertEqual(model["training_summary"]["label_source_counts"]["explicit"], 50)

    def test_sample_count_aware_maf_bins(self) -> None:
        scorer = DeNovoPathScorer(self.reference, self.gff, cds_fasta=self.cds, protein_fasta=self.pep)
        one_sample = scorer.score_record(
            ["chr1", "14", "one", "C", "G", "100", "PASS", ".", "GT", "0/1"],
            n_samples=1,
        )[0]
        small_cohort = scorer.score_record(
            ["chr1", "14", "small", "C", "G", "100", "PASS", ".", "GT", "0/1", "0/0", "0/0", "0/0"],
            n_samples=4,
        )[0]
        singleton = scorer.score_record(
            ["chr1", "14", "singleton", "C", "G", "100", "PASS", ".", "GT"]
            + ["0/1"]
            + ["0/0"] * 9,
            n_samples=10,
        )[0]
        common = scorer.score_record(
            ["chr1", "14", "common", "C", "G", "100", "PASS", ".", "GT"]
            + ["0/1"] * 5
            + ["0/0"] * 5,
            n_samples=10,
        )[0]
        fixed = scorer.score_record(
            ["chr1", "14", "fixed", "C", "G", "100", "PASS", ".", "GT"] + ["1/1"] * 10,
            n_samples=10,
        )[0]
        hwe_excess_het = scorer.score_record(
            ["chr1", "14", "hwe", "C", "G", "100", "PASS", ".", "GT"] + ["0/1"] * 10,
            n_samples=10,
        )[0]

        self.assertEqual(one_sample.maf_bin, "single_sample")
        self.assertEqual(small_cohort.maf_bin, "small_cohort")
        self.assertEqual(singleton.maf_bin, "singleton")
        self.assertEqual(common.maf_bin, "common")
        self.assertEqual(fixed.maf_bin, "fixed_or_near_fixed")
        self.assertGreater(hwe_excess_het.hwe_score, 0.4)
        self.assertAlmostEqual(hwe_excess_het.heterozygosity_observed, 1.0)
        self.assertAlmostEqual(hwe_excess_het.heterozygosity_expected, 0.5)
        self.assertAlmostEqual(hwe_excess_het.heterozygosity_deviation_score, 0.5)
        self.assertAlmostEqual(hwe_excess_het.inbreeding_coefficient, -1.0)
        self.assertEqual(singleton.ac, 1)
        self.assertEqual(singleton.an, 20)
        self.assertGreater(singleton.cohort_score, fixed.cohort_score)
        self.assertIn("annotation_stop_altering_subclasses", active_methods(1))
        self.assertIn("annotation_exon_boundary_spanning", active_methods(1))
        self.assertIn("protein_lm_kmer_delta_proxy", active_methods(1))
        self.assertIn("utr_rnafold_delta_g_heuristic", active_methods(1))
        self.assertIn("sequence_dna_lm_kmer_delta", active_methods(1))
        self.assertIn("cohort_small_sample_carrier_pattern", active_methods(4))
        self.assertIn("cohort_large_sample_maf_binning", active_methods(10))
        self.assertIn("cohort_large_sample_hwe_deviation", active_methods(10))
        self.assertIn("cohort_large_sample_heterozygosity", active_methods(10))
        self.assertIn("cohort_large_sample_inbreeding_coefficient", active_methods(10))
        self.assertIn("cohort_window_pi_theta", active_methods(10))
        self.assertIn("cohort_window_tajima_d", active_methods(10))
        self.assertIn("cohort_window_ld_haplotype_proxy", active_methods(10))
        self.assertIn("cohort_gene_constraint_proxy", active_methods(10))
        self.assertIn("protein_domain_annotation", active_methods(1, has_protein_domains=True))
        self.assertIn("protein_structure_alphafold_esmfold_annotation", active_methods(1, has_protein_structures=True))
        self.assertIn("protein_esm2_precomputed_delta", active_methods(1, has_protein_esm_scores=True))
        self.assertIn("regulatory_mirna_seed_disruption", active_methods(1, has_mirna_sites=True))
        self.assertIn("ml_json_model_inference", active_methods(1, has_ml_model=True))

        self.assertEqual(DeNovoPathScorer.maf_bin_label(0, 20000, 10000), "absent")
        self.assertEqual(DeNovoPathScorer.maf_bin_label(2, 20000, 10000), "ultra_rare")
        self.assertEqual(DeNovoPathScorer.maf_bin_label(50, 20000, 10000), "rare")
        self.assertEqual(DeNovoPathScorer.maf_bin_label(500, 20000, 10000), "low_frequency")
        self.assertEqual(DeNovoPathScorer.maf_bin_label(19999, 20000, 10000), "fixed_or_near_fixed")

    def test_sample_info_group_population_scores(self) -> None:
        scorer = DeNovoPathScorer(self.reference, self.gff, cds_fasta=self.cds, protein_fasta=self.pep)
        sample_info = SampleInfo(
            groups={"s1": "A", "s2": "A", "s3": "B", "s4": "B"},
            phenotypes={"s1": "case", "s2": "case", "s3": "control", "s4": "control"},
        )
        grouped = scorer.score_record(
            ["chr1", "14", "grouped", "C", "G", "100", "PASS", ".", "GT", "0/1", "0/1", "0/0", "0/0"],
            n_samples=4,
            sample_names=["s1", "s2", "s3", "s4"],
            sample_info=sample_info,
        )[0]

        self.assertEqual(grouped.group_af, "A:0.5000|B:0.0000")
        self.assertEqual(grouped.private_shared, "private:A")
        self.assertGreater(grouped.fst_score, 0.3)
        self.assertEqual(grouped.case_control_af, "case:0.5000|control:0.0000")
        self.assertAlmostEqual(grouped.case_control_score, 0.5)
        self.assertIn("cohort_group_subpopulation_af", active_methods(4, has_sample_info=True, has_phenotype=True))
        self.assertIn("cohort_group_private_shared", active_methods(4, has_sample_info=True, has_phenotype=True))
        self.assertIn("cohort_group_fst", active_methods(4, has_sample_info=True, has_phenotype=True))
        self.assertIn("cohort_case_control_enrichment", active_methods(4, has_sample_info=True, has_phenotype=True))

    def test_sample_info_case_control_with_single_group(self) -> None:
        scorer = DeNovoPathScorer(self.reference, self.gff, cds_fasta=self.cds, protein_fasta=self.pep)
        sample_info = SampleInfo(
            groups={"s1": "A", "s2": "A", "s3": "A", "s4": "A"},
            phenotypes={"s1": "case", "s2": "case", "s3": "control", "s4": "control"},
        )
        grouped = scorer.score_record(
            ["chr1", "14", "single_group", "C", "G", "100", "PASS", ".", "GT", "0/1", "0/1", "0/0", "0/0"],
            n_samples=4,
            sample_names=["s1", "s2", "s3", "s4"],
            sample_info=sample_info,
        )[0]

        self.assertEqual(grouped.group_af, ".")
        self.assertEqual(grouped.private_shared, ".")
        self.assertEqual(grouped.fst_score, 0.0)
        self.assertEqual(grouped.case_control_af, "case:0.5000|control:0.0000")
        self.assertAlmostEqual(grouped.case_control_score, 0.5)

    def test_multiallelic_genotype_classes_exclude_other_alts(self) -> None:
        counts = genotype_counts(["GT"], ["0/0", "0/1", "1/1", "2/2", "1/2"], 2)

        self.assertEqual(counts[0]["ac"], 4)
        self.assertEqual(counts[0]["an"], 10)
        self.assertEqual(counts[0]["carriers"], 3)
        self.assertEqual(counts[0]["n_called_diploid"], 3)
        self.assertEqual(counts[0]["n_hom_ref"], 1)
        self.assertEqual(counts[0]["n_het"], 1)
        self.assertEqual(counts[0]["n_hom_alt"], 1)

        self.assertEqual(counts[1]["ac"], 3)
        self.assertEqual(counts[1]["an"], 10)
        self.assertEqual(counts[1]["carriers"], 2)
        self.assertEqual(counts[1]["n_called_diploid"], 2)
        self.assertEqual(counts[1]["n_hom_ref"], 1)
        self.assertEqual(counts[1]["n_het"], 0)
        self.assertEqual(counts[1]["n_hom_alt"], 1)

    def test_large_cohort_window_population_statistics(self) -> None:
        vcf = os.path.join(self.work, "window.vcf")
        output = os.path.join(self.work, "window.scored.vcf.gz")
        summary = os.path.join(self.work, "window.summary.json")
        samples = [f"s{i}" for i in range(10)]
        with open(vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            genotypes = ["1/1"] * 5 + ["0/0"] * 5
            out.write("chr1\t14\tw1\tC\tG\t100\tPASS\t.\tGT\t" + "\t".join(genotypes) + "\n")
            out.write("chr1\t15\tw2\tA\tG\t100\tPASS\t.\tGT\t" + "\t".join(genotypes) + "\n")
        score_vcf(
            Args(
                vcf=vcf,
                reference=self.reference,
                gff=self.gff,
                cds=self.cds,
                pep=self.pep,
                output=output,
                window=10,
                pop_window=1000,
                config=None,
                sample_info=None,
                region=None,
                summary=summary,
                html_report=None,
                index_output="never",
                limit=0,
                progress=0,
            )
        )
        with gzip.open(output, "rt") as handle:
            text = handle.read()
        self.assertIn("cohort_window_pi_theta", text)
        self.assertIn("cohort_window_tajima_d", text)
        self.assertIn("cohort_window_ld_haplotype_proxy", text)
        self.assertIn("DNP_PI=1.0000", text)
        self.assertRegex(text, r"DNP_THETA=0\.[0-9]+")
        self.assertIn("DNP_TAJD=", text)
        self.assertIn("DNP_LD=1.0000", text)
        self.assertIn("DNP_HAP=0.0000", text)
        self.assertIn("DNP_GENELOF=", text)
        self.assertIn("DNP_GENEMIS=", text)
        self.assertIn("DNP_GENECON=", text)
        self.assertRegex(text, r"DNP_GENECON=0\.[0-9]+")
        with open(summary) as handle:
            summary_text = handle.read()
        self.assertIn('"population_window_size": 1000', summary_text)
        self.assertIn('"gene_constraint_genes": 1', summary_text)
        self.assertIn('"records_per_second":', summary_text)
        self.assertIn('"n_samples": 10', summary_text)

    def test_fast_cohort_mode_skips_optional_prescans(self) -> None:
        vcf = os.path.join(self.work, "fast.vcf")
        output = os.path.join(self.work, "fast.scored.vcf.gz")
        summary = os.path.join(self.work, "fast.summary.json")
        samples = [f"s{i}" for i in range(10)]
        with open(vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            out.write("chr1\t14\tfast1\tC\tG\t100\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 10) + "\n")
        score_vcf(
            Args(
                vcf=vcf,
                reference=self.reference,
                gff=self.gff,
                cds=self.cds,
                pep=self.pep,
                output=output,
                window=10,
                pop_window=0,
                skip_gene_constraint=True,
                config=None,
                sample_info=None,
                region=None,
                summary=summary,
                html_report=None,
                index_output="never",
                limit=0,
                progress=0,
            )
        )
        with gzip.open(output, "rt") as handle:
            text = handle.read()
        self.assertNotIn("cohort_window_pi_theta", text)
        self.assertNotIn("cohort_window_tajima_d", text)
        self.assertNotIn("cohort_window_ld_haplotype_proxy", text)
        self.assertNotIn("cohort_gene_constraint_proxy", text)
        self.assertIn("cohort_large_sample_hwe_deviation", text)
        self.assertIn("DNP_PI=0.0000", text)
        self.assertIn("DNP_THETA=0.0000", text)
        self.assertIn("DNP_GENECON=0.0000", text)
        with open(summary) as handle:
            summary_data = json.load(handle)
        self.assertFalse(summary_data["population_windows_enabled"])
        self.assertFalse(summary_data["gene_constraint_enabled"])
        self.assertEqual(summary_data["gene_constraint_genes"], 0)
        self.assertNotIn("cohort_window_pi_theta", summary_data["active_methods"])
        self.assertNotIn("cohort_gene_constraint_proxy", summary_data["active_methods"])

    def test_parallel_fast_mode_matches_single_process_records(self) -> None:
        vcf = os.path.join(self.work, "parallel.vcf")
        single_output = os.path.join(self.work, "parallel.single.vcf.gz")
        parallel_output = os.path.join(self.work, "parallel.sharded.vcf.gz")
        single_summary = os.path.join(self.work, "parallel.single.summary.json")
        parallel_summary = os.path.join(self.work, "parallel.sharded.summary.json")
        samples = [f"s{i}" for i in range(10)]
        with open(vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            out.write("chr1\t14\tp1\tC\tG\t100\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 10) + "\n")
            out.write("chr1\t15\tp2\tA\tG\t100\tPASS\t.\tGT\t" + "\t".join(["0/0"] * 5 + ["0/1"] * 5) + "\n")
            out.write("chr4\t4\tp3\tT\tC\t100\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 3 + ["0/0"] * 7) + "\n")

        common = dict(
            vcf=vcf,
            reference=self.reference,
            gff=self.gff,
            cds=self.cds,
            pep=self.pep,
            protein_domains=None,
            protein_structures=None,
            protein_lm_scores=None,
            mirna_sites=None,
            ml_model=None,
            sample_info=None,
            window=10,
            pop_window=0,
            skip_gene_constraint=True,
            config=None,
            region=None,
            html_report=None,
            index_output="never",
            limit=0,
            progress=0,
            strict_ref=False,
        )
        score_vcf(Args(**common, output=single_output, summary=single_summary))
        score_vcf_parallel(
            Args(
                **common,
                output=parallel_output,
                summary=parallel_summary,
                jobs=2,
                records_per_shard=1,
                temp_dir=self.work,
                keep_shards=False,
            )
        )

        def data_lines(path: str) -> List[str]:
            with gzip.open(path, "rt") as handle:
                return [line.rstrip("\n") for line in handle if not line.startswith("#")]

        self.assertEqual(data_lines(single_output), data_lines(parallel_output))
        with open(parallel_summary) as handle:
            summary_data = json.load(handle)
        self.assertTrue(summary_data["parallel"]["enabled"])
        self.assertEqual(summary_data["parallel"]["shards"], 3)
        self.assertEqual(summary_data["records_scored"], 3)
        self.assertEqual(summary_data["alt_alleles_scored"], 3)
        self.assertEqual(summary_data["validation"]["reference"]["mismatch_records"], 0)

    def test_parallel_fast_mode_region_merge_writes_single_header(self) -> None:
        vcf = os.path.join(self.work, "parallel_region.vcf")
        single_output = os.path.join(self.work, "parallel_region.single.vcf.gz")
        parallel_output = os.path.join(self.work, "parallel_region.sharded.vcf.gz")
        single_summary = os.path.join(self.work, "parallel_region.single.summary.json")
        parallel_summary = os.path.join(self.work, "parallel_region.sharded.summary.json")
        samples = [f"s{i}" for i in range(10)]
        with open(vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            out.write("chr4\t4\toutside\tT\tC\t100\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 10) + "\n")
            out.write("chr1\t14\tinside\tC\tG\t100\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 10) + "\n")

        common = dict(
            vcf=vcf,
            reference=self.reference,
            gff=self.gff,
            cds=self.cds,
            pep=self.pep,
            protein_domains=None,
            protein_structures=None,
            protein_lm_scores=None,
            mirna_sites=None,
            ml_model=None,
            sample_info=None,
            window=10,
            pop_window=0,
            skip_gene_constraint=True,
            config=None,
            region=["chr1:10-14"],
            html_report=None,
            index_output="never",
            limit=0,
            progress=0,
            strict_ref=False,
        )
        score_vcf(Args(**common, output=single_output, summary=single_summary))
        score_vcf_parallel(
            Args(
                **common,
                output=parallel_output,
                summary=parallel_summary,
                jobs=2,
                records_per_shard=1,
                temp_dir=self.work,
                keep_shards=False,
            )
        )

        def lines(path: str) -> Tuple[List[str], List[str]]:
            with gzip.open(path, "rt") as handle:
                header = []
                records = []
                for line in handle:
                    target = header if line.startswith("#") else records
                    target.append(line.rstrip("\n"))
                return header, records

        single_header, single_records = lines(single_output)
        parallel_header, parallel_records = lines(parallel_output)
        self.assertEqual(single_records, parallel_records)
        self.assertEqual(len(single_records), 1)
        self.assertEqual(sum(1 for line in parallel_header if line.startswith("#CHROM")), 1)
        self.assertEqual(sum(1 for line in parallel_header if line.startswith("##fileformat")), 1)
        with open(parallel_summary) as handle:
            summary_data = json.load(handle)
        self.assertEqual(summary_data["records_scored"], 1)
        self.assertEqual(summary_data["records_skipped_by_region"], 1)

    def test_region_parser(self) -> None:
        region = parse_region("chr1:10-20")
        self.assertTrue(region.contains("chr1", 10))
        self.assertTrue(region.contains("chr1", 20))
        self.assertFalse(region.contains("chr1", 21))
        self.assertFalse(region.contains("chr2", 15))


if __name__ == "__main__":
    unittest.main()
