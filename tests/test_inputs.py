"""Input-boundary tests; no biological inference is tested or implied."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from analysis.validate_inputs import GROUPS, validate_inputs


class InputPreflightTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.path = self.root / "results.csv"
        self.header = ["gene_id", "log2FoldChange", "stat", "pvalue", "padj"]
        self.rows = [["Gnai1", "0.3", "2", "0.02", "0.04"],
                     ["Gnai2", "-0.4", "-3", "0.005", "0.02"]]
        self.write(self.path, self.header, self.rows)

    @staticmethod
    def write(path, header, rows):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    def check(self, **overrides):
        options = dict(kind="deg", gene_id_type="symbol", numerator="LD SALINE",
                       denominator="NON-LD control", complete_tested_universe=True)
        options.update(overrides)
        return validate_inputs(self.path, **options)

    @staticmethod
    def codes(report, level="errors"):
        return {item["code"] for item in report[level]}

    def metadata(self, samples_per_group=2):
        path = self.root / "samples.csv"
        rows = [[f"sample_{i}_{j}", group] for i, group in enumerate(GROUPS)
                for j in range(samples_per_group)]
        self.write(path, ["sample_id", "group"], rows)
        return path, [row[0] for row in rows]

    def test_valid_deg_is_read_only_and_hashes_parsed_bytes(self):
        before = self.path.read_bytes()
        report = self.check()
        self.assertTrue(report["valid"])
        self.assertTrue(report["readiness"]["gsea_rank_input"])
        self.assertEqual(report["summary"]["data_rows"], 2)
        self.assertEqual(report["files"]["data"]["sha256"], hashlib.sha256(before).hexdigest())
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(sorted(p.name for p in self.root.iterdir()), ["results.csv"])

    def test_optional_raw_pvalue_column(self):
        self.write(self.path, self.header[:3] + ["padj"], [row[:3] + [row[4]] for row in self.rows])
        self.assertTrue(self.check()["valid"])

    def test_missing_padj_is_reported_but_not_used_to_filter_gsea(self):
        self.rows[0][-1] = "NA"
        self.write(self.path, self.header, self.rows)
        report = self.check()
        self.assertTrue(report["valid"])
        self.assertTrue(report["readiness"]["gsea_rank_input"])
        self.assertEqual(report["summary"]["rows_with_finite_signed_stat"], 2)
        self.assertEqual(report["summary"]["rows_with_effect_and_padj"], 1)
        self.assertEqual(report["exclusions"]["missing_padj"]["gene_ids"], ["Gnai1"])

    def test_missing_stat_requires_documented_reason(self):
        self.rows[0][2] = "NA"
        self.write(self.path, self.header, self.rows)
        self.assertIn("missing_stat", self.codes(self.check()))
        report = self.check(missing_stat_reason="Upstream all-zero untestable rows.")
        self.assertTrue(report["valid"])
        self.assertEqual(report["exclusions"]["missing_stat"]["gene_ids"], ["Gnai1"])
        self.assertFalse(report["readiness"]["gsea_rank_input"])

    def test_documented_untestable_rows_can_retain_na_de_fields(self):
        self.rows.append(["Gnai3", "NA", "NA", "NA", "NA"])
        self.write(self.path, self.header, self.rows)
        report = self.check(missing_stat_reason="All-zero genes were untestable upstream.")
        self.assertTrue(report["valid"])
        self.assertEqual(report["summary"]["data_rows"], 3)
        self.assertEqual(report["exclusions"]["missing_effect"]["gene_ids"], ["Gnai3"])
        self.assertTrue(report["readiness"]["gsea_rank_input"])

    def test_sign_mismatch_is_an_error_and_never_flipped(self):
        self.rows[1][2] = "3"
        self.write(self.path, self.header, self.rows)
        before = self.path.read_bytes()
        self.assertIn("statistic_direction", self.codes(self.check()))
        self.assertEqual(before, self.path.read_bytes())

    def test_nonfinite_and_out_of_range_values_fail(self):
        for column, value, code in [(1, "NaN", "non_finite"), (2, "Infinity", "non_finite"),
                                    (3, "1.01", "probability_range"), (4, "-0.1", "numeric_range"),
                                    (1, " 1", "invalid_numeric"), (2, "1_000", "invalid_numeric"),
                                    (1, "1e10000", "float_range"), (2, "-1e10000", "float_range"),
                                    (4, "1e-10000", "numeric_underflow")]:
            with self.subTest(column=column, value=value):
                rows = [row[:] for row in self.rows]
                rows[0][column] = value
                self.write(self.path, self.header, rows)
                self.assertIn(code, self.codes(self.check()))

    def test_duplicate_genes_headers_and_malformed_rows_fail(self):
        self.write(self.path, self.header, [self.rows[0], self.rows[0]])
        self.assertIn("duplicate_gene_id", self.codes(self.check()))
        self.write(self.path, self.header + ["padj"], [self.rows[0] + ["0.04"]])
        self.assertIn("duplicate_header", self.codes(self.check()))
        self.write(self.path, self.header, [self.rows[0][:-1]])
        self.assertIn("row_width", self.codes(self.check()))
        self.path.write_text('gene_id,log2FoldChange,stat,padj\n"unterminated,1,2,0.01\n')
        self.assertIn("malformed_csv", self.codes(self.check()))

    def test_empty_files_and_whitespace_identifiers_fail(self):
        self.path.write_text("")
        self.assertIn("empty_file", self.codes(self.check()))
        self.write(self.path, self.header, [])
        self.assertIn("no_data_rows", self.codes(self.check()))
        self.rows[0][0] = " Gnai1"
        self.write(self.path, self.header, self.rows)
        self.assertIn("invalid_gene_id", self.codes(self.check()))

    def test_enumerated_contrasts_and_id_declarations(self):
        self.assertIn("contrast_declaration", self.codes(self.check(numerator="S")))
        self.assertIn("contrast_declaration", self.codes(self.check(numerator="NON-LD control")))
        self.assertIn("gene_id_declaration", self.codes(self.check(gene_id_type=None)))
        self.rows[0][0] = "ENSMUSG00000000001.4"
        self.rows[1][0] = "ENSMUSG00000000002"
        self.write(self.path, self.header, self.rows)
        self.assertTrue(self.check(gene_id_type="ensembl_mouse")["valid"])
        self.assertIn(b"ENSMUSG00000000001.4", self.path.read_bytes())
        self.rows[0][0] = "ENSG00000000001"
        self.write(self.path, self.header, self.rows)
        self.assertIn("gene_id_type", self.codes(self.check(gene_id_type="ensembl_mouse")))

    def test_universe_is_not_inferred_from_gene_numbers_or_pvalues(self):
        report = self.check(complete_tested_universe=False)
        self.assertTrue(report["valid"])
        self.assertFalse(report["readiness"]["gsea_rank_input"])
        self.assertIn("universe_not_confirmed", self.codes(report, "warnings"))

    def test_inline_fpkm_uses_exact_sample_metadata_without_changing_de(self):
        metadata, samples = self.metadata()
        self.write(self.path, self.header + samples, [row + ["1.5"] * len(samples) for row in self.rows])
        before = self.path.read_bytes()
        report = self.check(metadata=metadata, sample_expression="fpkm")
        self.assertTrue(report["valid"])
        self.assertTrue(report["readiness"]["heatmap_expression_input"])
        self.assertEqual(report["summary"]["samples_per_group"], {g: 2 for g in GROUPS})
        self.assertEqual(self.path.read_bytes(), before)

    def test_inline_fpkm_missing_extra_or_invalid_sample_cells_fail(self):
        metadata, samples = self.metadata()
        self.write(self.path, self.header + samples[:-1], [row + ["1"] * (len(samples)-1) for row in self.rows])
        self.assertIn("sample_mismatch", self.codes(self.check(metadata=metadata, sample_expression="fpkm")))
        self.write(self.path, self.header + samples + ["unknown_sample"], [row + ["1"] * (len(samples)+1) for row in self.rows])
        self.assertIn("sample_mismatch", self.codes(self.check(metadata=metadata, sample_expression="fpkm")))
        self.write(self.path, self.header + samples, [row + ["-0.2"] * len(samples) for row in self.rows])
        self.assertIn("numeric_range", self.codes(self.check(metadata=metadata, sample_expression="fpkm")))
        self.assertIn("inline_expression_declaration", self.codes(self.check(sample_expression="fpkm")))

    def test_extra_gene_annotations_must_be_declared_for_inline_fpkm(self):
        metadata, samples = self.metadata()
        self.write(self.path, self.header + samples + ["description"], [row + ["1"] * len(samples) + ["gene annotation"] for row in self.rows])
        self.assertFalse(self.check(metadata=metadata, sample_expression="fpkm")["valid"])
        report = self.check(metadata=metadata, sample_expression="fpkm", annotation_columns=["description"])
        self.assertTrue(report["valid"])

    def test_inline_fpkm_accepts_two_compared_groups_but_requires_both(self):
        metadata = self.root / "samples.csv"
        rows = [["saline_1", "LD SALINE"], ["saline_2", "LD SALINE"],
                ["control_1", "NON-LD control"], ["control_2", "NON-LD control"]]
        self.write(metadata, ["sample_id", "group"], rows)
        samples = [row[0] for row in rows]
        self.write(self.path, self.header + samples, [row + ["1"] * 4 for row in self.rows])
        report = self.check(metadata=metadata, sample_expression="fpkm")
        self.assertTrue(report["valid"])
        self.assertEqual(report["summary"]["samples_per_group"]["LD LEPTIN"], 0)
        report = self.check(metadata=metadata, sample_expression="fpkm", numerator="LD LEPTIN")
        self.assertIn("insufficient_replicates", self.codes(report))

    def test_inline_sample_ids_cannot_reuse_a_de_statistic_column(self):
        metadata, samples = self.metadata()
        rows = [["stat" if sample == samples[0] else sample, GROUPS[i//2]] for i, sample in enumerate(samples)]
        self.write(metadata, ["sample_id", "group"], rows)
        self.write(self.path, self.header + samples[1:], [["Gnai1", "0.3", "2", "0.02", "0.04"] + ["1"] * 5])
        report = self.check(metadata=metadata, sample_expression="fpkm")
        self.assertIn("sample_column_collision", self.codes(report))

    def test_metadata_duplicate_or_repaired_group_labels_are_rejected(self):
        metadata, samples = self.metadata()
        self.write(self.path, ["gene_id"] + samples, [["Gnai1"] + ["1"] * len(samples)])
        metadata.write_text("sample_id,group\nsample_0_0,LD SALINE\nsample_0_0,LD SALINE\nsample_1_0,L\n")
        report = validate_inputs(self.path, kind="expression", gene_id_type="symbol", metadata=metadata)
        self.assertTrue({"duplicate_sample_id", "invalid_group", "insufficient_replicates"}.issubset(self.codes(report)))

    def test_raw_counts_are_nonnegative_integers_expression_can_be_fractional(self):
        metadata, samples = self.metadata()
        self.write(self.path, ["gene_id"] + samples, [["Gnai1"] + ["1.5"] * len(samples)])
        report = validate_inputs(self.path, kind="counts", gene_id_type="symbol", metadata=metadata)
        self.assertIn("non_integer_count", self.codes(report))
        report = validate_inputs(self.path, kind="expression", gene_id_type="symbol", metadata=metadata)
        self.assertTrue(report["valid"])
        self.assertTrue(report["readiness"]["heatmap_expression_input"])

    def test_utf8_bom_and_cli_json_exit_status(self):
        self.path.write_bytes(b"\xef\xbb\xbf" + self.path.read_bytes())
        self.assertTrue(self.check()["valid"])
        command = [sys.executable, str(Path(__file__).parents[1]/"analysis/validate_inputs.py"), str(self.path),
                   "--kind", "deg", "--gene-id-type", "symbol", "--numerator", "LD SALINE", "--denominator", "NON-LD control"]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0)
        self.assertFalse(json.loads(result.stdout)["readiness"]["gsea_rank_input"])
        self.path.write_text("gene_id,stat\nGnai1,1\n")
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 1)
        self.assertFalse(json.loads(result.stdout)["valid"])


if __name__ == "__main__":
    unittest.main()
