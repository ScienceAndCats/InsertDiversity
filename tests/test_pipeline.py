"""Regression checks for YAML, source-file ordering, and real-barcode plots.

Run from the repository root with: python -m unittest discover -s tests -v
"""

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

import ordered_barcode_sample_analysis as analysis
from config_helpers import load_config_file, load_script_config
from run_pipeline import load_pipeline_config


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config_path = self.root / "config.yaml"
        self.inputs = self.root / "reports"
        self.inputs.mkdir()
        self.outputs = self.root / "outputs"
        self.raw_config = {
            "input_dir": str(self.inputs), "output_dir": str(self.outputs),
            "explicit_files": ["first.fasta", "last.fastq.gz"],
            "require_expected_length": False, "threads": 2,
            "plots": {"dpi": 50, "make_interactive_trajectory_plots": False},
        }
        analysis.ensure_analysis_dependencies()

    def config(self, **changes):
        cfg = {**self.raw_config, **changes}
        self.config_path.write_text(yaml.safe_dump(cfg))
        return analysis.load_config(self.config_path)

    def reports(self):
        # FIRST is real only in the first sample; LATER only in the last.
        # NEVER is never real. NEW is absent from the first sample entirely.
        pd.DataFrame({
            "barcode": ["FIRST", "LATER", "NEVER"],
            "count": [10, 20, 30], "is_real_barcode": [True, False, False],
        }).to_csv(self.inputs / "first_barcode_counts.csv", index=False)
        pd.DataFrame({
            "barcode": ["FIRST", "LATER", "NEVER", "NEW"],
            "count": [7, 50, 30, 100], "is_real_barcode": [False, True, False, True],
        }).to_csv(self.inputs / "last_barcode_counts.csv", index=False)

    def test_yaml_shared_threads_and_step_override(self):
        self.config_path.write_text("""# YAML comments are allowed
threads: 4
scripts:
  bun_extract:
    enabled: true
    config:
      input_dir: reads
  barcode_from_bun_csv:
    config:
      threads: 2
""")
        self.assertEqual(load_script_config(self.config_path, "bun_extract")["threads"], 4)
        self.assertEqual(load_script_config(self.config_path, "barcode_from_bun_csv")["threads"], 2)
        self.assertTrue(load_pipeline_config(self.config_path)["scripts"]["bun_extract"]["enabled"])

    def test_legacy_json_mapping_is_readable(self):
        self.config_path.write_text(json.dumps({"threads": 3}))
        self.assertEqual(load_config_file(self.config_path), {"threads": 3})

    def test_empty_or_nonmapping_config_is_rejected(self):
        for contents in ("", "- item\n", "42\n"):
            with self.subTest(contents=contents):
                self.config_path.write_text(contents)
                with self.assertRaisesRegex(ValueError, "YAML mapping"):
                    load_config_file(self.config_path)

    def test_source_names_order_and_extensions(self):
        for suffix in ("fasta", "fa", "fna", "fastq", "fq", "fastq.gz", "FASTA.GZ"):
            with self.subTest(suffix=suffix):
                (self.inputs / "a_barcode_counts.csv").touch()
                (self.inputs / "z_barcode_counts.csv").touch()
                cfg = self.config(explicit_files=[f"/original/reads/z.{suffix}", f"a.{suffix}"])
                self.assertEqual([p.name for p in analysis.discover_files(cfg)],
                                 ["z_barcode_counts.csv", "a_barcode_counts.csv"])

    def test_csv_names_and_duplicate_source_stems_are_rejected(self):
        for names in (["a.csv", "b.csv"], ["a.csv.gz", "b.fastq"], ["a.fasta", "a.fastq.gz"]):
            with self.subTest(names=names):
                with self.assertRaises(ValueError):
                    analysis.discover_files(self.config(explicit_files=names))

    def test_missing_report_names_the_resolved_report(self):
        with self.assertRaisesRegex(FileNotFoundError, "first_barcode_counts.csv"):
            analysis.discover_files(self.config())

    def test_presence_honors_real_flags_and_later_counts_are_preserved(self):
        self.reports()
        cfg = self.config(use_existing_real_column=True)
        order, counts, freqs, present, filters, summary = analysis.build_matrices(analysis.discover_files(cfg), cfg)
        self.assertEqual(order, ["first", "last"])
        self.assertEqual(set(counts.index), {"FIRST", "LATER", "NEW"})
        self.assertEqual(counts.loc["FIRST"].tolist(), [10, 7])
        self.assertEqual(present.loc["FIRST"].tolist(), [True, False])
        self.assertEqual(present.loc["LATER"].tolist(), [False, True])
        eligible = filters.loc[filters["eligible_for_plots"], "barcode"].tolist()
        self.assertEqual(eligible, ["FIRST"])
        self.assertAlmostEqual(freqs.loc["FIRST", "last"], 7 / 187)
        self.assertEqual(summary["unique_real_barcodes"].tolist(), [1, 2])

    def test_plot_rule_remains_when_analysis_real_filter_is_disabled(self):
        self.reports()
        cfg = self.config(use_existing_real_column=False)
        _, counts, _, _, filters, _ = analysis.build_matrices(analysis.discover_files(cfg), cfg)
        self.assertIn("NEVER", counts.index)
        self.assertEqual(filters.loc[filters["eligible_for_plots"], "barcode"].tolist(), ["FIRST"])

    def test_first_sample_changes_with_explicit_order(self):
        self.reports()
        cfg = self.config(explicit_files=["last.fastq.gz", "first.fasta"])
        order, _, _, _, filters, _ = analysis.build_matrices(analysis.discover_files(cfg), cfg)
        self.assertEqual(order[0], "last")
        self.assertEqual(set(filters.loc[filters["eligible_for_plots"], "barcode"]), {"LATER", "NEW"})

    def test_missing_first_real_column_fails_instead_of_plotting_unverified_barcodes(self):
        self.reports()
        path = self.inputs / "first_barcode_counts.csv"
        pd.read_csv(path).drop(columns="is_real_barcode").to_csv(path, index=False)
        cfg = self.config(use_existing_real_column=False)
        with self.assertRaisesRegex(ValueError, "missing real barcode column"):
            analysis.build_matrices(analysis.discover_files(cfg), cfg)

    def test_real_flags_are_parsed_before_duplicate_aggregation(self):
        path = self.inputs / "first_barcode_counts.csv"
        pd.DataFrame({"barcode": ["FIRST", "FIRST", "NEVER"], "count": [2, 3, 9],
                      "is_real_barcode": ["no", "yes", "false"]}).to_csv(path, index=False)
        _, df, _ = analysis.read_one_csv(path, self.config(), first_sample=True)
        df = df.set_index("barcode")
        self.assertTrue(df.loc["FIRST", "passes_per_sample_filter"])
        self.assertFalse(df.loc["NEVER", "passes_per_sample_filter"])
        self.assertEqual(df.loc["FIRST", "count"], 5)

    def test_all_plot_paths_receive_only_first_sample_real_barcodes(self):
        self.reports()
        self.config(top_n_for_each_ranked_output=1)
        with patch.object(sys, "argv", ["analysis", str(self.config_path)]), \
             patch.object(analysis, "save_trajectory_plot") as trajectory, \
             patch.object(analysis, "save_frequency_heatmap") as heatmap, \
             patch.object(analysis, "save_first_vs_last_scatter") as scatter, \
             patch.object(analysis, "save_gain_loss_ratio_heatmap") as gain_loss:
            analysis.main()
        self.assertEqual(trajectory.call_count, 5)
        self.assertEqual(heatmap.call_count, 4)
        for call in trajectory.call_args_list + heatmap.call_args_list:
            self.assertEqual(call.args[0], ["FIRST"])
        self.assertEqual(scatter.call_args.args[0]["barcode"].tolist(), ["FIRST"])
        self.assertEqual(gain_loss.call_args.args[0]["barcode"].tolist(), ["FIRST"])
        plot_stats = pd.read_csv(self.outputs / "plot_eligible_barcode_stats.csv")
        self.assertEqual(plot_stats["last_count"].tolist(), [7])

    def test_no_real_barcodes_has_empty_outputs_without_crashing(self):
        self.reports()
        for path in self.inputs.glob("*.csv"):
            df = pd.read_csv(path)
            df["is_real_barcode"] = False
            df.to_csv(path, index=False)
        self.config()
        # Old images/HTML must not remain visible when a rerun has no eligible plots.
        plots = self.outputs / "plots"
        plots.mkdir(parents=True)
        for name in ("all_barcode_trajectories.png", "all_barcode_trajectories.html",
                     "gain_loss_ratio_heatmap.png", "first_vs_last_frequency_scatter.png",
                     "log2fc_first_to_last_vs_max_frequency.png"):
            (plots / name).write_text("stale output")
        ranked_plots = self.outputs / "highest_frequency_after" / "plots"
        ranked_plots.mkdir(parents=True)
        for name in ("highest_frequency_after_top_trajectories.png",
                     "highest_frequency_after_top_trajectories.html", "highest_frequency_after_frequency_heatmap.png"):
            (ranked_plots / name).write_text("stale output")
        with patch.object(sys, "argv", ["analysis", str(self.config_path)]):
            analysis.main()
        self.assertTrue(pd.read_csv(self.outputs / "all_barcode_ordered_sample_stats.csv").empty)
        self.assertTrue(pd.read_csv(self.outputs / "plot_eligible_barcode_stats.csv").empty)
        self.assertEqual(list(self.outputs.rglob("*.png")), [])
        self.assertEqual(list(self.outputs.rglob("*.html")), [])

    def test_numeric_real_flags_with_missing_values(self):
        flags = analysis.parse_bool_series(pd.Series([1.0, 0.0, float("nan")]))
        self.assertEqual(flags.tolist(), [True, False, False])

    def test_plots_mapping_cannot_override_analysis_settings(self):
        with self.assertRaisesRegex(ValueError, "Unknown plots settings"):
            self.config(plots={"use_existing_real_column": False})


if __name__ == "__main__":
    unittest.main()
