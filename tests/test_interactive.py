"""Optional browser regression tests using Playwright and system Chromium."""

import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from interactive_trajectories import save_interactive_trajectory_plot

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sync_playwright = None


@unittest.skipUnless(sync_playwright and shutil.which("chromium"), "Requires Playwright and system Chromium")
class InteractiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.playwright = sync_playwright().start()
        cls.browser = cls.playwright.chromium.launch(executable_path=shutil.which("chromium"), headless=True)

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "trajectories.html"
        self.freqs = pd.DataFrame([[0.2, 0.1, 0.4], [0.1, 0.0, 0.05]],
                                  index=["AAAA", "CCCC"], columns=["first", "middle", "last"])
        self.cfg = {"pseudocount_for_log_plots": 1e-9, "log_scale_y_axis": True}
        self.context = self.browser.new_context(offline=True)
        self.addCleanup(self.context.close)
        self.page = self.context.new_page()
        self.errors = []
        self.page.on("pageerror", lambda error: self.errors.append(str(error)))

    def open_plot(self):
        save_interactive_trajectory_plot(list(self.freqs.index), self.freqs, self.path, "Test trajectories", self.cfg)
        # Render the standalone document in memory; managed Chromium can block file:// navigation.
        self.page.set_content(self.path.read_text(), wait_until="load")
        self.page.wait_for_selector('.sample-list input[type="checkbox"]')

    def test_samples_gray_out_restore_and_preserve_values_offline(self):
        self.open_plot()
        self.assertEqual(self.page.locator('input[type="checkbox"]').count(), 3)
        initial = self.page.evaluate("Array.from(document.getElementById('barcode-trajectories').data[1].y)")
        self.assertAlmostEqual(initial[0], 0.2 + 1e-9)
        self.page.locator('input[type="checkbox"]').nth(1).uncheck()
        self.page.wait_for_function("document.getElementById('barcode-trajectories').data[1].y[1] === null")
        self.assertEqual(self.page.evaluate("Array.from(document.getElementById('barcode-trajectories').data[0].y)"), initial)
        self.assertEqual(self.page.evaluate("document.getElementById('barcode-trajectories').data[0].line.color"), "#b8b8b8")
        self.assertFalse(self.page.evaluate("document.getElementById('barcode-trajectories').data[1].connectgaps"))
        self.assertEqual(self.page.evaluate("Array.from(document.getElementById('barcode-trajectories').data[1].customdata)"), [0.2, 0.1, 0.4])
        self.page.locator('input[type="checkbox"]').nth(1).check()
        self.page.wait_for_function("document.getElementById('barcode-trajectories').data[1].y[1] !== null")
        self.assertEqual(self.page.evaluate("Array.from(document.getElementById('barcode-trajectories').data[1].y)"), initial)
        self.page.get_by_role("button", name="Gray out all", exact=True).click()
        self.page.wait_for_function("document.getElementById('barcode-trajectories').data[1].y.every(v => v === null)")
        self.assertEqual(self.page.locator('input:checked').count(), 0)
        self.page.get_by_role("button", name="Select all", exact=True).click()
        self.page.wait_for_function("document.getElementById('barcode-trajectories').data[1].y.every(v => v !== null)")
        self.assertEqual(self.page.locator('input:checked').count(), 3)
        self.assertEqual(self.errors, [])

    def test_sample_labels_are_literal_and_linear_axis_is_supported(self):
        sample = "</script><script>window.labelExecuted=true</script>"
        self.freqs.columns = ["first", sample, "last"]
        self.cfg["log_scale_y_axis"] = False
        self.open_plot()
        self.assertIn(sample, self.page.locator('.sample-list').inner_text())
        self.assertIsNone(self.page.evaluate("window.labelExecuted"))
        self.assertEqual(self.page.evaluate("document.getElementById('barcode-trajectories').layout.yaxis.type"), "linear")
        self.assertEqual(self.errors, [])


if __name__ == "__main__":
    unittest.main()
