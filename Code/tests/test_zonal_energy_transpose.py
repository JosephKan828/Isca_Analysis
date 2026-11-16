"""Check the frequency attribution, transport closure, and hydrostatic ordering."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "zonal_energy_transpose.py"
SPEC = importlib.util.spec_from_file_location("zonal_energy_transpose", SCRIPT)
zonal = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(zonal)


class FullRecordTransportTests(unittest.TestCase):
    def setUp(self):
        self.nt = 360
        t = np.arange(self.nt)
        self.u = np.full(self.nt, 4.0)
        self.s = np.full(self.nt, 300000.0)
        for period, amplitude_u, amplitude_s, phase in (
            (120, 2.0, 1000.0, 0.0),
            (20, 1.0, 600.0, np.pi),
            (6, 0.5, 200.0, np.pi / 3),
        ):
            self.u += amplitude_u * np.cos(2 * np.pi * t / period)
            self.s += amplitude_s * np.cos(2 * np.pi * t / period + phase)
        self.u = self.u[:, None, None, None]
        self.s = self.s[:, None, None, None]
        _, self.masks = zonal.frequency_bands(self.nt, 1.0, 30.0, 10.0)

    def test_signed_band_transport_and_all_sixteen_term_closure(self):
        result = zonal.transport_budget(self.u, self.s, self.masks)
        terms = result["transport_terms"][:, :, 0, 0, 0]
        expected = np.diag([1200000.0, 1000.0, -300.0, 25.0])
        np.testing.assert_allclose(terms, expected, rtol=1e-11, atol=1e-8)
        np.testing.assert_allclose(
            terms.sum(), np.mean(self.u * self.s), rtol=1e-12,
        )
        self.assertLess(result["reconstruction_errors"][0], 1e-12)
        self.assertLess(result["reconstruction_errors"][1], 1e-8)

    def test_odd_length_and_nyquist_reconstruction(self):
        for nt in (360, 367):
            t = np.arange(nt)
            field = 7 + np.cos(2 * np.pi * 3 * t / nt)
            field += np.cos(2 * np.pi * 20 * t / nt)
            field += np.cos(2 * np.pi * (nt // 2) * t / nt)
            _, masks = zonal.frequency_bands(nt, 1.0, 30.0, 10.0)
            mean, bands = zonal.split_time_bands(field, masks)
            self.assertFalse(any(mask[0] for mask in masks))
            np.testing.assert_allclose(mean + sum(bands), field, atol=1e-12)

    def test_height_is_referenced_to_bottom_in_either_pressure_order(self):
        for pressure in (np.array([250., 500., 1000.]), np.array([1000., 500., 250.])):
            temperature = np.full((2, 3, 1, 2), 300.)
            actual = zonal.geopotential_height(pressure, temperature)
            expected = zonal.RD / zonal.GRAVITY * 300 * np.log(1000 / pressure)
            np.testing.assert_allclose(actual[0, :, 0, 0], expected)
            self.assertEqual(actual[0, pressure.argmax(), 0, 0], 0)

    def test_meridional_widths_include_partial_cells_and_preserve_order(self):
        latitude = np.array([-30., -10., 10., 30.])
        expected = zonal.EARTH_RADIUS * np.deg2rad([5., 20., 20., 5.])
        widths = zonal.meridional_strip_widths(latitude, 25.)
        np.testing.assert_allclose(widths, expected)
        permutation = np.array([2, 0, 3, 1])
        np.testing.assert_allclose(
            zonal.meridional_strip_widths(latitude[permutation], 25.), expected[permutation],
        )

    def test_meridional_widths_cover_the_requested_arc_without_normalization(self):
        latitude = np.array([-85., -40., -5., 25., 70.])
        for limit in (20., 60., 90.):
            widths = zonal.meridional_strip_widths(latitude, limit)
            self.assertAlmostEqual(
                widths.sum() / zonal.EARTH_RADIUS, np.deg2rad(2 * limit),
            )
            self.assertTrue(np.all(widths >= 0))
        np.testing.assert_allclose(
            zonal.meridional_strip_widths(latitude, 20., radius=2 * zonal.EARTH_RADIUS),
            2 * zonal.meridional_strip_widths(latitude, 20.),
        )


if __name__ == "__main__":
    unittest.main()
