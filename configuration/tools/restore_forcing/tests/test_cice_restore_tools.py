from __future__ import annotations

from datetime import date, datetime
from pathlib import Path
import sys
import tempfile
import types
import unittest

import netCDF4
import numpy as np


TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
sys.modules.setdefault("xesmf", types.ModuleType("xesmf"))

import convert_cice_history_to_restore_forcing as history_converter
import convert_climate_ice_to_restore_forcing as climate_converter
from validate_cice_restore_forcing import RuntimeConfig, validate_restore_file


CONFIG = RuntimeConfig(nx=2, ny=2, ncat=2, nilyr=2, ntrcr=6, fyear_init=2017)


def make_forcing(path: Path, schema: str, nrecords: int = 12, data_type: str = "monthly_netcdf") -> None:
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("time", nrecords)
        ds.createDimension("nc", CONFIG.ncat)
        ds.createDimension("nj", CONFIG.ny)
        ds.createDimension("ni", CONFIG.nx)
        ds.createDimension("nilyr", CONFIG.nilyr)
        ds.createDimension("ntrcr", CONFIG.ntrcr)
        time = ds.createVariable("time", "f8", ("time",))
        time[:] = np.arange(nrecords)
        prefix = "days" if data_type == "daily_netcdf" else "months"
        time.units = f"{prefix} since 2017-01-01 00:00:00"
        for name in ("aicen", "vicen", "vsnon"):
            var = ds.createVariable(name, "f4", ("time", "nc", "nj", "ni"))
            var[:] = 0.0
        ds.variables["aicen"][:, 0, :, :] = 0.5
        ds.variables["vicen"][:, 0, :, :] = 0.75
        ds.variables["vsnon"][:, 0, :, :] = 0.05
        if schema == "full_trcrn":
            ds.createVariable("trcrn", "f4", ("time", "nc", "ntrcr", "nj", "ni"))[:] = 0.0
        elif schema == "bc_fields_v1":
            ds.createVariable("Tsfc", "f4", ("time", "nc", "nj", "ni"))[:] = -5.0
            ds.createVariable("Tinz", "f4", ("time", "nc", "nilyr", "nj", "ni"))[:] = -4.0
            ds.createVariable("Sinz", "f4", ("time", "nc", "nilyr", "nj", "ni"))[:] = 4.0
        ds.restore_format = schema


class ValidatorTests(unittest.TestCase):
    def test_accepts_all_supported_schemas(self):
        with tempfile.TemporaryDirectory() as directory:
            for schema in ("full_trcrn", "bc_fields_v1", "structure_only"):
                with self.subTest(schema=schema):
                    path = Path(directory) / f"{schema}.nc"
                    make_forcing(path, schema)
                    self.assertEqual(validate_restore_file(path, CONFIG, "monthly_netcdf", True), schema)

    def test_rejects_partial_thermodynamics(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "partial.nc"
            make_forcing(path, "structure_only")
            with netCDF4.Dataset(path, "a") as ds:
                ds.createVariable("Tsfc", "f4", ("time", "nc", "nj", "ni"))[:] = -5.0
            with self.assertRaisesRegex(ValueError, "incomplete thermodynamic schema"):
                validate_restore_file(path, CONFIG, "monthly_netcdf", True)

    def test_rejects_wrong_cyclic_record_count(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "short.nc"
            make_forcing(path, "structure_only", nrecords=11)
            with self.assertRaisesRegex(ValueError, "exactly 12"):
                validate_restore_file(path, CONFIG, "monthly_netcdf", True)

    def test_rejects_invalid_values_and_dimensions(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.nc"
            make_forcing(path, "structure_only")
            with netCDF4.Dataset(path, "a") as ds:
                ds.variables["aicen"][0, :, 0, 0] = (0.8, 0.7)
            with self.assertRaisesRegex(ValueError, "summed aicen"):
                validate_restore_file(path, CONFIG, "monthly_netcdf", True)
            with self.assertRaisesRegex(ValueError, "shape"):
                validate_restore_file(
                    path,
                    RuntimeConfig(nx=3, ny=2, ncat=2, nilyr=2, ntrcr=6, fyear_init=2017),
                    "monthly_netcdf",
                    True,
                )

    def test_accepts_small_positive_category_area_but_rejects_zero_area_with_mass(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "thin_tail.nc"
            make_forcing(path, "structure_only")
            with netCDF4.Dataset(path, "a") as ds:
                ds.variables["aicen"][0, 1, 0, 0] = 5.0e-7
                ds.variables["vicen"][0, 1, 0, 0] = 5.0e-6
            self.assertEqual(validate_restore_file(path, CONFIG, "monthly_netcdf", True), "structure_only")
            with netCDF4.Dataset(path, "a") as ds:
                ds.variables["aicen"][0, 1, 0, 0] = 0.0
            with self.assertRaisesRegex(ValueError, "where category concentration is zero"):
                validate_restore_file(path, CONFIG, "monthly_netcdf", True)

    def test_noncyclic_coverage_and_time_anchor(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "daily.nc"
            make_forcing(path, "structure_only", nrecords=3, data_type="daily_netcdf")
            self.assertEqual(
                validate_restore_file(
                    path,
                    CONFIG,
                    "daily_netcdf",
                    False,
                    datetime(2017, 1, 1),
                    datetime(2017, 1, 3),
                ),
                "structure_only",
            )
            with self.assertRaisesRegex(ValueError, "outside forcing interpolation coverage"):
                validate_restore_file(path, CONFIG, "daily_netcdf", False, simulation_end=datetime(2017, 1, 3, 12))
            with netCDF4.Dataset(path, "a") as ds:
                ds.variables["time"].units = "days since 2017-02-01 00:00:00"
            with self.assertRaisesRegex(ValueError, "time units"):
                validate_restore_file(path, CONFIG, "daily_netcdf", False)

    def test_legacy_one_based_time_requires_explicit_opt_in(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy_time.nc"
            make_forcing(path, "structure_only")
            with netCDF4.Dataset(path, "a") as ds:
                ds.variables["time"][:] = np.arange(1, 13)
            with self.assertRaisesRegex(ValueError, "zero-based"):
                validate_restore_file(path, CONFIG, "monthly_netcdf", True)
            with self.assertWarnsRegex(UserWarning, "legacy one-based"):
                self.assertEqual(
                    validate_restore_file(
                        path,
                        CONFIG,
                        "monthly_netcdf",
                        True,
                        allow_one_based_time=True,
                    ),
                    "structure_only",
                )


class ConverterTests(unittest.TestCase):
    def test_history_writer_uses_zero_based_index_time(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "history_restore.nc"
            shape = (3, CONFIG.ncat, CONFIG.ny, CONFIG.nx)
            aicen = np.zeros(shape, dtype=np.float32)
            vicen = np.zeros(shape, dtype=np.float32)
            vsnon = np.zeros(shape, dtype=np.float32)
            aicen[:, 0] = 0.5
            vicen[:, 0] = 0.75
            vsnon[:, 0] = 0.05
            history_converter.write_output(
                output_path=path,
                source_paths=[Path(f"iceh.2017-01-0{day}.nc") for day in range(1, 4)],
                day_offsets=np.arange(3, dtype=np.float64),
                aicen_all=aicen,
                vicen_all=vicen,
                vsnon_all=vsnon,
                trcrn_all=None,
                tsfc_all=None,
                tinz_all=None,
                sinz_all=None,
                compression_level=0,
                source_kind="history",
                restore_format="structure_only",
            )
            with netCDF4.Dataset(path) as ds:
                np.testing.assert_array_equal(ds.variables["time"][:], np.arange(3))
                self.assertEqual(ds.restore_format, "structure_only")
            self.assertEqual(validate_restore_file(path, CONFIG, "daily_netcdf", False), "structure_only")

    def test_climate_writer_declares_structure_only(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "climate_restore.nc"
            shape = (2, 5, 1, 1)
            fields = np.zeros(shape, dtype=np.float32)
            climate_converter.write_output(
                output_path=path,
                requested_months=[date(2017, 1, 1), date(2017, 2, 1)],
                aicen_all=fields,
                vicen_all=fields,
                vsnon_all=fields,
                compression_level=0,
                aice_file=Path("siconc.nc"),
                vice_file=Path("sivol.nc"),
                snow_source_desc="sisnthick.nc",
                regrid_method="bilinear",
                diagnostics={},
            )
            with netCDF4.Dataset(path) as ds:
                np.testing.assert_array_equal(ds.variables["time"][:], (0.0, 1.0))
                self.assertEqual(ds.restore_format, "structure_only")

    def test_target_grid_is_derived_from_ice_in(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            grid = root / "grids" / "grid.nc"
            grid.parent.mkdir()
            grid.touch()
            ice_in = root / "ice_in"
            ice_in.write_text(f"&grid_nml\n grid_file = '{grid}'\n/\n", encoding="utf-8")
            self.assertEqual(history_converter.resolve_target_grid_file(None, ice_in), grid)
            self.assertEqual(climate_converter.resolve_target_grid_file(None, ice_in), grid)

    def test_climate_units_and_category_assignment(self):
        np.testing.assert_allclose(climate_converter.normalize_aice(np.array([50.0]), "%"), [0.5])
        np.testing.assert_allclose(climate_converter.normalize_length(np.array([125.0]), "cm", "ice"), [1.25])
        aice = np.array([[[0.5]]], dtype=np.float32)
        vice = np.array([[[0.75]]], dtype=np.float32)
        vsno = np.array([[[0.02]]], dtype=np.float32)
        aicen, vicen, vsnon, _ = climate_converter.apply_category_template(aice, vice, vsno)
        self.assertEqual(np.flatnonzero(aicen[0, :, 0, 0]).tolist(), [2])
        self.assertAlmostEqual(float(vicen.sum()), 0.75)
        self.assertAlmostEqual(float(vsnon.sum()), 0.02)

    def test_bc_fields_mask_open_water(self):
        config = history_converter.IceInConfig(
            ncat=2,
            nilyr=2,
            nslyr=1,
            ktherm=2,
            tr_iage=False,
            tr_fy=False,
            tr_lvl=False,
            tr_pond_lvl=False,
            tr_pond_topo=False,
            tr_pond_sealvl=False,
            tr_snow=False,
            tr_iso=False,
            tr_aero=False,
            tr_fsd=False,
            tr_brine=False,
        )
        shape = (1, 2, 1, 1)
        fields = {
            "aicen": np.array([[[[0.5]], [[0.0]]]], dtype=np.float32),
            "vsnon": np.zeros(shape, dtype=np.float32),
            "tinz": np.full((1, 2, 2, 1, 1), -3.0, dtype=np.float32),
            "sinz": np.full((1, 2, 2, 1, 1), 4.0, dtype=np.float32),
            "tsnz": None,
            "tsfc": None,
        }
        tsfc, tinz, sinz = history_converter.build_bc_fields_from_history(fields, config)
        self.assertEqual(float(tsfc[0, 1, 0, 0]), 0.0)
        self.assertEqual(float(tinz[0, 1, 0, 0, 0]), 0.0)
        self.assertEqual(float(sinz[0, 0, 0, 0, 0]), 4.0)


if __name__ == "__main__":
    unittest.main()
