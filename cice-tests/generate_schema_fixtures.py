#!/usr/bin/env python3
"""Generate small PSKRIPS full-tracer and structure-only BC fixtures."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "configuration/tools/restore_forcing/convert_cice_history_to_restore_forcing.py"
ICE_IN = ROOT / "cice-tests/pskrips_bc_diag/ice_in"
SOURCE = Path(
    "/nfs/scratch/paulinan/cice-dirs/input/CICE_data/ic/pskrips/"
    "iced_pskrips_2017-02-01.onice_hi0p5_hs0p01.thin_tail_itd_minarea0p05.nc"
)
OUTPUT_DIR = ROOT / "cice-tests/fixtures"


def load_converter():
    spec = spec_from_file_location("restore_converter", CONVERTER)
    module = module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def stack_fields(base, count):
    return {
        name: None if value is None else np.repeat(value[None, ...], count, axis=0)
        for name, value in base.items()
    }


def scale_structure(fields, factors):
    scale = np.asarray(factors, dtype=np.float32)[:, None, None, None]
    for name in ("aicen", "vicen", "vsnon"):
        fields[name] = np.asarray(fields[name], dtype=np.float32) * scale


def main():
    conv = load_converter()
    config = conv.parse_ice_in(ICE_IN)
    layout = conv.build_tracer_layout(config)
    base = conv.read_restart_day(SOURCE, config)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    full_fields = stack_fields(base, 2)
    scale_structure(full_fields, [0.80, 0.60])
    full_trcrn = conv.build_trcrn_from_restart_fields(full_fields, config, layout)
    full_path = OUTPUT_DIR / "pskrips_full_trcrn_2day.nc"
    conv.write_output(
        output_path=full_path,
        source_paths=[Path("iced.2017-01-01-00000.nc"), Path("iced.2017-01-02-00000.nc")],
        day_offsets=np.asarray([0.0, 1.0]),
        aicen_all=full_fields["aicen"],
        vicen_all=full_fields["vicen"],
        vsnon_all=full_fields["vsnon"],
        trcrn_all=full_trcrn,
        tsfc_all=None,
        tinz_all=None,
        sinz_all=None,
        compression_level=2,
        source_kind="restart_fixture",
        restore_format="full_trcrn",
    )
    with conv.netCDF4.Dataset(full_path, "r+") as dataset:
        dataset.variables["time"].units = "days since 2017-01-01 00:00:00"

    structure_fields = stack_fields(base, 12)
    scale_structure(structure_fields, np.linspace(0.55, 0.77, 12))
    structure_path = OUTPUT_DIR / "pskrips_structure_only_12month.nc"
    conv.write_output(
        output_path=structure_path,
        source_paths=[Path("iced.2017-01-01-00000.nc"), Path("iced.2017-12-01-00000.nc")],
        day_offsets=np.arange(12, dtype=np.float64),
        aicen_all=structure_fields["aicen"],
        vicen_all=structure_fields["vicen"],
        vsnon_all=structure_fields["vsnon"],
        trcrn_all=None,
        tsfc_all=None,
        tinz_all=None,
        sinz_all=None,
        compression_level=2,
        source_kind="restart_fixture",
        restore_format="structure_only",
    )
    with conv.netCDF4.Dataset(structure_path, "r+") as dataset:
        dataset.variables["time"].units = "months since 2017-01-01 00:00:00"

    print(f"full_trcrn={full_path} ntrcr={layout.ntrcr}")
    print(f"structure_only={structure_path}")


if __name__ == "__main__":
    main()
