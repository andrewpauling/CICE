#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
RUNS = HERE / "runs"
INPUTS = HERE / "inputs"

CASE_INFO = {
    "pn_parent_remap": ("parent", "remap"),
    "pn_parent_upwind": ("parent", "upwind"),
    "pn_nest_remap": ("nest", "remap"),
    "pn_nest_upwind": ("nest", "upwind"),
}


def replace_namelist(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=\s*).*$")
    text, count = pattern.subn(rf"\g<1>{value}", text)
    if count != 1:
        raise RuntimeError(f"Expected exactly one {name} setting, found {count}")
    return text


def replace_setting(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^setenv\s+{re.escape(name)}\s+.*$")
    text, count = pattern.subn(f"setenv {name}  {value}", text)
    if count != 1:
        raise RuntimeError(f"Expected exactly one {name} setting, found {count}")
    return text


def configure(case_name: str, role: str, advection: str) -> None:
    case = CASES / case_name
    if not (case / "ice_in").exists():
        raise FileNotFoundError(f"Run cice.setup first: {case}")

    parent = role == "parent"
    paths = {
        "ice_ic": INPUTS / ("parent_initial_2017-01-01.nc" if parent else "nest_initial_2017-01-01.nc"),
        "grid_file": Path("/nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/grid_pskrips.nc") if parent else INPUTS / "nest_grid.nc",
        "kmt_file": Path("/nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/kmt_pskrips.nc") if parent else INPUTS / "nest_kmt.nc",
        "bathymetry_file": Path("/nfs/scratch/paulinan/cice-dirs/input/CICE_data/grid/pskrips/bath_pskrips.nc") if parent else INPUTS / "nest_bath.nc",
    }
    nx, ny = ((210, 240) if parent else (84, 80))

    ice_in_path = case / "ice_in"
    text = ice_in_path.read_text()
    settings = {
        "year_init": "2017",
        "month_init": "1",
        "day_init": "1",
        "use_leap_years": ".false.",
        "dt": "86400.0",
        "npt_unit": "'1'",
        "npt": "5",
        "ndtd": "1",
        "runtype": "'initial'",
        "restart_ext": ".false." if parent else ".true.",
        "use_restart_time": ".false.",
        "dumpfreq": "'1','x','x','x','x'",
        "histfreq": "'x','x','x','x','x'",
        "ice_ic": f"'{paths['ice_ic']}'",
        "grid_format": "'nc'",
        "grid_type": "'regional'",
        "grid_ice": "'C'",
        "grid_file": f"'{paths['grid_file']}'",
        "kmt_file": f"'{paths['kmt_file']}'",
        "bathymetry_file": f"'{paths['bathymetry_file']}'",
        "ktherm": "-1",
        "kdyn": "0",
        "kridge": "-1",
        "ktransport": "1",
        "advection": f"'{advection}'",
        "calc_strair": ".false.",
        "atm_data_type": "'default'",
        "ocn_data_type": "'default'",
        "restore_ice": ".false." if parent else ".true.",
        "restore_ice_data_type": "'daily_netcdf'",
        "restore_ice_data_file": f"'{INPUTS / ('parent_' + advection + '_nest_boundary.nc')}'",
        "restore_ice_cycle_year": ".false.",
        "restore_ice_use_west": ".true.",
        "restore_ice_use_east": ".true.",
        "restore_ice_use_south": ".true.",
        "restore_ice_use_north": ".true.",
        "trestore": "0",
        "nprocs": "30",
        "nx_global": str(nx),
        "ny_global": str(ny),
        "block_size_x": "7",
        "block_size_y": "8",
        "max_blocks": "-1",
        "distribution_type": "'roundrobin'",
        "distribution_wght": "'block'",
        "ew_boundary_type": "'open'",
        "ns_boundary_type": "'open'",
        "fyear_init": "2017",
    }
    for name, value in settings.items():
        text = replace_namelist(text, name, value)
    ice_in_path.write_text(text)

    settings_path = case / "cice.settings"
    text = settings_path.read_text()
    text = replace_setting(text, "ICE_RUNDIR", str(RUNS / case_name))
    text = replace_setting(text, "ICE_CPPDEFS", "-DRESTORE_DIAGNOSTICS")
    if not re.search(r"(?m)^setenv\s+ICE_QUIETMODE\s+", text):
        text = re.sub(
            r"(?m)^(setenv\s+ICE_CPPDEFS\s+.*)$",
            r"\1\nsetenv ICE_QUIETMODE  false",
            text,
            count=1,
        )
    settings_path.write_text(text)

    print(f"Configured {case_name}: {role}, {advection}, {nx}x{ny}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="*")
    args = parser.parse_args()
    for case_name in (args.cases or sorted(CASE_INFO)):
        if case_name not in CASE_INFO:
            parser.error(f"unknown case {case_name}; choose from {', '.join(sorted(CASE_INFO))}")
        configure(case_name, *CASE_INFO[case_name])


if __name__ == "__main__":
    main()
