#!/usr/bin/env python3
"""Read-only validation for index-timed CICE boundary-restore forcing."""

from __future__ import annotations

import argparse
import calendar
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
import warnings

import netCDF4
import numpy as np


STRUCTURE_DIMS = ("time", "nc", "nj", "ni")
THERMO_LAYER_DIMS = ("time", "nc", "nilyr", "nj", "ni")
TRACER_DIMS = ("time", "nc", "ntrcr", "nj", "ni")
TOLERANCE = 1.0e-6


@dataclass(frozen=True)
class RuntimeConfig:
    nx: int
    ny: int
    ncat: int
    nilyr: int
    ntrcr: int
    fyear_init: int


def _parse_scalar(raw: str) -> bool | int | str:
    value = raw.strip().rstrip(",")
    if value.lower() == ".true.":
        return True
    if value.lower() == ".false.":
        return False
    if value[:1] in {"'", '"'} and value[-1:] == value[:1]:
        return value[1:-1]
    return int(value)


def _namelist_value(text: str, name: str, default: bool | int | str | None = None):
    match = re.search(
        rf"(?im)^\s*{re.escape(name)}\s*=\s*('(?:[^']*)'|\"(?:[^\"]*)\"|[^!\n/,]+)",
        text,
    )
    return default if match is None else _parse_scalar(match.group(1))


def read_runtime_config(path: Path) -> RuntimeConfig:
    text = path.read_text(encoding="utf-8")
    required = {name: _namelist_value(text, name) for name in ("nx_global", "ny_global", "ncat", "nilyr", "nslyr")}
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise ValueError(f"{path} is missing required settings: {', '.join(missing)}")

    ntrcr = 1 + int(required["nilyr"]) + int(required["nslyr"]) + int(required["nilyr"])
    ntrcr += int(bool(_namelist_value(text, "tr_iage", False)))
    ntrcr += int(bool(_namelist_value(text, "tr_FY", False)))
    if bool(_namelist_value(text, "tr_lvl", False)):
        ntrcr += 2
    if any(bool(_namelist_value(text, name, False)) for name in ("tr_pond_lvl", "tr_pond_topo", "tr_pond_sealvl")):
        ntrcr += 3
    # CICE reserves the trailing brine slot even when tr_brine is disabled.
    ntrcr += 1

    fyear_init = _namelist_value(text, "fyear_init", _namelist_value(text, "year_init"))
    if fyear_init is None:
        raise ValueError(f"{path} does not define fyear_init or year_init")
    return RuntimeConfig(
        nx=int(required["nx_global"]),
        ny=int(required["ny_global"]),
        ncat=int(required["ncat"]),
        nilyr=int(required["nilyr"]),
        ntrcr=ntrcr,
        fyear_init=int(fyear_init),
    )


def _require_dims(ds: netCDF4.Dataset, name: str, dimensions: tuple[str, ...], shape: tuple[int, ...]) -> None:
    if name not in ds.variables:
        raise ValueError(f"missing variable {name}")
    variable = ds.variables[name]
    if variable.dimensions != dimensions:
        raise ValueError(f"{name} dimensions are {variable.dimensions}, expected {dimensions}")
    if variable.shape != shape:
        raise ValueError(f"{name} shape is {variable.shape}, expected {shape}")


def _values(ds: netCDF4.Dataset, name: str) -> np.ndarray:
    return np.ma.asarray(ds.variables[name][...], dtype=np.float64).filled(np.nan)


def _check_finite(ds: netCDF4.Dataset, names: tuple[str, ...]) -> None:
    for name in names:
        if not np.all(np.isfinite(_values(ds, name))):
            raise ValueError(f"{name} contains non-finite or missing values")


def _check_time(
    ds: netCDF4.Dataset,
    data_type: str,
    fyear_init: int,
    nrecords: int,
    allow_one_based_time: bool,
) -> None:
    if "time" not in ds.variables:
        raise ValueError("missing variable time")
    time = ds.variables["time"]
    if time.dimensions != ("time",):
        raise ValueError("time must have dimensions ('time',)")
    expected = np.arange(nrecords, dtype=np.float64)
    actual = np.ma.asarray(time[...], dtype=np.float64).filled(np.nan)
    if not np.array_equal(actual, expected):
        legacy_expected = np.arange(1, nrecords + 1, dtype=np.float64)
        if allow_one_based_time and np.array_equal(actual, legacy_expected):
            warnings.warn(
                "accepting legacy one-based time values; the CICE reader uses record position, "
                "but newly converted files must use zero-based indices",
                UserWarning,
                stacklevel=2,
            )
        else:
            raise ValueError("time values must be consecutive zero-based record indices")
    expected_prefix = "days since" if data_type == "daily_netcdf" else "months since"
    units = str(getattr(time, "units", ""))
    expected_anchor = f"{fyear_init:04d}-01-01"
    if not units.lower().startswith(expected_prefix) or expected_anchor not in units:
        raise ValueError(f"time units must be '{expected_prefix} {expected_anchor} ...'")


def _record_position(when: datetime, data_type: str, fyear_init: int) -> float:
    origin = datetime(fyear_init, 1, 1)
    if when < origin:
        raise ValueError("simulation coverage begins before fyear_init")
    if data_type == "daily_netcdf":
        return (when - origin).total_seconds() / 86400.0
    month_index = 12 * (when.year - fyear_init) + when.month - 1
    seconds_into_month = (when.day - 1) * 86400 + when.hour * 3600 + when.minute * 60 + when.second
    month_seconds = calendar.monthrange(when.year, when.month)[1] * 86400
    return month_index + seconds_into_month / month_seconds


def validate_restore_file(
    forcing_file: Path,
    config: RuntimeConfig,
    data_type: str,
    cycle_year: bool,
    simulation_start: datetime | None = None,
    simulation_end: datetime | None = None,
    allow_one_based_time: bool = False,
) -> str:
    if data_type not in {"daily_netcdf", "monthly_netcdf"}:
        raise ValueError(f"unsupported data type {data_type}")

    with netCDF4.Dataset(forcing_file, "r") as ds:
        for dimension in ("time", "nc", "nj", "ni"):
            if dimension not in ds.dimensions:
                raise ValueError(f"missing dimension {dimension}")
        nrecords = len(ds.dimensions["time"])
        structure_shape = (nrecords, config.ncat, config.ny, config.nx)
        for name in ("aicen", "vicen", "vsnon"):
            _require_dims(ds, name, STRUCTURE_DIMS, structure_shape)

        has_trcrn = "trcrn" in ds.variables
        thermo_present = tuple(name in ds.variables for name in ("Tsfc", "Tinz", "Sinz"))
        if any(thermo_present) and not all(thermo_present):
            raise ValueError("incomplete thermodynamic schema: Tsfc, Tinz, and Sinz must be supplied together")
        if has_trcrn and all(thermo_present):
            raise ValueError("ambiguous schema: supply trcrn or Tsfc/Tinz/Sinz, not both")

        if has_trcrn:
            schema = "full_trcrn"
            _require_dims(
                ds,
                "trcrn",
                TRACER_DIMS,
                (nrecords, config.ncat, config.ntrcr, config.ny, config.nx),
            )
            finite_names = ("aicen", "vicen", "vsnon", "trcrn")
        elif all(thermo_present):
            schema = "bc_fields_v1"
            _require_dims(ds, "Tsfc", STRUCTURE_DIMS, structure_shape)
            layer_shape = (nrecords, config.ncat, config.nilyr, config.ny, config.nx)
            _require_dims(ds, "Tinz", THERMO_LAYER_DIMS, layer_shape)
            _require_dims(ds, "Sinz", THERMO_LAYER_DIMS, layer_shape)
            finite_names = ("aicen", "vicen", "vsnon", "Tsfc", "Tinz", "Sinz")
        else:
            schema = "structure_only"
            finite_names = ("aicen", "vicen", "vsnon")

        declared_schema = getattr(ds, "restore_format", None)
        if declared_schema is not None and declared_schema != schema:
            raise ValueError(f"restore_format={declared_schema!r} does not match detected schema {schema!r}")

        _check_finite(ds, finite_names)
        aicen = _values(ds, "aicen")
        vicen = _values(ds, "vicen")
        vsnon = _values(ds, "vsnon")
        if np.any(aicen < -TOLERANCE) or np.any(aicen > 1.0 + TOLERANCE):
            raise ValueError("aicen must be between zero and one")
        if np.any(np.sum(aicen, axis=1) > 1.0 + TOLERANCE):
            raise ValueError("category-summed aicen exceeds one")
        if np.any(vicen < -TOLERANCE) or np.any(vsnon < -TOLERANCE):
            raise ValueError("vicen and vsnon must be nonnegative")
        if np.any((aicen <= 0.0) & ((vicen > TOLERANCE) | (vsnon > TOLERANCE))):
            raise ValueError("vicen/vsnon must be zero where category concentration is zero")
        if schema == "bc_fields_v1" and np.any(_values(ds, "Sinz") < -TOLERANCE):
            raise ValueError("Sinz must be nonnegative")

        expected_records = 365 if data_type == "daily_netcdf" else 12
        if cycle_year and nrecords != expected_records:
            raise ValueError(f"cyclic {data_type} forcing must contain exactly {expected_records} records")
        if nrecords < 1:
            raise ValueError("forcing contains no records")
        _check_time(ds, data_type, config.fyear_init, nrecords, allow_one_based_time)

        if not cycle_year:
            for label, when in (("start", simulation_start), ("end", simulation_end)):
                if when is None:
                    continue
                position = _record_position(when, data_type, config.fyear_init)
                if position > nrecords - 1 + TOLERANCE:
                    raise ValueError(f"simulation {label} is outside forcing interpolation coverage")

    return schema


def parse_datetime(text: str) -> datetime:
    return datetime.fromisoformat(text.replace("Z", "+00:00")).replace(tzinfo=None)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forcing-file", required=True, type=Path)
    parser.add_argument("--ice-in", required=True, type=Path)
    parser.add_argument("--data-type", required=True, choices=("daily_netcdf", "monthly_netcdf"))
    parser.add_argument("--cycle-year", action="store_true", help="Require one exact climatological year.")
    parser.add_argument("--simulation-start", type=parse_datetime, help="ISO date/time for non-cyclic coverage validation.")
    parser.add_argument("--simulation-end", type=parse_datetime, help="ISO date/time for non-cyclic coverage validation.")
    parser.add_argument(
        "--allow-one-based-time",
        action="store_true",
        help="Accept legacy 1..N record coordinates; newly converted files must use 0..N-1.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = read_runtime_config(args.ice_in)
    schema = validate_restore_file(
        args.forcing_file,
        config,
        args.data_type,
        args.cycle_year,
        args.simulation_start,
        args.simulation_end,
        args.allow_one_based_time,
    )
    print(f"Valid CICE restore forcing: {args.forcing_file}")
    print(f"Schema: {schema}")


if __name__ == "__main__":
    main()
