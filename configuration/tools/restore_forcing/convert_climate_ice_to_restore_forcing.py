#!/usr/bin/env python3

from __future__ import annotations

import argparse
from datetime import date
from datetime import datetime
from pathlib import Path
import re

import netCDF4
import numpy as np
import xarray as xr
import xesmf as xe


CAT_LIMITS = np.asarray([0.0, 0.6445072, 1.391433, 2.470179, 4.567288, 1.0e8], dtype=np.float64)
AICE_MIN = 0.01
HICE_MIN = 0.01
SNOW_MIN = 0.01
RHOS = 330.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert coarse monthly climate-model sea-ice fields into a category-resolved "
            "CICE boundary-restore NetCDF file without trcrn."
        )
    )
    parser.add_argument("--aice-file", required=True, help="Monthly sea-ice concentration NetCDF file.")
    parser.add_argument("--vice-file", required=True, help="Monthly sea-ice volume-per-area NetCDF file.")
    parser.add_argument("--vsno-file", help="Monthly snow volume-per-area NetCDF file.")
    parser.add_argument("--output", required=True, help="Path to the output restore-forcing NetCDF file.")
    parser.add_argument("--start-month", required=True, help="First month to include, in YYYY-MM format.")
    parser.add_argument("--end-month", required=True, help="Last month to include, in YYYY-MM format.")
    parser.add_argument(
        "--target-grid-file",
        help=(
            "Target CICE grid NetCDF file used for regridding. If omitted, grid_file is "
            "read from --ice-in."
        ),
    )
    parser.add_argument(
        "--ice-in",
        help="Active CICE ice_in namelist used to locate grid_file and confirm ncat=5.",
    )
    parser.add_argument(
        "--aice-var",
        default="siconc",
        help="Variable name for sea-ice concentration in the climate file.",
    )
    parser.add_argument(
        "--vice-var",
        default="sivol",
        help="Variable name for sea-ice volume-per-area in the climate file.",
    )
    parser.add_argument(
        "--vsno-var",
        default="sisnthick",
        help="Variable name for snow volume-per-area in the climate file.",
    )
    parser.add_argument(
        "--snow-mass-file",
        help="Monthly snow mass-per-sea-ice-area NetCDF file.",
    )
    parser.add_argument(
        "--snow-conc-file",
        help="Monthly snow concentration-on-sea-ice NetCDF file.",
    )
    parser.add_argument(
        "--snow-mass-var",
        default="sisnmass",
        help="Variable name for snow mass-per-sea-ice-area in the climate file.",
    )
    parser.add_argument(
        "--snow-conc-var",
        default="sisnconc",
        help="Variable name for snow concentration-on-sea-ice in the climate file.",
    )
    parser.add_argument(
        "--time-name",
        default="time",
        help="Time coordinate name in the climate files.",
    )
    parser.add_argument(
        "--lon-name",
        default="longitude",
        help="Longitude variable name in the climate files.",
    )
    parser.add_argument(
        "--lat-name",
        default="latitude",
        help="Latitude variable name in the climate files.",
    )
    parser.add_argument(
        "--regrid-method",
        choices=["bilinear", "nearest_s2d"],
        default="bilinear",
        help="xESMF regridding method for remapping climate fields to the coupled grid.",
    )
    parser.add_argument(
        "--reuse-weights",
        action="store_true",
        help="Reuse an existing xESMF weight file next to the output when available.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=1,
        help="NetCDF zlib compression level for output variables.",
    )
    args = parser.parse_args()
    if not args.vsno_file and not (args.snow_mass_file and args.snow_conc_file):
        parser.error("provide either --vsno-file or both --snow-mass-file and --snow-conc-file")
    if args.vsno_file and (args.snow_mass_file or args.snow_conc_file):
        parser.error("use either --vsno-file or the --snow-mass-file/--snow-conc-file pair, not both")
    if bool(args.snow_mass_file) != bool(args.snow_conc_file):
        parser.error("--snow-mass-file and --snow-conc-file must be provided together")
    return args


def parse_month(text: str) -> date:
    return datetime.strptime(text, "%Y-%m").date().replace(day=1)


def read_ice_in_value(path: Path, name: str) -> str | None:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf"(?im)^\s*{re.escape(name)}\s*=\s*('(?:[^']*)'|\"(?:[^\"]*)\"|[^!\n/,]+)",
        text,
    )
    if match is None:
        return None
    return match.group(1).strip().rstrip(",").strip("'\"")


def resolve_target_grid_file(target_grid_file: str | None, ice_in_path: Path | None) -> Path:
    if target_grid_file:
        return Path(target_grid_file)
    if ice_in_path is None:
        raise ValueError("provide --target-grid-file or --ice-in so the target CICE grid can be resolved")
    grid_file = read_ice_in_value(ice_in_path, "grid_file")
    if not grid_file:
        raise ValueError(f"{ice_in_path} does not define grid_file")
    path = Path(grid_file)
    return path if path.is_absolute() else ice_in_path.parent / path


def month_range(start_month: date, end_month: date) -> list[date]:
    if end_month < start_month:
        raise ValueError("end-month must not be earlier than start-month")

    months: list[date] = []
    current = start_month
    while current <= end_month:
        months.append(current)
        if current.month == 12:
            current = current.replace(year=current.year + 1, month=1)
        else:
            current = current.replace(month=current.month + 1)
    return months


def build_target_grid(path: Path) -> xr.Dataset:
    with xr.open_dataset(path) as ds:
        return xr.Dataset(
            data_vars={
                "lon": xr.DataArray(
                    np.mod(np.rad2deg(np.asarray(ds["ulon"][:, :], dtype=np.float64)), 360.0),
                    dims=("y", "x"),
                ),
                "lat": xr.DataArray(
                    np.rad2deg(np.asarray(ds["ulat"][:, :], dtype=np.float64)),
                    dims=("y", "x"),
                ),
            }
        )


def build_source_grid(path: Path, lon_name: str, lat_name: str) -> xr.Dataset:
    with xr.open_dataset(path) as ds:
        if lon_name not in ds.variables or lat_name not in ds.variables:
            raise ValueError(f"{path} does not contain {lon_name}/{lat_name}")

        lon = np.asarray(ds[lon_name][...], dtype=np.float64)
        lat = np.asarray(ds[lat_name][...], dtype=np.float64)

        if lon.ndim == 1 and lat.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        elif lon.ndim == 2 and lat.ndim == 2:
            lon2d, lat2d = lon, lat
        else:
            raise ValueError(f"Unsupported lon/lat rank in {path}: {lon.ndim} / {lat.ndim}")

        return xr.Dataset(
            data_vars={
                "lon": xr.DataArray(np.mod(lon2d, 360.0), dims=("nj", "ni")),
                "lat": xr.DataArray(lat2d, dims=("nj", "ni")),
            }
        )


def build_regridder(
    source_grid: xr.Dataset,
    target_grid: xr.Dataset,
    method: str,
    weight_file: Path,
    reuse_weights: bool,
) -> xe.Regridder:
    return xe.Regridder(
        source_grid,
        target_grid,
        method,
        periodic=False,
        filename=str(weight_file),
        reuse_weights=reuse_weights and weight_file.exists(),
    )


def maybe_regrid_records(records: np.ndarray, regridder: xe.Regridder) -> np.ndarray:
    ntime = records.shape[0]
    sample = np.asarray(regridder(xr.DataArray(records[0], dims=("nj", "ni"))), dtype=np.float32)
    output = np.zeros((ntime, sample.shape[0], sample.shape[1]), dtype=np.float32)

    for index in range(ntime):
        remapped = np.asarray(regridder(xr.DataArray(records[index], dims=("nj", "ni"))), dtype=np.float32)
        output[index] = np.where(np.isfinite(remapped), remapped, 0.0)

    return output


def _fill_missing(data: np.ndarray, variable: netCDF4.Variable) -> np.ndarray:
    result = np.ma.filled(data, np.nan).astype(np.float32)
    for attr_name in ("_FillValue", "missing_value"):
        fill_value = getattr(variable, attr_name, None)
        if fill_value is None:
            continue
        fill_value = np.asarray(fill_value, dtype=np.float32)
        result = np.where(np.isclose(result, fill_value, rtol=0.0, atol=0.0), np.nan, result)
    return np.where(np.isfinite(result), result, 0.0).astype(np.float32)


def read_monthly_field(
    path: Path,
    var_name: str,
    time_name: str,
    requested_months: list[date],
) -> tuple[list[date], np.ndarray, str, str]:
    with netCDF4.Dataset(path) as ds:
        if var_name not in ds.variables:
            raise ValueError(f"{path} does not contain variable {var_name}")
        if time_name not in ds.variables:
            raise ValueError(f"{path} does not contain time variable {time_name}")

        time_var = ds.variables[time_name]
        time_values = np.asarray(time_var[:], dtype=np.float64)
        time_units = getattr(time_var, "units", None)
        time_calendar = getattr(time_var, "calendar", "standard")
        if time_units is None:
            raise ValueError(f"{path} time variable {time_name} is missing units")

        decoded = netCDF4.num2date(time_values, time_units, calendar=time_calendar)
        month_lookup: dict[tuple[int, int], int] = {}
        for index, current_time in enumerate(decoded):
            month_lookup[(int(current_time.year), int(current_time.month))] = index

        selected_indices: list[int] = []
        selected_months: list[date] = []
        for month in requested_months:
            key = (month.year, month.month)
            if key not in month_lookup:
                raise ValueError(f"{path} is missing requested month {month:%Y-%m}")
            selected_indices.append(month_lookup[key])
            selected_months.append(month)

        variable = ds.variables[var_name]
        data = _fill_missing(variable[selected_indices, ...], variable)
        if data.ndim != 3:
            raise ValueError(f"Expected {var_name} in {path} to have dimensions (time,j,i), found {data.shape}")
        units = getattr(variable, "units", "")
        return selected_months, data, units, time_calendar


def normalize_aice(values: np.ndarray, units: str) -> np.ndarray:
    normalized_units = units.strip().lower()
    if normalized_units in {"%", "percent", "percentage"}:
        result = values / 100.0
    elif normalized_units in {"1", "", "fraction"}:
        result = values
    else:
        raise ValueError(f"Unsupported concentration units: {units!r}")
    return np.clip(result, 0.0, 1.0).astype(np.float32)


def normalize_length(values: np.ndarray, units: str, label: str) -> np.ndarray:
    normalized_units = units.strip().lower()
    if normalized_units in {"m", "meter", "meters", "metre", "metres", "1"}:
        result = values
    elif normalized_units in {"cm", "centimeter", "centimeters", "centimetre", "centimetres"}:
        result = values / 100.0
    elif normalized_units in {"mm", "millimeter", "millimeters", "millimetre", "millimetres"}:
        result = values / 1000.0
    else:
        raise ValueError(f"Unsupported {label} units: {units!r}")
    return np.maximum(result, 0.0).astype(np.float32)


def normalize_snow_mass(values: np.ndarray, units: str) -> np.ndarray:
    normalized_units = units.strip().lower().replace(" ", "")
    if normalized_units in {"kgm-2", "kgm^-2", "kg/m2", "kgm**-2", "kgm-2s"}:
        result = values
    elif normalized_units in {"gm-2", "gm^-2", "g/m2"}:
        result = values / 1000.0
    else:
        raise ValueError(f"Unsupported snow mass units: {units!r}")
    return np.maximum(result, 0.0).astype(np.float32)


def derive_vsno_from_snow_fields(aice: np.ndarray, snow_mass: np.ndarray, snow_conc: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    vsno = aice * (snow_mass / RHOS)
    snow_cover_mask = snow_conc > AICE_MIN
    vsno = np.where(snow_cover_mask, vsno, 0.0).astype(np.float32)
    diagnostics = {
        "snow_cover_cells": int(np.count_nonzero(snow_cover_mask)),
        "snow_masked_zero_cover_cells": int(np.count_nonzero((snow_conc <= AICE_MIN) & (snow_mass > 0.0))),
    }
    return vsno, diagnostics


def apply_category_template(
    aice: np.ndarray,
    vice: np.ndarray,
    vsno: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    ntime, ny, nx = aice.shape
    ncat = len(CAT_LIMITS) - 1

    aicen = np.zeros((ntime, ncat, ny, nx), dtype=np.float32)
    vicen = np.zeros_like(aicen)
    vsnon = np.zeros_like(aicen)

    hice = np.divide(vice, aice, out=np.zeros_like(vice), where=aice > AICE_MIN)
    valid_mask = (aice > AICE_MIN) & (hice > HICE_MIN)
    vsno_floor = np.where(valid_mask, np.maximum(vsno, SNOW_MIN * aice), 0.0).astype(np.float32)

    diagnostics: dict[str, int] = {
        "valid_ice_cells": int(np.count_nonzero(valid_mask)),
        "dropped_sparse_cells": int(np.count_nonzero((aice > 0.0) & (aice <= AICE_MIN))),
        "dropped_thin_cells": int(np.count_nonzero((aice > AICE_MIN) & (hice <= HICE_MIN))),
        "snow_floor_cells": int(np.count_nonzero(valid_mask & (vsno < SNOW_MIN * aice))),
    }

    for category_index in range(ncat):
        lower = CAT_LIMITS[category_index]
        upper = CAT_LIMITS[category_index + 1]
        category_mask = valid_mask & (hice > lower) & (hice <= upper)
        aicen[:, category_index, :, :] = np.where(category_mask, aice, 0.0)
        vicen[:, category_index, :, :] = np.where(category_mask, vice, 0.0)
        vsnon[:, category_index, :, :] = np.where(category_mask, vsno_floor, 0.0)
        diagnostics[f"category_{category_index + 1}_cells"] = int(np.count_nonzero(category_mask))

    return aicen, vicen, vsnon, diagnostics


def write_output(
    output_path: Path,
    requested_months: list[date],
    aicen_all: np.ndarray,
    vicen_all: np.ndarray,
    vsnon_all: np.ndarray,
    compression_level: int,
    aice_file: Path,
    vice_file: Path,
    snow_source_desc: str,
    regrid_method: str,
    diagnostics: dict[str, int],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ntime, ncat, ny, nx = aicen_all.shape
    month_offsets = np.asarray(
        [12 * (month.year - requested_months[0].year) + (month.month - requested_months[0].month) for month in requested_months],
        dtype=np.float64,
    )

    with netCDF4.Dataset(output_path, "w") as ds:
        ds.createDimension("time", ntime)
        ds.createDimension("nc", ncat)
        ds.createDimension("nj", ny)
        ds.createDimension("ni", nx)

        time_var = ds.createVariable("time", "f8", ("time",))
        time_var[:] = month_offsets
        time_var.long_name = "boundary restore record month offsets"
        time_var.units = f"months since {requested_months[0]:%Y-%m}-01 00:00:00"
        time_var.calendar = "proleptic_gregorian"

        nc_var = ds.createVariable("NCAT", "f4", ("nc",))
        nc_var[:] = np.arange(1, ncat + 1, dtype=np.float32)
        nc_var.long_name = "category index"
        nc_var.units = "1"

        var_kwargs = {
            "zlib": True,
            "complevel": compression_level,
            "shuffle": True,
        }

        aicen_var = ds.createVariable("aicen", "f4", ("time", "nc", "nj", "ni"), **var_kwargs)
        aicen_var[:, :, :, :] = aicen_all
        aicen_var.long_name = "ice area, categories"
        aicen_var.units = "1"

        vicen_var = ds.createVariable("vicen", "f4", ("time", "nc", "nj", "ni"), **var_kwargs)
        vicen_var[:, :, :, :] = vicen_all
        vicen_var.long_name = "ice volume, categories"
        vicen_var.units = "m"

        vsnon_var = ds.createVariable("vsnon", "f4", ("time", "nc", "nj", "ni"), **var_kwargs)
        vsnon_var[:, :, :, :] = vsnon_all
        vsnon_var.long_name = "snow volume on ice, categories"
        vsnon_var.units = "m"

        ds.setncattr("restore_source", "climate_monthly_template")
        ds.setncattr("restore_format", "structure_only")
        ds.setncattr("restore_trcrn_included", 0)
        ds.setncattr("restore_note", "Monthly climate restore file with METROMS-style fixed category template and no trcrn.")
        ds.setncattr("restore_category_template", "fixed_single_category_by_total_thickness")
        ds.setncattr("restore_category_limits_m", ",".join(str(limit) for limit in CAT_LIMITS))
        ds.setncattr("restore_threshold_aice", AICE_MIN)
        ds.setncattr("restore_threshold_hice_m", HICE_MIN)
        ds.setncattr("restore_snow_floor_m", SNOW_MIN)
        ds.setncattr("restore_regrid_method", regrid_method)
        ds.setncattr("restore_aice_file", str(aice_file))
        ds.setncattr("restore_vice_file", str(vice_file))
        ds.setncattr("restore_vsno_source", snow_source_desc)
        ds.setncattr("restore_start_month", f"{requested_months[0]:%Y-%m}")
        ds.setncattr("restore_end_month", f"{requested_months[-1]:%Y-%m}")
        for key, value in diagnostics.items():
            ds.setncattr(f"restore_diag_{key}", value)


def main() -> None:
    args = parse_args()

    aice_file = Path(args.aice_file)
    vice_file = Path(args.vice_file)
    vsno_file = Path(args.vsno_file) if args.vsno_file else None
    snow_mass_file = Path(args.snow_mass_file) if args.snow_mass_file else None
    snow_conc_file = Path(args.snow_conc_file) if args.snow_conc_file else None
    output_path = Path(args.output)
    ice_in_path = Path(args.ice_in) if args.ice_in else None

    if ice_in_path is not None and not ice_in_path.exists():
        raise FileNotFoundError(ice_in_path)
    target_grid_file = resolve_target_grid_file(args.target_grid_file, ice_in_path)
    if ice_in_path is not None:
        ncat = read_ice_in_value(ice_in_path, "ncat")
        if ncat is None or int(ncat) != len(CAT_LIMITS) - 1:
            raise ValueError(f"{ice_in_path} must define ncat={len(CAT_LIMITS) - 1} for the fixed category template")

    required_paths = [aice_file, vice_file, target_grid_file]
    if vsno_file is not None:
        required_paths.append(vsno_file)
    if snow_mass_file is not None:
        required_paths.append(snow_mass_file)
    if snow_conc_file is not None:
        required_paths.append(snow_conc_file)

    for path in required_paths:
        if not path.exists():
            raise FileNotFoundError(path)
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(output_path)

    start_month = parse_month(args.start_month)
    end_month = parse_month(args.end_month)
    if start_month.month != 1:
        raise ValueError("start-month must be January to match the reader's index-based time convention")
    requested_months = month_range(start_month, end_month)

    aice_months, aice_raw, aice_units, calendar = read_monthly_field(
        aice_file,
        args.aice_var,
        args.time_name,
        requested_months,
    )
    vice_months, vice_raw, vice_units, vice_calendar = read_monthly_field(
        vice_file,
        args.vice_var,
        args.time_name,
        requested_months,
    )
    if vsno_file is not None:
        vsno_months, vsno_raw, vsno_units, vsno_calendar = read_monthly_field(
            vsno_file,
            args.vsno_var,
            args.time_name,
            requested_months,
        )
        if aice_months != vice_months or aice_months != vsno_months:
            raise ValueError("Input climate files do not share the same selected monthly timestamps")
        if not (calendar == vice_calendar == vsno_calendar):
            raise ValueError("Input climate files do not share the same calendar")
    else:
        snow_conc_months, snow_conc_raw, snow_conc_units, snow_conc_calendar = read_monthly_field(
            snow_conc_file,
            args.snow_conc_var,
            args.time_name,
            requested_months,
        )
        snow_mass_months, snow_mass_raw, snow_mass_units, snow_mass_calendar = read_monthly_field(
            snow_mass_file,
            args.snow_mass_var,
            args.time_name,
            requested_months,
        )
        if aice_months != vice_months or aice_months != snow_conc_months or aice_months != snow_mass_months:
            raise ValueError("Input climate files do not share the same selected monthly timestamps")
        if not (calendar == vice_calendar == snow_conc_calendar == snow_mass_calendar):
            raise ValueError("Input climate files do not share the same calendar")

    source_grid = build_source_grid(aice_file, args.lon_name, args.lat_name)
    target_grid = build_target_grid(target_grid_file)
    weight_file = output_path.with_suffix(f".{args.regrid_method}.weights.nc")
    regridder = build_regridder(source_grid, target_grid, args.regrid_method, weight_file, args.reuse_weights)

    aice = normalize_aice(maybe_regrid_records(aice_raw, regridder), aice_units)
    vice = normalize_length(maybe_regrid_records(vice_raw, regridder), vice_units, "ice volume")
    if vsno_file is not None:
        vsno = normalize_length(maybe_regrid_records(vsno_raw, regridder), vsno_units, "snow volume")
        snow_diagnostics: dict[str, int] = {}
        snow_source_desc = str(vsno_file)
    else:
        snow_conc = normalize_aice(maybe_regrid_records(snow_conc_raw, regridder), snow_conc_units)
        snow_mass = normalize_snow_mass(maybe_regrid_records(snow_mass_raw, regridder), snow_mass_units)
        vsno, snow_diagnostics = derive_vsno_from_snow_fields(aice, snow_mass, snow_conc)
        snow_source_desc = f"snow_mass={snow_mass_file};snow_conc={snow_conc_file}"

    aicen_all, vicen_all, vsnon_all, diagnostics = apply_category_template(aice, vice, vsno)
    diagnostics.update(snow_diagnostics)

    write_output(
        output_path=output_path,
        requested_months=requested_months,
        aicen_all=aicen_all,
        vicen_all=vicen_all,
        vsnon_all=vsnon_all,
        compression_level=args.compression_level,
        aice_file=aice_file,
        vice_file=vice_file,
        snow_source_desc=snow_source_desc,
        regrid_method=args.regrid_method,
        diagnostics=diagnostics,
    )

    print(f"Wrote {output_path}")
    print(f"Records: {aicen_all.shape[0]}")
    print(f"Grid: ni={aicen_all.shape[3]} nj={aicen_all.shape[2]} nc={aicen_all.shape[1]}")
    print(f"Calendar: {calendar}")
    print(f"Regrid method: {args.regrid_method}")
    print(f"Weight file: {weight_file}")
    for key in sorted(diagnostics):
        print(f"{key}: {diagnostics[key]}")


if __name__ == "__main__":
    main()
