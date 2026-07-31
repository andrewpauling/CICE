#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import date
from datetime import datetime
from datetime import timedelta
from pathlib import Path
import re

import netCDF4
import numpy as np
import xarray as xr
import xesmf as xe


PUNY = 1.0e-11
RHOW = 1026.0
RHOI = 917.0
RHOS = 330.0
CP_OCN = 4218.0
CP_ICE = 2106.0
LFRESH = 2.835e6 - 2.501e6
AZ1_LIQ = -18.48
BZ1_LIQ = 0.0
AZ2_LIQ = -10.3085
BZ2_LIQ = 62.4
TB_LIQ = -7.6362968855167352
AZ1P_LIQ = AZ1_LIQ / 1000.0
BZ1P_LIQ = BZ1_LIQ / 1000.0
AZ2P_LIQ = AZ2_LIQ / 1000.0
BZ2P_LIQ = BZ2_LIQ / 1000.0


@dataclass(frozen=True)
class IceInConfig:
    ncat: int
    nilyr: int
    nslyr: int
    ktherm: int
    tr_iage: bool
    tr_fy: bool
    tr_lvl: bool
    tr_pond_lvl: bool
    tr_pond_topo: bool
    tr_pond_sealvl: bool
    tr_snow: bool
    tr_iso: bool
    tr_aero: bool
    tr_fsd: bool
    tr_brine: bool

    @property
    def tr_pond(self) -> bool:
        return self.tr_pond_lvl or self.tr_pond_topo or self.tr_pond_sealvl


@dataclass(frozen=True)
class TracerLayout:
    ntrcr: int
    nt_tsfc: int
    nt_qice: int
    nt_qsno: int
    nt_sice: int
    nt_iage: int | None
    nt_fy: int | None
    nt_alvl: int | None
    nt_vlvl: int | None
    nt_apnd: int | None
    nt_hpnd: int | None
    nt_ipnd: int | None
    nt_fbri: int | None


def _parse_fortran_value(text: str) -> bool | int | str:
    value = text.strip().rstrip(",")
    lowered = value.lower()
    if lowered == ".true.":
        return True
    if lowered == ".false.":
        return False
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1]
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1]
    return int(value)


def _read_ice_in_value(text: str, name: str, default: bool | int | str | None = None) -> bool | int | str | None:
    pattern = re.compile(
        rf"(?im)^\s*{re.escape(name)}\s*=\s*('(?:[^']*)'|\"(?:[^\"]*)\"|[^!\n/,]+)"
    )
    match = pattern.search(text)
    if match is None:
        return default
    return _parse_fortran_value(match.group(1))


def parse_ice_in(path: Path) -> IceInConfig:
    text = path.read_text(encoding="utf-8")
    config = IceInConfig(
        ncat=int(_read_ice_in_value(text, "ncat")),
        nilyr=int(_read_ice_in_value(text, "nilyr")),
        nslyr=int(_read_ice_in_value(text, "nslyr")),
        ktherm=int(_read_ice_in_value(text, "ktherm")),
        tr_iage=bool(_read_ice_in_value(text, "tr_iage", False)),
        tr_fy=bool(_read_ice_in_value(text, "tr_FY", False)),
        tr_lvl=bool(_read_ice_in_value(text, "tr_lvl", False)),
        tr_pond_lvl=bool(_read_ice_in_value(text, "tr_pond_lvl", False)),
        tr_pond_topo=bool(_read_ice_in_value(text, "tr_pond_topo", False)),
        tr_pond_sealvl=bool(_read_ice_in_value(text, "tr_pond_sealvl", False)),
        tr_snow=bool(_read_ice_in_value(text, "tr_snow", False)),
        tr_iso=bool(_read_ice_in_value(text, "tr_iso", False)),
        tr_aero=bool(_read_ice_in_value(text, "tr_aero", False)),
        tr_fsd=bool(_read_ice_in_value(text, "tr_fsd", False)),
        tr_brine=bool(_read_ice_in_value(text, "tr_brine", False)),
    )
    if config.ktherm != 2:
        raise ValueError(f"Only ktherm=2 mushy thermodynamics are supported for trcrn synthesis, found {config.ktherm}")
    unsupported = []
    if config.tr_snow:
        unsupported.append("tr_snow")
    if config.tr_iso:
        unsupported.append("tr_iso")
    if config.tr_aero:
        unsupported.append("tr_aero")
    if config.tr_fsd:
        unsupported.append("tr_fsd")
    if unsupported:
        raise ValueError(f"Unsupported tracer options for synthetic trcrn output: {', '.join(unsupported)}")
    return config


def build_tracer_layout(config: IceInConfig) -> TracerLayout:
    cursor = 1
    nt_tsfc = cursor
    cursor += 1
    nt_qice = cursor
    cursor += config.nilyr
    nt_qsno = cursor
    cursor += config.nslyr
    nt_sice = cursor
    cursor += config.nilyr

    nt_iage = None
    if config.tr_iage:
        nt_iage = cursor
        cursor += 1

    nt_fy = None
    if config.tr_fy:
        nt_fy = cursor
        cursor += 1

    nt_alvl = None
    nt_vlvl = None
    if config.tr_lvl:
        nt_alvl = cursor
        cursor += 1
        nt_vlvl = cursor
        cursor += 1

    nt_apnd = None
    nt_hpnd = None
    nt_ipnd = None
    if config.tr_pond:
        nt_apnd = cursor
        cursor += 1
        nt_hpnd = cursor
        cursor += 1
        nt_ipnd = cursor
        cursor += 1

    # CICE allocates trcrn with one reserved trailing slot for the brine tracer
    # even when tr_brine is disabled; the restore helper validates against that size.
    nt_fbri = cursor

    return TracerLayout(
        ntrcr=cursor,
        nt_tsfc=nt_tsfc,
        nt_qice=nt_qice,
        nt_qsno=nt_qsno,
        nt_sice=nt_sice,
        nt_iage=nt_iage,
        nt_fy=nt_fy,
        nt_alvl=nt_alvl,
        nt_vlvl=nt_vlvl,
        nt_apnd=nt_apnd,
        nt_hpnd=nt_hpnd,
        nt_ipnd=nt_ipnd,
        nt_fbri=nt_fbri,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert standalone daily CICE history or restart files into a boundary-restore "
            "NetCDF file for the coupled MITgcm+CICE restore helper."
        )
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--history-dir",
        help="Directory containing standalone daily iceh.YYYY-MM-DD.nc files.",
    )
    source_group.add_argument(
        "--restart-dir",
        help="Directory containing standalone daily iced.YYYY-MM-DD-00000.nc restart files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output restore-forcing NetCDF file.",
    )
    parser.add_argument(
        "--start-date",
        required=True,
        help="First day to include, in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--end-date",
        required=True,
        help="Last day to include, in YYYY-MM-DD format.",
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
    parser.add_argument(
        "--target-grid-file",
        help=(
            "Target CICE grid NetCDF file used for regridding. If omitted, grid_file is "
            "read from --ice-in."
        ),
    )
    parser.add_argument(
        "--source-grid-file",
        help=(
            "Optional NetCDF file carrying source-grid lon/lat coordinates. "
            "Useful for restart input, which usually omits TLON/TLAT."
        ),
    )
    parser.add_argument(
        "--regrid-method",
        choices=["bilinear", "nearest_s2d"],
        default="bilinear",
        help="xESMF regridding method for remapping standalone history fields to the coupled grid.",
    )
    parser.add_argument(
        "--reuse-weights",
        action="store_true",
        help="Reuse an existing xESMF weight file next to the output when available.",
    )
    parser.add_argument(
        "--ice-in",
        help="Path to the active CICE ice_in namelist used to reconstruct and write trcrn.",
    )
    parser.add_argument(
        "--restore-format",
        choices=["full_trcrn", "bc_fields_v1", "structure_only"],
        default="full_trcrn",
        help=(
            "Boundary thermo representation to write. full_trcrn writes packed tracers, "
            "bc_fields_v1 writes Tsfc/Tinz/Sinz, and structure_only writes only aicen/vicen/vsnon."
        ),
    )
    return parser.parse_args()


def parse_day(text: str) -> date:
    return datetime.strptime(text, "%Y-%m-%d").date()


def resolve_target_grid_file(target_grid_file: str | None, ice_in_path: Path | None) -> Path:
    if target_grid_file:
        return Path(target_grid_file)
    if ice_in_path is None:
        raise ValueError("provide --target-grid-file or --ice-in so the target CICE grid can be resolved")
    text = ice_in_path.read_text(encoding="utf-8")
    grid_file = _read_ice_in_value(text, "grid_file")
    if not isinstance(grid_file, str) or not grid_file:
        raise ValueError(f"{ice_in_path} does not define grid_file")
    path = Path(grid_file)
    return path if path.is_absolute() else ice_in_path.parent / path


def iter_days(start_day: date, end_day: date) -> list[date]:
    if end_day < start_day:
        raise ValueError("end-date must not be earlier than start-date")
    days: list[date] = []
    current = start_day
    while current <= end_day:
        days.append(current)
        current += timedelta(days=1)
    return days


def history_path(history_dir: Path, day: date) -> Path:
    return history_dir / f"iceh.{day.isoformat()}.nc"


def restart_path(restart_dir: Path, day: date) -> Path:
    return restart_dir / f"iced.{day.isoformat()}-00000.nc"


def build_source_grid(path: Path) -> xr.Dataset:
    with xr.open_dataset(path) as ds:
        if "TLON" in ds.variables and "TLAT" in ds.variables:
            lon = np.mod(np.asarray(ds["TLON"][:, :], dtype=np.float64), 360.0)
            lat = np.asarray(ds["TLAT"][:, :], dtype=np.float64)
        elif "ulon" in ds.variables and "ulat" in ds.variables:
            lon = np.mod(np.rad2deg(np.asarray(ds["ulon"][:, :], dtype=np.float64)), 360.0)
            lat = np.rad2deg(np.asarray(ds["ulat"][:, :], dtype=np.float64))
        else:
            raise ValueError(f"{path} does not contain source-grid lon/lat variables")
        source_grid = xr.Dataset(
            data_vars={
                "lon": xr.DataArray(
                    lon,
                    dims=("nj", "ni"),
                ),
                "lat": xr.DataArray(
                    lat,
                    dims=("nj", "ni"),
                ),
            }
        )
    return source_grid


def build_target_grid(path: Path) -> xr.Dataset:
    with xr.open_dataset(path) as ds:
        target_grid = xr.Dataset(
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
    return target_grid


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
    lead_shape = records.shape[:-2]
    sample_index = (0,) * len(lead_shape)
    sample = regridder(xr.DataArray(records[sample_index], dims=("nj", "ni")))
    target_ny, target_nx = sample.shape
    output = np.zeros((*lead_shape, target_ny, target_nx), dtype=np.float32)

    for index in np.ndindex(lead_shape):
        field = xr.DataArray(records[index], dims=("nj", "ni"))
        remapped = np.asarray(regridder(field), dtype=np.float32)
        output[index] = np.where(np.isfinite(remapped), remapped, 0.0)

    return output


def _read_history_var(ds: netCDF4.Dataset, name: str) -> np.ndarray:
    variable = ds.variables[name]
    data = np.ma.filled(variable[0, ...], np.nan).astype(np.float32)
    fill_candidates = [getattr(variable, "_FillValue", None), getattr(variable, "missing_value", None)]
    for fill_value in fill_candidates:
        if fill_value is None:
            continue
        fill_array = np.asarray(fill_value, dtype=np.float32)
        data = np.where(np.isclose(data, fill_array, rtol=0.0, atol=0.0), np.nan, data)
    return np.where(np.isfinite(data), data, 0.0).astype(np.float32)


def _read_optional_history_var(ds: netCDF4.Dataset, name: str) -> np.ndarray | None:
    if name not in ds.variables:
        return None
    return _read_history_var(ds, name)


def _read_restart_var(ds: netCDF4.Dataset, name: str) -> np.ndarray:
    variable = ds.variables[name]
    data = np.ma.filled(variable[:], np.nan).astype(np.float32)
    fill_candidates = [getattr(variable, "_FillValue", None), getattr(variable, "missing_value", None)]
    for fill_value in fill_candidates:
        if fill_value is None:
            continue
        fill_array = np.asarray(fill_value, dtype=np.float32)
        data = np.where(np.isclose(data, fill_array, rtol=0.0, atol=0.0), np.nan, data)
    return np.where(np.isfinite(data), data, 0.0).astype(np.float32)


def _read_optional_restart_var(ds: netCDF4.Dataset, *names: str) -> np.ndarray | None:
    for name in names:
        if name in ds.variables:
            return _read_restart_var(ds, name)
    return None


def _read_restart_layers(ds: netCDF4.Dataset, prefix: str, count: int) -> np.ndarray | None:
    if count == 0:
        return None
    layers: list[np.ndarray] = []
    for layer_index in range(1, count + 1):
        name = f"{prefix}{layer_index:03d}"
        if name not in ds.variables:
            raise ValueError(f"Restart file is missing required variable {name}")
        layers.append(_read_restart_var(ds, name))
    return np.stack(layers, axis=1).astype(np.float32)


def read_history_day(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray | None]]:
    with netCDF4.Dataset(path) as ds:
        time_var = np.asarray(ds.variables["time"][:], dtype=np.float64)
        if time_var.size != 1:
            raise ValueError(f"Expected one time sample in {path}, found {time_var.size}")

        fields: dict[str, np.ndarray | None] = {
            "aicen": _read_history_var(ds, "aicen_d"),
            "vicen": _read_history_var(ds, "vicen_d"),
            "vsnon": _read_history_var(ds, "vsnon_d"),
            "tsfc": _read_optional_history_var(ds, "Tsfc_d"),
            "iage": _read_optional_history_var(ds, "iage_d"),
            "fy": _read_optional_history_var(ds, "FY_d"),
            "alvl": _read_optional_history_var(ds, "alvl_d"),
            "vlvl": _read_optional_history_var(ds, "vlvl_d"),
            "apondn": _read_optional_history_var(ds, "apondn_d"),
            "hpondn": _read_optional_history_var(ds, "hpondn_d"),
            "ipondn": _read_optional_history_var(ds, "ipondn_d"),
            "tinz": _read_optional_history_var(ds, "Tinz_d"),
            "sinz": _read_optional_history_var(ds, "Sinz_d"),
            "tsnz": _read_optional_history_var(ds, "Tsnz_d"),
            "fbrine": _read_optional_history_var(ds, "fbrine_d"),
        }
        return time_var, fields


def read_restart_day(path: Path, config: IceInConfig | None) -> dict[str, np.ndarray | None]:
    with netCDF4.Dataset(path) as ds:
        vsnon = _read_restart_var(ds, "vsnon")
        fields: dict[str, np.ndarray | None] = {
            "aicen": _read_restart_var(ds, "aicen"),
            "vicen": _read_restart_var(ds, "vicen"),
            "vsnon": vsnon,
            "tsfc": None,
            "qice": None,
            "qsno": None,
            "tsno": None,
            "sice": None,
            "iage": None,
            "fy": None,
            "alvl": None,
            "vlvl": None,
            "apnd": None,
            "hpnd": None,
            "ipnd": None,
            "fbri": None,
        }

        if config is None:
            return fields

        fields["tsfc"] = _read_restart_var(ds, "Tsfcn")
        fields["qice"] = _read_restart_layers(ds, "qice", config.nilyr)
        fields["qsno"] = _read_restart_layers(ds, "qsno", config.nslyr)
        fields["sice"] = _read_restart_layers(ds, "sice", config.nilyr)
        if fields["qsno"] is not None:
            snow_mask = vsnon[:, None, :, :] > PUNY
            fields["tsno"] = np.where(
                snow_mask,
                np.minimum(icepack_snow_temperature_from_enthalpy(fields["qsno"]), 0.0),
                0.0,
            ).astype(np.float32)

        if config.tr_iage:
            fields["iage"] = _read_optional_restart_var(ds, "iage")
        if config.tr_fy:
            fields["fy"] = _read_optional_restart_var(ds, "FY")
        if config.tr_lvl:
            fields["alvl"] = _read_optional_restart_var(ds, "alvl")
            fields["vlvl"] = _read_optional_restart_var(ds, "vlvl")
        if config.tr_pond:
            fields["apnd"] = _read_optional_restart_var(ds, "apnd")
            fields["hpnd"] = _read_optional_restart_var(ds, "hpnd")
            fields["ipnd"] = _read_optional_restart_var(ds, "ipnd")
        fields["fbri"] = _read_optional_restart_var(ds, "fbri", "fbrine")

        return fields


def liquidus_brine_salinity_mush(z_tin: np.ndarray) -> np.ndarray:
    j1_liq = BZ1_LIQ / AZ1_LIQ
    k1_liq = 1.0 / 1000.0
    l1_liq = (1.0 + BZ1P_LIQ) / AZ1_LIQ
    j2_liq = BZ2_LIQ / AZ2_LIQ
    k2_liq = 1.0 / 1000.0
    l2_liq = (1.0 + BZ2P_LIQ) / AZ2_LIQ
    t_high = z_tin > TB_LIQ
    sbr = np.where(
        t_high,
        (z_tin + j1_liq) / (k1_liq * z_tin + l1_liq),
        (z_tin + j2_liq) / (k2_liq * z_tin + l2_liq),
    )
    return np.where(z_tin <= 0.0, sbr, 0.0)


def icepack_mushy_liquid_fraction(z_tin: np.ndarray, z_sin: np.ndarray) -> np.ndarray:
    sbr = np.maximum(liquidus_brine_salinity_mush(z_tin), PUNY)
    return z_sin / np.maximum(sbr, z_sin)


def icepack_enthalpy_mush(z_tin: np.ndarray, z_sin: np.ndarray) -> np.ndarray:
    phi = icepack_mushy_liquid_fraction(z_tin, z_sin)
    return (
        phi * (CP_OCN * RHOW - CP_ICE * RHOI) * z_tin
        + RHOI * CP_ICE * z_tin
        - (1.0 - phi) * RHOI * LFRESH
    )


def icepack_enthalpy_snow(z_tsn: np.ndarray) -> np.ndarray:
    return -RHOS * (-CP_ICE * z_tsn + LFRESH)


def icepack_snow_temperature_from_enthalpy(z_qsn: np.ndarray) -> np.ndarray:
    return (LFRESH + z_qsn / RHOS) / CP_ICE


def _broadcast_2d_to_categories(field: np.ndarray, ncat: int) -> np.ndarray:
    return np.repeat(field[:, None, :, :], ncat, axis=1)


def _broadcast_masked(field: np.ndarray, category_mask: np.ndarray) -> np.ndarray:
    return np.where(category_mask, _broadcast_2d_to_categories(field, category_mask.shape[1]), 0.0)


def _expand_area_weighted_field(field: np.ndarray, weight: np.ndarray, category_mask: np.ndarray) -> np.ndarray:
    expanded = np.divide(
        field[:, None, :, :],
        weight,
        out=np.zeros_like(weight, dtype=np.float32),
        where=weight > PUNY,
    )
    return np.where(category_mask, expanded, 0.0)


def _require_field(field_name: str, field_value: np.ndarray | None) -> np.ndarray:
    if field_value is None:
        raise ValueError(f"History files do not contain required field {field_name}")
    return field_value


def synthesize_trcrn(
    history_fields: dict[str, np.ndarray | None],
    config: IceInConfig,
    layout: TracerLayout,
) -> np.ndarray:
    aicen_all = _require_field("aicen_d", history_fields["aicen"])
    vsnon_all = _require_field("vsnon_d", history_fields["vsnon"])
    tinz_all = _require_field("Tinz_d", history_fields["tinz"])
    sinz_all = _require_field("Sinz_d", history_fields["sinz"])

    ntime, ncat, ny, nx = aicen_all.shape
    if config.ncat != ncat:
        raise ValueError(f"ice_in ncat={config.ncat} does not match history ncat={ncat}")
    if tinz_all.shape[2] != config.nilyr:
        raise ValueError(f"ice_in nilyr={config.nilyr} does not match history nkice={tinz_all.shape[2]}")
    if sinz_all.shape[2] != config.nilyr:
        raise ValueError(f"ice_in nilyr={config.nilyr} does not match history nkice={sinz_all.shape[2]}")

    trcrn_all = np.zeros((ntime, ncat, layout.ntrcr, ny, nx), dtype=np.float32)
    category_mask = aicen_all > PUNY
    layer_mask = category_mask[:, :, None, :, :]
    total_aice = np.sum(aicen_all, axis=1, keepdims=True)
    total_vice = np.sum(_require_field("vicen_d", history_fields["vicen"]), axis=1, keepdims=True)

    tsnz_all = history_fields["tsnz"]
    tsfc_all = history_fields["tsfc"]
    if tsnz_all is not None:
        tsfc_cat = np.where(vsnon_all > PUNY, tsnz_all[:, :, 0, :, :], tinz_all[:, :, 0, :, :])
    elif tsfc_all is not None:
        tsfc_cat = _broadcast_masked(tsfc_all, category_mask)
    else:
        tsfc_cat = np.where(category_mask, tinz_all[:, :, 0, :, :], 0.0)
    trcrn_all[:, :, layout.nt_tsfc - 1, :, :] = np.asarray(tsfc_cat, dtype=np.float32)

    qice_all = np.where(layer_mask, icepack_enthalpy_mush(tinz_all, sinz_all), 0.0)
    for layer_index in range(config.nilyr):
        tracer_index = layout.nt_qice + layer_index - 1
        trcrn_all[:, :, tracer_index, :, :] = np.asarray(qice_all[:, :, layer_index, :, :], dtype=np.float32)

    if config.nslyr > 0:
        if tsnz_all is not None:
            if tsnz_all.shape[2] != config.nslyr:
                raise ValueError(f"ice_in nslyr={config.nslyr} does not match history nksnow={tsnz_all.shape[2]}")
            qsno_all = np.where(vsnon_all[:, :, None, :, :] > PUNY, icepack_enthalpy_snow(tsnz_all), 0.0)
        else:
            qsno_surface = np.where(vsnon_all > PUNY, icepack_enthalpy_snow(tsfc_cat), 0.0)
            qsno_all = np.repeat(qsno_surface[:, :, None, :, :], config.nslyr, axis=2)
        for layer_index in range(config.nslyr):
            tracer_index = layout.nt_qsno + layer_index - 1
            trcrn_all[:, :, tracer_index, :, :] = np.asarray(qsno_all[:, :, layer_index, :, :], dtype=np.float32)

    for layer_index in range(config.nilyr):
        tracer_index = layout.nt_sice + layer_index - 1
        trcrn_all[:, :, tracer_index, :, :] = np.asarray(
            np.where(category_mask, sinz_all[:, :, layer_index, :, :], 0.0),
            dtype=np.float32,
        )

    if layout.nt_iage is not None:
        iage_all = history_fields["iage"]
        if iage_all is not None:
            trcrn_all[:, :, layout.nt_iage - 1, :, :] = np.asarray(
                _broadcast_masked(iage_all * (365.0 * 86400.0), category_mask),
                dtype=np.float32,
            )

    if layout.nt_fy is not None:
        trcrn_all[:, :, layout.nt_fy - 1, :, :] = np.asarray(
            np.where(category_mask, 0.0, 0.0),
            dtype=np.float32,
        )

    if layout.nt_alvl is not None:
        alvl_cat = np.where(category_mask, 1.0, 0.0)
        trcrn_all[:, :, layout.nt_alvl - 1, :, :] = np.asarray(alvl_cat, dtype=np.float32)

    if layout.nt_vlvl is not None:
        vlvl_cat = np.where(category_mask, 1.0, 0.0)
        trcrn_all[:, :, layout.nt_vlvl - 1, :, :] = np.asarray(vlvl_cat, dtype=np.float32)

    if layout.nt_apnd is not None:
        trcrn_all[:, :, layout.nt_apnd - 1, :, :] = np.asarray(
            np.where(category_mask, 0.0, 0.0),
            dtype=np.float32,
        )

    if layout.nt_hpnd is not None:
        trcrn_all[:, :, layout.nt_hpnd - 1, :, :] = np.asarray(
            np.where(category_mask, 0.0, 0.0),
            dtype=np.float32,
        )

    if layout.nt_ipnd is not None:
        trcrn_all[:, :, layout.nt_ipnd - 1, :, :] = np.asarray(
            np.where(category_mask, 0.0, 0.0),
            dtype=np.float32,
        )

    if layout.nt_fbri is not None:
        fbrine_all = history_fields["fbrine"]
        if fbrine_all is not None:
            fbri_cat = np.where(category_mask, np.maximum(fbrine_all, PUNY), 0.0)
        else:
            # CICE reserves this trailing tracer slot even when tr_brine is disabled,
            # and the thermo step still expects positive fbri values where ice exists.
            fbri_cat = np.where(category_mask, 1.0, 0.0)
        trcrn_all[:, :, layout.nt_fbri - 1, :, :] = np.asarray(fbri_cat, dtype=np.float32)

    return trcrn_all


def build_trcrn_from_restart_fields(
    restart_fields: dict[str, np.ndarray | None],
    config: IceInConfig,
    layout: TracerLayout,
) -> np.ndarray:
    aicen_all = _require_field("aicen", restart_fields["aicen"])
    ntime, ncat, ny, nx = aicen_all.shape
    if ncat != config.ncat:
        raise ValueError(f"ice_in ncat={config.ncat} does not match restart ncat={ncat}")

    category_mask = aicen_all > PUNY
    trcrn_all = np.zeros((ntime, ncat, layout.ntrcr, ny, nx), dtype=np.float32)

    tsfc_all = _require_field("Tsfcn", restart_fields["tsfc"])
    trcrn_all[:, :, layout.nt_tsfc - 1, :, :] = np.asarray(np.where(category_mask, tsfc_all, 0.0), dtype=np.float32)

    qice_all = _require_field("qice", restart_fields["qice"])
    if qice_all.shape[2] != config.nilyr:
        raise ValueError(f"Restart qice layers={qice_all.shape[2]} do not match ice_in nilyr={config.nilyr}")
    for layer_index in range(config.nilyr):
        tracer_index = layout.nt_qice + layer_index - 1
        trcrn_all[:, :, tracer_index, :, :] = np.asarray(
            np.where(category_mask, qice_all[:, :, layer_index, :, :], 0.0),
            dtype=np.float32,
        )

    if config.nslyr > 0:
        tsno_all = _require_field("tsno", restart_fields["tsno"])
        if tsno_all.shape[2] != config.nslyr:
            raise ValueError(f"Restart snow layers={tsno_all.shape[2]} do not match ice_in nslyr={config.nslyr}")
        vsnon_all = _require_field("vsnon", restart_fields["vsnon"])
        snow_mask = vsnon_all[:, :, None, :, :] > PUNY
        qsno_all = np.where(
            snow_mask,
            icepack_enthalpy_snow(np.minimum(tsno_all, 0.0)),
            -RHOS * LFRESH,
        )
        for layer_index in range(config.nslyr):
            tracer_index = layout.nt_qsno + layer_index - 1
            trcrn_all[:, :, tracer_index, :, :] = np.asarray(
                np.where(category_mask, qsno_all[:, :, layer_index, :, :], 0.0),
                dtype=np.float32,
            )

    sice_all = _require_field("sice", restart_fields["sice"])
    if sice_all.shape[2] != config.nilyr:
        raise ValueError(f"Restart sice layers={sice_all.shape[2]} do not match ice_in nilyr={config.nilyr}")
    for layer_index in range(config.nilyr):
        tracer_index = layout.nt_sice + layer_index - 1
        trcrn_all[:, :, tracer_index, :, :] = np.asarray(
            np.where(category_mask, sice_all[:, :, layer_index, :, :], 0.0),
            dtype=np.float32,
        )

    if layout.nt_iage is not None:
        iage_all = restart_fields["iage"]
        iage_cat = np.where(category_mask, iage_all, 0.0) if iage_all is not None else np.zeros_like(aicen_all, dtype=np.float32)
        trcrn_all[:, :, layout.nt_iage - 1, :, :] = np.asarray(iage_cat, dtype=np.float32)

    if layout.nt_fy is not None:
        fy_all = restart_fields["fy"]
        fy_cat = np.where(category_mask, fy_all, 0.0) if fy_all is not None else np.zeros_like(aicen_all, dtype=np.float32)
        trcrn_all[:, :, layout.nt_fy - 1, :, :] = np.asarray(fy_cat, dtype=np.float32)

    if layout.nt_alvl is not None:
        alvl_all = restart_fields["alvl"]
        alvl_cat = np.where(category_mask, alvl_all, 1.0) if alvl_all is not None else np.where(category_mask, 1.0, 0.0)
        trcrn_all[:, :, layout.nt_alvl - 1, :, :] = np.asarray(alvl_cat, dtype=np.float32)

    if layout.nt_vlvl is not None:
        vlvl_all = restart_fields["vlvl"]
        vlvl_cat = np.where(category_mask, vlvl_all, 1.0) if vlvl_all is not None else np.where(category_mask, 1.0, 0.0)
        trcrn_all[:, :, layout.nt_vlvl - 1, :, :] = np.asarray(vlvl_cat, dtype=np.float32)

    if layout.nt_apnd is not None:
        apnd_all = restart_fields["apnd"]
        apnd_cat = np.where(category_mask, apnd_all, 0.0) if apnd_all is not None else np.zeros_like(aicen_all, dtype=np.float32)
        trcrn_all[:, :, layout.nt_apnd - 1, :, :] = np.asarray(apnd_cat, dtype=np.float32)

    if layout.nt_hpnd is not None:
        hpnd_all = restart_fields["hpnd"]
        hpnd_cat = np.where(category_mask, hpnd_all, 0.0) if hpnd_all is not None else np.zeros_like(aicen_all, dtype=np.float32)
        trcrn_all[:, :, layout.nt_hpnd - 1, :, :] = np.asarray(hpnd_cat, dtype=np.float32)

    if layout.nt_ipnd is not None:
        ipnd_all = restart_fields["ipnd"]
        ipnd_cat = np.where(category_mask, ipnd_all, 0.0) if ipnd_all is not None else np.zeros_like(aicen_all, dtype=np.float32)
        trcrn_all[:, :, layout.nt_ipnd - 1, :, :] = np.asarray(ipnd_cat, dtype=np.float32)

    fbri_all = restart_fields["fbri"]
    if fbri_all is not None:
        fbri_cat = np.where(category_mask, np.maximum(fbri_all, PUNY), 0.0)
    else:
        fbri_cat = np.where(category_mask, 1.0, 0.0)
    trcrn_all[:, :, layout.nt_fbri - 1, :, :] = np.asarray(fbri_cat, dtype=np.float32)

    return trcrn_all


def build_bc_fields_from_history(
    history_fields: dict[str, np.ndarray | None],
    config: IceInConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aicen_all = _require_field("aicen_d", history_fields["aicen"])
    tinz_all = _require_field("Tinz_d", history_fields["tinz"])
    sinz_all = _require_field("Sinz_d", history_fields["sinz"])
    vsnon_all = _require_field("vsnon_d", history_fields["vsnon"])

    ntime, ncat, ny, nx = aicen_all.shape
    if config.ncat != ncat:
        raise ValueError(f"ice_in ncat={config.ncat} does not match history ncat={ncat}")
    if tinz_all.shape[2] != config.nilyr:
        raise ValueError(f"ice_in nilyr={config.nilyr} does not match history nkice={tinz_all.shape[2]}")
    if sinz_all.shape[2] != config.nilyr:
        raise ValueError(f"ice_in nilyr={config.nilyr} does not match history nkice={sinz_all.shape[2]}")

    category_mask = aicen_all > PUNY
    tsnz_all = history_fields["tsnz"]
    tsfc_all = history_fields["tsfc"]
    if tsnz_all is not None:
        if tsnz_all.shape[2] != config.nslyr:
            raise ValueError(f"ice_in nslyr={config.nslyr} does not match history nksnow={tsnz_all.shape[2]}")
        tsfc_cat = np.where(vsnon_all > PUNY, tsnz_all[:, :, 0, :, :], tinz_all[:, :, 0, :, :])
    elif tsfc_all is not None:
        tsfc_cat = _broadcast_masked(tsfc_all, category_mask)
    else:
        tsfc_cat = np.where(category_mask, tinz_all[:, :, 0, :, :], 0.0)

    tsfc_cat = np.asarray(np.where(category_mask, np.minimum(tsfc_cat, 0.0), 0.0), dtype=np.float32)
    tinz_bc = np.asarray(np.where(category_mask[:, :, None, :, :], tinz_all, 0.0), dtype=np.float32)
    sinz_bc = np.asarray(np.where(category_mask[:, :, None, :, :], sinz_all, 0.0), dtype=np.float32)
    return tsfc_cat, tinz_bc, sinz_bc


def sanitize_regridded_restart_fields(restart_fields: dict[str, np.ndarray | None]) -> None:
    aicen_all = restart_fields.get("aicen")
    if aicen_all is None:
        return

    category_mask = aicen_all > PUNY
    layer_mask = category_mask[:, :, None, :, :]

    for field_name in ("aicen", "vicen", "vsnon", "tsfc", "iage", "fy", "alvl", "vlvl", "apnd", "hpnd", "ipnd", "fbri"):
        field_value = restart_fields.get(field_name)
        if field_value is not None:
            restart_fields[field_name] = np.where(category_mask, field_value, 0.0).astype(np.float32)

    for field_name in ("qice", "qsno", "tsno", "sice"):
        field_value = restart_fields.get(field_name)
        if field_value is not None:
            restart_fields[field_name] = np.where(layer_mask, field_value, 0.0).astype(np.float32)


def write_output(
    output_path: Path,
    source_paths: list[Path],
    day_offsets: np.ndarray,
    aicen_all: np.ndarray,
    vicen_all: np.ndarray,
    vsnon_all: np.ndarray,
    trcrn_all: np.ndarray | None,
    tsfc_all: np.ndarray | None,
    tinz_all: np.ndarray | None,
    sinz_all: np.ndarray | None,
    compression_level: int,
    source_kind: str,
    restore_format: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ntime, ncat, ny, nx = aicen_all.shape

    with netCDF4.Dataset(output_path, "w") as ds:
        ds.createDimension("time", ntime)
        ds.createDimension("nc", ncat)
        ds.createDimension("nj", ny)
        ds.createDimension("ni", nx)
        if trcrn_all is not None:
            ds.createDimension("ntrcr", trcrn_all.shape[2])
        if tinz_all is not None or sinz_all is not None:
            ds.createDimension("nilyr", (tinz_all if tinz_all is not None else sinz_all).shape[2])

        time_var = ds.createVariable("time", "f8", ("time",))
        time_var[:] = day_offsets
        time_var.long_name = "boundary restore record day offsets"
        time_var.units = f"days since {source_paths[0].stem.split('.', 1)[1]} 00:00:00"
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

        # Match native CICE time-dependent file order so nf90_get_var with
        # start/count=(x,y,ncat,time) sees the expected Fortran dimensions.
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

        if trcrn_all is not None:
            trcrn_var = ds.createVariable("trcrn", "f4", ("time", "nc", "ntrcr", "nj", "ni"), **var_kwargs)
            trcrn_var[:, :, :, :, :] = trcrn_all
            trcrn_var.long_name = "category tracers"
            trcrn_var.units = "mixed"

        if tsfc_all is not None:
            tsfc_var = ds.createVariable("Tsfc", "f4", ("time", "nc", "nj", "ni"), **var_kwargs)
            tsfc_var[:, :, :, :] = tsfc_all
            tsfc_var.long_name = "category surface temperature"
            tsfc_var.units = "degC"

        if tinz_all is not None:
            tinz_var = ds.createVariable("Tinz", "f4", ("time", "nc", "nilyr", "nj", "ni"), **var_kwargs)
            tinz_var[:, :, :, :, :] = tinz_all
            tinz_var.long_name = "category internal ice temperature"
            tinz_var.units = "degC"

        if sinz_all is not None:
            sinz_var = ds.createVariable("Sinz", "f4", ("time", "nc", "nilyr", "nj", "ni"), **var_kwargs)
            sinz_var[:, :, :, :, :] = sinz_all
            sinz_var.long_name = "category internal ice salinity"
            sinz_var.units = "ppt"

        ds.setncattr("restore_source", f"standalone_cice_{source_kind}")
        ds.setncattr("restore_format", restore_format)
        ds.setncattr("restore_trcrn_included", 1 if trcrn_all is not None else 0)
        if restore_format == "bc_fields_v1":
            ds.setncattr(
                "restore_note",
                "This restore file carries aicen/vicen/vsnon plus Tsfc/Tinz/Sinz for thermo reconstruction inside the coupled restore helper.",
            )
        elif trcrn_all is None:
            ds.setncattr(
                "restore_note",
                "This restore file contains aicen/vicen/vsnon only; the coupled restore helper will leave trcrn unchanged.",
            )
        else:
            ds.setncattr(
                "restore_note",
                f"This restore file includes trcrn populated from standalone {source_kind} data and ice_in tracer settings.",
            )
        ds.setncattr("restore_input_dir", str(source_paths[0].parent))
        ds.setncattr("restore_start_file", str(source_paths[0]))
        ds.setncattr("restore_end_file", str(source_paths[-1]))


def main() -> None:
    args = parse_args()
    history_dir = Path(args.history_dir) if args.history_dir else None
    restart_dir = Path(args.restart_dir) if args.restart_dir else None
    output_path = Path(args.output)
    source_grid_file = Path(args.source_grid_file) if args.source_grid_file else None
    ice_in_path = Path(args.ice_in) if args.ice_in else None
    source_kind = "restart" if restart_dir is not None else "history"

    source_dir = restart_dir if restart_dir is not None else history_dir
    if source_dir is None or not source_dir.exists():
        raise FileNotFoundError(source_dir)
    if source_grid_file is not None and not source_grid_file.exists():
        raise FileNotFoundError(source_grid_file)
    if ice_in_path is not None and not ice_in_path.exists():
        raise FileNotFoundError(ice_in_path)
    target_grid_file = resolve_target_grid_file(args.target_grid_file, ice_in_path)
    if not target_grid_file.exists():
        raise FileNotFoundError(target_grid_file)
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(output_path)

    start_day = parse_day(args.start_date)
    end_day = parse_day(args.end_date)
    if (start_day.month, start_day.day) != (1, 1):
        raise ValueError("start-date must be 1 January to match the reader's index-based time convention")
    requested_days = iter_days(start_day, end_day)
    if restart_dir is not None:
        source_paths = [restart_path(restart_dir, day) for day in requested_days]
    else:
        source_paths = [history_path(history_dir, day) for day in requested_days]

    for path in source_paths:
        if not path.exists():
            raise FileNotFoundError(path)

    config = parse_ice_in(ice_in_path) if ice_in_path is not None and args.restore_format != "structure_only" else None
    layout = build_tracer_layout(config) if config is not None else None

    field_records: dict[str, list[np.ndarray]] = {}
    if restart_dir is not None:
        day_offsets = np.asarray([(day - start_day).days for day in requested_days], dtype=np.float64)
        restart_field_names = [
            "aicen",
            "vicen",
            "vsnon",
            "tsfc",
            "qice",
            "qsno",
            "tsno",
            "sice",
            "iage",
            "fy",
            "alvl",
            "vlvl",
            "apnd",
            "hpnd",
            "ipnd",
            "fbri",
        ]
        for path in source_paths:
            fields = read_restart_day(path, config)
            for field_name, field_value in fields.items():
                if field_value is None:
                    continue
                field_records.setdefault(field_name, []).append(field_value)
        stacked_fields: dict[str, np.ndarray | None] = {}
        for field_name in restart_field_names:
            records = field_records.get(field_name)
            stacked_fields[field_name] = np.stack(records, axis=0) if records else None
    else:
        for path in source_paths:
            _, fields = read_history_day(path)
            for field_name, field_value in fields.items():
                if field_value is None:
                    continue
                field_records.setdefault(field_name, []).append(field_value)
        stacked_fields = {}
        for field_name, records in field_records.items():
            stacked_fields[field_name] = np.stack(records, axis=0)
        for optional_name in [
            "aicen",
            "vicen",
            "vsnon",
            "tsfc",
            "iage",
            "fy",
            "alvl",
            "vlvl",
            "apondn",
            "hpondn",
            "ipondn",
            "tinz",
            "sinz",
            "tsnz",
            "fbrine",
        ]:
            stacked_fields.setdefault(optional_name, None)
        # CICE maps model time directly to zero-based record indices; source
        # history time units and epochs are intentionally not propagated.
        day_offsets = np.arange(len(requested_days), dtype=np.float64)

    source_grid = None
    if source_grid_file is not None:
        source_grid = build_source_grid(source_grid_file)
    else:
        try:
            source_grid = build_source_grid(source_paths[0])
        except ValueError:
            if restart_dir is None:
                raise
            history_candidate = history_path(restart_dir.parent / "history", start_day)
            if not history_candidate.exists():
                raise ValueError(
                    "Restart input requires --source-grid-file or an available sibling history file carrying TLON/TLAT"
                )
            source_grid = build_source_grid(history_candidate)

    target_grid = build_target_grid(target_grid_file)
    weight_file = output_path.with_suffix(f".{args.regrid_method}.weights.nc")
    regridder = build_regridder(
        source_grid=source_grid,
        target_grid=target_grid,
        method=args.regrid_method,
        weight_file=weight_file,
        reuse_weights=args.reuse_weights,
    )

    for field_name, field_value in list(stacked_fields.items()):
        if field_value is not None:
            stacked_fields[field_name] = maybe_regrid_records(field_value, regridder)

    if restart_dir is not None:
        sanitize_regridded_restart_fields(stacked_fields)

    if restart_dir is not None:
        aicen_all = _require_field("aicen", stacked_fields["aicen"])
        vicen_all = _require_field("vicen", stacked_fields["vicen"])
        vsnon_all = _require_field("vsnon", stacked_fields["vsnon"])
    else:
        aicen_all = _require_field("aicen_d", stacked_fields["aicen"])
        vicen_all = _require_field("vicen_d", stacked_fields["vicen"])
        vsnon_all = _require_field("vsnon_d", stacked_fields["vsnon"])

    trcrn_all = None
    tsfc_all = None
    tinz_all = None
    sinz_all = None
    if args.restore_format == "full_trcrn":
        if config is None or layout is None:
            raise ValueError("--restore-format full_trcrn requires --ice-in")
        if restart_dir is not None:
            trcrn_all = build_trcrn_from_restart_fields(stacked_fields, config, layout)
        else:
            trcrn_all = synthesize_trcrn(stacked_fields, config, layout)
    elif args.restore_format == "bc_fields_v1":
        if config is None:
            raise ValueError("--restore-format bc_fields_v1 requires --ice-in")
        if restart_dir is not None:
            raise ValueError("bc_fields_v1 export is currently implemented for history input only")
        tsfc_all, tinz_all, sinz_all = build_bc_fields_from_history(stacked_fields, config)

    write_output(
        output_path=output_path,
        source_paths=source_paths,
        day_offsets=day_offsets,
        aicen_all=aicen_all,
        vicen_all=vicen_all,
        vsnon_all=vsnon_all,
        trcrn_all=trcrn_all,
        tsfc_all=tsfc_all,
        tinz_all=tinz_all,
        sinz_all=sinz_all,
        compression_level=args.compression_level,
        source_kind=source_kind,
        restore_format=args.restore_format,
    )

    print(f"Wrote {output_path}")
    print(f"Records: {aicen_all.shape[0]}")
    print(f"Grid: ni={aicen_all.shape[3]} nj={aicen_all.shape[2]} nc={aicen_all.shape[1]}")
    print(f"Restore format: {args.restore_format}")
    if trcrn_all is not None:
        print(f"Tracers: ntrcr={trcrn_all.shape[2]}")
    if tinz_all is not None:
        print(f"Thermo layers: nilyr={tinz_all.shape[2]}")
    print(f"Regrid method: {args.regrid_method}")
    print(f"Weight file: {weight_file}")


if __name__ == "__main__":
    main()
