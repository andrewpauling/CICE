#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import netCDF4
import numpy as np


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
SOURCE_ROOT = Path("/nfs/scratch/paulinan/cice-dirs/input/CICE_data")
SOURCE_RESTART = SOURCE_ROOT / "ic/pskrips/iced_pskrips_2017-02-01.onice_hi0p5_hs0p01.thin_tail_itd_minarea0p05.nc"
SOURCE_GRID = SOURCE_ROOT / "grid/pskrips/grid_pskrips.nc"
SOURCE_KMT = SOURCE_ROOT / "grid/pskrips/kmt_pskrips.nc"
SOURCE_BATH = SOURCE_ROOT / "grid/pskrips/bath_pskrips.nc"

X0, X1 = 100, 184
Y0, Y1 = 140, 220
NX_PARENT, NY_PARENT = 210, 240


def create_like(source: Path, destination: Path, slices: dict[str, slice] | None = None) -> None:
    slices = slices or {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(source) as src, netCDF4.Dataset(destination, "w", format=src.data_model) as dst:
        dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
        for name, dim in src.dimensions.items():
            if name in slices:
                start, stop, step = slices[name].indices(len(dim))
                if step != 1:
                    raise ValueError(f"Only unit-stride slices are supported for {name}")
                size = stop - start
            else:
                size = None if dim.isunlimited() else len(dim)
            dst.createDimension(name, size)

        for name, var in src.variables.items():
            fill_value = getattr(var, "_FillValue", None)
            kwargs = {"fill_value": fill_value} if fill_value is not None else {}
            out = dst.createVariable(name, var.datatype, var.dimensions, **kwargs)
            attrs = {key: var.getncattr(key) for key in var.ncattrs() if key != "_FillValue"}
            if attrs:
                out.setncatts(attrs)
            selection = tuple(slices.get(dim, slice(None)) for dim in var.dimensions)
            out[:] = var[selection]


def build_synthetic_parent(path: Path) -> None:
    create_like(SOURCE_RESTART, path)
    with netCDF4.Dataset(SOURCE_KMT) as ds:
        wet = np.asarray(ds.variables["kmt"][:] > 0)

    yy, xx = np.meshgrid(np.arange(NY_PARENT), np.arange(NX_PARENT), indexing="ij")
    concentration = (
        0.42
        + 0.24 * np.sin(2.0 * np.pi * xx / 19.0)
        + 0.16 * np.cos(2.0 * np.pi * yy / 27.0)
        + 0.10 * np.sin(2.0 * np.pi * (xx + yy) / 37.0)
    )
    concentration = np.where(wet, np.clip(concentration, 0.02, 0.92), 0.0)
    fractions = np.asarray([0.10, 0.19, 0.29, 0.25, 0.17], dtype=np.float64)
    thickness = np.asarray([0.30, 0.90, 1.80, 3.40, 5.50], dtype=np.float64)

    with netCDF4.Dataset(path, "r+") as ds:
        aicen_original = np.asarray(ds.variables["aicen"][:])
        reference_points: list[tuple[int, int]] = []
        for category in range(aicen_original.shape[0]):
            candidates = np.argwhere(aicen_original[category] > 1.0e-5)
            if not len(candidates):
                raise RuntimeError(f"No reference state for category {category + 1}")
            reference_points.append(tuple(int(value) for value in candidates[len(candidates) // 2]))

        # Replicate a physically initialized tracer state for every category.
        # Structural fields are overwritten below.  This keeps the restart
        # internally valid while the disabled thermodynamics leave transport
        # as the only evolving process.
        for name, var in ds.variables.items():
            if var.dimensions != ("ncat", "y", "x"):
                continue
            if name in {"aicen", "vicen", "vsnon"}:
                continue
            values = np.asarray(var[:])
            for category, (ref_y, ref_x) in enumerate(reference_points):
                values[category, :, :] = values[category, ref_y, ref_x]
            if np.issubdtype(values.dtype, np.floating):
                values[:, ~wet] = 0.0
            var[:] = values

        aicen = fractions[:, None, None] * concentration[None, :, :]
        vicen = aicen * thickness[:, None, None]
        vsnon = aicen * 0.08
        ds.variables["aicen"][:] = aicen
        ds.variables["vicen"][:] = vicen
        ds.variables["vsnon"][:] = vsnon

        # Uniform grid-relative flow makes the velocity halo continuation
        # exact, so state/technical transport halos are what the test probes.
        ds.variables["uvel"][:] = np.where(wet, 0.025, 0.0)
        ds.variables["vvel"][:] = np.where(wet, 0.015, 0.0)
        # A C-grid restart carries the normal face components separately.
        # Upwind transport consumes these fields directly, whereas remap also
        # retains the legacy U-point velocity representation.  The production
        # initial-condition file predates the C-grid face fields, so add them
        # explicitly for a matched transport-only experiment.
        if "uvelE" not in ds.variables:
            ds.createVariable("uvelE", "f8", ("y", "x"), fill_value=1.0e30)
        if "vvelN" not in ds.variables:
            ds.createVariable("vvelN", "f8", ("y", "x"), fill_value=1.0e30)
        ds.variables["uvelE"][:] = np.where(wet, 0.025, 0.0)
        ds.variables["vvelN"][:] = np.where(wet, 0.015, 0.0)
        ds.setncattr("parent_nest_test", "smooth transport-only P-SKRIPS parent state")
        ds.setncattr("parent_nest_uvel_m_per_s", 0.025)
        ds.setncattr("parent_nest_vvel_m_per_s", 0.015)


def main() -> None:
    INPUTS.mkdir(parents=True, exist_ok=True)
    parent_restart = INPUTS / "parent_initial_2017-01-01.nc"
    build_synthetic_parent(parent_restart)

    create_like(SOURCE_GRID, INPUTS / "nest_grid.nc", {"x": slice(X0, X1), "y": slice(Y0, Y1)})
    create_like(SOURCE_KMT, INPUTS / "nest_kmt.nc", {"x": slice(X0, X1), "y": slice(Y0, Y1)})
    create_like(SOURCE_BATH, INPUTS / "nest_bath.nc", {"ni": slice(X0, X1), "nj": slice(Y0, Y1)})
    # The transport-only nest has no dynamics step to populate velocity
    # halos.  Read a one-cell extended restart so its C-grid west/south face
    # velocities match the corresponding parent faces.  State halos are
    # subsequently replaced by the time-varying prescribed-halo forcing.
    create_like(
        parent_restart,
        INPUTS / "nest_initial_2017-01-01.nc",
        {"x": slice(X0 - 1, X1 + 1), "y": slice(Y0 - 1, Y1 + 1)},
    )

    with netCDF4.Dataset(INPUTS / "nest_kmt.nc") as ds:
        wet = np.asarray(ds.variables["kmt"][:] > 0)
        if not np.all(wet):
            raise RuntimeError("Selected nested test region is not completely wet")

    print(f"Created parent {NX_PARENT}x{NY_PARENT} and nest {X1-X0}x{Y1-Y0} inputs in {INPUTS}")


if __name__ == "__main__":
    main()
