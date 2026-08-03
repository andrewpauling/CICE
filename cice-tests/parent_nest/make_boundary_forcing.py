#!/usr/bin/env python3

from __future__ import annotations

import argparse
from datetime import date, timedelta
from pathlib import Path

import netCDF4
import numpy as np


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
RUNS = HERE / "runs"
X0, X1 = 100, 184
Y0, Y1 = 140, 220


def restart_paths(advection: str) -> list[Path]:
    paths = [INPUTS / "parent_initial_2017-01-01.nc"]
    restart_dir = RUNS / f"pn_parent_{advection}" / "restart"
    for offset in range(1, 6):
        day = date(2017, 1, 1) + timedelta(days=offset)
        candidates = sorted(restart_dir.glob(f"iced.{day.isoformat()}-*.nc"))
        if len(candidates) != 1:
            raise RuntimeError(f"Expected one restart for {day} in {restart_dir}, found {candidates}")
        paths.append(candidates[0])
    return paths


def make_file(advection: str) -> Path:
    paths = restart_paths(advection)
    output = INPUTS / f"parent_{advection}_nest_boundary.nc"
    with netCDF4.Dataset(paths[0]) as first:
        ncat = len(first.dimensions["ncat"])

    with netCDF4.Dataset(output, "w", format="NETCDF4_CLASSIC") as dst:
        dst.createDimension("time", len(paths))
        dst.createDimension("nc", ncat)
        dst.createDimension("nj", Y1 - Y0)
        dst.createDimension("ni", X1 - X0)

        time = dst.createVariable("time", "f8", ("time",))
        time[:] = np.arange(len(paths), dtype=np.float64)
        time.units = "days since 2017-01-01 00:00:00"
        time.calendar = "proleptic_gregorian"
        ncat_var = dst.createVariable("NCAT", "f4", ("nc",))
        ncat_var[:] = np.arange(1, ncat + 1, dtype=np.float32)

        outputs = {}
        for name, units in (("aicen", "1"), ("vicen", "m"), ("vsnon", "m")):
            variable = dst.createVariable(name, "f8", ("time", "nc", "nj", "ni"), zlib=True, complevel=1)
            variable.units = units
            outputs[name] = variable

        for record, path in enumerate(paths):
            with netCDF4.Dataset(path) as src:
                for name, variable in outputs.items():
                    parent = np.asarray(src.variables[name][:], dtype=np.float64)
                    target = parent[:, Y0:Y1, X0:X1].copy()

                    # structure_only forcing does not restore physical-edge
                    # cells.  CICE uses each edge value as the target for the
                    # adjacent exterior halo, so populate those file edges
                    # from the corresponding parent cells just outside the
                    # nested domain.  This makes the one-cell transport
                    # stencil comparable to the same region in the parent.
                    target[:, :, 0] = parent[:, Y0:Y1, X0 - 1]
                    target[:, :, -1] = parent[:, Y0:Y1, X1]
                    target[:, 0, :] = parent[:, Y0 - 1, X0:X1]
                    target[:, -1, :] = parent[:, Y1, X0:X1]
                    variable[record, :, :, :] = target

        dst.restore_source = "parent_cice_restart_exterior_edge_cells"
        dst.restore_format = "structure_only"
        dst.edge_semantics = "file physical edges contain adjacent parent exterior cells for prescribed halos"
        dst.parent_advection = advection
        dst.parent_x_slice = f"{X0}:{X1}"
        dst.parent_y_slice = f"{Y0}:{Y1}"
        dst.parent_restart_files = "\n".join(str(path) for path in paths)
    print(f"Created {output}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("advection", nargs="*")
    args = parser.parse_args()
    for advection in (args.advection or ("remap", "upwind")):
        if advection not in ("remap", "upwind"):
            parser.error("advection must be remap or upwind")
        make_file(advection)


if __name__ == "__main__":
    main()
