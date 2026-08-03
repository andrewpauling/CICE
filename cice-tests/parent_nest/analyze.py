#!/usr/bin/env python3

from __future__ import annotations

from datetime import date
import json
from pathlib import Path

import netCDF4
import numpy as np


HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
INPUTS = HERE / "inputs"
ANALYSIS = HERE / "analysis"
X0, X1 = 100, 184
Y0, Y1 = 140, 220
INITIAL_DAY = date(2017, 1, 1)
FINAL_DAY = date(2017, 1, 6)


def one_restart(case: str, day: date = FINAL_DAY) -> Path:
    paths = sorted((RUNS / case / "restart").glob(f"iced.{day.isoformat()}-*.nc"))
    if len(paths) != 1:
        raise RuntimeError(f"Expected one restart for {case} on {day}, found {paths}")
    return paths[0]


def state(path: Path, subset: bool) -> dict[str, np.ndarray]:
    selection = (slice(None), slice(Y0, Y1), slice(X0, X1)) if subset else np.s_[:, :, :]
    with netCDF4.Dataset(path) as ds:
        category = {name: np.asarray(ds.variables[name][selection], dtype=np.float64) for name in ("aicen", "vicen", "vsnon")}
    if not subset:
        # Nested runs use restart_ext so fixed face velocities are available
        # on the one-cell exterior stencil in this transport-only test.
        expected = (Y1 - Y0, X1 - X0)
        shape = category["aicen"].shape[-2:]
        if shape == (expected[0] + 2, expected[1] + 2):
            category = {name: values[:, 1:-1, 1:-1] for name, values in category.items()}
        elif shape != expected:
            raise RuntimeError(f"Unexpected nested restart shape {shape}; expected {expected} or extended grid")
    return {
        "aice": category["aicen"].sum(axis=0),
        "vice": category["vicen"].sum(axis=0),
        "vsno": category["vsnon"].sum(axis=0),
        "aicen": category["aicen"],
        "vicen": category["vicen"],
        "vsnon": category["vsnon"],
    }


def area_weights() -> np.ndarray:
    with netCDF4.Dataset(INPUTS / "nest_grid.nc") as ds:
        hte = np.asarray(ds.variables["hte"][:], dtype=np.float64) / 100.0
        htn = np.asarray(ds.variables["htn"][:], dtype=np.float64) / 100.0
    return hte * htn


def metrics(error: np.ndarray) -> dict[str, float]:
    return {
        "max_abs": float(np.max(np.abs(error))),
        "mean_abs": float(np.mean(np.abs(error))),
        "rms": float(np.sqrt(np.mean(error * error))),
        "bias": float(np.mean(error)),
    }


def distance_profiles(error: np.ndarray) -> dict[str, dict[str, float]]:
    ny, nx = error.shape
    yy, xx = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    distance = np.minimum.reduce((xx, nx - 1 - xx, yy, ny - 1 - yy))
    profile = {}
    for ring in range(0, 21):
        values = error[distance == ring]
        if values.size:
            profile[str(ring)] = metrics(values)
    values = error[distance >= 20]
    if values.size:
        profile[">=20"] = metrics(values)
    return profile


def edge_profiles(error: np.ndarray) -> dict[str, list[dict[str, float]]]:
    result = {}
    for edge in ("west", "east", "south", "north"):
        rows = []
        for distance in range(0, 16):
            if edge == "west":
                values = error[:, distance]
            elif edge == "east":
                values = error[:, -(distance + 1)]
            elif edge == "south":
                values = error[distance, :]
            else:
                values = error[-(distance + 1), :]
            rows.append({"distance": distance, **metrics(values)})
        result[edge] = rows
    return result


def edge_core_metrics(error: np.ndarray, corner_trim: int = 10) -> dict[str, object]:
    """Metrics on the outermost row/column away from corner interactions."""
    return {
        "corner_trim": corner_trim,
        "west": metrics(error[corner_trim:-corner_trim, 0]),
        "east": metrics(error[corner_trim:-corner_trim, -1]),
        "south": metrics(error[0, corner_trim:-corner_trim]),
        "north": metrics(error[-1, corner_trim:-corner_trim]),
    }


def time_series(advection: str, weights: np.ndarray) -> list[dict[str, object]]:
    rows = []
    for offset in range(1, (FINAL_DAY - INITIAL_DAY).days + 1):
        day = date.fromordinal(INITIAL_DAY.toordinal() + offset)
        parent = state(one_restart(f"pn_parent_{advection}", day), subset=True)["aice"]
        nest = state(one_restart(f"pn_nest_{advection}", day), subset=False)["aice"]
        error = nest - parent
        parent_integral = float(np.sum(parent * weights))
        nest_integral = float(np.sum(nest * weights))
        rows.append({
            "day": day.isoformat(),
            "global": metrics(error),
            "relative_integral_error": (nest_integral - parent_integral) / parent_integral,
            "edge_core": edge_core_metrics(error),
        })
    return rows


def compare(advection: str, weights: np.ndarray) -> dict[str, object]:
    parent_path = one_restart(f"pn_parent_{advection}")
    nest_path = one_restart(f"pn_nest_{advection}")
    parent = state(parent_path, subset=True)
    nest = state(nest_path, subset=False)
    output: dict[str, object] = {
        "parent_restart": str(parent_path),
        "nest_restart": str(nest_path),
        "fields": {},
    }
    for field in ("aice", "vice", "vsno"):
        error = nest[field] - parent[field]
        parent_integral = float(np.sum(parent[field] * weights))
        nest_integral = float(np.sum(nest[field] * weights))
        output["fields"][field] = {
            "global": metrics(error),
            "distance_from_nearest_edge": distance_profiles(error),
            "by_edge": edge_profiles(error),
            "edge_core": edge_core_metrics(error),
            "parent_integral": parent_integral,
            "nest_integral": nest_integral,
            "relative_integral_error": (nest_integral - parent_integral) / parent_integral,
        }

    threshold = 0.15
    parent_edge = parent["aice"] >= threshold
    nest_edge = nest["aice"] >= threshold
    output["aice_threshold"] = threshold
    output["aice_threshold_disagreement_fraction"] = float(np.mean(parent_edge != nest_edge))
    output["aice_time_series"] = time_series(advection, weights)
    return output


def text_summary(summary: dict[str, object]) -> str:
    lines = ["P-SKRIPS standalone parent-to-nest transport comparison", ""]
    for advection in ("remap", "upwind"):
        result = summary[advection]
        lines.append(advection.upper())
        for field in ("aice", "vice", "vsno"):
            data = result["fields"][field]
            global_data = data["global"]
            interior = data["distance_from_nearest_edge"].get(">=20", {})
            lines.append(
                f"  {field}: max={global_data['max_abs']:.8e} "
                f"rms={global_data['rms']:.8e} "
                f"interior_rms={interior.get('rms', float('nan')):.8e} "
                f"integral_rel={data['relative_integral_error']:.8e}"
            )
        lines.append(
            f"  aice>=0.15 mask disagreement: "
            f"{result['aice_threshold_disagreement_fraction']:.8e}"
        )
        edge = result["fields"]["aice"]["edge_core"]
        lines.append(
            "  aice outer-ring core RMS (W/E/S/N): "
            + "/".join(f"{edge[name]['rms']:.8e}" for name in ("west", "east", "south", "north"))
        )
        lines.append(
            "  daily aice RMS: "
            + ", ".join(f"{row['global']['rms']:.8e}" for row in result["aice_time_series"])
        )
        lines.append("")

    remap = summary["remap"]["fields"]["aice"]["global"]["rms"]
    upwind = summary["upwind"]["fields"]["aice"]["global"]["rms"]
    lines.append(f"AICE RMS remap/upwind ratio: {remap / upwind if upwind else float('inf'):.8e}")
    return "\n".join(lines) + "\n"


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    weights = area_weights()
    summary = {advection: compare(advection, weights) for advection in ("remap", "upwind")}
    (ANALYSIS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    report = text_summary(summary)
    (ANALYSIS / "summary.txt").write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
