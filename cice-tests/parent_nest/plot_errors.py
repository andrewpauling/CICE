#!/usr/bin/env python3

"""Render spatial parent-to-nest error evolution without cartographic dependencies."""

from __future__ import annotations

from datetime import date, timedelta
from pathlib import Path

import netCDF4
import numpy as np
from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
OUTPUT = HERE / "analysis" / "figures"
X0, X1 = 100, 184
Y0, Y1 = 140, 220
START = date(2017, 1, 1)
DAYS = [START + timedelta(days=offset) for offset in range(1, 6)]

PANEL_W = 380
PANEL_H = 390
MAP_W = (X1 - X0) * 4
MAP_H = (Y1 - Y0) * 4
MAP_X = 38
MAP_Y = 50


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "LiberationSans-Bold.ttf" if bold else "LiberationSans-Regular.ttf"
    return ImageFont.truetype(f"/usr/share/fonts/liberation-sans/{name}", size=size)


FONT_SMALL = font(13)
FONT_TICK = font(12)
FONT_PANEL = font(16, bold=True)
FONT_ROW = font(19, bold=True)
FONT_TITLE = font(25, bold=True)
FONT_SUBTITLE = font(15)


def restart(case: str, day: date) -> Path:
    candidates = sorted((RUNS / case / "restart").glob(f"iced.{day.isoformat()}-*.nc"))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one restart for {case} on {day}, found {candidates}")
    return candidates[0]


def aice(path: Path, parent: bool) -> np.ndarray:
    with netCDF4.Dataset(path) as ds:
        values = np.asarray(ds.variables["aicen"][:], dtype=np.float64).sum(axis=0)
    if parent:
        return values[Y0:Y1, X0:X1]
    expected = (Y1 - Y0, X1 - X0)
    if values.shape == (expected[0] + 2, expected[1] + 2):
        return values[1:-1, 1:-1]
    if values.shape != expected:
        raise RuntimeError(f"Unexpected nested restart shape {values.shape}")
    return values


def errors() -> dict[str, list[np.ndarray]]:
    result: dict[str, list[np.ndarray]] = {}
    for scheme in ("remap", "upwind"):
        result[scheme] = []
        for day in DAYS:
            parent = aice(restart(f"pn_parent_{scheme}", day), parent=True)
            nest = aice(restart(f"pn_nest_{scheme}", day), parent=False)
            result[scheme].append(nest - parent)
    return result


def diverging_rgb(values: np.ndarray, limit: float) -> np.ndarray:
    """Blue-white-red map with white fixed at zero."""
    scaled = np.clip(values / limit, -1.0, 1.0)
    white = np.asarray([247.0, 247.0, 247.0])
    blue = np.asarray([33.0, 102.0, 172.0])
    red = np.asarray([178.0, 24.0, 43.0])
    magnitude = np.abs(scaled)[..., None]
    endpoint = np.where((scaled >= 0.0)[..., None], red, blue)
    return np.rint(white + magnitude * (endpoint - white)).astype(np.uint8)


def magnitude_rgb(error: np.ndarray, log_min: float = -8.0, log_max: float = 0.0) -> np.ndarray:
    """Sequential map of log10 absolute error, with exact zeros at the floor."""
    logged = np.full(error.shape, log_min, dtype=np.float64)
    nonzero = np.abs(error) > 0.0
    logged[nonzero] = np.log10(np.abs(error[nonzero]))
    scaled = np.clip((logged - log_min) / (log_max - log_min), 0.0, 1.0)
    positions = np.asarray([0.0, 0.32, 0.65, 1.0])
    colors = np.asarray([
        [247.0, 247.0, 247.0],
        [116.0, 173.0, 209.0],
        [254.0, 224.0, 139.0],
        [178.0, 24.0, 43.0],
    ])
    output = np.empty(error.shape + (3,), dtype=np.float64)
    for channel in range(3):
        output[..., channel] = np.interp(scaled, positions, colors[:, channel])
    return np.rint(output).astype(np.uint8)


def text_center(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, text_font: ImageFont.FreeTypeFont, fill=(20, 20, 20)) -> None:
    box = draw.textbbox((0, 0), value, font=text_font)
    width = box[2] - box[0]
    draw.text((xy[0] - width // 2, xy[1]), value, font=text_font, fill=fill)


def panel(
    error: np.ndarray,
    day: date,
    limit: float,
    show_x: bool,
    show_y: bool,
    title: str | None = None,
) -> Image.Image:
    canvas = Image.new("RGB", (PANEL_W, PANEL_H), "white")
    draw = ImageDraw.Draw(canvas)
    image = Image.fromarray(diverging_rgb(np.flipud(error), limit), mode="RGB")
    image = image.resize((MAP_W, MAP_H), resample=Image.Resampling.NEAREST)
    canvas.paste(image, (MAP_X, MAP_Y))
    draw.rectangle((MAP_X, MAP_Y, MAP_X + MAP_W - 1, MAP_Y + MAP_H - 1), outline=(25, 25, 25), width=1)

    rms = float(np.sqrt(np.mean(error * error)))
    max_abs = float(np.max(np.abs(error)))
    text_center(draw, (MAP_X + MAP_W // 2, 4), title or day.strftime("%d Jan 2017"), FONT_PANEL)
    text_center(draw, (MAP_X + MAP_W // 2, 27), f"RMS {rms:.3g}   max |e| {max_abs:.3g}", FONT_SMALL, (55, 55, 55))

    if show_x:
        for value in (0, 21, 42, 63, 83):
            x = MAP_X + round(value / 83 * (MAP_W - 1))
            draw.line((x, MAP_Y + MAP_H, x, MAP_Y + MAP_H + 4), fill=(30, 30, 30), width=1)
            text_center(draw, (x, MAP_Y + MAP_H + 5), str(value), FONT_TICK)
    if show_y:
        for value in (0, 20, 40, 60, 79):
            y = MAP_Y + MAP_H - 1 - round(value / 79 * (MAP_H - 1))
            draw.line((MAP_X - 4, y, MAP_X, y), fill=(30, 30, 30), width=1)
            label = str(value)
            box = draw.textbbox((0, 0), label, font=FONT_TICK)
            draw.text((MAP_X - 8 - (box[2] - box[0]), y - 7), label, font=FONT_TICK, fill=(20, 20, 20))
    return canvas


def magnitude_panel(error: np.ndarray, day: date, show_x: bool, show_y: bool) -> Image.Image:
    canvas = Image.new("RGB", (PANEL_W, PANEL_H), "white")
    draw = ImageDraw.Draw(canvas)
    image = Image.fromarray(magnitude_rgb(np.flipud(error)), mode="RGB")
    image = image.resize((MAP_W, MAP_H), resample=Image.Resampling.NEAREST)
    canvas.paste(image, (MAP_X, MAP_Y))
    draw.rectangle((MAP_X, MAP_Y, MAP_X + MAP_W - 1, MAP_Y + MAP_H - 1), outline=(25, 25, 25), width=1)
    rms = float(np.sqrt(np.mean(error * error)))
    max_abs = float(np.max(np.abs(error)))
    text_center(draw, (MAP_X + MAP_W // 2, 4), day.strftime("%d Jan 2017"), FONT_PANEL)
    text_center(draw, (MAP_X + MAP_W // 2, 27), f"RMS {rms:.3g}   max |e| {max_abs:.3g}", FONT_SMALL, (55, 55, 55))
    if show_x:
        for value in (0, 21, 42, 63, 83):
            x = MAP_X + round(value / 83 * (MAP_W - 1))
            draw.line((x, MAP_Y + MAP_H, x, MAP_Y + MAP_H + 4), fill=(30, 30, 30), width=1)
            text_center(draw, (x, MAP_Y + MAP_H + 5), str(value), FONT_TICK)
    if show_y:
        for value in (0, 20, 40, 60, 79):
            y = MAP_Y + MAP_H - 1 - round(value / 79 * (MAP_H - 1))
            draw.line((MAP_X - 4, y, MAP_X, y), fill=(30, 30, 30), width=1)
            label = str(value)
            box = draw.textbbox((0, 0), label, font=FONT_TICK)
            draw.text((MAP_X - 8 - (box[2] - box[0]), y - 7), label, font=FONT_TICK, fill=(20, 20, 20))
    return canvas


def colorbar(canvas: Image.Image, x: int, y: int, height: int, limit: float, label: str) -> None:
    draw = ImageDraw.Draw(canvas)
    gradient = np.linspace(limit, -limit, height)[:, None]
    bar = Image.fromarray(diverging_rgb(gradient, limit), mode="RGB").resize((24, height))
    canvas.paste(bar, (x, y))
    draw.rectangle((x, y, x + 23, y + height - 1), outline=(30, 30, 30), width=1)
    for fraction, value in ((0.0, limit), (0.25, limit / 2), (0.5, 0.0), (0.75, -limit / 2), (1.0, -limit)):
        yy = y + round(fraction * (height - 1))
        draw.line((x + 24, yy, x + 29, yy), fill=(30, 30, 30), width=1)
        draw.text((x + 33, yy - 7), f"{value:.3g}", font=FONT_TICK, fill=(20, 20, 20))
    rotated = Image.new("RGBA", (height, 28), (255, 255, 255, 0))
    rdraw = ImageDraw.Draw(rotated)
    text_center(rdraw, (height // 2, 2), label, FONT_SMALL)
    rotated = rotated.rotate(90, expand=True)
    canvas.paste(rotated, (x + 67, y), rotated)


def magnitude_colorbar(canvas: Image.Image, x: int, y: int, height: int) -> None:
    draw = ImageDraw.Draw(canvas)
    samples = np.power(10.0, np.linspace(0.0, -8.0, height))[:, None]
    bar = Image.fromarray(magnitude_rgb(samples), mode="RGB").resize((24, height))
    canvas.paste(bar, (x, y))
    draw.rectangle((x, y, x + 23, y + height - 1), outline=(30, 30, 30), width=1)
    for exponent in (0, -2, -4, -6, -8):
        fraction = -exponent / 8.0
        yy = y + round(fraction * (height - 1))
        draw.line((x + 24, yy, x + 29, yy), fill=(30, 30, 30), width=1)
        draw.text((x + 33, yy - 7), f"1e{exponent}", font=FONT_TICK, fill=(20, 20, 20))
    rotated = Image.new("RGBA", (height, 28), (255, 255, 255, 0))
    rdraw = ImageDraw.Draw(rotated)
    text_center(rdraw, (height // 2, 2), "absolute AICE error", FONT_SMALL)
    rotated = rotated.rotate(90, expand=True)
    canvas.paste(rotated, (x + 67, y), rotated)


def row_label(canvas: Image.Image, name: str, top: int) -> None:
    label = Image.new("RGBA", (PANEL_H, 38), (255, 255, 255, 0))
    draw = ImageDraw.Draw(label)
    text_center(draw, (PANEL_H // 2, 5), name, FONT_ROW)
    label = label.rotate(90, expand=True)
    canvas.paste(label, (4, top), label)


def grid_figure(data: dict[str, list[np.ndarray]], limits: dict[str, float], output: Path, title: str) -> None:
    left = 52
    top = 72
    gap = 4
    right = 115
    width = left + 5 * PANEL_W + 4 * gap + right
    height = top + 2 * PANEL_H + 52
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    text_center(draw, (width // 2, 12), title, FONT_TITLE)
    text_center(draw, (width // 2, 44), "signed error = nested AICE - matching parent AICE", FONT_SUBTITLE, (55, 55, 55))

    for row, scheme in enumerate(("remap", "upwind")):
        y = top + row * PANEL_H
        row_label(canvas, scheme.upper(), y)
        for column, (day, error) in enumerate(zip(DAYS, data[scheme])):
            rendered = panel(error, day, limits[scheme], show_x=(row == 1), show_y=(column == 0))
            canvas.paste(rendered, (left + column * (PANEL_W + gap), y))
        colorbar(canvas, left + 5 * PANEL_W + 4 * gap + 13, y + MAP_Y, MAP_H, limits[scheme], "AICE error")

    text_center(draw, (left + (5 * PANEL_W + 4 * gap) // 2, height - 30), "nested-grid x index (west to east); y increases south to north", FONT_SUBTITLE)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, optimize=True)


def difference_figure(data: dict[str, list[np.ndarray]], limit: float, output: Path) -> None:
    differences = [remap - upwind for remap, upwind in zip(data["remap"], data["upwind"])]
    left = 52
    top = 72
    gap = 4
    right = 115
    width = left + 5 * PANEL_W + 4 * gap + right
    height = top + PANEL_H + 52
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    text_center(draw, (width // 2, 12), "Difference between remap and upwind parent-to-nest errors", FONT_TITLE)
    text_center(draw, (width // 2, 44), "(remap nest - remap parent) - (upwind nest - upwind parent)", FONT_SUBTITLE, (55, 55, 55))
    row_label(canvas, "REMAP - UPWIND", top)
    for column, (day, error) in enumerate(zip(DAYS, differences)):
        rendered = panel(error, day, limit, show_x=True, show_y=(column == 0))
        canvas.paste(rendered, (left + column * (PANEL_W + gap), top))
    colorbar(canvas, left + 5 * PANEL_W + 4 * gap + 13, top + MAP_Y, MAP_H, limit, "difference in AICE error")
    text_center(draw, (left + (5 * PANEL_W + 4 * gap) // 2, height - 30), "nested-grid x index (west to east); y increases south to north", FONT_SUBTITLE)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, optimize=True)


def magnitude_grid_figure(data: dict[str, list[np.ndarray]], output: Path) -> None:
    left = 52
    top = 72
    gap = 4
    right = 115
    width = left + 5 * PANEL_W + 4 * gap + right
    height = top + 2 * PANEL_H + 52
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    text_center(draw, (width // 2, 12), "Magnitude and inward growth of AICE parent-to-nest error", FONT_TITLE)
    text_center(draw, (width // 2, 44), "common logarithmic scale; white <= 1e-8; west and south are inflow boundaries", FONT_SUBTITLE, (55, 55, 55))
    for row, scheme in enumerate(("remap", "upwind")):
        y = top + row * PANEL_H
        row_label(canvas, scheme.upper(), y)
        for column, (day, error) in enumerate(zip(DAYS, data[scheme])):
            rendered = magnitude_panel(error, day, show_x=(row == 1), show_y=(column == 0))
            canvas.paste(rendered, (left + column * (PANEL_W + gap), y))
        magnitude_colorbar(canvas, left + 5 * PANEL_W + 4 * gap + 13, y + MAP_Y, MAP_H)
    text_center(draw, (left + (5 * PANEL_W + 4 * gap) // 2, height - 30), "nested-grid x index (west to east); y increases south to north", FONT_SUBTITLE)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, optimize=True)


def animation(data: dict[str, list[np.ndarray]], limits: dict[str, float], output: Path) -> None:
    frames = []
    panel_x = {"remap": 42, "upwind": 526}
    width = 1020
    height = PANEL_H + 92
    for index, day in enumerate(DAYS):
        canvas = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(canvas)
        text_center(draw, (width // 2, 10), f"AICE parent-to-nest error: {day.strftime('%d %B %Y')}", FONT_TITLE)
        for scheme in ("remap", "upwind"):
            x = panel_x[scheme]
            rendered = panel(
                data[scheme][index], day, limits[scheme], show_x=True, show_y=True,
                title=f"{scheme.upper()} (scale +/-{limits[scheme]:.3g})",
            )
            canvas.paste(rendered, (x, 52))
            colorbar(canvas, x + PANEL_W + 4, 52 + MAP_Y, MAP_H, limits[scheme], f"{scheme} AICE error")
        frames.append(canvas)
    output.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(output, save_all=True, append_images=frames[1:], duration=1200, loop=0, optimize=True)


def main() -> None:
    data = errors()
    common_limit = 0.70
    scheme_limits = {"remap": 0.70, "upwind": 0.035}
    difference_limit = 0.70
    grid_figure(
        data,
        {"remap": common_limit, "upwind": common_limit},
        OUTPUT / "aice_error_evolution_common_scale.png",
        "AICE parent-to-nest error evolution: common color scale",
    )
    grid_figure(
        data,
        scheme_limits,
        OUTPUT / "aice_error_evolution_scheme_scales.png",
        "AICE parent-to-nest error evolution: detail at each scheme's scale",
    )
    difference_figure(data, difference_limit, OUTPUT / "aice_error_remap_minus_upwind.png")
    magnitude_grid_figure(data, OUTPUT / "aice_error_magnitude_logscale.png")
    animation(data, scheme_limits, OUTPUT / "aice_error_evolution.gif")
    print(f"Wrote figures to {OUTPUT}")


if __name__ == "__main__":
    main()
