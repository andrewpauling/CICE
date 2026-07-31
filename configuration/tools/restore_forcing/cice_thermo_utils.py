#!/usr/bin/env python3

from __future__ import annotations

import numpy as np


N_ICE_LAYERS = 7
MIN_ACTIVE_AREA = 1.0e-6
TMIN = -100.0

RHO_SNOW = 330.0
RHO_ICE = 917.0
RHO_WATER = 1026.0
CP_ICE = 2106.0
CP_OCN = 4218.0
L_FRESH = 3.34e5
PUNY = 1.0e-11

AZ1_LIQ = -18.48
BZ1_LIQ = 0.0
AZ2_LIQ = -10.3085
BZ2_LIQ = 62.4
TB_LIQ = -7.6362968855167352
SB_LIQ = 123.66702800276086
AZ1P_LIQ = AZ1_LIQ / 1000.0
BZ1P_LIQ = BZ1_LIQ / 1000.0
AZ2P_LIQ = AZ2_LIQ / 1000.0
BZ2P_LIQ = BZ2_LIQ / 1000.0


def snow_enthalpy_from_temperature(snow_temperature_c: np.ndarray) -> np.ndarray:
    return -RHO_SNOW * (L_FRESH - CP_ICE * snow_temperature_c)


def liquidus_brine_salinity_mush(ice_temperature_c: np.ndarray) -> np.ndarray:
    j1_liq = BZ1_LIQ / AZ1_LIQ
    k1_liq = 1.0 / 1000.0
    l1_liq = (1.0 + BZ1P_LIQ) / AZ1_LIQ
    j2_liq = BZ2_LIQ / AZ2_LIQ
    k2_liq = 1.0 / 1000.0
    l2_liq = (1.0 + BZ2P_LIQ) / AZ2_LIQ

    t_high = (ice_temperature_c > TB_LIQ).astype(np.float64)
    subzero = (ice_temperature_c <= 0.0).astype(np.float64)

    salinity = ((ice_temperature_c + j1_liq) / (k1_liq * ice_temperature_c + l1_liq)) * t_high
    salinity += ((ice_temperature_c + j2_liq) / (k2_liq * ice_temperature_c + l2_liq)) * (1.0 - t_high)
    return salinity * subzero


def mushy_liquid_fraction(ice_temperature_c: np.ndarray, bulk_salinity: np.ndarray) -> np.ndarray:
    brine_salinity = np.maximum(liquidus_brine_salinity_mush(ice_temperature_c), PUNY)
    return bulk_salinity / np.maximum(brine_salinity, bulk_salinity)


def ice_enthalpy_mush(ice_temperature_c: np.ndarray, bulk_salinity: np.ndarray) -> np.ndarray:
    liquid_fraction = mushy_liquid_fraction(ice_temperature_c, bulk_salinity)
    return (
        liquid_fraction * (CP_OCN * RHO_WATER - CP_ICE * RHO_ICE) * ice_temperature_c
        + RHO_ICE * CP_ICE * ice_temperature_c
        - (1.0 - liquid_fraction) * RHO_ICE * L_FRESH
    )


def liquidus_temperature_mush(brine_salinity: np.ndarray) -> np.ndarray:
    m1_liq = AZ1_LIQ
    n1_liq = -AZ1P_LIQ
    o1_liq = -BZ1_LIQ / AZ1_LIQ
    m2_liq = AZ2_LIQ
    n2_liq = -AZ2P_LIQ
    o2_liq = -BZ2_LIQ / AZ2_LIQ

    t_high = (brine_salinity <= SB_LIQ).astype(np.float64)
    temperature = ((brine_salinity / (m1_liq + n1_liq * brine_salinity)) + o1_liq) * t_high
    temperature += ((brine_salinity / (m2_liq + n2_liq * brine_salinity)) + o2_liq) * (1.0 - t_high)
    return temperature


def mushy_temperature_from_enthalpy(ice_enthalpy: np.ndarray, bulk_salinity: np.ndarray) -> np.ndarray:
    as1_liq = AZ1P_LIQ * (RHO_WATER * CP_OCN - RHO_ICE * CP_ICE)
    ac1_liq = RHO_ICE * CP_ICE * AZ1_LIQ
    bs1_liq = (1.0 + BZ1P_LIQ) * (RHO_WATER * CP_OCN - RHO_ICE * CP_ICE) + RHO_ICE * L_FRESH * AZ1P_LIQ
    bq1_liq = -AZ1_LIQ
    bc1_liq = RHO_ICE * CP_ICE * BZ1_LIQ - RHO_ICE * L_FRESH * AZ1_LIQ
    cs1_liq = RHO_ICE * L_FRESH * (1.0 + BZ1P_LIQ)
    cq1_liq = -BZ1_LIQ
    cc1_liq = -RHO_ICE * L_FRESH * BZ1_LIQ

    as2_liq = AZ2P_LIQ * (RHO_WATER * CP_OCN - RHO_ICE * CP_ICE)
    ac2_liq = RHO_ICE * CP_ICE * AZ2_LIQ
    bs2_liq = (1.0 + BZ2P_LIQ) * (RHO_WATER * CP_OCN - RHO_ICE * CP_ICE) + RHO_ICE * L_FRESH * AZ2P_LIQ
    bq2_liq = -AZ2_LIQ
    bc2_liq = RHO_ICE * CP_ICE * BZ2_LIQ - RHO_ICE * L_FRESH * AZ2_LIQ
    cs2_liq = RHO_ICE * L_FRESH * (1.0 + BZ2P_LIQ)
    cq2_liq = -BZ2_LIQ
    cc2_liq = -RHO_ICE * L_FRESH * BZ2_LIQ

    d_liq = ((1.0 + AZ1P_LIQ * TB_LIQ + BZ1P_LIQ) / (AZ1_LIQ * TB_LIQ + BZ1_LIQ)) * (
        (CP_OCN * RHO_WATER - CP_ICE * RHO_ICE) * TB_LIQ + L_FRESH * RHO_ICE
    )
    e_liq = CP_ICE * RHO_ICE * TB_LIQ - L_FRESH * RHO_ICE

    f1_liq = (-1000.0 * CP_OCN * RHO_WATER) / AZ1_LIQ
    g1_liq = -1000.0
    h1_liq = (-BZ1_LIQ * CP_OCN * RHO_WATER) / AZ1_LIQ
    f2_liq = (-1000.0 * CP_OCN * RHO_WATER) / AZ2_LIQ
    g2_liq = -1000.0
    h2_liq = (-BZ2_LIQ * CP_OCN * RHO_WATER) / AZ2_LIQ
    i_liq = 1.0 / (CP_OCN * RHO_WATER)

    s_low = (bulk_salinity < SB_LIQ).astype(np.float64)
    q0 = ((f1_liq * bulk_salinity) / (g1_liq + bulk_salinity) + h1_liq) * s_low
    q0 += ((f2_liq * bulk_salinity) / (g2_liq + bulk_salinity) + h2_liq) * (1.0 - s_low)
    q_melt = (ice_enthalpy > q0).astype(np.float64)

    qb = d_liq * bulk_salinity + e_liq
    t_high = (ice_enthalpy > qb).astype(np.float64)
    t_low = 1.0 - t_high

    a_term = (as1_liq * bulk_salinity + ac1_liq) * t_high + (as2_liq * bulk_salinity + ac2_liq) * t_low
    b_term = (bs1_liq * bulk_salinity + bq1_liq * ice_enthalpy + bc1_liq) * t_high
    b_term += (bs2_liq * bulk_salinity + bq2_liq * ice_enthalpy + bc2_liq) * t_low
    c_term = (cs1_liq * bulk_salinity + cq1_liq * ice_enthalpy + cc1_liq) * t_high
    c_term += (cs2_liq * bulk_salinity + cq2_liq * ice_enthalpy + cc2_liq) * t_low

    discriminant = np.maximum(b_term**2 - 4.0 * a_term * c_term, PUNY)
    temperature = (-b_term + np.sqrt(discriminant)) / (2.0 * a_term)
    return q_melt * ice_enthalpy * i_liq + (1.0 - q_melt) * temperature


def build_linear_ice_temperature_profile(
    surface_temperature_c: np.ndarray,
    freezing_temperature_c: np.ndarray,
    layer_count: int = N_ICE_LAYERS,
) -> np.ndarray:
    layers = []
    slope = freezing_temperature_c - surface_temperature_c
    for layer_index in range(layer_count):
        fraction = (layer_index + 0.5) / layer_count
        layers.append(surface_temperature_c + slope * fraction)
    return np.stack(layers, axis=0)


def compute_qice_layers(
    surface_temperature_c: np.ndarray,
    freezing_temperature_c: np.ndarray,
    salinity_layers: np.ndarray,
    active_category_mask: np.ndarray,
) -> list[np.ndarray]:
    layer_count = salinity_layers.shape[1]
    qice_layers = np.zeros((layer_count,) + surface_temperature_c.shape, dtype=np.float64)
    for category in range(salinity_layers.shape[0]):
        category_temperature = build_linear_ice_temperature_profile(
            surface_temperature_c[category],
            freezing_temperature_c,
            layer_count=layer_count,
        )
        category_active_mask = active_category_mask[category]
        category_salinity = salinity_layers[category][:, category_active_mask]
        if category_salinity.shape != category_temperature[:, category_active_mask].shape:
            category_salinity = category_salinity.T
        qice_layers[:, category, category_active_mask] = ice_enthalpy_mush(
            category_temperature[:, category_active_mask],
            category_salinity,
        )
    return [qice_layers[layer_index] for layer_index in range(layer_count)]


def validate_thermo_state(
    qice_layers: np.ndarray,
    salinity_layers: np.ndarray,
    active_category_mask: np.ndarray,
    min_temperature: float = TMIN,
) -> None:
    if qice_layers.shape != (salinity_layers.shape[1], salinity_layers.shape[0], salinity_layers.shape[2], salinity_layers.shape[3]):
        raise ValueError("Unexpected thermo-array shape mismatch.")

    salinity_for_inversion = np.transpose(salinity_layers, (1, 0, 2, 3))
    reconstructed_temperature = mushy_temperature_from_enthalpy(qice_layers, salinity_for_inversion)
    bad_mask = active_category_mask[None, :, :, :] & (reconstructed_temperature < min_temperature)
    if np.any(bad_mask):
        first_bad = np.argwhere(bad_mask)[0]
        layer_index, category, y_index, x_index = [int(value) for value in first_bad]
        raise ValueError(
            "Reconstructed ice temperature below Tmin after thermo generation: "
            f"layer={layer_index + 1}, category={category + 1}, y={y_index + 1}, x={x_index + 1}, "
            f"temperature={reconstructed_temperature[layer_index, category, y_index, x_index]:.6f}"
        )
