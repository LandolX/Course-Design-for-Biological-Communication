import numpy as np
from scipy.constants import c
from typing import Dict, Tuple, Optional
from dataclasses import dataclass


TISSUE_PROPERTIES_CN: Dict[str, Dict[str, Tuple[float, float, float]]] = {
    "\u8102\u80aa (Fat)": {
        "MICS": (5.28, 0.20, 0.45),
        "ISM":  (5.28, 0.14, 1.12),
    },
    "\u808c\u8089 (Muscle)": {
        "MICS": (52.73, 0.74, 2.50),
        "ISM":  (52.73, 1.74, 4.80),
    },
    "\u76ae\u80a4 (Skin)": {
        "MICS": (38.00, 0.46, 1.80),
        "ISM":  (38.00, 1.46, 4.20),
    },
    "\u9aa8\u9abc (Bone)": {
        "MICS": (11.38, 0.09, 1.20),
        "ISM":  (11.38, 0.39, 3.10),
    },
    "\u8840\u6db2 (Blood)": {
        "MICS": (58.30, 1.38, 2.80),
        "ISM":  (58.30, 2.54, 5.60),
    },
}

FREQ_BANDS = {
    "MICS": (402e6, 405e6),
    "ISM":  (2.4e9, 2.48e9),
}

TISSUE_PROPERTIES = {
    "skin": {
        "relative_permittivity": {100e6: 46.7, 433e6: 46.0, 915e6: 41.0, 2.4e9: 38.0, 3e9: 37.0},
        "conductivity": {100e6: 0.69, 433e6: 0.70, 915e6: 0.87, 2.4e9: 1.46, 3e9: 1.60},
    },
    "fat": {
        "relative_permittivity": {100e6: 5.46, 433e6: 5.40, 915e6: 5.35, 2.4e9: 5.28, 3e9: 5.25},
        "conductivity": {100e6: 0.05, 433e6: 0.05, 915e6: 0.06, 2.4e9: 0.10, 3e9: 0.15},
    },
    "muscle": {
        "relative_permittivity": {100e6: 66.2, 433e6: 56.0, 915e6: 55.0, 2.4e9: 52.7, 3e9: 50.8},
        "conductivity": {100e6: 0.73, 433e6: 0.80, 915e6: 0.95, 2.4e9: 1.74, 3e9: 2.20},
    },
    "bone": {
        "relative_permittivity": {100e6: 12.6, 433e6: 11.5, 915e6: 11.0, 2.4e9: 10.5, 3e9: 10.0},
        "conductivity": {100e6: 0.08, 433e6: 0.10, 915e6: 0.15, 2.4e9: 0.30, 3e9: 0.40},
    },
}

TISSUE_THICKNESS = {
    "skin": 1.5e-3,
    "fat": 10.0e-3,
    "muscle": 20.0e-3,
    "bone": float("inf"),
}

TISSUE_LAYER_ORDER = ["skin", "fat", "muscle", "bone"]


def _interp_property(tissue: str, prop: str, freq: float) -> float:
    data = TISSUE_PROPERTIES[tissue][prop]
    freqs = sorted(data.keys())
    if freq <= freqs[0]:
        return data[freqs[0]]
    if freq >= freqs[-1]:
        return data[freqs[-1]]
    values = [data[f] for f in freqs]
    return float(np.interp(freq, freqs, values))


def get_permittivity(tissue: str, freq: float) -> float:
    return _interp_property(tissue, "relative_permittivity", freq)


def get_conductivity(tissue: str, freq: float) -> float:
    return _interp_property(tissue, "conductivity", freq)


def calc_attenuation_coefficient(tissue: str, freq: float) -> float:
    mu0 = 4 * np.pi * 1e-7
    eps0 = 8.854e-12
    eps_r = get_permittivity(tissue, freq)
    sigma = get_conductivity(tissue, freq)
    eps_c = eps_r - 1j * sigma / (2 * np.pi * freq * eps0)
    gamma = 1j * 2 * np.pi * freq * np.sqrt(mu0 * eps0 * eps_c)
    alpha = np.real(gamma)
    return alpha


def calc_tissue_path_loss(tissue: str, freq: float, depth: float) -> float:
    alpha = calc_attenuation_coefficient(tissue, freq)
    path_loss_tissue = 20 * np.log10(np.exp(alpha * depth))
    return path_loss_tissue


def calc_free_space_path_loss(freq: float, distance: float) -> float:
    wavelength = c / freq
    if distance <= 0 or wavelength <= 0:
        return 0.0
    return 20 * np.log10(4 * np.pi * distance / wavelength)


def calc_multilayer_path_loss(
    freq: float,
    depths: list[float],
    tissue_types: list[str],
    extra_distance: float = 0.01,
) -> float:
    pl_total = 0.0
    total_depth = sum(depths)
    for tissue, depth in zip(tissue_types, depths):
        if depth <= 0:
            continue
        pl_total += calc_tissue_path_loss(tissue, freq, depth)
    pl_total += calc_free_space_path_loss(freq, total_depth + extra_distance)
    return pl_total


def calc_snr(
    freq: float,
    depth: float,
    tissue_type: str = "muscle",
    tx_power_dbm: float = 0.0,
    noise_figure_db: float = 5.0,
    bandwidth_hz: float = 1e6,
) -> float:
    pl = calc_tissue_path_loss(tissue_type, freq, depth) + calc_free_space_path_loss(
        freq, depth + 0.01
    )
    thermal_noise_dbm = -174 + 10 * np.log10(bandwidth_hz)
    total_noise_dbm = thermal_noise_dbm + noise_figure_db
    snr_db = tx_power_dbm - pl - total_noise_dbm
    return snr_db


def estimate_path_loss_implant(
    freq: float,
    implant_depth_m: float,
    tissue_layers: list[str] | None = None,
    layer_thicknesses: list[float] | None = None,
) -> float:
    if tissue_layers is None:
        tissue_layers = TISSUE_LAYER_ORDER
    if layer_thicknesses is None:
        layer_thicknesses = [TISSUE_THICKNESS[t] for t in tissue_layers]
    remaining = implant_depth_m
    used_layers = []
    used_depths = []
    for tissue, thick in zip(tissue_layers, layer_thicknesses):
        if remaining <= 0:
            break
        d = min(thick if np.isfinite(thick) else remaining, remaining)
        used_layers.append(tissue)
        used_depths.append(d)
        remaining -= d
    if remaining > 0:
        used_layers.append("muscle")
        used_depths.append(remaining)
    return calc_multilayer_path_loss(freq, used_depths, used_layers)


@dataclass
class PathLossResult:
    distance: np.ndarray
    path_loss: np.ndarray
    shadowing: np.ndarray
    multipath_loss: np.ndarray
    total_loss: np.ndarray


class ImplantPathLossModel:
    def __init__(
        self,
        frequency: float = 403e6,
        band: str = "MICS",
        tissue_type: str = "\u808c\u8089 (Muscle)",
    ):
        self.frequency = frequency
        self.band = band
        self.tissue_type = tissue_type

        band_key = "MICS" if frequency < 1e9 else "ISM"
        props = TISSUE_PROPERTIES_CN[tissue_type][band_key]
        self.epsilon_r = props[0]
        self.sigma = props[1]
        self.attenuation_coeff = props[2]

        self._calc_wave_properties()

    def _calc_wave_properties(self):
        epsilon_0 = 8.854e-12
        mu_0 = 4 * np.pi * 1e-7
        omega = 2 * np.pi * self.frequency

        epsilon_c = epsilon_0 * (self.epsilon_r - 1j * self.sigma / (omega * epsilon_0))
        self.gamma = 1j * omega * np.sqrt(mu_0 * epsilon_c)
        self.alpha = np.real(self.gamma)
        self.beta = np.imag(self.gamma)
        self.wavelength = 2 * np.pi / self.beta

    def get_ref_path_loss(self, d0: float = 0.05) -> float:
        epsilon_0 = 8.854e-12
        mu_0 = 4 * np.pi * 1e-7
        eta_0 = np.sqrt(mu_0 / epsilon_0)
        eta_tissue = np.sqrt(mu_0 / (epsilon_0 * self.epsilon_r))
        gamma_t = (eta_tissue - eta_0) / (eta_tissue + eta_0)

        transmit_gain = 10 ** (-2.0 / 10)
        receive_gain = 10 ** (-2.0 / 10)

        pl0 = (
            -20 * np.log10(self.wavelength / (4 * np.pi * d0))
            - transmit_gain
            - receive_gain
            + 20 * np.log10(1 - gamma_t ** 2)
        )
        return pl0

    def compute_path_loss(
        self,
        distance: np.ndarray,
        n: float = 4.5,
        d0: float = 0.05,
        sigma: float = 3.5,
        include_shadowing: bool = True,
        include_multipath: bool = True,
    ) -> PathLossResult:
        distance = np.asarray(distance, dtype=np.float64)
        distance = np.maximum(distance, 1e-6)

        pl0 = self.get_ref_path_loss(d0)

        path_loss = pl0 + 10 * n * np.log10(np.maximum(distance / d0, 1e-6))

        tissue_factor = self.attenuation_coeff * 100
        tissue_loss = tissue_factor * distance
        path_loss = path_loss + tissue_loss

        distance_mask = distance >= d0
        path_loss = np.where(distance_mask, path_loss, pl0 * np.ones_like(distance))

        if include_shadowing:
            shadowing = np.random.normal(0, sigma, size=distance.shape)
        else:
            shadowing = np.zeros_like(distance)

        if include_multipath:
            multipath_loss = self._calc_multipath_effect(distance)
        else:
            multipath_loss = np.zeros_like(distance)

        total_loss = path_loss + shadowing + multipath_loss

        return PathLossResult(
            distance=distance,
            path_loss=path_loss,
            shadowing=shadowing,
            multipath_loss=multipath_loss,
            total_loss=total_loss,
        )

    def _calc_multipath_effect(self, distance: np.ndarray) -> np.ndarray:
        tissue_list = list(TISSUE_PROPERTIES_CN.keys())
        num_layers = len(tissue_list)

        reflection_coeffs = 0.3 * np.random.uniform(0.5, 1.0, num_layers)

        mp_loss = np.zeros_like(distance)
        for i, rc in enumerate(reflection_coeffs):
            path_diff = 2 * (i + 1) * 0.01
            phase_diff = 2 * np.pi * path_diff / self.wavelength
            mp_factor = np.sqrt(
                1 + rc ** 2 + 2 * rc * np.cos(phase_diff + 2 * np.pi * distance / self.wavelength)
            )
            mp_loss += -20 * np.log10(np.maximum(mp_factor, 1e-6))

        return mp_loss

    def calc_freq_dependent_attenuation(self, frequencies: np.ndarray) -> np.ndarray:
        f_ghz = np.asarray(frequencies, dtype=np.float64) / 1e9
        a0 = self.attenuation_coeff
        attenuation = a0 * (f_ghz / (self.frequency / 1e9)) ** 0.8
        return attenuation

    @staticmethod
    def compare_with_literature(
        distances: np.ndarray,
        freq: float = 403e6,
    ) -> Dict[str, np.ndarray]:
        model_muscle = ImplantPathLossModel(frequency=freq, band="MICS", tissue_type="\u808c\u8089 (Muscle)")
        model_fat = ImplantPathLossModel(frequency=freq, band="MICS", tissue_type="\u8102\u80aa (Fat)")
        model_skin = ImplantPathLossModel(frequency=freq, band="MICS", tissue_type="\u76ae\u80a4 (Skin)")

        result_muscle = model_muscle.compute_path_loss(
            distances, n=4.5, include_shadowing=False, include_multipath=False,
        )
        result_fat = model_fat.compute_path_loss(
            distances, n=4.5, include_shadowing=False, include_multipath=False,
        )
        result_skin = model_skin.compute_path_loss(
            distances, n=4.5, include_shadowing=False, include_multipath=False,
        )

        pl_ieee_deep = 47.14 + 10 * 4.26 * np.log10(np.maximum(distances / 0.05, 1e-6))
        pl_ieee_surface = 40.01 + 10 * 5.84 * np.log10(np.maximum(distances / 0.05, 1e-6))

        return {
            "distance": distances,
            "muscle_model": result_muscle.path_loss,
            "fat_model": result_fat.path_loss,
            "skin_model": result_skin.path_loss,
            "ieee_deep_implant": pl_ieee_deep,
            "ieee_surface_implant": pl_ieee_surface,
        }
