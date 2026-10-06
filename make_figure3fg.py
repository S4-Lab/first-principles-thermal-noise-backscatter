"""
Figure 3(f,g): mismatched high bandwidth case — measured vs numerical vs theoretical.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# ------------------------- configuration -------------------------

@dataclass(frozen=True)
class Config:
    seed: int = 42

    # system parameters
    c: float = 3e8
    k: float = 1.38e-23
    B: float = 982_000_000
    center_f: float = 11e9
    T_amb: float = 290
    T_fw: float = 99
    T_bw: float = 55
    T_s: float = 290
    int_time: float = 0.01
    phase_center_mm: float = 20.6
    electrical_length_m: float = 0.42
    ant_gain_dB: float = 28.37
    phi: float = np.pi / 4

    # plotting
    xlim_mm: tuple[float, float] = (90, 400)
    ylim_dbmhz: tuple[float, float] = (-172.87, -172.66)
    xticks_mm: tuple[float, ...] = (1e2, 2e2, 3e2, 4e2)
    yticks_dbmhz: tuple[float, ...] = (-172.70, -172.75, -172.80, -172.85)

    ylim_diff_W: tuple[float, float] = (0.0, 2e-13)
    yticks_diff_W: tuple[float, ...] = (0.0, 0.5e-13, 1.0e-13, 1.5e-13)


# ------------------------- style -------------------------

def set_plot_style():
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica"],
            "axes.labelsize": 25,
            "xtick.labelsize": 25,
            "ytick.labelsize": 25,
            "legend.fontsize": 23,
            "mathtext.default": "regular",
        }
    )


# ------------------------- signal/noise helpers -------------------------

def fractional_delay(signal: np.ndarray, delay_samples: float) -> np.ndarray:
    N = len(signal)
    S = np.fft.fft(signal)
    f = np.fft.fftfreq(N)
    return np.fft.ifft(S * np.exp(-2j * np.pi * f * delay_samples))


def generate_thermal_noise_voltage(
    rng: np.random.Generator,
    k: float,
    num_samples: int,
    bandwidth: float,
    temperature: float,
    resistance: float = 50,
) -> np.ndarray:
    sigma = np.sqrt(k * temperature * bandwidth * resistance)
    a = rng.rayleigh(scale=sigma / np.sqrt(2), size=num_samples)
    phase = rng.uniform(0, 2 * np.pi, size=num_samples)
    return a * np.exp(1j * phase)


def power_dbm_from_v(v: np.ndarray, resistance: float = 50) -> float:
    p_W = np.mean((np.abs(v) ** 2) / resistance)
    return 10 * np.log10(p_W * 1000)


# ------------------------- models -------------------------

def numerical_matched_dbm(
    rng: np.random.Generator,
    cfg: Config,
    ant_gain_lin: float,
    lam: float,
    r_m: float,
    num_samples: int,
) -> float:
    h = (np.sqrt(ant_gain_lin) * lam / (4 * np.pi * r_m)) * np.exp(1j * (-2 * np.pi * r_m / lam))

    v_s = generate_thermal_noise_voltage(rng, cfg.k, num_samples, cfg.B, cfg.T_s)
    v_s_rx = h * v_s

    v_amb = generate_thermal_noise_voltage(rng, cfg.k, num_samples, cfg.B, cfg.T_amb)
    v_amb_rx = v_amb * np.sqrt(1 - np.abs(h) ** 2)

    v_fw = generate_thermal_noise_voltage(rng, cfg.k, num_samples, cfg.B, cfg.T_fw)

    return power_dbm_from_v(v_s_rx + v_amb_rx + v_fw)


def theory_mismatched_power_W(cfg: Config, ant_gain_lin: float, lam: float, r_m: float, phi: float) -> float:
    tau = 2 * (r_m + cfg.electrical_length_m) / cfg.c
    loss = np.sqrt(ant_gain_lin) * lam / (4 * np.pi)

    first = cfg.T_amb + cfg.T_fw
    second = (cfg.T_bw - cfg.T_amb) * (loss**4) / (r_m**4)
    third = (2 * np.sqrt(cfg.T_bw * cfg.T_fw) * (loss**2) * np.sinc(cfg.B * tau) * np.cos(phi + (-4 * np.pi * r_m / lam))) / (r_m**2)

    return cfg.k * cfg.B * (first + second + third)


def numerical_mismatched_dbm(
    rng: np.random.Generator,
    cfg: Config,
    ant_gain_lin: float,
    lam: float,
    r_m: float,
    num_samples: int,
    phi: float,
) -> float:
    magnitude = (np.sqrt(ant_gain_lin) * lam) / (4 * np.pi * r_m)
    phase = -2 * np.pi * r_m / lam
    h = magnitude * np.exp(1j * (phase + phi))

    delay_samples = 2 * (r_m + cfg.electrical_length_m) * cfg.B / cfg.c  # fractional
    alpha = np.sqrt(cfg.T_bw / cfg.T_fw)

    v_fw = generate_thermal_noise_voltage(rng, cfg.k, num_samples, cfg.B, cfg.T_fw)
    v_fw_delayed = fractional_delay(v_fw, delay_samples)

    v_src_back = alpha * (h**2) * v_fw  # matches your original structure

    v_amb = generate_thermal_noise_voltage(rng, cfg.k, num_samples, cfg.B, cfg.T_amb)
    v_amb_rx = v_amb * np.sqrt(1 - np.abs(h) ** 4)

    return power_dbm_from_v(v_fw_delayed + v_src_back + v_amb_rx)


# ------------------------- IO helpers -------------------------

def load_and_calibrate_measured(cfg: Config, path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.load(path, allow_pickle=True).item()
    ranges_mm = np.array(d["Position [mm]"], dtype=float) + cfg.phase_center_mm
    p_on = np.array(d["50 Ohm"], dtype=float)
    p_off = np.array(d["Short"], dtype=float)

    slope = (cfg.T_amb + cfg.T_fw) / np.mean(p_on)
    meas_50ohm_W = cfg.k * cfg.B * slope * p_on
    meas_short_W = cfg.k * cfg.B * slope * p_off
    return ranges_mm, meas_50ohm_W, meas_short_W


# ------------------------- main -------------------------

def main():
    cfg = Config()
    set_plot_style()

    here = Path(__file__).resolve().parent
    raw_dir = here / "rawmeasurements"
    results_dir = here / "data/Figure3"
    out_dir = here / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    meas_path = raw_dir / "rawFig3fg.npy"

    rng = np.random.default_rng(cfg.seed)

    lam = cfg.c / cfg.center_f
    ant_gain_lin = 10 ** (cfg.ant_gain_dB / 10)
    num_samples = int(cfg.B * cfg.int_time)
    offset_db = 10 * np.log10(cfg.B)

    # ---- measured ----
    ranges_mm, meas_50ohm_W, meas_short_W = load_and_calibrate_measured(cfg, meas_path)
    theo_50ohm_W = np.mean(meas_50ohm_W) * np.ones_like(meas_50ohm_W)  # same as your original

    # ---- numerical/theory ----
    num_50ohm_dbm = np.zeros(len(ranges_mm))
    num_short_dbm = np.zeros(len(ranges_mm))
    theo_short_W = np.zeros(len(ranges_mm))

    for i, r_m in enumerate(ranges_mm / 1000.0):
        print(str(i + 1) + "/" + str(len(ranges_mm)))
        num_50ohm_dbm[i] = numerical_matched_dbm(rng, cfg, ant_gain_lin, lam, r_m, num_samples)
        num_short_dbm[i] = numerical_mismatched_dbm(rng, cfg, ant_gain_lin, lam, r_m, num_samples, phi=cfg.phi)
        theo_short_W[i] = theory_mismatched_power_W(cfg, ant_gain_lin, lam, r_m, phi=cfg.phi)

    # ---- plot: absolute (dBm/Hz) ----
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(ranges_mm, num_50ohm_dbm - offset_db, color="#ff4d4d", linewidth=5, label=r"Numerical $50\,\Omega$")
    ax.plot(ranges_mm, 10 * np.log10(meas_50ohm_W * 1000) - offset_db, color="#ff9999", linewidth=5, label=r"Measured $50\,\Omega$")
    ax.plot(ranges_mm, 10 * np.log10(theo_50ohm_W * 1000) - offset_db, color="#b30000", linewidth=4, linestyle="--", label=r"Theoretical $50\,\Omega$")

    ax.plot(ranges_mm, num_short_dbm - offset_db, color="#4d79ff", linewidth=5, label="Numerical Short Circuit")
    ax.plot(ranges_mm, 10 * np.log10(meas_short_W * 1000) - offset_db, color="#99ccff", linewidth=5, label="Measured Short Circuit")
    ax.plot(ranges_mm, 10 * np.log10(theo_short_W * 1000) - offset_db, color="#003399", linewidth=4, linestyle="--", label="Theoretical Short Circuit")

    ax.text(0.45, 0.92, "982 MHz Noise Bandwidth", color="black", fontsize=25, fontweight="bold", transform=ax.transAxes)
    ax.set_xlim(cfg.xlim_mm)
    ax.set_ylim(cfg.ylim_dbmhz)
    ax.set_xlabel("Range (mm)")
    ax.set_ylabel("Received Power (dBm/Hz)")
    ax.grid(True)
    ax.legend()
    ax.set_xticks(list(cfg.xticks_mm))
    ax.set_yticks(list(cfg.yticks_dbmhz))
    ax.set_xticklabels([f"{t:g}" for t in cfg.xticks_mm])
    ax.set_yticklabels([f"{t:g}" for t in cfg.yticks_dbmhz])
    fig.tight_layout(pad=2.0)
    fig.savefig(out_dir / "Fig3f.png", dpi=900)
    plt.show()

    # ---- plot: difference (W) ----
    fig, ax = plt.subplots(figsize=(10, 6))
    c_num, c_theo, c_meas = "#D95319", "#0072BD", "#EDB120"

    diff_meas_W = meas_50ohm_W - meas_short_W
    diff_num_W = np.mean(meas_50ohm_W) - (10 ** (num_short_dbm / 10)) / 1000.0
    diff_theo_W = np.mean(meas_50ohm_W) - theo_short_W

    ax.plot(ranges_mm, diff_meas_W, color=c_meas, linewidth=5, label="Measured")
    ax.plot(ranges_mm, diff_num_W, color=c_num, linewidth=5, label="Numerical")
    ax.plot(ranges_mm, diff_theo_W, color=c_theo, linestyle="--", linewidth=4, label="Theoretical")

    ax.text(0.45, 0.92, "982 MHz Noise Bandwidth", color="black", fontsize=25, fontweight="bold", transform=ax.transAxes)
    ax.set_xlim(cfg.xlim_mm)
    ax.set_ylim(cfg.ylim_diff_W)
    ax.set_xlabel("Range (mm)")
    ax.set_ylabel("Received Power Difference (W)")
    ax.grid(True, which="both", linestyle="-", alpha=0.4)
    ax.legend(loc="center right")
    ax.set_xticks(list(cfg.xticks_mm))
    ax.set_yticks(list(cfg.yticks_diff_W))
    ax.set_xticklabels([f"{t:g}" for t in cfg.xticks_mm])
    ax.set_yticklabels([f"{t:g}" for t in cfg.yticks_diff_W])
    fig.tight_layout(pad=2.0)
    fig.savefig(out_dir / "Fig3g.png", dpi=900)
    plt.show()

    # ---- save processed vectors (single file) ----
    results_vector = {
        "seed": cfg.seed,
        "B_Hz": cfg.B,
        "offset_dB": offset_db,
        "range_mm": ranges_mm,
        "measured_50ohm_dBm_per_Hz": 10 * np.log10(meas_50ohm_W * 1000) - offset_db,
        "theoretical_50ohm_dBm_per_Hz": 10 * np.log10(theo_50ohm_W * 1000) - offset_db,
        "numerical_50ohm_dBm_per_Hz": num_50ohm_dbm - offset_db,
        "measured_short_dBm_per_Hz": 10 * np.log10(meas_short_W * 1000) - offset_db,
        "theoretical_short_dBm_per_Hz": 10 * np.log10(theo_short_W * 1000) - offset_db,
        "numerical_short_dBm_per_Hz": num_short_dbm - offset_db,
        "diff_measured_W": diff_meas_W,
        "diff_numerical_W": diff_num_W,
        "diff_theoretical_W": diff_theo_W,
    }
    np.save(results_dir / "Fig3fg.npy", results_vector)


if __name__ == "__main__":
    main()