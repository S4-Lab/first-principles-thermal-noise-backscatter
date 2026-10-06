"""
Figure 2(c,d): matched case — measured vs numerical vs theoretical.
"""

from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


# -------------------------- configuration --------------------------

@dataclass(frozen=True)
class Config:
    # physical/system params
    c: float = 3e8
    k: float = 1.38e-23
    B: float = 982_000_000
    center_f: float = 11.2e9
    T_AMB: float = 290
    T_REC: float = 99
    T_S_values: tuple[float, float] = (55, 290)
    duration: float = 0.01
    phase_center_mm: float = 20.6
    ant_gain_dB: float = 28.37

    # behavior
    seed: int = 42

# -------------------------- helpers (math/sim) --------------------------

def generate_thermal_noise_voltage(
    rng: np.random.Generator,
    num_samples: int,
    bandwidth: float,
    temperature: float = 290,
    resistance: float = 50,
) -> np.ndarray:
    k = 1.38064852e-23  # Boltzmann constant [J/K]
    voltage_sigma = np.sqrt(k * temperature * bandwidth * resistance)

    amplitude = rng.rayleigh(scale=voltage_sigma / np.sqrt(2), size=num_samples)
    phase = rng.uniform(0, 2 * np.pi, size=num_samples)
    return amplitude * np.exp(1j * phase)


def vsamples_power_dbm(signal: np.ndarray, resistance: float = 50) -> float:
    p_w = np.mean((np.abs(signal) ** 2) / resistance)
    return 10 * np.log10(p_w * 1000)


def theory_matched_temperature_K(
    r_m: float, Ts: float, G_TX_G_RX: float, lam: float, Tamb: float, Trec: float
) -> float:
    L = np.sqrt(G_TX_G_RX) * lam / (4 * np.pi)
    return Ts * (L**2) / (r_m**2) + Tamb * (1 - (L**2) / (r_m**2)) + Trec


def numerical_matched_power_dbm(
    rng: np.random.Generator,
    B: float,
    ant_gain_lin: float,
    lam: float,
    r_m: float,
    num_samples: int,
    T_S: float,
    T_AMB: float,
    T_REC: float,
) -> float:
    magnitude = (np.sqrt(ant_gain_lin) * lam) / (4 * np.pi * r_m)
    phase = -2 * np.pi * r_m / lam
    h = magnitude * np.exp(1j * phase)

    # Source noise
    V_S = generate_thermal_noise_voltage(rng, num_samples, B, T_S)
    V_S_REC = h * V_S

    # Ambient noise
    V_AMB = generate_thermal_noise_voltage(rng, num_samples, B, T_AMB)
    V_AMB_REC = V_AMB * np.sqrt(1 - np.abs(h) ** 2)

    # Receiver noise
    V_REC = generate_thermal_noise_voltage(rng, num_samples, B, T_REC)

    return vsamples_power_dbm(V_S_REC + V_AMB_REC + V_REC)


def theo_num_matched(
    rng: np.random.Generator,
    cfg: Config,
    distances_mm: np.ndarray,
    results_dir: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lam = cfg.c / cfg.center_f
    ant_gain_lin = 10 ** (cfg.ant_gain_dB / 10)
    num_samples = int(cfg.B * cfg.duration)

    distances_m = distances_mm / 1000.0

    P_num = np.zeros((len(distances_m), len(cfg.T_S_values)))

    P_theory = np.zeros_like(P_num)

    for j, T_S in enumerate(cfg.T_S_values):
        for i, r_m in enumerate(distances_m):
            print(str((i+1)+(j*len(distances_m)))+"/"+str(len(distances_m)*len(cfg.T_S_values)))
            P_num[i, j] = numerical_matched_power_dbm(
                rng, cfg.B, ant_gain_lin, lam, r_m, num_samples, T_S, cfg.T_AMB, cfg.T_REC)

            T_eq = theory_matched_temperature_K(r_m, T_S, ant_gain_lin, lam, cfg.T_AMB, cfg.T_REC)
            P_theory[i, j] = 10 * np.log10(T_eq * (cfg.k * cfg.B) * 1000)

    return P_num, P_theory, distances_mm


# -------------------------- plotting --------------------------

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


def plot_raw(
    out_path: Path,
    distances_mm: np.ndarray,
    y_amb_w: np.ndarray,
    y_lna_w: np.ndarray,
    P_num_dbm: np.ndarray,
    P_theory_dbm: np.ndarray,
    B: float,
):
    set_plot_style()
    fig, ax = plt.subplots(figsize=(10, 6))
    offset = 10 * np.log10(B)

    ax.plot(distances_mm, 10 * np.log10(y_amb_w * 1000) - offset, color="#ff9999", linewidth=5, label=r"Measured 50 $\Omega$")
    ax.plot(distances_mm, P_num_dbm[:, 1] - offset, color="#ff4d4d", linewidth=5, label=r"Numerical 50 $\Omega$")
    ax.plot(distances_mm, P_theory_dbm[:, 1] - offset, color="#b30000", linewidth=4, linestyle="--", label=r"Theoretical 50 $\Omega$")

    ax.plot(distances_mm, 10 * np.log10(y_lna_w * 1000) - offset, color="#99ccff", linewidth=5, label="Measured LNA Input")
    ax.plot(distances_mm, P_num_dbm[:, 0] - offset, color="#4d79ff", linewidth=5, label="Numerical LNA Input")
    ax.plot(distances_mm, P_theory_dbm[:, 0] - offset, color="#003399", linewidth=4, linestyle="--", label="Theoretical LNA Input")

    ax.set_ylim([-173.35, -172.55])
    ax.set_xlabel("Range (mm)")
    ax.set_ylabel("Received Power (dBm/Hz)")
    ax.grid(True)
    ax.legend()
    ax.text(0.48, 0.92, "982 MHz Noise Bandwidth", color="black", fontsize=25, fontweight="bold", transform=ax.transAxes)

    x_ticks = [1e2, 2e2, 3e2, 4e2, 5e2, 6e2]
    y_ticks = [-172.6, -172.7, -172.8, -172.9, -173.0, -173.1, -173.2, -173.3]
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticklabels([f"{t:g}" for t in x_ticks])
    ax.set_yticklabels([f"{t:g}" for t in y_ticks])

    fig.tight_layout(pad=2.0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=900, pad_inches=0.25)
    plt.show()


def plot_difference(
    out_path: Path,
    distances_mm: np.ndarray,
    y_amb_w: np.ndarray,
    y_lna_w: np.ndarray,
    P_num_dbm: np.ndarray,
    P_theory_dbm: np.ndarray,
):
    set_plot_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    c1, c2, c3 = "#D95319", "#0072BD", "#EDB120"
    ax.loglog(distances_mm, y_amb_w - y_lna_w, color=c1, linewidth=5, label="Measured")
    ax.loglog(
        distances_mm,
        ((10 ** (P_num_dbm[:, 1] / 10)) - (10 ** (P_num_dbm[:, 0] / 10))) / 1000,
        color=c2,
        linewidth=5,
        label="Numerical",
    )
    ax.loglog(
        distances_mm,
        ((10 ** (P_theory_dbm[:, 1] / 10)) - (10 ** (P_theory_dbm[:, 0] / 10))) / 1000,
        color=c3,
        linestyle="--",
        linewidth=4,
        label="Theoretical",
    )

    # manual ticks (match original)
    x_ticks = [1e2, 2e2, 3e2, 4e2, 5e2, 6e2]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f"{t:g}" for t in x_ticks])

    ax.set_yticks([1e-13])
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.set_yticklabels([r"10$^{-13}$"], fontname="Helvetica")

    ax.grid(True, which="both", linestyle="-", alpha=0.4)
    ax.set_xlabel("Range (mm)")
    ax.set_ylabel("Received Power Difference (W)")
    ax.legend()

    # slope annotation
    x = np.log10(distances_mm / 1000)
    y = np.log10((y_amb_w - y_lna_w) * 1000)
    ax.set_ylim([7e-15, 1.5e-12])
    slope, _logA = np.polyfit(x, y, 1)
    ax.text(0.22, 0.5, f"slope = {slope:.2f}", color="black", fontsize=25, fontweight="bold", transform=ax.transAxes)
    ax.text(0.48, 0.92, "982 MHz Noise Bandwidth", color="black", fontsize=25, fontweight="bold", transform=ax.transAxes)

    fig.tight_layout(pad=2.0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=900, pad_inches=0.25)
    plt.show()


# -------------------------- main --------------------------

def main():
    cfg = Config()
    rng = np.random.default_rng(cfg.seed)

    # paths relative to *this* file
    here = Path(__file__).resolve().parent
    data_dir = here / "rawmeasurements"
    out_dir = here / "out"
    results_dir = here / "data/Figure2"

    measured_path = data_dir / "rawFig2cd.npy"
    raw_fig_path = out_dir / "Fig2c.png"
    dif_fig_path = out_dir / "Fig2d.png"

    # --- load measured data ---
    data = np.load(measured_path, allow_pickle=True).item()
    distances_mm = np.array(data["Position [mm]"], dtype=float) + cfg.phase_center_mm

    y_AMB = np.array(data["50 Ohm"], dtype=float)   # calibration point
    sdr_slope = (cfg.T_AMB + cfg.T_REC) / np.mean(y_AMB)
    y_LNA = cfg.k * cfg.B * sdr_slope * np.array(data["LNA Input"], dtype=float)
    y_AMB = cfg.k * cfg.B * sdr_slope * y_AMB

    # --- cleaning ---
    min_diff = 1e-14
    mask_switch = (y_AMB - y_LNA) > min_diff
    distances_mm = distances_mm[mask_switch]
    y_LNA_clean = y_LNA[mask_switch]
    y_AMB_clean = y_AMB[mask_switch]

    # --- compute numerical + theory ---
    P_num_dbm, P_theory_dbm, distances_mm = theo_num_matched(rng, cfg, distances_mm, results_dir)

    # --- plots ---
    plot_raw(raw_fig_path, distances_mm, y_AMB_clean, y_LNA_clean, P_num_dbm, P_theory_dbm, cfg.B)
    plot_difference(dif_fig_path, distances_mm, y_AMB_clean, y_LNA_clean, P_num_dbm, P_theory_dbm)

    # --- save processed vectors for sharing/replotting ---
    offset = 10 * np.log10(cfg.B)
    results_vector = {
        "distances_mm": distances_mm,
        "measured_50_dBm_per_Hz": 10 * np.log10(y_AMB_clean * 1000) - offset,
        "numerical_50_dBm_per_Hz": P_num_dbm[:, 1] - offset,
        "theoretical_50_dBm_per_Hz": P_theory_dbm[:, 1] - offset,
        "measured_LNAInput_dBm_per_Hz": 10 * np.log10(y_LNA_clean * 1000) - offset,
        "numerical_LNAInput_dBm_per_Hz": P_num_dbm[:, 0] - offset,
        "theoretical_LNAInput_dBm_per_Hz": P_theory_dbm[:, 0] - offset,
    }
    np.save(results_dir / "Fig2cd.npy", results_vector)


if __name__ == "__main__":
    main()