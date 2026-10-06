"""
Figure 3(d,e): mismatched low bandwidth case — measured vs numerical vs theoretical.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# ------------------------- configuration -------------------------

@dataclass(frozen=True)
class Config:
    # RNG
    seed: int = 42

    # system parameters
    c: float = 3e8
    k: float = 1.38e-23
    B: float = 50_000_000
    center_f: float = 11e9
    T_amb: float = 290
    T_fw: float = 99
    T_bw: float = 55
    int_time: float = 0.1
    phase_center_mm: float = 20.6
    electrical_length_m: float = 0.42
    ant_gain_dB: float = 28.37

    # plotting
    ylim: tuple[float, float] = (-173.1, -172.5)
    x_ticks: tuple[float, ...] = (1e2, 2e2, 3e2, 4e2, 5e2)
    y_ticks: tuple[float, ...] = (-172.5, -172.6, -172.7, -172.8, -172.9, -173.0, -173.1)


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

def calibrate_measured(cfg: Config, meas_on: np.ndarray, meas_off: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    slope = (cfg.T_amb + cfg.T_fw) / np.mean(meas_on)
    cal_50ohm_W = cfg.k * cfg.B * slope * meas_on
    meas_mismatch_W = cfg.k * cfg.B * slope * meas_off
    return cal_50ohm_W, meas_mismatch_W


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

# ------------------------- models -------------------------

def theory_mismatched_power_W(cfg: Config, r_m: float, phi: float = 0.0) -> float:
    lam = cfg.c / cfg.center_f
    ant_gain = 10 ** (cfg.ant_gain_dB / 10)

    tau = 2 * (r_m + cfg.electrical_length_m) / cfg.c
    L = np.sqrt(ant_gain) * lam / (4 * np.pi)

    return (cfg.k * cfg.B) * (cfg.T_amb + cfg.T_fw + ((cfg.T_bw - cfg.T_amb) * L**4 / r_m**4)
        + (2 * np.sqrt(cfg.T_bw * cfg.T_fw) * L**2 * np.sinc(cfg.B * tau) / r_m**2
            * np.cos(-4 * np.pi / lam * r_m + phi)))


def numerical_mismatched_power_dbm(rng: np.random.Generator, cfg: Config, r_m: float, phi: float = 0.0) -> float:
    lam = cfg.c / cfg.center_f
    ant_gain = 10 ** (cfg.ant_gain_dB / 10)
    num_samples = int(cfg.B * cfg.int_time)

    magnitude = (np.sqrt(ant_gain) * lam) / (4 * np.pi * r_m)
    phase = -2 * np.pi * r_m / lam
    h = magnitude * np.exp(1j * (phase + phi))

    sample_delay = int(np.round(2 * (r_m + cfg.electrical_length_m) * cfg.B / cfg.c))
    alpha = np.sqrt(cfg.T_bw / cfg.T_fw)

    # forwared propagating noise
    v_rec = generate_thermal_noise_voltage(rng, cfg.k, num_samples + sample_delay, cfg.B, cfg.T_fw)

    # ambient noise
    v_amb = generate_thermal_noise_voltage(rng, cfg.k, num_samples, cfg.B, cfg.T_amb)
    v_amb_rx = v_amb * np.sqrt(1 - np.abs(h) ** 4)

    # backscattered noise
    v_src_rx = alpha * (h**2) * v_rec[sample_delay : num_samples + sample_delay]

    v_rec_curr = v_rec[:num_samples]
    p_W = np.mean((np.abs(v_src_rx + v_amb_rx + v_rec_curr) ** 2) / 50.0)
    return 10 * np.log10(p_W * 1000)  # dBm


def simulate_curves(cfg: Config, ranges_mm: np.ndarray, phi: float, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    num_dbm = np.zeros(len(ranges_mm))
    theo_W = np.zeros(len(ranges_mm))
    for i, r_m in enumerate(ranges_mm / 1000.0):
        print(str(i+1)+"/"+str(len(ranges_mm)))
        num_dbm[i] = numerical_mismatched_power_dbm(rng, cfg, r_m, phi=phi)
        theo_W[i] = theory_mismatched_power_W(cfg, r_m, phi=phi)
    return num_dbm, theo_W

# ------------------------- plot helpers -------------------------

def plot_case(
    cfg: Config,
    ranges_mm: np.ndarray,
    meas_mismatch_W: np.ndarray,
    num_mismatch_dbm: np.ndarray,
    theo_mismatch_W: np.ndarray,
    label_case: str,
    out_path: Path,
):
    set_plot_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    offset_db = 10 * np.log10(cfg.B)

    ax.plot(ranges_mm, num_mismatch_dbm - offset_db, color="#D95319", linewidth=5, label=f"Numerical {label_case}")
    ax.plot(ranges_mm, 10 * np.log10(meas_mismatch_W * 1000) - offset_db, color="#0072BD", linewidth=5, label=f"Measured {label_case}")
    ax.plot(ranges_mm, 10 * np.log10(theo_mismatch_W * 1000) - offset_db, color="#EDB120", linewidth=4, linestyle="--", label=f"Theoretical {label_case}")

    ax.text(0.48, 0.92, "50 MHz Noise Bandwidth", color="black", fontsize=25, fontweight="bold", transform=ax.transAxes)

    ax.set_xlabel("Range (mm)")
    ax.set_ylabel("Received Power (dBm/Hz)")
    ax.grid(True)
    ax.legend()
    ax.set_ylim(cfg.ylim)

    ax.set_xticks(list(cfg.x_ticks))
    ax.set_yticks(list(cfg.y_ticks))
    ax.set_xticklabels([f"{t:g}" for t in cfg.x_ticks])
    ax.set_yticklabels([f"{t:g}" for t in cfg.y_ticks])

    fig.tight_layout(pad=2.0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=900)
    plt.show()


# ------------------------- main -------------------------

def main():
    cfg = Config()

    here = Path(__file__).resolve().parent
    raw_dir = here / "rawmeasurements"
    results_dir = here / "data/Figure3"
    out_dir = here / "out"

    short_file = raw_dir / "rawFig3e.npy"
    open_file = raw_dir / "rawFig3d.npy"

    rng = np.random.default_rng(cfg.seed)

    # ---------------- SHORT ----------------
    d = np.load(short_file, allow_pickle=True).item()
    r_mm = np.array(d["Position [mm]"], dtype=float)
    meas_on = np.array(d["50 Ohm"], dtype=float)
    meas_off = np.array(d["Short"], dtype=float)
    range_mm_short = r_mm + cfg.phase_center_mm
    cal_50_W_short, meas_short_W = calibrate_measured(cfg, meas_on, meas_off)
    num_short_dbm, theo_short_W = simulate_curves(cfg, range_mm_short, phi=np.pi / 8, rng=rng)

    plot_case(cfg, range_mm_short, meas_short_W, num_short_dbm, theo_short_W, "Short Circuit", out_dir / "Fig3e.png")

    # ---------------- OPEN ----------------
    d = np.load(open_file, allow_pickle=True).item()
    r_mm = np.array(d["Position [mm]"], dtype=float)
    meas_50 = np.array(d["50 Ohm"], dtype=float)
    meas_Open = np.array(d["Open"], dtype=float)
    range_mm_open = r_mm + cfg.phase_center_mm
    cal_50_W_open, meas_open_W = calibrate_measured(cfg, meas_50, meas_Open)
    num_open_dbm, theo_open_W = simulate_curves(cfg, range_mm_open, phi=0.0, rng=rng)

    plot_case(cfg, range_mm_open, meas_open_W, num_open_dbm, theo_open_W, "Open Circuit", out_dir / "Fig3d.png")

    # ---------------- SAVE BOTH TOGETHER ----------------
    offset = 10 * np.log10(cfg.B)
    results_vector = {
        "B_Hz": cfg.B,
        "offset_dB": offset,
        "seed": cfg.seed,

        "short": {
            "range_mm": range_mm_short,
            "measured_50ohm_cal_dBm_per_Hz": 10 * np.log10(cal_50_W_short * 1000) - offset,
            "measured_short_dBm_per_Hz": 10 * np.log10(meas_short_W * 1000) - offset,
            "numerical_short_dBm_per_Hz": num_short_dbm - offset,
            "theoretical_short_dBm_per_Hz": 10 * np.log10(theo_short_W * 1000) - offset,
        },
        "open": {
            "range_mm": range_mm_open,
            "measured_50ohm_cal_dBm_per_Hz": 10 * np.log10(cal_50_W_open * 1000) - offset,
            "measured_open_dBm_per_Hz": 10 * np.log10(meas_open_W * 1000) - offset,
            "numerical_open_dBm_per_Hz": num_open_dbm - offset,
            "theoretical_open_dBm_per_Hz": 10 * np.log10(theo_open_W * 1000) - offset,
        },
    }

    np.save(results_dir / "Fig3de.npy", results_vector)


if __name__ == "__main__":
    main()