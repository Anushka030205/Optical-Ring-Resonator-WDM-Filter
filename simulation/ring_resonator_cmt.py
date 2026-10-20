"""Coupled-mode-theory (CMT) simulation of dual- and triple-ring add-drop resonators.

Python port of simulation/matlab/dual_ring_add_drop.m and
simulation/matlab/triple_ring_add_drop.m. Solves M * A = B at every wavelength
point, then extracts resonance, FWHM, Q factor, FSR and five-channel DWDM
crosstalk at the drop port.

Usage (from the repository root):
    python simulation/ring_resonator_cmt.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["axes.formatter.useoffset"] = False

RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS.mkdir(exist_ok=True)

# ---------------- Shared parameters (Table 1 of the report) ----------------
C = 3e8
LAMBDA = np.linspace(1511.2e-9, 1511.9e-9, 200_000)
LAMBDA0 = 1511.54e-9
NG = 4.2          # group index
L = 100e-6        # circumference of every ring (m)
K_BUS = K_DROP = 0.01
K12 = K23 = 0.02
CHANNEL_SPACING = 0.4e-9   # ~50 GHz
NUM_CHANNELS = 5
SIGMA = 0.08e-9            # Gaussian channel width

COLORS = {"dual": "#2a6fbb", "triple": "#d4572a"}


def simulate(n_rings, loss_dB_cm, s_add, add_in_wdm):
    """Run the CMT model for a chain of `n_rings` (2 or 3) rings."""
    omega = 2 * np.pi * C / LAMBDA
    omega0 = 2 * np.pi * C / LAMBDA0
    delta = omega - omega0

    alpha = loss_dB_cm / 4.343 * 100          # Np/m
    a_loss = np.exp(-alpha * L)
    tau = C / (NG * L)
    g_i = (1 - a_loss) * tau
    g_e1, g_e2 = K_BUS * tau, K_DROP * tau
    k12, k23 = K12 * tau, K23 * tau

    n = len(LAMBDA)
    M = np.zeros((n, n_rings, n_rings), dtype=complex)
    M[:, 0, 0] = 1j * delta + g_i + g_e1 + g_e2
    for r in range(1, n_rings):
        M[:, r, r] = 1j * delta + g_i
    M[:, 0, 1] = M[:, 1, 0] = -1j * k12
    if n_rings == 3:
        M[:, 1, 2] = M[:, 2, 1] = -1j * k23

    B = np.zeros((n, n_rings, 1), dtype=complex)
    B[:, 0, 0] = np.sqrt(2 * g_e1) * 1.0 + np.sqrt(2 * g_e2) * s_add
    a1 = np.linalg.solve(M, B)[:, 0, 0]

    t_through = np.abs(1.0 - np.sqrt(2 * g_e1) * a1) ** 2
    t_drop = np.abs(np.sqrt(2 * g_e2) * a1 - s_add) ** 2

    # --- Q factor from drop-port FWHM ---
    idx = int(np.argmax(t_drop))
    lam_res = LAMBDA[idx]
    half = t_drop[idx] / 2
    left = np.where(t_drop[: idx + 1] <= half)[0][-1]
    right = np.where(t_drop[idx:] <= half)[0][0] + idx
    fwhm = LAMBDA[right] - LAMBDA[left]
    q = lam_res / fwhm
    fsr = lam_res**2 / (NG * L)

    # --- WDM filtering and crosstalk ---
    channels = lam_res + np.arange(-2, 3) * CHANNEL_SPACING
    spectrum = sum(np.exp(-((LAMBDA - c) ** 2) / (2 * SIGMA**2)) for c in channels)
    spectrum = spectrum / spectrum.max()
    drive = spectrum
    if add_in_wdm:
        drive = spectrum + np.exp(-((LAMBDA - lam_res) ** 2) / (2 * SIGMA**2))
    out = drive * t_drop

    power = np.array([out[np.argmin(np.abs(LAMBDA - c))] for c in channels])
    power = np.maximum(power, 1e-12)
    crosstalk = 10 * np.log10(power / power[NUM_CHANNELS // 2])

    return dict(t_through=t_through, t_drop=t_drop, lam_res=lam_res, fwhm=fwhm, q=q,
                fsr=fsr, channels=channels, spectrum=spectrum, out=out, crosstalk=crosstalk)


def plot_spectra(res, name, title, fname):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(LAMBDA * 1e9, res["t_through"], color="#555555", lw=1.4, label="Through port")
    ax.plot(LAMBDA * 1e9, res["t_drop"], color=COLORS[name], lw=1.6, label="Drop port")
    ax.set(xlabel="Wavelength (nm)", ylabel="Transmission", title=title)
    ax.grid(alpha=0.35)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / fname, dpi=150)
    plt.close(fig)


def plot_wdm(res, name, title, fname):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(LAMBDA * 1e9, res["spectrum"], "--", color="#555555", lw=1.1, label="Input WDM")
    ax.plot(LAMBDA * 1e9, res["out"], color=COLORS[name], lw=1.6, label="After filter")
    ax.set(xlabel="Wavelength (nm)", ylabel="Normalized power", title=title)
    ax.grid(alpha=0.35)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / fname, dpi=150)
    plt.close(fig)


def main():
    # Parameters exactly as used in the report: dual 2 dB/cm, add 0.05; triple 0.5 dB/cm, add 0.02
    dual = simulate(2, loss_dB_cm=2.0, s_add=0.05, add_in_wdm=False)
    triple = simulate(3, loss_dB_cm=0.5, s_add=0.02, add_in_wdm=True)

    plot_spectra(dual, "dual", "Dual Ring Add-Drop: Through and Drop Port", "dual_ring_spectrum.png")
    plot_wdm(dual, "dual", "WDM Filtering with Dual Ring", "dual_ring_wdm.png")
    plot_spectra(triple, "triple", "Triple Ring Add-Drop: Through and Drop Port", "triple_ring_spectrum.png")
    plot_wdm(triple, "triple", "WDM Filtering with Triple Ring", "triple_ring_wdm.png")

    # Drop-port overlay (normalised) and crosstalk comparison
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for name, res, lab in (("dual", dual, "Dual ring"), ("triple", triple, "Triple ring")):
        ax.plot(LAMBDA * 1e9, res["t_drop"] / res["t_drop"].max(), color=COLORS[name], lw=1.6, label=lab)
    ax.set(xlabel="Wavelength (nm)", ylabel="Normalized drop-port transmission",
           title="Drop-Port Response: Dual vs Triple Ring")
    ax.grid(alpha=0.35)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "drop_port_comparison.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = np.arange(1, 6)
    w = 0.38
    ax.bar(x - w / 2, dual["crosstalk"], w, color=COLORS["dual"], label="Dual ring")
    ax.bar(x + w / 2, triple["crosstalk"], w, color=COLORS["triple"], label="Triple ring")
    ax.set(xlabel="Channel", ylabel="Crosstalk (dB)", title="WDM Crosstalk by Channel", xticks=x)
    ax.grid(axis="y", alpha=0.35)
    ax.legend(loc="lower center")
    fig.tight_layout()
    fig.savefig(RESULTS / "crosstalk_comparison.png", dpi=150)
    plt.close(fig)

    # Tables
    lines = ["config,resonance_nm,FWHM_pm,Q,FSR_nm"]
    for name, r in (("dual", dual), ("triple", triple)):
        lines.append(f"{name},{r['lam_res']*1e9:.4f},{r['fwhm']*1e12:.4f},{r['q']:.4e},{r['fsr']*1e9:.4f}")
    (RESULTS / "metrics.csv").write_text("\n".join(lines) + "\n")

    lines = ["channel,dual_wavelength_nm,dual_crosstalk_dB,triple_wavelength_nm,triple_crosstalk_dB"]
    for i in range(NUM_CHANNELS):
        lines.append(f"{i+1},{dual['channels'][i]*1e9:.2f},{dual['crosstalk'][i]:.2f},{triple['channels'][i]*1e9:.2f},{triple['crosstalk'][i]:.2f}")
    (RESULTS / "crosstalk.csv").write_text("\n".join(lines) + "\n")

    for name, r in (("DUAL RING", dual), ("TRIPLE RING", triple)):
        print(f"\n===== {name} =====")
        print(f"Resonance = {r['lam_res']*1e9:.4f} nm | FWHM = {r['fwhm']*1e12:.4f} pm | "
              f"Q = {r['q']:.3e} | FSR = {r['fsr']*1e9:.4f} nm")
        for i, ct in enumerate(r["crosstalk"], 1):
            print(f"  Channel {i} crosstalk = {ct:.2f} dB")


if __name__ == "__main__":
    main()
