#!/usr/bin/env python3
"""Reproduce the GHZ-256 manuscript's Supplementary Figure S1.

Dependencies: Python 3, numpy, matplotlib
    python -m pip install numpy matplotlib

Usage:
    python ghz_figure.py ghz_counts_fidelity.json
    python ghz_figure.py "ghz_counts_fidelity(2).json" --output ghz_figure.png
    python ghz_figure.py ghz_counts_fidelity.json --output ghz_figure.pdf

Input schema:
    {"qubits": 256, "records": [
        {"phi": null, "shots": 15000, "counts": {"<256-bit string>": count}},
        {"phi": 0.0, "shots": 59, "counts": {"<256-bit string>": count}},
        ...
    ]}

The input must contain one Z-basis record and all 257 equatorial settings
phi_k = k*pi/257, k = 0,...,256, in acquisition order. Full-precision phases
are read from the JSON, not from the rounded acquisition log. Strict JSON
is required; use the corrected supplementary JSON without trailing commas.

This standalone extraction preserves the plotting style used in the paper.
All plotted quantities are recomputed from counts, not copied from result
fields. Panel (b) reorders all 257 settings by theta = (256*phi) mod (2*pi)
without binning or averaging. Error bars use the independent-shot variance
estimate (1-parity**2)/(shots-1). No readout correction or filtering is
applied. The printed Hoeffding interval covers counting uncertainty under
independent trials and stable measurement settings, not systematic errors.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def analyze(path):
    """Validate the supplied GHZ-256 data and recompute its estimators."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("qubits") != 256:
        raise ValueError("This script reproduces the 256-qubit GHZ experiment.")
    records = data["records"]
    if len(records) != 258 or records[0].get("phi") is not None:
        raise ValueError("Expected one Z-basis record followed by 257 phase records.")
    for index, record in enumerate(records):
        shots, counts = record["shots"], record["counts"]
        if type(shots) is not int or shots <= 1:
            raise ValueError(f"Record {index}: shots must be an integer greater than 1.")
        if not counts or any(
            len(bits) != 256 or set(bits) - {"0", "1"}
            or type(count) is not int or count < 0
            for bits, count in counts.items()
        ):
            raise ValueError(f"Record {index}: invalid bitstrings or counts.")
        if sum(counts.values()) != shots:
            raise ValueError(f"Record {index}: counts do not sum to shots.")

    z = records[0]
    n_zero = z["counts"].get("0" * 256, 0)
    n_one = z["counts"].get("1" * 256, 0)
    z_counts = np.array([n_zero, n_one, z["shots"] - n_zero - n_one])
    population = (n_zero + n_one) / z["shots"]
    phase_records = records[1:]
    phi = np.array([record["phi"] for record in phase_records], dtype=float)
    shots = np.array([record["shots"] for record in phase_records])
    if not np.allclose(phi, np.arange(257) * np.pi / 257, rtol=0, atol=1e-12):
        raise ValueError("Expected the full-precision grid phi_k = k*pi/257.")
    even = np.array([
        sum(count for bits, count in record["counts"].items()
            if bits.count("1") % 2 == 0)
        for record in phase_records
    ])
    parity = (2 * even - shots) / shots
    k = len(phi)
    cosine, sine = np.cos(256 * phi), np.sin(256 * phi)
    weights = np.array([2 * cosine / k, 2 * sine / k])
    a, b = weights @ parity
    magnitude = np.hypot(a, b)
    offset = np.arctan2(-b, a)
    variance = (1 - parity**2) / (shots - 1)
    covariance = (weights * variance) @ weights.T
    population_variance = population * (1 - population) / (z["shots"] - 1)
    fidelity = (population + a) / 2
    fidelity_se = np.sqrt(population_variance + covariance[0, 0]) / 2
    half_width = np.sqrt(np.log(40) / 2 * (
        1 / (4 * z["shots"]) + 4 / k**2 * np.sum(cosine**2 / shots)
    ))
    predicted = a * cosine + b * sine
    total_variation = np.sum((parity - parity.mean())**2)
    r_squared = (
        1 - np.sum((parity - predicted)**2) / total_variation
        if total_variation > 0 else float("nan")
    )
    summary = {
        "qubits": 256,
        "phase_settings": k,
        "z_shots": z["shots"],
        "rotated_shots": int(shots.sum()),
        "total_shots": z["shots"] + int(shots.sum()),
        "population": float(population),
        "cosine_component": float(a),
        "sine_component": float(b),
        "coherence_magnitude_raw": float(magnitude),
        "phase_offset_radians": float(offset),
        "fidelity_ghz_plus_raw": float(fidelity),
        "standard_error_estimate": float(fidelity_se),
        "confidence_95_hoeffding_raw": [
            float(fidelity - half_width), float(fidelity + half_width)
        ],
        "confidence_95_hoeffding_physical": [
            float(max(0, fidelity - half_width)),
            float(min(1, fidelity + half_width)),
        ],
        "phase_optimized_fidelity_raw": float((population + magnitude) / 2),
        "unweighted_r_squared": float(r_squared),
    }
    return summary, z_counts, phi, parity, variance


def plot(summary, z_counts, phi, parity, variance, output, dpi=300):
    """Save the same two-panel layout used in the manuscript."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.labelsize": 10,
    })
    fig, axes = plt.subplots(
        1, 2, figsize=(8.3, 3.45), gridspec_kw={"width_ratios": [1, 2.1]}
    )
    ax = axes[0]
    probabilities = z_counts / summary["z_shots"]
    ax.bar(np.arange(3), probabilities,
           color=["#244969", "#54849b", "#b9c2c9"], width=.62)
    for index, (value, count) in enumerate(zip(probabilities, z_counts)):
        ax.text(index, value + .015, f"{count:,}", ha="center", fontsize=9)
    ax.set_xticks([0, 1, 2], [r"$0^{256}$", r"$1^{256}$", "Other"])
    ax.set_ylim(0, .64)
    ax.set_ylabel("Measured probability")
    ax.set_title("a  Computational basis", loc="left", fontsize=11,
                 fontweight="bold", pad=12)
    ax.text(.5, .89,
            f'{summary["z_shots"]:,} shots; P = {summary["population"]:.4f}',
            transform=ax.transAxes, ha="center", fontsize=9)

    ax = axes[1]
    theta = np.mod(256 * phi, 2 * np.pi)
    order = np.argsort(theta)
    ax.errorbar(
        theta[order], parity[order], yerr=np.sqrt(variance)[order],
        fmt="o", ms=2.6, elinewidth=.55, alpha=.60, color="#315d82",
        label="Measured parity ± 1 SE",
    )
    t = np.linspace(0, 2 * np.pi, 1200)
    ax.plot(
        t, summary["coherence_magnitude_raw"]
        * np.cos(t + summary["phase_offset_radians"]),
        color="#b44632", lw=1.6, label="256th Fourier harmonic",
    )
    ax.set_xlim(-.05, 2 * np.pi + .05)
    ax.set_ylim(-1.23, 1.26)
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
                  ["0", "π/2", "π", "3π/2", "2π"])
    ax.set_xlabel("Effective phase θ = (256φ) mod 2π (rad)")
    ax.set_ylabel("Parity")
    ax.set_title("b  Phase-resolved coherence", loc="left", fontsize=11,
                 fontweight="bold", pad=12)
    ax.legend(fontsize=8.3, loc="upper center", frameon=False,
              bbox_to_anchor=(.5, -.23), ncol=1)
    ax.text(
        .5, .94,
        f'{summary["phase_settings"]} settings; '
        f'{summary["rotated_shots"]:,} shots; '
        f'R² = {summary["unweighted_r_squared"]:.4f}',
        ha="center", transform=ax.transAxes, fontsize=9,
    )
    fig.subplots_adjust(left=.085, right=.985, top=.86, bottom=.28, wspace=.39)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", default="ghz_counts_fidelity.json",
                        type=Path, help="Input count JSON (default: %(default)s)")
    parser.add_argument("-o", "--output", type=Path,
                        default=Path("ghz_supplementary_figure_s1.png"),
                        help="Output PNG, PDF, or SVG (default: %(default)s)")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Raster resolution (default: %(default)s)")
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    try:
        analysis = analyze(args.input)
        plot(*analysis, output=args.output, dpi=args.dpi)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(json.dumps(analysis[0], indent=2, allow_nan=False))
    print(f"Saved figure: {args.output.resolve()}")


if __name__ == "__main__":
    main()
