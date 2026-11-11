"""Reproduce the completed finite Fourier checks, not an advantage benchmark.

The asymptotic argument is in ../science/gap-and-signal.md. These checks
independently form the primitive Fourier multiplier, restrict its real-space
kernel to the electronic cell grid, and compare the aliased spectrum.
"""

import json
from pathlib import Path

import numpy as np


def check_size(length):
    k = 2 * np.pi * np.fft.fftfreq(length)
    q = 2 * np.pi * np.arange(length // 2) / (length // 2)
    shape = (length, length, length // 2)
    a = np.sin(k[:, None, None] / 2) ** 2 + np.sin(k[None, :, None] / 2) ** 2
    u = np.sin(q[None, None, :] / 4) ** 2
    first = np.divide(u, a + u, out=np.zeros(shape), where=(a + u) > 0)
    second = np.divide(1 - u, a + 1 - u, out=np.zeros(shape), where=(a + 1 - u) > 0)
    symbol = (first + second) / 2
    low = (
        (np.abs(k[:, None, None]) <= np.pi / 4 + 1e-14)
        & (np.abs(k[None, :, None]) <= np.pi / 4 + 1e-14)
        & (q[None, None, :] >= np.pi / 2 - 1e-14)
        & (q[None, None, :] <= 3 * np.pi / 2 + 1e-14)
    )
    high = np.broadcast_to(
        (np.abs(k[:, None, None]) >= np.pi / 2 - 1e-14)
        & (np.abs(k[None, :, None]) >= np.pi / 2 - 1e-14),
        shape,
    )
    lam = 4 * (
        np.sin(k[:, None, None] / 2) ** 2
        + np.sin(k[None, :, None] / 2) ** 2
        + np.sin(k[None, None, :] / 2) ** 2
    )
    primitive_symbol = np.divide(
        4 * np.sin(k[None, None, :] / 2) ** 2,
        lam,
        out=np.zeros((length, length, length)),
        where=lam > 0,
    )
    cell_kernel = np.fft.ifftn(primitive_symbol).real[:, :, ::2]
    aliased_symbol = np.fft.fftn(cell_kernel)
    result = {
        "L": length,
        "cell_frequency_count": int(symbol.size),
        "min_on_low_set": float(symbol[low].min()),
        "max_on_high_set": float(symbol[high].max()),
        "low_fraction": float(low.mean()),
        "high_fraction": float(high.mean()),
        "variance": float(symbol.var()),
        "claimed_variance_lower": 1 / 900,
        "independent_fft_alias_error": float(np.max(np.abs(aliased_symbol - symbol))),
        "mean_vs_diagonal_identity_error": float(abs(symbol.mean() - (1 - length**-3) / 3)),
    }
    assert result["min_on_low_set"] >= 8 / 15 - 1e-12
    assert result["max_on_high_set"] <= 1 / 3 + 1e-12
    assert result["low_fraction"] >= 1 / 32
    assert result["high_fraction"] >= 1 / 4
    assert result["variance"] >= 1 / 900
    assert result["independent_fft_alias_error"] < 1e-12
    assert result["mean_vs_diagonal_identity_error"] < 1e-12
    return result


def main():
    result = {
        "scope": "finite checks of analytic Fourier alias and uniform-signal lemma, not an asymptotic proof or advantage test",
        "results": [check_size(length) for length in (8, 16, 32)],
    }
    destination = (Path(__file__).parent / "results" / "signal_certificate.json")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
