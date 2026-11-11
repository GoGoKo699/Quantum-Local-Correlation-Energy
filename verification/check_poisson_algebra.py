#!/usr/bin/env python3
"""Independent finite checks of the procedural Poisson LCCSD construction.

The tiny L=2 torus is an algebra check, not a material or scaling benchmark.
Uses an independent exact fermion operator builder; all Poisson, orbital,
Gram, source and rank-changing formulas below are constructed independently.
"""
from __future__ import annotations
import itertools
import json
from pathlib import Path
import numpy as np

import fermion_reference as fermion


def main():
    L, g = 2, 0.005
    xyz = np.array(list(itertools.product(range(L), repeat=3)))
    momenta = 2 * np.pi * xyz / L
    N, V = len(xyz), len(xyz) // 2
    fourier = np.exp(-1j * momenta @ xyz.T) / np.sqrt(N)
    lam = 4 * np.sum(np.sin(momenta / 2) ** 2, axis=1)
    inv_lam = np.zeros(N)
    inv_lam[lam > 1e-12] = 1 / lam[lam > 1e-12]
    coulomb = (fourier.conj().T @ np.diag(4 * np.pi * g * inv_lam) @ fourier).real
    Z = 7 * np.pi * (1 - 1 / L)
    G0 = float(np.diag(coulomb).max() / g)
    Gbar = 3 * np.pi + np.pi**3 / 24
    green_columns = np.diag(np.sqrt(4 * np.pi * inv_lam / Z)) @ fourier

    p = np.array([0.20, 0.35, 0.60, 0.80])
    gamma = np.sqrt(p * (1 - p))
    orbital = np.zeros((N, N))
    transition = np.zeros((N, V))
    for x in range(V):
        a, b = 2 * x, 2 * x + 1
        orbital[a, a], orbital[b, a] = np.sqrt(p[x]), np.sqrt(1 - p[x])
        orbital[a, b], orbital[b, b] = -np.sqrt(1 - p[x]), np.sqrt(p[x])
        transition[a, x], transition[b, x] = -gamma[x], gamma[x]
    occupied, virtual = tuple(range(0, N, 2)), tuple(range(1, N, 2))
    fock = np.zeros((N, N))
    for x in range(V):
        fock[2*x, 2*x] = -0.5 - 0.01*x
        fock[2*x+1, 2*x+1] = 0.5 + 0.007*x
    for x in range(V):
        for y in range(x+1, V):
            if np.count_nonzero(xyz[2*x] != xyz[2*y]) == 1:
                fock[2*x, 2*y] = fock[2*y, 2*x] = 0.013
                fock[2*x+1, 2*y+1] = fock[2*y+1, 2*x+1] = -0.017
    density = np.einsum("xp,xq->xpq", orbital, orbital)
    eri = np.einsum("xpq,xy,yrs->pqrs", density, coulomb, density, optimize=True)
    one_body = fock.copy()
    for i in occupied:
        one_body -= eri[:, :, i, i] - eri[:, i, i, :]
    model = {"norb": N, "nelec": V, "one_body": one_body, "eri": eri}
    states, position, H = fermion.hamiltonian(model)
    model_f = {"norb": N, "nelec": V, "one_body": fock, "eri": np.zeros_like(eri)}
    _, _, HF = fermion.hamiltonian(model_f)
    reference = sum(1 << i for i in occupied)
    ref = position[reference]
    excitations = []
    for rank in (1, 2):
        for holes in itertools.combinations(occupied, rank):
            for particles in itertools.combinations(virtual, rank):
                state, sign = fermion.act(reference, fermion.excitation_operations(holes, particles))
                excitations.append((holes, particles, state, sign))
    indices = np.array([position[e[2]] for e in excitations])
    signs = np.array([e[3] for e in excitations])
    A = signs[:, None] * H[np.ix_(indices, indices)] * signs[None, :] - H[ref, ref] * np.eye(len(indices))
    F = signs[:, None] * HF[np.ix_(indices, indices)] * signs[None, :] - HF[ref, ref] * np.eye(len(indices))
    b = signs * H[indices, ref]
    ns = V*V
    K = A[ns:, :ns]
    J = transition.T @ coulomb @ transition
    W = np.diag(np.sqrt(4 * np.pi * g * inv_lam)) @ fourier @ transition
    alpha_J = 4 * np.pi * g * float(gamma.max())**2
    B0 = alpha_J / np.sqrt(2)
    alpha_K = 4 * g * float(gamma.max()) * np.sqrt(2 * np.pi * Z)
    b_formula = np.zeros_like(b)
    K_formula = np.zeros_like(K)
    target_by_state = {e[2]: (r, e[3]) for r, e in enumerate(excitations[ns:])}
    for x, y in itertools.combinations(range(V), 2):
        for r, (holes, particles, _, _) in enumerate(excitations):
            if holes == (2*x, 2*y) and particles == (2*x+1, 2*y+1):
                b_formula[r] = J[x, y]
    for c, (holes, particles, state, sign) in enumerate(excitations[:ns]):
        charge = density[:, particles[0], particles[0]] - density[:, holes[0], holes[0]]
        column = transition.T @ coulomb @ charge
        for z in range(V):
            result = fermion.act(state, [("ann", 2*z), ("cre", 2*z+1)])
            if result is None or result[0] not in target_by_state:
                continue
            r, target_sign = target_by_state[result[0]]
            K_formula[r, c] += target_sign * result[1] * sign * column[z]
    interaction_norm = float(np.linalg.norm(A-F, ord=2))
    Gamma = g * (6*Gbar + 16*np.pi*float(gamma.max())**2 + 4*float(gamma.max())*np.sqrt(2*np.pi*Gbar))
    f_o = fock[np.ix_(occupied, occupied)]
    f_v = fock[np.ix_(virtual, virtual)]
    Eg = float(np.min(np.diag(f_v) - (np.sum(np.abs(f_v), axis=1)-np.abs(np.diag(f_v)))) - np.max(np.diag(f_o) + (np.sum(np.abs(f_o), axis=1)-np.abs(np.diag(f_o)))))
    gap = float(np.linalg.eigvalsh(A).min())
    errors = {
        "Fourier_unitarity": float(np.linalg.norm(fourier.conj().T@fourier-np.eye(N))),
        "orbital_orthogonality": float(np.linalg.norm(orbital.T@orbital-np.eye(N))),
        "charged_overlap": float(np.linalg.norm(g*Z*green_columns.conj().T@green_columns-coulomb)),
        "neutral_Gram": float(np.linalg.norm(W.conj().T@W-J)),
        "source_formula": float(np.linalg.norm(b-b_formula)),
        "single_double_formula": float(np.linalg.norm(K-K_formula)),
        "HF_singles_residual": float(np.linalg.norm(b[:ns])),
    }
    result = {
        "status": "finite_algebra_check_only",
        "L": L, "primitive_sites": N, "cells": V, "full_fixed_electron_dimension": len(states),
        "SD_dimension": len(indices), "interaction_g": g, "inhomogeneous_p": p.tolist(),
        "errors": errors,
        "normalizations": {"G0_exact_finite_sum": G0, "G0_uniform_bound": Gbar, "Z_known_envelope": Z,
            "J_norm": float(np.linalg.norm(J, 2)), "alpha_J": alpha_J,
            "source_norm_squared_per_cell": float(b@b/V), "B0_squared": B0**2,
            "K_norm": float(np.linalg.norm(K, 2)), "alpha_K": alpha_K,
            "max_K_target_predecessors": int(np.sum(np.abs(K)>1e-11, axis=1).max()),
            "interaction_norm": interaction_norm, "Gamma_uniform_bound": Gamma,
            "Eg_Gershgorin": Eg, "certified_Delta": Eg-Gamma, "actual_minimum_eigenvalue": gap},
        "correlation_energy_per_cell": float(-b@np.linalg.solve(A,b)/V),
        "scope": "Checks finite Poisson/Gram/source/K algebra and conservative bound constants. No asymptotic, material, runtime or quantum advantage claim."
    }
    assert max(errors.values()) < 1e-10
    assert np.linalg.norm(J,2) <= alpha_J + 1e-12
    assert b@b/V <= B0**2 + 1e-12
    assert np.linalg.norm(K,2) <= alpha_K + 1e-12
    assert result["normalizations"]["max_K_target_predecessors"] <= 2
    assert interaction_norm <= Gamma and gap >= Eg-Gamma > 0
    (Path(__file__).parent / "results" / "poisson_algebra.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
