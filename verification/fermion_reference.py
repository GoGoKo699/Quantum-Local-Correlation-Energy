#!/usr/bin/env python3
"""Independent exact second-quantized fermion reference builder.

Enumerates the full fixed-electron-number sector. Used only for tiny algebra
checks; it is not the classical benchmark algorithm.
"""
from __future__ import annotations
import itertools
import numpy as np

def act(state: int, operations: list[tuple[str, int]]) -> tuple[int, int] | None:
    """Apply operations in right-to-left operator order; return state and sign."""
    sign = 1
    for kind, orbital in operations:
        mask = 1 << orbital
        occupied = bool(state & mask)
        if (kind == "ann" and not occupied) or (kind == "cre" and occupied):
            return None
        if (state & (mask - 1)).bit_count() & 1:
            sign = -sign
        state = state ^ mask
    return state, sign


def excitation_operations(holes: tuple[int, ...], particles: tuple[int, ...]):
    return [("ann", i) for i in holes] + [("cre", a) for a in reversed(particles)]


def hamiltonian(model):
    norb, nelec = model["norb"], model["nelec"]
    states = [sum(1 << i for i in c) for c in itertools.combinations(range(norb), nelec)]
    position = {state: j for j, state in enumerate(states)}
    h, eri = model["one_body"], model["eri"]
    matrix = np.zeros((len(states), len(states)))
    for col, state in enumerate(states):
        occ = [i for i in range(norb) if state & (1 << i)]
        for q in occ:
            after = state ^ (1 << q)
            for p in range(norb):
                if after & (1 << p):
                    continue
                transformed = act(state, [("ann", q), ("cre", p)])
                assert transformed is not None
                output, sign = transformed
                matrix[position[output], col] += h[p, q] * sign
        for q, s in itertools.combinations(occ, 2):
            after = state ^ (1 << q) ^ (1 << s)
            empty = [p for p in range(norb) if not after & (1 << p)]
            for p, r in itertools.combinations(empty, 2):
                transformed = act(state, [("ann", q), ("ann", s), ("cre", r), ("cre", p)])
                assert transformed is not None
                output, sign = transformed
                element = eri[p, q, r, s] - eri[p, s, r, q]
                matrix[position[output], col] += element * sign
    return states, position, matrix
