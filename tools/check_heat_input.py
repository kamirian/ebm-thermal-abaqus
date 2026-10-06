"""Integrate the DFLUX beam flux over the mesh and compare with the nominal absorbed power.

    python tools/check_heat_input.py

For first-order heat-transfer elements Abaqus evaluates a body flux at the nodes and
lumps it with nodal volumes. This script does the same for the POWDER_LAYER elements,
with the flux formula of model/dflux.for, and compares the total with
ETA * U * IB / 2 (half model).
"""

from pathlib import Path

import numpy as np

MODEL = Path(__file__).resolve().parents[1] / "model"

# Parameters of model/dflux.for (mm, s, mW)
PI, PHI, ETA, U, IB, V, S, Z0 = 3.1415926, 0.55, 0.9, 60.0e3, 6.7, 632.6, 0.062, 10.07


def read_block(lines, header):
    i = next(k for k, l in enumerate(lines) if l.startswith(header)) + 1
    rows = []
    while i < len(lines) and not lines[i].startswith("*"):
        rows.append([v for v in lines[i].split(",") if v.strip()])
        i += 1
    return rows


def main():
    mesh = (MODEL / "ebm_mesh.inp").read_text().splitlines()
    sets = (MODEL / "ebm_sets.inp").read_text().splitlines()
    nodes = np.array(read_block(mesh, "*Node"), dtype=float)
    X = np.zeros((int(nodes[:, 0].max()) + 1, 3))
    X[nodes[:, 0].astype(int)] = nodes[:, 1:]
    elems = np.array(read_block(mesh, "*Element"), dtype=int)
    first, last, step = (int(v) for v in read_block(sets, "*Elset, elset=POWDER_LAYER")[0])
    layer = elems[np.isin(elems[:, 0], np.arange(first, last + 1, step))]
    conn = layer[:, 1:]

    # Nodal volumes: integral of each shape function over the element (2x2x2 Gauss)
    nat = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    g = 1 / np.sqrt(3)
    vol = np.zeros(len(X))
    xe = X[conn]
    for p in np.array([[a, b, c] for a in (-g, g) for b in (-g, g) for c in (-g, g)]):
        N = 0.125 * np.prod(1 + nat * p, axis=1)
        dN = np.empty((8, 3))
        for k in range(3):
            o = [j for j in range(3) if j != k]
            dN[:, k] = 0.125 * nat[:, k] * (1 + nat[:, o[0]] * p[o[0]]) * (1 + nat[:, o[1]] * p[o[1]])
        detJ = np.abs(np.linalg.det(np.einsum("ak,eaj->ekj", dN, xe)))
        np.add.at(vol, conn, detJ[:, None] * N[None, :])

    used = np.unique(conn)
    x, y, z = X[used].T
    w = vol[used]
    d = Z0 - z
    inside = (d >= 0) & (d <= S)
    iz = (-2.25 * (d / S) ** 2 + 1.5 * (d / S) + 0.75) / 0.75
    nominal = 0.5 * ETA * U * IB
    print(f"powder-layer elements: {len(layer)}, nodes: {len(used)}")
    print(f"nominal absorbed power (half model): {nominal / 1e3:.2f} W")
    for t in (0.005, 0.01, 0.02):
        x0 = -PHI + V * t
        hs = 2 * U * IB / (PI * PHI ** 2) * np.exp(-2 * ((x - x0) ** 2 + y ** 2) / PHI ** 2)
        power = np.sum(np.where(inside, ETA * hs * iz / S, 0.0) * w)
        print(f"t = {t:5.3f} s, beam at x = {x0:6.3f} mm: absorbed {power / 1e3:7.2f} W "
              f"({100 * power / nominal:.2f}% of nominal)")


if __name__ == "__main__":
    main()
