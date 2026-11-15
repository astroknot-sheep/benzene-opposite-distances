"""Checks on synthetic crystals with a known answer (no CSD data needed).

Run:  python tests/test_synthetic.py        (or: pytest tests)

A regular benzene ring with C-C = 1.39 A has para distance 2 * 1.39 = 2.78 A.
"""
import math
import sys
import tempfile
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benzene_opposite_distances import rows_for_cif  # noqa: E402

CC = 1.39
CELL = (7.1, 9.2, 6.8, 88.0, 95.0, 101.0)
P1 = ["x,y,z"]
P21C = ["x,y,z", "-x,y+1/2,-z+1/2", "-x,-y,-z", "x,-y+1/2,z+1/2"]
PM1 = ["x,y,z", "-x,-y,-z"]


def ring_atoms(centre_frac, cell, seed, n=6, elements="CCCCCC", radius=CC, dz=0.0):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    u, v, w = q[:, 0], q[:, 1], q[:, 2]
    c = np.array(cell.orthogonalize(gemmi.Fractional(*centre_frac)).tolist())
    out = []
    for k in range(n):
        th = k * math.pi / 3
        base = c + radius * (math.cos(th) * u + math.sin(th) * v)
        base = base + (dz if k % 2 else -dz) * w
        hyd = c + (radius + 1.08) * (math.cos(th) * u + math.sin(th) * v)
        out.append((elements[k], base, hyd))
    return out


def write_cif(path, code, ops, atoms, cell_params=CELL):
    lines = [f"data_{code}", f"_database_code_CSD {code}"]
    for tag, val in zip(["a", "b", "c"], cell_params[:3]):
        lines.append(f"_cell_length_{tag} {val}")
    for tag, val in zip(["alpha", "beta", "gamma"], cell_params[3:]):
        lines.append(f"_cell_angle_{tag} {val}")
    lines += ["loop_", "_space_group_symop_operation_xyz"] + [f"'{o}'" for o in ops]
    lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol",
              "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
    cell = gemmi.UnitCell(*cell_params)
    for label, el, xyz in atoms:
        f = cell.fractionalize(gemmi.Position(*xyz))
        lines.append(f"{label} {el} {f.x:.6f} {f.y:.6f} {f.z:.6f}")
    Path(path).write_text("\n".join(lines) + "\n")


def build(tmp, name, ops, centre, keep, seed, first=1, **kw):
    cell = gemmi.UnitCell(*CELL)
    ring = ring_atoms(centre, cell, seed, **kw)
    atoms = []
    for k in keep:
        el, c_xyz, h_xyz = ring[k]
        atoms.append((f"{el}{first + k}", el, c_xyz))
        atoms.append((f"H{first + k}", "H", h_xyz))
    p = Path(tmp) / f"{name}.cif"
    write_cif(p, name, ops, atoms)
    return p


def check(rows, expected_pairs):
    assert len(rows) == expected_pairs, rows
    for _, _, _, d in rows:
        assert abs(d - 2 * CC) < 1e-3, d


def test_half_molecule_on_inversion_centre():
    with tempfile.TemporaryDirectory() as tmp:
        p = build(tmp, "HALF01", PM1, (0.0, 0.5, 0.0), keep=[0, 1, 2], seed=1)
        rows = rows_for_cif(p)
        check(rows, 3)
        assert [(r[1], r[2]) for r in rows] == [("C1", "C1'"), ("C2", "C2'"), ("C3", "C3'")]
        assert all(r[0] == "HALF01" for r in rows)


def test_full_ring_straddling_cell_boundary():
    with tempfile.TemporaryDirectory() as tmp:
        p = build(tmp, "FULL", P1, (0.98, 0.03, 0.99), keep=range(6), seed=2)
        rows = rows_for_cif(p)
        check(rows, 3)
        assert [(r[1], r[2]) for r in rows] == [("C1", "C4"), ("C2", "C5"), ("C3", "C6")]


def test_general_position_in_P21c_gives_one_ring():
    with tempfile.TemporaryDirectory() as tmp:
        p = build(tmp, "P21C", P21C, (0.3, 0.2, 0.15), keep=range(6), seed=3)
        check(rows_for_cif(p), 3)


def test_two_independent_rings():
    with tempfile.TemporaryDirectory() as tmp:
        cell = gemmi.UnitCell(*CELL)
        atoms = []
        for first, centre, seed in ((1, (0.2, 0.2, 0.2), 4), (7, (0.7, 0.7, 0.7), 5)):
            for k, (el, c_xyz, h_xyz) in enumerate(ring_atoms(centre, cell, seed)):
                atoms.append((f"C{first + k}", "C", c_xyz))
                atoms.append((f"H{first + k}", "H", h_xyz))
        p = Path(tmp) / "TWO.cif"
        write_cif(p, "TWO", P1, atoms)
        check(rows_for_cif(p), 6)


def test_chair_cyclohexane_is_ignored():
    with tempfile.TemporaryDirectory() as tmp:
        p = build(tmp, "CHAIR", P1, (0.4, 0.4, 0.4), keep=range(6), seed=6,
                  radius=1.46, dz=0.25)
        assert rows_for_cif(p) == []


def test_pyridine_is_ignored():
    with tempfile.TemporaryDirectory() as tmp:
        p = build(tmp, "PYR", P1, (0.4, 0.4, 0.4), keep=range(6), seed=7,
                  elements="CCCCCN")
        assert rows_for_cif(p) == []


def test_concatenated_cifs_resolve_each_deposition():
    with tempfile.TemporaryDirectory() as tmp:
        blocks = []
        for code, deposition in [("BENZEN11", 298305), ("BENZEN12", 298306)]:
            path = build(tmp, code, P1, (0.4, 0.4, 0.4), keep=range(6), seed=7)
            blocks.append(path.read_text().replace(
                f"_database_code_CSD {code}",
                f"_database_code_depnum_ccdc_archive 'CCDC {deposition}'",
            ).replace(f"data_{code}", f"data_deposition_{deposition}"))
        batch = Path(tmp) / "bundle.cif"
        batch.write_text("\n".join(blocks))
        rows = rows_for_cif(batch)
        check(rows, 6)
        assert [row[0] for row in rows] == ["BENZEN11"] * 3 + ["BENZEN12"] * 3


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            try:
                fn()
                print("PASS", name)
            except Exception as exc:
                failed += 1
                print("FAIL", name, repr(exc))
    sys.exit(1 if failed else 0)
