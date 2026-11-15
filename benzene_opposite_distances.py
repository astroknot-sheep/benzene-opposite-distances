#!/usr/bin/env python3
"""Distances between diametrically opposite carbon atoms of benzene rings in crystals.

Input : individual or concatenated CIF files from the BENZEN refcode family.
Output: CSV with columns  CSD_id, atom_id1, atom_id2, distance  (Angstrom).

How it works
------------
1. Read the cell, space group and the atoms of the asymmetric unit (gemmi).
2. Apply every symmetry operator and lattice translation, so molecules that
   sit on a symmetry element (benzene lies on an inversion centre, so the
   asymmetric unit holds only half a ring) are completed.
3. Build the C-C bond graph (bond if distance < 1.75 A) and find all 6-membered
   carbon rings that contain an atom of the asymmetric unit.
4. Keep rings that look like benzene: every ring atom has exactly two ring
   neighbours, mean ring bond length 1.30-1.50 A, ring planar to within 0.15 A.
5. For each ring write its 3 pairs of opposite (para) atoms and their
   Euclidean distance.

Atom labels come from the CIF. A symmetry-generated copy of an atom gets a
trailing apostrophe, e.g. C1 and C1' for the two ends of a diameter of a ring
on an inversion centre.
"""
import argparse
import csv
import itertools
import re
import sys
from pathlib import Path

import gemmi
import numpy as np

BOND_MAX = 1.75          # A, C-C bond cutoff
BOND_MIN = 0.40          # A, ignore coincident images
RING_BOND_RANGE = (1.30, 1.50)   # A, mean bond length of a benzene-like ring
PLANAR_TOL = 0.15        # A, max distance of a ring atom from the mean plane
MERGE_TOL = 0.15         # A, images closer than this are the same atom

# Deposition numbers from CCDC Access Structures for BENZEN through BENZEN28.
BENZEN_REFCODES = {
    "1108749": "BENZEN",
    "1108750": "BENZEN01",
    "1108751": "BENZEN02",
    "1108752": "BENZEN03",
    "1108753": "BENZEN04",
    "1108754": "BENZEN05",
    "1108755": "BENZEN06",
    "1108756": "BENZEN07",
    "1108757": "BENZEN08",
    "1108758": "BENZEN09",
    "1108759": "BENZEN10",
    "298305": "BENZEN11",
    "298306": "BENZEN12",
    "298307": "BENZEN13",
    "682617": "BENZEN14",
    "757059": "BENZEN15",
    "757060": "BENZEN16",
    "757061": "BENZEN17",
    "725244": "BENZEN18",
    "1454032": "BENZEN19",
    "1423904": "BENZEN20",
    "1581798": "BENZEN21",
    "1581824": "BENZEN22",
    "1579553": "BENZEN23",
    "1579564": "BENZEN24",
    "1843296": "BENZEN25",
    "1843297": "BENZEN26",
    "1913527": "BENZEN27",
    "2201163": "BENZEN28",
}

_OFFSETS = np.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float)


def natural_key(label):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", label)]


class CarbonGraph:
    """Carbon atoms of a crystal, generated lazily from the asymmetric unit."""

    def __init__(self, sm):
        self.M = np.array(sm.cell.orth.mat.tolist(), dtype=float)
        sg = sm.spacegroup or gemmi.SpaceGroup("P 1")
        ops = list(sg.operations())
        sites = [s for s in sm.sites
                 if s.element.name == "C" and s.disorder_group <= 1]
        self.labels = [s.label for s in sites]
        frac = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in sites]).reshape(-1, 3)

        # every symmetry image of every asymmetric-unit carbon (fractional coords)
        self.img, self.img_j = [], []
        for op in ops:
            for j, f in enumerate(frac):
                self.img.append(op.apply_to_xyz([float(x) for x in f]))
                self.img_j.append(j)
        self.img = np.array(self.img).reshape(-1, 3)

        # nodes = distinct atoms in the (unwrapped) crystal; the first len(sites)
        # are the asymmetric-unit atoms exactly as listed in the CIF ("home" atoms)
        self.pos = [f.copy() for f in frac]
        self.label = list(self.labels)
        self.home = [True] * len(sites)
        self.n_home = len(sites)
        self._nb = {}

    def _node_at(self, p, j):
        d = (np.array(self.pos) - p) @ self.M.T
        dist = np.linalg.norm(d, axis=1)
        i = int(np.argmin(dist))
        if dist[i] < MERGE_TOL:
            return i
        self.pos.append(p)
        self.label.append(self.labels[j])
        self.home.append(False)
        return len(self.pos) - 1

    def neighbours(self, i):
        if i in self._nb:
            return self._nb[i]
        p = self.pos[i]
        d = self.img - p
        d = d - np.round(d)                       # nearest lattice image
        cand = d[None, :, :] + _OFFSETS[:, None, :]
        dist = np.linalg.norm(cand @ self.M.T, axis=2)
        o_idx, n_idx = np.nonzero((dist < BOND_MAX) & (dist > BOND_MIN))
        out = {self._node_at(p + cand[o, n], self.img_j[n]) for o, n in zip(o_idx, n_idx)}
        out.discard(i)
        self._nb[i] = sorted(out)
        return self._nb[i]

    def six_rings(self, start):
        found = {}
        path = [start]

        def dfs():
            last = path[-1]
            if len(path) == 6:
                if start in self.neighbours(last):
                    found.setdefault(frozenset(path), list(path))
                return
            for m in self.neighbours(last):
                if m not in path:
                    path.append(m)
                    dfs()
                    path.pop()

        dfs()
        return list(found.values())

    def cart(self, i):
        return self.M @ self.pos[i]

    def is_benzene_like(self, cyc):
        members = set(cyc)
        if any(len(members & set(self.neighbours(n))) != 2 for n in cyc):
            return False                           # chords / fused or bridged
        pts = np.array([self.cart(n) for n in cyc])
        bonds = [np.linalg.norm(pts[k] - pts[(k + 1) % 6]) for k in range(6)]
        if not RING_BOND_RANGE[0] <= float(np.mean(bonds)) <= RING_BOND_RANGE[1]:
            return False
        centred = pts - pts.mean(axis=0)
        normal = np.linalg.svd(centred)[2][2]
        return float(np.abs(centred @ normal).max()) <= PLANAR_TOL

    def name(self, i):
        return self.label[i] if self.home[i] else self.label[i] + "'"


def refcode_of(block, path, single_block):
    value = block.find_value("_database_code_CSD")
    if value and not gemmi.cif.is_null(value):
        return gemmi.cif.as_string(value).strip().upper()

    for tag in ("_database_code_depnum_ccdc_archive", "_database_code_depnum_ccdc"):
        value = block.find_value(tag) or ""
        match = re.fullmatch(r"(?:CCDC\s+)?(\d+)", value.strip("'\" "), re.IGNORECASE)
        if match and match[1] in BENZEN_REFCODES:
            return BENZEN_REFCODES[match[1]]

    candidates = [block.name, Path(path).stem] if single_block else [block.name]
    for candidate in candidates:
        if re.fullmatch(r"BENZEN(?:[0-9]{2})?", candidate, re.IGNORECASE):
            return candidate.upper()
    raise ValueError(f"data_{block.name}: unable to identify the CSD refcode")


def rows_for_block(block, csd_id):
    sm = gemmi.make_small_structure_from_block(block)
    g = CarbonGraph(sm)
    rings = {}
    for start in range(g.n_home):
        for cyc in g.six_rings(start):
            if g.is_benzene_like(cyc):
                rings.setdefault(frozenset(cyc), cyc)
    rows = []
    for cyc in rings.values():
        for k in range(3):
            a, b = cyc[k], cyc[k + 3]
            dist = float(np.linalg.norm(g.cart(a) - g.cart(b)))
            n1, n2 = sorted((g.name(a), g.name(b)), key=natural_key)
            rows.append((csd_id, n1, n2, dist))
    rows.sort(key=lambda r: (natural_key(r[1]), natural_key(r[2])))
    return rows


def rows_for_cif(path, warn=None):
    """All (CSD_id, atom_id1, atom_id2, distance) rows for one CIF file."""
    doc = gemmi.cif.read_file(str(path))
    rows = []
    for block in doc:
        csd_id = refcode_of(block, path, single_block=len(doc) == 1)
        block_rows = rows_for_block(block, csd_id)
        if not block_rows and warn is not None:
            warn(f"{csd_id}: no measurable benzene ring")
        rows.extend(block_rows)
    return rows


def collect_files(inputs):
    files = []
    for item in inputs:
        p = Path(item)
        files += sorted(p.glob("*.cif")) if p.is_dir() else [p]
    return files


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("inputs", nargs="+", help="CIF files and/or folders of CIF files")
    ap.add_argument("-o", "--output", default="benzene_distances.csv")
    ap.add_argument("--decimals", type=int, default=4)
    args = ap.parse_args(argv)

    all_rows, problems = [], []
    files = collect_files(args.inputs)
    if not files:
        ap.error("no CIF files found")
    if not 0 <= args.decimals <= 15:
        ap.error("--decimals must be between 0 and 15")

    read_failed = False
    for path in files:
        try:
            rows = rows_for_cif(path, warn=problems.append)
        except Exception as exc:                   # keep going; report at the end
            read_failed = True
            problems.append(f"{path.name}: could not read ({exc})")
            continue
        all_rows += rows
        print(f"{path.name}: {len(rows)} pairs", file=sys.stderr)

    all_rows.sort(key=lambda r: (natural_key(r[0]), natural_key(r[1]), natural_key(r[2])))
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["CSD_id", "atom_id1", "atom_id2", "distance"])
        for csd_id, a1, a2, d in all_rows:
            w.writerow([csd_id, a1, a2, f"{d:.{args.decimals}f}"])

    print(f"\nWrote {len(all_rows)} rows from {len(files)} files to {args.output}", file=sys.stderr)
    for msg in problems:
        print("WARNING:", msg, file=sys.stderr)
    return int(read_failed or not all_rows)


if __name__ == "__main__":
    sys.exit(main())
