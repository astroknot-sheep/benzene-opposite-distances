# Opposite-atom distances in benzene rings (CSD refcode family BENZEN)

Reads CIF files of the BENZEN family and writes a CSV with one row per pair of
diametrically opposite ring carbons:

| column     | meaning                                              |
|------------|------------------------------------------------------|
| `CSD_id`   | refcode, e.g. `BENZEN01`                             |
| `atom_id1` | label of one carbon, e.g. `C1`                       |
| `atom_id2` | label of the diametrically opposite carbon           |
| `distance` | Euclidean distance between them, in Angstrom         |

A benzene ring has three such pairs (1-4, 2-5, 3-6), so a ring gives three rows.
A regular ring with C-C = 1.39 A gives about 2.78 A.

## Usage

```bash
python -m pip install -r requirements.txt
# put the downloaded CIFs (BENZEN.cif, BENZEN01.cif, ...) in ./cifs
python benzene_opposite_distances.py cifs -o benzene_distances.csv
# Or pass a combined download directly
python benzene_opposite_distances.py BENZEN-family.cif -o benzene_distances.csv
python tests/test_synthetic.py          # self-check on synthetic crystals
```

The script reports CIF blocks with no measurable ring on standard error.
It uses the CIF's CSD refcode when available, then checks the deposition number
against the BENZEN family mapping. In a combined CIF, each block is identified
independently.

## Getting the CIFs

Download the BENZEN family from [CCDC Access Structures](https://www.ccdc.cam.ac.uk/structures/)
as individual CIFs or use **Download Selected** for a combined CIF.
Keep downloads local under the [Access Structures terms](https://www.ccdc.cam.ac.uk/access-structures-terms).

## How it works

1. Read cell, space group and asymmetric-unit atoms with
   [gemmi](https://gemmi.readthedocs.io/).
2. Apply all symmetry operators and lattice translations. This matters for
   benzene: the molecule lies on an inversion centre, so the asymmetric unit has
   only half a ring (C1, C2, C3) and the opposite atoms are symmetry copies.
3. Build the C-C bond graph (bond if closer than 1.75 A) and find 6-membered
   carbon rings through the asymmetric-unit atoms.
4. Keep benzene-like rings: each atom has exactly two ring neighbours, mean ring
   bond 1.30-1.50 A, ring planar to within 0.15 A. This drops chair cyclohexane,
   and rings containing N, O, etc. are ignored because only carbon is used.
5. Report the three opposite pairs of each ring.

### Label convention

Labels are those of the CIF. A symmetry-generated copy gets a trailing
apostrophe, so for a ring on an inversion centre the rows read `C1, C1'`,
`C2, C2'`, `C3, C3'`. For a ring whose six atoms are all listed, the rows read
`C1, C4`, `C2, C5`, `C3, C6`.

## Limits

* Atoms in minor disorder components (disorder group 2 and higher) are skipped.
* Only the cutoffs in the constants at the top of the script decide what counts as
  a bond and as benzene-like; adjust them if a structure is missed.
## Validation

The synthetic tests cover inversion symmetry, a ring spanning a unit-cell
boundary, P2₁/c, two independent rings, two non-benzene controls, and refcode
resolution in a combined CIF.

Tested end to end on 29 downloaded CCDC records with Python 3.11, Gemmi 0.7.5,
and NumPy 2.4.6. The run yielded 69 distances across 23 entries. BENZEN05,
BENZEN08–BENZEN10, BENZEN27, and BENZEN28 lacked atomic coordinates. The three
BENZEN distances agreed with Mercury to its displayed precision (2.782,
2.789, and 2.768 Å).

## Next step

The CCDC CSD Python API (<https://downloads.ccdc.cam.ac.uk/documentation/API/index.html>)
can read entries straight from an installed CSD by refcode, which avoids
downloading CIFs. It needs a CSD installation and licence.
