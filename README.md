# UMDFP — Universal Molecular Description File Processor

UMDFP is a cheminformatics pipeline that converts **SMILES strings** into a compact, binary-searchable **Universal Molecular Format (UMF)** database that can be instantly exported to common molecular file formats (Mol2, SDF, PDBQT, SMILES).

---

## Overview

The pipeline has two main stages:

```
SMILES file  ──(smiles2umd)──►  .umd file  ──(umfgen)──►  .umf + .umfp files
                                                                   │
                          mol2 / SDF / PDBQT / SMILES  ◄──(umfdump / umfparser)──┘
```

| Stage | Tool | Description |
|-------|------|-------------|
| 1 | `smiles2umd` (Python) | Converts a `.smi` file into a human-readable UMD text file. Embeds 3-D coordinates, minimises geometry, and assigns partial charges. |
| 2 | `umfgen` (C++) | Compresses one or more UMD files into a single binary UMF file with fast look-ahead indices. |
| Export | `umfdump` / `umfparser` (C++) | Convert UMF molecules back to standard cheminformatics formats. |

---

## File Formats

### UMD — Universal Molecular Description (text)

A human-readable, semicolon-commented text format.  Each molecule is enclosed in `START` / `END` markers:

```
START
; Name
Aspirin
; SMILES String
CC(=O)Oc1ccccc1C(=O)O
; # Atoms
21
; # Bonds
21
; Atoms
; Index Element X Y Z Q Hyb Aromatic FormalQ
0 C  1.2345 -0.6789 0.0000 0.123 3 0 0
...
; Bonds
; Index At1 At2 2*Order
0 0 1 2
...
END
```

Atom hybridisation is encoded as an integer (0 = S, 1 = SP, 2 = SP2, 3 = SP3, 4 = SP3D, 5 = SP3D2).  
Bond order is stored as `2 × order` (single = 2, double = 4, triple = 6, aromatic = 3).

### UMF — Universal Molecular Format (binary)

A compact binary format designed for databases of **over one billion molecules**.

* Prefixed with a 32-byte file signature (`UMF_SIGNATURE`).
* Each molecule block carries a per-molecule header with:
  * the block's byte size (for O(1) skipping),
  * a 256-molecule look-ahead size (for fast batch seeking), and
  * a 65 536-molecule look-ahead size (for very large jumps).
* A companion `.umfp` pointer file enables O(log *n*) lookup by molecule name.

---

## Tools

### `smiles2umd` — SMILES → UMD

Python script located at `src/smiles2umd`.

**Requirements:** Python 3, RDKit, OpenFF Toolkit, NumPy.

```
Usage: smiles2umd <input.smi> [options]

Options:
  -o, --output FILE       Output UMD file (default: output.umd)
  -c, --cpus N            Number of parallel worker threads (default: 1)
  -b, --batch_size N      SMILES to load per batch (default: 4096)
  -e, --error FILE        Log file for failed molecules (default: error.log)
  -s, --segment N         Re-run a specific CPU segment (advanced use)
  -f, --force-overwrite   Overwrite existing output files
```

**Example:**
```bash
python src/smiles2umd library.smi -o library.umd -c 8
```

The script:
1. Parses SMILES with **RDKit**.
2. Adds explicit hydrogens and embeds a 3-D conformer using **ETKDGv2**.
3. Minimises geometry with the **MMFF94** force field.
4. Assigns **Gasteiger** partial charges via OpenFF.
5. Writes each molecule to the UMD file.

Molecules that fail (undefined stereochemistry, unsupported radicals, valency errors) are silently skipped and their names are written to the error log.

---

### `umfgen` — UMD → UMF

Compresses a UMD file into a binary UMF database and immediately validates it.

```bash
# Build
cd src && make umfgen

# Run
./umfgen <input.umd> [output.umf]   # default output: test.umf
```

After writing the UMF file it automatically builds the `.umfp` pointer file and performs a read-back validation pass.

---

### `umfparser` — Query UMF by name or index

```bash
./umfparser <file.umf> <name|index> [-f mol2|sdf|pdbqt|smi] [-o output_file] [-s]
```

| Flag | Description |
|------|-------------|
| `-f` | Output format: `mol2` (default), `sdf`, `pdbqt`, `smi` |
| `-o` | Write output to a file instead of stdout |
| `-s` | Write each molecule to a separate file (requires `-o`) |

The second argument can be:
* A **molecule name** (string) — resolved via the `.umfp` pointer file.
* A **0-based integer index** — resolved by seeking through look-ahead indices.
* A **path to a text file** containing one molecule name per line.

---

### `umfparser-live` — Interactive query shell

An interactive REPL that keeps the UMF reader open between queries.

```bash
./umfparser-live <file.umf>
# Enter query (Ligand name or index, or 'exit' to quit):
```

---

### `umfdump` — Export all molecules

Dumps every molecule in a UMF file to a single output file (or one file per molecule with `-s`).

```bash
./umfdump <file.umf> [output.mol2] [-f mol2|sdf|pdbqt|smi] [-s]
```

---

### `umfchecker` — Validate a UMF file

Reads through every molecule in the file and reports any corruption.

```bash
./umfchecker <file.umf>
```

---

## Building the C++ Tools

Requires **g++** with C++17 support.

```bash
cd src
make          # builds all tools: umfgen, umfparser, umfparser-live, umfchecker, umfdump
make clean    # removes object files and executables
```

Pre-built binaries for Linux are also available in the `bin/` directory.

---

## Supported Output Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Tripos Mol2 | `.mol2` | Full atom-type and charge annotation |
| SDF / MDL Molfile | `.sdf` | Includes formal charge (`M  CHG`) and SMILES tag |
| PDBQT | `.pdbqt` | AutoDock-compatible; includes `ROOT`/`BRANCH`/`TORSDOF` records |
| SMILES | `.smi` | One line per molecule: `<SMILES> <name>` |

---

## Supported Elements

H, Li, Be, B, C, N, O, F, Mg, Si, P, S, Cl, Fe, Se, Br, Zn, Pd, I, Pt.

---

## Repository Layout

```
UMDFP/
├── src/
│   ├── smiles2umd          # Python: SMILES → UMD converter
│   ├── General.py          # Python utilities (SequentialSMILESLoader, Mol2FileLoaded, …)
│   ├── UMDFP.h             # Core C++ data structures (UMDAtom, UMDBond, UMDMolecule,
│   │                       #   UMDReader, UMFWriter, UMFReader, …)
│   ├── UMDBabel.h          # Format converters (Mol2Format, SDFFormat, PDBQTFormat, …)
│   ├── UMFHelpers.h        # Topology utilities (DFS traversal, ring detection, …)
│   ├── SearchableFP.h      # Fingerprint I/O and Tanimoto / cosine similarity
│   ├── UMDcompressor.cpp   # → umfgen
│   ├── UMFParser.cpp       # → umfparser-live
│   ├── UMFQuery.cpp        # → umfparser
│   ├── UMFChecker.cpp      # → umfchecker
│   ├── UMFDumper.cpp       # → umfdump
│   └── Makefile
├── bin/                    # Pre-built Linux binaries
├── UMFTests.ipynb          # Jupyter notebook demonstrating batch SMILES→UMD conversion
└── LICENSE
```

---

## License

See [LICENSE](LICENSE) for details.
