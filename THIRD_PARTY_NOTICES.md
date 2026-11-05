# Third-party notices

OpenGenEdit is Copyright (C) 2026 Antoine Champie and is licensed under the GNU Affero General Public License, version 3
([LICENSE](LICENSE)). It contains, links to, or is built with the third-party components below, each under its own licence.
Those licences are, as far as we can tell, compatible with distributing OpenGenEdit under the AGPL-3.0; this file is a good-faith
list of notices, not legal advice. If you find a component that is missing or wrong here, please open an issue.

How to read the "where" column: *bundled* = its code or data is inside this repository or inside the Windows `OpenGenEdit.exe`;
*linked* = comes from a package installed at build time and ends up in the exe; *build* = used to build, not shipped;
*test* = used only by the tests, not shipped.

## Main components

| Component | Version | Licence | URL | Where |
|---|---|---|---|---|
| PySide6-Essentials / shiboken6 (Qt for Python), with the Qt 6 libraries | 6.12.0 | LGPL-3.0 (also offered under GPL-2.0 / GPL-3.0; OpenGenEdit uses it under the LGPL-3.0) | https://doc.qt.io/qtforpython-6/ | linked, in the exe |
| PyO3 | 0.23.5 | MIT or Apache-2.0 | https://github.com/pyo3/pyo3 | linked (Rust to Python glue) |
| maturin | 1.15.0 (1.7 or newer, below 2) | MIT or Apache-2.0 | https://github.com/PyO3/maturin | build |
| gb-io (patched copy) | 0.9.0 | MIT, Copyright (c) 2018 David Leslie | https://github.com/dlesl/gb-io | bundled in `core/vendor/gb-io/` (licence: `core/vendor/gb-io/LICENSE`) |
| noodles (and its sub-crates `-bam`, `-sam`, `-core`, `-bgzf`, `-csi`, `-vcf`, `-tabix`) | 0.118.0 | MIT | https://github.com/zaeleus/noodles | linked |
| quick-xml | 0.42.0 | MIT | https://github.com/tafia/quick-xml | linked |
| 3Dmol.js | 2.4.2 | BSD-3-Clause (it incorporates code from GLmol, Three.js and jQuery) | https://3dmol.org/ | bundled in `app/vendor/3Dmol-min.js` (licence: `app/vendor/3Dmol.LICENSE`) |
| PyInstaller | 6.22.3 | GPL-2.0-or-later with a special exception that allows the programs it builds to be distributed under other licences | https://pyinstaller.org/ | build (its bootloader is inside the exe) |
| Python (CPython) and the libraries it ships with (OpenSSL, zlib, libffi, expat, ...) | 3.13 | PSF-2.0 and the licences of those libraries | https://www.python.org/ | inside the exe |
| DM Sans (font) | static Medium and Bold instances of the Google Fonts release | SIL Open Font License 1.1, Copyright 2014 The DM Sans Project Authors | https://github.com/googlefonts/dm-fonts | build: in `packaging/fonts/` (licence: `packaging/fonts/OFL.txt`), used by `packaging/make_icon.py`; the wordmark images carry outlines of its letters, the font itself is not in the exe |
| pytest, Biopython, pyflakes | current | MIT; Biopython licence (BSD-3-Clause style); MIT | https://pytest.org/ , https://biopython.org/ , https://github.com/PyCQA/pyflakes | test |

Versions are those used when this file was written (`core/Cargo.lock` and the build environment); `core/Cargo.lock` is the exact list for the
Rust side. Check the licence text of a component in its own repository or package for the exact wording.

### 3Dmol.js (BSD-3-Clause)

The 3D structure page embeds 3Dmol.js. Its full copyright and licence notice, which must accompany it, is in
[`app/vendor/3Dmol.LICENSE`](app/vendor/3Dmol.LICENSE) (inside the exe: `vendor\3Dmol.LICENSE`). Copyright (c) 2014, University of Pittsburgh and
contributors. Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. Bioinformatics 31, 1322-1324 (2015),
doi:10.1093/bioinformatics/btu829.

### gb-io (MIT), patched

`core/vendor/gb-io` is a copy of gb-io 0.9.0 with a few changes of ours (the GenBank writer's line wrapping, an optional blank after `ORIGIN`
so an empty record can be re-read, the gap letter `-`, and the end of a location). Every change is marked "OpenGenEdit patch" in the source.
The original MIT licence and copyright (c) 2018 David Leslie stay in `core/vendor/gb-io/LICENSE`. The patch is built in through
`[patch.crates-io]` in `core/Cargo.toml`.

## Qt / PySide6 (LGPL): what you can do, and how to get the source

The Windows exe is built with PyInstaller and contains Qt 6 and PySide6, which are free software under the **GNU Lesser General Public
License version 3** (Qt for Python is also available under GPL-2.0 and GPL-3.0; the Qt libraries it bundles carry further notices for
the third-party code inside Qt, such as FreeType, HarfBuzz, libpng, zlib and ICU, listed on the Qt page below).

- Licence texts: LGPL-3.0 https://www.gnu.org/licenses/lgpl-3.0.html , GPL-3.0 https://www.gnu.org/licenses/gpl-3.0.html , and the list of
  licences of what Qt contains: https://doc.qt.io/qt-6/licenses-used-in-qt.html . The `PySide6` and `shiboken6` Python packages carry
  their own copies in their `dist-info` folders when installed with pip.
- Source code of Qt and PySide6 (we do not modify them): Qt https://download.qt.io/official_releases/qt/ ; Qt for Python
  https://code.qt.io/cgit/pyside/pyside-setup.git/ (tag `v6.12.0`) and https://pypi.org/project/PySide6-Essentials/#files (source and wheels).
- Using a different Qt: the LGPL lets you replace the Qt libraries with your own build of Qt/PySide6 and run OpenGenEdit with them. Each release
  has two Windows downloads. **`OpenGenEdit-folder.zip`** is the program as a folder: Qt is plain files in `_internal\PySide6` and
  `_internal\shiboken6`, and the `README.txt` inside says how to swap them for another PySide6 version. `OpenGenEdit.exe` is a PyInstaller *one-file*
  exe that packs those libraries inside itself, so you cannot swap a file in it; for that one, **rebuild OpenGenEdit from its source** (public,
  AGPL-3.0) against any PySide6 version you like: `pip install PySide6-Essentials==<version>` and follow the build steps in the README.
  Help > About shows the Qt version in use.

## Rust crates

OpenGenEdit's core (`core/`) is Rust and uses the crates below (the full dependency tree for the Windows build, read from `cargo metadata`
on `core/Cargo.lock`: 71 crates besides our own). Many are only used while compiling (procedural macros, build scripts) and are not in the
final binary. MIT, Apache-2.0 and the like ask that copyright and licence notices are kept; each crate's repository (last column) holds its
licence file and copyright line, and `cargo` stores them in `~/.cargo/registry` after a build. "A OR B" means the crate may be used under either.

| Crate | Version | Licence | Repository |
|---|---|---|---|
| autocfg | 1.5.1 | Apache-2.0 OR MIT | https://github.com/cuviper/autocfg |
| bit-vec | 0.10.1 | Apache-2.0 OR MIT | https://github.com/contain-rs/bit-vec |
| bitflags | 2.13.2 | MIT OR Apache-2.0 | https://github.com/bitflags/bitflags |
| borsh | 1.8.1 | MIT OR Apache-2.0 | https://github.com/near/borsh-rs |
| borsh-derive | 1.8.1 | Apache-2.0 | https://github.com/near/borsh-rs |
| bstr | 1.13.1 | MIT OR Apache-2.0 | https://github.com/BurntSushi/bstr |
| bytes | 1.12.1 | MIT | https://github.com/tokio-rs/bytes |
| cfg-if | 1.0.5 | MIT OR Apache-2.0 | https://github.com/rust-lang/cfg-if |
| cfg_aliases | 0.2.2 | MIT | https://github.com/katharostech/cfg_aliases |
| circular | 0.3.0 | MIT | https://github.com/sozu-proxy/circular |
| crossbeam-channel | 0.5.17 | MIT OR Apache-2.0 | https://github.com/crossbeam-rs/crossbeam |
| crossbeam-deque | 0.8.8 | MIT OR Apache-2.0 | https://github.com/crossbeam-rs/crossbeam |
| crossbeam-epoch | 0.9.21 | MIT OR Apache-2.0 | https://github.com/crossbeam-rs/crossbeam |
| crossbeam-utils | 0.8.23 | MIT OR Apache-2.0 | https://github.com/crossbeam-rs/crossbeam |
| either | 1.19.0 | MIT OR Apache-2.0 | https://github.com/rayon-rs/either |
| equivalent | 1.0.3 | Apache-2.0 OR MIT | https://github.com/indexmap-rs/equivalent |
| gb-io | 0.9.0 | MIT | https://github.com/dlesl/gb-io |
| hashbrown | 0.17.1 | MIT OR Apache-2.0 | https://github.com/rust-lang/hashbrown |
| heck | 0.5.0 | MIT OR Apache-2.0 | https://github.com/withoutboats/heck |
| indexmap | 2.14.2 | Apache-2.0 OR MIT | https://github.com/indexmap-rs/indexmap |
| indoc | 2.0.8 | MIT OR Apache-2.0 | https://github.com/dtolnay/indoc |
| itertools | 0.14.0 | MIT OR Apache-2.0 | https://github.com/rust-itertools/itertools |
| lexical-core | 1.0.6 | MIT / Apache-2.0 | https://github.com/Alexhuszagh/rust-lexical |
| lexical-parse-float | 1.0.6 | MIT / Apache-2.0 | https://github.com/Alexhuszagh/rust-lexical |
| lexical-parse-integer | 1.0.6 | MIT / Apache-2.0 | https://github.com/Alexhuszagh/rust-lexical |
| lexical-util | 1.0.7 | MIT / Apache-2.0 | https://github.com/Alexhuszagh/rust-lexical |
| lexical-write-float | 1.0.6 | MIT / Apache-2.0 | https://github.com/Alexhuszagh/rust-lexical |
| lexical-write-integer | 1.0.6 | MIT / Apache-2.0 | https://github.com/Alexhuszagh/rust-lexical |
| libc | 0.2.191 | MIT OR Apache-2.0 | https://github.com/rust-lang/libc |
| log | 0.4.34 | MIT OR Apache-2.0 | https://github.com/rust-lang/log |
| memchr | 2.8.3 | Unlicense OR MIT | https://github.com/BurntSushi/memchr |
| memoffset | 0.9.1 | MIT | https://github.com/Gilnaa/memoffset |
| nom | 8.0.0 | MIT | https://github.com/rust-bakery/nom |
| noodles | 0.118.0 | MIT | https://github.com/zaeleus/noodles |
| noodles-bam | 0.97.0 | MIT | https://github.com/zaeleus/noodles |
| noodles-bgzf | 0.52.0 | MIT | https://github.com/zaeleus/noodles |
| noodles-core | 0.21.0 | MIT | https://github.com/zaeleus/noodles |
| noodles-csi | 0.62.0 | MIT | https://github.com/zaeleus/noodles |
| noodles-sam | 0.92.0 | MIT | https://github.com/zaeleus/noodles |
| noodles-tabix | 0.68.0 | MIT | https://github.com/zaeleus/noodles |
| noodles-vcf | 0.95.0 | MIT | https://github.com/zaeleus/noodles |
| once_cell | 1.21.4 | MIT OR Apache-2.0 | https://github.com/matklad/once_cell |
| percent-encoding | 2.3.2 | MIT OR Apache-2.0 | https://github.com/servo/rust-url/ |
| portable-atomic | 1.15.0 | Apache-2.0 OR MIT | https://github.com/taiki-e/portable-atomic |
| proc-macro-crate | 3.5.0 | MIT OR Apache-2.0 | https://github.com/bkchr/proc-macro-crate |
| proc-macro2 | 1.0.107 | MIT OR Apache-2.0 | https://github.com/dtolnay/proc-macro2 |
| pyo3 | 0.23.5 | MIT OR Apache-2.0 | https://github.com/pyo3/pyo3 |
| pyo3-build-config | 0.23.5 | MIT OR Apache-2.0 | https://github.com/pyo3/pyo3 |
| pyo3-ffi | 0.23.5 | MIT OR Apache-2.0 | https://github.com/pyo3/pyo3 |
| pyo3-macros | 0.23.5 | MIT OR Apache-2.0 | https://github.com/pyo3/pyo3 |
| pyo3-macros-backend | 0.23.5 | MIT OR Apache-2.0 | https://github.com/pyo3/pyo3 |
| quick-xml | 0.42.0 | MIT | https://github.com/tafia/quick-xml |
| quote | 1.0.47 | MIT OR Apache-2.0 | https://github.com/dtolnay/quote |
| rayon | 1.12.0 | MIT OR Apache-2.0 | https://github.com/rayon-rs/rayon |
| rayon-core | 1.13.0 | MIT OR Apache-2.0 | https://github.com/rayon-rs/rayon |
| rustversion | 1.0.23 | MIT OR Apache-2.0 | https://github.com/dtolnay/rustversion |
| serde | 1.0.229 | MIT OR Apache-2.0 | https://github.com/serde-rs/serde |
| serde_core | 1.0.229 | MIT OR Apache-2.0 | https://github.com/serde-rs/serde |
| serde_derive | 1.0.229 | MIT OR Apache-2.0 | https://github.com/serde-rs/serde |
| syn | 2.0.119 | MIT OR Apache-2.0 | https://github.com/dtolnay/syn |
| syn | 3.0.9 | MIT OR Apache-2.0 | https://github.com/dtolnay/syn |
| target-lexicon | 0.12.16 | Apache-2.0 WITH LLVM-exception | https://github.com/bytecodealliance/target-lexicon |
| thiserror | 2.0.21 | MIT OR Apache-2.0 | https://github.com/dtolnay/thiserror |
| thiserror-impl | 2.0.21 | MIT OR Apache-2.0 | https://github.com/dtolnay/thiserror |
| toml_datetime | 1.1.2 | MIT OR Apache-2.0 | https://github.com/toml-rs/toml |
| toml_edit | 0.25.17 | MIT OR Apache-2.0 | https://github.com/toml-rs/toml |
| toml_parser | 1.1.5 | MIT OR Apache-2.0 | https://github.com/toml-rs/toml |
| unicode-ident | 1.0.27 | (MIT OR Apache-2.0) AND Unicode-3.0 | https://github.com/dtolnay/unicode-ident |
| unindent | 0.2.4 | MIT OR Apache-2.0 | https://github.com/dtolnay/indoc |
| winnow | 1.0.4 | MIT | https://github.com/winnow-rs/winnow |
| zlib-rs | 0.6.8 | Zlib | https://github.com/trifectatechfoundation/zlib-rs |

## Data

- **`app/parts.tsv`** (the "Detect Common Features" library, 94 sequences): built by `scripts/build_parts.py` from annotated features of public
  GenBank / INSDC records of well-known vectors, such as pBR322 (J01749), pUC19 (L09137), pKD46 (AY048746), pEGFP-N1 (U55762), pCas9 (KM099231) and
  pCAMBIA-1300 (AF234296); every row names its source accession in the `source` column. The records are public database entries
  (policies: https://www.ncbi.nlm.nih.gov/home/about/policies/). NCBI places no restrictions of its own on their use, but the people who submitted a
  record may claim rights in it, so check the record itself if you want to reuse a sequence outside this program. Short functional elements
  (promoters, origins, resistance genes) are, as far as we can tell, facts about biology; we have not reviewed each record's submitter terms or patents one by one.
- **Restriction enzymes** (`app/enzymes.py`, 47 common enzymes): recognition sites and cut positions, which are facts about the enzymes. Not taken from
  REBASE or any supplier's catalogue.
- **Gel ladders** (`app/tools.py`): "typical" band sets of common ladder types, written for OpenGenEdit; not any manufacturer's exact list.
- **Golden Gate overhangs** (`app/overhangs.py`): the 25 high-fidelity overhangs and the junction sets of the 12- and 24-fragment cassettes from
  Potapov et al., ACS Synth Biol 7, 2665-2674 (2018), doi:10.1021/acssynbio.8b00333 (supplementary tables). Facts published in a paper, cited in Help > Citations.
- **Demo files** (`demo/demo_*.gb`): random sequences generated by `scripts/make_demo_files.py`. They are not real parts.
- **Example sequences** (`demo/pBAD30_eilR.gb`, `demo/pFG051_Chlor.gb`, `demo/pFG051_Chlor.dna`, `demo/BW25113_2026.gb`): real plasmids and an
  *E. coli* BW25113 genome built by the author, who owns them and releases them with the program under the same licence. They are shipped so that
  the tools can be tried on real constructs rather than on random DNA; the `.dna` file is there so the .dna reader can be exercised on a real
  file. No third party holds rights in them.
- **Data fetched while the program runs** (only when you use the feature; see [PRIVACY.md](PRIVACY.md)): NCBI records (public; per-record terms),
  AlphaFold DB models (CC BY 4.0, EMBL-EBI and Google DeepMind; please cite Jumper et al. 2021 and Varadi et al. 2022), RCSB PDB files (CC0), ESM Atlas
  predictions (Meta; see its terms). Help > Citations lists the references for each.

## Names

SnapGene, BLAST, AlphaFold, ESMFold, Addgene, NEB, Qt, GitHub, Windows and other names in this repository are trademarks of their owners and are used only
to say what a file format, a method or a service is. OpenGenEdit is an independent project: it is not affiliated with, endorsed by or sponsored by any of them.

The project's own names, "OpenGenEdit" and "oDNA", and its wordmark and icon: how they may be used, and why a modified version
takes another name, is in [TRADEMARKS.md](TRADEMARKS.md).
