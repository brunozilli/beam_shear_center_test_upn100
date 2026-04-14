# Open Section Shear Centre Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Fortran](https://img.shields.io/badge/Fortran-F77%2FF90-blue.svg)](https://fortran-lang.org)

A robust, industrial-grade solver for computing the **shear centre** of 
open thin-walled sections (UPN, C-channels, L-sections) from a 1D midline 
mesh or 2D triangular mesh.

## Features

- ✅ Supports both **1D edge meshes** (UNV format) and **2D triangular meshes**
- ✅ Automatic **web/flange classification** based on local tangent
- ✅ **Open chain reconstruction** from leftmost to rightmost flange tip
- ✅ **Correct shear flow integration** with projected force components
- ✅ **Area-weighted centroid** calculation for 2D meshes
- ✅ **Explicit X/Y axes** output (flange axis vs. web axis)
- ✅ **Pure Fortran 77** (fixed-form) compatible with `gfortran`
- ✅ **No external dependencies** beyond LAPACK/BLAS

## Quick Start

### Prerequisites

- `gfortran` (or any Fortran compiler)
- LAPACK and BLAS libraries

### Compilation

```bash
make clean && make
Running a Test
bash
./test/test_shear_center meshes/upn100_1D.unv
Expected Output (UPN100)
text
Shear Center Results:
  X_sc = -26.94 mm  (flange axis - out of web)
  Y_sc =   1.00 mm  (web axis - vertical)
Mesh Requirements
1D Mesh (UNV Format)
Elements of type 11 (linear edges) or 41 (triangles)

Nodes should lie on the midline of the section

Open chain topology (two free ends)

2D Mesh (UNV Format)
Triangular elements (type 41)

The solver automatically extracts the boundary contour

File Structure
text
.
├── src/
│   ├── shear_center.f              # Main solver
│   ├── compute_section_properties.f # Section properties (2D)
│   ├── read_section_mesh_unv.f     # UNV mesh reader
│   ├── mesh_checker.f              # Mesh validation
│   ├── build_D_matrix.f            # Utility routines
│   └── section_database.f          # Section database
├── test/
│   └── test_shear_center.f         # Test driver
├── meshes/
│   └── upn100_1D.unv               # Example UPN100 midline mesh
├── Makefile
├── LICENSE.md
├── DISCLAIMER.md
├── THEORY.md
└── README.md
API Reference
Main Subroutine
fortran
subroutine compute_shear_center(nn, ne, nodes, elements, 
                                Ixx, Iyy, Ixy, x_sc, y_sc)
Parameter	Type	Description
nn	integer	Number of nodes
ne	integer	Number of elements
nodes	double precision(3,nn)	Node coordinates (X, Y, Z)
elements	integer(3,ne)	Element connectivity
Ixx, Iyy, Ixy	double precision	Second moments of area
x_sc, y_sc	double precision	Shear centre coordinates
Limitations
Open sections only – closed sections require Bredt's constant

Constant thickness assumed – variable thickness not modelled

Coarse meshes may introduce integration errors

Contributing
Contributions are welcome! Please open an issue or submit a pull request.

Authors
Bruno Zilli

DeepSeek

License
This project is licensed under the MIT License – see the LICENSE.md file for details.

Disclaimer
This software is provided for educational and research purposes only.
See DISCLAIMER.md for full details.

References
Vlasov, V.Z. (1961). Thin-Walled Elastic Beams

Timoshenko, S.P. & Goodier, J.N. (1970). Theory of Elasticity

Pilkey, W.D. (2002). Analysis and Design of Elastic Beams
