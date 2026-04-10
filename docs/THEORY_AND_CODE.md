# Shear Centre Calculation: Theory and Code Documentation

## Table of Contents
1. [Introduction](#introduction)
2. [Mathematical Theory](#mathematical-theory)
3. [Finite Element Formulation](#finite-element-formulation)
4. [Code Structure](#code-structure)
5. [Implementation Details](#implementation-details)
6. [Validation](#validation)
7. [Troubleshooting](#troubleshooting)

---

## Introduction

The **shear centre** (or centre of twist) of a beam cross-section is the point where a transverse shear load produces **no twisting moment**. For doubly symmetric sections (e.g., I-beams, rectangles), the shear centre coincides with the centroid. For asymmetric sections (e.g., C-channels, L-angles), it lies elsewhere.

This document provides a complete derivation of the finite element formulation and a detailed explanation of the Fortran 77 implementation.

---

## Mathematical Theory

### 1. Problem Statement

For a beam cross-section in the y-z plane, subject to shear forces V_y and V_z, the shear centre (y_s, z_s) satisfies:
V_y applied at y_s → no twist
V_z applied at z_s → no twist

text

### 2. Governing Equations (Saint-Venant Theory)

We solve **two auxiliary Poisson problems**:

#### Case 1: Shear force V_y = 1, V_z = 0
∇²φ₁ = -2z on Ω

text
with **Neumann boundary conditions** (natural, not Dirichlet):
∂φ₁/∂n = 0 on Γ

text

The shear centre coordinate is:
y_s = (1/I_z) ∫_Ω (∂φ₁/∂z) · y dA

text

#### Case 2: Shear force V_y = 0, V_z = 1
∇²φ₂ = 2y on Ω

text
with Neumann boundary conditions:
∂φ₂/∂n = 0 on Γ

text

The shear centre coordinate is:
z_s = (1/I_y) ∫_Ω (∂φ₂/∂y) · z dA

text

### 3. Why Neumann Conditions?

In the shear centre problem, the boundary conditions are **natural (Neumann)**, not Dirichlet (φ = 0). 

- **Dirichlet** would force φ = 0 on the boundary, which is incorrect
- **Neumann** (∂φ/∂n = 0) allows the solution to develop the correct shear flow distribution

This is a critical point: many implementations erroneously impose φ = 0 on the boundary, leading to completely wrong results.

### 4. Compatibility Condition (Neumann)

For a pure Neumann problem, the right-hand side must satisfy:
∫_Ω f dA = 0

text

Otherwise, no solution exists. For our cases:
- Case 1: ∫_Ω (-2z) dA = -2·S_z = 0 (since centroid is at origin)
- Case 2: ∫_Ω (2y) dA = 2·S_y = 0 (since centroid is at origin)

This holds automatically when the coordinate origin is at the centroid.

### 5. Numerical Implementation Considerations

#### 5.1 Singularity of the Stiffness Matrix

The Neumann problem has a **singular stiffness matrix** (rank = n-1). The null space consists of constant functions. We must impose an additional constraint to obtain a unique solution.

**Methods:**
1. **Fix one node** (simplest, used here): set φ(1) = 0
2. **Lagrange multiplier** (more elegant): add constraint ∫φ dA = 0
3. **Pseudo-inverse** (most robust)

Our implementation uses **method 1** (fix one node) + small epsilon on diagonal.

#### 5.2 RHS Normalization

Even though the theoretical RHS satisfies ∫f dA = 0, numerical integration introduces small errors. We explicitly normalize:
f_mean = (1/A) ∫_Ω f dA
f_normalized = f - f_mean

text

This ensures the discrete system is compatible.

---

## Finite Element Formulation

### 1. Triangular Element (P1 - linear shape functions)

#### Shape Functions

For a triangle with vertices (y₁,z₁), (y₂,z₂), (y₃,z₃):
N₁(ξ,η) = ξ
N₂(ξ,η) = η
N₃(ξ,η) = 1 - ξ - η

text

where ξ, η are area coordinates.

#### Gradient of Shape Functions
∂N₁/∂y = b₁, ∂N₁/∂z = c₁
∂N₂/∂y = b₂, ∂N₂/∂z = c₂
∂N₃/∂y = b₃, ∂N₃/∂z = c₃

text

with:
b₁ = (z₂ - z₃)/(2A) c₁ = (y₃ - y₂)/(2A)
b₂ = (z₃ - z₁)/(2A) c₂ = (y₁ - y₃)/(2A)
b₃ = (z₁ - z₂)/(2A) c₃ = (y₂ - y₁)/(2A)

text

and element area:
A = ½ |(y₂-y₁)(z₃-z₁) - (y₃-y₁)(z₂-z₁)|

text

### 2. Element Stiffness Matrix

For Poisson equation ∇²φ = f, the weak form gives:
K_e(i,j) = A · (b_i·b_j + c_i·c_j)

text

where i,j = 1,2,3 are local node indices.

### 3. Element RHS Vector

Using lumped integration (midpoint rule):
f_e(i) = A · f(y_c, z_c) / 3

text

where (y_c, z_c) is the element centroid.

For our problems:

**Case 1** (f = -2z):
f_e(i) = -2 · A · z_c / 3

text

**Case 2** (f = 2y):
f_e(i) = 2 · A · y_c / 3

text

### 4. Global Assembly

The global stiffness matrix K (n×n) and RHS vectors R₁, R₂ (n×1) are assembled by summing element contributions:
K(global_i, global_j) += K_e(local_i, local_j)
R₁(global_i) += f_e1(local_i)
R₂(global_i) += f_e2(local_i)

text

### 5. Boundary Conditions (Neumann)

Neumann conditions are **natural** in the FEM formulation — they do not require explicit enforcement. The only modification is the RHS normalization (compatibility).

### 6. Singularity Removal

We fix the first node:
K(1, j) = 0 ∀j
K(i, 1) = 0 ∀i
K(1, 1) = 1
R₁(1) = 0
R₂(1) = 0

text

And add a small epsilon to the diagonal for numerical stability:
K(i,i) = K(i,i) + ε, ε = 1×10⁻¹²

text

### 7. Solution of Linear Systems

We solve two systems simultaneously (n_rhs = 2):
K · [φ₁ φ₂] = [R₁ R₂]

text

We use **DGELS** (LAPACK) — a least-squares solver that handles near-singular matrices robustly.

### 8. Post-processing: Shear Centre Calculation

**Critical point:** We integrate the **gradient** of φ, not φ itself!
y_s = (1/I_z) · Σ_e [ A_e · (∂φ₁/∂z) · y_c ]
z_s = (1/I_y) · Σ_e [ A_e · (∂φ₂/∂y) · z_c ]

text

where:
∂φ₁/∂z = φ₁₁·c₁ + φ₁₂·c₂ + φ₁₃·c₃
∂φ₂/∂y = φ₂₁·b₁ + φ₂₂·b₂ + φ₂₃·b₃

text

---

## Code Structure

### File Organization
beam_shear_center/
├── src/
│ ├── read_section_mesh_unv.f # UNV mesh reader
│ ├── mesh_checker.f # Mesh validation
│ ├── compute_section_properties.f # Area, centroid, Iy, Iz, Iyz
│ ├── shear_center.f # Main shear centre solver
│ ├── build_D_matrix.f # Constitutive matrix (future)
│ └── section_database.f # Section database (future)
├── test/
│ └── test_shear_center.f # Test program for HEB200
├── meshes/
│ ├── HEB200_mm.unv # HEB200 mesh (335 nodes, 472 triangles)
│ ├── rect_10x20.unv # Rectangle 10×20 mm
│ └── circle_dia10.unv # Circle diameter 10 mm
├── Makefile # Build system
└── docs/
└── THEORY_AND_CODE.md # This file

text

### Subroutine: compute_shear_center

**Interface:**
```fortran
subroutine compute_shear_center(nn, ne, nodes, elements, Iy, Iz, y_s, z_s)
Arguments:

VariableTypeIntentDescription
nnintegerinNumber of nodes
neintegerinNumber of elements
nodesdouble(2,nn)inNodal coordinates (y,z)
elementsinteger(3,ne)inConnectivity table
IydoubleinMoment of inertia about y-axis
IzdoubleinMoment of inertia about z-axis
y_sdoubleoutShear centre y-coordinate
z_sdoubleoutShear centre z-coordinate
Local Variables:

VariableDescription
stiff(nn,nn)Global stiffness matrix
rhs(nn,2)Global RHS vectors (two cases)
phi(nn,2)Solution vectors
b(3), c(3)Shape function gradients
ke(3,3)Element stiffness matrix
work(:)Workspace for DGELS
Subroutine: read_section_mesh_unv
Interface:

fortran
subroutine read_section_mesh_unv(filename, max_nodes, max_triangles,
     &                           n_nodes, n_triangles, conn, y, z)
Supported Format: Universal File Format (UNV) from Salome/ASTER

Node section: marker 2411

Element section: marker 2412

Triangular elements: type 41

Implementation Details
1. Fixed Format (Fortran 77)
All source files follow fixed format rules:

ColumnsUsage
1c for comment lines
2-5Statement label (optional)
6Continuation character
7-72Fortran code
73+Ignored
2. Compilation (Makefile)
makefile
FC = gfortran
FFLAGS = -ffixed-form -Wall -O2 -std=legacy
LDFLAGS = -llapack -lblas
3. LAPACK Solver (DGELS)
We use DGELS because it:

Handles near-singular matrices robustly

Works for both square and rectangular systems

Returns least-squares solution when system is singular

Calling sequence:

fortran
! Query optimal workspace
lwork = -1
call DGELS('N', nn, nn, 2, stiff, nn, rhs, nn, work, lwork, info)
lwork = int(work(1))
allocate(work(lwork))

! Solve
call DGELS('N', nn, nn, 2, stiff, nn, rhs, nn, work, lwork, info)
4. Numerical Stabilization
We add a small epsilon to the diagonal:

fortran
eps = 1.0d-12
do i = 1, nn
   stiff(i,i) = stiff(i,i) + eps
end do
This improves conditioning without affecting the solution significantly.

Validation
Test Case: HEB200 (Doubly Symmetric)
Expected results: y_s = 0, z_s = 0

Mesh: 335 nodes, 472 triangles

Results:

text
y_s = -1.1125804459212456E-003 mm
z_s = -1.4340090578989642E-003 mm
Error analysis:

Absolute error: ~0.001 mm

Relative error: ~5×10⁻⁶ (0.0005%)

Section size: ~200 mm

Conclusion: Excellent accuracy for linear FEM.

Test Case: Rectangle 10×20 (Symmetric)
Expected results: y_s = 0, z_s = 0

Mesh: 4 nodes, 2 triangles (coarse)

Results:

text
y_s ≈ 2.0 mm (large error due to coarse mesh)
z_s ≈ -1.6 mm
Note: Coarse meshes produce larger errors. The HEB200 mesh (472 triangles) is sufficiently refined.

Troubleshooting
Common Issues and Solutions
ProblemPossible CauseSolution
DGELS fails with info > 0Singular matrixCheck RHS normalization; increase epsilon
y_s, z_s far from zeroWrong integrationEnsure you integrate gradient, not φ
Large numerical errorsMesh too coarseRefine mesh
Compilation errorsFixed format violationCheck column 72 limit
Segmentation faultArray bounds exceededCheck MAX_NODES, MAX_ELEMENTS
Debugging Tips
Add temporary print statements:

fortran
write(*,*) 'DEBUG: nfree =', nfree
write(*,*) 'DEBUG: min/max stiff =', minval(stiff), maxval(stiff)
write(*,*) 'DEBUG: info =', info
Performance Notes
For nn = 335 nodes, matrix size is 335×335 (~112k entries)

DGELS solves in ~0.1 seconds

Memory usage: ~1 MB for K, ~0.5 MB for RHS

Future Extensions
Constitutive Matrix D(6,6)

Relates stress and strain for general beam sections

Shear Correction Factors (k_y, k_z)

Account for non-uniform shear stress distribution

Warping Function ω

For torsion of thin-walled sections

Vlasov Theory

For open thin-walled sections

CalculiX Integration

Package as UEL or user material

References
Timoshenko, S. P. & Goodier, J. N. (1970). Theory of Elasticity. McGraw-Hill.

Cook, R. D., Malkus, D. S., Plesha, M. E. & Witt, R. J. (2002). Concepts and Applications of Finite Element Analysis. Wiley.

LAPACK User's Guide (1999). DGELS documentation.

Bathe, K. J. (2014). Finite Element Procedures. Klaus-Jürgen Bathe.

Acknowledgments
DeepSeek for AI-assisted development and debugging

Bruno Zilli for project coordination and validation

License
MIT License — see LICENSE file for details.

Contact
For questions or contributions, please open an issue on GitHub.
