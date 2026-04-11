# Theory and Implementation Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Geometric Properties](#geometric-properties)
3. [Shear Centre Theory](#shear-centre-theory)
4. [Finite Element Formulation](#finite-element-formulation)
5. [Code Structure](#code-structure)
6. [Validation](#validation)
7. [References](#references)

---

## Introduction

This software computes the **shear centre** (also called centre of twist or
flexural centre) for arbitrary thin-walled beam sections using the Finite
Element Method (FEM). The shear centre is the point in the cross-section
where a transverse load produces no twisting moment.

### Why is the shear centre important?

For asymmetric sections (L-channels, C-channels, etc.), the shear centre
does NOT coincide with the centroid. Loading applied at the shear centre
produces bending without torsion. Loading applied elsewhere produces
coupled bending-torsion response.

---

## Geometric Properties

### Section Properties

For a planar cross-section in the `(y,z)` plane, the following properties
are computed:

| Property | Symbol | Formula |
|----------|--------|---------|
| Area | `A` | `∫ dA` |
| First moments | `S_y, S_z` | `∫ y dA, ∫ z dA` |
| Centroid | `(ȳ, z̄)` | `(S_y/A, S_z/A)` |
| Moments of inertia | `I_y, I_z` | `∫ (z-z̄)² dA, ∫ (y-ȳ)² dA` |
| Product of inertia | `I_yz` | `∫ (y-ȳ)(z-z̄) dA` |

### Finite Element Computation

For triangular elements:
Area = ½ |(y₂-y₁)(z₃-z₁) - (y₃-y₁)(z₂-z₁)|

Centroid = (y₁+y₂+y₃, z₁+z₂+z₃) / 3

I_y = Σ area * [ (z̄ - z_c)² + (z̄ - z_c)²? ] [full formula in code]

text

---

## Shear Centre Theory

### Fundamental Equations

The shear centre coordinates `(y_s, z_s)` are given by:
y_s = (I_z · M_y - I_yz · M_z) / (I_y I_z - I_yz²)

z_s = (I_y · M_z - I_yz · M_y) / (I_y I_z - I_yz²)

text

where `M_y` and `M_z` are integrals of the shear flow:
M_y = ∫ (∂φ₁/∂z) · y dA
M_z = ∫ (∂φ₂/∂y) · z dA

text

### Prandtl Stress Function Analogy

The functions `φ₁` and `φ₂` satisfy Poisson's equation:
∇²φ₁ = -2·(z - z̄)
∇²φ₂ = +2·(y - ȳ)

text

with Neumann boundary conditions on free edges:
∂φ/∂n = 0

text

### Variational Formulation

The weak form (Galerkin) is:
∫ (∇φ · ∇v) dA = ∫ f v dA ∀v ∈ H¹(Ω)

text

This leads to the linear system:
K · φ = f

text

where `K` is the stiffness matrix (Laplacian), `f` is the forcing vector.

---

## Finite Element Formulation

### Element Type

**Linear triangular element (T3)** with 3 nodes and linear shape functions:
N₁ = 1 - ξ - η
N₂ = ξ
N₃ = η

text

### Shape Function Gradients

For an element with nodes `(y₁,z₁)`, `(y₂,z₂)`, `(y₃,z₃)`:
Area = ½ |(y₂-y₁)(z₃-z₁) - (y₃-y₁)(z₂-z₁)|

b₁ = (z₂ - z₃) / (2·Area)
b₂ = (z₃ - z₁) / (2·Area)
b₃ = (z₁ - z₂) / (2·Area)

c₁ = (y₃ - y₂) / (2·Area)
c₂ = (y₁ - y₃) / (2·Area)
c₃ = (y₂ - y₁) / (2·Area)

text

### Element Stiffness Matrix
Kᵉ_{ij} = Area · (b_i·b_j + c_i·c_j)

text

### Element RHS Vector

For shear in z-direction (case 1):
fᵉ_i = -Area · (z̄_element - z̄_global) / 3

text

For shear in y-direction (case 2):
fᵉ_i = +Area · (ȳ_element - ȳ_global) / 3

text

### Global Assembly
K = Σ Kᵉ
f = Σ fᵉ

text

---

## Code Structure

### Main Subroutine: `compute_shear_center`

```fortran
subroutine compute_shear_center(nn, ne, nodes, elements,
     +                          Iy, Iz, Iyz, y_s, z_s)
Input:

nn : number of nodes

ne : number of elements

nodes(2, nn) : nodal coordinates (y, z)

elements(3, ne) : element connectivity

Iy, Iz, Iyz : moments of inertia (from compute_section_properties)

Output:

y_s, z_s : shear centre coordinates (relative to centroid)

Step-by-Step Implementation
1. Centroid Computation
fortran
y_bar = 0.0d0; z_bar = 0.0d0; total_area = 0.0d0
do e = 1, ne
    ! compute element area and centroid
    ! accumulate total_area and weighted centroids
end do
y_bar = y_bar / total_area
z_bar = z_bar / total_area
2. System Assembly
fortran
do e = 1, ne
    ! get element geometry
    ! compute area, b(), c()
    
    ! element stiffness
    do i = 1,3; do j = 1,3
        ke(i,j) = area*(b(i)*b(j)+c(i)*c(j))
    end do; end do
    
    ! assemble into global K
    do i = 1,3; do j = 1,3
        K(elem(i,e), elem(j,e)) += ke(i,j)
    end do; end do
    
    ! RHS for both load cases
    do i = 1,3
        rhs(elem(i,e),1) -= area * (zc - z_bar) / 3
        rhs(elem(i,e),2) += area * (yc - y_bar) / 3
    end do
end do
3. Neumann Compatibility
For Poisson problems with pure Neumann BCs, the RHS must satisfy:

text
∫ f dA = 0
fortran
fmean = sum(rhs) / nn
rhs = rhs - fmean
4. Singularity Removal
The stiffness matrix is singular (nullspace of constants). We fix by:

fortran
! Pin node 1 to zero
do j = 1, nn
    K(1,j) = 0.0d0; K(j,1) = 0.0d0
end do
K(1,1) = 1.0d0
rhs(1,:) = 0.0d0

! Add small diagonal perturbation
eps = 1.0d-12
do i = 1, nn
    K(i,i) = K(i,i) + eps
end do
5. Linear Solver (DGELS)
We use LAPACK's DGELS which solves:

text
min ||K·φ - rhs||₂
fortran
call DGELS('N', nn, nn, 2, K, nn, phi, nn, work, lwork, info)
6. Shear Flow Integration
fortran
do e = 1, ne
    ! recompute b(), c() for element
    
    ! shear flow from gradients
    qy = phi(n1,1)*b(1) + phi(n2,1)*b(2) + phi(n3,1)*b(3)
    qz = phi(n1,2)*c(1) + phi(n2,2)*c(2) + phi(n3,2)*c(3)
    
    ! moment integrals
    M_y = M_y + area * qz * yc
    M_z = M_z + area * qy * zc
end do
7. Tensor Inversion
fortran
detI = Iy * Iz - Iyz * Iyz

y_s = (Iz * M_y - Iyz * M_z) / detI
z_s = (Iy * M_z - Iyz * M_y) / detI
Validation
Test Case 1: HEB200 (Doubly Symmetric)
Property	Theoretical	FEM	Error
Area [mm²]	11340	11340	0%
y_s [mm]	0	<1e-6	✓
z_s [mm]	0	<1e-6	✓
Expected: Shear centre at centroid.

Test Case 2: L100×100×10 (Angle Section)
Property	Theoretical	FEM	Error
Area [mm²]	1900	1915	0.8%
Distance from corner [mm]	35.7	34.6	3.1%
Expected: Shear centre on line y=z, at ~35.7 mm from corner.

References
Timoshenko, S.P. (1970). Theory of Elasticity. McGraw-Hill.

Cook, R.D. (2002). Concepts and Applications of Finite Element
Analysis. Wiley.

Vlasov, V.Z. (1961). Thin-Walled Elastic Beams. Israel Program for
Scientific Translations.

LAPACK User's Guide (1999). SIAM.

Salome Platform - Open-source pre/post-processor.
https://www.salome-platform.org

Appendix: Notation
Symbol	Description	Units
A	Cross-sectional area	mm²
I_y, I_z	Moments of inertia	mm⁴
I_yz	Product of inertia	mm⁴
y_s, z_s	Shear centre coordinates	mm
φ₁, φ₂	Stress functions	mm³
K	Stiffness matrix	-
f	RHS vector	mm
∇²	Laplacian operator	-
Document Version: 1.0
*Last Updated: 2024-04-11*
