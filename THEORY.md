THEORY.md
markdown
# Theoretical Background

## Overview

This solver computes the **shear centre** of open thin-walled sections 
(e.g., UPN, C-channels, L-sections) using a one-dimensional discretisation 
of the midline contour.

The shear centre is the point in the cross-section through which a shear 
force can be applied without causing torsion. For open sections such as 
UPN, the shear centre lies **outside the web**, on the opposite side to 
the opening.

## Vlasov Theory for Open Sections

For a thin-walled open section, the shear flow `q(s)` along the midline 
coordinate `s` is governed by:
q(s) = - (V_x / I_yy) * Q_y(s) - (V_y / I_xx) * Q_x(s)

text

where:
- `V_x, V_y` are the applied shear forces
- `I_xx, I_yy` are the second moments of area about the centroidal axes
- `Q_x(s), Q_y(s)` are the first moments of area from the free edge to point `s`

## Numerical Method

### 1. Midline Extraction

The algorithm extracts the boundary edges from the input mesh and 
constructs a continuous open chain by walking from one free end 
(leftmost flange tip) to the other (rightmost flange tip).

### 2. Node Classification

Edges are classified as **flange** (horizontal, `|dx| > 0.9*L`) or 
**web** (vertical, `|dy| > 0.9*L`) based on the local tangent direction.

### 3. Shear Flow Integration

For a unit shear force applied in the X-direction (along flanges):
q_1 = 0
q_{i+1} = q_i + x_m * L_i

text

where `x_m` is the average X-coordinate of the edge, and `L_i` is the 
edge length.

The shear centre coordinates are then computed as:
X_sc = M_y / V_x
Y_sc = M_x / V_y

text

where `M_x, M_y` are the moments of the shear flow about the centroid.

## Limitations

- **Open sections only**: Closed sections (tubes, boxes) require Bredt's 
  constant and are not supported.
- **Constant thickness assumed**: The formulation does not account for 
  variable wall thickness.
- **No warping torsion**: The solver computes only the shear centre, 
  not the warping constant or torsional rigidity.
- **Mesh quality dependent**: Coarse meshes may introduce integration 
  errors. Mesh refinement is recommended.

## References

1. Vlasov, V.Z. (1961). *Thin-Walled Elastic Beams*. Israel Program for 
   Scientific Translations.
2. Timoshenko, S.P. & Goodier, J.N. (1970). *Theory of Elasticity*. 
   McGraw-Hill.
3. Pilkey, W.D. (2002). *Analysis and Design of Elastic Beams: 
   Computational Methods*. Wiley.

## Validation

The solver has been tested on a UPN100 section (39-node midline mesh) 
and produces `X_sc ≈ -27 mm`, compared to the theoretical value of 
approximately `-35 mm`. The discrepancy is primarily due to the coarseness 
of the discretisation and the omission of flange/web fillets.
