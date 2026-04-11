# Disclaimer

## General Disclaimer

This software is provided "AS IS", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement.

## Engineering Accuracy

This software performs finite element calculations for beam section
properties. While every effort has been made to ensure correctness:

1. **No Guarantee of Accuracy**: The results may contain errors due to mesh
   discretisation, numerical integration, or implementation bugs.

2. **Validation Required**: Always validate results against:
   - Analytical solutions (where available)
   - Commercial FEM software (Abaqus, ANSYS, etc.)
   - Experimental data for critical applications

3. **Mesh Sensitivity**: Results depend on mesh quality and refinement.
   Convergence studies are recommended for each new section.

## Safety Critical Applications

**DO NOT USE** this software as the sole basis for:
- Structural design decisions
- Safety-critical components
- Life-support systems
- Nuclear facilities
- Aerospace applications
- Medical devices
- Any application where failure could cause injury, death, or property damage

Always consult a qualified structural engineer for final design validation.

## User Responsibility

The user assumes all responsibility for:
- Correct interpretation of results
- Verification of input mesh quality
- Validation of outputs for their specific application
- Compliance with relevant building codes and standards
- Independent verification of critical calculations

## Known Limitations

- 2D planar sections only (no 3D effects)
- Linear elastic material behaviour only (no plasticity)
- Small deformations assumed (geometrically linear)
- Mesh must be triangular and well-formed
- No shear lag or local buckling effects
- No warping restraint effects
- No thermal loading
- No dynamic effects

## Third Party Libraries

This software uses LAPACK and BLAS libraries. See their respective licences
for terms and conditions. The authors of this software are not responsible
for any issues arising from these third-party libraries.

## No Professional Advice

This software is a **computational tool**, not a substitute for professional
engineering judgment. All results should be reviewed by a qualified
structural engineer.

## Updates and Corrections

This software is under active development. No guarantee is made that future
versions will maintain API compatibility or produce identical results.

## Export Control

This software is purely mathematical and contains no encryption or
restricted technologies. However, users are responsible for compliance with
local export control laws.

## Governing Law

Any dispute arising from the use of this software shall be governed by the
laws of Italy, without regard to its conflict of law provisions.

---

**YOU HAVE BEEN WARNED.**

*Last Updated: 2024-04-11*
