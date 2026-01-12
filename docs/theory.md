# Cartesian Decomposition of Vibrational Intensities

This document summarizes the analytical background underlying the XYZ decomposition of IR and Raman intensities implemented in VibXYZ.

---

## IR intensity

For a vibrational normal mode \( Q \), the infrared intensity is proportional to the square of the dipole derivative:

\[
I_{\mathrm{IR}} \propto \left|\frac{\partial \boldsymbol{\mu}}{\partial Q}\right|^2
\]

Expanding into Cartesian components:

\[
\left|\frac{\partial \boldsymbol{\mu}}{\partial Q}\right|^2
=
\left(\frac{\partial \mu_x}{\partial Q}\right)^2
+
\left(\frac{\partial \mu_y}{\partial Q}\right)^2
+
\left(\frac{\partial \mu_z}{\partial Q}\right)^2
\]

Gaussian reports only the total activity. VibXYZ extracts each Cartesian contribution independently.

---

## Raman activity

Raman activity is computed from invariants of the polarizability derivative tensor \( \boldsymbol{\alpha} \):

### Isotropic invariant
\[
\alpha = \frac{1}{3}(\alpha_{xx} + \alpha_{yy} + \alpha_{zz})
\]

### Anisotropic invariant
\[
\beta^2 =
\frac{1}{2}
\left[
(\alpha_{xx} - \alpha_{yy})^2
+ (\alpha_{xx} - \alpha_{zz})^2
+ (\alpha_{yy} - \alpha_{zz})^2
+ 6(\alpha_{xy}^2 + \alpha_{xz}^2 + \alpha_{yz}^2)
\right]
\]

### Raman activity
\[
I_{\mathrm{Raman}} \propto 45\alpha^2 + 7\beta^2
\]

---

## Cartesian-resolved Raman contributions

In VibXYZ, Cartesian-resolved Raman contributions are computed by isolating diagonal tensor components and associated anisotropic terms. This enables directional analysis of Raman activity under anisotropic excitation or confinement.

---

## Notes

- This decomposition is intended for qualitative and comparative analysis
- Total Raman activity is recovered by summing Cartesian contributions
- Orientation-dependent selection rules can be analyzed directly