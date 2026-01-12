# VibXYZ

**Cartesian (X, Y, Z) decomposition of IR and Raman vibrational intensities from Gaussian frequency calculations.**

Gaussian reports total IR and Raman activities per vibrational mode, but does not provide their individual Cartesian contributions.  
VibXYZ enables extraction of the **X, Y, and Z components** of vibrational intensities by post-processing Gaussian frequency outputs.

This tool is useful for orientation-dependent vibrational analysis, anisotropic spectroscopy, and surface–molecule interaction studies.

---

## 🔬 What this code does

VibXYZ parses Gaussian `freq` calculations to:
- Extract vibrational eigenvectors
- Read dipole derivatives (IR)
- Read polarizability derivatives (Raman)
- Compute:
  - Total IR / Raman activities
  - Cartesian-resolved contributions (**X, Y, Z**)

---

## ⚙️ Gaussian requirements

Your Gaussian frequency calculation **must include**:

```text
#P functional/basis nosymm freq=(raman,savenormalmodes) iop(7/33=1)
```
### Important notes:

- #P is required to print dipole and polarizability derivatives
- iop(7/33=1) enables polarizability derivative output
- No atoms may be frozen (-1 constraints will cause the parser to fail)

## ▶️ How to run

- Run a Gaussian frequency calculation
- Place gaus_vib_xyz_raman_hacky_v2.py in the same directory as the Gaussian output
- Rename the Gaussian output file to:
```text
gaus.out.log
```
- Run:
```text
python gaus_vib_xyz_raman_hacky_v2.py
```

## 📊 Example notebook

The notebook below demonstrates how to load and visualize XYZ-decomposed spectra:
```text
Plot_XYZ-decomp-Gaussian_example.ipynb
```

## 🧠 Theory background

Brief analytical background is provided in:
```text
docs/theory.md
```

This includes:
- IR intensity decomposition
- Raman activity invariants
- Cartesian-resolved contributions

## 🧾 Credits

This code is adapted and extended from earlier Raman post-processing tools developed by:
- Alexandr Fonari (Georgia Tech)
- Shannon Stauffer (UT Austin)

Gaussian-specific parsing, XYZ decomposition, and visualization workflows were implemented by the authors of VibXYZ.