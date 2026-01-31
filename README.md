# 🧬 Molecular Dynamics: Systems Building & Simulation

We have extensive experience in building and simulating all-atom molecular dynamics systems. This repository shares tutorials at different difficulty levels, execution scripts for each step, MD evaluations, and various analyses relevant to each stage of development.

---

## 📑 Table of Contents
* [Peptides Construction](#peptides-construction-using-namd-193)
* [Membrane Preparation](#membrane-using-charmm-gui)
* [Contact](#contact)

---

## 🧪 Peptides (Construction Using NAMD 1.9.3)

### 🔹 Peptide Construction Steps in VMD
Detailed workflow for assembling the primary structure and initial coordinates.

### 🔹 Saving the PDB File
Best practices for exporting coordinates ensuring compatibility with NAMD.

### 🔹 Generating the PSF File and Adjusting N- and C-Termini
* **Step 1:** Initial topology generation.
* **Step 2:** Patching N-terminus (e.g., Acetylation).
* **Step 3:** Patching C-terminus (e.g., Amidation).
* **Step 4:** Final verification and charge neutrality.

> 📚 [References - Peptides](./path/to/peptide-refs.md)

---

## 🦠 Membrane (Using CHARMM-GUI)

### 🔹 System Preparation in CHARMM-GUI
Setting up lipids, hydration, and ion concentration (0.15 M NaCl).

### 🔹 Downloading and Choosing the System Files
Which files are essential for NAMD and how to organize the `toppar` directory.

### 🔹 Energy Minimization and Equilibration
* **Minimization:** Removing steric clashes.
* **Equilibration:** NVT and NPT ensembles with area-constant constraints.

> 📚 [References - Membrane](./path/to/membrane-refs.md)

---

## 📬 Contact

Any questions or suggestions, please contact us!

| Alexandre Suman de Araújo | Responsible Researcher | alexandre.suman@unesp.br |

| **SSiLiq Team** | Research & Development | [GitHub](https://github.com/SSiLiq) | [Instragram](https://www.instagram.com/ssiliq_ibilce/) |

---
*Developed with ❤️ for the Scientific Community.*
