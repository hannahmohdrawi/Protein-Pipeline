# 🧬 Protein Analysis Pipeline

## 📌 Overview

This project is a Python-based bioinformatics pipeline that reads protein sequences from FASTA files and computes key biochemical properties such as molecular weight, amino acid composition, and isoelectric point (pI). It also visualizes protein charge behavior across pH values.

---

## ⚙️ Installation Requirements

This project requires Python 3.8+ and the following Python libraries:

### 📦 Required Packages

* `numpy` – numerical computations
* `matplotlib` – plotting graphs

---

### 💻 Install Dependencies

Run the following command in your terminal or command prompt:

```bash id="q1p9ab"
pip install numpy matplotlib
```

If that does not work, try:

```bash id="w8k2nc"
python -m pip install numpy matplotlib
```

or on Windows:

```bash id="r3m9xd"
py -m pip install numpy matplotlib
```

---

## 🚀 How to Run

1. Ensure your FASTA file (e.g. `haemoglobin.fasta`) is in the same folder as the script.

2. Run the program:

```bash id="t6v9qp"
python fasta_reader.py
```

---

## 🧪 Features Implemented

### 🔹 FASTA Parsing

* Reads protein sequences from `.fasta` files
* Extracts header and amino acid sequence

### 🔹 Molecular Weight Calculation

* Computes protein molecular weight using amino acid composition
* Accounts for peptide bond formation (water loss correction)

### 🔹 Amino Acid Composition

* Calculates frequency of each amino acid in the sequence

### 🔹 Isoelectric Point (pI)

* Estimates pH at which the protein has net zero charge
* Uses pKa-based charge modelling and numerical search

### 🔹 Charge vs pH Visualization

* Plots net charge across pH 0–14
* Marks estimated pI on the graph
* Saves plot as `charge_vs_ph.png`

---

## 📊 Output Example

```text id="v2k8sd"
Header: >sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens
Sequence length: 142
Molecular Weight: 15126.0
Isoelectric point (pI): 8.72
```

A graph image file will also be generated:

```
charge_vs_ph.png
```

---

## 📂 Project Structure

```text id="m4n8qa"
protein_pipeline/
│
├── fasta_reader.py
├── haemoglobin.fasta
├── charge_vs_ph.png (generated)
└── README.md
```

---

## 🔬 Data Source

Protein sequences can be obtained from:

* UniProt

---

## ⚠️ Notes

* Ensure all dependencies are installed before running the script
* The pI calculation is an approximation based on pKa values
* Designed for educational and learning purposes

---

## 🚀 Future Improvements

* Add hydrophobicity analysis
* Improve pI calculation accuracy
* Add command-line arguments for input FASTA files
* Export results to CSV or JSON
* Build ML-based protein classification module

---
