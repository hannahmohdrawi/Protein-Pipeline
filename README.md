# 🧬 Protein Analysis Pipeline

## 📌 Description

This project is a Python-based bioinformatics pipeline that reads protein sequences from FASTA files and computes key biochemical properties such as molecular weight, amino acid composition, isoelectric point (pI), and hydrophobicity. It also visualizes protein charge behavior across pH values.

---

## 🔄 Pipeline Overview

```text
[FASTA File]
     │
     ▼
[read_fasta()]
     │
     ▼
[Protein Sequence]
     │
     ├───────────────┬────────────────────┬────────────────────┐
     ▼               ▼                    ▼                    ▼
[Molecular Weight] [AA Composition] [Hydrophobicity] [Charge Calculation]
                                                           │
                                                           ▼
                                                   [pI Estimation]
                                                           │
                                                           ▼
                                               [Charge vs pH Plot]
                                                           │
                                                           ▼
                           [Saved Graph (PNG)] + [Saved Features (JSON)]
```

---

## ⚙️ Installation Requirements

This project requires Python 3.8+ and the following Python libraries:

### 📦 Required Packages

* `numpy` – numerical computations
* `matplotlib` – plotting graphs

---

### 💻 Install Dependencies

Run the following command:

```bash
pip install numpy matplotlib
```

Alternative options:

```bash
python -m pip install numpy matplotlib
```

```bash
py -m pip install numpy matplotlib
```

---

## 🚀 How to Run

1. Place your FASTA file inside:

```text
Protein FASTA/
```

2. Run the program:

```bash
python fasta_reader.py
```

---

## 🧪 Features Implemented

### 🔹 FASTA Parsing

* Reads protein sequences from `.fasta` files
* Extracts header and amino acid sequence

---

### 🔹 Molecular Weight Calculation

* Computes protein molecular weight using amino acid composition
* Accounts for peptide bond formation (water loss correction)

---

### 🔹 Amino Acid Composition

* Calculates frequency of each amino acid in the sequence

---

### 🔹 Isoelectric Point (pI)

* Estimates pH at which the protein has net zero charge
* Uses pKa-based charge modelling and numerical search

---

### 🔹 Hydrophobicity (NEW)

* Calculates average hydrophobicity using the **Kyte–Doolittle scale**
* Indicates whether a protein is:

  * **Hydrophobic** (positive values → membrane-associated tendency)
  * **Hydrophilic** (negative values → soluble proteins)

---

### 🔹 Charge vs pH Visualization

* Plots net charge across pH 0–14
* Marks estimated pI on the graph
* Saves plot as a PNG file in `Protein Graphs/`

---

### 🔹 Feature Extraction (ML-Ready)

* Converts protein sequences into structured numerical features:

  * Length
  * Molecular weight
  * pI
  * Hydrophobicity
  * Amino acid frequencies (20 features)

* Saves features as JSON files in `Features/`

---

## 📊 Output Example

```text
Header: >sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens
Sequence length: 142
Molecular Weight: 15126.0
Isoelectric point (pI): 8.72
hydrophobicity: -0.45
```

---

### 📈 Generated Files

**Graph:**

```text
Protein Graphs/HBA_HUMAN_charge_vs_ph.png
```

**Features:**

```text
Features/HBA_HUMAN_features.json
```

Example JSON:

```json
{
    "length": 142,
    "molecular_weight": 15126.0,
    "pI": 8.72,
    "hydrophobicity": -0.45,
    "freq_A": 0.07,
    "freq_R": 0.03
}
```

---

## 📂 Project Structure

```text
protein_pipeline/
│
├── fasta_reader.py              # Main analysis script
│
├── Protein FASTA/               # Input folder (raw protein sequences)
│   ├── haemoglobin.fasta
│   ├── example_protein.fasta
│   └── ...
│
├── Protein Graphs/              # Output folder (generated plots)
│   ├── HBA_HUMAN_charge_vs_ph.png
│   └── ...
│
├── Features/                    # Extracted ML features
│   ├── HBA_HUMAN_features.json
│   └── ...
│
├── README.md                    # Project documentation
```

---

## 🔬 Data Source

Protein sequences can be obtained from:

* UniProt

---

## ⚠️ Notes

* Ensure all dependencies are installed before running the script
* The pI calculation is an approximation based on pKa values
* Hydrophobicity is based on the Kyte–Doolittle scale
* Designed for educational and learning purposes

---

## 🚀 Future Improvements

### 🔹 Machine Learning Integration (Planned)

Extend the pipeline to include machine learning for protein classification and prediction tasks.

**Proposed workflow:**

```text
Protein FASTA → Feature Extraction → Dataset → ML Model → Prediction
```

**Planned additions:**

* 📊 **Dataset generation**

  * Process multiple FASTA files
  * Combine extracted features into a structured dataset (CSV)

* 🧮 **Model training**

  * Use models such as:

    * Random Forest
    * Logistic Regression
  * Train on labelled protein data (e.g. function, location)

* 🔮 **Prediction capability**

  * Predict properties such as:

    * Protein function
    * Cellular location (cytosolic vs membrane)
    * Protein class/type

* 📈 **Model evaluation**

  * Accuracy, precision, recall
  * Train/test split validation

* 💾 **Model persistence**

  * Save trained models for reuse (`.pkl` files)

---

### 🔹 Additional Enhancements

* Batch processing of multiple FASTA files
* Export features to CSV for dataset creation
* Add more biochemical properties:

  * Aromaticity
  * Instability index
  * Secondary structure prediction
* Command-line interface (CLI) for flexible input

---

### 🧠 Long-Term Vision

Transform this project into a complete:

```text
Protein Analysis → Feature Engineering → Machine Learning Platform
```

capable of supporting real-world bioinformatics workflows.


---

## 📌 Summary

This project demonstrates a complete:

**FASTA → Analysis → Visualization → Feature Extraction pipeline**

and provides a strong foundation for bioinformatics and machine learning applications.
