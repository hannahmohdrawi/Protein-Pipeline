# 🧬 Protein Analysis Pipeline

## 📌 Overview

This project is a simple bioinformatics pipeline written in Python that reads a protein sequence from a FASTA file and computes key biochemical properties.

It is designed as a learning project to build foundational skills in:

* Python programming
* File handling
* Bioinformatics data formats
* Protein property calculations

---

## 🚀 Features Implemented

### 1. FASTA File Parsing

* Reads a protein sequence from a `.fasta` file
* Extracts:

  * Header (metadata line starting with `>`)
  * Amino acid sequence

---

### 2. Sequence Analysis

* Computes **sequence length**
* Handles multi-line FASTA sequences correctly

---

### 3. Molecular Weight Calculation

* Calculates protein molecular weight based on amino acid composition
* Adjusts for peptide bond formation by subtracting water mass

---

### 4. Isoelectric Point (pI) Calculation

* Estimates the pH at which the protein has **net zero charge**
* Uses:

  * pKa values for ionizable groups
  * Henderson–Hasselbalch-based charge calculations
* Scans pH range (0–14) to find approximate pI

---

## 🧪 Example Output

```text
Header: >sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens
Sequence length: 142
Molecular Weight: 15126.0
Isoelectric point (pI): 8.72
```

---

## 📂 Project Structure

```text
protein_pipeline/
│
├── main.py                # Main script
├── haemoglobin.fasta     # Example input file
└── README.md             # Project documentation
```

---

## ⚙️ How to Run

### 1. Install Python

Make sure Python is installed and added to PATH.

### 2. Place FASTA File

Ensure your FASTA file (e.g., `haemoglobin.fasta`) is in the same directory as the script.

### 3. Run the Script

```bash
python main.py
```

---

## 🧠 How It Works

### Pipeline Flow

```text
FASTA file
   ↓
read_fasta()
   ↓
sequence
   ↓
molecular_weight()
calculate_pI()
   ↓
results printed
```

---

## 🔬 Key Concepts Used

* File I/O in Python (`with open`)
* String manipulation
* Dictionaries for biochemical data
* Numerical approximation methods
* Basic protein chemistry:

  * Amino acid properties
  * pKa values
  * Net charge calculations

---

## ⚠️ Limitations

* Assumes a **single sequence** per FASTA file
* pI calculation is an **approximation** (not as precise as tools like ExPASy ProtParam)
* No error handling for invalid sequences yet

---

## 🔧 Future Improvements

* Add amino acid composition analysis
* Improve FASTA parser for multiple sequences
* Optimize pI calculation (binary search instead of brute force)
* Add hydrophobicity calculations
* Build a command-line interface
* Validate results against biological databases

---

## 📚 Data Source

Protein sequences can be obtained from:

* UniProt

---

## 👨‍💻 Author

This project was created as part of learning bioinformatics and Python by Hannah Mohd-Rawi.

---
