# 🧬 Protein Analysis Pipeline

## 📌 Description

This project is a Python-based bioinformatics pipeline that reads protein sequences from FASTA files and computes key biochemical properties such as molecular weight, amino acid composition, isoelectric point (pI), and hydrophobicity. It also visualizes protein charge behavior across pH values and includes a machine learning module for protein classification and prediction.

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
                                                           │
                                                           ▼
                                                [Feature Dataset (CSV)]
                                                           │
                                                           ▼
                                                   [ML Model Training]
                                                           │
                                                           ▼
                                              [Prediction + Evaluation]
                                                           │
                                                           ▼
                                               [Saved Model (.pkl file)]
```

---

## ⚙️ Installation Requirements

This project requires Python 3.8+ and the following Python libraries:

### 📦 Required Packages

* `numpy` – numerical computations
* `matplotlib` – plotting graphs
* `scikit-learn` – machine learning models and evaluation
* `pandas` – dataset creation and CSV handling
* `joblib` – model persistence (saving/loading `.pkl` files)

---

### 💻 Install Dependencies

Run the following command:

```bash
pip install numpy matplotlib scikit-learn pandas joblib
```

Alternative options:

```bash
python -m pip install numpy matplotlib scikit-learn pandas joblib
```

```bash
py -m pip install numpy matplotlib scikit-learn pandas joblib
```

---

## 🚀 How to Run

1. Place your FASTA file inside:

```text
Protein FASTA/
```

2. Run the main analysis pipeline:

```bash
python fasta_reader.py
```

3. Generate a dataset from multiple FASTA files:

```bash
python build_dataset.py
```

4. Train and evaluate a machine learning model:

```bash
python train_model.py
```

5. Predict protein properties for a new sequence:

```bash
python predict.py --fasta Protein\ FASTA/new_protein.fasta
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

### 🔹 Hydrophobicity

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
* Exports combined dataset as CSV in `Dataset/`

---

### 🔹 Machine Learning Module

#### 📊 Dataset Generation

* Processes multiple FASTA files in batch
* Combines extracted features into a structured dataset (`dataset.csv`)
* Supports labelled data for supervised learning (e.g. protein function, cellular location)

#### 🧮 Model Training

Supports the following classification models:

* **Random Forest** – robust ensemble method, handles non-linear relationships
* **Logistic Regression** – fast linear baseline, interpretable coefficients

```python
# Example usage
from train_model import train
model, report = train(dataset="Dataset/dataset.csv", model_type="random_forest", label="location")
```

#### 🔮 Prediction

* Accepts a new FASTA file as input
* Extracts features automatically
* Returns predicted class and confidence score

```python
from predict import predict_protein
result = predict_protein("Protein FASTA/new_protein.fasta", model_path="Models/rf_model.pkl")
print(result)
# {'prediction': 'membrane', 'confidence': 0.87}
```

#### 📈 Model Evaluation

* Train/test split (default: 80/20)
* Reports accuracy, precision, recall, and F1-score
* Outputs a confusion matrix plot to `Protein Graphs/`

```text
Classification Report:
              precision    recall  f1-score   support

    membrane       0.91      0.88      0.89        50
    cytosolic      0.87      0.90      0.88        48

    accuracy                           0.89        98
```

#### 💾 Model Persistence

* Trained models are saved as `.pkl` files in `Models/`
* Reload at any time for inference without retraining

```python
import joblib
model = joblib.load("Models/rf_model.pkl")
```

---

## 📊 Output Example

```text
Header: >sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens
Sequence length: 142
Molecular Weight: 15126.0
Isoelectric point (pI): 8.72
Hydrophobicity: -0.45
Predicted location: cytosolic (confidence: 0.92)
```

---

### 📈 Generated Files

**Graph:**

```text
Protein Graphs/HBA_HUMAN_charge_vs_ph.png
Protein Graphs/confusion_matrix.png
```

**Features:**

```text
Features/HBA_HUMAN_features.json
```

**Dataset:**

```text
Dataset/dataset.csv
```

**Model:**

```text
Models/rf_model.pkl
Models/lr_model.pkl
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
├── build_dataset.py             # Batch processing + CSV dataset builder
├── train_model.py               # ML model training and evaluation
├── predict.py                   # Predict properties for new sequences
│
├── Protein FASTA/               # Input folder (raw protein sequences)
│   ├── haemoglobin.fasta
│   ├── example_protein.fasta
│   └── ...
│
├── Protein Graphs/              # Output folder (generated plots)
│   ├── HBA_HUMAN_charge_vs_ph.png
│   ├── confusion_matrix.png
│   └── ...
│
├── Features/                    # Extracted ML features (per protein)
│   ├── HBA_HUMAN_features.json
│   └── ...
│
├── Dataset/                     # Combined feature datasets
│   └── dataset.csv
│
├── Models/                      # Saved trained models
│   ├── rf_model.pkl
│   └── lr_model.pkl
│
├── README.md                    # Project documentation
```

---

## 🔬 Data Source

Protein sequences can be obtained from:

* UniProt

Labelled datasets for ML training can be sourced from:

* UniProt (functional annotations)
* SwissProt (curated protein function and location data)

---

## ⚠️ Notes

* Ensure all dependencies are installed before running the script
* The pI calculation is an approximation based on pKa values
* Hydrophobicity is based on the Kyte–Doolittle scale
* ML model performance depends on dataset size and label quality
* Designed for educational and learning purposes

---

## 🚀 Future Improvements

### 🔹 Additional Biochemical Features

* Aromaticity
* Instability index
* Secondary structure prediction

### 🔹 Extended ML Capabilities

* Deep learning models (e.g. protein embeddings with ESM or ProtBERT)
* Multi-label classification (predict multiple properties simultaneously)
* Cross-validation and hyperparameter tuning

### 🔹 Tooling

* Command-line interface (CLI) for flexible input/output options
* Web interface for non-technical users
* Export reports as PDF

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

**FASTA → Analysis → Visualization → Feature Extraction → Machine Learning pipeline**

and provides a strong foundation for bioinformatics and machine learning applications.