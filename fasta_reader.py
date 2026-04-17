import re
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter


# Reading protein file and returning sequence
def read_fasta(file_path):
    header = ""
    sequence = ""

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                header = line
            else:
                sequence += line 

    return header, sequence


# Extract clean protein name from FASTA header
def get_protein_name(header):
    name = header.split()[1] if " " in header else "protein"
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    return name


# Amino acid weights
aa_weights = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1,
    'C': 121.2, 'E': 147.1, 'Q': 146.2, 'G': 75.1,
    'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
    'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1,
    'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
}


# Calculating molecular weight
def molecular_weight(seq):
    total = sum(aa_weights.get(aa, 0) for aa in seq)
    water_mass = 18.015
    return total - (len(seq) - 1) * water_mass


# pKa Values
pKa = {
    'Cterm': 2.34,
    'Nterm': 9.69,
    'C': 8.33,
    'D': 3.86,
    'E': 4.25,
    'H': 6.00,
    'K': 10.53,
    'R': 12.48,
    'Y': 10.07
}


# Count charged residues once
def get_counts(seq):
    return {aa: seq.count(aa) for aa in "CDEHKRY"}


# Net charge calculation
def net_charge_from_counts(counts, pH):

    pos = (
        (10**pKa['Nterm'] / (10**pKa['Nterm'] + 10**pH)) +
        counts['K'] * (10**pKa['K'] / (10**pKa['K'] + 10**pH)) +
        counts['R'] * (10**pKa['R'] / (10**pKa['R'] + 10**pH)) +
        counts['H'] * (10**pKa['H'] / (10**pKa['H'] + 10**pH))
    )

    neg = (
        (10**pH / (10**pKa['Cterm'] + 10**pH)) +
        counts['D'] * (10**pH / (10**pKa['D'] + 10**pH)) +
        counts['E'] * (10**pH / (10**pKa['E'] + 10**pH)) +
        counts['C'] * (10**pH / (10**pKa['C'] + 10**pH)) +
        counts['Y'] * (10**pH / (10**pKa['Y'] + 10**pH))
    )

    return pos - neg


# Find pI
def calculate_pI(counts):
    best_pH = 0
    min_charge = float("inf")

    for pH in [x * 0.01 for x in range(0, 1400)]:
        charge = net_charge_from_counts(counts, pH)

        if abs(charge) < min_charge:
            min_charge = abs(charge)
            best_pH = pH

    return best_pH


# Amino acid composition
def aa_composition(seq):
    counts = Counter(seq)
    total = len(seq)
    return {aa: counts[aa]/total for aa in counts}


# Plot charge curve
def plot_charge_curve(counts, pI, protein_name):

    pH_values = np.linspace(0, 14, 200)
    charges = [net_charge_from_counts(counts, pH) for pH in pH_values]

    plt.plot(pH_values, charges)
    plt.axhline(0)
    plt.axvline(pI)

    plt.xlabel("pH")
    plt.ylabel("Net Charge")
    plt.title(f"Charge vs pH Curve: {protein_name}")

    plt.text(pI, 0, f" pI ≈ {pI:.2f}")

    filename = f"{protein_name}_charge_vs_ph.png"
    plt.savefig(filename)
    plt.show()

    print(f"Graph saved as: {filename}")

# Execution
if __name__ == "__main__":

    header, seq = read_fasta("haemoglobin.fasta")
    counts = get_counts(seq)

    protein_name = get_protein_name(header)

    print("Header:", header)
    print("Sequence length:", len(seq))
    print("Molecular Weight:", molecular_weight(seq))

    pI = calculate_pI(counts)
    print("Isoelectric point (pI):", pI)

    plot_charge_curve(counts, pI, protein_name)