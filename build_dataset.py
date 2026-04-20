import os
import csv

def build_dataset(fasta_folder, output_file):
    data = []

    for filename in os.listdir(fasta_folder):
        if filename.endswith(".fasta"):

            file_path = os.path.join(fasta_folder, filename)

            header, seq = read_fasta(file_path)
            protein_name = get_protein_name(header)

            features = extract_features(seq)
            features["protein_name"] = protein_name

            data.append(features)

    # Get all feature keys
    keys = data[0].keys()

    # Write to CSV
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(data)

    print(f"Dataset saved to {output_file}")

if __name__ == "__main__":
build_dataset("Protein FASTA", "protein_dataset.csv")