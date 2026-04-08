
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


# Haemoglobin check 
header, seq = read_fasta("haemoglobin.fasta")

print("Header", header)
print("Sequence length:", len(seq))
print("First 50 aa:", seq[:50])