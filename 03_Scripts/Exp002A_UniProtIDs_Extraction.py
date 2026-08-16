# Extract UniProt IDs

"""
    Experiment No.: 002A
    Title: Extraction of UniProt IDs from MerR Family Regulatory Proteins
    Author: Anirban Majumder
    Supervisor: Dr. Rudra Prasad Saha
    Date: 4th August, 2026

    Objective: To extract UniProt IDs from the MerR FASTA dataset generated after CD-HIT

    Input: MerR_NR90_Galaxy.fasta

    Output: EXP002A_UniProtIDs_Output.txt
"""

# Import Required Libraries
from Bio import SeqIO
from pathlib import Path

# Configuration
input_fasta = Path(r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work\02_Processed_Data\MerR_NR90_Galaxy.fasta")
output_file = Path(r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work\05_Metadata\EXP002A_UniProt_IDs.txt")

# Function to Extract UniProt IDs from a FASTA File
def extract_Uniprot_IDs (input_file):
    
    uniprot_IDs = []
    
    for record in SeqIO.parse(input_file, "fasta"):

        accession = record.id.split("|")[0]

        uniprot_IDs.append(accession)

    return sorted(set(uniprot_IDs))

# Function to Save UniProt IDs to a Text File
def save_Uniprot_IDs (ids, output_file):
    
    with open(output_file, "w") as file:

        for accession in ids:

            file.write(accession + "\n")

# Main Function
def main ():
    
    print("=" * 60)
    print("Experiment 002A")
    print("UniProt ID Extraction")
    print("=" * 60)

    ids = extract_Uniprot_IDs(input_fasta)

    save_Uniprot_IDs(ids, output_file)

    print(f"\nInput FASTA File : {input_fasta.name}")
    print(f"Total IDs        : {len(ids)}")
    print(f"Output File      : {output_file.name}")

    print("\nStatus : Completed Successfully")

    print("=" * 60)

# Execute the Main Function
if __name__ == "__main__":
    main()

# Output:
"""
============================================================
Experiment 002A
UniProt ID Extraction
============================================================

Input FASTA File : MerR_NR90_Galaxy.fasta
Total IDs        : 14299
Output File      : EXP002A_UniProt_IDs.txt

Status : Completed Successfully
============================================================
"""

# UniProt IDs extracted from the MerR FASTA dataset have been saved to the output file "EXP002A_UniProt_IDs.txt".
# UniProt IDs of all 14299 sequences have been successfully extracted and saved to the output file.