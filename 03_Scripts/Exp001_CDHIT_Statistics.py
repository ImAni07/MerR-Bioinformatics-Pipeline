# CD HIT Statistics Script

"""
    Experiment No.: 001
    Title: Generation of a Non - Redundant Dataset of MerR Family Regulatory Proteins using CD-HIT
    Author: Anirban Majumder
    Supervisor: Dr. Rudra Prasad Saha
    Date: 4th August, 2026
    
    Objective: To generate a non-redundant protein sequence dataset of the MerR family regulatory proteins by clustering highly similar sequences using CD-HIT at a sequence identity threshold of 90%.
    
    Input:
        protein-matching-PF00376.fasta
        MerR_NR90_Galaxy.fasta
        MerR_NR90_Galaxy.clstr.txt
    
    Output: EXP001_CDHIT_Statistics_Output.txt
"""

# Import Required Libraries
from Bio import SeqIO
from pathlib import Path

"""
    Experiment 001
    Generation of a Non - Redundant Dataset of MerR Family Regulatory Proteins using CD-HIT
"""

# Load Files

# Input FASTA File (Raw Sequences)
raw_fasta = r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work\01_RawData\protein-matching-PF00376.fasta"

# Output FASTA File (Non-Redundant Sequences)
nr_fasta = r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work\02_Processed_Data\MerR_NR90_Galaxy.fasta"

# Output Cluster File (CD-HIT Cluster Information)
cluster_file = r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work\02_Processed_Data\MerR_NR90_Galaxy.clstr.txt"

# Function to Count the Number of Sequences in a FASTA File
def count_fasta_sequences(fasta_file):
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

# Function to Count the Number of Clusters in a CD-HIT Cluster File
def count_clusters(clstr_file):

    clusters = 0

    with open(clstr_file, "r") as file:

        for line in file:

            if line.startswith(">Cluster"):

                clusters += 1

    return clusters

# Main Function
def main ():
    
    # Count the Number of Sequences in the Input FASTA File
    input_sequences = count_fasta_sequences(raw_fasta)

    # Count the Number of Sequences in the Output FASTA File
    output_sequences = count_fasta_sequences(nr_fasta)

    # Count the Number of Clusters in the CD-HIT Cluster File
    number_of_clusters = count_clusters(cluster_file)

    # Calculate the Percentage of Redundancy Removed
    redundancy_removed = ((input_sequences - output_sequences) / input_sequences) * 100

    # Print the Statistics
    print("=" * 50)
    print("CD-HIT STATISTICS")
    print("=" * 50)

    print(f"Input Sequences           : {input_sequences}")
    print(f"Output Sequences          : {output_sequences}")
    print(f"Number of Clusters        : {number_of_clusters}")
    print(f"Redundancy Removed (%)    : {redundancy_removed:.2f}%")

    print("=" * 50)
    
    # Save the Output
    output_file = Path(r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work\04_Results\EXP001_CDHIT_Statistics_Output.txt")

    with open(output_file, "w") as file:

        file.write("=" * 50 + "\n")
        file.write("CD-HIT STATISTICS\n")
        file.write("=" * 50 + "\n")

        file.write(f"Input Sequences          : {input_sequences}\n")
        file.write(f"Representative Sequences : {output_sequences}\n")
        file.write(f"Number of Clusters       : {number_of_clusters}\n")
        file.write(f"Redundancy Removed (%)   : {redundancy_removed:.2f}%\n")
        
        file.write("=" * 50 + "\n")

    print(f"\nStatistics saved to:\n{output_file}")

# Execute the Main Function
if __name__ == "__main__":
    main()

# Output / Result
"""
==================================================
CD-HIT STATISTICS
==================================================
Input Sequences           : 22945
Output Sequences          : 14299
Number of Clusters        : 14299
Redundancy Removed (%)    : 37.68%
==================================================
"""