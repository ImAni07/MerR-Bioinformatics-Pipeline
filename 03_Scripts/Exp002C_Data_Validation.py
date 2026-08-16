# Experiment 002C

# Metadata Cleaning & Validation

# Import Required Libraries
import os
import re
import time
import logging
import datetime
import numpy as np
import pandas as pd
from Bio import SeqIO

# Configuration
PROJECT_DIR = r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work"
INPUT_DIR = os.path.join(PROJECT_DIR, "05_Metadata")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "05_Metadata")
RESULT_DIR = os.path.join(PROJECT_DIR, "04_Results")
LOG_DIR = os.path.join(PROJECT_DIR, "06_Logs")

input_file = os.path.join(INPUT_DIR, "EXP002B_UniProt_Metadata.csv")
fasta_file = os.path.join(PROJECT_DIR, "02_Processed_Data", "MerR_NR90_Galaxy.fasta")
output_file = os.path.join(OUTPUT_DIR, "EXP002C_Clean_Metadata.csv")
validation_report = os.path.join(RESULT_DIR, "EXP002C_Validation_Report.txt")
report_file = os.path.join(RESULT_DIR, "EXP002C_Metadata_Report.txt")
log_file = os.path.join(LOG_DIR, "EXP002C_Log.txt")

# Logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    filemode="w"
)

logging.info("=" * 60)
logging.info("EXP002C Started")
logging.info("=" * 60)

print("=" * 50)
print("Experiment 002C")
print("Metadata Integration, Cleaning & Validation")
print("=" * 50)

print(f"Input File  : {input_file}")
print(f"FASTA File  : {fasta_file}")
print(f"Output File : {output_file}")
print(f"Log File    : {log_file}")

logging.info("Configuration completed successfully.")
logging.info(f"Metadata File : {input_file}")
logging.info(f"FASTA File    : {fasta_file}")
logging.info(f"Output File   : {output_file}")

# Function to Load Metadata
def load_metadata():

    print("\nLoading metadata...")
    print("-" * 50)

    df = pd.read_csv(input_file)

    print(f"Rows    : {len(df)}")
    print(f"Columns : {len(df.columns)}")

    logging.info(f"Metadata loaded successfully.")
    logging.info(f"Rows: {len(df)}")
    logging.info(f"Columns: {len(df.columns)}")

    return df

# Function to Load FASTA Sequences
def load_fasta_sequences():

    print("\nLoading FASTA sequences...")
    print("-" * 50)

    fasta_sequences = {}
    total_sequences = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        total_sequences += 1
        header = record.id

        # UniProt ID is assumed to be the first token
        uniprot_id = header.split("|")[0]
        fasta_sequences[uniprot_id] = str(record.seq)

    print(f"FASTA Sequences Loaded : {total_sequences}")

    logging.info(f"FASTA sequences loaded successfully.")
    logging.info(f"Total FASTA sequences : {total_sequences}")

    return fasta_sequences


# Function to Merge FASTA Sequences
def merge_fasta_sequences(metadata_df, fasta_sequences):

    print("\nIntegrating FASTA sequences...")
    print("-" * 50)

    merged_df = metadata_df.copy()
    merged_df["FASTA_Sequence"] = (merged_df["UniProt_ID"].map(fasta_sequences))
    mapped = merged_df["FASTA_Sequence"].notna().sum()
    
    missing = merged_df["FASTA_Sequence"].isna().sum()
    
    if missing > 0:
        logging.warning(f"{missing} metadata records have no corresponding FASTA sequence.")
    else:
        logging.info("All metadata records successfully matched with FASTA sequences.")

    print(f"Sequences Mapped : {mapped}")
    print(f"Missing Sequences : {missing}")

    logging.info(f"Mapped FASTA sequences : {mapped}")
    logging.info(f"Missing FASTA sequences : {missing}")

    print("FASTA sequences integrated successfully.")

    return merged_df

# Function to Clean Metadata
def clean_metadata(df):

    print("\nCleaning metadata...")
    print("-" * 50)

    cleaned_df = df.copy()

    # Remove leading/trailing whitespace
    object_columns = cleaned_df.select_dtypes(include="object").columns

    for column in object_columns:
        cleaned_df[column] = (
            cleaned_df[column]
                .apply(
                    lambda x: x.strip()
                    if isinstance(x, str)
                    else x
                )
        )

    # Collapse multiple spaces
    for column in object_columns:
        cleaned_df[column] = (
            cleaned_df[column]
                .apply(
                    lambda x: re.sub(
                        r"\s+",
                        " ",
                        x
                    )
                    if isinstance(x, str)
                    else x
                )
        )

    # Convert blank strings to NaN
    cleaned_df.replace(r"^\s*$", np.nan, regex=True, inplace=True)

    # Convert numeric columns
    cleaned_df["Taxonomy_ID"] = pd.to_numeric(cleaned_df["Taxonomy_ID"],errors="coerce")
    cleaned_df["Sequence_Length"] = pd.to_numeric(cleaned_df["Sequence_Length"],errors="coerce")

    # Sort by UniProt ID
    cleaned_df.sort_values(by="UniProt_ID",inplace=True)
    cleaned_df.reset_index(drop=True,inplace=True)

    print("Metadata cleaned successfully.")

    logging.info("Metadata cleaned successfully.")

    return cleaned_df

# Function for Metadata Validation
def validate_metadata(metadata_df):

    print("\nValidating metadata...")
    print("-" * 50)

    validation = {}

    # Basic statistics
    validation["Total_Records"] = len(metadata_df)
    validation["Total_Columns"] = metadata_df.shape[1]

    # Duplicate IDs
    validation["Duplicate_UniProt_IDs"] = metadata_df["UniProt_ID"].duplicated().sum()

    # Missing values
    validation["Missing_UniProt_IDs"] = metadata_df["UniProt_ID"].isna().sum()
    validation["Missing_Protein_Name"] = metadata_df["Protein_Name"].isna().sum()
    validation["Missing_Gene_Name"] = metadata_df["Gene_Name"].isna().sum()
    validation["Missing_Organism"] = metadata_df["Organism"].isna().sum()
    validation["Missing_Taxonomy_ID"] = metadata_df["Taxonomy_ID"].isna().sum()
    validation["Missing_Taxonomic_Lineage"] = metadata_df["Taxonomic_Lineage"].isna().sum()
    validation["Missing_Sequence_Length"] = metadata_df["Sequence_Length"].isna().sum()
    validation["Missing_FASTA_Sequence"] = metadata_df["FASTA_Sequence"].isna().sum()
    validation["FASTA_Sequences_Available"] = (metadata_df["FASTA_Sequence"].notna().sum())
    
    validation["Sequence_Length_Mismatch"] = (
        metadata_df[
            metadata_df["FASTA_Sequence"].notna()
            ]
        .apply(
            lambda row:
                len(row["FASTA_Sequence"]) != row["Sequence_Length"],
                axis=1
                )
        .sum())
    
    validation["Missing_Function"] = metadata_df["Function"].isna().sum()

    # Reviewed / Unreviewed
    validation["Reviewed"] = (metadata_df["Reviewed"].astype(str).str.lower().eq("reviewed").sum())
    validation["Unreviewed"] = (metadata_df["Reviewed"].astype(str).str.lower().ne("reviewed").sum())

    print(f"FASTA Sequences Available : {validation['FASTA_Sequences_Available']}")
    print(f"Missing FASTA Sequences   : {validation['Missing_FASTA_Sequence']}")
    print(f"Length Mismatches         : {validation['Sequence_Length_Mismatch']}")
    
    print("Validation completed successfully.")

    logging.info(f"FASTA sequences available : {validation['FASTA_Sequences_Available']}")
    logging.info(f"Missing FASTA sequences : {validation['Missing_FASTA_Sequence']}")
    logging.info(f"Sequence length mismatches : {validation['Sequence_Length_Mismatch']}")
    
    logging.info("Metadata validation completed successfully.")

    return validation

# Function to Export Clean Metadata
def export_clean_metadata(metadata_df):

    print("\nExporting Clean Metadata...")
    print("-" * 50)

    metadata_df.to_csv(output_file,index=False)
    exported_df = pd.read_csv(output_file)
    
    # Verify FASTA_Sequence column exists
    if "FASTA_Sequence" not in exported_df.columns:
        raise RuntimeError("FASTA_Sequence column missing in exported CSV.")
    
    logging.info("FASTA_Sequence column verified in exported metadata.")

    print(f"Records Exported : {len(exported_df)}")

    if len(exported_df) == len(metadata_df):
        print("✓ CSV export verified.")
        logging.info("Clean metadata export verified.")

    else:
        raise RuntimeError("CSV export verification failed.")

    print("\nClean metadata successfully exported to:")
    print(output_file)

    logging.info(f"Clean metadata exported: {output_file}")

# Function to Export Validation Report
def export_validation_report(validation_summary):

    print("\nExporting Validation Report...")
    print("-" * 50)

    with open(validation_report, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write("EXP002C Validation Report\n")
        f.write("=" * 60 + "\n\n")

        for key, value in validation_summary.items():
            f.write(f"{key:<30}: {value}\n")

    print("Validation report generated successfully.")
    print(validation_report)

    logging.info("Validation report generated successfully.")
    logging.info("Validation report includes FASTA integration statistics.")
    logging.info(f"Validation report saved to: {validation_report}")

# Function to Generate Experiment Report
def generate_experiment_report(metadata_df, validation_summary, execution_time):

    print("\nGenerating Experiment Report...")
    print("-" * 50)

    report_lines = []

    report_lines.append("=" * 60)
    report_lines.append("Experiment : EXP002C")
    report_lines.append("Title      : Metadata Integration, Cleaning & Validation")
    report_lines.append("=" * 60)

    report_lines.append("")
    report_lines.append(f"Date : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")
    report_lines.append("INPUT")
    report_lines.append("-" * 60)
    report_lines.append(f"Input File           : {input_file}")
    report_lines.append(f"FASTA File           : {fasta_file}")
    report_lines.append(f"Input Records        : {validation_summary['Total_Records']}")
    report_lines.append("")
    report_lines.append("PROCESSING")
    report_lines.append("-" * 60)
    report_lines.append(f"Records Processed    : {len(metadata_df)}")
    report_lines.append(f"Metadata Columns     : {metadata_df.shape[1]}")
    report_lines.append("New Column Added      : FASTA_Sequence")
    report_lines.append(f"Duplicate IDs        : {validation_summary['Duplicate_UniProt_IDs']}")
    report_lines.append(f"Missing Gene Names   : {validation_summary['Missing_Gene_Name']}")
    report_lines.append(f"Missing Functions    : {validation_summary['Missing_Function']}")
    report_lines.append("")
    report_lines.append("FASTA INTEGRATION")
    report_lines.append("-" * 60)
    report_lines.append(f"FASTA Sequences Available : " f"{validation_summary['FASTA_Sequences_Available']}")
    report_lines.append(f"Missing FASTA Sequences   : "f"{validation_summary['Missing_FASTA_Sequence']}")
    report_lines.append(f"Length Mismatches         : " f"{validation_summary['Sequence_Length_Mismatch']}")
    report_lines.append("")
    report_lines.append("OUTPUT")
    report_lines.append("-" * 60)
    report_lines.append(f"Clean Metadata CSV   : {output_file}")
    report_lines.append("Metadata Status       : FASTA sequences successfully integrated")
    report_lines.append(f"Validation Report    : {validation_report}")
    report_lines.append(f"Log File             : {log_file}")
    report_lines.append("")
    report_lines.append("PIPELINE")
    report_lines.append("-" * 60)
    report_lines.append(f"Execution Time       : {round(execution_time,2)} seconds")
    report_lines.append("Overall Status        : SUCCESS")
    report_lines.append("")
    report_lines.append("=" * 60)
    report_lines.append("End of Report")
    report_lines.append("=" * 60)

    with open(report_file, "w", encoding="utf-8") as f:
        for line in report_lines:
            f.write(line + "\n")

    print("Experiment report generated successfully.")
    print(report_file)

    logging.info("Experiment report generated successfully.")
    logging.info("Experiment report includes FASTA integration summary.")
    logging.info(f"Experiment report saved to: {report_file}")

# Main Function
def main():

    start_time = time.time()

    metadata_df = load_metadata()
    fasta_sequences = load_fasta_sequences()
    metadata_df = merge_fasta_sequences(metadata_df, fasta_sequences)
    metadata_df = clean_metadata(metadata_df)
    validation_summary = validate_metadata(metadata_df)
    
    export_clean_metadata(metadata_df)
    export_validation_report(validation_summary)

    execution_time = time.time() - start_time

    generate_experiment_report(metadata_df, validation_summary, execution_time)

    print("\nExecution Summary")
    print("-" * 50)
    print(f"Metadata Records       : {len(metadata_df)}")
    print(f"Metadata Columns       : {metadata_df.shape[1]}")
    print(f"FASTA Sequences Added  : "f"{validation_summary['FASTA_Sequences_Available']}")
    print(f"Missing FASTA Sequence : "f"{validation_summary['Missing_FASTA_Sequence']}")
    
    print("\n")
    print("=" * 50)
    print("EXP002C COMPLETED SUCCESSFULLY")
    print("=" * 50)

    logging.info("Execution Summary")
    logging.info(f"Metadata Records             : {len(metadata_df)}")
    logging.info(f"Metadata Columns             : {metadata_df.shape[1]}")
    logging.info(f"FASTA Sequences Added        : "f"{validation_summary['FASTA_Sequences_Available']}")
    logging.info(f"Missing FASTA Sequences      : "f"{validation_summary['Missing_FASTA_Sequence']}")
    logging.info(f"Sequence Length Mismatches   : "f"{validation_summary['Sequence_Length_Mismatch']}")
    
    logging.info("=" * 50)
    logging.info("EXP002C COMPLETED SUCCESSFULLY")
    logging.info("=" * 50)

# Execute the Main Function
if __name__ == "__main__":

    try:
        main()

    except Exception as e:
        print("\n")
        print("=" * 50)
        print("EXP002C pipeline terminated due to an unexpected error.")
        print("=" * 50)
        print(e)

        logging.exception("Pipeline Failed")

        raise