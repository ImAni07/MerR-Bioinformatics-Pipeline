# Retrieve UniProt Metadata

"""
Experiment: EXP002B
Title: Download UniProt Metadata for Non-Redundant MerR Dataset
Author: Anirban Majumder
Supervisor: Dr. Rudra Prasad Saha

Objective: Download metadata for all representative MerR proteins from the UniProt REST API using batch retrieval.

Input: EXP002A_UniProt_IDs.txt

Output:
    EXP002B_UniProt_Metadata.csv
    EXP002B_Metadata_Report.txt
    EXP002B_API_Log.txt
"""

# 1. Configuration and Setup

# Import Libraries
from pathlib import Path
import logging
import time
import datetime
import re
from urllib import response
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from urllib.parse import urlparse, parse_qs
from io import StringIO

# Project Directories
project_root = Path(r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work")

raw_data = project_root / "01_RawData"
processed_data = project_root / "02_Processed_Data"
scripts = project_root / "03_Scripts"
results = project_root / "04_Results"
metadata = project_root / "05_Metadata"
logs = project_root / "06_Logs"

# Input Files
input_file = metadata / "EXP002A_UniProt_IDs.txt"

# Output Files
output_file = metadata / "EXP002B_UniProt_Metadata.csv"
failed = metadata / "EXP002B_Failed_IDs.txt"
report_file = results / "EXP002B_Metadata_Report.txt"
log_file = logs / "EXP002B_API_Log.txt"

# Create Output Directories

metadata.mkdir(exist_ok=True)
results.mkdir(exist_ok=True)
logs.mkdir(exist_ok=True)

# Configure Logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logging.info("=" * 50)
logging.info("EXP002B Started")
logging.info("=" * 50)

# UniProt REST API Configuration
BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 200
TIMEOUT = 60
MAX_RETRIES = 5
BACKOFF_FACTOR = 1

# Configure HTTP Session
retry_strategy = Retry(
    total=MAX_RETRIES,
    backoff_factor=BACKOFF_FACTOR,
    status_forcelist=[500, 502, 503, 504],
)

session = requests.Session()
adapter = HTTPAdapter(max_retries=retry_strategy)
session.mount("https://", adapter)

# Display Project Information
print("=" * 50)
print("Experiment 002B")
print("UniProt Metadata Retrieval")
print("=" * 50)

print(f"Input IDs      : {input_file.name}")
print(f"Output CSV     : {output_file.name}")
print(f"Log File       : {log_file.name}")
print(f"Batch Size     : {BATCH_SIZE}")
print("=" * 50)

logging.info("Configuration completed successfully.")


# 2. Read & Validate UniProt IDs

# Function to Validate UniProt Accession Numbers
def validate_uniprot_id(accession):

    pattern = r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-Z0-9]{10}$"

    return bool(re.match(pattern, accession))

# Function to Read UniProt IDs from a Text File
def read_uniprot_ids(input_file):

    logging.info("Reading UniProt accession IDs...")

    valid_ids = []
    invalid_ids = []

    with open(input_file, "r") as file:

        for line in file:
            accession = line.strip()

            if accession == "":
                continue

            if validate_uniprot_id(accession):
                valid_ids.append(accession)

            else:
                invalid_ids.append(accession)

    # Remove duplicates
    total_read = len(valid_ids)
    valid_ids = sorted(set(valid_ids))
    duplicates_removed = total_read - len(valid_ids)

    logging.info(f"Total IDs Read          : {total_read}")
    logging.info(f"Valid IDs              : {len(valid_ids)}")
    logging.info(f"Duplicate IDs Removed  : {duplicates_removed}")
    logging.info(f"Invalid IDs            : {len(invalid_ids)}")

    print("\nModule 2 Summary")
    print("-" * 50)
    print(f"Total IDs Read          : {total_read}")
    print(f"Valid IDs              : {len(valid_ids)}")
    print(f"Duplicate IDs Removed  : {duplicates_removed}")
    print(f"Invalid IDs            : {len(invalid_ids)}")
    print("-" * 50)

    return valid_ids, invalid_ids

# 3. Retrieve UniProt Metadata

POLLING_INTERVAL = 3

# Function to check HTTP Response
def check_response(response):

    try:
        response.raise_for_status()

    except requests.HTTPError:

        try:
            print(response.json())

        except Exception:
            print(response.text)

        raise

# Function to submit ID Mapping Job
def submit_mapping_job(batch):

    print(f"Submitting batch containing {len(batch)} IDs...")

    response = session.post(
        "https://rest.uniprot.org/idmapping/run",
        data={
            "from": "UniProtKB_AC-ID",
            "to": "UniProtKB",
            "ids": ",".join(batch)
        },
        timeout=TIMEOUT
    )

    check_response(response)
    job_id = response.json()["jobId"]
    logging.info(f"Job submitted successfully: {job_id}")
    
    print(f"Job ID : {job_id}")

    return job_id

# Function to poll the UniProt ID Mapping API until the submitted job has completed successfully
def wait_until_ready(job_id):

    print("\nWaiting for UniProt to complete the job...")

    while True:
        response = session.get(f"https://rest.uniprot.org/idmapping/status/{job_id}", timeout=TIMEOUT)

        check_response(response)
        status = response.json()

        # Job Finished
        if "results" in status or "failedIds" in status:
            print("✓ Job completed successfully.")
            logging.info(f"Job {job_id} completed successfully.")

            return {
                "job_id": job_id,
                "status": "COMPLETED"
                }

        # Job Still Running
        elif status.get("jobStatus") in ["NEW", "RUNNING"]:

            print(f"Current Status : {status['jobStatus']} ... waiting...")
            time.sleep(POLLING_INTERVAL)

        # Unexpected Status
        else:
            raise RuntimeError(f"Unexpected UniProt job status: {status}")

# Function to Retrieve the Redirect URL for a Completed Job 
def get_results_url(job_id): 
    
    response = session.get(f"https://rest.uniprot.org/idmapping/details/{job_id}", timeout=TIMEOUT) 
    check_response(response) 
    redirect_url = response.json()["redirectURL"] 
    
    print("\nRedirect URL:") 
    print(redirect_url) 
    
    logging.info(f"Redirect URL: {redirect_url}") 
    
    return redirect_url

# Function to Get Next Page URL from HTTP Link Header
def get_next_link(headers):

    link_header = headers.get("Link")

    if not link_header:
        return None

    match = re.search(r'<(.+?)>;\s*rel="next"', link_header)

    if match:
        return match.group(1)

    return None

# Function to Download One Page of Results
def fetch_results(url):

    response = session.get(url, timeout=TIMEOUT)
    check_response(response)

    return response.json(), response.headers

# Function to Retrieve Complete Metadata from UniProt
def retrieve_metadata_json(job_id):

    print("\nRetrieving UniProt metadata...")

    redirect_url = get_results_url(job_id)
    url = redirect_url + "?format=json"

    all_results = []
    page = 1

    while url:

        print(f"Downloading page {page}...")

        data, headers = fetch_results(url)
        results = data.get("results", [])

        print(f"Records on page {page}: {len(results)}")

        logging.info(f"Page {page}: {len(results)} records.")

        all_results.extend(results)
        url = get_next_link(headers)
        page += 1

    print(f"\nTotal Records Retrieved : {len(all_results)}")

    logging.info(f"Retrieved {len(all_results)} protein records.")
    
    return all_results

# Function to Validate Retrieved Records
def validate_batch(expected_ids, results):

    expected = len(expected_ids)
    observed = len(results)

    print("\nBatch Validation")
    print("-" * 50)
    print(f"Expected : {expected}")
    print(f"Retrieved : {observed}")

    if expected == observed:
        print("Status : OK")

        logging.info("Batch validation passed.")

    else:
        missing = expected - observed

        print(f"Missing : {missing}")

        logging.warning(f"Expected {expected}, Retrieved {observed}")

# 4. 

# Function to Parse a Single UniProt JSON Record
def parse_protein_record(record):

    protein = record.get("to", {})

    # UniProt Accession
    accession = protein.get("primaryAccession", "")

    # Entry Name
    entry_name = protein.get("uniProtkbId", "")

    # Organism
    organism = protein.get("organism", {})
    organism_name = organism.get("scientificName", "")
    taxonomy_id = organism.get("taxonId", "")
    lineage = "; ".join(organism.get("lineage", []))

    # Protein Length
    sequence_length = (protein.get("sequence", {}).get("length", ""))

    # Reviewed / Unreviewed
    entry_type = protein.get("entryType", "")
    reviewed = ("Reviewed" if "reviewed" in entry_type.lower() else "Unreviewed")

    # Protein Existence
    protein_existence = protein.get("proteinExistence","")

    # Protein Name
    protein_name = ""
    description = protein.get("proteinDescription",{})

    if "recommendedName" in description:
        protein_name = (description["recommendedName"].get("fullName", {}).get("value", ""))

    elif "submissionNames" in description:
        names = description["submissionNames"]

        if len(names):
            protein_name = (names[0].get("fullName", {}).get("value", ""))

    # Gene Name
    gene_name = ""
    genes = protein.get("genes", [])

    if genes:
        gene = genes[0]

        if "geneName" in gene:
            gene_name = (gene["geneName"].get("value", ""))

        elif "orfNames" in gene:
            orfs = gene["orfNames"]

            if len(orfs):
                gene_name = orfs[0].get("value","")

    # Function Annotation
    function = ""

    comments = protein.get("comments", [])

    for comment in comments:
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])

            if texts:
                function = texts[0].get("value","")

                break

    return {
        "UniProt_ID": accession,
        "Entry_Name": entry_name,
        "Protein_Name": protein_name,
        "Gene_Name": gene_name,
        "Organism": organism_name,
        "Taxonomy_ID": taxonomy_id,
        "Taxonomic_Lineage": lineage,
        "Reviewed": reviewed,
        "Sequence_Length": sequence_length,
        "Protein_Existence": protein_existence,
        "Function": function
    }

# Function to Parse All UniProt Records
def parse_metadata(metadata_json):

    print("\nParsing metadata...")

    parsed_records = []

    for record in metadata_json:
        parsed_records.append(parse_protein_record(record))

    metadata_df = pd.DataFrame(parsed_records)

    print(f"Parsed {len(metadata_df)} protein records.")

    logging.info(f"Parsed {len(metadata_df)} protein records.")

    return metadata_df

# 5.

# Function to Process All UniProt IDs in Batches
def process_all_batches(valid_ids):
    
    pipeline_start = time.time()

    print("\nStarting metadata retrieval for all UniProt IDs...")

    total_ids = len(valid_ids)
    total_batches = (total_ids + BATCH_SIZE - 1) // BATCH_SIZE

    print(f"Total IDs      : {total_ids}")
    print(f"Batch Size     : {BATCH_SIZE}")
    print(f"Total Batches  : {total_batches}")

    logging.info(f"Processing {total_ids} UniProt IDs in {total_batches} batches.")

    all_dataframes = []
    failed_ids = []
    total_records = 0

    for batch_number in range(total_batches):
        start = batch_number * BATCH_SIZE
        end = min(start + BATCH_SIZE, total_ids)
        batch = valid_ids[start:end]

        print("\n" + "=" * 50)
        print(f"Batch {batch_number + 1} of {total_batches}")
        print("=" * 50)
        print(f"First ID : {batch[0]}")
        print(f"Last ID  : {batch[-1]}")
        
        batch_start_time = time.time()

        try:
            job_id = submit_mapping_job(batch)
            job_info = wait_until_ready(job_id)
            
            # Retry retrieval if UniProt returns an empty result
            MAX_DOWNLOAD_RETRIES = 3
            metadata_json = []

            for attempt in range(1, MAX_DOWNLOAD_RETRIES + 1):
                metadata_json = retrieve_metadata_json(job_info["job_id"])

                if len(metadata_json) > 0:
                    print(f"Metadata successfully retrieved on attempt {attempt}.")
                    break

                print(f"No metadata returned. Retrying ({attempt}/{MAX_DOWNLOAD_RETRIES})...")

                logging.warning(f"Batch {batch_number + 1}: Empty metadata. Retry {attempt}.")
                
                time.sleep(10)

            validate_batch(batch,metadata_json)
            
            batch_df = parse_metadata(metadata_json)
            all_dataframes.append(batch_df)
            
            total_records += len(batch_df)
            
            print(f"Running Total : {total_records} proteins")

            logging.info(f"Running Total : {total_records}")
            logging.info(f"Batch {batch_number + 1} completed successfully.")

        except Exception as error:

            logging.error(f"Batch {batch_number + 1} failed : {error}")
            failed_ids.extend(batch)

            print(f"Batch {batch_number + 1} failed.")
        
        batch_end_time = time.time()

        elapsed = batch_end_time - batch_start_time

        print(f"Batch Execution Time : {elapsed:.2f} seconds")

        logging.info(f"Batch Execution Time : {elapsed:.2f} seconds")

    pipeline_end = time.time()
    pipeline_time = pipeline_end - pipeline_start

    print("\n" + "=" * 50)
    print("Pipeline Summary")
    print("=" * 50)
    print(f"Total Records : {total_records}")
    print(f"Failed IDs    : {len(failed_ids)}")
    print(f"Total Time    : {pipeline_time/60:.2f} minutes")

    logging.info(f"Pipeline completed in {pipeline_time/60:.2f} minutes.")
    
    return all_dataframes, failed_ids

# Function to Merge All Batch DataFrames
def merge_dataframes(all_dataframes):

    print("\nMerging all batch DataFrames...")

    if len(all_dataframes) == 0:
        raise ValueError("No DataFrames available for merging.")

    master_df = pd.concat(
        all_dataframes,
        ignore_index=True
    )

    expected_rows = sum(len(df) for df in all_dataframes)
    actual_rows = len(master_df)

    if expected_rows != actual_rows:
        raise ValueError(f"Merge Error: Expected {expected_rows} rows but found {actual_rows}.")

    logging.info("Merge validation passed.")

    print(f"Master DataFrame contains {actual_rows} records.")

    logging.info(f"Merged {len(all_dataframes)} DataFrames.")
    logging.info(f"Master DataFrame contains {actual_rows} records.")

    return master_df

# Function to Prepare Master DataFrame
def prepare_master_dataframe(master_df):

    print("\nPreparing Master DataFrame...")
    print("-" * 50)

    # Sort by UniProt ID
    master_df.sort_values(by="UniProt_ID", inplace=True)

    # Reset Index
    master_df.reset_index(drop=True, inplace=True)

    # Column Order
    master_df = master_df[
        [
            "UniProt_ID",
            "Entry_Name",
            "Protein_Name",
            "Gene_Name",
            "Organism",
            "Taxonomy_ID",
            "Taxonomic_Lineage",
            "Reviewed",
            "Sequence_Length",
            "Protein_Existence",
            "Function"
        ]
    ]

    # Duplicate Check
    duplicate_count = master_df["UniProt_ID"].duplicated().sum()

    if duplicate_count > 0:
        raise ValueError(f"Duplicate UniProt IDs detected: {duplicate_count}")

    # Missing ID Check
    missing_count = master_df["UniProt_ID"].isna().sum()

    if missing_count > 0:
        raise ValueError(f"Missing UniProt IDs detected: {missing_count}")

    # Final Summary
    print("Master DataFrame successfully prepared.\n")
    print(f"Rows                : {len(master_df)}")
    print(f"Columns             : {len(master_df.columns)}")
    print(f"Duplicate IDs       : {duplicate_count}")
    print(f"Missing IDs         : {missing_count}")

    logging.info("Master DataFrame prepared successfully.")
    logging.info(f"Rows: {len(master_df)}")
    logging.info(f"Columns: {len(master_df.columns)}")
    logging.info(f"Duplicate IDs: {duplicate_count}")
    logging.info(f"Missing IDs: {missing_count}")

    return master_df

# Function to Export Metadata CSV
def export_metadata_csv(master_df):

    print("\nExporting metadata CSV...")
    print("-" * 50)

    master_df.to_csv(output_file, index=False)
    
    print(f"Records Exported : {len(master_df)}")
    
    logging.info(f"Exported {len(master_df)} metadata records.")
    
    if output_file.exists():
        print("✓ CSV export verified.")

        logging.info("CSV export verified.")

    else:
        raise FileNotFoundError(f"Failed to create:\n{output_file}")
    
    print(f"Metadata successfully exported to:\n{output_file}")

    logging.info(f"Metadata CSV exported: {output_file}")

    return

# Function to Export Failed UniProt IDs
def export_failed_ids(failed_ids):

    print("\nExporting failed UniProt IDs...")
    print("-" * 60)

    if len(failed_ids) == 0:
        print("No failed UniProt IDs.")

        logging.info("No failed UniProt IDs.")

        return

    with open(failed, "w") as file:
        for accession in failed_ids:
            file.write(f"{accession}\n")

    print(f"{len(failed_ids)} failed IDs exported to:\n{failed}")

    logging.info(f"Failed IDs exported: {failed}")

# Function to Generate Experiment Report
def generate_experiment_report(valid_ids,
    master_df,
    failed_ids,
    total_batches,
    total_time,
    input_file,
    output_file,
    report_file,
    log_file):

    print("\nGenerating Experiment Report...")
    print("-" * 50)

    report_lines = []

    # Report Header
    report_lines.append("=" * 50)
    report_lines.append("Experiment : EXP002B")
    report_lines.append("Title      : UniProt Metadata Retrieval")
    report_lines.append("=" * 50)
    report_lines.append("")
    report_lines.append(f"Date : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Input Information
    report_lines.append("")
    report_lines.append("INPUT")
    report_lines.append("-" * 50)
    report_lines.append(f"Input File            : {input_file}")
    report_lines.append(f"Total UniProt IDs     : {len(valid_ids)}")

    # Processing Information
    report_lines.append("")
    report_lines.append("PROCESSING")
    report_lines.append("-" * 50)
    report_lines.append(f"Batch Size            : {BATCH_SIZE}")
    report_lines.append(f"Total Batches         : {total_batches}")
    report_lines.append(f"Records Exported      : {len(master_df)}")
    report_lines.append(f"Metadata Columns      : {len(master_df.columns)}")
    report_lines.append(f"Failed IDs            : {len(failed_ids)}")

    # Output Files
    report_lines.append("")
    report_lines.append("OUTPUT")
    report_lines.append("-" * 50)
    report_lines.append(f"Metadata CSV          : {output_file}")
    report_lines.append(f"Failed IDs File       : {failed}")
    report_lines.append(f"Log File              : {log_file}")

    # Pipeline Summary
    report_lines.append("")
    report_lines.append("PIPELINE")
    report_lines.append("-" * 50)
    report_lines.append(f"Execution Time        : {round(total_time, 2)} minutes ({round(total_time * 60, 2)} seconds)")

    if len(failed_ids) == 0:
        status = "SUCCESS"

    elif len(master_df) > 0:
        status = "PARTIAL SUCCESS"

    else:
        status = "FAILED"

    report_lines.append(f"Overall Status        : {status}")

    # Report Footer
    report_lines.append("")
    report_lines.append("=" * 50)
    report_lines.append("End of Report")
    report_lines.append("=" * 50)

    # Write Report to Disk
    with open(report_file, "w", encoding="utf-8") as file:

        for line in report_lines:
            file.write(line + "\n")

    print("Experiment report generated successfully.")
    print(report_file)

    logging.info("Experiment report generated successfully.")

# Function to Check Pipeline Integrity
def pipeline_integrity_check(valid_ids, master_df, failed_ids):

    print("\nPerforming Pipeline Integrity Check...")
    print("-" * 50)

    input_ids = len(valid_ids)
    exported_records = len(master_df)
    failed_records = len(failed_ids)
    duplicate_records = master_df["UniProt_ID"].duplicated().sum()
    exported_ids = set(master_df["UniProt_ID"])
    input_id_set = set(valid_ids)
    missing_ids = sorted(input_id_set - exported_ids)

    print(f"Input UniProt IDs        : {input_ids}")
    print(f"Exported Records         : {exported_records}")
    print(f"Duplicate Records        : {duplicate_records}")
    print(f"Failed IDs              : {failed_records}")
    print(f"Missing UniProt IDs      : {len(missing_ids)}")

    logging.info("=" * 50)
    logging.info("PIPELINE INTEGRITY CHECK")
    logging.info("=" * 50)
    logging.info(f"Input IDs        : {input_ids}")
    logging.info(f"Exported Records : {exported_records}")
    logging.info(f"Duplicate Records: {duplicate_records}")
    logging.info(f"Failed IDs       : {failed_records}")
    logging.info(f"Missing IDs      : {len(missing_ids)}")

    integrity_pass = (
        duplicate_records == 0 and
        failed_records == 0
    )

    if integrity_pass:
        print("\nIntegrity Status : PASSED")

        if len(missing_ids) > 0:
            print(f"Warning : {len(missing_ids)} input IDs did not produce metadata.")

            logging.warning(f"{len(missing_ids)} IDs missing from exported metadata.")

        else:
            print("All input IDs accounted for.")

            logging.info("All input IDs accounted for.")

    else:
        print("\nIntegrity Status : FAILED")

        logging.error("Pipeline integrity check failed.")

    return missing_ids

# Main Function
def main():

    try:
        valid_ids, invalid_ids = read_uniprot_ids(input_file)

        pipeline_start = time.time()
        all_dataframes, failed_ids = process_all_batches(valid_ids)
        pipeline_end = time.time()
        total_time = (pipeline_end - pipeline_start) / 60
        total_batches = (len(valid_ids) + BATCH_SIZE - 1) // BATCH_SIZE

        master_df = merge_dataframes(all_dataframes)

        master_df = prepare_master_dataframe(master_df)

        export_metadata_csv(master_df)
        export_failed_ids(failed_ids)

        generate_experiment_report(
            valid_ids=valid_ids,
            master_df=master_df,
            failed_ids=failed_ids,
            total_batches=total_batches,
            total_time=total_time,
            input_file=input_file,
            output_file=output_file,
            report_file=report_file,
            log_file=log_file
        )
        
        missing_ids = pipeline_integrity_check(
            valid_ids,
            master_df,
            failed_ids
        )
        
        # Export Missing IDs (if any)
        if len(missing_ids) > 0:
            missing_file = metadata / "EXP002B_Missing_Metadata_IDs.txt"

            with open(missing_file, "w") as f:
                for accession in missing_ids:
                    f.write(accession + "\n")
    
            print("\nMissing IDs exported to:")
            print(missing_file)

            logging.info(f"Missing metadata IDs exported to {missing_file}")

        print("\n")
        print("=" * 50)
        print("EXP002B COMPLETED SUCCESSFULLY")
        print("=" * 50)

        logging.info("EXP002B completed successfully.")

    except Exception as error:

        print("\n")
        print("=" * 50)
        print("PIPELINE FAILED")
        print("=" * 50)
        print(error)

        logging.exception(error)

        raise

# Execute the Main Function
if __name__ == "__main__":
    main()