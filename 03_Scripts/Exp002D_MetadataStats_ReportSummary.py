# Experiment 002D

# Metadata Statistics & Report Summary

# Import Required Libraries
import os
import time
import logging
import datetime
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# Configuration
PROJECT_DIR = r"C:\Users\Anirban Majumder\OneDrive - RICE Group\Desktop\Academics\New folder\PhD_Research_Work"
INPUT_DIR = os.path.join(PROJECT_DIR, "05_Metadata")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "05_Metadata")
RESULT_DIR = os.path.join(PROJECT_DIR, "04_Results")
FIGURE_DIR = os.path.join(PROJECT_DIR, "07_Figures")
TABLE_DIR = os.path.join(PROJECT_DIR, "08_Tables")
LOG_DIR = os.path.join(PROJECT_DIR, "06_Logs")

# Files
input_file = os.path.join(INPUT_DIR, "EXP002C_Clean_Metadata.csv")
output_file = os.path.join(OUTPUT_DIR, "EXP002D_Metadata_Summary.csv")
statistics_report = os.path.join(RESULT_DIR, "EXP002D_Metadata_Statistics.txt")
experiment_report = os.path.join(RESULT_DIR, "EXP002D_Metadata_Report.txt")
log_file = os.path.join(LOG_DIR, "EXP002D_Log.txt")

# Table Files
TABLE_FILES = {
    "organisms":
        os.path.join(TABLE_DIR, "Organism_Frequency.csv"),
    "protein_names":
        os.path.join(TABLE_DIR, "Protein_Name_Frequency.csv"),
    "protein_existence":
        os.path.join(TABLE_DIR, "Protein_Existence_Counts.csv"),
    "missing_values":
        os.path.join(TABLE_DIR, "Missing_Value_Summary.csv"),
    "sequence_statistics":
        os.path.join(TABLE_DIR, "Sequence_Length_Statistics.csv")
}

# Figure Files
FIGURE_FILES = {
    "histogram":
        os.path.join(FIGURE_DIR, "Sequence_Length_Histogram.png"),
    "boxplot":
        os.path.join(FIGURE_DIR,"Sequence_Length_Boxplot.png"),
    "protein_existence":
        os.path.join(FIGURE_DIR, "Protein_Existence_Pie.png"),
    "organisms":
        os.path.join(FIGURE_DIR, "Top20_Organisms.png"),
    "protein_names":
        os.path.join(FIGURE_DIR, "Top20_Protein_Names.png")
}

# Logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    filemode="w"
)

logging.info("=" * 50)
logging.info("EXP002D Started")
logging.info("=" * 50)
logging.info("Configuration completed successfully.")

print("=" * 50)
print("Experiment 002D")
print("Metadata Statistics & Summary")
print("=" * 50)
print(f"Input File  : {input_file}")
print(f"Output File : {output_file}")
print(f"Log File    : {log_file}")

# Experiment Summary
print("\nExperiment Summary")
print("-" * 50)
print("Input Dataset      : EXP002C_Clean_Metadata.csv")
print("Analysis Modules   :")
print("  • Dataset Statistics")
print("  • Biological Statistics")
print("  • Sequence Statistics")
print("  • Frequency Tables")
print("  • Visualizations")
print("  • Reports")

print("\nOutput Folders")
print(f"  • {RESULT_DIR}")
print(f"  • {OUTPUT_DIR}")
print(f"  • {FIGURE_DIR}")
print(f"  • {TABLE_DIR}")

logging.info("Experiment summary displayed.")
logging.info(f"Input File : {input_file}")
logging.info(f"Output Directory : {OUTPUT_DIR}")
logging.info(f"Figure Directory : {FIGURE_DIR}")
logging.info(f"Table Directory : {TABLE_DIR}")

# Load & Verify Metadata
def load_metadata():

    print("\nLoading Clean Metadata...")
    print("-" * 50)

    logging.info("Checking input metadata file...")

    # Check whether input file exists
    if not os.path.exists(input_file):
        logging.error("Input metadata file not found.")
        raise FileNotFoundError(f"Input metadata file not found:\n{input_file}")

    logging.info("Input metadata file found.")

    # Load metadata
    metadata_df = pd.read_csv(input_file)

    print(f"Rows          : {len(metadata_df)}")
    print(f"Columns       : {metadata_df.shape[1]}")

    memory_mb = metadata_df.memory_usage(deep=True).sum() / (1024 ** 2)

    print(f"Memory Usage  : {memory_mb:.2f} MB")

    logging.info(f"Rows: {len(metadata_df)}")
    logging.info(f"Columns: {metadata_df.shape[1]}")
    logging.info(f"Memory Usage: {memory_mb:.2f} MB")

    # Verify required columns
    required_columns = [
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

    missing_columns = [
        column
        for column in required_columns
        if column not in metadata_df.columns
    ]

    if len(missing_columns) > 0:
        logging.error(f"Missing columns detected: {missing_columns}")

        raise ValueError(
            "Required metadata columns are missing:\n"
            + "\n".join(missing_columns)
        )

    print("\nRequired column verification : PASSED")

    logging.info("Required column verification passed.")

    # Display column names
    print("\nMetadata Columns")
    print("-" * 50)

    for column in metadata_df.columns:
        print(column)

    logging.info(
        "Metadata loaded successfully with %d records and %d columns.",
        len(metadata_df),
        metadata_df.shape[1]
    )

    return metadata_df

# Dataset Statistics
def generate_dataset_statistics(metadata_df):

    print("\nGenerating Dataset Statistics...")
    print("-" * 50)

    logging.info("Generating dataset statistics...")

    statistics = {
        "Dataset Summary": {
            "Total Records":
                len(metadata_df),
            "Total Columns":
                metadata_df.shape[1],
            "Memory Usage (MB)":
                round(
                    metadata_df.memory_usage(
                        deep=True
                    ).sum() / (1024 ** 2),
                    2
                )
        },
        "Uniqueness": {
            "Duplicate UniProt IDs":
                metadata_df["UniProt_ID"].duplicated().sum(),
            "Unique Organisms":
                metadata_df["Organism"].nunique(),
            "Unique Taxonomy IDs":
                metadata_df["Taxonomy_ID"].nunique()
        },
        "Review Status": {
            "Reviewed":
                metadata_df["Reviewed"]
                .astype(str)
                .str.lower()
                .eq("reviewed")
                .sum(),
            "Unreviewed":
                metadata_df["Reviewed"]
                .astype(str)
                .str.lower()
                .ne("reviewed")
                .sum()
        },
        "Missing Values": {
            column: {
                "Count":
                    int(metadata_df[column].isna().sum()),
                "Percentage":
                    round(
                        metadata_df[column].isna().mean() * 100,
                        2
                    )
            }
        for column in metadata_df.columns
        }
    }

    print("Dataset statistics generated successfully.")

    logging.info("Dataset statistics generated successfully.")

    return statistics

# Biological Statistics
def generate_biological_statistics(metadata_df):

    print("\nGenerating Biological Statistics...")
    print("-" * 50)

    logging.info("Generating biological statistics...")

    biological_statistics = {}

    # Protein Statistics
    protein_stats = {
        "Unique Protein Names":
            metadata_df["Protein_Name"].nunique(),
        "Top 20 Most Frequent Protein Names":
            metadata_df["Protein_Name"]
            .value_counts()
            .head(20)
            .to_dict()
    }

    # Gene Statistics
    gene_stats = {
        "Unique Gene Names":
            metadata_df["Gene_Name"].nunique(),
        "Missing Gene Names":
            int(metadata_df["Gene_Name"].isna().sum()),
        "Top 20 Most Frequent Gene Names":
            metadata_df["Gene_Name"]
            .dropna()
            .value_counts()
            .head(20)
            .to_dict()
    }

    # Protein Existence Statistics
    existence_stats = (
        metadata_df["Protein_Existence"]
        .value_counts(dropna=False)
    )

    protein_existence = {}
    total_records = len(metadata_df)

    for label, count in existence_stats.items():
        protein_existence[label] = {
            "Count": int(count),
            "Percentage": round(
                count / total_records * 100,
                2
            )
        }

    # Organism Statistics
    organism_counts = metadata_df["Organism"].value_counts()
    organism_stats = {
        "Unique Organisms":
            metadata_df["Organism"].nunique(),
        "Average Proteins per Organism":
            round(
                len(metadata_df)
                /
                metadata_df["Organism"].nunique(),
                2
            ),
        "Median Proteins per Organism":
            int(organism_counts.median()),
        "Maximum Proteins in One Organism":
            int(organism_counts.max()),
        "Top 20 Most Represented Organisms":
            organism_counts
            .head(20)
            .to_dict(),
        "Top 20 Most Frequent Taxonomy IDs":
            metadata_df["Taxonomy_ID"]
            .value_counts()
            .head(20)
            .to_dict()
    }

    biological_statistics["Protein Statistics"] = protein_stats
    biological_statistics["Gene Statistics"] = gene_stats
    biological_statistics["Protein Existence"] = protein_existence
    biological_statistics["Organism Statistics"] = organism_stats

    print("Biological statistics generated successfully.")

    logging.info("Biological statistics generated successfully.")

    return biological_statistics

# Sequence Statistics
def generate_sequence_statistics(metadata_df):

    print("\nGenerating Sequence Statistics...")
    print("-" * 50)

    logging.info("Generating sequence statistics...")

    sequence = metadata_df["Sequence_Length"].dropna()
    
    sequence_statistics = {
        "Total Sequences":
            len(sequence),
        "Minimum Length":
            int(sequence.min()),
        "Maximum Length":
            int(sequence.max()),
        "Range":
            int(sequence.max() - sequence.min()),
        "Mean Length":
            round(sequence.mean(), 2),
        "Median Length":
            round(sequence.median(), 2),
        "Variance":
            round(sequence.var(), 2),
        "Standard Deviation":
            round(sequence.std(), 2),
        "Q1":
            round(sequence.quantile(0.25), 2),
        "Q2":
            round(sequence.quantile(0.50), 2),
        "Q3":
            round(sequence.quantile(0.75), 2),
        "IQR":
            round(sequence.quantile(0.75) - sequence.quantile(0.25), 2)
    }

    print("Sequence statistics generated successfully.")

    logging.info("Sequence statistics generated successfully.")

    return sequence_statistics

# Frequency Tables
def generate_frequency_tables(metadata_df):

    print("\nGenerating Frequency Tables...")
    print("-" * 50)

    logging.info("Generating frequency tables...")

    frequency_tables = {}

    frequency_tables["Protein_Names"] = (
        metadata_df["Protein_Name"]
        .value_counts()
        .head(20)
        .reset_index(name="Count")
        .rename(columns={"Protein_Name":"Protein_Name"})
    )

    frequency_tables["Gene_Names"] = (
        metadata_df["Gene_Name"]
        .dropna()
        .value_counts()
        .head(20)
        .reset_index(name="Count")
        .rename(columns={"Gene_Name":"Gene_Name"})
    )

    frequency_tables["Organisms"] = (
        metadata_df["Organism"]
        .value_counts()
        .head(20)
        .reset_index(name="Count")
        .rename(columns={"Organism":"Organism"})
    )

    frequency_tables["Taxonomy_IDs"] = (
        metadata_df["Taxonomy_ID"]
        .value_counts()
        .head(20)
        .reset_index(name="Count")
        .rename(columns={"Taxonomy_ID":"Taxonomy_ID"})
    )

    frequency_tables["Protein_Existence"] = (
        metadata_df["Protein_Existence"]
        .value_counts()
        .reset_index(name="Count")
        .rename(columns={"Protein_Existence":"Protein_Existence"})
    )

    missing = metadata_df.isna().sum()

    frequency_tables["Missing_Values"] = pd.DataFrame({
        "Column": missing.index,
        "Missing_Count": missing.values,
        "Missing_Percentage":
            (missing.values / len(metadata_df) * 100).round(2)
    })

    print("Frequency tables generated successfully.")

    logging.info("Frequency tables generated successfully.")

    return frequency_tables

# Data Visualization

PLOT_TEMPLATE = "plotly_white"

TITLE_FONT = dict(family="Arial", size=22)
AXIS_FONT = dict(family="Arial", size=16)

# Figure 1: Sequence Length Histogram
def create_sequence_length_histogram(metadata_df):

    print("\nCreating Sequence Length Histogram...")

    logging.info("Creating sequence length histogram.")

    fig = px.histogram(
        metadata_df,
        x="Sequence_Length",
        nbins=50,
        title="Distribution of Protein Sequence Lengths",
        template=PLOT_TEMPLATE
    )

    fig.update_layout(
        title_font=TITLE_FONT,
        xaxis_title="Sequence Length (Amino Acids)",
        yaxis_title="Number of Proteins",
        font=AXIS_FONT
    )

    output = os.path.join(FIGURE_DIR, "Figure_01_Sequence_Length_Histogram.png")

    fig.write_image(
        output,
        width=2200,
        height=1600,
        scale=4
    )

    logging.info("Histogram saved.")

    print("✓ Histogram created.")

    return fig

# Sequence Length Box Plot
def create_sequence_length_boxplot(metadata_df):

    print("Creating Sequence Length Boxplot...")

    logging.info("Creating sequence length boxplot.")

    fig = px.box(
        metadata_df,
        y="Sequence_Length",
        points="outliers",
        title="Protein Sequence Length Distribution",
        template=PLOT_TEMPLATE
    )

    fig.update_layout(
        title_font=TITLE_FONT,
        yaxis_title="Sequence Length (Amino Acids)",
        font=AXIS_FONT
    )

    output = os.path.join(FIGURE_DIR, "Figure_02_Sequence_Length_Boxplot.png")

    fig.write_image(
        output,
        width=1600,
        height=1200,
        scale=3
    )

    logging.info("Boxplot saved.")

    print("✓ Boxplot created.")

    return fig

# Figure 3: Protein Existence Pie Chart
def create_protein_existence_piechart(metadata_df):

    print("Creating Protein Existence Pie Chart...")

    logging.info("Creating protein existence pie chart.")

    counts = (
        metadata_df["Protein_Existence"]
        .value_counts()
        .reset_index()
    )

    counts.columns = [
        "Protein_Existence",
        "Count"
    ]

    fig = px.pie(
        counts,
        names="Protein_Existence",
        values="Count",
        hole=0.40,
        title="Protein Existence Evidence",
        template=PLOT_TEMPLATE
    )

    fig.update_traces(
        textposition="inside",
        textinfo="percent+label"
    )

    fig.update_layout(
        title_font=TITLE_FONT,
        font=AXIS_FONT
    )

    output = os.path.join(FIGURE_DIR, "Figure_03_Protein_Existence_PieChart.png")

    fig.write_image(
        output,
        width=1600,
        height=1200,
        scale=3
    )

    logging.info("Protein existence pie chart saved.")

    print("✓ Pie chart created.")

    return fig

# Figure 4: Top 20 Protein Names
def create_top20_protein_bar(frequency_tables):

    print("Creating Top 20 Protein Names Bar Chart...")

    logging.info("Creating Top 20 Protein Names Bar Chart.")

    protein_df = frequency_tables["Protein_Names"]

    fig = px.bar(
        protein_df,
        x="Count",
        y="Protein_Name",
        orientation="h",
        title="Top 20 Most Frequent Protein Names",
        template=PLOT_TEMPLATE
    )

    fig.update_layout(
        title_font=TITLE_FONT,
        xaxis_title="Protein Count",
        yaxis_title="Protein Name",
        font=AXIS_FONT,
        yaxis=dict(categoryorder="total ascending")
    )

    output = os.path.join(FIGURE_DIR, "Figure_04_Top20_Protein_Names.png")

    fig.write_image(
        output,
        width=2200,
        height=1600,
        scale=3
    )

    print("✓ Protein Name chart created.")

    logging.info("Protein Name chart saved.")

    return fig

# Figure 5: Top 20 Gene Names
def create_top20_gene_bar(frequency_tables):

    print("Creating Top 20 Gene Names Bar Chart...")

    logging.info("Creating Top 20 Gene Names Bar Chart.")

    gene_df = frequency_tables["Gene_Names"]

    fig = px.bar(
        gene_df,
        x="Count",
        y="Gene_Name",
        orientation="h",
        title="Top 20 Most Frequent Gene Names",
        template=PLOT_TEMPLATE
    )

    fig.update_layout(
        title_font=TITLE_FONT,
        xaxis_title="Frequency",
        yaxis_title="Gene Name",
        font=AXIS_FONT,
        yaxis=dict(categoryorder="total ascending")
    )

    output = os.path.join(FIGURE_DIR, "Figure_05_Top20_Gene_Names.png")

    fig.write_image(
        output,
        width=1800,
        height=1400,
        scale=3
    )

    print("✓ Gene Name chart created.")

    logging.info("Gene Name chart saved.")

    return fig

# Figure 6: Top 20 Organisms
def create_top20_organism_bar(frequency_tables):

    print("Creating Top 20 Organisms Bar Chart...")

    logging.info("Creating Top 20 Organisms Bar Chart.")

    organism_df = frequency_tables["Organisms"]

    fig = px.bar(
        organism_df,
        x="Count",
        y="Organism",
        orientation="h",
        title="Top 20 Most Represented Organisms",
        template=PLOT_TEMPLATE
    )

    fig.update_layout(
        title_font=TITLE_FONT,
        xaxis_title="Protein Count",
        yaxis_title="Organism",
        font=AXIS_FONT,
        yaxis=dict(categoryorder="total ascending")
    )

    output = os.path.join(FIGURE_DIR, "Figure_06_Top20_Organisms.png")

    fig.write_image(
        output,
        width=2200,
        height=1600,
        scale=3
    )

    print("✓ Organism chart created.")

    logging.info("Organism chart saved.")

    return fig

# Figure 7: Missing Values
def create_missing_values_bar(frequency_tables):

    print("Creating Missing Values Chart...")

    logging.info("Creating Missing Values Chart.")

    missing_df = frequency_tables["Missing_Values"]

    fig = px.bar(
        missing_df,
        x="Missing_Count",
        y="Column",
        orientation="h",
        title="Missing Values Across Metadata Fields",
        template=PLOT_TEMPLATE
    )

    fig.update_layout(
        title_font=TITLE_FONT,
        xaxis_title="Missing Records",
        yaxis_title="Metadata Column",
        font=AXIS_FONT,
        yaxis=dict(categoryorder="total ascending")
    )

    output = os.path.join(FIGURE_DIR, "Figure_07_Missing_Values.png")

    fig.write_image(
        output,
        width=1800,
        height=1400,
        scale=3
    )

    print("✓ Missing Value chart created.")

    logging.info("Missing Value chart saved.")

    return fig

# Generate All Figures
def generate_all_visualizations(metadata_df, frequency_tables):

    print("\nGenerating Plotly Visualizations...")
    print("-" * 50)

    logging.info("Generating all Plotly visualizations.")

    create_sequence_length_histogram(metadata_df)
    create_sequence_length_boxplot(metadata_df)
    create_protein_existence_piechart(metadata_df)
    create_top20_protein_bar(frequency_tables)
    create_top20_gene_bar(frequency_tables)
    create_top20_organism_bar(frequency_tables)
    create_missing_values_bar(frequency_tables)

    print("\nAll visualizations generated successfully.")

    logging.info("All visualizations generated successfully.")

# Export Results
def export_results(metadata_df, dataset_statistics, biological_statistics, sequence_statistics, frequency_tables):

    print("\nExporting Results...")
    print("-" * 50)

    logging.info("Exporting experiment results...")

    # Create output directories if they do not exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(TABLE_DIR, exist_ok=True)
    os.makedirs(RESULT_DIR, exist_ok=True)
    os.makedirs(FIGURE_DIR, exist_ok=True)

    # Export Metadata Summary
    try:
        metadata_df.to_csv(
            output_file,
            index=False
        )

        print("✓ Metadata summary exported.")

        logging.info(f"Metadata summary exported successfully:\n{output_file}")

    except Exception as error:
        logging.exception("Failed to export metadata summary.")

        print(f"Error exporting metadata summary:\n{error}")

    # Export Frequency Tables
    export_mapping = {
        "Protein_Names":
            os.path.join(TABLE_DIR, "Protein_Name_Frequency.csv"),
        "Gene_Names":
            os.path.join(TABLE_DIR, "Gene_Name_Frequency.csv"),
        "Organisms":
            os.path.join(TABLE_DIR, "Organism_Frequency.csv"),
        "Taxonomy_IDs":
            os.path.join(TABLE_DIR, "Taxonomy_ID_Frequency.csv"),
        "Protein_Existence":
            os.path.join(TABLE_DIR, "Protein_Existence_Counts.csv"),
        "Missing_Values":
            os.path.join(TABLE_DIR, "Missing_Value_Summary.csv")
    }

    for table_name, dataframe in frequency_tables.items():
        try:
            if table_name in export_mapping:
                dataframe.to_csv(
                    export_mapping[table_name],
                    index=False
                )

                print(f"✓ {table_name} exported.")

                logging.info(f"{table_name} exported successfully.")

        except Exception as error:
            logging.exception(f"Failed to export {table_name}.")

            print(f"Error exporting {table_name}:\n{error}")

    # Export Sequence Statistics
    try:
        sequence_df = pd.DataFrame(
            sequence_statistics.items(),
            columns=[
                "Statistic",
                "Value"
            ]
        )

        sequence_df.to_csv(
            os.path.join(TABLE_DIR, "Sequence_Length_Statistics.csv"),
            index=False
        )

        print("✓ Sequence statistics exported.")

        logging.info("Sequence statistics exported successfully.")

    except Exception as error:
        logging.exception("Failed to export sequence statistics.")

        print(f"Error exporting sequence statistics:\n{error}")

    # Export Dataset Statistics
    try:
        dataset_rows = []

        for section, values in dataset_statistics.items():
            if isinstance(values, dict):
                for key, value in values.items():
                    if isinstance(value, dict):
                        
                        for sub_key, sub_value in value.items():
                            dataset_rows.append({
                                "Section": section,
                                "Metric":f"{key} ({sub_key})",
                                "Value": sub_value
                            })

                    else:
                        dataset_rows.append({
                            "Section": section,
                            "Metric": key,
                            "Value": value
                        })

        pd.DataFrame(dataset_rows).to_csv(
            os.path.join(TABLE_DIR, "Dataset_Statistics.csv"),
            index=False
        )

        print("✓ Dataset statistics exported.")

        logging.info("Dataset statistics exported successfully.")

    except Exception as error:
        logging.exception("Failed to export dataset statistics.")

        print(f"Error exporting dataset statistics:\n{error}")

    # Export Biological Statistics
    try:
        biological_rows = []

        for section, values in biological_statistics.items():
            for key, value in values.items():
                if isinstance(value, dict):
                    for sub_key, sub_value in value.items():
                        if isinstance(sub_value, dict):
                            
                            for x, y in sub_value.items():
                                biological_rows.append({
                                    "Section": section,
                                    "Metric": f"{key} | {sub_key} | {x}",
                                    "Value": y
                                })

                        else:
                            biological_rows.append({
                                "Section": section,
                                "Metric": f"{key} | {sub_key}",
                                "Value": sub_value
                            })

                else:
                    biological_rows.append({
                        "Section": section,
                        "Metric": key,
                        "Value": value
                    })

        pd.DataFrame(biological_rows).to_csv(os.path.join(TABLE_DIR, "Biological_Statistics.csv"),index=False)

        print("✓ Biological statistics exported.")

        logging.info("Biological statistics exported successfully.")

    except Exception as error:
        logging.exception("Failed to export biological statistics.")

        print(f"Error exporting biological statistics:\n{error}")

    print("\nAll results exported successfully.")

    logging.info("=" * 50)
    logging.info("Completed successfully.")
    logging.info("=" * 50)

    return

# Statistics & Experiment Reports
def generate_reports(metadata_df, dataset_statistics, biological_statistics, sequence_statistics, frequency_tables):

    print("\nGenerating Reports...")
    print("-" * 50)

    logging.info("Generating statistics and experiment reports...")

    # Statistics Report
    try:
        with open(statistics_report, "w", encoding="utf-8") as report:
            report.write("=" * 50 + "\n")
            report.write("EXP002D : METADATA STATISTICS REPORT\n")
            report.write("=" * 50 + "\n\n")

            # Dataset Statistics
            report.write("1. DATASET STATISTICS\n")
            report.write("-" * 50 + "\n")

            for section, values in dataset_statistics.items():
                report.write(f"\n{section}\n")

                if isinstance(values, dict):
                    for key, value in values.items():
                        if isinstance(value, dict):
                            report.write(f"\n{key}\n")

                            for sub_key, sub_value in value.items():
                                report.write(f"   {sub_key:<30}: {sub_value}\n")

                        else:
                            report.write(f"{key:<35}: {value}\n")

            # Biological Statistics
            report.write("\n\n")
            report.write("2. BIOLOGICAL STATISTICS\n")
            report.write("-" * 50 + "\n")

            for section, values in biological_statistics.items():
                report.write(f"\n{section}\n")

                for key, value in values.items():
                    if isinstance(value, dict):
                        report.write(f"\n{key}\n")

                        for sub_key, sub_value in value.items():
                            if isinstance(sub_value, dict):
                                report.write(f"\n{sub_key}\n")

                                for x, y in sub_value.items():
                                    report.write(f"   {x:<20}: {y}\n")

                            else:
                                report.write(f"   {sub_key:<30}: {sub_value}\n")

                    else:
                        report.write(f"{key:<35}: {value}\n")

            # Sequence Statistics
            report.write("\n\n")
            report.write("3. SEQUENCE STATISTICS\n")
            report.write("-" * 50 + "\n")

            for key, value in sequence_statistics.items():
                report.write(f"{key:<35}: {value}\n")

            report.write("\n\n")
            report.write("=" * 50 + "\n")
            report.write("END OF REPORT\n")
            report.write("=" * 50 + "\n")

        print("✓ Statistics report generated.")

        logging.info("Statistics report generated successfully.")

    except Exception as error:
        logging.exception("Failed to generate statistics report.")

        print(error)

    # Experiment Report
    try:
        with open(experiment_report, "w", encoding="utf-8") as report:

            report.write("=" * 50 + "\n")
            report.write("EXPERIMENT 002D REPORT\n")
            report.write("=" * 50 + "\n\n")
            
            report.write("Experiment Name\n")
            report.write("------------------------------\n")
            report.write("Metadata Statistics & Summary\n\n")
            report.write("Objective\n")
            report.write("------------------------------\n")
            
            report.write(
                "Generate descriptive statistics, "
                "frequency tables, visualizations "
                "and summary reports from the "
                "clean metadata generated in "
                "Experiment 002C.\n\n"
            )

            report.write("Input Dataset\n")
            report.write("------------------------------\n")
            report.write(f"{input_file}\n\n")

            report.write("Records Analysed\n")
            report.write("------------------------------\n")
            report.write(f"{len(metadata_df)}\n\n")

            report.write("Modules Executed\n")
            report.write("------------------------------\n")

            report.write("Module 1 : Configuration\n")
            report.write("Module 2 : Load Metadata\n")
            report.write("Module 3 : Dataset Statistics\n")
            report.write("Module 4 : Biological Statistics\n")
            report.write("Module 5 : Sequence Statistics\n")
            report.write("Module 6 : Frequency Tables\n")
            report.write("Module 7 : Visualizations\n")
            report.write("Module 8 : Export Results\n")
            report.write("Module 9 : Report Generation\n\n")

            report.write("Output Directories\n")
            report.write("------------------------------\n")

            report.write(f"Metadata : {OUTPUT_DIR}\n")
            report.write(f"Results  : {RESULT_DIR}\n")
            report.write(f"Figures  : {FIGURE_DIR}\n")
            report.write(f"Tables   : {TABLE_DIR}\n")
            report.write(f"Logs     : {LOG_DIR}\n\n")

            report.write("Generated Figures\n")
            report.write("------------------------------\n")

            for figure in sorted(os.listdir(FIGURE_DIR)):
                report.write(f"{figure}\n")

            report.write("\nGenerated Tables\n")
            report.write("------------------------------\n")

            for table in sorted(os.listdir(TABLE_DIR)):
                report.write(f"{table}\n")

            report.write("\nExperiment Status\n")
            report.write("------------------------------\n")
            report.write("SUCCESS\n\n")

            report.write("Completion Time : " + datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") + "\n")

            report.write("\n")
            report.write("=" * 80 + "\n")
            report.write("END OF EXPERIMENT REPORT\n")
            report.write("=" * 80 + "\n")

        print("✓ Experiment report generated.")

        logging.info("Experiment report generated successfully.")

    except Exception as error:
        logging.exception("Failed to generate experiment report.")

        print(error)

    print("\nAll reports generated successfully.")

    logging.info("=" * 50)
    logging.info("Experiment (Statistics & Experiment Reports) completed successfully.")
    logging.info("=" * 60)

# Main Function
def main():

    start_time = time.time()

    print("\nStarting Experiment 002D...")
    print("=" * 50)

    logging.info("=" * 50)
    logging.info("Experiment 002D execution started.")
    logging.info("=" * 50)

    try:
        metadata_df = load_metadata()
        dataset_statistics = generate_dataset_statistics(metadata_df)
        biological_statistics = generate_biological_statistics(metadata_df)
        sequence_statistics = generate_sequence_statistics(metadata_df)
        frequency_tables = generate_frequency_tables(metadata_df)

        generate_all_visualizations(metadata_df, frequency_tables)

        export_results(metadata_df, dataset_statistics, biological_statistics, sequence_statistics, frequency_tables)

        generate_reports(metadata_df, dataset_statistics, biological_statistics, sequence_statistics, frequency_tables)

        # Execution Summary
        execution_time = time.time() - start_time

        print("\n")
        print("=" * 50)
        print("Experiment 002D Completed Successfully")
        print("=" * 50)

        print(f"Execution Time : {execution_time:.2f} seconds")
        print(f"Records Processed : {len(metadata_df)}")

        print("\nOutput Directories")

        print(f"Metadata : {OUTPUT_DIR}")
        print(f"Results  : {RESULT_DIR}")
        print(f"Figures  : {FIGURE_DIR}")
        print(f"Tables   : {TABLE_DIR}")
        print(f"Logs     : {LOG_DIR}")

        print("\nAll outputs generated successfully.")

        logging.info("=" * 50)
        logging.info("Experiment 002D completed successfully.")
        logging.info(f"Execution Time : {execution_time:.2f} seconds")
        logging.info(f"Records Processed : {len(metadata_df)}")
        logging.info("=" * 50)

    except Exception as error:

        print("\nExperiment failed.")
        print(error)

        logging.exception("Experiment 002D terminated due to an unexpected error.")

        raise

# Execute the Main Function
if __name__ == "__main__":
    main()