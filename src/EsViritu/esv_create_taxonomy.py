#!/usr/bin/env python
"""
esv_create_taxonomy.py

Generates an esviritu-style taxonomy table by:
  1. Joining with a user-provided accession table
  2. Loading an accession-to-taxid mapping
  3. Retrieving full lineage info via pytaxonkit
  4. Outputting a TSV with taxonomic ranks and metadata

Usage:
  python esv_create_taxonomy.py -i input.tsv -o output.tsv -a /path/to/acc2taxid.tsv.gz -t /path/to/taxdump
"""

import argparse
import sys
import polars as pl
import gzip
import os


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate esviritu-style taxonomy table from accession list.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use organism name as subspecies label
  python esv_create_taxonomy.py -i accessions.tsv -s organism -o output.tsv -a /path/to/acc2taxid.tsv.gz -t /path/to/taxdump

  # Use genotype as subspecies label
  python esv_create_taxonomy.py -i accessions.tsv -s genotype -o output.tsv -a /path/to/acc2taxid.tsv.gz -t /path/to/taxdump

  # Use NCBI default subspecies column (no -s flag)
  python esv_create_taxonomy.py -i accessions.tsv -o output.tsv -a /path/to/acc2taxid.tsv.gz -t /path/to/taxdump
"""
    )
    parser.add_argument(
        "-i", "--input",
        dest="select_tab",
        required=True,
        help="Input TSV file containing accessions (must have 'Accession', 'Organism_Name', and 'Length' columns). \
            Please include the 'Assembly' column if possible."
    )
    parser.add_argument(
        "-s", "--subspecies-label",
        dest="subspecies_lab",
        default = "",
        choices=["organism", "subspecies", "genotype", "species"],
        help="Default: The script will use pytaxonkit to label subspecies name (recommended). "
             "To use custom subspecies label: 'organism' (Organism_Name), "
             "'subspecies' (existing column), 'genotype' (Genotype), or 'species' (Species). "
    )
    parser.add_argument(
        "-o", "--output",
        dest="outfile",
        required=True,
        help="Output TSV file for the taxonomy table"
    )
    parser.add_argument(
        "-a", "--accfile",
        dest="accfile",
        required=True,
        help="path to pre-formatted acc2taxid.tsv.gz file"
    )
    parser.add_argument(
        "-t", "--taxonkit_dir",
        dest="taxonkit_dir",
        required=True,
        help="path to pre-formatted taxdump directory"
    )

    return parser.parse_args()


def main():
    # Parse command-line arguments
    args = parse_args()
    select_tab = args.select_tab
    subspecies_lab = args.subspecies_lab
    outfile = args.outfile
    taxonkit_dir = args.taxonkit_dir
    accfile = args.accfile

    try:
        import pytaxonkit as pytax
    except Exception as e:
        sys.exit(f'{e}\npytaxonkit can be installed via:\nconda install bioconda::pytaxonkit">=0.10"')



    ## Load table of selected accessions, handling subspecies label based on user choice
    # Common schema overrides for consistent column types
    schema_overrides = {"Genotype": str, "Segment": str, "TaxID": str}

    try:
        print("Loading accession table file...")
        acc_df = pl.read_csv(
            select_tab,
            has_header=True,
            separator="\t",
            schema_overrides=schema_overrides
        )    
    except FileNotFoundError:
        sys.exit(f"Error: Input file not found: {select_tab}")
    except pl.exceptions.ComputeError as e:
        sys.exit(f"Error parsing input TSV: {e}")
    except Exception as e:
        sys.exit(f"Error loading accession table: {e}")
    
    required_fields = ["Accession", "Organism_Name", "Length"]
    missing_cols = [f for f in required_fields if f not in acc_df.columns]

    if missing_cols:
        sys.exit(f"Error: Missing required columns {missing_cols} in accession table")


    if subspecies_lab == "organism":
        # Use Organism_Name column, prefixed with "t__"
        acc_df = acc_df.with_columns(
            pl.concat_str(
                [pl.lit("t"), pl.col("Organism_Name")],
                separator="__"
            ).alias("subspecies")
        )

    elif subspecies_lab == "subspecies":
        # Use existing subspecies column as-is
        if "subspecies" not in acc_df.columns:
            sys.exit("Error: -s subspecies specified but 'subspecies' column not found in input")
        print("Using existing subspecies column")

    elif subspecies_lab == "genotype":
        # Use Genotype column, prefixed with "t__"
        if "Genotype" not in acc_df.columns:
            sys.exit("Error: -s genotype specified but 'Genotype' column not found in input")
        acc_df = acc_df.with_columns(
            pl.concat_str(
                [pl.lit("t"), pl.col("Genotype")],
                separator="__"
            ).alias("subspecies")
        )

    elif subspecies_lab == "species":
        # use Species column, prefixed with "t__"
        if "Species" not in acc_df.columns:
            sys.exit("Error: -s species specified but 'Species' column not found in input")
        acc_df = acc_df.with_columns(
            pl.concat_str(
                [pl.lit("t"), pl.col("Species")],
                separator="__"
            ).alias("subspecies")
        )

    else:
        # no subspecies field yet - remove if present to avoid conflicts
        if "subspecies" in acc_df.columns:
            acc_df = acc_df.select(pl.exclude("subspecies"))

    ## Add Segment column as null if missing from input
    if 'Segment' not in acc_df.columns:
        acc_df = acc_df.with_columns(
            pl.lit(None).alias("Segment")
        )
    ## Add Assembly column as set:{Accession} if missing from input
    if 'Assembly' not in acc_df.columns:
        acc_df = acc_df.with_columns(
            pl.concat_str(
                [pl.lit("set"), pl.col("Accession")],
                separator=":"
            ).alias("Assembly")
        )

    ## fix NA values in assembly, rename Organism_Name
    acc_df = acc_df.with_columns(
        pl.when(
            pl.col("Assembly").cast(pl.String) == "NA"
        ).then(
            pl.concat_str(
                [pl.lit("set"), pl.col("Accession")],
                separator=":"
            )
        ).otherwise(
            pl.col("Assembly").cast(pl.String)
        ).alias("Assembly")
    ).rename(
        {"Organism_Name": "Name"}
    )

    ## Add description column as {Name} if missing from input
    if 'description' not in acc_df.columns:
        acc_df = acc_df.with_columns(pl.col("Name").alias("description"))

    ## remove pre-existing taxonomic levels
    acc_df = acc_df.select(pl.exclude(
        ["kingdom", "phylum", "tclass",
        "order", "family", "genus", "species"]
        )
    )

    ## Load accession-to-taxid mapping from gzip file

    acc2taxid_message = """
    INSTRUCTIONS TO DOWNLOAD AND FORMAT acc2taxid.tsv.gz:


    wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

    zcat nucl_gb.accession2taxid.gz | cut -f 2,3 | gzip -c > acc2taxid.tsv.gz
    ## need cols 2 (accession.version) and 3 (taxid)
    ## estimated file size = 900 MB
    """
    try:
        print("Loading acc2taxid file...")
        with gzip.open(accfile, 'rb') as f:
            acc2taxid_df = pl.read_csv(f, separator="\t")
    except FileNotFoundError:
        sys.exit(f"Error: acc2taxid file not found: {accfile}\n{acc2taxid_message}")
    except gzip.BadGzipFile:
        sys.exit(f"Error: Invalid gzip file: {accfile}")
    except Exception as e:
        sys.exit(f"Error loading acc2taxid table: {e}\n{acc2taxid_message}")

    ## Join accession table with acc2taxid to get taxonomy IDs
    taxid_df = acc_df.join(
        acc2taxid_df,
        how="left",
        left_on="Accession",
        right_on="accession.version"
    ).with_columns(
        pl.col("taxid").cast(pl.Int64).alias("taxid")
    )

    taxondir_message = """
    INSTRUCTIONS TO DOWNLOAD AND FORMAT taxdump directory:


    mkdir taxdump && cd taxdump

    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

    tar -xvf taxdump.tar.gz
    ## estimated directory size = 3 GB
    """

    ## Retrieve full lineage info using pytaxonkit
    if subspecies_lab == "":
        pytax_format_string = "k__{kingdom};p__{phylum};c__{class};o__{order};f__{family};g__{genus};s__{species};t__{subspecies|strain|no rank}"
        field_list = [
            "kingdom", "phylum", "tclass", "order", 
            "family", "genus", "species", "subspecies"
            ]
    else:
        pytax_format_string = "k__{kingdom};p__{phylum};c__{class};o__{order};f__{family};g__{genus};s__{species}"
        field_list = [
            "kingdom", "phylum", "tclass", "order", 
            "family", "genus", "species"
            ]
    
    try:
        print("Retrieving lineage info from pytaxonkit...")
        result_pl = pl.from_pandas(
            pytax.lineage(
                taxid_df["taxid"].to_list(),
                data_dir=taxonkit_dir,
                threads=os.cpu_count(),
                formatstr=pytax_format_string
            )
        )
    except pytax.TaxonKitCLIError as e:
        sys.exit(f"Error running taxonkit lineage: {e}")
    except FileNotFoundError:
        sys.exit(f"Error: taxonkit data directory not found: {taxonkit_dir}\n{taxondir_message}")
    except Exception as e:
        sys.exit(f"Error retrieving lineage info: {e}\n{taxondir_message}")

    result_pl = result_pl.with_columns(
        pl.col("TaxID").cast(pl.Int64).alias("TaxID")
    )

    ## Split lineage string into taxonomic rank columns
    try:
        full_tax_df = result_pl.join(
            taxid_df,
            how="left",
            left_on="TaxID",
            right_on="taxid"
        ).with_columns(
            pl.col("Lineage")
            .str.split_exact(";", 8)
            .struct.rename_fields(field_list)
            .alias("fields")
        ).unnest(
            "fields"
        ).unique(
        ).with_columns(
            # Sum Length for all records sharing the same Assembly identifier
            pl.col("Length").sum().over("Assembly").alias("Asm_length")
        ).select(
            pl.exclude(["Lineage"])
        )
    except pl.exceptions.ColumnNotFoundError as e:
        sys.exit(f"Error: Missing expected column during processing: {e}")
    except Exception as e:
        sys.exit(f"Error processing lineage data: {e}")

    ordered_fields = [
        "Accession", "description", "Name", "Segment",
        "kingdom", "phylum", "tclass", "order", 
        "family", "genus", "species", "subspecies", 
        "Length", "TaxID", "Assembly", "Asm_length"
        ]



    ## Reorder columns: ordered_fields first, then any remaining fields
    existing_ordered = [f for f in ordered_fields if f in full_tax_df.columns]
    remaining = [c for c in full_tax_df.columns if c not in existing_ordered]
    full_tax_df = full_tax_df.select(existing_ordered + remaining)

    ## Write output TSV
    try:
        full_tax_df.write_csv(
            outfile,
            separator="\t",
            quote_style='never'
        )
        print(f"Output written to: {outfile}")
    except PermissionError:
        sys.exit(f"Error: Permission denied writing to: {outfile}")
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")


if __name__ == "__main__":
    main()
