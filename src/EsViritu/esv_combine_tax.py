#!/usr/bin/env python
"""
esv_combine_tax.py

Combine several EsViritu-style taxonomy tables (.tsv)

Usage:
  python esv_combine_tax.py taxon1.tsv taxon2.tsv -o outdir/virus_pathogen_database.all_metadata.tsv
"""

import argparse
import sys
import polars as pl
import os


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate esviritu-style taxonomy table from accession list.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('tables', nargs="+", help='tables (.tsv) to combine.')
    parser.add_argument(
        "-o", "--output",
        dest="outfile",
        default="virus_pathogen_database.all_metadata.tsv",
        help="Output TSV file for the combined taxonomy table."
    )

    return parser.parse_args()


def main():
    # Parse command-line arguments
    args = parse_args()

    schema_overrides = {"Genotype": str, "Segment": str, "TaxID": str}
    dfs = []
    for fp in args.tables:
        print(f'processing {fp}')
        try:
            
            df = pl.read_csv(
                fp, separator="\t",
                infer_schema_length=10000,
                schema_overrides=schema_overrides
                )
            ordered_fields = [
                "Accession", "description", "Name", "Segment",
                "kingdom", "phylum", "tclass", "order", 
                "family", "genus", "species", "subspecies", 
                "Length", "TaxID", "Assembly", "Asm_length"
                ]
            existing_ordered = [f for f in ordered_fields if f in df.columns]
            if existing_ordered != ordered_fields:
                sys.exit(
                    f'Error: missing one or more required fields in:\n{fp}\n\nrequired fields:\n{ordered_fields}'
                )
            remaining = [c for c in df.columns if c not in existing_ordered]
            df = df.select(existing_ordered + remaining)
            dfs.append(df)
        except Exception as e:
            print(f"Failed to read {fp}: {e}")
    if dfs:
        merged = pl.concat(dfs, how="vertical_relaxed")
        merged = merged.unique(maintain_order=True)

        merged.write_csv(args.outfile, separator="\t")
        print(f"Wrote merged file: {args.outfile}")

if __name__ == "__main__":
    main()