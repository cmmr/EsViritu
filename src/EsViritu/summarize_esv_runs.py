import argparse
import os
import glob
import polars as pl
import subprocess

def summarize_esv_runs(directory):
    # Define file types to search for
    file_patterns = [
        "*.detected_virus.info.tsv",
        "*.detected_virus.assembly_summary.tsv",
        "*.virus_coverage_windows.tsv"
    ]
    
    output_files = {}
    for pattern in file_patterns:
        file_paths = glob.glob(os.path.join(directory, pattern))
        if not file_paths:
            print(f"No files found for pattern: {pattern}")
            continue
        print(f"Found {len(file_paths)} files for pattern: {pattern}")
        dfs = []
        for fp in file_paths:
            try:
                df = pl.read_csv(
                    fp, separator="\t",
                    infer_schema_length=10000
                    )
                if pattern == "*.virus_coverage_windows.tsv":
                    base = os.path.basename(fp)
                    sample_name = base.removesuffix(".virus_coverage_windows.tsv")
                    df = df.with_columns(pl.lit(sample_name).alias("sample_ID"))
                dfs.append(df)
            except Exception as e:
                print(f"Failed to read {fp}: {e}")
        if dfs:
            merged = pl.concat(dfs, how="vertical_relaxed")
            dir_basename = os.path.basename(os.path.abspath(directory))
            out_name = f"{dir_basename}.{pattern.replace('*.', '')}"
            out_path = os.path.join(os.getcwd(), out_name)
            merged.write_csv(out_path, separator="\t")
            print(f"Wrote merged file: {out_path}")
            output_files[pattern] = out_path

    # If all three summary tables were made, run the R script
    if all(p in output_files for p in file_patterns):
        
        r_script = os.path.join(os.path.dirname(__file__), "EsViritu_project_reactable.R")
        coverage_tsv = output_files["*.virus_coverage_windows.tsv"]
        main_tsv = output_files["*.detected_virus.info.tsv"]
        # Add two empty string arguments for the 4th and 5th args if needed
        cmd = [
            "Rscript", r_script,
            coverage_tsv,
            main_tsv,
            os.path.basename(os.path.abspath(directory))
        ]
        try:
            print(f"Generating report: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
        except Exception as e:
            print(f"Failed to run R script: {e}")

def main():
    """CLI entry point for the `summarize_esv_runs` console script."""
    parser = argparse.ArgumentParser(
        description="Summarize EsViritu run outputs across a directory."
    )
    parser.add_argument(
        "directory",
        type=str,
        help="Directory containing EsViritu .tsv outputs",
    )
    args = parser.parse_args()
    summarize_esv_runs(args.directory)


if __name__ == "__main__":
    main()
