import argparse
import os
import glob
import polars as pl

def summarize_esv_runs(directory):
    # Define file types to search for
    file_patterns = [
        "*.detected_virus.info.tsv",
        "*.detected_virus.assembly_summary.tsv",
        "*.virus_coverage_windows.tsv"
    ]
    
    for pattern in file_patterns:
        file_paths = glob.glob(os.path.join(directory, pattern))
        if not file_paths:
            print(f"No files found for pattern: {pattern}")
            continue
        print(f"Found {len(file_paths)} files for pattern: {pattern}")
        dfs = []
        for fp in file_paths:
            try:
                df = pl.read_csv(fp, separator="\t")
                dfs.append(df)
            except Exception as e:
                print(f"Failed to read {fp}: {e}")
        if dfs:
            merged = pl.concat(dfs, how="vertical_relaxed")
            dir_basename = os.path.basename(os.path.abspath(directory))
            out_name = f"{dir_basename}.{pattern.replace('*.', '')}"
            out_path = os.path.join(directory, out_name)
            merged.write_csv(out_path, separator="\t")
            print(f"Wrote merged file: {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize EsViritu run outputs across a directory.")
    parser.add_argument("directory", type=str, help="Directory containing EsViritu .tsv outputs")
    args = parser.parse_args()
    summarize_esv_runs(args.directory)
