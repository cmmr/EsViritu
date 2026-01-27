import argparse
import os
import glob
import polars as pl
import subprocess
import yaml

def summarize_esv_runs():
    """CLI entry point for the `summarize_esv_runs` console script."""
    parser = argparse.ArgumentParser(
        description="Summarize EsViritu run outputs across a directory."
    )
    parser.add_argument(
        "directory",
        type=str,
        help="Directory containing EsViritu .tsv outputs",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Directory to save output files (default: current working directory)",
    )
    args = parser.parse_args()
    
    outdir = args.outdir if args.outdir else os.getcwd()
    os.makedirs(outdir, exist_ok=True)
    # Define tsv file patterns to search for
    tsv_patterns = [
        "*.detected_virus.info.tsv",
        "*.detected_virus.assembly_summary.tsv",
        "*.virus_coverage_windows.tsv",
        "*.tax_profile.tsv"
    ]
    
    output_files = {}
    for pattern in tsv_patterns:
        file_paths = glob.glob(os.path.join(args.directory, pattern))
        if not file_paths:
            print(f"No files found for pattern: {pattern}")
            continue
        print(f"Found {len(file_paths)} files for pattern: {pattern}")
        if pattern == "*.tax_profile.tsv":
            spthresh_list = []
            subsp_thresh_list = []
            yaml_paths = glob.glob(os.path.join(args.directory, "*_esviritu.params.yaml"))
            for yaml_file_path in yaml_paths:
                with open(yaml_file_path, 'r') as yaml_file:
                    data = yaml.safe_load(yaml_file)
                    spthresh_list.append(data["spthresh"])
                    subsp_thresh_list.append(data["subspthresh"])

            if len(set(spthresh_list)) > 1:
                print("ERROR: disparate species-level thresholds detected in tax_profile tables, not combining.")
                print(f"values detected {set(spthresh_list)}")
                continue
            if len(set(subsp_thresh_list)) > 1:
                print("ERROR: disparate subspecies-level thresholds detected in tax_profile tables, not combining.")
                print(f"values detected {set(subsp_thresh_list)}")
                continue
                
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
            dir_basename = os.path.basename(os.path.abspath(args.directory))
            out_name = f"{dir_basename}.{pattern.replace('*.', '')}"
            out_path = os.path.join(outdir, out_name)
            merged.write_csv(out_path, separator="\t")
            print(f"Wrote merged file: {out_path}")
            output_files[pattern] = out_path

    # make input read summary from readstats.yaml
    read_dfs = []
    pattern = "*_esviritu.readstats.yaml"
    file_paths = glob.glob(os.path.join(args.directory, pattern))
    if not file_paths:
        print(f"No files found for pattern: {pattern}")
    else:
        print(f"Found {len(file_paths)} files for pattern: {pattern}")
        for yaml_file_path in file_paths:
            base = os.path.basename(yaml_file_path)
            sample_name = base.removesuffix("_esviritu.readstats.yaml")
            with open(yaml_file_path, 'r') as yaml_file:
                data = yaml.safe_load(yaml_file)
                data["sample_ID"] = sample_name

                read_df = pl.DataFrame(data)

                read_dfs.append(read_df)
    if read_dfs:
        merged_read_df = pl.concat(read_dfs, how="diagonal_relaxed")
        merged_read_df = merged_read_df.select(["sample_ID", pl.all().exclude("sample_ID")])
        dir_basename = os.path.basename(os.path.abspath(args.directory))
        out_name = f"{dir_basename}.readstats.tsv"
        out_path = os.path.join(outdir, out_name)
        merged_read_df.write_csv(out_path, separator="\t")
        print(f"Wrote merged file: {out_path}")
        output_files[pattern] = out_path
    # If both required summary tables were made, run the R script

    req_patterns = [
        "*.detected_virus.info.tsv",
        "*.virus_coverage_windows.tsv",
    ]
    if all(p in output_files for p in req_patterns):
        
        r_script = os.path.join(os.path.dirname(__file__), "EsViritu_project_reactable.R")
        coverage_tsv = output_files["*.virus_coverage_windows.tsv"]
        main_tsv = output_files["*.detected_virus.info.tsv"]
        # Add two empty string arguments for the 4th and 5th args if needed
        cmd = [
            "Rscript", r_script,
            coverage_tsv,
            main_tsv,
            os.path.basename(os.path.abspath(args.directory)),
            outdir
        ]
        try:
            print(f"Generating report: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
        except Exception as e:
            print(f"Failed to run R script: {e}")

if __name__ == "__main__":
    summarize_esv_runs()
