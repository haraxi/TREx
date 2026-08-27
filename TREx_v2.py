#!/usr/bin/env python3
"""
TREx: Taxonomic‑based Relative‑Abundance Extractor 
TREx_v2.py 

Designed based on TREx_v1.py, with additional features and improvements. 
Specifically, this version allows computation of RPKMF averages for output from EsViritu 1.x.x onwards.

"""

import argparse
import os
import glob
import re
from datetime import datetime
from dateutil import parser as date_parser
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend for plotting
import matplotlib.pyplot as plt
import pandas as pd


## ----------------- Helper Functions ----------------- ##

# ----------------- Date Parsing ----------------- #
def parse_date_input(date_args):
    if not date_args:
        return None, None
    try:
        start_date = date_parser.parse(date_args[0])
    except Exception:
        raise ValueError(f"Unable to parse start date: {date_args[0]}")
    if len(date_args) == 1:
        end_date = datetime.now()
    else:
        try:
            end_date = date_parser.parse(date_args[1])
        except Exception:
            raise ValueError(f"Unable to parse end date: {date_args[1]}")
    return start_date, end_date

#----------------- Accession Parsing ----------------- #
def accession_targets(acc_string):
    if not acc_string:
        return []
    targets = []
    for tok in acc_string.split("+"):
        tok = tok.strip().lower()
        if not tok:
            continue
        if ":" in tok:
            tok = tok.split(":", 1)[1]
        targets.append(tok)
    return targets


def parse_assembly_tokens(cell):
    if not cell or str(cell).lower() == "nan":
        return set()
    tokens = set()
    for tok in str(cell).split(","):
        tok = tok.strip().strip('"')
        if not tok:
            continue
        if ":" in tok:
            tok = tok.split(":", 1)[1]
        tokens.add(tok.strip().lower())
    return tokens


#----------------- Row Matching ----------------- #
def matches_search(row, args):
    # ---------------- accession filter ---------------- #
    if args.accession:
        tokens = parse_assembly_tokens(row.get("assembly_list", ""))
        if not (tokens & set(args._acc_list)):
            return False

    # ---------------- sample ID filter ---------------- #
    if args.sample_id:
        if str(row.get("Sample_ID", "")).strip().lower() != args.sample_id.strip().lower():
            return False
    return True

#----------------- File Processing ----------------- #
def process_file(tsv_file, args, date_range, metadata_folder):
    base = os.path.basename(tsv_file)
    m = re.match(r"esviritu_(p\d+).*\.tsv", base, re.IGNORECASE)
    if not m:
        return []
    pool_id = m.group(1)

    try:
        df = pd.read_csv(tsv_file, sep="\t")
    except Exception as e:
        print(f"Error reading {tsv_file}: {e}")
        return []

    meta_path = os.path.join(metadata_folder, f"{pool_id}_metadata.xlsx")
    meta_df = pd.read_excel(meta_path) if os.path.exists(meta_path) else None

    hits = []
    for _, row in df.iterrows():
        if not matches_search(row, args):
            continue

        raw_sid = str(row.get("sample_ID", ""))
        sample_id = raw_sid.split(".", 1)[0] if "." in raw_sid else raw_sid

        if meta_df is None:
            continue
        mdf = meta_df[meta_df["Sample_ID"].astype(str) == sample_id]
        if mdf.empty:
            continue
        meta = mdf.iloc[0]

        # ---------------- metadata filters ---------------- #
        if args.city and meta.get("City", "") not in args.city:
            continue
        if date_range[0] and date_range[1]:
            try:
                mdate = date_parser.parse(str(meta["Date"]))
            except Exception:
                continue
            if not (date_range[0] <= mdate <= date_range[1]):
                continue
        if args.site and meta.get("Site", "") != args.site:
            continue

        hits.append({
            "Site": meta.get("Site", ""),
            "City": meta.get("City", ""),
            "Date": meta.get("Date", ""),
            "Flow": meta.get("Flow", ""),
            "assembly_list": row.get("assembly_list", ""),
            "sequence_name": row.get("sequence_name", ""),
            "taxid": row.get("taxid", ""),
            "kingdom": row.get("kingdom", ""),
            "phylum": row.get("phylum", ""),
            "class": row.get("class", ""),
            "order": row.get("order", ""),
            "family": row.get("family", ""),
            "genus": row.get("genus", ""),
            "species": row.get("species", ""),
            "subspecies": row.get("subspecies", ""),
            "strain": row.get("strain", ""),
            "RPKMF": row.get("RPKMF", 0),
            "reference_length": row.get("reference_length", ""),
            "covered_bases": row.get("covered_bases", ""),
            "reads_aligned": row.get("reads_aligned", ""),
            "mean_coverage": row.get("mean_coverage", ""),
            "total_filtered_reads_in_sample": row.get("total_filtered_reads_in_sample", ""),
            "PoolID": pool_id,
            "Sample_ID": sample_id,
        })
    return hits

#----------------- Total Reads Lookup ----------------- #
def build_total_reads_lookup(tax_folder):
    lookup = {}
    pattern = os.path.join(tax_folder, "esviritu_p*.tax_profiles.tsv")
    for tsv in glob.glob(pattern):
        base = os.path.basename(tsv)
        m = re.match(r"esviritu_(p\d+)\.tax_profiles\.tsv", base, re.IGNORECASE)
        pool = m.group(1) if m else ""
        try:
            df = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        for _, row in df.iterrows():
            raw = str(row.get("sample_ID", ""))
            sample, pool_from_dot = (raw.split(".", 1)[0], raw.split(".", 1)[1]) if "." in raw else (raw, pool)
            key = (pool_from_dot, sample)
            val = row.get("filtered_reads_in_sample", None)
            if pd.notnull(val):
                lookup[key] = val
    return lookup


#----------------- RPKMF Averaging ----------------- #
def calculate_RPKMF_avg(df, metadata_folder, output_folder,
                        group_dim, time_res, allowed_cities=None):
    if df.empty:
        print("No data for RPKMF averaging.")
        return

    det = df.groupby(["PoolID","Sample_ID"], as_index=False)["RPKMF"].sum()

    # load metadata
    metas = []
    for mf in glob.glob(os.path.join(metadata_folder,"*_metadata.xlsx")):
        try:
            mdf = pd.read_excel(mf)
            name = os.path.basename(mf)
            pm = re.match(r"(p\d+)_metadata\.xlsx", name, re.IGNORECASE)
            mdf["PoolID"] = pm.group(1) if pm else ""
            metas.append(mdf)
        except Exception:
            pass
    if not metas:
        print("No metadata found for averaging.")
        return
    big = pd.concat(metas, ignore_index=True)

    merged = pd.merge(big, det, on=["PoolID","Sample_ID"], how="left")
    merged["RPKMF"] = merged["RPKMF"].fillna(0)
    merged["Date"] = pd.to_datetime(merged["Date"], errors="coerce")

    # time grouping
    tr = time_res.lower()
    if tr == "week":
        merged["Time_Group"] = merged["Date"].dt.to_period("W").apply(lambda r: r.start_time)
    elif tr == "month":
        merged["Time_Group"] = merged["Date"].dt.to_period("M").apply(lambda r: r.start_time)
    elif tr == "year":
        merged["Time_Group"] = merged["Date"].dt.to_period("Y").apply(lambda r: r.start_time)
    else:
        raise ValueError("Time must be week, month, or year")

    gd = group_dim.lower()
    if gd == "all":
        merged["All"] = "All"
        agg = merged.groupby(["All","Time_Group"], as_index=False)["RPKMF"].mean()
    elif gd == "city":
        if allowed_cities and len(allowed_cities)>1:
            merged["CityGroup"] = "All"
            group_col = "CityGroup"
        else:
            group_col = "City"
        agg = merged.groupby([group_col,"Time_Group"], as_index=False)["RPKMF"].mean()
        if allowed_cities:
            agg["Cities_Included"] = ", ".join(allowed_cities)
    elif gd == "site":
        agg = merged.groupby(["Site","Time_Group"], as_index=False)["RPKMF"].mean()
    else:
        raise ValueError("Group must be all, city, or site")

    out = os.path.join(output_folder, "RPKMF_averages.xlsx")
    agg.to_excel(out, index=False)
    print(f"Wrote RPKMF averages to {out}")


def plot_RPKMF_averages(output_folder, output_name):
    avg_path = os.path.join(output_folder, "RPKMF_averages.xlsx")
    if not os.path.exists(avg_path):
        print(f"No RPKMF_averages.xlsx found at {avg_path}")
        return

    df = pd.read_excel(avg_path)

    if "Time_Group" not in df.columns or "RPKMF" not in df.columns:
        print("Missing expected columns in RPKMF_averages.xlsx")
        return

    # Convert time to datetime
    df["Time_Group"] = pd.to_datetime(df["Time_Group"])

    # Determine group column
    group_col = None
    for col in ["All", "CityGroup", "City", "Site"]:
        if col in df.columns:
            group_col = col
            break

    if not group_col:
        print("No group column found to separate lines.")
        return

    # Plot setup
    plt.figure(figsize=(10, 6))
    for key, grp in df.groupby(group_col):
        plt.plot(grp["Time_Group"], grp["RPKMF"], marker="o", label=str(key))

    plt.xlabel("Time (MM/YYYY)")
    plt.ylabel("Average RPKMF")
    plt.title(f"RPKMF Averages Over Time - {output_name}")
    plt.legend()
    plt.tight_layout()

    # Format time axis
    plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%m/%Y'))
    plt.xticks(rotation=45)

    # Save as PNG and SVG
    png_path = os.path.join(output_folder, f"{output_name}_RPKMF_plot.png")
    svg_path = os.path.join(output_folder, f"{output_name}_RPKMF_plot.svg")
    plt.savefig(png_path)
    plt.savefig(svg_path)
    plt.close()
    print(f"Saved RPKMF plot to {png_path} and {svg_path}")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Search metagenomic TSVs, merge metadata, and optionally compute "
            "RPKMF averages (zeros for non-detects).  Multiple accession IDs "
            "can be supplied with '+' separators, e.g. 'NC_001802.1+AF052428.1'."
        )
    )
    group = parser.add_mutually_exclusive_group(required=True)
    #group.add_argument("--taxid", help="Taxonomic ID (non-sensitive)")
    #group.add_argument("--tax_name", help="Taxonomic name (non-sensitive)")
    group.add_argument("--accession", help="Accession(s) separated by '+', both modes")

    parser.add_argument("--sample_id", help="Filter by Sample ID")
    parser.add_argument("--date", nargs="+", help="Date or range: YYYY, MM-YYYY, DD-MM-YYYY")
    parser.add_argument("--city", nargs="+", help="City filter(s), e.g. 'Houston, TX'")
    parser.add_argument("--site", help="Site filter")
    parser.add_argument("--RPKMF_avg", nargs=2, metavar=("GROUP", "TIME"),
                        help="Compute RPKMF avg: GROUP=all|city|site, TIME=week|month|year")
    parser.add_argument("--output", help="Optional: custom output folder name (no spaces recommended)")


    args = parser.parse_args()

    # pre-compute list of accession targets for speed
    if args.accession:
        args._acc_list = accession_targets(args.accession)
        if not args._acc_list:
            parser.error("--accession provided but no valid IDs parsed.")
    else:
        args._acc_list = []

    date_range = parse_date_input(args.date) if args.date else (None, None)

    # Output folder logic (prioritizes --output if provided)
    if args.output:
        folder = args.output.strip().replace(" ", "_")
    elif args.tax_name:
        folder = args.tax_name.lower()
    elif args.accession:
        folder = args.accession.lower().replace("+", "_")
    elif args.taxid:
        folder = str(args.taxid).lower()
    else:
        folder = "output"
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.abspath(os.path.join(script_dir, ".."))
    input_folder = os.path.join(base_dir, "taxonomy") #Change as required to point to the correct input folder
    metadata_folder = os.path.join(base_dir, "metadata") #Change as required to point to the correct input folder
    output_folder = os.path.join(base_dir, folder)
    os.makedirs(output_folder, exist_ok=True)

    # process TSVs
    combined = []
    
    for f in glob.glob(os.path.join(input_folder, "esviritu_*.tax_profile.tsv")):
        combined.extend(process_file(f, args, date_range, metadata_folder))

    if not combined:
        print("No matching records found.")
        return

    # write combined Excel
    cols_norm = [
        "Site", "City", "Date", "assembly_list", "sequence_name", "taxid", "kingdom", "phylum",
        "class", "order", "family", "genus", "species", "subspecies", "strain", "RPKMF",
        "reference_length", "covered_bases", "reads_aligned", "mean_coverage",
        "total_filtered_reads_in_sample", "PoolID", "Sample_ID", "Flow",
    ]
    df_out = pd.DataFrame(combined, columns=(cols_norm))
    out_file = os.path.join(output_folder, "combined_data.xlsx")
    df_out.to_excel(out_file, index=False)
    print(f"Wrote combined data to {out_file}")

    # optional RPKMF averaging
    if args.RPKMF_avg:
        group_dim, time_res = args.RPKMF_avg
        allowed = args.city if group_dim.lower() == "city" else None
        calculate_RPKMF_avg(df_out, metadata_folder, output_folder, group_dim, time_res, allowed)
        plot_RPKMF_averages(output_folder, folder)



if __name__ == "__main__":
    main()
