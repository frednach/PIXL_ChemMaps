import py_functions.webOperations as webOp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import PurePosixPath,Path
import re

#create oxide data dataFrame
#path to folder of oxides and corresponding xml files
oxide_dir = Path("pixl_cache/oxides")
xml_dir = Path("pixl_cache/oxide_xml")

all_oxide_data = []

for filepath in oxide_dir.glob("*.csv"):

    try:
        sol, rtt, product_type = webOp.parse_oxide_filename(filepath)

        # Find corresponding XML file
        xml_path = xml_dir / filepath.with_suffix(".xml").name
        
        if xml_path.exists():
            target = webOp.get_target_from_xml(xml_path)
        else:
            target = None
            print(f"No XML found for {filepath.name}")
            print("XML:",xml_path)

        # Read oxide data
        df = pd.read_csv(filepath)
        
        df.columns = df.columns.str.strip()
        
        # Add metadata
        df["SOL"] = sol
        df["RTT"] = rtt
        df["TARGET"] = target
        df["PRODUCT_TYPE"] = product_type
        df["SOURCE_FILE"] = filepath.name

        all_oxide_data.append(df)

    except Exception as e:
        print(f"Problem with {filepath.name}: {e}")

oxide_data = pd.concat(all_oxide_data, ignore_index=True)
#reorder dataframe structure to have identifiers at the front
identifier_cols = [
    "SOL",
    "TARGET",
    "RTT",
    "PMC",
    "PRODUCT_TYPE",
    "SOURCE_FILE"
]
other_cols = [
    col for col in oxide_data.columns
    if col not in identifier_cols
]
oxide_data = oxide_data[
    identifier_cols + other_cols
]
#inspect the data
# #print("SHAPE:",oxide_data.shape)
# #print("COLUMNS:",oxide_data.columns)
# #print("HEAD:",oxide_data.head())

#define output directory
output_dir = Path("pixl_cache/processed_data")
#create name for output file
output_file = output_dir / "oxide_data_master.csv"
#create .csv and save to specified location
oxide_data.to_csv(
    output_file,
    index=False
)
#print confirmation
print(f"Saved master oxide table to: {output_file}")
print("OxideData has",len(oxide_data),"rows")

#
oxide_keys = set(
    zip(
        oxide_data["SOL"],
        oxide_data["RTT"],
        oxide_data["PMC"]
    )
)

#Do the same for the rxl files
#rxl directory
rxl_dir = Path("pixl_cache/rxl")

oxide_key_df = (
    oxide_data[
        ["SOL", "RTT", "PMC"]
    ]
    .drop_duplicates()
)

all_rxl_data = []

for filepath in rxl_dir.glob("*.csv"):

    try:
        sol, rtt = webOp.parse_rxl_filename(filepath)

        df = pd.read_csv(filepath)
        df.columns = df.columns.str.strip()

        df["SOL"] = sol
        df["RTT"] = rtt
        df["SOURCE_FILE_RXL"] = filepath.name

        # Retain only points with oxide data
        df = df.merge(
            oxide_key_df,
            on=["SOL", "RTT", "PMC"],
            how="inner"
        )

        if not df.empty:
            all_rxl_data.append(df)

    except Exception as e:
        print(f"Problem with {filepath.name}: {e}")


#combine the data into one dataframe
rxl_data = pd.concat(
    all_rxl_data,
    ignore_index=True
)
#rearrange dataFrame to put identifiers first
identifier_cols = [
    "SOL",
    "RTT",
    "PMC",
    "SOURCE_FILE_RXL"
]
other_cols = [
    col for col in rxl_data.columns
    if col not in identifier_cols]

rxl_data = rxl_data[identifier_cols + other_cols]


#save file as csv
output_dir = Path("pixl_cache/processed_data")
rxl_output = output_dir / "rxl_data_master.csv"
rxl_data.to_csv(rxl_output,index=False)
#print confirmation
print(f"RXL Data has {len(rxl_data)} rows in {rxl_output}")

oxide_unique = (
    oxide_data[
        ["SOL", "RTT", "PMC"]
    ]
    .drop_duplicates()
)

rxl_unique = (
    rxl_data[
        ["SOL", "RTT", "PMC"]
    ]
    .drop_duplicates()
)

print("Unique oxide points:", len(oxide_unique))
print("Unique RXL points:", len(rxl_unique))

#just use data observed by both detectors = rqc
oxide_data_rqc = oxide_data[
    oxide_data["PRODUCT_TYPE"].str.lower() == "rqc"
].copy()


#create dataframe of unique oxide data only
output_dir = Path("pixl_cache/processed_data")
ox_output = output_dir / "oxide_unique.csv"
oxide_unique.to_csv(ox_output,index=False)

#join the two files into a single one
pixl_merged_c = oxide_data_rqc.merge(
    rxl_data,
    on=["SOL", "RTT", "PMC"],
    how="inner",
    suffixes=("_OXIDE", "_RXL")
)
# save to .csv
output_dir = Path("pixl_cache/processed_data")
output_file = output_dir / "master_merged.csv"
pixl_merged_c.to_csv(
    output_file,
    index=False
)
print(f"Saved merged dataset to: {output_file}")
#save to parquet for faster loading later
pixl_merged_c.to_parquet(
    output_dir / "master_merged_c.parquet",
    index=False
)
#check size after merge
print(
    "Unique merged PIXL points:",
    len(
        pixl_merged_c[
            ["SOL", "RTT", "PMC"]
        ].drop_duplicates()
    )
)

