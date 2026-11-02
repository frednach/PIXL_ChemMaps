import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import requests as rq
from bs4 import BeautifulSoup as bs
from urllib.parse import urljoin, urlparse
import pds4_tools
import os
from py_functions import webOperations as webOp
from pathlib import PurePosixPath,Path
from concurrent.futures import ThreadPoolExecutor, as_completed


#import all .csv files that contain chemical composition data
#static url to the pds-database folder with the oxide data files
base_url_oxides = "https://pds-geosciences.wustl.edu/m2020/urn-nasa-pds-mars2020_pixl/data_oxides_pmc/"
pixl_oxides_urls, pixl_oxides_xmls = webOp.get_pds_csv_urls(base_url_oxides)
print(f"Found {len(pixl_oxides_urls)} oxide csv files")

#Download all relevant xml files
xml_cache = Path("pixl_cache/oxide_xml")
for url in pixl_oxides_xmls:
    webOp.download_file(url,xml_cache)
    
#Download all oxide data files
cache_dir = Path("pixl_cache/oxides")
with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:

    futures = {
        executor.submit(webOp.download_file, url, cache_dir): url
        for url in pixl_oxides_urls
    }

    for i, future in enumerate(as_completed(futures), 1):
        url = futures[future]

        try:
            future.result()
        except Exception as e:
            print(f"Failed: {url}: {e}")

        if i % 20 == 0:
            print(f"Downloaded {i}/{len(pixl_oxides_urls)}")

analysis_dir = f"pixl_cache/target_analysis"
oxide_dir = f"pixl_cache/oxide_xml"
oxides_dir = f"pixl_cache/oxides"
plots_dir = f"pixl_cache/plots"
processed_dir = f"pixl_cache/processed_data"
rxl_dir = f"pixl_cache/rxl"
target_dir = f"pixl_cache/target_analysis"
targets_dir = f"pixl_cache/targets"

directories = [analysis_dir, oxide_dir, oxides_dir, plots_dir, processed_dir, rxl_dir, target_dir, targets_dir]
for dir in directories:
    # Create Directory Structure if needed
    directory_path = Path(dir)
    try:
        directory_path.mkdir()
    except FileExistsError:
        continue
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_path}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

#static url to the pds-database folder with the rxl files
base_url_rxl = (
    "https://pds-geosciences.wustl.edu/m2020/"
    "urn-nasa-pds-mars2020_pixl/data_processed/"
)

# get number of subfolders (division by sol of observation)
sol_urls = webOp.get_sol_urls(base_url_rxl)
print(f"Found {len(sol_urls)} sol folders.")

# get a list of all available and relevant rxl files
all_rxl_urls = []
with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:

    futures = {
        executor.submit(webOp.get_rxl_urls_from_sol, sol_url): sol_url
        for sol_url in sol_urls
    }

    for i, future in enumerate(as_completed(futures), 1):
        sol_url = futures[future]

        try:
            rxl_urls = future.result()
            all_rxl_urls.extend(rxl_urls)

        except Exception as e:
            print(f"Problem reading {sol_url}: {e}")

        if i % 20 == 0:
            print(
                f"Processed {i}/{len(sol_urls)} folders | "
                f"RXL files found: {len(all_rxl_urls)}"
            )

print(f"Total RXL files: {len(all_rxl_urls)}")

# download all available rxl files
cache_dir = Path("pixl_cache/rxl")
with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:

    futures = {
        executor.submit(webOp.download_file, url, cache_dir): url
        for url in all_rxl_urls
    }

    for i, future in enumerate(as_completed(futures), 1):
        url = futures[future]

        try:
            future.result()
        except Exception as e:
            print(f"Failed: {url}: {e}")

        if i % 20 == 0:
            print(f"Downloaded {i}/{len(all_rxl_urls)}")

