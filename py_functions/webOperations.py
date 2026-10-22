import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import requests as rq
from bs4 import BeautifulSoup as bs
from urllib.parse import urljoin, urlparse
import pds4_tools
import re
from pathlib import PurePosixPath
import xml.etree.ElementTree as ET

def get_pds_csv_urls(base_url):

    # ----- Find sol directories -----
    response = rq.get(base_url)
    response.raise_for_status()

    soup = bs(response.text, "html.parser")

    sol_urls = []

    for link in soup.find_all("a", href=True):

        full_url = urljoin(base_url, link["href"])

        path = urlparse(full_url).path.rstrip("/")
        folder_name = path.split("/")[-1]

        if folder_name.startswith("sol_"):
            sol_urls.append(full_url)

    print(f"Found {len(sol_urls)} sol folders.")

    # ----- Find CSV files -----
    csv_urls = []
    xml_urls = []

    for sol_url in sol_urls:

        response = rq.get(sol_url)
        response.raise_for_status()

        soup = bs(response.text, "html.parser")

        for link in soup.find_all("a", href=True):

            href = link["href"]

            if href.lower().endswith(".csv"):
                csv_urls.append(urljoin(sol_url, href))
                xml_urls.append(urljoin(sol_url, href).rsplit(".", 1)[0] + ".xml")

    print(f"Found {len(csv_urls)} CSV files.")
    print(f"Found {len(xml_urls)} XML files.")

    return csv_urls,xml_urls




def get_pixl_metadata(xml_url):

    structures = pds4_tools.read(xml_url)
    label_dict = structures.label.to_dict()

    product = label_dict["Product_Observational"]

    # Product title
    title = product["Identification_Area"]["title"]

    # Sol number
    sol = (
        product["Observation_Area"]
        ["Mission_Area"]
        ["mars2020:Mars2020_Parameters"]
        ["mars2020:Observation_Information"]
        ["mars2020:sol_number"]
    )

    # Target name from title
    match = re.search(r"Product:\s*(.*?)\s*\(Sol", title)

    target_name = match.group(1) if match else None

    return {
        "target": target_name,
        "sol": int(sol),
        "title": title
    }



def get_rtt(url_or_filename):
    filename = PurePosixPath(url_or_filename).name
    return filename[35:44]



def get_all_urls(parent_url):
    visited_folders = set()
    file_urls = []

    def crawl(folder_url):
        if folder_url in visited_folders:
            return

        visited_folders.add(folder_url)

        response = rq.get(folder_url)
        response.raise_for_status()

        soup = bs(response.text, "html.parser")

        for link in soup.find_all("a", href=True):
            href = link["href"]

            # Ignore parent-directory links, anchors, etc.
            if href in ("../", "./", "/"):
                continue

            full_url = urljoin(folder_url, href)

            # Stay inside the original parent directory
            if not full_url.startswith(parent_url):
                continue

            parsed = urlparse(full_url)

            # Ignore query strings/fragments
            full_url = parsed._replace(
                query="",
                fragment=""
            ).geturl()

            # Directory links normally end with "/"
            if full_url.endswith("/"):
                crawl(full_url)
            else:
                file_urls.append(full_url)

    crawl(parent_url)

    return file_urls




def load_pixl_scan(oxide_url, rxl_url):

    oxide_filename = PurePosixPath(oxide_url).name
    rxl_filename = PurePosixPath(rxl_url).name

    oxide_rtt = get_rtt(oxide_filename)
    rxl_rtt = get_rtt(rxl_filename)

    if oxide_rtt != rxl_rtt:
        raise ValueError(
            f"RTT mismatch: "
            f"{oxide_rtt} != {rxl_rtt}"
        )

    oxide = pd.read_csv(oxide_url)
    rxl = pd.read_csv(rxl_url)

    oxide.columns = oxide.columns.str.strip().str.upper()
    rxl.columns = rxl.columns.str.strip().str.upper()

    oxide["RTT"] = oxide_rtt
    rxl["RTT"] = rxl_rtt

    combined = oxide.merge(
        rxl,
        on=["RTT", "PMC"],
        how="left",
        validate="one_to_one"
    )

    return combined


def download_file(url,cache_dir):
    filename = url.split("/")[-1]
    local_path = cache_dir / filename

    if not local_path.exists():
        response = rq.get(url)
        response.raise_for_status()

        local_path.write_bytes(response.content)

    return local_path


def parse_oxide_filename(filepath):
    try:
        filename = filepath.name

        # Sol: immediately after ps__
        sol_match = re.search(r"ps__(\d{4})_", filename, re.IGNORECASE)

        # Product type, e.g. rqa, rqb, rqc
        product_match = re.search(r"_000(rq[a-z])__", filename, re.IGNORECASE)

        if sol_match is None:
            raise ValueError(f"Could not determine SOL from {filename}")

        sol = int(sol_match.group(1))

        # RTT according to PIXL filename position
        rtt = filename[35:44]

        product_type = (
            product_match.group(1).lower()
            if product_match
            else None
        )

    except Exception as e:
        print(f"Problem in funct for {filepath.name}: {e}")

    return sol, rtt, product_type



def get_target_from_xml(xml_path):

    tree = ET.parse(xml_path)
    root = tree.getroot()

    # PDS4 uses XML namespaces, so search by the end of the tag name
    title = None

    for element in root.iter():
        if element.tag.endswith("title"):
            title = element.text
            break

    if title is None:
        return None

    match = re.search(
        r"Product:\s*(.*?)\s*\(Sol",
        title
    )

    return match.group(1) if match else None



def parse_rxl_filename(filepath):
    filename = filepath.name

    # Sol
    sol_match = re.search(r"__(\d{4})_", filename)

    if sol_match is None:
        raise ValueError(f"Could not determine SOL from {filename}")

    sol = int(sol_match.group(1))
    
    # Same RTT extraction convention as before
    rtt = filename[35:44]

    return sol, rtt


def get_sol_urls(base_url):
    response = rq.get(base_url, timeout=30)
    response.raise_for_status()

    soup = bs(response.text, "html.parser")

    sol_urls = []

    for link in soup.find_all("a", href=True):
        full_url = urljoin(base_url, link["href"])

        folder_name = (
            urlparse(full_url)
            .path
            .rstrip("/")
            .split("/")[-1]
        )

        if folder_name.startswith("sol_"):
            sol_urls.append(full_url)

    return sol_urls



def get_rxl_urls_from_sol(sol_url):
    response = rq.get(sol_url, timeout=30)
    response.raise_for_status()

    soup = bs(response.text, "html.parser")

    rxl_urls = []

    for link in soup.find_all("a", href=True):
        href = link["href"]

        if (
            "rxl" in href.lower()
            and href.lower().endswith(".csv")
        ):
            rxl_urls.append(urljoin(sol_url, href))

    return rxl_urls