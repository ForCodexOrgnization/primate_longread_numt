#!/usr/bin/env python3
"""Print official ENA FASTQ metadata for a run and 1-based FASTQ index.

Output: tab-separated expected MD5, expected bytes, URL (no header).
Fail closed if the accession, URL, MD5 or byte count cannot be verified.
"""
import argparse
import json
import re
from urllib.parse import urlencode
from urllib.request import urlopen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--url", required=True, help="manifest HTTPS or FTP URL")
    args = ap.parse_args()
    query = urlencode({
        "accession": args.run,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
        "format": "json",
    })
    with urlopen("https://www.ebi.ac.uk/ena/portal/api/filereport?" + query,
                 timeout=60) as response:
        records = json.load(response)
    if len(records) != 1 or records[0].get("run_accession") != args.run:
        raise SystemExit("ENA metadata: unexpected run or record count")
    record = records[0]
    urls = record.get("fastq_ftp", "").split(";")
    md5s = record.get("fastq_md5", "").split(";")
    sizes = record.get("fastq_bytes", "").split(";")
    if not (len(urls) == len(md5s) == len(sizes)):
        raise SystemExit("ENA metadata: unequal URL/MD5/size list lengths")
    index = args.index - 1
    if index < 0 or index >= len(urls):
        raise SystemExit("ENA metadata: FASTQ index out of range")
    official_url = urls[index].strip()
    normalized = args.url.removeprefix("https://").removeprefix("http://")
    if official_url != normalized:
        raise SystemExit(f"ENA URL mismatch: {normalized!r} != {official_url!r}")
    md5 = md5s[index].strip().lower()
    size = sizes[index].strip()
    if not re.fullmatch(r"[0-9a-f]{32}", md5) or not size.isdigit() or int(size) <= 0:
        raise SystemExit("ENA metadata: invalid checksum or file size")
    print(f"{md5}\t{size}\t{official_url}")


if __name__ == "__main__":
    main()
