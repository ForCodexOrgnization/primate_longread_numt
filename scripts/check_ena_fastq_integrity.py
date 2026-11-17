#!/usr/bin/env python3
"""Compare an ENA FASTQ (including quarantined files) against ENA fastq_md5.

Usage: python3 scripts/check_ena_fastq_integrity.py --run SRR27194480
       --file data/.../read_1.fastq.gz.corrupt.<timestamp> --index 1
"""
import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", required=True)
    p.add_argument("--file", required=True)
    p.add_argument("--index", type=int, default=1, help="1-based FASTQ index")
    args = p.parse_args()
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?"
           + urlencode({"accession": args.run, "result": "read_run",
                        "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
                        "format": "json"}))
    with urlopen(url, timeout=60) as response:
        records = json.load(response)
    if len(records) != 1 or records[0].get("run_accession") != args.run:
        raise SystemExit("Unexpected ENA response: " + repr(records))
    record = records[0]
    urls = record.get("fastq_ftp", "").split(";")
    checksums = record.get("fastq_md5", "").split(";")
    sizes = record.get("fastq_bytes", "").split(";")
    if not (len(urls) == len(checksums) == len(sizes)):
        raise SystemExit("ENA URL/MD5/size list lengths differ")
    i = args.index - 1
    if i < 0 or i >= len(urls):
        raise SystemExit("FASTQ index outside ENA metadata")
    expected = checksums[i].strip().lower()
    if not re.fullmatch("[0-9a-f]{32}", expected):
        raise SystemExit("Missing or invalid ENA MD5: " + expected)
    path = Path(args.file)
    h = hashlib.md5()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    actual = h.hexdigest()
    expected_size = int(sizes[i])
    actual_size = path.stat().st_size
    print("Run:", args.run)
    print("ENA FASTQ URL:", urls[i])
    print("File:", path)
    print("ENA size:", expected_size, "local size:", actual_size)
    print("ENA MD5:", expected)
    print("Local MD5:", actual)
    print("SIZE_MATCH:", actual_size == expected_size)
    print("MD5_MATCH:", actual == expected)
    if actual != expected or actual_size != expected_size:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
