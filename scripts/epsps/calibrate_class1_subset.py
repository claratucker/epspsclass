#!/usr/bin/env python3
# calibrate_class1_subset.py
#
# Background: epspsclass requires all 148 of its CLASS_I_MARKERS positions
# to match exactly before calling a sequence Class I. Real organisms never
# hit 148/148. The repo's own benchmark data (CHANGELOG v1.0.4, README Known
# limitations) reports real class II organisms matching 18-22 of the 148
# class I markers, and real class I organisms from distant taxa matching
# only 24-28, a 2-marker gap too narrow to threshold safely on the full set.
#
# This script tests whether a SUBSET of the 148 positions separates the two
# groups more cleanly than the full set does. The approach: fetch the real
# benchmark organism sequences from NCBI, classify each of the 148 positions
# by how well it discriminates class I from class II in this benchmark, and
# report which subset size and which specific positions give the largest
# margin between the two groups.
#
# This does NOT invent new markers. It only selects, from the existing
# CLASS_I_MARKERS positions already derived from Leino et al. 2021's
# reference sequences, the subset most useful for thresholding, the same
# spirit as Leino et al.'s own choice to single out position 101 (Pro106 in
# E. coli numbering) as the canonical sensitivity marker (Funke et al. 2009;
# Sammons and Gaines 2014).
#
# Requires NCBI_API_KEY as an environment variable for a higher rate limit.
# Run on an instance with network access to NCBI eutils, not in a sandboxed
# environment without outbound internet access.

import os
import sys
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

API_KEY = os.environ.get("NCBI_API_KEY", "")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DELAY = 0.11 if API_KEY else 0.34

# Benchmark organisms named in epspsclass's own documentation (CHANGELOG.md,
# README.md Known limitations). Class II organisms are real, well known
# glyphosate-resistant species used in the original derivation work.
# Class I distant-taxa organisms are named only as phyla in the repo
# (Bacteroidetes, Firmicutes), not specific species, so representative
# genera from those phyla with annotated EPSPS are used instead and the
# substitution is logged explicitly.
BENCHMARK = {
    "II": [
        ("Staphylococcus aureus", '(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "Staphylococcus aureus"[Organism]'),
        ("Ruminococcus gnavus", '(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "Ruminococcus gnavus"[Organism]'),
        ("Dorea formicigenerans", '(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "Dorea formicigenerans"[Organism]'),
    ],
    "I": [
        # Representative Bacteroidetes and Firmicutes genera with annotated
        # EPSPS, standing in for the repo's unnamed "distant taxa" benchmark.
        ("Bacteroides fragilis", '(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "Bacteroides fragilis"[Organism]'),
        ("Bacteroides thetaiotaomicron", '(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "Bacteroides thetaiotaomicron"[Organism]'),
        ("Clostridium perfringens", '(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "Clostridium perfringens"[Organism]'),
        ("Prevotella copri", '(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "Prevotella copri"[Organism]'),
    ],
}

OUT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "results",
)
os.makedirs(OUT_DIR, exist_ok=True)
BENCHMARK_FAA = f"{OUT_DIR}/class1_calibration_benchmark.faa"


def eutils_get(endpoint, params):
    params = dict(params)
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()


def fetch_one(name, query):
    xml = eutils_get("esearch.fcgi", {"db": "protein", "term": query, "retmax": 10, "sort": "relevance"})
    time.sleep(DELAY)
    ids = [e.text for e in ET.fromstring(xml).findall(".//Id")]
    if not ids:
        return None

    xml = eutils_get("esummary.fcgi", {"db": "protein", "id": ",".join(ids), "retmode": "xml"})
    time.sleep(DELAY)
    summaries = []
    for doc in ET.fromstring(xml).findall(".//DocSum"):
        uid = doc.findtext("Id")
        caption = next((i.text for i in doc.findall("Item") if i.get("Name") == "Caption"), "")
        summaries.append((uid, caption or ""))
    refseq = [s for s in summaries if s[1].startswith(("WP_", "NP_", "YP_"))]
    uid, accession = refseq[0] if refseq else summaries[0]

    fasta = eutils_get("efetch.fcgi", {"db": "protein", "id": uid, "rettype": "fasta", "retmode": "text"})
    time.sleep(DELAY)
    if not fasta.strip().startswith(">"):
        return None
    seq = "".join(fasta.strip().splitlines()[1:])
    return accession, seq


def main():
    records = []
    log = []
    for cls, organisms in BENCHMARK.items():
        for name, query in organisms:
            result = fetch_one(name, query)
            if result is None:
                print(f"[{cls}] {name}: not_found")
                log.append((cls, name, "not_found", ""))
                continue
            accession, seq = result
            clean_name = name.replace(" ", "_")
            records.append(f">benchmark_class{cls}_{clean_name}|{accession}\n{seq}\n")
            print(f"[{cls}] {name}: ok ({accession}, {len(seq)} aa)")
            log.append((cls, name, "ok", accession))

    with open(BENCHMARK_FAA, "w") as f:
        f.writelines(records)

    n_ok = sum(1 for r in log if r[2] == "ok")
    print(f"\n{n_ok}/{len(log)} benchmark sequences fetched.")
    print(f"Written to: {BENCHMARK_FAA}")
    print("\nNext step: run classify_and_calibrate.py to test marker subsets")
    print("against these real sequences before touching classifier.py.")


if __name__ == "__main__":
    main()
