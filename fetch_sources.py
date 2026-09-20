import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from comparability import config

SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML"
SUPPFILES = "https://www.ebi.ac.uk/europepmc/webservices/rest/%s/supplementaryFiles"

SOURCES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sources.tsv")
PAUSE = 0.5
TIMEOUT = 60
RETRIES = 3


def get(url):
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:
                return fh.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt == RETRIES - 1:
                raise
        except urllib.error.URLError:
            if attempt == RETRIES - 1:
                raise
        time.sleep(PAUSE * (attempt + 2))
    return None


def resolve_pmcid(doi):
    query = urllib.parse.urlencode(
        {"query": 'DOI:"%s"' % doi, "format": "json", "pageSize": "1"})
    body = get("%s?%s" % (SEARCH, query))
    if not body:
        return ""
    hits = json.loads(body.decode("utf-8")).get("resultList", {}).get("result", [])
    return hits[0].get("pmcid", "") if hits else ""


def main():
    out_xml = config.FULLTEXT_DIR
    out_supp = config.data_path("supplementary")
    os.makedirs(out_xml, exist_ok=True)
    os.makedirs(out_supp, exist_ok=True)

    with open(SOURCES, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    have, fetched, missing = 0, 0, []
    for row in rows:
        doi, pmcid = row["doi"], row["pmcid"]
        if not pmcid:
            pmcid = resolve_pmcid(doi)
            time.sleep(PAUSE)
        if not pmcid:
            missing.append(doi)
            continue

        target = os.path.join(out_xml, pmcid + ".xml")
        if os.path.exists(target):
            have += 1
        else:
            body = get(FULLTEXT % pmcid)
            if body is None:
                missing.append(doi)
                continue
            with open(target, "wb") as fh:
                fh.write(body)
            fetched += 1
            time.sleep(PAUSE)

        archive = os.path.join(out_supp, pmcid + ".zip")
        if not os.path.exists(archive):
            body = get(SUPPFILES % pmcid)
            if body:
                with open(archive, "wb") as fh:
                    fh.write(body)
            time.sleep(PAUSE)

    print("full text already present: %d" % have)
    print("full text downloaded: %d" % fetched)
    print("no open-access full text at Europe PMC: %d" % len(missing))
    for doi in missing:
        print("  %s  https://doi.org/%s" % (doi, doi))
    print("full text in ROSE_DATA/fulltext, supplementary archives in "
          "ROSE_DATA/supplementary")


if __name__ == "__main__":
    main()
