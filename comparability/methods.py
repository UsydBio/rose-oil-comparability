import re

from . import tables
from .samples import SPECIES


_METHODS_TITLE = re.compile(
    r"(material|method|experimental|sample\s*(collection|preparation)|"
    r"plant\s*material|reagent)", re.I)


_CULTIVAR_PATTERNS = [
    re.compile(r"[‘'\"“]([A-Z][\w .\-]{2,30})[’'\"”]"),
    re.compile(r"\bcv\.?\s+([A-Z][\w .\-]{2,30})"),
    re.compile(r"\bcultivars?\s+([A-Z][\w .\-]{2,30})"),
    re.compile(r"\bvar\.\s+([a-z][\w\-]{2,30})"),
]


_CULTIVAR_BLACKLIST = re.compile(
    r"^(the|this|that|and|or|not|see|table|figure|fig|supplementary|note|"
    r"control|blank|standard|method|methods|results|discussion|introduction|"
    r"n\.?d\.?|trace|total|others?|"

    r"hmisc|performanceanalytics|ggplot2?|tidyverse|dplyr|vegan|corrplot|"
    r"pheatmap|factoextra|mixomics|scipy|numpy|pandas|sklearn|"

    r"sour|sweet|fruity|green|fatty|woody|floral|spicy|citrus|herbal|minty|"
    r"honey|rosy|waxy|balsamic|earthy|musty|creamy|nutty)\b", re.I)


_NON_CULTIVAR_CONTEXT = re.compile(
    r"(package|software|library|version\s*\d|R\s*Core|descriptor|attribute|"
    r"odou?r\s*(note|term)|sensory\s*(term|attribute|descriptor))", re.I)

COUNTRIES = [
    "Bulgaria", "Turkey", "Türkiye", "Iran", "India", "China", "Morocco",
    "France", "Italy", "Spain", "Greece", "Egypt", "Saudi Arabia", "Japan",
    "Korea", "Republic of Korea", "Russia", "Ukraine", "Poland", "Germany",
    "Romania", "Serbia", "Albania", "Pakistan", "Afghanistan", "Taiwan",
    "Portugal", "Lithuania", "Estonia", "Tunisia", "Algeria", "Lebanon",
    "Syria", "Iraq", "Mexico", "Brazil", "Chile", "Australia", "Netherlands",
    "Switzerland", "Austria", "Hungary", "Bulgaria",
]
_COUNTRY_RE = re.compile(r"\b(%s)\b" % "|".join(sorted(set(COUNTRIES), key=len, reverse=True)))


REGIONS = {
    "kazanlak": "Kazanlak, Bulgaria", "kazanlik": "Kazanlak, Bulgaria",
    "isparta": "Isparta, Türkiye", "burdur": "Burdur, Türkiye",
    "kashan": "Kashan, Iran", "taif": "Taif, Saudi Arabia",
    "grasse": "Grasse, France", "pingyin": "Pingyin, China",
    "kushui": "Kushui (Yongdeng), China", "yunnan": "Yunnan, China",
    "hotan": "Hotan/Yutian, China", "khotan": "Hotan/Yutian, China",
}
_REGION_RE = re.compile(r"\b(%s)\b" % "|".join(REGIONS), re.I)

STAGES = [
    ("bud", re.compile(r"\b(bud|budding|unopened|closed\s*flower)\b", re.I)),
    ("half_open", re.compile(r"\b(half[\s-]*open|semi[\s-]*open|partially\s*open)\b", re.I)),
    ("full_bloom", re.compile(r"\b(full\s*(bloom|flower|open)|fully\s*open|anthesis)\b", re.I)),
    ("senescent", re.compile(r"\b(senescen|withered|fading|post[\s-]*anthesis)\b", re.I)),
]

_MONTH = re.compile(
    r"\b(January|February|March|April|May|June|July|August|September|"
    r"October|November|December)\b", re.I)


_TIME_OF_DAY = re.compile(
    r"\b([01]?\d|2[0-3])(?::([0-5]\d))\s*(?:a\.?m\.?|p\.?m\.?|h(?:ours?)?|o.clock)?\b"
    r"|\b([01]?\d|2[0-3])[:.]([0-5]\d)\s*(?:a\.?m\.?|p\.?m\.?|o.clock)\b", re.I)

_GPS = re.compile(r"(\d{1,3}[°º]\s*\d{1,2}['′]?\s*[\d.]*\s*[\"″]?\s*[NSEW])")


def methods_text(root):

    chunks = []
    for sec in root.iter("sec"):
        title = sec.find("title")
        label = tables.cell_text(title) if title is not None else ""
        if _METHODS_TITLE.search(label):
            chunks.append(tables.cell_text(sec))
    if not chunks:

        body = root.find(".//body")
        if body is not None:
            chunks.append(tables.cell_text(body))
    return "\n".join(chunks)


def cultivar_candidates(text):

    found = []
    for pattern in _CULTIVAR_PATTERNS:
        for m in pattern.finditer(text):
            name = m.group(1).strip(" .")
            if len(name) < 3 or _CULTIVAR_BLACKLIST.match(name):
                continue
            window = text[max(0, m.start() - 120):m.end() + 120]
            if _NON_CULTIVAR_CONTEXT.search(window):
                continue
            if name.lower() in SPECIES:
                continue
            if name not in found:
                found.append(name)
    return found



FIELDS = ["methods_species", "methods_species_all", "methods_species_n",
          "methods_cultivars", "methods_cultivar_n", "methods_country",
          "methods_region", "methods_stage", "methods_harvest_month",
          "methods_harvest_time", "methods_gps", "methods_processing",
          "methods_chars"]
