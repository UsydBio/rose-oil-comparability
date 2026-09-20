import re

import xml.etree.ElementTree as ET

from . import analytical


CENSUS_IS_UPPER_BOUND = (
    "Every percentage in this census counts a study as reporting an item when "
    "the item is named anywhere in the searched text. It does not check that "
    "the value is correct, that it applies to the quantitative table, or that "
    "it is complete. Reporting rates are therefore upper bounds.")


NON_VOLATILE_STUDIES = {
    "PMC10572640",
    "PMC11557480",
    "PMC11554761",
}


_REVIEW_SELF = re.compile(
    r"th(?:is|e\s+present)\s+review\s+(?:provides|presents|summari|aims|"
    r"critically|examines|covers|discusses)|"
    r"in\s+this\s+review\s*,?\s*we", re.I)


def is_review(root):

    if root is None:
        return False, ""
    if root.get("article-type", "") == "review-article":
        return True, 'article-type="review-article"'
    ab = root.find(".//abstract")
    if ab is not None:
        text = " ".join(x.strip() for x in ab.itertext())
        m = _REVIEW_SELF.search(text)
        if m:
            return True, m.group(0)
    return False, ""


def wrap_text(text, title="Materials and methods"):

    art = ET.Element("article")
    body = ET.SubElement(art, "body")
    sec = ET.SubElement(body, "sec")
    ET.SubElement(sec, "title").text = title
    ET.SubElement(sec, "p").text = text
    return art


EXTRA_DETECTORS = [
    ("GC-MS", re.compile(
        r"(?<![A-Za-z])gc\s*[-–—/×x]?\s*(?:tof|q\s*-?\s*tof|q|qq?q|it|"
        r"time\s*-?\s*of\s*-?\s*flight|quadrupole|ion\s*trap)\s*[-–—/]?\s*ms"
        r"(?![A-Za-z])|"
        r"gas\s+chromatograph\w*[^.]{0,40}?time\s*-?\s*of\s*-?\s*flight", re.I)),
    ("GC-IMS", re.compile(
        r"(?<![A-Za-z])gc\s*[-–—/]\s*ims(?![A-Za-z])|"
        r"ion\s+mobility\s+spectrometr", re.I)),
]


def detectors(text):

    found = [p for p in analytical._platforms_in(text)
             if analytical._PLATFORM_BY_NAME[p][0] == "detection"]
    for name, pat in EXTRA_DETECTORS:
        if pat.search(text) and name not in found:
            found.append(name)
    return found


def _window(text, start, end, span):
    return text[max(0, start - span):end + span]


def _hit(pattern, text, context=None, span=200, veto=None, veto_span=None):

    for m in pattern.finditer(text):
        if context is not None and not context.search(
                _window(text, m.start(), m.end(), span)):
            continue
        if veto is not None and veto.search(
                _window(text, m.start(), m.end(), veto_span or span)):
            continue
        return re.sub(r"\s+", " ",
                      _window(text, m.start(), m.end(), 90)).strip()
    return ""


_RF_ANY = re.compile(r"response\s+factors?|correction\s+factors?", re.I)
_RF_DECLINED = re.compile(

    r"without\s+(?:\w+ing\s+)?(?:the\s+use\s+of\s+)?(?:any\s+)?"
    r"(?:response\s+factor|correction)|no\s+response\s+factors?|"
    r"response\s+factors?\s+(?:were|was)\s+not|"
    r"(?:assum\w+|taken|set|regard\w*)[^.]{0,45}response\s+factors?"
    r"[^.]{0,25}(?:of\s+)?(?:1(?:\.0+)?|one|unity|equal)|"
    r"(?:equal|identical|unity)\s+response\s+factors?|"
    r"uncorrected\s+(?:peak\s+)?areas?|no\s+correction\s+factors?", re.I)
_RF_APPLIED = re.compile(
    r"response\s+factors?[^.]{0,60}(?:appli|us(?:ed|ing)|correct|calculat|"
    r"determin|obtain|from\s+the\s+literature|taken\s+from)|"
    r"(?:appli|us(?:ed|ing)|correct\w*|calculat\w*|multipl\w*)[^.]{0,60}"
    r"response\s+factors?", re.I)


def response_factor(text):

    ev = _hit(_RF_DECLINED, text)
    if ev:
        return "declined", ev
    ev = _hit(_RF_APPLIED, text)
    if ev:
        return "applied", ev
    return "", ""


_REPLICATE = re.compile(
    r"in\s+(?:triplicate|duplicate|quadruplicate)|"
    r"(?:triplicate|duplicate|quadruplicate)s?\s+(?:analys|measur|determin|"
    r"extract|injection|sample|run)|"
    r"(?:three|two|four|five|six|\d{1,2})\s+(?:biological\s+|technical\s+|"
    r"independent\s+|parallel\s+)?(?:replicat\w+|repetitions?|repeats?|"
    r"determinations?)|"
    r"(?:replicat\w+|repetitions?)\s*(?:\(\s*n\s*=\s*\d+\s*\)|n\s*=\s*\d+)|"
    r"(?<![A-Za-z])n\s*=\s*[2-9]\d?(?![\d%])", re.I)


_REPLICATE_CHEM = re.compile(
    r"chromatograph|(?<![A-Za-z])GC(?![A-Za-z])|essential\s+oil|volatil|"
    r"distill|extract\w*|aroma|injection|SPME|headspace|compound", re.I)


_REPLICATE_NOT = re.compile(
    r"antioxidant|DPPH|ABTS|FRAP|cytotox|(?<![A-Za-z])cells?(?![A-Za-z])|"
    r"antimicrob|antibacter|antifung|inhibit\w*\s+zone|(?<![A-Za-z])MIC"
    r"(?![A-Za-z])|disc\s+diffusion|sensory|panel\w*|viability|"
    r"(?<![A-Za-z])qRT?\s*-?\s*PCR(?![A-Za-z])|expression|"

    r"disease|plantlet|seedling|inocul|infect|lesion|germinat|"
    r"pot\s+experiment|greenhouse|field\s+trial", re.I)


def replicates(text):

    loose = _hit(_REPLICATE, text)
    strict = _hit(_REPLICATE, text, context=_REPLICATE_CHEM, span=300,
                  veto=_REPLICATE_NOT, veto_span=220)
    return loose, strict


_HARVEST_CTX = re.compile(
    r"harvest|collect|pick(?:ed|ing)?|gather|plucked|"
    r"flowers?\s+(?:was|were)|petals?\s+(?:was|were)|"
    r"sampl\w+\s+(?:was|were)\s+(?:taken|obtained|collected)|"
    r"time\s+points?|during\s+the\s+(?:flowering|blooming)", re.I)


_MONTH = re.compile(
    r"(?<![A-Za-z])(January|February|March|April|June|July|August|September|"
    r"October|November|December)(?![A-Za-z])", re.I)
_MONTH_MAY = re.compile(r"(?<![A-Za-z])May(?![A-Za-z])")
_YEAR = re.compile(r"(?<![\d.])(19[5-9]\d|20[0-4]\d)(?![\d.])")
_DAY = re.compile(
    r"(?<![\d.])([1-3]?\d)(?:st|nd|rd|th)?\s+(?:of\s+)?"
    r"(?:January|February|March|April|May|June|July|August|September|"
    r"October|November|December)|"
    r"(?:January|February|March|April|May|June|July|August|September|"
    r"October|November|December)\s+([1-3]?\d)(?:st|nd|rd|th)?(?![\d])", re.I)


_TOD = re.compile(
    r"(?<![\d.:])([01]?\d|2[0-3])\s*[:.]\s*([0-5]\d)\s*"
    r"(?:a\.?m\.?|p\.?m\.?|h(?:ours?)?|o.clock)|"
    r"(?<![\d.:])([01]?\d|2[0-3])\s*(?:a\.?m\.?|p\.?m\.?)(?![A-Za-z])|"
    r"(?<![\d.:])([01]?\d|2[0-3])\s*[:.]\s*([0-5]\d)\s*(?=[-–—]|\s+(?:to|and)\s+)"
    r"\s*(?:[01]?\d|2[0-3])\s*[:.]\s*[0-5]\d", re.I)
_CITATION = re.compile(r"\(\s*(?:19|20)\d\d\s*\)|Page\s+\d|vol\.|pp\.|doi", re.I)


def harvest_date(text):

    month = (_hit(_MONTH, text, _HARVEST_CTX, 250)
             or _hit(_MONTH_MAY, text, _HARVEST_CTX, 250))
    year = _hit(_YEAR, text, _HARVEST_CTX, 250, veto=_CITATION, veto_span=40)
    day = _hit(_DAY, text, _HARVEST_CTX, 250)
    return month, year, day


def harvest_time(text):

    return _hit(_TOD, text, _HARVEST_CTX, 140, veto=_CITATION, veto_span=45)


_SUPPLIER = re.compile(
    r"sigma|aldrich|merck|supelco|restek|agilent|shimadzu|thermo\s|fisher|"
    r"waters\s+corp|bruker|perkin|varian|hewlett|packard|gmbh|operon|"
    r"co\.\,?\s*ltd|(?<![A-Za-z])inc(?![A-Za-z])|s\.\s*a\.|b\.\s*v\.|"
    r"purchas|suppli(?:ed|er)|manufactur|reagent|chemicals?\s+(?:were|was)|"
    r"solvents?\s+(?:were|was)|standards?\s+(?:were|was)\s+(?:purchas|obtain)|"
    r"software|instrument|apparatus|incubator|centrifug|spectrometer|"
    r"lyophili|freeze[\s-]*dr\w*\s*\(|shaker|oven\s*\(|balance\s*\(|"
    r"(?<![A-Za-z])kits?(?![A-Za-z])|cell\s+lines?|(?<![A-Za-z])strains?"
    r"(?![A-Za-z])|(?<![A-Za-z])medium(?![A-Za-z])|"

    r"(?<![A-Za-z])[A-Z]{2,}[\s-]?\d{3,}(?![\d])", re.I)


_COLLECT_CTX = re.compile(
    r"harvest|collect(?:ed|ion)|grown|cultivat|plantation|orchard|"
    r"originat|native\s+to|growing\s+(?:region|area)|picked|gathered|"
    r"sampling\s+site|experimental\s+(?:station|plot|area|field)|"
    r"planted|nursery|(?:rose|rosa|flower|petal|plant|hip|fruit|leaf|"
    r"cultivar|accession|landrace|germplasm)\w*", re.I)
_GPS = re.compile(
    r"\d{1,3}\s*[°º]\s*\d{1,2}\s*['′]?\s*[\d.]*\s*[\"″]?\s*[NSEW]|"
    r"(?<![\d.])[-+]?\d{1,2}\.\d{3,7}\s*[°º]?\s*[NS]\s*[,;]\s*"
    r"[-+]?\d{1,3}\.\d{3,7}", re.I)


_ORIGIN_ACQ = re.compile(
    r"(?:rose|rosa|flower|petal|blossom|hip|fruit|leaf|leaves|plant\s+material|"
    r"sample|material|specimen|cultivar|accession|genotype|landrace|"
    r"germplasm|oil)s?[^.]{0,140}?(?:were|was|are|is)\s+(?:kindly\s+)?"
    r"(?:being\s+)?(?:collect|harvest|obtain|grow|cultivat|provid|pick|"
    r"gather|suppl|purchas|plant|harvest)\w*|"
    r"(?:collected|harvested|obtained|grown|cultivated|provided|picked|"
    r"gathered|originating|sourced)\s+(?:from|in|at|by)(?![A-Za-z])|"
    r"(?:grow|cultivat)\w*\s+(?:in|at|near)(?![A-Za-z])", re.I)
_QUOTED_RESULT = re.compile(
    r"\[\s*\d|et\s+al\.|according\s+to|reported\s+(?:by|that|in)|"
    r"previous\s+stud|other\s+(?:studies|authors)|"
    r"(?:University|Faculty|Department|Institute\s+of)\s+of\s+[A-Z]|"
    r"^\s*\d\s*(?:Department|Faculty|School|College)", re.I | re.M)


def geographic_origin(text, country_re, region_re):

    country = ""
    for m in country_re.finditer(text):
        near = _window(text, m.start(), m.end(), 200)
        if not _ORIGIN_ACQ.search(near):
            continue
        if _SUPPLIER.search(_window(text, m.start(), m.end(), 110)):
            continue
        if _QUOTED_RESULT.search(_window(text, m.start(), m.end(), 130)):
            continue
        country = re.sub(r"\s+", " ",
                         _window(text, m.start(), m.end(), 90)).strip()
        break
    region = _hit(region_re, text, _COLLECT_CTX, 250,
                  veto=_QUOTED_RESULT, veto_span=130)
    gps = _hit(_GPS, text)
    return country, region, gps


_EXTRACTION_CTX = re.compile(
    r"distill|extract\w*|clevenger|macerat|steep|soak|percolat|"
    r"(?<![A-Za-z])SPME(?![A-Za-z])|headspace|equilibrat|absorb\w*\s+onto", re.I)
_EXTRACTION_TIME = re.compile(
    r"(?:for|during|over|of)\s+(\d{1,3}(?:\.\d+)?)\s*"
    r"(h(?:ours?|rs?)?|min(?:utes?)?)(?![A-Za-z])", re.I)
_EXTRACTION_TEMP = re.compile(
    r"(?:at|to|of|and)\s+(\d{1,3})\s*(?:±\s*\d+\s*)?[°º℃]\s*C(?![a-z])|"
    r"room\s+temperature|ambient\s+temperature", re.I)


_INSTRUMENT_TEMP = re.compile(
    r"injector|inlet|injection\s+port|oven|transfer\s+line|ion\s+source|"
    r"detector|desorb|interface|held\s+at\s+\d{2,3}\s*[°º]?\s*C\s+for\s+\d+\s*min",
    re.I)
_RATIO = re.compile(
    r"(?<![\d.])(\d{1,3})\s*[:：]\s*(\d{1,3})(?![\d:])"
    r"\s*(?:\(?\s*(?:w\s*/\s*v|v\s*/\s*w|w\s*/\s*w|g\s*/\s*m?L|m?L\s*/\s*g)"
    r"\s*\)?)?", re.I)


_RATIO_PLANT = re.compile(
    r"flowers?|petals?|plant\s+material|herb|biomass|sample|solid|material", re.I)
_RATIO_LIQUID = re.compile(
    r"water|liquid|solvent|distilled\s+water|hydro", re.I)
_RATIO_WORD = re.compile(r"ratio|module|proportion", re.I)


def extraction_params(text):

    time_ev = _hit(_EXTRACTION_TIME, text, _EXTRACTION_CTX, 160)
    temp_ev = _hit(_EXTRACTION_TEMP, text, _EXTRACTION_CTX, 130,
                   veto=_INSTRUMENT_TEMP, veto_span=90)
    ratio_ev = ""
    for m in _RATIO.finditer(text):
        win = _window(text, m.start(), m.end(), 120)
        if (_RATIO_WORD.search(win) and _RATIO_PLANT.search(win)
                and _RATIO_LIQUID.search(win)):
            ratio_ev = re.sub(r"\s+", " ",
                              _window(text, m.start(), m.end(), 90)).strip()
            break
    return time_ev, temp_ev, ratio_ev


_MATERIAL_CTX = re.compile(
    r"harvest|collect|pick(?:ed|ing)?|gather|plucked|plant\s+material|"
    r"sampl\w+|subjected\s+to|distill|extract\w*|air[\s-]*dr|fresh(?:ly)?|"
    r"(?:was|were)\s+(?:used|taken|obtained|ground|crushed|weighed)", re.I)
ORGANS = [
    ("flower", re.compile(r"(?<![A-Za-z])flowers?(?![A-Za-z])|blossoms?|"
                          r"(?<![A-Za-z])florets?(?![A-Za-z])", re.I)),
    ("petal", re.compile(r"(?<![A-Za-z])petals?(?![A-Za-z])", re.I)),
    ("bud", re.compile(r"(?<![A-Za-z])(?:flower\s+)?buds?(?![A-Za-z])", re.I)),
    ("leaf", re.compile(r"(?<![A-Za-z])leaf|leaves(?![A-Za-z])", re.I)),
    ("hip_fruit", re.compile(r"(?<![A-Za-z])(?:rose\s*)?hips?(?![A-Za-z])|"
                             r"(?<![A-Za-z])fruits?(?![A-Za-z])", re.I)),
    ("seed", re.compile(r"(?<![A-Za-z])seeds?(?![A-Za-z])", re.I)),
    ("stem", re.compile(r"(?<![A-Za-z])stems?(?![A-Za-z])|"
                        r"(?<![A-Za-z])shoots?(?![A-Za-z])", re.I)),
    ("root", re.compile(r"(?<![A-Za-z])roots?(?![A-Za-z])", re.I)),
    ("sepal", re.compile(r"(?<![A-Za-z])sepals?(?![A-Za-z])|calyx", re.I)),
    ("callus", re.compile(r"(?<![A-Za-z])callus(?:es)?(?![A-Za-z])", re.I)),
]


STAGE_PATTERNS = [
    ("bud", re.compile(r"(?<![A-Za-z])bud(?:ding)?\s+stage|unopened\s+(?:flower|bud)|"
                       r"closed\s+(?:flower|bud)", re.I)),
    ("half_open", re.compile(r"half[\s-]*open|semi[\s-]*open|partially\s+open|"
                             r"initial[\s-]*open", re.I)),
    ("full_bloom", re.compile(r"full(?:y)?[\s-]*(?:bloom|open|flower)|anthesis|"
                              r"blooming\s+stage", re.I)),
    ("senescent", re.compile(r"senescen\w*|withered|fading|post[\s-]*anthesis", re.I)),
    ("days_after_anthesis", re.compile(
        r"\d{1,3}\s*(?:days?|d)\s+after\s+(?:anthesis|flowering|full\s+bloom)", re.I)),
    ("named_phase", re.compile(
        r"(?:developmental|flowering|growth|maturity|ripening)\s+"
        r"(?:stages?|phases?)|stages?\s+(?:S\s*_?\s*1|I(?:\s*[-–,]\s*V)?)|"
        r"phase\s+(?:IV|V|III|II)(?![A-Za-z])", re.I)),
]


def organ_stated(text):

    found, ev = [], ""
    for name, pat in ORGANS:
        hit = _hit(pat, text, _MATERIAL_CTX, 160)
        if hit:
            found.append(name)
            ev = ev or hit
    return found, ev


def stage_stated(text):

    found, ev = [], ""
    for name, pat in STAGE_PATTERNS:
        hit = _hit(pat, text, _MATERIAL_CTX, 200)
        if hit:
            found.append(name)
            ev = ev or hit
    return found, ev


CULTIVAR_JUNK = re.compile(
    r"^(?:peak\s*area|values?|conc\w*|product\s+on\s+mass|content|amount|"
    r"mean|median|total|sum|ratio|percent\w*|relative|area|concentration|"
    r"mass|weight|yield|sd|se|std\w*|rsd|n\.?d\.?|trace|sample|samples|"
    r"control|blank|standard|std|rt|ri|no\.?|number|index|code|"
    r"[\d.,\s%±/-]*)$", re.I)


CULTIVAR_NOT_A_NAME = re.compile(
    r"(?<![A-Za-z])(?:μ|µ|u|m|n|k)?(?:g|L|l|mol)\s*[/·]|"
    r"(?<![A-Za-z])(?:conc|area|mass|weight|yield|content|ratio|mean|median|"
    r"stdev|st\.?\s*dev|rsd|sd|se|total|sum|value|percent|product\s+on)"
    r"(?![A-Za-z])|[%±]", re.I)
_CULT_TRAILER = re.compile(
    r"\s*[(,].*$|\s*\bcv\.?\b\s*|\s*\bvar\.?\b\s*", re.I)
_CULT_SPECIES = re.compile(
    r"(?<![A-Za-z])(?:rosa|r\.)\s*(?:x\s*)?[a-z]+(?![A-Za-z])", re.I)


def normalize_cultivar(value):

    v = (value or "").strip()
    if not v or CULTIVAR_JUNK.match(v) or CULTIVAR_NOT_A_NAME.search(v):
        return ""
    v = _CULT_TRAILER.sub("", v)
    v = _CULT_SPECIES.sub("", v)
    v = re.sub(r"[‘’'\"“”]", "", v)
    v = re.sub(r"\s+", " ", v).strip(" -,;:.'\u2018\u2019")
    if len(v) < 2 or CULTIVAR_JUNK.match(v):
        return ""
    return v.lower()


ITEMS = [

    ("platform_detector", "Analytical platform and detector named",
     "instrument", "essential", "cost:platform_detector"),
    ("column_phase", "Column stationary phase named",
     "instrument", "essential", "cost:column_polarity"),
    ("column_dimensions", "Column dimensions (length x i.d. x film)",
     "instrument", "recommended", ""),
    ("carrier_gas", "Carrier gas named",
     "instrument", "recommended", ""),
    ("oven_program", "Oven temperature programme given",
     "instrument", "recommended", ""),
    ("injector_temp", "Injector temperature given",
     "instrument", "desirable", ""),
    ("split_mode", "Split ratio or splitless stated",
     "instrument", "desirable", ""),

    ("ri_reported", "Retention indices reported at all",
     "identification", "essential", "cost:ri_coverage"),
    ("ri_basis_named", "RI basis named (Kovats / linear / arithmetic)",
     "identification", "essential", "cost:ri_basis"),
    ("ri_alkane_series", "n-alkane series named",
     "identification", "recommended", ""),
    ("ident_basis", "Identification basis stated (library / RI / standard)",
     "identification", "essential", "cost:ident_basis"),

    ("quant_mode", "Quantification mode stated",
     "quantification", "essential", "cost:quant_mode"),
    ("internal_standard", "Internal standard named",
     "quantification", "recommended", "cost:internal_standard"),
    ("response_factor", "Response-factor correction applied or declined",
     "quantification", "essential", "cost:response_factor"),

    ("species", "Species stated",
     "provenance", "essential", "cost:species"),
    ("cultivar", "Cultivar / accession stated",
     "provenance", "recommended", "yield:cultivar_bridges"),
    ("organ", "Plant organ stated",
     "provenance", "recommended", ""),
    ("stage", "Developmental stage stated",
     "provenance", "recommended", "cost:stage"),
    ("harvest_year", "Harvest year given",
     "provenance", "recommended", ""),
    ("harvest_month", "Harvest month given",
     "provenance", "recommended", ""),
    ("harvest_day", "Harvest day given",
     "provenance", "desirable", ""),
    ("harvest_time", "Harvest time of day given",
     "provenance", "desirable", ""),
    ("origin_country", "Country of origin given",
     "provenance", "recommended", "cost:origin_country"),
    ("origin_region", "Growing region given",
     "provenance", "desirable", ""),
    ("origin_gps", "Coordinates given",
     "provenance", "desirable", ""),

    ("material_type", "Material / product form stated (oil, absolute, "
     "hydrosol, extract, headspace)",
     "processing", "essential", "cost:material_type"),
    ("extraction_method", "Extraction / sampling method named",
     "processing", "essential", "cost:extraction_method"),
    ("extraction_time", "Extraction time given",
     "processing", "recommended", ""),
    ("extraction_temp", "Extraction temperature given",
     "processing", "recommended", ""),
    ("extraction_ratio", "Plant-to-liquid ratio given",
     "processing", "desirable", ""),

    ("replicates_strict", "Replication of the chemical measurement stated",
     "uncertainty", "recommended", "cost:dispersion"),
    ("dispersion_reported", "Dispersion attached to the reported values",
     "uncertainty", "essential", "cost:dispersion"),
]

ITEM_KEYS = [k for k, _l, _g, _t, _j in ITEMS]
ITEM_LABEL = {k: l for k, l, _g, _t, _j in ITEMS}
ITEM_GROUP = {k: g for k, _l, g, _t, _j in ITEMS}
ITEM_TIER = {k: t for k, _l, _g, t, _j in ITEMS}
ITEM_JUSTIFIED_BY = {k: j for k, _l, _g, _t, j in ITEMS}


POOLABILITY_ITEMS = ["platform_detector", "column_phase", "quant_mode",
                     "response_factor", "ri_basis_named", "species",
                     "material_type", "extraction_method"]
