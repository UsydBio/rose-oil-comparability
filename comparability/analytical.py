import re

from . import tables


_ANALYTICAL_TITLE = re.compile(
    r"(material|method|experimental|instrument|apparatus|equipment|"
    r"chromatograph|condition|analys[ie]s|analytical|determination|"
    r"identification|quantif|volatile|aroma|essential\s*oil|gc|ms\b|"
    r"sample\s*(collection|preparation)|reagent)", re.I)


def analytical_text(root):

    chunks = []

    def walk(node):

        for sec in node.findall("sec"):
            title = sec.find("title")
            label = tables.cell_text(title) if title is not None else ""
            if _ANALYTICAL_TITLE.search(label):
                chunks.append(tables.cell_text(sec))
            else:
                walk(sec)

    body = root.find(".//body")
    walk(body if body is not None else root)
    if chunks:
        return "\n".join(chunks), "methods_sections"

    if body is None:
        return "", ""
    parts = [tables.cell_text(sec) for sec in body.findall("sec")]
    if not parts:
        parts = [tables.cell_text(body)]
    return "\n".join(parts), "body_fallback"


def body_text(root):

    body = root.find(".//body")
    if body is None:
        return ""
    for ref in list(body.iter("ref-list")):
        for child in list(ref):
            ref.remove(child)
    return tables.cell_text(body)


def table_text(root):

    parts = []
    for tw in tables.iter_table_wraps(root):
        parts.append(" ".join([tw["label"], tw["caption"], tw["footer"]]))
        parts.extend(tw["headers"])
    return "\n".join(parts)


PLATFORMS = [

    ("GCxGC", "detection", re.compile(
        r"(?<![A-Za-z])gc\s*[×x×*]\s*gc(?![A-Za-z])|"
        r"comprehensive\s+two[\s-]*dimensional\s+gas\s+chromatograph", re.I)),

    ("GC-MS", "detection", re.compile(
        r"(?<![A-Za-z])gc\s*[-–—/−]?\s*"
        r"(?:q?[-–—/]?(?:tof|qq?q|it|q|sq|hr)\s*[-–—/]?\s*)?"
        r"ms(?![A-Za-z])|"
        r"(?<![A-Za-z])gc\s*[-–—/]\s*ms\s*/\s*ms|"
        r"gas\s+chromatograph\w*[\s-]*(?:–|-|/|and\s+)?\s*"
        r"(?:time[\s-]*of[\s-]*flight\s+)?mass\s+spectro", re.I)),
    ("GC-FID", "detection", re.compile(
        r"(?<![A-Za-z])gc\s*[-–—/]?\s*fid(?![A-Za-z])|"
        r"flame[\s-]*ionization\s+detect", re.I)),
    ("GC-O", "detection", re.compile(
        r"(?<![A-Za-z])gc\s*[-–—/]\s*o(?![A-Za-z])|"
        r"olfactometr|sniffing\s*port|gas\s+chromatography[\s-]*olfactometr", re.I)),
    ("LC-MS", "detection", re.compile(
        r"(?<![A-Za-z])(?:u?h?plc|lc)\s*[-–—/]?\s*ms(?![A-Za-z])|"
        r"liquid\s+chromatograph\w*[\s-]*(?:–|-|/|and\s+)?\s*mass\s+spectro|"
        r"(?<![A-Za-z])q[\s-]?tof(?![A-Za-z])", re.I)),
    ("LC", "separation", re.compile(
        r"(?<![A-Za-z])(uplc|uhplc|hplc)(?![A-Za-z])|"
        r"(ultra[\s-]*)?high[\s-]*performance\s+liquid\s+chromatograph", re.I)),
    ("PTR-MS", "detection", re.compile(
        r"(?<![A-Za-z])ptr\s*[-–]?\s*(?:tof\s*[-–]?\s*)?ms(?![A-Za-z])|"
        r"proton[\s-]*transfer[\s-]*reaction", re.I)),
    ("e-nose", "detection", re.compile(
        r"electronic\s+nose|(?<![A-Za-z])e\s*[-–]\s*nose(?![A-Za-z])", re.I)),
    ("NMR", "detection", re.compile(
        r"(?<![A-Za-z])nmr(?![A-Za-z])|nuclear\s+magnetic\s+resonance", re.I)),
    ("FTIR", "detection", re.compile(
        r"(?<![A-Za-z])ft\s*[-–]?\s*ir(?![A-Za-z])|"
        r"fourier[\s-]*transform\s+infrared", re.I)),

    ("HS-SPME", "introduction", re.compile(
        r"(?<![A-Za-z])hs\s*[-–—/]?\s*spme(?![A-Za-z])|"
        r"head\s*-?\s*space\s+solid[\s-]*phase\s+micro\s*-?extraction", re.I)),
    ("SPME", "introduction", re.compile(
        r"(?<![A-Za-z])spme(?![A-Za-z])|solid[\s-]*phase\s+micro\s*-?extraction", re.I)),
    ("headspace", "introduction", re.compile(
        r"head\s*-?\s*space", re.I)),
    ("SAFE", "introduction", re.compile(

        r"(?<![A-Za-z])SAFE(?![A-Za-z])|"
        r"(?i:solvent[\s-]*assisted\s+flavou?r\s+evaporation)")),
    ("SDE", "introduction", re.compile(
        r"(?<![A-Za-z])sde(?![A-Za-z])|simultaneous\s+distillation[\s-]*extraction", re.I)),
    ("hydrodistillation", "introduction", re.compile(
        r"hydro\s*-?\s*distillation|steam\s+distillation|Clevenger", re.I)),
    ("TD", "introduction", re.compile(
        r"thermal\s+desorption|(?<![A-Za-z])td\s*[-–]\s*gc(?![A-Za-z])", re.I)),
]

_PLATFORM_BY_NAME = {name: (axis, pat) for name, axis, pat in PLATFORMS}


_PLATFORM_IMPLIES = {"HS-SPME": ("SPME", "headspace")}


def _platforms_in(text):

    hits = [name for name, _axis, pat in PLATFORMS if pat.search(text)]
    drop = set()
    for name in hits:
        drop.update(_PLATFORM_IMPLIES.get(name, ()))
    return [h for h in hits if h not in drop]


_COLUMN_STRONG = (r"DB|HP|ZB|Rtx|RTX|Innowax|InertCap|Stabilwax|Supelcowax|"
                  r"Carbowax|EconoCap|BPX|SPB|TRB|SolGel|Restek")


_COLUMN_WEAK = (r"CP|TR|TG|BP|VF|SE|OV|EC|PE|AT|MEGA|NST|Elite|Optima|"
                r"Equity|CBP|TC|MDN")

_COLUMN_STRONG_RE = re.compile(
    r"(?<![A-Za-z])(?P<fam>" + _COLUMN_STRONG + r")(?:TM|™|®)?[\s‐‑‒–-]?"
    r"(?P<spec>(?:INNO|Inno)?WAX|(?:INNO|Inno)?Wax|wax|\d{1,3})(?!\d)"
    r"(?P<suffix>MS|ms|Ms|HT|ht|M(?![A-Za-z])|Sil|sil)?",
    re.I)
_COLUMN_WEAK_RE = re.compile(
    r"(?<![A-Za-z])(?P<fam>" + _COLUMN_WEAK + r")(?:TM|™|®)?[‐‑‒–-]"
    r"(?P<spec>WAX|Wax|\d{1,3})(?!\d)"
    r"(?P<suffix>MS|ms|Ms|HT|Sil)?")


_FAMILY_CASE = {
    "db": "DB", "hp": "HP", "zb": "ZB", "rtx": "Rtx", "cp": "CP", "tr": "TR",
    "tg": "TG", "bp": "BP", "bpx": "BPX", "vf": "VF", "spb": "SPB",
    "se": "SE", "ov": "OV", "ec": "EC", "pe": "PE", "at": "AT",
    "mega": "MEGA", "nst": "NST", "cbp": "CBP", "tc": "TC", "mdn": "MDN",
    "elite": "Elite", "optima": "Optima", "equity": "Equity",
    "innowax": "Innowax", "inertcap": "InertCap", "stabilwax": "Stabilwax",
    "supelcowax": "Supelcowax", "carbowax": "Carbowax",
    "econocap": "EconoCap", "trb": "TRB", "solgel": "SolGel",
    "restek": "Restek",
}


_COLUMN_CONTEXT = re.compile(
    r"column|capillar|stationary\s+phase|chromatograph|(?<![A-Za-z])GC(?![A-Za-z])", re.I)


_COLUMN_ANTICONTEXT = re.compile(
    r"resin|adsorb|adsorp|desorption\s+curve|silica\s+gel|sephadex|"
    r"macroporous|ion[\s-]*exchange|flash\s+chromatograph|"
    r"solid[\s-]*phase\s+extraction|vacuum\s+liquid", re.I)


_COLUMN_PROSE = [
    ("PEG/wax", re.compile(r"polyethylene\s*glycol|(?<![A-Za-z])PEG(?![A-Za-z])|"
                           r"(?<![A-Za-z])FFAP(?![A-Za-z])", re.I)),
    ("5%-phenyl-methylpolysiloxane", re.compile(
        r"5\s*%\s*(?:-)?\s*phenyl[\s-]*(?:95\s*%\s*(?:-)?\s*)?"
        r"(?:methyl|dimethyl)[\s-]*(?:poly)?siloxane", re.I)),
    ("50%-phenyl-methylpolysiloxane", re.compile(
        r"(?:50|35)\s*%\s*(?:-)?\s*phenyl", re.I)),
    ("cyclodextrin", re.compile(r"cyclodextrin|chirasil|lipodex", re.I)),
]


_POLARITY_RULES = [
    ("polar", re.compile(r"wax|peg|ffap|20M|innowax|carbowax|stabilwax|supelcowax", re.I)),
    ("chiral", re.compile(r"cyclodextrin|BGB|chirasil|lipodex", re.I)),
    ("medium-polar", re.compile(r"(?<![A-Za-z\d])(17|35|1301|624|1701|225)(?!\d)", re.I)),
    ("non-polar", re.compile(r"(?<![A-Za-z\d])(1|5|8|30|54|101)(?!\d)", re.I)),
]


_POLARITY_WORD = (r"non[\s-]*polar|apolar|weakly\s+polar|"
                  r"(?:medium|mid|middle|semi)[\s-]*polar|"
                  r"strong\s+polarity|(?<![a-z-])polar")
_POLARITY_ASSERTED = re.compile(
    r"(" + _POLARITY_WORD + r")\s*"
    r"(?:capillary\s+)?(?:\d?D?\s*)?(?:column|phase|polarity)|"
    r"column(?:[^.]|\.\d){0,40}?(" + _POLARITY_WORD + r")",
    re.I)
_POLARITY_NORM = {
    "nonpolar": "non-polar", "non polar": "non-polar", "apolar": "non-polar",
    "weakly polar": "medium-polar", "medium polar": "medium-polar",
    "mediumpolar": "medium-polar", "mid polar": "medium-polar",
    "midpolar": "medium-polar", "middle polar": "medium-polar",
    "semi polar": "medium-polar", "semipolar": "medium-polar",
    "strong polarity": "polar",
}


_DIM_RE = re.compile(
    r"(?P<len>\d{1,3}(?:\.\d+)?)\s*m\b[^,;()]{0,18}?"
    r"[×x×*£X]\s*"
    r"(?P<id>\d{1,3}(?:\.\d+)?)\s*(?P<idu>mm|μm|µm|um)\b"
    r"(?:[^,;()]{0,40}?[×x×*£X,;]\s*"
    r"(?P<film>\d{1,2}(?:\.\d+)?)\s*(?P<filmu>μm|µm|um|mm)?)?",
    re.I)


def _near(text, start, end, pattern, span=200):
    return bool(pattern.search(text[max(0, start - span):end + span]))


def _column_names(text):

    found = []

    def add(name):
        if name not in found:
            found.append(name)

    for pat in (_COLUMN_STRONG_RE, _COLUMN_WEAK_RE):
        for m in pat.finditer(text):
            if not _near(text, m.start(), m.end(), _COLUMN_CONTEXT):
                continue
            if _near(text, m.start(), m.end(), _COLUMN_ANTICONTEXT, 80):
                continue
            fam = _FAMILY_CASE.get(m.group("fam").lower(), m.group("fam"))
            spec = m.group("spec")
            suffix = m.group("suffix") or ""
            if spec.lower().endswith("wax"):
                spec = "INNOWax" if spec.lower().startswith("inno") else "Wax"
                suffix = ""
            add("%s-%s%s" % (fam, spec, suffix.upper() if suffix else ""))

    for m in re.finditer(r"(?<![A-Za-z])BGB[\s-]?(\d{2,3})(?!\d)", text):
        if _near(text, m.start(), m.end(), _COLUMN_CONTEXT):
            add("BGB-%s" % m.group(1))
    found = [n for n in found if not _COLUMN_ANTICONTEXT.search(n)]

    for label, pat in _COLUMN_PROSE:
        for m in pat.finditer(text):
            if _near(text, m.start(), m.end(), _COLUMN_CONTEXT, 120):
                add(label)
                break
    return found


def _polarity_of(names):

    labels = []
    for name in names:
        got = ""
        for label, pat in _POLARITY_RULES:
            if pat.search(name):
                got = label
                break
        if got and got not in labels:
            labels.append(got)
    if not labels:
        return "", ""
    if len(labels) == 1:
        return labels[0], ";".join(labels)
    return "mixed", ";".join(labels)


def _dimensions(text):

    for m in _DIM_RE.finditer(text):
        length = float(m.group("len"))

        if not (5 <= length <= 120):
            continue
        idv = float(m.group("id"))
        idu = m.group("idu").lower()
        if idu in ("μm", "µm", "um"):
            idv = idv / 1000.0
        elif idv > 5:
            idv = idv / 1000.0
        if not (0.05 <= idv <= 1.0):
            continue
        film = ""
        if m.group("film"):
            fv = float(m.group("film"))
            fu = (m.group("filmu") or "").lower()
            if fu == "mm":
                fv = fv * 1000.0
            if 0.05 <= fv <= 60:
                film = "%g" % fv
        return m.group(0).strip(), "%g" % length, "%g" % idv, film
    return "", "", "", ""


_RI_BASIS = [
    ("kovats", re.compile(r"kov[aá]ts|kovat", re.I)),
    ("linear", re.compile(r"linear\s+retention\s+ind|"
                          r"(?<![A-Za-z])lri(?![A-Za-z])|"
                          r"van\s+den\s+dool", re.I)),
    ("arithmetic", re.compile(r"arithmetic\s+(retention\s+)?ind|"
                              r"adams.{0,30}arithmetic", re.I)),
]
_RI_ANY = re.compile(
    r"retention\s+ind|(?<![A-Za-z])(ri|lri|kri|rri)\s*(?:value|s)?(?![A-Za-z])", re.I)


_ALKANE_RE = re.compile(
    r"(?<![A-Za-z])C\s*_?\s*(\d{1,2})\s*"
    r"(?:[-‐‑‒–—−~]|to)\s*"
    r"C\s*_?\s*(\d{1,2})(?!\d)", re.I)


_ALKANE_CONTEXT = re.compile(
    r"n\s*-?\s*alkane|paraffin|homologous\s+series|"
    r"retention\s+ind\w*\s*(?:\(?(?:ri|lri|rri|ki)\)?)?[^A-Za-z]{0,20}"
    r"(?:was|were|are|is|calculat|determin|comput|based)|"
    r"calculated\s+(?:against|using|from)|kov[aá]ts", re.I)


def _alkane_series(text):
    for m in _ALKANE_RE.finditer(text):
        window = text[max(0, m.start() - 200):m.end() + 200]
        if not _ALKANE_CONTEXT.search(window):
            continue
        lo, hi = int(m.group(1)), int(m.group(2))
        if not (3 <= lo < hi <= 44):
            continue
        return "C%d-C%d" % (lo, hi)
    return ""


_QUANT_MODES = [
    ("relative_area_pct", re.compile(
        r"(peak[\s-]*)?area\s+normali[sz]|normali[sz]ation\s+method|"
        r"relative\s+(peak\s+)?area|percent(age)?\s+(relative\s+)?peak\s+area|"
        r"area\s+percent|internal\s+normali[sz]ation|"
        r"without\s+(any\s+)?(the\s+use\s+of\s+)?(response\s+factor\s+)?correction", re.I)),
    ("internal_standard", re.compile(r"internal\s+standard", re.I)),
    ("external_calibration", re.compile(
        r"external\s+standard|external\s+calibration|"
        r"(standard|calibration)\s+curve", re.I)),
    ("absolute_conc", re.compile(
        r"(mass\s+)?concentration\s+of\s+each|"
        r"content\s+of\s+each\s+(component|compound)|"
        r"(μ|µ|m|n)g\s*[·/∙]\s*(g|kg|l|ml|m?L)(?![A-Za-z])", re.I)),
]


_CHEM = (r"(?:[\dA-Za-z()\[\],'’αβγ−–.+-]+\s+)?"
         r"[\dA-Za-z(][\dA-Za-z()\[\],'’αβγ−–.-]{3,}")


_IS_NAME_PATTERNS = [
    re.compile(r"(" + _CHEM + r")\s*\((?:the\s+|an?\s+)?internal\s+standard", re.I),
    re.compile(r"(?:using|with|of|added|adding|to)\s+(" + _CHEM + r")\s+as\s+"
               r"(?:the|an|a)\s+internal\s+standard", re.I),
    re.compile(r"(" + _CHEM + r")\s+(?:was|were)\s+(?:used\s+)?as\s+"
               r"(?:the|an|a)\s+internal\s+standard", re.I),
    re.compile(r"(" + _CHEM + r")\s+internal\s+standard\s+solution", re.I),

    re.compile(r"internal\s+standard\s*[(:,]\s*(" + _CHEM + r")", re.I),
    re.compile(r"(?:with|using)\s+internal\s+standard\s+(" + _CHEM + r")", re.I),
]


_IS_STOP = re.compile(
    r"(?<![A-Za-z])(the|and|or|with|using|used|was|were|as|to|for|in|on|by|"
    r"that|which|this|internal|external|solution|standard|standards|method|"
    r"sample|samples|area|peak|concentration|ratio|volume|each|its|their|"
    r"chromatogram|mixture|added|respectively|of|a|an|at|from|our|"
    r"procedure|protocol|approach|previous|described|above|below|"
    r"contain\w*)(?![A-Za-z])", re.I)


_IS_LEAD_NOISE = re.compile(
    r"^(?:[\d.,]+\s*)?(?:μ|µ|u|m|n)?[LlgG](?![A-Za-z-])\s*(?:of\s+)?|"
    r"^[\d.,]+\s+(?:of\s+)?|^of\s+", re.I)


def _internal_standard(text):

    names = []
    for pat in _IS_NAME_PATTERNS:
        for m in pat.finditer(text):
            cand = m.group(1).strip(" ,;:.")
            for _ in range(3):
                new = _IS_LEAD_NOISE.sub("", cand).strip(" ,;:.")
                if new == cand:
                    break
                cand = new

            words = cand.split()
            while len(words) > 1 and not re.search(r"[A-Za-z]{3}", words[-1]):
                words.pop()
            cand = " ".join(words).strip(" ,;:.([")
            if len(cand) < 5 or _IS_STOP.search(cand):
                continue
            if not re.search(r"[A-Za-z]{4}", cand):
                continue
            if cand.lower() not in [n.lower() for n in names]:
                names.append(cand)
    return names


_CARRIER_MARK = re.compile(
    r"carrier\s*gas|as\s+(?:a|the)\s+carrier|carrier\s+flow", re.I)
_GAS_RE = re.compile(
    r"(?i:helium|hydrogen|nitrogen|argon)|"
    r"(?<![A-Za-z])(?:He|H2|N2|Ar)(?![A-Za-z0-9])")
_CARRIER_NORM = {"he": "helium", "h2": "hydrogen", "n2": "nitrogen",
                 "ar": "argon"}


def _carrier_gas(text):

    for mark in _CARRIER_MARK.finditer(text):
        lo, hi = max(0, mark.start() - 200), mark.end() + 200
        window = text[lo:hi]
        best, bestd = "", 10 ** 9
        for g in _GAS_RE.finditer(window):
            d = abs((lo + g.start()) - mark.start())
            if d < bestd:
                best, bestd = g.group(0), d
        if best:
            return _CARRIER_NORM.get(best.lower(), best.lower())
    return ""


_FLOW_RE = re.compile(
    r"(?:flow(?:\s+rate)?|velocity)[^.]{0,40}?(\d{1,2}(?:\.\d+)?)\s*"
    r"(?:mL|ml|cm)\s*[·/\s]?\s*(?:min|s)(?![A-Za-z])|"
    r"(\d{1,2}(?:\.\d+)?)\s*(?:mL|ml)\s*[·/]\s*min(?![A-Za-z])", re.I)
_FLOW_CONTEXT = re.compile(
    r"carrier|column\s+flow|constant\s+(flow|pressure)|linear\s+velocity", re.I)

_INJ_TEMP_RE = re.compile(
    r"(?:injector|inlet|injection\s+port)(?:[^.]|\.\d){0,70}?"
    r"(\d{2,3})\s*[°℃℉]?\s*C(?![a-z])", re.I)

_SPLIT_RATIO_RE = re.compile(
    r"split(?:\s+(?:ratio|mode))?(?:[^.]|\.\d){0,30}?"
    r"(?:\(\s*)?(\d{1,3}\s*[:：]\s*\d{1,3})", re.I)
_SPLITLESS_RE = re.compile(r"split\s*-?\s*less|non\s*-?\s*split|without\s+split", re.I)

_OVEN_RE = re.compile(
    r"[^.]{0,40}(?:oven|column)\s+temperature[^.]{0,400}\.|"
    r"[^.]{0,40}temperature\s+program[^.]{0,400}\.|"
    r"[^.]{0,40}(?:initially|initial)\s+(?:held|maintained|set)\s+at\s+\d{1,3}\s*"
    r"[°℃]?\s*C[^.]{0,400}\.", re.I)

_EI_RE = re.compile(r"(\d{2,3})\s*eV(?![A-Za-z])|electron\s+(impact|ioni[sz]ation)", re.I)


_LIB_PATTERNS = [
    ("NIST", re.compile(r"(?<![A-Za-z])NIST\s*[-.]?\s*(\d{1,2}(?:\.\d)?)?(?![\d.])")),

    ("NIST", re.compile(r"National\s+Institute\s+of\s+Standards\s+and\s+Technology"
                        r"[^.]{0,80}?\(?\s*(?:Version|v\.?)\s*(\d{1,2}(?:\.\d)?)", re.I)),
    ("NIST", re.compile(r"National\s+Institute\s+of\s+Standards\s+and\s+Technology()", re.I)),
    ("Wiley", re.compile(r"(?<![A-Za-z])Wiley\s*[-]?\s*(\d{1,2})?(?:th|n\.?l)?", re.I)),
    ("FFNSC", re.compile(r"(?<![A-Za-z])FFNSC\s*(\d(?:\.\d)?)?", re.I)),
    ("Adams", re.compile(r"(?<![A-Za-z])Adams(?![A-Za-z])\s*(\d{4})?")),
    ("MassFinder", re.compile(r"MassFinder\s*(\d)?", re.I)),
]

_IDENT_BASIS = [
    ("ms_library", re.compile(
        r"(?<![A-Za-z])NIST(?![A-Za-z])|(?i:wiley|ffnsc|massfinder|"
        r"mass\s+spectr\w+\s+(library|database)|spectral\s+library|"
        r"match(?:ing)?\s+factor|similarity\s+(index|of|greater))")),
    ("ri_comparison", re.compile(
        r"retention\s+ind\w*(?:[^.]|\.\d){0,120}?(literature|reported|published|"
        r"compar|reference|database)|"
        r"(compar|match)\w*(?:[^.]|\.\d){0,60}?retention\s+ind", re.I)),
    ("authentic_standard", re.compile(
        r"authentic\s+(standard|compound|sample)|co[\s-]*injection|"
        r"reference\s+standard|standard\s+(compound|substance)s?\s+"
        r"(were|was)\s+(used|purchased|run)", re.I)),
    ("odor_description", re.compile(
        r"odou?r\s+descript|aroma\s+descript|sniffing", re.I)),
]


def _first(pattern, text, group=1):
    m = pattern.search(text)
    if not m:
        return ""
    if isinstance(group, int):
        return (m.group(group) or "").strip()
    for g in group:
        if m.group(g):
            return m.group(g).strip()
    return ""


def parse_analytical(root):

    text, text_source = analytical_text(root)
    whole = body_text(root)
    tabtext = table_text(root)

    out = {name: "" for name in FIELDS}
    out["study_analytical_level"] = "study"
    out["study_analytical_text_source"] = text_source
    out["study_analytical_chars"] = len(text)
    out["study_platform_n"] = 0
    out["study_column_n"] = 0
    if not text and not whole:
        return out

    in_methods = _platforms_in(text)
    in_body = _platforms_in(whole)
    out["study_platform"] = ";".join(in_methods)
    out["study_platform_n"] = len(in_methods)
    out["study_platform_source"] = text_source if in_methods else ""
    out["study_platform_mentioned_only"] = ";".join(
        p for p in in_body if p not in in_methods)
    out["study_platform_detection"] = ";".join(
        p for p in in_methods if _PLATFORM_BY_NAME[p][0] == "detection")
    out["study_platform_introduction"] = ";".join(
        p for p in in_methods if _PLATFORM_BY_NAME[p][0] == "introduction")

    in_tables = _platforms_in(tabtext) if tabtext else []

    tbl_det = [p for p in in_tables if _PLATFORM_BY_NAME[p][0] == "detection"]
    det = [p for p in in_methods if _PLATFORM_BY_NAME[p][0] == "detection"]
    quant_sentences = " ".join(
        m.group(0) for m in re.finditer(
            r"[^.]{0,200}(?:quantif|relative\s+content|peak\s+area\s+normali|"
            r"relative\s+(?:peak\s+)?area|percentage\s+of\s+each)[^.]{0,200}\.",
            text, re.I))
    quant_plat = [p for p in _platforms_in(quant_sentences)
                  if _PLATFORM_BY_NAME[p][0] == "detection"]
    if tbl_det:
        out["study_platform_primary"] = ";".join(tbl_det)
        out["study_platform_primary_source"] = "table_caption"
    elif len(quant_plat) == 1:
        out["study_platform_primary"] = quant_plat[0]
        out["study_platform_primary_source"] = "quantification_sentence"
    elif len(det) == 1:
        out["study_platform_primary"] = det[0]
        out["study_platform_primary_source"] = "sole_detection_platform"

    names = _column_names(text)
    if not names and text_source == "methods_sections":

        names = _column_names(whole)
        if names:
            out["study_column_source"] = "body_fallback"
    if names and not out["study_column_source"]:
        out["study_column_source"] = text_source
    out["study_column_phase"] = names[0] if names else ""
    out["study_column_phase_all"] = ";".join(names)
    out["study_column_n"] = len(names)

    polarity, detail = _polarity_of(names)
    asserted = []
    for m in _POLARITY_ASSERTED.finditer(text):
        word = (m.group(1) or m.group(2) or "").lower()
        word = re.sub(r"[\s-]+", " ", word).strip()
        word = _POLARITY_NORM.get(word, word)
        if word and word not in asserted:
            asserted.append(word)
    if asserted:
        out["study_column_polarity"] = (
            asserted[0] if len(asserted) == 1 else "mixed")
        out["study_column_polarity_source"] = "text_assertion"
        out["study_column_polarity_inferred"] = ";".join(asserted)
    elif polarity:
        out["study_column_polarity"] = polarity
        out["study_column_polarity_source"] = "inferred_from_phase_name"
        out["study_column_polarity_inferred"] = detail

    raw, length, idmm, film = _dimensions(text)
    if not raw:
        raw, length, idmm, film = _dimensions(whole)
        if raw:
            out["study_column_dim_source"] = "body_fallback"
    elif raw:
        out["study_column_dim_source"] = text_source
    out["study_column_dimensions"] = raw
    out["study_column_length_m"] = length
    out["study_column_id_mm"] = idmm
    out["study_column_film_um"] = film

    basis = [name for name, pat in _RI_BASIS if pat.search(text)]
    if basis:
        out["study_ri_basis"] = ";".join(basis)

        out["study_ri_basis_named"] = ";".join(basis)
        out["study_ri_basis_source"] = text_source
    elif _RI_ANY.search(text):
        out["study_ri_basis"] = "unspecified"
        out["study_ri_basis_source"] = text_source
    out["study_ri_alkane_series"] = _alkane_series(text) or _alkane_series(whole)
    if out["study_ri_alkane_series"]:
        out["study_ri_series_source"] = (
            text_source if _alkane_series(text) else "body_fallback")

    modes = [name for name, pat in _QUANT_MODES if pat.search(text)]
    out["study_quant_mode"] = ";".join(modes)
    out["study_quant_mode_n"] = len(modes)
    if modes:
        out["study_quant_mode_source"] = text_source
    istds = _internal_standard(text)
    out["study_internal_standard"] = ";".join(istds[:3])
    if istds:
        out["study_internal_standard_source"] = text_source

    out["study_carrier_gas"] = _carrier_gas(text)
    out["study_carrier_flow"] = ""
    for m in _FLOW_RE.finditer(text):
        if _near(text, m.start(), m.end(), _FLOW_CONTEXT, 120):
            out["study_carrier_flow"] = (m.group(1) or m.group(2) or "").strip()
            break
    out["study_injector_temp_c"] = _first(_INJ_TEMP_RE, text)
    ratio = _first(_SPLIT_RATIO_RE, text)
    if ratio:
        out["study_split_ratio"] = re.sub(r"\s+", "", ratio)
        out["study_split_mode"] = "split"
    elif _SPLITLESS_RE.search(text):
        out["study_split_mode"] = "splitless"
    oven = _first(_OVEN_RE, text, 0)
    out["study_oven_program"] = re.sub(r"\s+", " ", oven)[:300]
    ev = _first(_EI_RE, text)
    out["study_ms_ionization"] = ("EI %s eV" % ev) if ev else (
        "EI" if _EI_RE.search(text) else "")
    if any(out[f] for f in ("study_carrier_gas", "study_injector_temp_c",
                            "study_split_mode", "study_oven_program")):
        out["study_conditions_source"] = text_source

    libs = []
    for label, pat in _LIB_PATTERNS:
        for m in pat.finditer(text):
            ver = (m.group(1) or "").strip()
            entry = ("%s %s" % (label, ver)).strip()
            if entry not in libs:
                libs.append(entry)

    libs = [l for l in libs
            if " " in l or not any(o.startswith(l + " ") for o in libs)]
    out["study_ident_library"] = ";".join(libs[:5])

    if out["study_platform_primary"]:
        det_key = out["study_platform_primary"].replace(";", "+")
        out["study_platform_group_basis"] = "primary"
    elif det:
        det_key = "+".join(det)
        out["study_platform_group_basis"] = "detector_set"
    else:
        det_key = ""
        out["study_platform_group_basis"] = ""
    if det_key:
        out["study_platform_group"] = "%s|%s" % (
            det_key, out["study_column_polarity"] or "unknown-polarity")
    ident = [name for name, pat in _IDENT_BASIS if pat.search(text)]
    out["study_ident_basis"] = ";".join(ident)
    out["study_ident_basis_n"] = len(ident)
    if ident or libs:
        out["study_ident_source"] = text_source

    return out


FIELDS = [
    "study_analytical_level", "study_analytical_text_source",
    "study_analytical_chars",

    "study_platform", "study_platform_n", "study_platform_detection",
    "study_platform_introduction", "study_platform_primary",
    "study_platform_primary_source", "study_platform_mentioned_only",
    "study_platform_source", "study_platform_group",
    "study_platform_group_basis",

    "study_column_phase", "study_column_phase_all", "study_column_n",
    "study_column_polarity", "study_column_polarity_source",
    "study_column_polarity_inferred", "study_column_source",
    "study_column_dimensions", "study_column_length_m", "study_column_id_mm",
    "study_column_film_um", "study_column_dim_source",

    "study_ri_basis", "study_ri_basis_named", "study_ri_basis_source",
    "study_ri_alkane_series",
    "study_ri_series_source",

    "study_quant_mode", "study_quant_mode_n", "study_quant_mode_source",
    "study_internal_standard", "study_internal_standard_source",

    "study_carrier_gas", "study_carrier_flow", "study_injector_temp_c",
    "study_split_ratio", "study_split_mode", "study_oven_program",
    "study_ms_ionization", "study_conditions_source",

    "study_ident_library", "study_ident_basis", "study_ident_basis_n",
    "study_ident_source",
]
