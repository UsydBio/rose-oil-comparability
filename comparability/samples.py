import re


SPECIES = {
    "damascena": "Rosa x damascena",
    "damascene": "Rosa x damascena",
    "damascenea": "Rosa x damascena",
    "centifolia": "Rosa x centifolia",
    "gallica": "Rosa gallica",
    "galica": "Rosa gallica",
    "alba": "Rosa alba",
    "rugosa": "Rosa rugosa",
    "chinensis": "Rosa chinensis",
    "hybrida": "Rosa hybrida",
    "canina": "Rosa canina",
    "moschata": "Rosa moschata",
    "multiflora": "Rosa multiflora",
    "odorata": "Rosa odorata",
    "laevigata": "Rosa laevigata",
    "banksiae": "Rosa banksiae",
    "willmottiae": "Rosa willmottiae",
    "sertata": "Rosa sertata",
    "majalis": "Rosa majalis",
    "roxburghii": "Rosa roxburghii",
}
_SPECIES_RE = re.compile(
    r"(?:\bR(?:osa)?\.?\s*)?(?:×|x\s)?\s*\b(%s)\b" % "|".join(SPECIES), re.I)


NAMED_CULTIVARS = {
    "pingyin": "Pingyin",
    "kushui": "Kushui",
    "jinbian": "Jinbian",
    "mohong": "Mohong",
    "dianhong": "Dianhong",
    "kazanlak": "Kazanlak",
    "kazanlik": "Kazanlak",
    "isparta": "Isparta",
    "taif": "Taif",
}
_CULTIVAR_NAMED_RE = re.compile(r"\b(%s)\b" % "|".join(NAMED_CULTIVARS), re.I)

ORGANS = {
    "petal": "petal", "petals": "petal",
    "flower": "flower", "flowers": "flower", "floral": "flower",
    "bud": "bud", "buds": "bud",
    "fruit": "fruit", "fruits": "fruit", "hip": "hip", "hips": "hip",
    "leaf": "leaf", "leaves": "leaf",
    "stem": "stem", "stems": "stem", "root": "root", "roots": "root",
    "seed": "seed", "seeds": "seed", "pollen": "pollen",
    "calyx": "calyx", "stamen": "stamen",
}
_ORGAN_RE = re.compile(r"\b(%s)\b" % "|".join(ORGANS), re.I)

PROCESSING = [
    ("hydrodistillation", re.compile(r"\b(hydro\s*-?\s*distillation|\bhd\b)\b", re.I)),
    ("steam_distillation", re.compile(r"\bsteam\s*-?\s*distill", re.I)),
    ("supercritical_co2", re.compile(r"\b(sfe|supercritical|sc\s*-?\s*co2|co2\s*extract)\b", re.I)),
    ("solvent_extraction", re.compile(r"\b(solvent\s*extract|hexane|petroleum\s*ether)\b", re.I)),
    ("absolute", re.compile(r"\babsolute\b", re.I)),
    ("concrete", re.compile(r"\bconcrete\b", re.I)),
    ("hydrosol", re.compile(r"\b(hydrosol|rose\s*water|distillation\s*water)\b", re.I)),
    ("microwave", re.compile(r"\bmicrowave\b", re.I)),
    ("ultrasound", re.compile(r"\b(ultrasound|ultrasonic)\b", re.I)),
    ("headspace", re.compile(r"\b(headspace|hs\s*-?\s*spme|spme)\b", re.I)),
    ("dried", re.compile(r"\b(dried|drying|air\s*-?\s*dry|freeze\s*-?\s*dr)\b", re.I)),
    ("fresh", re.compile(r"\bfresh\b", re.I)),
]


SOLVENTS = {
    "hexane": "hexane", "n-hexane": "hexane",
    "dichloromethane": "dichloromethane", "dcm": "dichloromethane",
    "methanol": "methanol", "meoh": "methanol",
    "ethanol": "ethanol", "etoh": "ethanol",
    "ethyl acetate": "ethyl_acetate", "acetone": "acetone",
    "petroleum ether": "petroleum_ether", "diethyl ether": "diethyl_ether",
    "water": "water", "pentane": "pentane", "chloroform": "chloroform",
}
_SOLVENT_RE = re.compile(r"\b(%s)\b" % "|".join(sorted(SOLVENTS, key=len, reverse=True)), re.I)


_LABEL_BLACKLIST = re.compile(
    r"^(gc\s*-?\s*ms|gc\s*-?\s*fid|ms|fid|hplc|nmr|ctrl|control|blank|"
    r"n\s*a|na|nd|dilution\s*stock|stock|standard|std|sample|samples|"
    r"p\s*-?\s*value|range|mean|median|min|max|lit\.?|exp\.?|"
    r"code|serial\s*number|serial|no\.?|number|id|index|"
    r"r\s*\d?|rri.*|rt.*|t\s*_?\s*r.*)$", re.I)

PLATFORM = [
    ("GC-FID", re.compile(r"\bgc\s*-?\s*fid\b", re.I)),
    ("GC-MS", re.compile(r"\bgc\s*[-x×/]?\s*ms\b", re.I)),
    ("GCxGC", re.compile(r"\bgc\s*[x×]\s*gc\b", re.I)),
    ("LC-MS", re.compile(r"\b(lc|hplc|uplc|uhplc)\s*-?\s*ms\b", re.I)),
]


_MEASURE_NOISE = re.compile(
    r"(relative\s*content|rel\.?\s*%?,?\s*as\s*determined(\s*(by|using))?|"
    r"area\s*%\s*(in)?|content\s*%?\s*[-–]?|rate\s*\(\s*%\s*\)|"
    r"concentration|\bcontent\b|\brel\.?\b|\brelative\b|\bmean\b|\bsd\b|"
    r"\(\s*(%|mg\s*/\s*kg|µg\s*/\s*l|ug\s*/\s*l|mg\s*/\s*l|ppm|ppb)[^)]*\)|"
    r"\b(µg|ug|mg|ng)\s*[·/]\s*(kg|g|l|ml)\b|%)",
    re.I,
)


_CODE_RE = re.compile(r"^(?:[A-Z]{1,4}[-_]?\d{1,3}|[IVX]{1,5}|[A-Z]{1,3})$")


def species_in(text):

    if not text:
        return []
    seen = []
    for m in _SPECIES_RE.finditer(text):
        name = SPECIES[m.group(1).lower()]
        if name not in seen:
            seen.append(name)
    return seen



FIELDS = ["sample_species", "sample_cultivar", "sample_organ",
          "sample_processing", "sample_platform", "sample_label",
          "sample_parse_confidence",

          "sample_stage", "sample_origin", "sample_treatment",
          "sample_code_axis", "sample_code_source"]
