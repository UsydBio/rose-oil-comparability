import math
from collections import Counter, defaultdict

from . import stats

__all__ = [
    "design_separability", "restricted_conditional", "bridging_levels",
    "log_replace", "SEPARABILITY_FIELDS",
]

SEPARABILITY_FIELDS = [
    "matrix", "block", "factor", "class", "n", "n_levels", "n_levels_named",
    "sep_code",
    "n_studies", "n_studies_multilevel", "n_samples_in_multilevel",
    "n_levels_in_ge2_studies", "n_samples_in_bridged_levels",
    "cramers_v_vs_study", "separable", "separability_verdict",
    "restricted_n", "restricted_n_studies", "restricted_df_factor",
    "restricted_r2_partial", "restricted_r2_of_within_study",
    "restricted_n_cells", "restricted_n_cells_ge2", "restricted_df_residual",
    "restricted_F", "restricted_p", "restricted_testable", "restricted_reason",
]


def design_separability(strata, labels, unstated="(unstated)"):

    strata = list(strata)
    labels = list(labels)
    n = len(labels)
    by_stratum = defaultdict(list)
    for i, s in enumerate(strata):
        by_stratum[s].append(i)

    multilevel = []
    for s, idx in by_stratum.items():
        if len({labels[i] for i in idx}) >= 2:
            multilevel.append(s)
    n_in_multi = sum(len(by_stratum[s]) for s in multilevel)

    by_level = defaultdict(set)
    level_n = Counter(labels)
    for i, lab in enumerate(labels):
        by_level[lab].add(strata[i])
    bridged = [lab for lab, ss in by_level.items()
               if lab != unstated and len(ss) >= 2]
    n_bridged_samples = sum(level_n[lab] for lab in bridged)

    named = [lab for lab in by_level if lab != unstated]
    v = stats.cramers_v(labels, strata) if n else float("nan")

    if len(by_level) < 2:
        code = "no-contrast"
        verdict = "only 1 level -- no contrast; R^2 cannot be reported"
    elif multilevel and bridged:
        code = "within+bridge"
        verdict = ("separable within study: %d studies contain a within-study "
                   "contrast (%d samples), and %d levels occur in >= 2 studies "
                   "(both routes are available)"
                   % (len(multilevel), n_in_multi, len(bridged)))
    elif multilevel:
        code = "within-only"
        verdict = ("separable within study only: %d studies contain a "
                   "within-study contrast (%d samples); no level is repeated "
                   "across studies, so cross-laboratory bridging does not hold "
                   "and the conclusion cannot be extrapolated to other "
                   "laboratories" % (len(multilevel), n_in_multi))
    elif bridged:
        code = "study-level-only"
        verdict = ("study-level only, nested within study: %d levels occur "
                   "in >= 2 studies (%d samples), but no single study contains "
                   ">= 2 levels: the factor can only be tested with study as the "
                   "unit of replication, and it is indistinguishable from every "
                   "other difference between those studies. "
                   "This is not separability." % (len(bridged), n_bridged_samples))
    else:
        code = "structurally-inseparable"
        verdict = ("structurally inseparable: no single study contains "
                   ">= 2 levels, and no level occurs in >= 2 studies; the factor "
                   "encodes the same partition as study")

    sep = bool(multilevel)

    return {
        "n": n,
        "n_levels": len(by_level),
        "n_levels_named": len(named),
        "n_studies": len(by_stratum),
        "n_studies_multilevel": len(multilevel),
        "multilevel_studies": sorted(multilevel),
        "n_samples_in_multilevel": n_in_multi,
        "n_levels_in_ge2_studies": len(bridged),
        "bridged_levels": sorted(bridged),
        "n_samples_in_bridged_levels": n_bridged_samples,
        "cramers_v_vs_study": v,
        "separable": sep,
        "sep_code": code,
        "separability_verdict": verdict,
    }


def restricted_conditional(d2, strata, labels, n_perm=9999, seed=20260912):

    strata = list(strata)
    labels = list(labels)
    by_stratum = defaultdict(list)
    for i, s in enumerate(strata):
        by_stratum[s].append(i)
    keep = []
    for s, idx in by_stratum.items():
        if len({labels[i] for i in idx}) >= 2:
            keep.extend(idx)
    keep.sort()
    if not keep:
        return {"testable": False, "n": 0, "n_studies": 0,
                "reason": "no single study contains >= 2 levels -- "
                          "structurally inseparable",
                "r2_partial": None, "r2_of_within_study": None,
                "F": None, "p": None, "df_factor": 0}
    sub = [[d2[i][j] for j in keep] for i in keep]
    res = stats.permanova_conditional(sub, [strata[i] for i in keep],
                                      [labels[i] for i in keep], n_perm, seed)
    res["n"] = len(keep)
    res["n_studies"] = len({strata[i] for i in keep})

    cellc = Counter((strata[i], labels[i]) for i in keep)
    res["n_cells"] = len(cellc)
    res["n_cells_ge2"] = sum(1 for v in cellc.values() if v >= 2)
    return res


def bridging_levels(strata, labels, unstated="(unstated)"):

    by_level = defaultdict(set)
    for s, lab in zip(strata, labels):
        if lab != unstated:
            by_level[lab].add(s)
    out = [(lab, sorted(ss)) for lab, ss in by_level.items() if len(ss) >= 2]
    out.sort(key=lambda x: (-len(x[1]), x[0]))
    return out


def log_replace(vec, deltas, base=10.0):

    if len(vec) != len(deltas):
        raise ValueError("log_replace: length mismatch")
    lb = math.log(base)
    out = []
    for v, d in zip(vec, deltas):
        x = v if v > 0 else d
        if x <= 0:
            raise ValueError("log_replace: still non-positive after replacement")
        out.append(math.log(x) / lb)
    return out
