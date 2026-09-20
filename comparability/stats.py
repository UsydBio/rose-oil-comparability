import math
import random
from collections import Counter, defaultdict

__all__ = [
    "closure", "multiplicative_replacement", "clr",
    "euclidean", "bray_curtis", "jaccard",
    "dist_matrix", "sq_dist_matrix",
    "ss_total", "ss_within", "permanova", "permanova_conditional",
    "mean", "sd", "median", "r2_rmse", "cramers_v",
]


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def sd(xs, ddof=1):
    xs = list(xs)
    if len(xs) - ddof <= 0:
        return float("nan")
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - ddof))


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return float("nan")
    if n % 2:
        return xs[n // 2]
    return 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def closure(vec, total=1.0):

    s = sum(vec)
    if s <= 0:
        raise ValueError("closure: non-positive total")
    return [v * total / s for v in vec]


def multiplicative_replacement(vec, deltas):

    if len(vec) != len(deltas):
        raise ValueError("multiplicative_replacement: length mismatch")
    s = sum(vec)
    if s <= 0:
        raise ValueError("multiplicative_replacement: all components are zero; CLR is undefined")
    zero_mass = sum(d for v, d in zip(vec, deltas) if v <= 0)
    if zero_mass >= s:
        raise ValueError("multiplicative_replacement: replacement mass exceeds the total")
    scale = (s - zero_mass) / s
    out = []
    for v, d in zip(vec, deltas):
        out.append(d if v <= 0 else v * scale)
    return out


def clr(vec):

    if any(v <= 0 for v in vec):
        raise ValueError("clr: input contains non-positive values; apply zero replacement first")
    logs = [math.log(v) for v in vec]
    g = sum(logs) / len(logs)
    return [x - g for x in logs]


def euclidean(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def bray_curtis(a, b):

    num = sum(abs(x - y) for x, y in zip(a, b))
    den = sum(a) + sum(b)
    return num / den if den > 0 else 0.0


def jaccard(sa, sb):

    u = len(sa | sb)
    return 1.0 - (len(sa & sb) / u) if u else 0.0


def dist_matrix(items, fn):
    n = len(items)
    d = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            v = fn(items[i], items[j])
            d[i][j] = d[j][i] = v
    return d


def sq_dist_matrix(items, fn):
    d = dist_matrix(items, fn)
    return [[v * v for v in row] for row in d]


def ss_total(d2):
    n = len(d2)
    s = 0.0
    for i in range(n):
        row = d2[i]
        for j in range(i + 1, n):
            s += row[j]
    return s / n


def ss_within(d2, groups):

    s = 0.0
    for idx in groups.values():
        m = len(idx)
        if m < 2:
            continue
        acc = 0.0
        for a in range(m):
            row = d2[idx[a]]
            for b in range(a + 1, m):
                acc += row[idx[b]]
        s += acc / m
    return s


def _group_index(labels):
    g = defaultdict(list)
    for i, lab in enumerate(labels):
        g[lab].append(i)
    return g


def permanova(d2, labels, n_perm=9999, seed=20260912):

    n = len(labels)
    sizes = Counter(labels)
    a = len(sizes)
    res = {
        "n": n, "n_levels": a,
        "level_sizes": dict(sizes),
        "n_levels_ge2": sum(1 for v in sizes.values() if v >= 2),
        "n_perm": n_perm, "testable": False, "reason": "",
        "ss_total": None, "ss_among": None, "ss_within": None,
        "r2": None, "r2_adj": None, "r2_null_mean": None, "r2_excess": None,
        "F": None, "p": None, "df_among": a - 1, "df_within": n - a,
    }
    if a < 2:
        res["reason"] = "only 1 level; not testable"
        return res
    if n - a <= 0:
        res["reason"] = "every level holds exactly 1 sample (df_within = 0); not testable"
        return res
    sst = ss_total(d2)
    groups = _group_index(labels)
    ssw = ss_within(d2, groups)
    ssa = sst - ssw
    res.update(ss_total=sst, ss_among=ssa, ss_within=ssw,
               r2=(ssa / sst if sst > 0 else float("nan")))
    if sst <= 0 or ssw <= 0:
        res["reason"] = "within-group sum of squares is 0; F is undefined (samples coincide exactly)"
        res["F"] = float("inf")
        return res
    f_obs = (ssa / (a - 1)) / (ssw / (n - a))
    rng = random.Random(seed)
    perm = list(labels)
    ge = 0
    r2_null = 0.0
    for _ in range(n_perm):
        rng.shuffle(perm)
        w = ss_within(d2, _group_index(perm))
        f = ((sst - w) / (a - 1)) / (w / (n - a)) if w > 0 else float("inf")
        r2_null += (sst - w) / sst
        if f >= f_obs:
            ge += 1

    ms_t = sst / (n - 1)
    ms_w = ssw / (n - a)
    r2_null /= n_perm
    res.update(testable=True, F=f_obs, p=(ge + 1) / (n_perm + 1),
               r2_adj=(1 - ms_w / ms_t) if ms_t > 0 else float("nan"),
               r2_null_mean=r2_null, r2_excess=res["r2"] - r2_null)
    return res


def permanova_conditional(d2, strata, labels, n_perm=9999, seed=20260912):

    n = len(labels)
    strata = list(strata)
    labels = list(labels)
    by_stratum = defaultdict(list)
    for i, s in enumerate(strata):
        by_stratum[s].append(i)

    df_fac = 0
    informative = 0
    for s, idx in by_stratum.items():
        lv = Counter(labels[i] for i in idx)
        if len(lv) >= 2:
            df_fac += len(lv) - 1
            informative += len(idx)
    res = {
        "n": n, "df_factor": df_fac, "n_in_varying_strata": informative,
        "n_perm": n_perm, "testable": False, "reason": "",
        "ss_total": None, "ss_study": None, "ss_factor_given_study": None,
        "ss_residual": None, "r2_partial": None,
        "r2_of_within_study": None, "F": None, "p": None,
    }
    if df_fac == 0:
        res["reason"] = ("the factor does not vary within any single study -- "
                         "completely confounded with study; no conditional "
                         "test is possible")
        return res
    sst = ss_total(d2)
    ssw_study = ss_within(d2, _group_index(strata))
    cells = _group_index(list(zip(strata, labels)))
    ssw_cell = ss_within(d2, cells)
    ss_fac = ssw_study - ssw_cell
    df_res = n - len(cells)
    res.update(ss_total=sst, ss_study=sst - ssw_study,
               ss_factor_given_study=ss_fac, ss_residual=ssw_cell,
               r2_partial=(ss_fac / sst if sst > 0 else float("nan")),
               r2_of_within_study=(ss_fac / ssw_study if ssw_study > 0 else float("nan")),
               df_residual=df_res)
    if df_res <= 0 or ssw_cell <= 0:
        res["reason"] = ("residual df = 0 (each study x level cell holds "
                         "exactly 1 sample); not testable")
        return res
    f_obs = (ss_fac / df_fac) / (ssw_cell / df_res)
    rng = random.Random(seed)
    cur = list(labels)
    ge = 0
    for _ in range(n_perm):
        for s, idx in by_stratum.items():
            vals = [cur[i] for i in idx]
            rng.shuffle(vals)
            for i, v in zip(idx, vals):
                cur[i] = v
        w = ss_within(d2, _group_index(list(zip(strata, cur))))
        fac = ssw_study - w
        f = (fac / df_fac) / (w / df_res) if w > 0 else float("inf")
        if f >= f_obs:
            ge += 1
    res.update(testable=True, F=f_obs, p=(ge + 1) / (n_perm + 1))
    return res


def r2_rmse(obs, pred, ref):

    ss_p = sum((o - p) ** 2 for o, p in zip(obs, pred))
    ss_r = sum((o - r) ** 2 for o, r in zip(obs, ref))
    n = len(obs)
    return {
        "n": n,
        "rmse": math.sqrt(ss_p / n) if n else float("nan"),
        "rmse_ref": math.sqrt(ss_r / n) if n else float("nan"),
        "r2": (1 - ss_p / ss_r) if ss_r > 0 else float("nan"),
    }


def cramers_v(xs, ys):

    n = len(xs)
    if n == 0:
        return float("nan")
    tab = Counter(zip(xs, ys))
    rx = Counter(xs)
    ry = Counter(ys)
    chi2 = 0.0
    for i in rx:
        for j in ry:
            e = rx[i] * ry[j] / n
            o = tab.get((i, j), 0)
            chi2 += (o - e) ** 2 / e
    k = min(len(rx), len(ry))
    if k < 2:
        return float("nan")
    return math.sqrt((chi2 / n) / (k - 1))
