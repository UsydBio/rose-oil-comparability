import collections
import itertools
import json
import math
import os
import re
import zlib

from .. import config
from .. import contrast_analysis as CA

SEED = 20260913
MIN_SHARED = 3

MEASUREMENTS = config.data_path("measurements_long.tsv")
SMILES_CACHE = config.data_path("smiles_cache.json")


CONTROL_RE = re.compile(
    r"(?<![A-Za-z])(control|untreated|un-?treated|ctrl|unstressed|"
    r"check\s+variety|normal\s+conditions)(?![A-Za-z])", re.I)


UNORIENTABLE_AXES = {
    "origin": "level names are each paper's own villages, plots or compass "
              "aspects; no two studies share one",
    "colour_morph": "one study only, and the two morphs have no cross-study "
                    "canonical labels",
    "sampling_mode": "analytical, not biological (HS vs SPME), and one study only",
    "mixed_stage_diel": "superseded by three derived contrasts",
}


def load_blocks(table_path=MEASUREMENTS, drop_aggregate_rows=True):

    from .. import contrasts as K

    samples, kept, dropped = K.load_samples(table_path)
    cells, names, _ = CA.load_value_cells(table_path)
    specs, genuine = CA.build_specs(samples)

    superseded = set()
    for s in specs:
        if s.parent:
            superseded.add(s.parent)

    blocks, info = {}, {}
    for spec in specs:
        if spec.cid in superseded:
            info[spec.cid] = {"status": "superseded_parent"}
            continue
        b, meta = _block(spec, samples, cells, names, drop_aggregate_rows)
        info[spec.cid] = meta
        if b is not None:
            blocks[spec.cid] = b
    info["_corpus"] = {
        "n_sample_columns": len(samples), "rows_in_scope": kept,
        "rows_out_of_scope": dropped, "n_value_cells": len(cells),
        "n_specs": len(specs), "n_genuine": len(genuine),
        "n_superseded_parents": len(superseded),
    }
    return blocks, names, info


def _block(spec, samples, cells, names, drop_aggregate_rows):
    all_mats = CA.materials_of(spec, samples)
    mats, dropped_ct, _ = CA.drop_cross_table_pseudoreplicates(all_mats, spec.level_of)
    ids = [m[0] if isinstance(m[0], str) else "|".join(str(x) for x in m[0])
           for m in mats]
    cols = [m[1] for m in mats]
    labels = [spec.level_of[c.key] for c in cols]

    basis = sorted(CA.shared_basis(cols))
    tabs = set(c.table_id for c in cols)
    meta = {
        "axis": spec.axis,
        "study_id": spec.study_id,
        "measure_kind": spec.measure_kind,
        "ordered": bool(spec.ordered),
        "n_basis_before_aggregate_drop": len(basis),
        "n_aggregate_rows": sum(1 for c in basis if CA.looks_aggregate(names.get(c, ""))),
        "basis_preselected": any(
            any(t in tb.lower() for t in CA.PRESELECTED_TABLE_TOKENS) for tb in tabs),
        "caveat": CA.CAVEATS.get(spec.cid, ""),
        "confounded_with": spec.enum.get("confounded_with", ""),
        "n_cross_table_pseudoreps_dropped": len(dropped_ct),
    }
    if drop_aggregate_rows:
        basis = [c for c in basis if not CA.looks_aggregate(names.get(c, ""))]
    meta["n_basis"] = len(basis)
    if len(basis) < MIN_SHARED:
        meta["status"] = "basis_below_%d" % MIN_SHARED
        return None, meta

    raw = CA.value_matrix(cols, basis, cells)
    keep, dropped_dup = CA.drop_identical(raw, labels, ids)
    raw = [raw[i] for i in keep]
    labels = [labels[i] for i in keep]
    ids = [ids[i] for i in keep]
    kz, closed, clrs, n_capped = CA.clr_block(raw)
    labels = [labels[i] for i in kz]
    ids = [ids[i] for i in kz]

    meta.update({"status": "ok", "n_duplicate_materials_dropped": len(dropped_dup),
                 "n_zero_mass_capped": n_capped, "n_analysed": len(clrs),
                 "n_levels_analysed": len(set(labels)),
                 "n_basis_inchikey": sum(1 for c in basis if c.startswith("sk:"))})
    return {"spec": spec, "cid": spec.cid, "study": spec.study_id, "axis": spec.axis,
            "clr": clrs, "closed": closed, "labels": labels, "basis": basis,
            "ids": ids, "meta": meta}, meta


def oriented_pairs(block):

    spec = block["spec"]
    axis = block["axis"]
    levels = sorted(set(block["labels"]))
    out = []

    if spec.ordered:
        ls = sorted(set(block["labels"]), key=lambda l: spec.order_of[l])
        if len(ls) >= 2:
            out.append((ls[0], ls[-1], "earliest->latest", "strict"))
        return out

    if axis == "treatment":
        ctl = tuple(l for l in levels if CONTROL_RE.search(l))
        per = tuple(l for l in levels if l not in ctl)
        if ctl and per:
            out.append((ctl, per, "control->perturbed", "oriented"))
        return out

    canon = {}
    for l in levels:
        c = CA.canon_level(axis, l)
        if c:
            canon.setdefault(c, l)
    for u, v in itertools.combinations(sorted(canon), 2):
        out.append((canon[u], canon[v], "%s->%s" % (u, v), "strict"))
    return out


def subcompositional_clr(block, keys):

    idx = [j for j, k in enumerate(block["basis"]) if k in keys]
    out = []
    for row in block["clr"]:
        vals = [row[j] for j in idx]
        m = sum(vals) / len(vals)
        out.append([v - m for v in vals])
    return [block["basis"][j] for j in idx], out


def level_deltas(block, level_a, level_b, restrict_to=None):

    a = level_a if isinstance(level_a, tuple) else (level_a,)
    b = level_b if isinstance(level_b, tuple) else (level_b,)
    if restrict_to is not None:
        keys, clr = subcompositional_clr(block, restrict_to)
    else:
        keys, clr = block["basis"], block["clr"]
    out = {}
    for j, key in enumerate(keys):
        ca = [r[j] for r, l in zip(clr, block["labels"]) if l in a]
        cb = [r[j] for r, l in zip(clr, block["labels"]) if l in b]
        if not ca or not cb:
            continue
        jj = block["basis"].index(key)
        ra = [r[jj] for r, l in zip(block["closed"], block["labels"]) if l in a]
        rb = [r[jj] for r, l in zip(block["closed"], block["labels"]) if l in b]
        out[key] = (CA.S.mean(cb) - CA.S.mean(ca), CA.S.mean(rb) - CA.S.mean(ra))
    return out


def effect_rows(blocks, names, inchikey_only=True, dedupe_measure_kind=True):

    rows = []
    for cid in sorted(blocks):
        b = blocks[cid]
        for la, lb, key, tier in oriented_pairs(b):
            d = level_deltas(b, la, lb)
            for ck, (clr, raw) in sorted(d.items()):
                if inchikey_only and not ck.startswith("sk:"):
                    continue
                if clr == 0.0:
                    continue
                rows.append({
                    "contrast_id": cid, "study_id": b["study"], "axis": b["axis"],
                    "measure_kind": b["spec"].measure_kind,
                    "orient_key": key, "tier": tier,
                    "level_a": la if isinstance(la, str) else " + ".join(la),
                    "level_b": lb if isinstance(lb, str) else " + ".join(lb),
                    "compound_key": ck, "compound_name": names.get(ck, ck),
                    "effect_clr": clr, "effect_raw_pct": raw,
                    "sign": 1 if clr > 0 else -1,
                    "sign_raw_pct": 0 if raw == 0 else (1 if raw > 0 else -1),
                    "abs_effect_clr": abs(clr),
                    "n_basis": len(b["basis"]),
                    "basis_preselected": b["meta"].get("basis_preselected", False),
                    "caveat": b["meta"].get("caveat", ""),
                })
    if dedupe_measure_kind:
        levels_of = {}
        for r in rows:
            levels_of.setdefault(r["contrast_id"], set()).update(
                (r["level_a"], r["level_b"]))
        best = {}
        for r in rows:
            k = (r["study_id"], r["axis"], r["orient_key"],
                 frozenset(levels_of[r["contrast_id"]]))
            if k not in best or r["n_basis"] > best[k][1]:
                best[k] = (r["contrast_id"], r["n_basis"])
        keep = {v[0] for v in best.values()}
        rows = [r for r in rows if r["contrast_id"] in keep]
    return rows


def binomial_two_sided(agree, total):

    if total <= 0:
        return None
    k = max(agree, total - agree)

    log_half_n = -total * math.log(2.0)
    terms = []
    for i in range(k, total + 1):
        terms.append(math.lgamma(total + 1) - math.lgamma(i + 1)
                     - math.lgamma(total - i + 1) + log_half_n)
    top = max(terms)
    tail = top + math.log(sum(math.exp(t - top) for t in terms))
    return min(1.0, 2.0 * math.exp(tail))


def cross_study_pairs(blocks, min_shared=MIN_SHARED, same_study=False,
                      subcompositional=False):

    items = []
    for cid in sorted(blocks):
        b = blocks[cid]
        for la, lb, key, tier in oriented_pairs(b):
            items.append((cid, b, la, lb, key, tier))

    out = []
    for (ca, A, laa, lab, ka, ta), (cb, B, lba, lbb, kb, tb) in \
            itertools.combinations(items, 2):
        if ka != kb or A["axis"] != B["axis"]:

            continue
        same = A["study"] == B["study"]
        if same != same_study:
            continue
        if same and ca == cb:
            continue
        common = sorted(set(k for k in A["basis"] if k.startswith("sk:")) &
                        set(k for k in B["basis"] if k.startswith("sk:")))
        if subcompositional:
            if len(common) < 2:
                continue
            da = level_deltas(A, laa, lab, restrict_to=set(common))
            db = level_deltas(B, lba, lbb, restrict_to=set(common))
        else:
            da = level_deltas(A, laa, lab)
            db = level_deltas(B, lba, lbb)
        common = [k for k in common if k in da and k in db]
        zb = (_zero_backed(A, laa, lab, common) |
              _zero_backed(B, lba, lbb, common))
        agree = dis = tie = 0
        agree_raw = dis_raw = 0
        agree_q = dis_q = 0
        detail = []
        for c in common:
            x, xr = da[c]
            y, yr = db[c]
            if x == 0 or y == 0:
                tie += 1
            elif (x > 0) == (y > 0):
                agree += 1
                if c not in zb:
                    agree_q += 1
            else:
                dis += 1
                if c not in zb:
                    dis_q += 1
            if xr != 0 and yr != 0:
                if (xr > 0) == (yr > 0):
                    agree_raw += 1
                else:
                    dis_raw += 1
            detail.append((c, x, y))
        n = agree + dis
        out.append({
            "axis": A["axis"], "orient_key": ka, "tier": ta,
            "basis": "shared_subcomposition" if subcompositional else "own_full_basis",
            "status": "ok" if len(common) >= min_shared else "too_few_shared",
            "study_a": A["study"], "contrast_a": ca,
            "level_pair_a": _pair_text(laa, lab),
            "study_b": B["study"], "contrast_b": cb,
            "level_pair_b": _pair_text(lba, lbb),
            "n_shared_inchikey": len(common),
            "n_agree": agree, "n_disagree": dis, "n_tie": tie,
            "frac_agree": (agree / float(n)) if n else None,
            "n_agree_raw_pct": agree_raw, "n_disagree_raw_pct": dis_raw,
            "n_zero_backed": len(zb),
            "n_agree_quantified": agree_q, "n_disagree_quantified": dis_q,
            "frac_agree_quantified": (agree_q / float(agree_q + dis_q))
                                     if (agree_q + dis_q) else None,
            "binomial_p_quantified": binomial_two_sided(agree_q, agree_q + dis_q)
                                     if (agree_q + dis_q) >= min_shared else None,
            "binomial_p": binomial_two_sided(agree, n) if len(common) >= min_shared else None,
            "detail": detail,
        })
    return out


def _zero_backed(block, level_a, level_b, keys):

    a = level_a if isinstance(level_a, tuple) else (level_a,)
    b = level_b if isinstance(level_b, tuple) else (level_b,)
    rows = [i for i, l in enumerate(block["labels"]) if l in a or l in b]
    idx = {k: j for j, k in enumerate(block["basis"])}
    out = set()
    for k in keys:
        j = idx.get(k)
        if j is None:
            continue
        if any(block["closed"][i][j] <= 0 for i in rows):
            out.add(k)
    return out


def _pair_text(a, b):
    a = a if isinstance(a, str) else " + ".join(a)
    b = b if isinstance(b, str) else " + ".join(b)
    return "%s -> %s" % (a, b)


def agreement_by_effect_threshold(pairs, thresholds=(0.0, 0.25, 0.5, 1.0)):

    out = []
    for p in pairs:
        if p["status"] != "ok":
            continue
        for t in thresholds:
            ag = dis = 0
            for c, x, y in p["detail"]:
                if x == 0 or y == 0 or abs(x) < t or abs(y) < t:
                    continue
                if (x > 0) == (y > 0):
                    ag += 1
                else:
                    dis += 1
            out.append({
                "relation": p.get("relation", ""), "tier": p["tier"],
                "axis": p["axis"], "orient_key": p["orient_key"],
                "study_a": p["study_a"], "study_b": p["study_b"],
                "contrast_a": p["contrast_a"], "contrast_b": p["contrast_b"],
                "min_abs_effect_clr": t, "n_agree": ag, "n_disagree": dis,
                "n_compared": ag + dis,
                "frac_agree": (ag / float(ag + dis)) if (ag + dis) else None,
                "binomial_p": binomial_two_sided(ag, ag + dis) if (ag + dis) else None,
            })
    return out



_ELEMENT = re.compile(r"Cl|Br|[BCNOPSFIbcnops]")
_ATOMIC_WEIGHT = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007, "S": 32.06,
                  "F": 18.998, "Cl": 35.45, "Br": 79.904, "I": 126.904,
                  "P": 30.974, "Si": 28.085, "B": 10.81}


def smiles_descriptors(smiles, formula=""):

    s = smiles or ""
    el = collections.Counter(_ELEMENT.findall(s))
    n_c = el["C"] + el["c"]
    n_o = el["O"] + el["o"]
    n_n = el["N"] + el["n"]
    n_s = el["S"] + el["s"]
    n_hal = el["F"] + el["Cl"] + el["Br"] + el["I"]
    aromatic = el["c"] + el["n"] + el["o"] + el["s"]
    heavy = sum(el.values())

    rings = len(re.findall(r"(?<!%)\d", re.sub(r"\[[^\]]*\]", "", s))) // 2
    n_double, n_triple, n_branch = s.count("="), s.count("#"), s.count("(")
    mw = 0.0
    for sym, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        if sym in _ATOMIC_WEIGHT:
            mw += _ATOMIC_WEIGHT[sym] * (int(cnt) if cnt else 1)
    hetero = n_o + n_n + n_s + n_hal
    return {
        "n_heavy": float(heavy), "n_C": float(n_c), "n_O": float(n_o),
        "n_N": float(n_n), "n_S": float(n_s), "n_halogen": float(n_hal),
        "n_aromatic_atoms": float(aromatic),
        "frac_aromatic": aromatic / heavy if heavy else 0.0,
        "n_rings": float(rings), "n_double": float(n_double),
        "n_triple": float(n_triple), "n_branch": float(n_branch),
        "n_stereocentres": float(s.count("@")), "mw": mw,
        "dbe": float(n_double + n_triple + rings),
        "o_over_c": n_o / n_c if n_c else 0.0,
        "heteroatom_frac": hetero / heavy if heavy else 0.0,
        "is_hydrocarbon": 1.0 if hetero == 0 else 0.0,
        "fg_acid": 1.0 if re.search(r"C\(=O\)O(?![A-Za-z0-9@\[])", s) else 0.0,
        "fg_ester": 1.0 if re.search(r"C\(=O\)O[C@\[c]", s) else 0.0,
        "fg_carbonyl": 1.0 if "=O" in s else 0.0,
        "fg_hydroxyl": 1.0 if re.search(r"(?<![=\)])O(?![A-Za-z0-9@\[=\)])", s) else 0.0,
        "fg_ether": 1.0 if re.search(r"[Cc]O[Cc]", s) else 0.0,
        "terpene_c10": 1.0 if n_c == 10 else 0.0,
        "terpene_c15": 1.0 if n_c == 15 else 0.0,
        "long_chain": 1.0 if n_c >= 17 else 0.0,
    }


DESCRIPTOR_FIELDS = sorted(smiles_descriptors("CCO", "C2H6O").keys())


def compound_structures(compound_keys, table_path=MEASUREMENTS, cache_path=None):

    import csv
    if cache_path is None:
        cache_path = os.path.abspath(SMILES_CACHE)
    want = {k[3:] for k in compound_keys if k.startswith("sk:")}
    cid_votes = collections.defaultdict(collections.Counter)
    name_votes = collections.defaultdict(collections.Counter)
    formula_votes = collections.defaultdict(collections.Counter)
    with open(table_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sk = r.get("inchikey_skeleton") or ""
            if sk not in want:
                continue
            if r.get("pubchem_cid"):
                cid_votes[sk][r["pubchem_cid"]] += 1
            if r.get("compound_clean_name"):
                name_votes[sk][r["compound_clean_name"]] += 1
            if r.get("compound_formula"):
                formula_votes[sk][r["compound_formula"]] += 1
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8") as fh:
            cache = json.load(fh)
    out = {}
    for sk in sorted(want):
        cid = cid_votes[sk].most_common(1)[0][0] if cid_votes.get(sk) else ""
        hit = cache.get(cid, {}) if cid else {}
        formula = hit.get("formula") or (
            formula_votes[sk].most_common(1)[0][0] if formula_votes.get(sk) else "")
        smiles = hit.get("smiles", "")
        out["sk:" + sk] = {
            "pubchem_cid": cid, "smiles": smiles, "formula": formula,
            "name": name_votes[sk].most_common(1)[0][0] if name_votes.get(sk) else sk,
            "descriptors": smiles_descriptors(smiles, formula),
            "has_smiles": bool(smiles),
        }
    return out


def _hashed_substrings(smiles, n_bits=64, k=4):

    v = [0.0] * n_bits
    s = smiles or ""
    for i in range(max(0, len(s) - k + 1)):

        v[zlib.crc32(s[i:i + k].encode("utf-8")) % n_bits] = 1.0
    return v


def feature_matrix(rows, structures, mode="struct", axes=None, compounds=None,
                   n_bits=64):

    import numpy as np

    if axes is None:
        axes = sorted({r["axis"] for r in rows})
    if compounds is None:
        compounds = sorted({r["compound_key"] for r in rows})
    ax_ix = {a: i for i, a in enumerate(axes)}
    c_ix = {c: i for i, c in enumerate(compounds)}

    cols, feats = [], []
    if mode not in ("compound_onehot", "axis_onehot"):
        feats += ["desc:" + f for f in DESCRIPTOR_FIELDS]
    if mode in ("struct+axis", "axis_onehot"):
        feats += ["axis:" + a for a in axes]
    if mode == "struct+ngram":
        feats += ["ngram:%d" % i for i in range(n_bits)]
    if mode == "compound_onehot":
        feats += ["cmp:" + c for c in compounds]

    for r in rows:
        st = structures.get(r["compound_key"], {})
        d = st.get("descriptors") or smiles_descriptors("", "")
        v = []
        if mode not in ("compound_onehot", "axis_onehot"):
            v += [d[f] for f in DESCRIPTOR_FIELDS]
        if mode in ("struct+axis", "axis_onehot"):
            a = [0.0] * len(axes)
            if r["axis"] in ax_ix:
                a[ax_ix[r["axis"]]] = 1.0
            v += a
        if mode == "struct+ngram":
            v += _hashed_substrings(st.get("smiles", ""), n_bits=n_bits)
        if mode == "compound_onehot":
            c = [0.0] * len(compounds)
            if r["compound_key"] in c_ix:
                c[c_ix[r["compound_key"]]] = 1.0
            v += c
        cols.append(v)
    return np.asarray(cols, dtype=float), feats


def leave_study_out(rows, structures, seed=SEED, min_test=10, min_train_studies=3):

    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    studies = sorted({r["study_id"] for r in rows})
    fold_rows, skipped = [], []

    for held in studies:
        tr = [r for r in rows if r["study_id"] != held]
        te = [r for r in rows if r["study_id"] == held]
        if len(te) < min_test:
            skipped.append((held, len(te), "held-out study has < %d instances" % min_test))
            continue
        if len({r["study_id"] for r in tr}) < min_train_studies:
            skipped.append((held, len(te), "fewer than %d training studies" % min_train_studies))
            continue

        y_te = np.array([r["sign"] for r in te])
        y_tr = np.array([r["sign"] for r in tr])

        maj_sign = 1 if (y_tr > 0).sum() >= (y_tr <= 0).sum() else -1

        cmp_key = collections.defaultdict(lambda: [0, 0])
        pair_key = collections.defaultdict(lambda: [0, 0])
        for r in tr:
            cmp_key[(r["axis"], r["compound_key"])][0 if r["sign"] > 0 else 1] += 1
            pair_key[(r["axis"], r["orient_key"], r["compound_key"])][
                0 if r["sign"] > 0 else 1] += 1

        def lookup(table, key):
            if key not in table:
                return None
            up, down = table[key]
            if up == down:
                return None
            return 1 if up > down else -1

        preds = {"coin_flip": [None] * len(te),
                 "train_majority": [maj_sign] * len(te),
                 "compound_majority": [lookup(cmp_key, (r["axis"], r["compound_key"]))
                                       for r in te],
                 "pair_majority": [lookup(pair_key, (r["axis"], r["orient_key"],
                                                     r["compound_key"])) for r in te]}

        axes = sorted({r["axis"] for r in tr})
        Xtr, _ = feature_matrix(tr, structures, "struct+axis", axes=axes)
        Xte, _ = feature_matrix(te, structures, "struct+axis", axes=axes)
        sc = StandardScaler().fit(Xtr)
        clf = LogisticRegression(max_iter=4000, random_state=seed)
        try:
            clf.fit(sc.transform(Xtr), y_tr)
            preds["structure_logreg"] = list(clf.predict(sc.transform(Xte)))
        except Exception as exc:
            preds["structure_logreg"] = [None] * len(te)
            skipped.append((held, len(te), "logreg failed: %s" % exc))

        answered = np.array([p is not None for p in preds["compound_majority"]])
        for name, p in preds.items():
            ok = np.array([q is not None for q in p])
            hit = np.array([1 if (q is not None and q == t) else 0
                            for q, t in zip(p, y_te)])
            both = ok & answered
            fold_rows.append({
                "held_out_study": held, "model": name,
                "n_test": len(te), "n_answered": int(ok.sum()),
                "coverage": round(float(ok.mean()), 4),
                "accuracy": round(float(hit[ok].mean()), 4) if ok.any() else None,
                "n_correct": int(hit[ok].sum()),
                "accuracy_on_shared_subset":
                    round(float(hit[both].mean()), 4) if both.any() else None,
                "n_shared_subset": int(both.sum()),
                "binomial_p": binomial_two_sided(int(hit[ok].sum()), int(ok.sum()))
                              if ok.any() else None,
                "test_axes": ";".join(sorted({r["axis"] for r in te})),
                "axes_seen_in_training": ";".join(
                    sorted({r["axis"] for r in te} & set(axes))),
            })

        for fr in fold_rows[-len(preds):]:
            if fr["model"] == "coin_flip":
                fr.update({"n_answered": len(te), "coverage": 1.0, "accuracy": 0.5,
                           "n_correct": None, "binomial_p": None,
                           "accuracy_on_shared_subset": 0.5,
                           "n_shared_subset": int(answered.sum())})
    return fold_rows, skipped


def summarise_folds(fold_rows):

    import numpy as np
    by = collections.defaultdict(list)
    for r in fold_rows:
        by[r["model"]].append(r)
    out = []
    for model in sorted(by):
        rs = [r for r in by[model] if r["accuracy"] is not None]
        if not rs:
            continue
        n_ans = sum(r["n_answered"] for r in rs)
        n_hit = sum(r["n_correct"] or 0 for r in rs)
        sub = [r for r in by[model] if r["accuracy_on_shared_subset"] is not None]
        out.append({
            "model": model, "n_folds": len(rs),
            "macro_accuracy": round(float(np.mean([r["accuracy"] for r in rs])), 4),
            "macro_accuracy_sd": round(float(np.std([r["accuracy"] for r in rs])), 4),
            "pooled_accuracy": round(n_hit / n_ans, 4) if n_ans and model != "coin_flip" else (0.5 if model == "coin_flip" else None),
            "n_answered": n_ans,
            "mean_coverage": round(float(np.mean([r["coverage"] for r in by[model]])), 4),
            "macro_accuracy_shared_subset": round(
                float(np.mean([r["accuracy_on_shared_subset"] for r in sub])), 4) if sub else None,
            "pooled_binomial_p": binomial_two_sided(n_hit, n_ans)
                                 if model != "coin_flip" and n_ans else None,
        })
    return out


def accuracy_by_effect_size(rows, n_bins=4):

    import numpy as np
    vals = np.array([r["abs_effect_clr"] for r in rows])
    if len(vals) == 0:
        return []
    edges = np.quantile(vals, np.linspace(0, 1, n_bins + 1))
    studies = sorted({r["study_id"] for r in rows})
    out = []
    for b in range(n_bins):
        lo, hi = edges[b], edges[b + 1]
        hit = tot = 0
        for held in studies:
            tr = [r for r in rows if r["study_id"] != held]
            table = collections.defaultdict(lambda: [0, 0])
            for r in tr:
                table[(r["axis"], r["compound_key"])][0 if r["sign"] > 0 else 1] += 1
            for r in rows:
                if r["study_id"] != held:
                    continue
                if not (lo <= r["abs_effect_clr"] <= hi if b == n_bins - 1
                        else lo <= r["abs_effect_clr"] < hi):
                    continue
                k = (r["axis"], r["compound_key"])
                if k not in table:
                    continue
                up, down = table[k]
                if up == down:
                    continue
                pred = 1 if up > down else -1
                tot += 1
                hit += int(pred == r["sign"])
        out.append({"bin": b + 1, "abs_clr_low": round(float(lo), 4),
                    "abs_clr_high": round(float(hi), 4), "n_answered": tot,
                    "accuracy": round(hit / tot, 4) if tot else None,
                    "binomial_p": binomial_two_sided(hit, tot) if tot else None})
    return out


def probe_representations(rows, structures, seed=SEED):

    from . import evaluate

    seen, uniq = set(), []
    for r in rows:
        k = (r["study_id"], r["compound_key"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(r)

    codes = {s_: i for i, s_ in enumerate(sorted({r["study_id"] for r in uniq}))}
    labels = [codes[r["study_id"]] for r in uniq]

    out = {"n_probe_rows": len(uniq), "n_studies": len(codes)}
    for mode in ("struct", "axis_onehot", "struct+axis", "struct+ngram",
                 "compound_onehot"):
        X, feats = feature_matrix(uniq, structures, mode)
        res = evaluate.study_identity_probe(X, labels, seed=seed)
        res["n_features"] = len(feats)
        out[mode] = res
    return out
