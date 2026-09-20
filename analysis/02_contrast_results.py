import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from comparability import config, contrasts as K
from comparability import contrast_analysis as CA

MEASUREMENTS = config.data_path("measurements_long.tsv")
OUT_TSV = config.results_path("contrast_results.tsv")
OUT_COMP = config.results_path("contrast_compound_effects.tsv")
OUT_DIR = config.results_path("contrast_direction_agreement.tsv")

ALPHA = 0.05
TOP_DRIVERS_PER_CONTRAST = 12

FIELDS = [
    "contrast_id", "kind", "parent", "superseded_by", "study_id", "doi", "pmcid",
    "tier", "axis", "axis_biological", "axis_note", "factor_field", "level_source",
    "measure_kind", "values_comparable", "identity_status", "platform",
    "column_polarity", "n_columns", "n_materials_enumerated", "n_materials",
    "n_cross_table_pseudoreps_dropped", "n_levels_with_cross_table_pseudoreps",
    "cross_table_pseudoreps", "pseudorep_pair_aitchison",
    "between_level_aitchison_predrop", "n_duplicate_materials_dropped",
    "duplicate_pairs", "n_allzero_materials_dropped", "n_analysed",
    "n_levels", "n_levels_analysed", "level_sizes", "min_n_per_level",
    "max_n_per_level", "n_levels_replicated", "n_shared_enumerated", "n_basis",
    "n_basis_inchikey", "n_basis_aggregate_rows", "n_tables", "basis_preselected", "n_zero_mass_capped",
    "table_confounded", "confounded_with", "status", "inference_class",
    "df_among", "df_within", "n_perm", "p_min_achievable",
    "r2_aitchison", "r2_excess_aitchison", "r2_adj_aitchison", "F_aitchison",
    "p_aitchison", "q_aitchison",
    "r2_bray", "F_bray", "p_bray", "q_bray",
    "mean_between_level_aitchison", "mean_within_level_aitchison",
    "separation_ratio",
    "mantel_r", "mantel_p", "mantel_q", "mantel_n_orders", "mantel_exact",
    "n_trend_q05", "trend_q_floor", "n_drivers_q05", "top_drivers", "verdict", "verdict_reason",
    "caveat", "derivation", "study_title",
]

COMP_FIELDS = ["contrast_id", "axis", "study_id", "rank", "compound_key",
               "compound_name", "statistic", "effect_clr", "effect_raw_pct",
               "p", "q", "trend_rho", "trend_p", "trend_q"]

DIR_FIELDS = ["axis", "comparison", "status", "study_a", "contrast_a", "level_pair_a",
              "study_b", "contrast_b", "level_pair_b", "n_shared_inchikey",
              "n_agree", "n_disagree", "n_tie", "frac_agree",
              "n_agree_raw_pct", "n_disagree_raw_pct", "sign_test_p",
              "shared_compound_signs"]


def verdict_of(r):

    pmin = r.get("p_min_achievable")
    if r.get("status") == "tested":
        q = r.get("q_aitchison")
        hard = []
        if r.get("basis_preselected"):
            hard.append("basis pre-filtered for significance (DAMs sheet)")
        if r.get("confounded_with"):
            hard.append("factor fully covaries with " + r["confounded_with"])
        if r.get("caveat"):
            hard.append("see caveat")
        if q is not None and q <= ALPHA:
            if r.get("inference_class") == "inferential" and not hard \
                    and (r.get("n_basis") or 0) >= 10:
                return "publishable", "FDR q=%.4f on %d shared compounds, every level n>=2" % (
                    q, r.get("n_basis") or 0)
            return "suggestive", "FDR q=%.4f but: %s" % (q, "; ".join(hard) or
                                                         "only one level replicated")
        if pmin is not None and pmin >= ALPHA:
            return ("design_cannot_reject",
                    "permutation p floor is %.3f >= %.2f: this design cannot reach "
                    "significance whatever the data" % (pmin, ALPHA))
        return "underpowered_null", "tested, q=%s, no effect detected" % CA.fmt(q)
    if r.get("mantel_q") is not None and r["mantel_q"] <= ALPHA:
        return ("suggestive", "no within-level replication, but the ordered-level "
                "Mantel trend survives FDR (q=%.4f)" % r["mantel_q"])
    return ("description_only",
            "every level n=1 after removing cross-table re-reports: no group-variance "
            "test exists (df_within=0)")


def main():
    print("[1/5] Loading the measurement table ...")
    samples, kept, dropped = K.load_samples(MEASUREMENTS)
    cells, names, n_cell_rows = CA.load_value_cells(MEASUREMENTS)
    print("      %d measurement columns in scope; %d rows in scope,"
          " %d rows out of scope; %d value cells"
          % (len(samples), kept, dropped, len(cells)))

    print("[2/5] Rebuilding the contrast inventory, splitting out derived contrasts ...")
    specs, genuine = CA.build_specs(samples)
    derived_parents = set(s.parent for s in specs if s.kind == "derived" and s.parent)
    children = {}
    for s in specs:
        if s.parent:
            children.setdefault(s.parent, []).append(s.cid)
    print("      %d genuine contrasts; %d derived; %d analysed in total"
          % (len(genuine), len(specs) - len(genuine), len(specs)))

    print("[3/5] Per-contrast tests (PERMANOVA, %d permutations, seed %d) ..." % (CA.N_PERM, CA.SEED))
    results, comp_rows, blocks = [], [], {}
    for spec in specs:
        res, comp, trend, block = CA.analyse_contrast(
            spec, samples, cells, names, n_perm=CA.N_PERM, seed=CA.SEED)
        res["axis_biological"] = spec.axis not in CA.ANALYTICAL_AXES
        res["superseded_by"] = ";".join(sorted(children.get(spec.cid, [])))
        results.append(res)
        if block is not None:
            blocks[spec.cid] = block
        trend_by_key = {}
        if trend:
            trend_by_key = dict((d["compound_key"], d) for d in trend["compounds"])
        for i, d in enumerate(comp[:TOP_DRIVERS_PER_CONTRAST]):
            t = trend_by_key.get(d["compound_key"], {})
            comp_rows.append({
                "contrast_id": res["contrast_id"], "axis": res["axis"],
                "study_id": res["study_id"], "rank": i + 1,
                "compound_key": d["compound_key"], "compound_name": d["compound_name"],
                "statistic": d["statistic"], "effect_clr": d["effect_clr"],
                "effect_raw_pct": d["effect_raw"], "p": d["p"], "q": d["q"],
                "trend_rho": t.get("rho"), "trend_p": t.get("p"), "trend_q": t.get("q"),
            })
        print("      %-56s %-17s basis=%-4d %s"
              % (spec.cid[:56], res.get("inference_class", ""), res.get("n_basis", 0),
                 res.get("status", "")))

    print("[4/5] Multiple comparisons: BH correction (three independent test families) ...")

    for src, dst in (("p_aitchison", "q_aitchison"), ("p_bray", "q_bray"),
                     ("mantel_p", "mantel_q")):
        ps = [r.get(src) for r in results]
        for r, q in zip(results, CA.benjamini_hochberg(ps)):
            r[dst] = q
    for r in results:
        v, why = verdict_of(r)
        r["verdict"], r["verdict_reason"] = v, why

    print("[5/5] Cross-study direction agreement (counts, not a pooled p) ...")
    by_axis = {}
    for r in results:
        by_axis.setdefault(r["axis"], []).append(r)
    dir_rows = []
    for axis in sorted(by_axis):
        sub = dict((r["contrast_id"], blocks[r["contrast_id"]])
                   for r in by_axis[axis]
                   if r["contrast_id"] in blocks and not r.get("superseded_by"))
        for d in CA.direction_agreement(sub, axis):
            dir_rows.append({
                "axis": d["axis"], "comparison": d["tag"], "status": d["status"],
                "study_a": d["study_a"], "contrast_a": d["contrast_a"],
                "level_pair_a": d["level_pair_a"], "study_b": d["study_b"],
                "contrast_b": d["contrast_b"], "level_pair_b": d["level_pair_b"],
                "n_shared_inchikey": d["n_shared_inchikey"], "n_agree": d["n_agree"],
                "n_disagree": d["n_disagree"], "n_tie": d["n_tie"],
                "frac_agree": d["frac_agree"], "n_agree_raw_pct": d["n_agree_raw"],
                "n_disagree_raw_pct": d["n_disagree_raw"],
                "sign_test_p": d["sign_test_p"],
                "shared_compound_signs": "; ".join(
                    "%s(%s%s)" % (names.get(c, c), "+" if x > 0 else "-",
                                  "+" if y > 0 else "-")
                    for c, x, y in d["detail"][:24]),
                "detail": [(names.get(c, c), x, y) for c, x, y in d["detail"]],
            })

    CA.write_tsv(OUT_TSV, FIELDS, results)
    CA.write_tsv(OUT_COMP, COMP_FIELDS, comp_rows)
    CA.write_tsv(OUT_DIR, DIR_FIELDS, dir_rows)
    print()
    print("Written:")
    for p in (OUT_TSV, OUT_COMP, OUT_DIR):
        print("  ", os.path.basename(p))
    counts = {}
    for r in results:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    print("Verdict distribution:", "; ".join("%s=%d" % kv for kv in sorted(counts.items())))


AXIS_TITLES = {
    "cultivar": "cultivar",
    "species": "species",
    "colour_morph": "colour_morph: flower colour morph (recorded in species_level)",
    "stage_floral": "stage: floral opening stage",
    "stage_fruit": "stage: fruit ripening (days after anthesis)",
    "diel": "diel: time of day (recorded in developmental_stage)",
    "mixed_stage_diel": "(not analysable) stage and diel combined in one field",
    "organ": "organ",
    "origin": "origin: geographic origin",
    "processing": "processing: processing / extraction method",
    "treatment": "treatment",
    "sampling_mode": "sampling_mode: sampling method (analytical, not biological)",
}

VERDICT_ORDER = ["publishable", "suggestive", "design_cannot_reject",
                 "underpowered_null", "description_only"]


if __name__ == "__main__":
    main()
