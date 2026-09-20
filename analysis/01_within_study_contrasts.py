import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from comparability import config, contrasts as C

MEASUREMENTS = config.data_path("measurements_long.tsv")
OUT_CONTRASTS = config.results_path("within_study_contrasts.tsv")

MIN_SHARED = 3


def main():
    samples, n_kept, n_dropped = C.load_samples(MEASUREMENTS)
    print("%d rows in composition scope, %d rows excluded" % (n_kept, n_dropped))
    print("In-scope measurement columns: %d, studies: %d"
          % (len(samples), len(set(s.study_id for s in samples.values()))))

    all_c = C.build_contrasts(samples)
    for c in all_c:
        c["genuine"] = C.is_genuine(c, MIN_SHARED)
        c["score"] = C.score_contrast(c)
    all_c.sort(key=lambda c: (-c["score"], c["study_id"], c["factor"]))

    genuine = [c for c in all_c if c["genuine"]]

    in_any = set()
    in_confirmed = set()
    mats_any = set()
    for c in genuine:
        for k in c["_samples"]:
            in_any.add(k)
            if c["identity_status"] == "confirmed":
                in_confirmed.add(k)
    per_study_in = {}
    for k in in_any:
        per_study_in.setdefault(k[0], []).append(samples[k])
    for st, ss_ in per_study_in.items():
        for mk in C._dedupe_materials(ss_):
            mats_any.add((st, mk))

    cols = ["study_id", "doi", "pmcid", "tier", "factor", "factor_kind", "level_source",
            "measure_kind", "n_levels", "n_columns", "n_materials",
            "min_materials_per_level", "max_materials_per_level",
            "n_levels_with_replication", "has_replication", "n_shared_compounds",
            "n_shared_quantified", "n_shared_any_sample", "n_tables", "table_confounded",
            "platform", "platform_constant", "column_polarity", "polarity_constant",
            "column_class", "identity_status", "values_comparable", "confounded_with",
            "suspect_levels", "genuine", "score", "levels"]
    C.write_tsv(OUT_CONTRASTS, cols, all_c)
    print("Wrote %s (%d candidate contrasts, %d of them genuine)"
          % (os.path.basename(OUT_CONTRASTS), len(all_c), len(genuine)))


if __name__ == "__main__":
    main()
