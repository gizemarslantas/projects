#!/usr/bin/env python3
import argparse
import json
import os
import time

import joblib
import numpy as np
import pandas as pd
import pgenlib
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, roc_auc_score


def read_table(path):
    return pd.read_csv(path, sep=r"\s+", dtype=str)


def iid_column(df, label):
    for col in ("#IID", "IID"):
        if col in df.columns:
            return col
    raise ValueError(f"{label}: IID sütunu bulunamadı. Sütunlar: {list(df.columns)}")


def impute_missing_inplace(g):
    missing = g < 0
    if not missing.any():
        return 0
    valid = ~missing
    valid_n = valid.sum(axis=1)
    sums = np.where(valid, g, 0.0).sum(axis=1, dtype=np.float64)
    means = np.divide(sums, valid_n, out=np.zeros_like(sums), where=valid_n > 0)
    row_idx = np.nonzero(missing)[0]
    g[missing] = means[row_idx]
    return int(missing.sum())


def make_x(g_variant_major, pcs):
    variant_n, sample_n = g_variant_major.shape
    x = np.empty((sample_n, variant_n + pcs.shape[1]), dtype=np.float32)
    x[:, :variant_n] = g_variant_major.T
    x[:, variant_n:] = pcs
    return x


def make_rf(n_estimators, n_jobs, seed, min_leaf, oob=False):
    return RandomForestClassifier(
        n_estimators=n_estimators,
        criterion="gini",
        max_features="sqrt",
        max_depth=None,
        min_samples_leaf=min_leaf,
        class_weight="balanced_subsample",
        bootstrap=True,
        max_samples=0.80,
        oob_score=oob,
        n_jobs=n_jobs,
        random_state=seed,
    )


def load_inputs(base, chrom):
    prefix = os.path.join(base, "per_chr", f"chr{chrom}", f"SLD_EA_chr{chrom}_merged")
    pgen_path = prefix + ".pgen"
    pvar_path = prefix + ".pvar"
    psam_path = prefix + ".psam"
    pheno_path = os.path.join(base, "SLD_EA_c2.pheno.tsv")

    for path in (pgen_path, pvar_path, psam_path, pheno_path):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    psam = read_table(psam_path)
    pheno = read_table(pheno_path)
    psam_iid = iid_column(psam, "PSAM")
    pheno_iid = iid_column(pheno, "PHENO")

    iids = psam[psam_iid].astype(str)
    if iids.duplicated().any():
        raise ValueError("PSAM içinde tekrarlanan IID var.")
    if pheno[pheno_iid].duplicated().any():
        raise ValueError("Fenotip dosyasında tekrarlanan IID var.")

    pheno = pheno.set_index(pheno_iid)
    missing_pheno = iids[~iids.isin(pheno.index)]
    if len(missing_pheno):
        raise ValueError(f"Fenotipi olmayan {len(missing_pheno)} örnek var.")

    raw_y = pd.to_numeric(pheno.loc[iids, "PHENO"], errors="raise").to_numpy()
    observed = set(np.unique(raw_y).tolist())
    if observed != {1, 2}:
        raise ValueError(f"PHENO değerleri 1/2 olmalı; bulunan: {sorted(observed)}")
    y = (raw_y == 2).astype(np.int8)

    pc_cols = []
    pcs = np.empty((len(y), 0), dtype=np.float32)
    if not np.isfinite(pcs).all():
        raise ValueError("PC dosyasında eksik/sonsuz değer var.")

    pvar = pd.read_csv(
        pvar_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2, 3, 4],
        names=["CHR", "POS", "ID", "REF", "ALT"],
        dtype={"CHR": str, "POS": np.int64, "ID": str, "REF": str, "ALT": str},
    )
    return pgen_path, pvar, y, pcs


def main():
    parser = argparse.ArgumentParser(description="Bloklu iki aşamalı kromozom RF SNP sıralaması")
    parser.add_argument("--base", required=True)
    parser.add_argument("--chr", required=True, type=int, choices=range(1, 23))
    parser.add_argument("--block-size", type=int, default=10000)
    parser.add_argument("--block-trees", type=int, default=50)
    parser.add_argument("--top-per-block", type=int, default=50)
    parser.add_argument("--final-trees", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=40417)
    parser.add_argument("--min-leaf", type=int, default=5)
    parser.add_argument("--rf-root", default="random_forest_noPC")
    args = parser.parse_args()

    started = time.time()
    chrom = args.chr
    out_dir = os.path.join(args.base, args.rf_root, f"chr{chrom}")
    os.makedirs(out_dir, exist_ok=True)

    pgen_path, pvar, y, pcs = load_inputs(args.base, chrom)
    sample_n = len(y)
    variant_n = len(pvar)
    print(f"chr{chrom}: {sample_n} örnek, {variant_n} SNP", flush=True)
    print(f"Kontrol (PHENO=1): {(y == 0).sum()}", flush=True)
    print(f"Olgu    (PHENO=2): {(y == 1).sum()}", flush=True)

    reader = pgenlib.PgenReader(pgen_path.encode())
    if reader.get_raw_sample_ct() != sample_n:
        raise ValueError("PGEN ve PSAM örnek sayıları farklı.")
    if reader.get_variant_ct() != variant_n:
        raise ValueError("PGEN ve PVAR varyant sayıları farklı.")

    candidate_rows = []
    block_summary = []
    block_no = 0
    try:
        for start in range(0, variant_n, args.block_size):
            block_no += 1
            end = min(start + args.block_size, variant_n)
            block_n = end - start
            g = np.empty((block_n, sample_n), dtype=np.float32)
            reader.read_dosages_range(start, end, g)
            missing_n = impute_missing_inplace(g)
            x = make_x(g, pcs)
            del g

            rf = make_rf(
                args.block_trees,
                args.threads,
                args.seed + block_no,
                args.min_leaf,
                oob=False,
            )
            rf.fit(x, y)
            snp_imp = rf.feature_importances_[:block_n].astype(np.float64)
            pc_imp = rf.feature_importances_[block_n:].astype(np.float64)
            keep_n = min(args.top_per_block, block_n)
            local_order = np.lexsort((np.arange(block_n), -snp_imp))[:keep_n]

            for local_idx in local_order:
                global_idx = start + int(local_idx)
                candidate_rows.append(
                    {
                        "VARIANT_INDEX": global_idx,
                        "BLOCK": block_no,
                        "BLOCK_RF_IMPORTANCE": float(snp_imp[local_idx]),
                    }
                )
            block_summary.append(
                {
                    "BLOCK": block_no,
                    "START_INDEX": start,
                    "END_INDEX_EXCLUSIVE": end,
                    "SNP_COUNT": block_n,
                    "MISSING_DOSAGE": missing_n,
                    "SNP_IMPORTANCE_SUM": float(snp_imp.sum()),
                    "PC_IMPORTANCE_SUM": float(pc_imp.sum()),
                }
            )
            print(
                f"chr{chrom} blok {block_no}: {start}-{end - 1} tamam; "
                f"SNP önem toplamı={snp_imp.sum():.6f}, PC önem toplamı={pc_imp.sum():.6f}",
                flush=True,
            )
            del x, rf, snp_imp, pc_imp

        candidates = pd.DataFrame(candidate_rows).drop_duplicates("VARIANT_INDEX")
        candidate_idx = np.sort(candidates["VARIANT_INDEX"].to_numpy(np.uint32))
        print(f"chr{chrom}: final RF adayı {len(candidate_idx)} SNP", flush=True)

        g_final = np.empty((len(candidate_idx), sample_n), dtype=np.float32)
        reader.read_dosages_list(candidate_idx, g_final)
        final_missing = impute_missing_inplace(g_final)
    finally:
        reader.close()

    x_final = make_x(g_final, pcs)
    del g_final
    final_rf = make_rf(
        args.final_trees,
        args.threads,
        args.seed + 100000 + chrom,
        args.min_leaf,
        oob=True,
    )
    final_rf.fit(x_final, y)
    snp_imp = final_rf.feature_importances_[:len(candidate_idx)].astype(np.float64)
    pc_imp = final_rf.feature_importances_[len(candidate_idx):].astype(np.float64)

    oob_prob = final_rf.oob_decision_function_[:, 1]
    oob_pred = (oob_prob >= 0.5).astype(np.int8)
    oob_auc = roc_auc_score(y, oob_prob)
    oob_bal_acc = balanced_accuracy_score(y, oob_pred)

    candidate_block_imp = candidates.set_index("VARIANT_INDEX")["BLOCK_RF_IMPORTANCE"]
    ranking = pvar.iloc[candidate_idx].copy().reset_index(drop=True)
    ranking.insert(0, "VARIANT_INDEX", candidate_idx.astype(np.int64))
    ranking["RF_IMPORTANCE"] = snp_imp
    ranking["BLOCK_RF_IMPORTANCE"] = candidate_block_imp.loc[candidate_idx].to_numpy()
    ranking = ranking.sort_values(
        ["RF_IMPORTANCE", "BLOCK_RF_IMPORTANCE", "POS"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    ranking.insert(0, "RF_RANK", np.arange(1, len(ranking) + 1))

    ranking_path = os.path.join(out_dir, f"RF_chr{chrom}_ranking.tsv")
    top1000_path = os.path.join(out_dir, f"RF_chr{chrom}_top1000.tsv")
    blocks_path = os.path.join(out_dir, f"RF_chr{chrom}_block_summary.tsv")
    model_path = os.path.join(out_dir, f"RF_chr{chrom}_model.joblib")
    summary_path = os.path.join(out_dir, f"RF_chr{chrom}_summary.json")

    ranking.to_csv(ranking_path, sep="\t", index=False)
    ranking.head(1000).to_csv(top1000_path, sep="\t", index=False)
    pd.DataFrame(block_summary).to_csv(blocks_path, sep="\t", index=False)
    joblib.dump(
        {
            "model": final_rf,
            "variant_indices": candidate_idx,
            "pc_columns": [],
        },
        model_path,
        compress=3,
    )

    summary = {
        "chromosome": chrom,
        "samples": sample_n,
        "controls": int((y == 0).sum()),
        "cases": int((y == 1).sum()),
        "all_snps_screened": variant_n,
        "block_size": args.block_size,
        "block_trees": args.block_trees,
        "top_per_block": args.top_per_block,
        "candidate_snps": int(len(candidate_idx)),
        "final_trees": args.final_trees,
        "min_samples_leaf": args.min_leaf,
        "final_missing_dosages_imputed": final_missing,
        "oob_auc": float(oob_auc),
        "oob_balanced_accuracy": float(oob_bal_acc),
        "final_snp_importance_sum": float(snp_imp.sum()),
        "final_pc_importance_sum": float(pc_imp.sum()),
        "elapsed_minutes": float((time.time() - started) / 60.0),
    }
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    print("\n===== RF TAMAMLANDI =====", flush=True)
    for key, value in summary.items():
        print(f"{key}: {value}", flush=True)
    print(f"Sıralama: {ranking_path}", flush=True)
    print(f"Top 1000: {top1000_path}", flush=True)


if __name__ == "__main__":
    main()
