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


def read_ws(path):
    return pd.read_csv(path, sep=r"\s+", dtype=str)


def iid_col(df, label):
    if "#IID" in df.columns:
        return "#IID"
    if "IID" in df.columns:
        return "IID"
    raise ValueError(f"{label}: IID sütunu bulunamadı: {list(df.columns)}")


def impute_missing_inplace(g):
    missing = g < 0
    if not missing.any():
        return 0
    valid = ~missing
    valid_n = valid.sum(axis=1)
    sums = np.where(valid, g, 0.0).sum(axis=1, dtype=np.float64)
    means = np.divide(sums, valid_n, out=np.zeros_like(sums), where=valid_n > 0)
    rows = np.nonzero(missing)[0]
    g[missing] = means[rows]
    return int(missing.sum())


def load_sample_data(base):
    psam_path = os.path.join(base, "per_chr", "chr1", "SLD_EA_chr1_merged.psam")
    pheno_path = os.path.join(base, "SLD_EA_c2.pheno.tsv")

    psam = read_ws(psam_path)
    pheno = read_ws(pheno_path)
    pi = iid_col(psam, "PSAM")
    yi = iid_col(pheno, "PHENO")
    iids = psam[pi].astype(str).reset_index(drop=True)

    if iids.duplicated().any():
        raise ValueError("PSAM içinde tekrarlanan IID var.")
    if pheno[yi].duplicated().any():
        raise ValueError("Fenotip dosyasında tekrarlanan IID var.")

    pheno = pheno.set_index(yi)
    if not iids.isin(pheno.index).all():
        raise ValueError("Fenotipi eksik örnek var.")

    raw_y = pd.to_numeric(pheno.loc[iids, "PHENO"], errors="raise").to_numpy()
    if set(np.unique(raw_y).tolist()) != {1, 2}:
        raise ValueError("PHENO yalnızca 1 ve 2 içermeli.")
    y = (raw_y == 2).astype(np.int8)
    pc_names = []
    pcs = np.empty((len(y), 0), dtype=np.float32)
    if not np.isfinite(pcs).all():
        raise ValueError("PC matrisinde eksik/sonsuz değer var.")
    return iids, y, pcs, pc_names


def main():
    parser = argparse.ArgumentParser(description="Kromozom RF finalistlerini genom çapında yeniden sıralar")
    parser.add_argument("--base", required=True)
    parser.add_argument("--trees", type=int, default=2000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=4042201)
    parser.add_argument("--min-leaf", type=int, default=5)
    parser.add_argument("--rf-root", default="random_forest_noPC")
    args = parser.parse_args()

    started = time.time()
    out_dir = os.path.join(args.base, args.rf_root, "combined")
    os.makedirs(out_dir, exist_ok=True)
    iids, y, pcs, pc_names = load_sample_data(args.base)
    sample_n = len(iids)

    per_chr = []
    total_candidates = 0
    total_original_snps = 0
    for chrom in range(1, 23):
        ranking_path = os.path.join(
            args.base, args.rf_root, f"chr{chrom}", f"RF_chr{chrom}_ranking.tsv"
        )
        if not os.path.isfile(ranking_path):
            raise FileNotFoundError(ranking_path)
        chr_summary_path = os.path.join(
            args.base, args.rf_root, f"chr{chrom}", f"RF_chr{chrom}_summary.json"
        )
        if not os.path.isfile(chr_summary_path):
            raise FileNotFoundError(chr_summary_path)
        with open(chr_summary_path, "r", encoding="utf-8") as handle:
            chr_summary = json.load(handle)
        total_original_snps += int(chr_summary["all_snps_screened"])
        rank = pd.read_csv(ranking_path, sep="\t")
        required = {
            "VARIANT_INDEX", "CHR", "POS", "ID", "REF", "ALT",
            "RF_IMPORTANCE", "BLOCK_RF_IMPORTANCE"
        }
        missing = required.difference(rank.columns)
        if missing:
            raise ValueError(f"chr{chrom} sıralamasında eksik sütunlar: {sorted(missing)}")
        if rank["VARIANT_INDEX"].duplicated().any():
            raise ValueError(f"chr{chrom} içinde tekrarlanan VARIANT_INDEX var.")
        rank = rank.sort_values("VARIANT_INDEX", kind="mergesort").reset_index(drop=True)
        per_chr.append((chrom, rank))
        total_candidates += len(rank)
        print(f"chr{chrom}: {len(rank)} blok finalisti", flush=True)

    print(f"Toplam genomik RF adayı: {total_candidates}", flush=True)
    estimated_mb = sample_n * (total_candidates + len(pc_names)) * 4 / (1024 ** 2)
    print(f"Ana X matrisi yaklaşık: {estimated_mb:.1f} MiB", flush=True)

    x = np.empty((sample_n, total_candidates + len(pc_names)), dtype=np.float32)
    metadata_parts = []
    offset = 0
    total_missing = 0

    reference_iids = iids.to_numpy()
    for chrom, rank in per_chr:
        prefix = os.path.join(
            args.base, "per_chr", f"chr{chrom}", f"SLD_EA_chr{chrom}_merged"
        )
        pgen_path = prefix + ".pgen"
        psam_path = prefix + ".psam"
        psam = read_ws(psam_path)
        current_iids = psam[iid_col(psam, f"chr{chrom} PSAM")].astype(str).to_numpy()
        if not np.array_equal(current_iids, reference_iids):
            raise ValueError(f"chr{chrom} örnek sırası chr1 ile aynı değil.")

        variant_idx = np.ascontiguousarray(rank["VARIANT_INDEX"].to_numpy(np.uint32))
        g = np.empty((len(variant_idx), sample_n), dtype=np.float32)
        reader = pgenlib.PgenReader(pgen_path.encode())
        try:
            if reader.get_raw_sample_ct() != sample_n:
                raise ValueError(f"chr{chrom} PGEN örnek sayısı uyumsuz.")
            if len(variant_idx) and int(variant_idx[-1]) >= reader.get_variant_ct():
                raise ValueError(f"chr{chrom} VARIANT_INDEX PGEN sınırını aşıyor.")
            reader.read_dosages_list(variant_idx, g)
        finally:
            reader.close()

        missing_n = impute_missing_inplace(g)
        total_missing += missing_n
        end = offset + len(rank)
        x[:, offset:end] = g.T
        del g

        meta = rank[[
            "VARIANT_INDEX", "CHR", "POS", "ID", "REF", "ALT",
            "RF_IMPORTANCE", "BLOCK_RF_IMPORTANCE"
        ]].copy()
        meta = meta.rename(columns={"RF_IMPORTANCE": "CHR_RF_IMPORTANCE"})
        metadata_parts.append(meta)
        offset = end
        print(f"chr{chrom}: dozajlar yüklendi; eksik={missing_n}", flush=True)

    x[:, total_candidates:] = pcs
    metadata = pd.concat(metadata_parts, ignore_index=True)
    if len(metadata) != total_candidates or offset != total_candidates:
        raise RuntimeError("Aday veya matris sütun sayısı uyuşmuyor.")

    print(f"Genom çapında RF başlıyor: {args.trees} ağaç", flush=True)
    rf = RandomForestClassifier(
        n_estimators=args.trees,
        criterion="gini",
        max_features="sqrt",
        max_depth=None,
        min_samples_leaf=args.min_leaf,
        class_weight="balanced_subsample",
        bootstrap=True,
        max_samples=0.80,
        n_jobs=args.threads,
        random_state=args.seed,
    )
    rf.fit(x, y)

    snp_imp = rf.feature_importances_[:total_candidates].astype(np.float64)
    pc_imp = rf.feature_importances_[total_candidates:].astype(np.float64)
    feature_order = metadata[["CHR", "POS", "ID", "REF", "ALT"]].copy()
    metadata["GENOME_RF_IMPORTANCE"] = snp_imp
    metadata = metadata.sort_values(
        ["GENOME_RF_IMPORTANCE", "CHR_RF_IMPORTANCE", "BLOCK_RF_IMPORTANCE", "CHR", "POS"],
        ascending=[False, False, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    metadata.insert(0, "GENOME_RF_RANK", np.arange(1, len(metadata) + 1))

    ranking_path = os.path.join(out_dir, "RF_allchr_genomewide_ranking.tsv")
    top1000_path = os.path.join(out_dir, "RF_allchr_genomewide_top1000.tsv")
    model_path = os.path.join(out_dir, "RF_allchr_genomewide_model.joblib")
    summary_path = os.path.join(out_dir, "RF_allchr_genomewide_summary.json")

    metadata.to_csv(ranking_path, sep="\t", index=False)
    metadata.head(1000).to_csv(top1000_path, sep="\t", index=False)
    joblib.dump(
        {
            "model": rf,
            "features_in_model_order": feature_order,
            "pc_columns": pc_names,
        },
        model_path,
        compress=3,
    )

    summary = {
        "samples": sample_n,
        "controls": int((y == 0).sum()),
        "cases": int((y == 1).sum()),
        "all_original_snps_screened": total_original_snps,
        "genome_rf_candidate_snps": total_candidates,
        "trees": args.trees,
        "min_samples_leaf": args.min_leaf,
        "missing_dosages_imputed": total_missing,
        "genome_snp_importance_sum": float(snp_imp.sum()),
        "pc_importance_sum": float(pc_imp.sum()),
        "pc_importances": {pc_names[i]: float(pc_imp[i]) for i in range(len(pc_names))},
        "elapsed_minutes": float((time.time() - started) / 60.0),
    }
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    print("\n===== GENOM ÇAPINDA RF TAMAMLANDI =====", flush=True)
    for key, value in summary.items():
        print(f"{key}: {value}", flush=True)
    print(f"Tüm sıralama: {ranking_path}", flush=True)
    print(f"Top 1000: {top1000_path}", flush=True)


if __name__ == "__main__":
    main()
