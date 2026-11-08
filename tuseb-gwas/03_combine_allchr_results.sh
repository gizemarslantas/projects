#!/usr/bin/env bash
set -euo pipefail

OUT="/volume1/homes/tuseb/new_data_to_be_used/GWAS_404_vs_SLD_allchr"
GWAS_DIR="${OUT}/gwas_results"
COMBINED_DIR="${OUT}/combined_results"

mkdir -p "${COMBINED_DIR}"

UNSORTED="${COMBINED_DIR}/allchr_ADD_valid.unsorted.tsv"
SORTED="${COMBINED_DIR}/allchr_ADD_valid.sorted_by_P.tsv"
COUNT_FILE="${COMBINED_DIR}/valid_test_count.txt"

ALL_GZ="${COMBINED_DIR}/GWAS_allchr_corrected_PC5_all_results.tsv.gz"
TOP1000="${COMBINED_DIR}/GWAS_allchr_corrected_PC5_top1000_basic.tsv"
STATS="${COMBINED_DIR}/GWAS_allchr_corrected_PC5_summary.tsv"

echo "===== Sonuçlar birleştiriliyor ====="

python3 - "${GWAS_DIR}" "${UNSORTED}" "${COUNT_FILE}" <<'PYTHON'
import math
import os
import sys

gwas_dir, output_path, count_path = sys.argv[1:4]

columns = [
    "CHROM", "POS", "ID", "REF", "ALT", "A1",
    "OR", "LOG(OR)_SE", "Z_STAT", "P",
    "FIRTH?", "OBS_CT", "ERRCODE"
]

total = 0

with open(output_path, "w", encoding="utf-8") as out:
    for chrom in range(1, 23):
        path = os.path.join(
            gwas_dir,
            f"GWAS_chr{chrom}_corrected_PC5.PHENO.glm.logistic.hybrid"
        )

        if not os.path.isfile(path):
            raise FileNotFoundError(path)

        with open(path, "r", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            header = [name.lstrip("#") for name in header]
            index = {name: i for i, name in enumerate(header)}

            required = [
                "CHROM", "POS", "ID", "REF", "ALT",
                "A1", "OR", "P", "TEST"
            ]
            missing = [name for name in required if name not in index]

            if missing:
                raise RuntimeError(
                    f"{path}: eksik sütunlar: {','.join(missing)}"
                )

            for line in handle:
                fields = line.rstrip("\n").split("\t")

                if fields[index["TEST"]] != "ADD":
                    continue

                p_text = fields[index["P"]]

                try:
                    p_value = float(p_text)
                except (TypeError, ValueError):
                    continue

                if not math.isfinite(p_value) or not 0 <= p_value <= 1:
                    continue

                row = []
                for name in columns:
                    if name in index and index[name] < len(fields):
                        row.append(fields[index[name]])
                    else:
                        row.append(".")

                out.write("\t".join(row) + "\n")
                total += 1

        print(f"chr{chrom}: tamamlandı", flush=True)

with open(count_path, "w", encoding="utf-8") as handle:
    handle.write(str(total) + "\n")

print(f"Toplam geçerli test: {total}")
PYTHON

M=$(cat "${COUNT_FILE}")

if [[ "${M}" -ne 6003487 ]]; then
    echo "HATA: Beklenen 6003487 geçerli test, bulunan ${M}."
    exit 1
fi

echo "===== ${M} p-değeri sıralanıyor ====="

LC_ALL=C sort \
    -t $'\t' \
    -k10,10g \
    "${UNSORTED}" \
    > "${SORTED}"

SORTED_N=$(wc -l < "${SORTED}")
SORTED_N=$(echo "${SORTED_N}" | tr -d ' ')

if [[ "${SORTED_N}" -ne "${M}" ]]; then
    echo "HATA: Sıralama sonrasında satır sayısı değişti."
    exit 1
fi

echo "===== Global Bonferroni, FDR ve lambda hesaplanıyor ====="

python3 - \
    "${SORTED}" \
    "${M}" \
    "${ALL_GZ}" \
    "${TOP1000}" \
    "${STATS}" <<'PYTHON'
import gzip
import math
import sys
from array import array
from statistics import NormalDist

sorted_path, m_text, all_gz, top1000_path, stats_path = sys.argv[1:6]
m = int(m_text)

input_columns = [
    "CHROM", "POS", "ID", "REF", "ALT", "A1",
    "OR", "LOG(OR)_SE", "Z_STAT", "P",
    "FIRTH", "OBS_CT", "ERRCODE"
]

output_columns = input_columns + [
    "RANK", "BONFERRONI", "FDR_BH"
]

raw_q = array("d")
median_ranks = {(m + 1) // 2} if m % 2 else {m // 2, m // 2 + 1}
median_p_values = []

observed = 0

with open(sorted_path, "r", encoding="utf-8") as handle:
    for rank, line in enumerate(handle, start=1):
        fields = line.rstrip("\n").split("\t")
        p_value = float(fields[9])

        raw_q.append(min(1.0, p_value * m / rank))

        if rank in median_ranks:
            median_p_values.append(p_value)

        observed += 1

if observed != m:
    raise RuntimeError(f"Beklenen {m} satır, okunan {observed}")

running_min = 1.0

for i in range(m - 1, -1, -1):
    if raw_q[i] < running_min:
        running_min = raw_q[i]
    raw_q[i] = running_min

bonf_threshold = 0.05 / m
p5_count = 0
bonf_count = 0
fdr_count = 0

top_columns = [
    "CHR", "POS", "ID", "REF", "ALT", "A1",
    "OR", "P_VALUE", "BONFERRONI", "FDR_BH",
    "LOG_OR_SE", "Z_STAT", "FIRTH", "OBS_CT", "ERRCODE"
]

with open(sorted_path, "r", encoding="utf-8") as source, \
     gzip.open(all_gz, "wt", encoding="utf-8", compresslevel=6) as all_out, \
     open(top1000_path, "w", encoding="utf-8") as top_out:

    all_out.write("\t".join(output_columns) + "\n")
    top_out.write("\t".join(top_columns) + "\n")

    for rank, line in enumerate(source, start=1):
        fields = line.rstrip("\n").split("\t")
        p_value = float(fields[9])
        bonferroni = min(1.0, p_value * m)
        fdr_bh = raw_q[rank - 1]

        if p_value < 5e-8:
            p5_count += 1
        if p_value < bonf_threshold:
            bonf_count += 1
        if fdr_bh < 0.05:
            fdr_count += 1

        all_out.write(
            "\t".join(
                fields + [
                    str(rank),
                    f"{bonferroni:.12g}",
                    f"{fdr_bh:.12g}"
                ]
            ) + "\n"
        )

        if rank <= 1000:
            top_row = [
                fields[0],   # CHROM
                fields[1],   # POS
                fields[2],   # ID
                fields[3],   # REF
                fields[4],   # ALT
                fields[5],   # A1
                fields[6],   # OR
                fields[9],   # P
                f"{bonferroni:.12g}",
                f"{fdr_bh:.12g}",
                fields[7],   # LOG(OR)_SE
                fields[8],   # Z_STAT
                fields[10],  # FIRTH
                fields[11],  # OBS_CT
                fields[12],  # ERRCODE
            ]
            top_out.write("\t".join(top_row) + "\n")

median_chi_values = []

for p_value in median_p_values:
    z_value = NormalDist().inv_cdf(1.0 - p_value / 2.0)
    median_chi_values.append(z_value * z_value)

median_chi_square = sum(median_chi_values) / len(median_chi_values)
lambda_gc = median_chi_square / 0.4549364231195727

with open(stats_path, "w", encoding="utf-8") as stats:
    stats.write("METRIC\tVALUE\n")
    stats.write(f"VALID_TESTS\t{m}\n")
    stats.write(f"BONFERRONI_THRESHOLD\t{bonf_threshold:.12g}\n")
    stats.write("GENOME_WIDE_THRESHOLD\t5e-8\n")
    stats.write(f"LAMBDA_GC\t{lambda_gc:.6f}\n")
    stats.write(f"P_LT_5E-8\t{p5_count}\n")
    stats.write(f"P_LT_BONFERRONI\t{bonf_count}\n")
    stats.write(f"FDR_LT_0.05\t{fdr_count}\n")

print(f"Geçerli test: {m}")
print(f"Bonferroni eşiği: {bonf_threshold:.12g}")
print(f"Lambda GC: {lambda_gc:.6f}")
print(f"P < 5e-8: {p5_count}")
print(f"P < Bonferroni: {bonf_count}")
print(f"FDR < 0.05: {fdr_count}")
PYTHON

echo
echo "===== ANALİZ TAMAMLANDI ====="
cat "${STATS}"

echo
echo "Tüm sonuçlar: ${ALL_GZ}"
echo "Top 1000: ${TOP1000}"
echo "Özet: ${STATS}"
