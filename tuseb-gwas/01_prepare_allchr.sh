#!/usr/bin/env bash
set -euo pipefail

BASE="/volume1/homes/tuseb/new_data_to_be_used"
OUT="${BASE}/GWAS_404_vs_SLD_allchr"

EA="${BASE}/QC/EA_c2_qc_dedup"
EA_KEEP="${OUT}/EA_c2_1580.keep"
SLD_DIR="${BASE}/imputation_for_SLD_data/output/plink_final"

START_CHR="${1:-1}"
END_CHR="${2:-22}"

THREADS=4

mkdir -p \
  "${OUT}/per_chr" \
  "${OUT}/harmonization" \
  "${OUT}/logs"

LOG="${OUT}/logs/prepare_allchr_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

command -v plink2 >/dev/null
command -v python3 >/dev/null

if [ ! -s "${EA_KEEP}" ]; then
  echo "HATA: EA keep dosyası bulunamadı: ${EA_KEEP}"
  exit 1
fi

count_psam() {
  awk '!/^#/ && NF {n++} END {print n+0}' "$1"
}

count_pvar() {
  awk '!/^#/ && NF {n++} END {print n+0}' "$1"
}

for CHR in $(seq "${START_CHR}" "${END_CHR}"); do
  echo
  echo "========================================"
  echo "CHR${CHR} BAŞLADI: $(date)"
  echo "========================================"

  CHR_OUT="${OUT}/per_chr/chr${CHR}"
  HARM="${OUT}/harmonization/chr${CHR}"

  EA_PREFIX="${CHR_OUT}/EA_chr${CHR}_normalized"
  SLD_PREFIX="${CHR_OUT}/SLD_chr${CHR}_normalized"
  SLD_SOURCE="${SLD_DIR}/SLD_chr${CHR}_qc"

  mkdir -p "${CHR_OUT}" "${HARM}"

  for EXT in pgen pvar psam; do
    if [ ! -f "${SLD_SOURCE}.${EXT}" ]; then
      echo "HATA: chr${CHR} SLD ${EXT} dosyası eksik."
      exit 1
    fi
  done

  echo "chr${CHR}: EA normalize ediliyor."

  plink2 \
    --pfile "${EA}" \
    --chr "${CHR}" \
    --keep "${EA_KEEP}" \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --rm-dup exclude-all \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${EA_PREFIX}"

  echo "chr${CHR}: SLD normalize ediliyor."

  plink2 \
    --pfile "${SLD_SOURCE}" \
    --chr "${CHR}" \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --rm-dup exclude-all \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${SLD_PREFIX}"

  EA_SAMPLES=$(count_psam "${EA_PREFIX}.psam")
  SLD_SAMPLES=$(count_psam "${SLD_PREFIX}.psam")
  EA_VARIANTS=$(count_pvar "${EA_PREFIX}.pvar")
  SLD_VARIANTS=$(count_pvar "${SLD_PREFIX}.pvar")

  if [ "${EA_SAMPLES}" -ne 1580 ]; then
    echo "HATA: chr${CHR} EA örnek sayısı ${EA_SAMPLES}."
    exit 1
  fi

  if [ "${SLD_SAMPLES}" -ne 771 ]; then
    echo "HATA: chr${CHR} SLD örnek sayısı ${SLD_SAMPLES}."
    exit 1
  fi

  echo "EA varyant: ${EA_VARIANTS}"
  echo "SLD varyant: ${SLD_VARIANTS}"

  awk '!/^#/ && NF {print $3}' "${SLD_PREFIX}.pvar" |
    LC_ALL=C sort -u \
    > "${HARM}/SLD_unique.ids"

  : > "${HARM}/EA_exact.ids"
  : > "${HARM}/EA_swapped_rename.tsv"

  awk \
    -v exact="${HARM}/EA_exact.ids" \
    -v swapped="${HARM}/EA_swapped_rename.tsv" '
  BEGIN {OFS="\t"}

  NR==FNR {
    sld[$1]=1
    next
  }

  !/^#/ && NF {
    exact_id=$1 ":" $2 ":" $4 ":" $5
    swapped_id=$1 ":" $2 ":" $5 ":" $4

    if (exact_id in sld) {
      print $3 > exact
    }
    else if (swapped_id in sld) {
      print $3, swapped_id > swapped
    }
  }
  ' "${HARM}/SLD_unique.ids" "${EA_PREFIX}.pvar"

  awk '{print $1}' "${HARM}/EA_swapped_rename.tsv" \
    > "${HARM}/EA_swapped_old.ids"

  awk '{print $2}' "${HARM}/EA_swapped_rename.tsv" \
    > "${HARM}/SLD_swapped.ids"

  cat \
    "${HARM}/EA_exact.ids" \
    "${HARM}/EA_swapped_old.ids" |
    LC_ALL=C sort -u \
    > "${HARM}/EA_all_testable_old.ids"

  cat \
    "${HARM}/EA_exact.ids" \
    "${HARM}/SLD_swapped.ids" |
    LC_ALL=C sort -u \
    > "${HARM}/SLD_all_testable.ids"

  EXACT=$(wc -l < "${HARM}/EA_exact.ids" | tr -d ' ')
  SWAPPED=$(wc -l < "${HARM}/EA_swapped_rename.tsv" | tr -d ' ')
  TESTABLE_EA=$(wc -l < "${HARM}/EA_all_testable_old.ids" | tr -d ' ')
  TESTABLE_SLD=$(wc -l < "${HARM}/SLD_all_testable.ids" | tr -d ' ')

  if [ "${TESTABLE_EA}" -eq 0 ] ||
     [ "${TESTABLE_EA}" -ne "${TESTABLE_SLD}" ]; then
    echo "HATA: chr${CHR} test edilebilir SNP listeleri uyumsuz."
    exit 1
  fi

  echo "Exact: ${EXACT}"
  echo "Swapped: ${SWAPPED}"
  echo "Test edilebilir: ${TESTABLE_EA}"

  awk '
  BEGIN {OFS="\t"}

  NR==FNR {
    keep[$1]=1
    next
  }

  !/^#/ && ($3 in keep) {
    print $3, $4, $5
  }
  ' \
    "${HARM}/SLD_all_testable.ids" \
    "${SLD_PREFIX}.pvar" \
    > "${HARM}/SLD_target_alleles.tsv"

  TARGET_ALLELES=$(wc -l < "${HARM}/SLD_target_alleles.tsv" | tr -d ' ')

  if [ "${TARGET_ALLELES}" -ne "${TESTABLE_EA}" ]; then
    echo "HATA: chr${CHR} target allel tablosu eksik."
    exit 1
  fi

  echo "chr${CHR}: EA eski ID’lerle çıkarılıyor."

  plink2 \
    --pfile "${EA_PREFIX}" \
    --extract "${HARM}/EA_all_testable_old.ids" \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/EA_chr${CHR}_target_oldids"

  echo "chr${CHR}: EA swapped ID’leri değiştiriliyor."

  plink2 \
    --pfile "${CHR_OUT}/EA_chr${CHR}_target_oldids" \
    --update-name "${HARM}/EA_swapped_rename.tsv" 2 1 \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/EA_chr${CHR}_target_renamed"

  RENAMED_COUNT=$(count_pvar \
    "${CHR_OUT}/EA_chr${CHR}_target_renamed.pvar")

  if [ "${RENAMED_COUNT}" -ne "${TESTABLE_EA}" ]; then
    echo "HATA: chr${CHR} rename sonrası SNP kaybı var."
    exit 1
  fi

  echo "chr${CHR}: EA allelleri SLD’ye hizalanıyor."

  plink2 \
    --pfile "${CHR_OUT}/EA_chr${CHR}_target_renamed" \
    --ref-allele force "${HARM}/SLD_target_alleles.tsv" 2 1 \
    --alt1-allele force "${HARM}/SLD_target_alleles.tsv" 3 1 \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/EA_chr${CHR}_target_aligned"

  echo "chr${CHR}: SLD hedef SNP’leri çıkarılıyor."

  plink2 \
    --pfile "${SLD_PREFIX}" \
    --extract "${HARM}/SLD_all_testable.ids" \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/SLD_chr${CHR}_target_raw"

  awk '!/^#/ && NF {print $1,$2,$3,$4,$5}' \
    "${CHR_OUT}/EA_chr${CHR}_target_aligned.pvar" \
    > "${HARM}/EA_target.key"

  awk '!/^#/ && NF {print $1,$2,$3,$4,$5}' \
    "${CHR_OUT}/SLD_chr${CHR}_target_raw.pvar" \
    > "${HARM}/SLD_target.key"

  if ! diff -q \
    "${HARM}/EA_target.key" \
    "${HARM}/SLD_target.key"; then
    echo "HATA: chr${CHR} EA-SLD varyant yapısı uyumsuz."
    exit 1
  fi

  echo "chr${CHR}: SLD doz yönü düzeltiliyor."

  awk '
  BEGIN {OFS="\t"}
  !/^#/ && NF {
    print $3, $5, $4
  }
  ' "${CHR_OUT}/SLD_chr${CHR}_target_raw.pvar" \
    > "${HARM}/SLD_swap_order.tsv"

  awk '
  BEGIN {OFS="\t"}
  !/^#/ && NF {
    print $3, $5, $4, $4, $5
  }
  ' "${CHR_OUT}/SLD_chr${CHR}_target_raw.pvar" \
    > "${HARM}/SLD_relabel.tsv"

  plink2 \
    --pfile "${CHR_OUT}/SLD_chr${CHR}_target_raw" \
    --ref-allele force "${HARM}/SLD_swap_order.tsv" 2 1 \
    --alt1-allele force "${HARM}/SLD_swap_order.tsv" 3 1 \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/SLD_chr${CHR}_order_swapped"

  plink2 \
    --pfile "${CHR_OUT}/SLD_chr${CHR}_order_swapped" \
    --update-alleles "${HARM}/SLD_relabel.tsv" \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/SLD_chr${CHR}_target_corrected"

  awk '!/^#/ && NF {print $1,$2,$3,$4,$5}' \
    "${CHR_OUT}/SLD_chr${CHR}_target_corrected.pvar" \
    > "${HARM}/SLD_corrected_target.key"

  if ! diff -q \
    "${HARM}/EA_target.key" \
    "${HARM}/SLD_corrected_target.key"; then
    echo "HATA: chr${CHR} doz düzeltmesi sonrası varyant yapısı bozuldu."
    exit 1
  fi

  echo "chr${CHR}: frekans doğrulaması yapılıyor."

  plink2 \
    --pfile "${CHR_OUT}/EA_chr${CHR}_target_aligned" \
    --freq \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/EA_chr${CHR}_freq"

  plink2 \
    --pfile "${CHR_OUT}/SLD_chr${CHR}_target_corrected" \
    --freq \
    --threads "${THREADS}" \
    --out "${CHR_OUT}/SLD_chr${CHR}_freq_corrected"

  METRICS=$(python3 - \
    "${CHR_OUT}/SLD_chr${CHR}_freq_corrected.afreq" \
    "${CHR_OUT}/EA_chr${CHR}_freq.afreq" \
    "${TESTABLE_EA}" <<'PY'
import sys
import pandas as pd

def read_freq(path, label):
    df = pd.read_csv(path, sep=r"\s+")

    freq_col = next(
        (
            col for col in
            ["ALT_FREQS", "ALT_FREQ", "ALT1_FREQ"]
            if col in df.columns
        ),
        None
    )

    if freq_col is None:
        candidates = [
            col for col in df.columns
            if "FREQ" in col.upper()
        ]
        if len(candidates) != 1:
            raise RuntimeError(
                f"Frekans sütunu belirlenemedi: {list(df.columns)}"
            )
        freq_col = candidates[0]

    return df[["ID", freq_col]].rename(
        columns={freq_col: label}
    )

sld = read_freq(sys.argv[1], "AF_SLD")
ea = read_freq(sys.argv[2], "AF_EA")
expected = int(sys.argv[3])

merged = sld.merge(
    ea,
    on="ID",
    how="inner",
    validate="one_to_one"
).dropna()

count = len(merged)
corr = merged["AF_SLD"].corr(merged["AF_EA"])
corr_complement = (1-merged["AF_SLD"]).corr(merged["AF_EA"])
median_difference = (
    merged["AF_SLD"]-merged["AF_EA"]
).abs().median()

if count != expected:
    raise RuntimeError(
        f"Frekans eşleşmesi eksik: {count}/{expected}"
    )

if not corr > 0.90:
    raise RuntimeError(
        f"Düzeltilmiş frekans korelasyonu düşük: {corr}"
    )

if not median_difference < 0.10:
    raise RuntimeError(
        f"Medyan frekans farkı yüksek: {median_difference}"
    )

print(
    f"{count}\t"
    f"{corr:.8f}\t"
    f"{corr_complement:.8f}\t"
    f"{median_difference:.8f}"
)
PY
)

  IFS=$'\t' read -r \
    MATCHED CORR CORR_COMPLEMENT MEDIAN_DIFF \
    <<< "${METRICS}"

  printf \
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${CHR}" \
    "${EA_VARIANTS}" \
    "${SLD_VARIANTS}" \
    "${EXACT}" \
    "${SWAPPED}" \
    "${TESTABLE_EA}" \
    "${MATCHED}" \
    "${CORR}" \
    "${CORR_COMPLEMENT}" \
    "${MEDIAN_DIFF}" \
    > "${OUT}/logs/chr${CHR}_prepare_summary.tsv"

  echo "chr${CHR}: BAŞARILI"
  echo "Korelasyon: ${CORR}"
  echo "Medyan AF farkı: ${MEDIAN_DIFF}"
done

MASTER="${OUT}/harmonization/allchr_prepare_summary.tsv"
TEMP_MASTER="${MASTER}.tmp"

printf \
  "CHR\tEA_VARIANTS\tSLD_VARIANTS\tEXACT\tSWAPPED\tTESTABLE\tFREQ_MATCHED\tCORR_CORRECTED\tCORR_COMPLEMENT\tMEDIAN_AF_DIFF\n" \
  > "${TEMP_MASTER}"

for CHR in $(seq 1 22); do
  SUMMARY="${OUT}/logs/chr${CHR}_prepare_summary.tsv"

  if [ -s "${SUMMARY}" ]; then
    cat "${SUMMARY}" >> "${TEMP_MASTER}"
  fi
done

mv "${TEMP_MASTER}" "${MASTER}"

echo
echo "Tamamlanan kromozom özeti:"
column -t -s $'\t' "${MASTER}" || cat "${MASTER}"
echo
echo "BİTTİ: $(date)"
