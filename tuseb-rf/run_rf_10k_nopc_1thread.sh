#!/usr/bin/env bash
set -euo pipefail

SOURCE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE="${1:-$PWD}"
BASE="$(cd -- "${BASE}" && pwd)"
RF_ROOT="RF_noPC_10k_$(date +%Y%m%d_%H%M%S)_$$"
OUT="${BASE}/${RF_ROOT}"

# Sayısal kütüphanelerin paralel bellek kullanımını sınırla.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Başlamadan önce bağımlılıkları ve dosyaları kontrol et.
python3 -c 'import numpy, pandas, sklearn, joblib, pgenlib'

for SCRIPT in rf_rank_chr_nopc.py rf_combine_allchr_nopc.py; do
    [[ -s "${SOURCE}/${SCRIPT}" ]] || {
        echo "Eksik betik: ${SOURCE}/${SCRIPT}" >&2
        exit 1
    }
done

[[ -s "${BASE}/SLD_EA_c2.pheno.tsv" ]] || {
    echo "Fenotip dosyası eksik." >&2
    exit 1
}

for CHR in $(seq 1 22); do
    for EXT in pgen pvar psam; do
        FILE="${BASE}/per_chr/chr${CHR}/SLD_EA_chr${CHR}_merged.${EXT}"
        [[ -s "${FILE}" ]] || {
            echo "Eksik genotip dosyası: ${FILE}" >&2
            exit 1
        }
    done
done

mkdir "${OUT}"
cp "${SOURCE}/rf_rank_chr_nopc.py" "${OUT}/"
cp "${SOURCE}/rf_combine_allchr_nopc.py" "${OUT}/"

printf 'Sonuç klasörü: %s\n' "${OUT}"
printf 'PC yok | Blok: 10000 SNP | İş parçacığı: 1\n'

for CHR in $(seq 1 22); do
    printf '\n===== chr%s BAŞLIYOR =====\n' "${CHR}"
    date

    python3 -u "${OUT}/rf_rank_chr_nopc.py" \
        --base "${BASE}" \
        --chr "${CHR}" \
        --rf-root "${RF_ROOT}" \
        --block-size 10000 \
        --block-trees 50 \
        --top-per-block 50 \
        --final-trees 1000 \
        --threads 1 \
        --seed 40417 \
        --min-leaf 5
done

printf '\n===== GENOM ÇAPINDA SON RF =====\n'

python3 -u "${OUT}/rf_combine_allchr_nopc.py" \
    --base "${BASE}" \
    --rf-root "${RF_ROOT}" \
    --trees 2000 \
    --threads 1 \
    --seed 4042201 \
    --min-leaf 5

printf '\nTAMAMLANDI\n'
printf 'Tüm sıralama: %s/combined/RF_allchr_genomewide_ranking.tsv\n' "${OUT}"
printf 'İlk 1000 SNP: %s/combined/RF_allchr_genomewide_top1000.tsv\n' "${OUT}"
