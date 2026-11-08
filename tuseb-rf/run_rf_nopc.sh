#!/usr/bin/env bash
set -euo pipefail
BASE="${1:-/volume1/homes/tuseb/new_data_to_be_used/GWAS_404_vs_SLD_allchr}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RF_ROOT="random_forest_noPC_$(date +%Y%m%d_%H%M%S)_$$"
python3 -c 'import numpy, pandas, sklearn, joblib, pgenlib'
for CHR in $(seq 1 22); do
  for EXT in pgen pvar psam; do
    P="$BASE/per_chr/chr${CHR}/SLD_EA_chr${CHR}_merged.$EXT"
    [[ -s "$P" ]] || { echo "Eksik/boş: $P" >&2; exit 1; }
  done
done
[[ -s "$BASE/SLD_EA_c2.pheno.tsv" ]] || { echo 'Fenotip dosyası eksik' >&2; exit 1; }
mkdir "$BASE/$RF_ROOT"
cp "$SCRIPT_DIR/rf_rank_chr_nopc.py" "$SCRIPT_DIR/rf_combine_allchr_nopc.py" "$BASE/$RF_ROOT/"
exec > >(tee "$BASE/$RF_ROOT/driver.log") 2>&1
printf 'PC içermeyen RF. Girdi: %s\nÇıktı: %s\n' "$BASE" "$BASE/$RF_ROOT"
for CHR in $(seq 1 22); do
  python3 -u "$BASE/$RF_ROOT/rf_rank_chr_nopc.py" \
    --base "$BASE" --chr "$CHR" --rf-root "$RF_ROOT" \
    --block-size 10000 --block-trees 50 --top-per-block 50 \
    --final-trees 1000 --threads 4 --seed 40417 --min-leaf 5
done
python3 -u "$BASE/$RF_ROOT/rf_combine_allchr_nopc.py" \
  --base "$BASE" --rf-root "$RF_ROOT" \
  --trees 2000 --threads 4 --seed 4042201 --min-leaf 5
printf '\nTamamlandı: %s/combined/RF_allchr_genomewide_ranking.tsv\n' "$BASE/$RF_ROOT"
