#!/usr/bin/env bash
set -euo pipefail

OUT="/volume1/homes/tuseb/new_data_to_be_used/GWAS_404_vs_SLD_allchr"
BCFTOOLS="/volume1/homes/tuseb/new_data_to_be_used/tools/envs/bcftools/bin/bcftools"

PHENO="${OUT}/SLD_EA_c2.pheno.tsv"
COVAR="${OUT}/SLD_EA_c2_PC5.covar.tsv"
SUMMARY="${OUT}/harmonization/allchr_prepare_summary.tsv"

GWAS_DIR="${OUT}/gwas_results"
LOG_DIR="${OUT}/logs"

mkdir -p "${GWAS_DIR}" "${LOG_DIR}"

command -v plink2 >/dev/null
[[ -x "${BCFTOOLS}" ]]
[[ -s "${PHENO}" ]]
[[ -s "${COVAR}" ]]
[[ -s "${SUMMARY}" ]]

PHENO_N=$(awk 'NR > 1 && NF {n++} END {print n+0}' "${PHENO}")
COVAR_N=$(awk 'NR > 1 && NF {n++} END {print n+0}' "${COVAR}")

echo "Fenotip örneği: ${PHENO_N}"
echo "Kovaryat örneği: ${COVAR_N}"

[[ "${PHENO_N}" -eq 2351 ]]
[[ "${COVAR_N}" -eq 2351 ]]

for CHR in $(seq 1 22); do
    echo
    echo "===== chr${CHR} BAŞLIYOR ====="
    date

    CHR_DIR="${OUT}/per_chr/chr${CHR}"

    EA="${CHR_DIR}/EA_chr${CHR}_target_aligned"
    SLD="${CHR_DIR}/SLD_chr${CHR}_target_corrected"

    EA_VCF="${CHR_DIR}/EA_chr${CHR}_target_aligned.vcf.gz"
    SLD_VCF="${CHR_DIR}/SLD_chr${CHR}_target_corrected.vcf.gz"
    MERGED_VCF="${CHR_DIR}/SLD_EA_chr${CHR}_merged.vcf.gz"
    MERGED="${CHR_DIR}/SLD_EA_chr${CHR}_merged"

    GWAS="${GWAS_DIR}/GWAS_chr${CHR}_corrected_PC5"
    RESULT="${GWAS}.PHENO.glm.logistic.hybrid"

    DONE="${LOG_DIR}/chr${CHR}_stage2.done"
    CHR_SUMMARY="${LOG_DIR}/chr${CHR}_gwas_summary.tsv"

    EXPECTED=$(awk -F '\t' -v chr="${CHR}" '
        NR > 1 && $1 == chr {print $6}
    ' "${SUMMARY}")

    [[ -n "${EXPECTED}" ]]
    [[ -s "${EA}.pgen" ]]
    [[ -s "${EA}.pvar" ]]
    [[ -s "${EA}.psam" ]]
    [[ -s "${SLD}.pgen" ]]
    [[ -s "${SLD}.pvar" ]]
    [[ -s "${SLD}.psam" ]]

    if [[ -s "${DONE}" && -s "${RESULT}" && -s "${CHR_SUMMARY}" ]]; then
        echo "chr${CHR}: daha önce tamamlanmış, atlanıyor."
        continue
    fi

    echo "Beklenen SNP: ${EXPECTED}"

    plink2 \
        --pfile "${SLD}" \
        --export vcf bgz vcf-dosage=DS-force \
        --out "${CHR_DIR}/SLD_chr${CHR}_target_corrected"

    plink2 \
        --pfile "${EA}" \
        --export vcf bgz vcf-dosage=DS-force \
        --out "${CHR_DIR}/EA_chr${CHR}_target_aligned"

    "${BCFTOOLS}" index -f -t "${SLD_VCF}"
    "${BCFTOOLS}" index -f -t "${EA_VCF}"

    "${BCFTOOLS}" merge \
        -m none \
        -i - \
        --threads 2 \
        -Oz \
        -o "${MERGED_VCF}" \
        "${SLD_VCF}" \
        "${EA_VCF}"

    "${BCFTOOLS}" index -f -t "${MERGED_VCF}"

    plink2 \
        --vcf "${MERGED_VCF}" dosage=DS \
        --make-pgen \
        --sort-vars \
        --threads 4 \
        --out "${MERGED}"

    SAMPLE_N=$(awk '!/^#/ && NF {n++} END {print n+0}' "${MERGED}.psam")
    VARIANT_N=$(awk '!/^#/ && NF {n++} END {print n+0}' "${MERGED}.pvar")

    echo "Birleşmiş örnek: ${SAMPLE_N}"
    echo "Birleşmiş varyant: ${VARIANT_N}"

    [[ "${SAMPLE_N}" -eq 2351 ]]
    [[ "${VARIANT_N}" -eq "${EXPECTED}" ]]

    plink2 \
        --pfile "${MERGED}" \
        --pheno "${PHENO}" \
        --pheno-name PHENO \
        --covar "${COVAR}" \
        --covar-name PC1,PC2,PC3,PC4,PC5 \
        --glm firth-fallback hide-covar \
        --threads 4 \
        --out "${GWAS}"

    [[ -s "${RESULT}" ]]

    read -r ADD_N VALID_N < <(
        awk -F '\t' '
        NR == 1 {
            for (i=1; i<=NF; i++) {
                if ($i == "TEST") test_col=i
                if ($i == "P") p_col=i
            }
            next
        }
        $test_col == "ADD" {
            add++
            if ($p_col != "NA" && $p_col != ".")
                valid++
        }
        END {
            print add+0, valid+0
        }' "${RESULT}"
    )

    echo "ADD satırı: ${ADD_N}"
    echo "Geçerli p-değeri: ${VALID_N}"

    [[ "${ADD_N}" -eq "${EXPECTED}" ]]

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${CHR}" "${SAMPLE_N}" "${VARIANT_N}" "${ADD_N}" "${VALID_N}" \
        > "${CHR_SUMMARY}"

    date > "${DONE}"

    echo "===== chr${CHR} TAMAMLANDI ====="
done

{
    printf 'CHR\tSAMPLES\tVARIANTS\tADD_ROWS\tVALID_P\n'
    for CHR in $(seq 1 22); do
        cat "${LOG_DIR}/chr${CHR}_gwas_summary.tsv"
    done
} > "${GWAS_DIR}/allchr_gwas_summary.tsv"

echo
echo "===== TÜM KROMOZOMLAR TAMAMLANDI ====="
cat "${GWAS_DIR}/allchr_gwas_summary.tsv"
date
