#!/usr/bin/env python3

import argparse
import gzip
import os
import re

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo


def parse_attributes(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def load_genes(gtf_path):
    records = {}
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chrom = fields[0][3:] if fields[0].startswith("chr") else fields[0]
            if chrom == "MT":
                chrom = "M"
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            records.setdefault(chrom, []).append(
                {
                    "START": int(fields[3]),
                    "END": int(fields[4]),
                    "GENE": attrs.get("gene_name", gene_id),
                    "GENE_BIOTYPE": attrs.get(
                        "gene_biotype", attrs.get("gene_type", "")
                    ),
                    "ENSEMBL_GENE_ID": gene_id,
                }
            )

    genes = {}
    for chrom, chrom_records in records.items():
        genes[chrom] = (
            pd.DataFrame(chrom_records)
            .drop_duplicates(["START", "END", "ENSEMBL_GENE_ID"])
            .sort_values(["START", "END", "GENE"])
            .reset_index(drop=True)
        )
    return genes


def annotate_position(chrom, pos, genes):
    chrom = str(chrom)
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    chrom = re.sub(r"\.0$", "", chrom)
    pos = int(pos)
    if chrom not in genes or genes[chrom].empty:
        return {
            "GENE": "",
            "GENE_RELATION": "unmapped",
            "DISTANCE_BP": None,
            "GENE_BIOTYPE": "",
            "ENSEMBL_GENE_ID": "",
        }

    chrom_genes = genes[chrom]
    starts = chrom_genes["START"].to_numpy(np.int64)
    ends = chrom_genes["END"].to_numpy(np.int64)
    overlap = (starts <= pos) & (ends >= pos)

    if overlap.any():
        selected = chrom_genes.loc[overlap].copy()
        relation = "overlap"
        distance = 0
    else:
        distances = np.where(
            pos < starts, starts - pos, np.where(pos > ends, pos - ends, 0)
        )
        distance = int(distances.min())
        selected = chrom_genes.loc[distances == distance].copy()
        relation = "nearest"

    selected = (
        selected.drop_duplicates("ENSEMBL_GENE_ID")
        .sort_values(["GENE", "ENSEMBL_GENE_ID"])
    )
    return {
        "GENE": "; ".join(selected["GENE"].fillna("").astype(str).unique()),
        "GENE_RELATION": relation,
        "DISTANCE_BP": distance,
        "GENE_BIOTYPE": "; ".join(
            selected["GENE_BIOTYPE"].fillna("").astype(str).unique()
        ),
        "ENSEMBL_GENE_ID": "; ".join(
            selected["ENSEMBL_GENE_ID"].fillna("").astype(str).unique()
        ),
    }


def annotate_dataframe(df, chr_col, pos_col, genes):
    annotations = [
        annotate_position(chrom, pos, genes)
        for chrom, pos in zip(df[chr_col], df[pos_col])
    ]
    return pd.concat(
        [df.reset_index(drop=True), pd.DataFrame(annotations)], axis=1
    )


def python_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_excel(df, path, sheet_name, title, table_name):
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = sheet_name
    sheet.sheet_view.showGridLines = False

    last_col = get_column_letter(len(df.columns))
    sheet.merge_cells(f"A1:{last_col}1")
    sheet["A1"] = title
    sheet["A1"].font = Font(name="Arial", size=14, bold=True, color="1F1F1F")
    sheet["A1"].alignment = Alignment(horizontal="left", vertical="center")
    sheet.row_dimensions[1].height = 24

    sheet.merge_cells(f"A2:{last_col}2")
    sheet["A2"] = (
        "Gen anotasyonu: Ensembl GRCh37 release 87 (hg19). "
        "Gen dışındaki SNP'ler en yakın gene atanmıştır."
    )
    sheet["A2"].font = Font(
        name="Arial", size=9, italic=True, color="666666"
    )

    header_row = 4
    for col_idx, col_name in enumerate(df.columns, 1):
        cell = sheet.cell(header_row, col_idx, col_name)
        cell.font = Font(name="Arial", size=10, bold=True, color="FFFFFF")
        cell.fill = PatternFill(fill_type="solid", fgColor="1F4E78")
        cell.alignment = Alignment(horizontal="center", vertical="center")

    for row_idx, row in enumerate(df.itertuples(index=False, name=None), 5):
        for col_idx, value in enumerate(row, 1):
            value = python_value(value)
            cell = sheet.cell(row_idx, col_idx, value)
            cell.font = Font(name="Arial", size=10)
            cell.alignment = Alignment(
                horizontal="left" if isinstance(value, str) else "right",
                vertical="center",
            )

    last_row = header_row + len(df)
    table = Table(displayName=table_name, ref=f"A4:{last_col}{last_row}")
    table.tableStyleInfo = TableStyleInfo(
        name="TableStyleMedium2",
        showFirstColumn=False,
        showLastColumn=False,
        showRowStripes=True,
        showColumnStripes=False,
    )
    sheet.add_table(table)
    sheet.freeze_panes = "A5"

    scientific = {
        "P_VALUE",
        "BONFERRONI",
        "FDR_BH",
        "GENOME_RF_IMPORTANCE",
        "CHR_RF_IMPORTANCE",
        "BLOCK_RF_IMPORTANCE",
    }
    integers = {"GWAS_RANK", "GENOME_RF_RANK", "CHR", "POS", "DISTANCE_BP"}
    widths = {
        "GWAS_RANK": 12,
        "GENOME_RF_RANK": 16,
        "CHR": 8,
        "POS": 14,
        "ID": 24,
        "REF": 8,
        "ALT": 8,
        "A1": 8,
        "OR": 14,
        "P_VALUE": 14,
        "BONFERRONI": 14,
        "FDR_BH": 14,
        "GENOME_RF_IMPORTANCE": 22,
        "CHR_RF_IMPORTANCE": 20,
        "BLOCK_RF_IMPORTANCE": 22,
        "GENE": 30,
        "GENE_RELATION": 15,
        "DISTANCE_BP": 15,
        "GENE_BIOTYPE": 24,
        "ENSEMBL_GENE_ID": 26,
    }

    for col_idx, col_name in enumerate(df.columns, 1):
        letter = get_column_letter(col_idx)
        sheet.column_dimensions[letter].width = widths.get(col_name, 14)
        if col_name in scientific:
            fmt = "0.000E+00"
        elif col_name in integers:
            fmt = "0"
        elif col_name == "OR":
            fmt = "0.0000E+00"
        else:
            continue
        for row_idx in range(5, last_row + 1):
            sheet.cell(row_idx, col_idx).number_format = fmt

    workbook.save(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base",
        default="/volume1/homes/tuseb/new_data_to_be_used/GWAS_404_vs_SLD_allchr",
    )
    args = parser.parse_args()
    base = args.base

    gwas_file = os.path.join(
        base, "combined_results", "GWAS_allchr_corrected_PC5_top1000_basic.tsv"
    )
    overlap_file = os.path.join(
        base,
        "random_forest",
        "combined",
        "GWAS_RF_top1000_overlap.tsv",
    )
    gtf_file = os.path.join(
        base, "reference", "Homo_sapiens.GRCh37.87.gtf.gz"
    )
    output_dir = os.path.join(base, "excel")
    os.makedirs(output_dir, exist_ok=True)

    for path in (gwas_file, overlap_file, gtf_file):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    print("GRCh37 gen anotasyonları yükleniyor...", flush=True)
    genes = load_genes(gtf_file)

    gwas = pd.read_csv(gwas_file, sep="\t")
    if len(gwas) != 1000:
        raise ValueError(f"GWAS satırı 1000 yerine {len(gwas)}")
    gwas.insert(0, "GWAS_RANK", np.arange(1, 1001))
    gwas["_KEY"] = (
        gwas["CHR"].astype(str)
        + ":"
        + gwas["POS"].astype(str)
        + ":"
        + gwas["REF"].astype(str)
        + ":"
        + gwas["ALT"].astype(str)
    )
    gwas = annotate_dataframe(gwas, "CHR", "POS", genes)
    gwas_cols = [
        "GWAS_RANK", "CHR", "POS", "ID", "REF", "ALT", "A1", "OR",
        "P_VALUE", "BONFERRONI", "FDR_BH", "GENE", "GENE_RELATION",
        "DISTANCE_BP", "GENE_BIOTYPE", "ENSEMBL_GENE_ID",
    ]
    gwas_out = gwas[gwas_cols].copy()

    overlap = pd.read_csv(overlap_file, sep="\t")
    if len(overlap) != 22:
        raise ValueError(f"Ortak SNP satırı 22 yerine {len(overlap)}")
    rank_map = dict(zip(gwas["_KEY"], gwas["GWAS_RANK"]))
    overlap["GWAS_RANK"] = overlap["_VARIANT_KEY"].map(rank_map)
    if overlap["GWAS_RANK"].isna().any():
        raise ValueError("Bazı ortak SNP'ler GWAS Top 1000 içinde bulunamadı")

    overlap_basic = pd.DataFrame(
        {
            "GWAS_RANK": overlap["GWAS_RANK"].astype(int),
            "GENOME_RF_RANK": overlap["GENOME_RF_RANK"],
            "CHR": overlap["CHR_GWAS"],
            "POS": overlap["POS_GWAS"],
            "ID": overlap["ID_GWAS"],
            "REF": overlap["REF_GWAS"],
            "ALT": overlap["ALT_GWAS"],
            "A1": overlap["A1"],
            "OR": overlap["OR"],
            "P_VALUE": overlap["P_VALUE"],
            "BONFERRONI": overlap["BONFERRONI"],
            "FDR_BH": overlap["FDR_BH"],
            "GENOME_RF_IMPORTANCE": overlap["GENOME_RF_IMPORTANCE"],
            "CHR_RF_IMPORTANCE": overlap["CHR_RF_IMPORTANCE"],
            "BLOCK_RF_IMPORTANCE": overlap["BLOCK_RF_IMPORTANCE"],
        }
    ).sort_values(["GWAS_RANK", "GENOME_RF_RANK"])
    overlap_out = annotate_dataframe(overlap_basic, "CHR", "POS", genes)

    gwas_xlsx = os.path.join(output_dir, "GWAS_Top1000_hg19_Gene.xlsx")
    overlap_xlsx = os.path.join(output_dir, "GWAS_RF_Ortak22_hg19_Gene.xlsx")
    write_excel(
        gwas_out,
        gwas_xlsx,
        "GWAS Top 1000",
        "GWAS Top 1000 SNP ve hg19 gen anotasyonu",
        "GWASTop1000Table",
    )
    write_excel(
        overlap_out,
        overlap_xlsx,
        "Ortak 22 SNP",
        "GWAS ve Random Forest Top 1000 ortak SNP'leri",
        "Overlap22Table",
    )

    print("\n===== TAMAMLANDI =====")
    print("GWAS satırı:", len(gwas_out))
    print("Ortak SNP satırı:", len(overlap_out))
    print("GWAS Excel:", gwas_xlsx)
    print("Ortak 22 Excel:", overlap_xlsx)


if __name__ == "__main__":
    main()
