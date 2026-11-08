#!/usr/bin/env python3
from pathlib import Path
import argparse
from datetime import datetime
import numpy as np
import pandas as pd

try:
    import openpyxl
    from openpyxl.styles import Font, PatternFill
    from openpyxl.worksheet.table import Table, TableStyleInfo
except ImportError:
    raise SystemExit(
        "openpyxl eksik. Önce şu komutu çalıştırıp tekrar dene:\n"
        "python3 -m pip install --user openpyxl"
    )

parser = argparse.ArgumentParser(description="RF importance ortalaması üzerindeki adayları Excel'e aktarır.")
parser.add_argument("--combined", type=Path, default=Path("RF_noPC_10k_20261108_181230_4597/combined"), help="RF son model sonuç klasörü")
args = parser.parse_args()
folder = args.combined.resolve()
source = folder / "RF_allchr_genomewide_ranking.tsv"

if not source.is_file():
    raise SystemExit(f"Sonuç dosyası bulunamadı: {source}")

df = pd.read_csv(
    source,
    sep="\t",
    dtype={"ID": str, "CHR": str, "REF": str, "ALT": str}
)

column = "GENOME_RF_IMPORTANCE"
if column not in df.columns:
    raise SystemExit(f"Gerekli sütun bulunamadı: {column}")

importance = pd.to_numeric(df[column], errors="raise")
if not np.isfinite(importance).all() or (importance < 0).any():
    raise SystemExit("Importance sütununda geçersiz değer var.")

if len(df) != 30550:
    raise SystemExit(
        f"Beklenen aday sayısı 30550; dosyada {len(df)} satır var."
    )

total = float(importance.sum())
if not np.isclose(total, 1.0, rtol=0, atol=1e-6):
    raise SystemExit(f"Importance toplamı 1 değil: {total}")

df[column] = importance

# Yuvarlanmış sayı yerine dosyadaki gerçek ortalamayı kullan.
threshold = float(importance.mean())
selected = df.loc[importance > threshold].copy()
selected = selected.sort_values(
    [column, "GENOME_RF_RANK"],
    ascending=[False, True]
)

summary = pd.DataFrame([
    ("Kaynak dosya", source.name),
    ("Filtre sütunu", column),
    ("Eşik yöntemi", "Son modeldeki adayların ortalama importance değeri"),
    ("Eşik değeri", threshold),
    ("Karşılaştırma", "Importance > eşik"),
    ("Son modeldeki aday SNP sayısı", len(df)),
    ("Eşik üzerindeki SNP sayısı", len(selected)),
    ("Importance toplamı", total),
    ("PC kullanımı", "Yok"),
    ("Yorum", "Aday önceliklendirme eşiğidir; istatistiksel anlamlılık değildir."),
], columns=["Alan", "Değer"])

stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
output = folder / f"RF_ortalama_ustu_SNPler_{stamp}.xlsx"

with pd.ExcelWriter(output, engine="openpyxl") as writer:
    selected.to_excel(writer, sheet_name="SNPler", index=False)
    summary.to_excel(writer, sheet_name="Ozet", index=False)

    for ws in writer.book.worksheets:
        ws.freeze_panes = "A2"
        for cell in ws[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="17365D")
        ws.row_dimensions[1].height = 28

    ws = writer.book["SNPler"]
    for idx, name in enumerate(selected.columns, 1):
        letter = openpyxl.utils.get_column_letter(idx)
        width = 28 if name == "ID" or "IMPORTANCE" in name else 18
        ws.column_dimensions[letter].width = max(width, len(name) + 2)

        for cells in ws.iter_rows(min_row=2, min_col=idx, max_col=idx):
            cell = cells[0]
            if "IMPORTANCE" in name:
                cell.number_format = "0.0000000000"
            elif name in {"POS", "VARIANT_INDEX", "GENOME_RF_RANK"}:
                cell.number_format = "0"
            elif name in {"ID", "CHR", "REF", "ALT"}:
                cell.number_format = "@"

    if len(selected):
        table = Table(displayName="RFEsikUstu", ref=ws.dimensions)
        table.tableStyleInfo = TableStyleInfo(
            name="TableStyleMedium2",
            showRowStripes=True
        )
        ws.add_table(table)

    overview = writer.book["Ozet"]
    overview.column_dimensions["A"].width = 38
    overview.column_dimensions["B"].width = 100
    overview["B5"].number_format = "0.000000000000"

# Kaydedilen dosyanın satır sayısını kontrol et.
check = openpyxl.load_workbook(output, read_only=True, data_only=True)
saved_count = check["SNPler"].max_row - 1
check.close()
if saved_count != len(selected):
    raise RuntimeError("Excel satır sayısı doğrulanamadı.")

print(f"Eşik: {threshold:.12f}")
print(f"Toplam aday SNP: {len(df)}")
print(f"Eşik üzerindeki SNP: {len(selected)}")
print(f"Excel doğrulandı: {output.resolve()}")
