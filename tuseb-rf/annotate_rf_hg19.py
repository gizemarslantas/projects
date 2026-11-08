#!/usr/bin/env python3
from pathlib import Path
from datetime import datetime
from collections import defaultdict, Counter
import argparse
import gzip
import re
import shutil
import urllib.request

import numpy as np
import pandas as pd

try:
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.worksheet.table import Table, TableStyleInfo
    from openpyxl.utils import get_column_letter
except ImportError:
    raise SystemExit(
        "openpyxl eksik. Önce çalıştır:\n"
        "python3 -m pip install --user openpyxl"
    )


def normalize_chr(value):
    chrom = str(value).strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    return chrom


parser = argparse.ArgumentParser(description="RF eşik üstü Excel listesine GRCh37/hg19 gen anotasyonu ekler.")
parser.add_argument("--base", type=Path, default=Path.cwd(), help="NAS analiz klasörü")
parser.add_argument("--input", type=Path, help="SNPler sayfasını içeren RF Excel dosyası")
parser.add_argument("--reference", type=Path, help="Ensembl GRCh37 release 87 GTF.gz dosyası")
args = parser.parse_args()
base = args.base.resolve()
original = Path("/volume1/homes/tuseb/new_data_to_be_used")

folder = (
    base
    / "RF_noPC_10k_20261108_181230_4597"
    / "combined"
)

source = (
    folder
    / "RF_ortalama_ustu_SNPler_20261108_221824_581056.xlsx"
)

if args.input is not None:
    source = args.input.resolve()
    folder = source.parent

if not source.is_file():
    raise SystemExit(
        "RF Excel dosyası bulunamadı:\n{}".format(source)
    )

# Kaynak dosyayı oku. Özgün sıra ve importance değerlerini koru.
df = pd.read_excel(
    source,
    sheet_name="SNPler",
    dtype={
        "ID": str,
        "CHR": str,
        "REF": str,
        "ALT": str
    }
)

required = {
    "ID", "CHR", "POS", "REF", "ALT",
    "GENOME_RF_RANK", "GENOME_RF_IMPORTANCE"
}

missing = required - set(df.columns)
if missing:
    raise SystemExit(
        "Eksik sütunlar: {}".format(sorted(missing))
    )

if df.empty:
    raise SystemExit("RF listesi boş.")

annotation_columns = [
    "GENE",
    "GENE_RELATION",
    "DISTANCE_BP",
    "GENE_BIOTYPE",
    "ENSEMBL_GENE_ID"
]

if set(annotation_columns).intersection(df.columns):
    raise SystemExit(
        "Kaynak dosyada zaten gen anotasyon sütunları var. "
        "Orijinal RF listesini kullan."
    )

positions = pd.to_numeric(df["POS"], errors="raise")

if (
    not np.isfinite(positions.to_numpy(dtype=float)).all()
    or (positions <= 0).any()
    or (positions % 1 != 0).any()
):
    raise SystemExit("Geçersiz SNP pozisyonu var.")

normalized_chr = df["CHR"].map(normalize_chr)
allowed_chr = set(str(i) for i in range(1, 23))

if not normalized_chr.isin(allowed_chr).all():
    raise SystemExit(
        "Listede 1–22 dışında veya geçersiz kromozom var."
    )

keys = pd.DataFrame({
    "CHR": normalized_chr,
    "POS": positions.astype("int64"),
    "REF": df["REF"],
    "ALT": df["ALT"]
})

if keys.duplicated().any():
    raise SystemExit("RF listesinde tekrarlanan varyant var.")

# Gen referansını bul.
gtf_name = "Homo_sapiens.GRCh37.87.gtf.gz"

candidates = [
    base / "reference" / gtf_name,
    original / "GWAS_404_vs_SLD_allchr" / "reference" / gtf_name,
    original / "reference" / gtf_name
]

if args.reference is not None:
    if not args.reference.is_file():
        raise SystemExit("Belirtilen gen referansı bulunamadı: {}".format(args.reference))
    candidates.insert(0, args.reference.resolve())

gtf = next(
    (p for p in candidates if p.is_file() and p.stat().st_size > 0),
    None
)

url = (
    "https://ftp.ensembl.org/pub/grch37/release-87/"
    "gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
)

if gtf is None:
    gtf = candidates[0]
    gtf.parent.mkdir(parents=True, exist_ok=True)

    temporary = gtf.with_name(
        gtf.name + ".{}.part".format(
            datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        )
    )

    print(
        "GRCh37/hg19 gen referansı indiriliyor...",
        flush=True
    )

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "TUSEB-gene-annotation"}
    )

    with urllib.request.urlopen(request, timeout=120) as response:
        with temporary.open("wb") as writer:
            shutil.copyfileobj(
                response, writer, length=1024 * 1024
            )

    with gzip.open(temporary, "rt") as reader:
        first = reader.readline()
        if not first.startswith("#"):
            raise SystemExit(
                "İndirilen referans beklenen GTF biçiminde değil."
            )

    temporary.replace(gtf)

print("Gen referansı: {}".format(gtf), flush=True)

# Gen gövdelerini yükle. İntronlar ve tüm gen biyotipleri kapsamda.
genes = defaultdict(list)
assembly_verified = False

with gzip.open(gtf, "rt") as reader:
    for line in reader:
        if line.startswith("#"):
            if "GRCh37" in line:
                assembly_verified = True
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9 or fields[2] != "gene":
            continue

        chrom = normalize_chr(fields[0])
        if chrom not in allowed_chr:
            continue

        attributes = dict(
            re.findall(r'(\S+)\s+"([^"]*)"', fields[8])
        )

        gene_id = attributes.get("gene_id")
        if not gene_id:
            continue

        genes[chrom].append({
            "start": int(fields[3]),
            "end": int(fields[4]),
            "name": attributes.get("gene_name", gene_id),
            "biotype": attributes.get("gene_biotype", ""),
            "id": gene_id
        })

if not assembly_verified:
    raise SystemExit(
        "GTF başlığında GRCh37 doğrulanamadı."
    )

if not genes:
    raise SystemExit(
        "Referanstan gen aralıkları okunamadı."
    )

print("Gen anotasyonu başlıyor...", flush=True)

annotations = []

for index, row in enumerate(df.itertuples(index=False), 1):
    chrom = normalize_chr(row.CHR)
    pos = int(row.POS)
    chrom_genes = genes.get(chrom)

    if not chrom_genes:
        raise SystemExit(
            "Referansta kromozom bulunamadı: {}".format(chrom)
        )

    # Gen içindeyse mesafe 0.
    # Gen dışındaysa en yakın gen sınırına uzaklık.
    distances = [
        (
            max(g["start"] - pos, pos - g["end"], 0),
            g
        )
        for g in chrom_genes
    ]

    distance = min(d for d, g in distances)

    selected = sorted(
        [g for d, g in distances if d == distance],
        key=lambda g: (g["name"], g["id"])
    )

    annotations.append({
        "GENE": "; ".join(g["name"] for g in selected),
        "GENE_RELATION": (
            "overlap" if distance == 0 else "nearest"
        ),
        "DISTANCE_BP": distance,
        "GENE_BIOTYPE": "; ".join(
            g["biotype"] for g in selected
        ),
        "ENSEMBL_GENE_ID": "; ".join(
            g["id"] for g in selected
        )
    })

    if index % 100 == 0 or index == len(df):
        print(
            "{}/{} SNP anotasyonlandı.".format(index, len(df)),
            flush=True
        )

result = pd.concat(
    [
        df.reset_index(drop=True),
        pd.DataFrame(annotations)
    ],
    axis=1
)

assert len(result) == len(df)
assert result["GENOME_RF_RANK"].tolist() == (
    df["GENOME_RF_RANK"].tolist()
)
assert np.array_equal(
    result["GENOME_RF_IMPORTANCE"].to_numpy(),
    df["GENOME_RF_IMPORTANCE"].to_numpy()
)

stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
output = folder / (
    "RF_ortalama_ustu_hg19_Gen_{}.xlsx".format(stamp)
)

counts = Counter(result["GENE_RELATION"])

with pd.ExcelWriter(output, engine="openpyxl") as writer:
    result.to_excel(
        writer,
        sheet_name="RF SNP Gen",
        index=False,
        startrow=3
    )

    ws = writer.book["RF SNP Gen"]
    last_col = get_column_letter(len(result.columns))
    last_row = len(result) + 4

    ws.sheet_view.showGridLines = False
    ws.freeze_panes = "A5"

    ws.merge_cells("A1:{}1".format(last_col))
    ws["A1"] = (
        "Random Forest ortalama importance üzerindeki SNP'ler"
    )
    ws["A1"].font = Font(
        name="Arial",
        size=14,
        bold=True,
        color="1F1F1F"
    )
    ws.row_dimensions[1].height = 25

    ws.merge_cells("A2:{}2".format(last_col))
    ws["A2"] = (
        "Ensembl GRCh37 release 87 (hg19). "
        "Gen dışındaki SNP'ler en yakın gene atanmıştır."
    )
    ws["A2"].font = Font(
        name="Arial", size=10, color="595959"
    )

    table = Table(
        displayName="RFAnnotated",
        ref="A4:{}{}".format(last_col, last_row)
    )
    table.tableStyleInfo = TableStyleInfo(
        name="TableStyleMedium2",
        showRowStripes=True,
        showColumnStripes=False
    )
    ws.add_table(table)

    for cell in ws[4]:
        cell.fill = PatternFill(
            "solid", fgColor="1F4E78"
        )
        cell.font = Font(
            name="Arial",
            size=10,
            bold=True,
            color="FFFFFF"
        )
        cell.alignment = Alignment(
            horizontal="center",
            vertical="center",
            wrap_text=True
        )

    ws.row_dimensions[4].height = 34

    for column_index, name in enumerate(result.columns, 1):
        letter = get_column_letter(column_index)

        width = {
            "ID": 27,
            "GENE": 32,
            "GENE_RELATION": 18,
            "GENE_BIOTYPE": 26,
            "ENSEMBL_GENE_ID": 30,
            "REF": 8,
            "ALT": 8,
            "CHR": 8
        }.get(name, max(16, len(name) + 2))

        ws.column_dimensions[letter].width = width

        for cells in ws.iter_rows(
            min_row=5,
            max_row=last_row,
            min_col=column_index,
            max_col=column_index
        ):
            cell = cells[0]
            cell.font = Font(name="Arial", size=10)
            cell.alignment = Alignment(
                vertical="center",
                wrap_text=name in {
                    "GENE",
                    "GENE_BIOTYPE",
                    "ENSEMBL_GENE_ID"
                }
            )

            if "IMPORTANCE" in name:
                cell.number_format = "0.0000000000"
            elif name in {
                "GENOME_RF_RANK",
                "VARIANT_INDEX",
                "POS",
                "DISTANCE_BP"
            }:
                cell.number_format = "0"
            elif name in {
                "ID", "CHR", "REF", "ALT"
            }:
                cell.number_format = "@"

    for excel_row, names in enumerate(result["GENE"], 5):
        gene_count = names.count(";") + 1
        ws.row_dimensions[excel_row].height = max(
            25, 18 * gene_count
        )

    notes = [
        (
            "overlap: SNP gen gövdesi içinde; intronlar dahil. "
            "DISTANCE_BP=0."
        ),
        (
            "nearest: SNP gen dışında; DISTANCE_BP en yakın "
            "gen sınırına baz çifti uzaklığıdır."
        ),
        (
            "Birden çok örtüşen veya eşit uzaklıktaki gen, "
            "aynı sırayla noktalı virgülle ayrılmıştır."
        ),
        (
            "Gen yakınlığı, SNP'nin o geni etkilediğini "
            "kanıtlamaz. Gen adları Ensembl 87 sürümündendir."
        ),
        (
            "Toplam SNP: {}. Genle örtüşen: {}. "
            "En yakın gen atanan: {}."
        ).format(
            len(result),
            counts["overlap"],
            counts["nearest"]
        ),
        "Kaynak RF: {}".format(source.name),
        "Gen referansı: {}".format(gtf_name),
        "Referans adresi: {}".format(url)
    ]

    for row_number, text in enumerate(notes, last_row + 3):
        ws.merge_cells(
            "A{}:{}{}".format(
                row_number, last_col, row_number
            )
        )
        cell = ws.cell(row_number, 1, text)
        cell.font = Font(
            name="Arial", size=10, color="595959"
        )
        cell.alignment = Alignment(
            wrap_text=True, vertical="center"
        )
        ws.row_dimensions[row_number].height = 30

# Yalnızca veri satırlarını yeniden okuyarak Excel'i doğrula.
saved = pd.read_excel(
    output,
    sheet_name="RF SNP Gen",
    header=3,
    nrows=len(result),
    dtype={"ID": str}
)

assert len(saved) == len(result)
assert saved["ID"].tolist() == result["ID"].tolist()
assert saved["GENOME_RF_RANK"].tolist() == (
    result["GENOME_RF_RANK"].tolist()
)
assert saved["GENE"].tolist() == result["GENE"].tolist()
assert saved["GENE_RELATION"].tolist() == (
    result["GENE_RELATION"].tolist()
)
assert saved["DISTANCE_BP"].tolist() == (
    result["DISTANCE_BP"].tolist()
)
assert saved["ENSEMBL_GENE_ID"].tolist() == (
    result["ENSEMBL_GENE_ID"].tolist()
)

for column in result.columns:
    if "IMPORTANCE" in column:
        assert np.allclose(
            saved[column].to_numpy(dtype=float),
            result[column].to_numpy(dtype=float),
            rtol=1e-12,
            atol=1e-18
        )

print("\nToplam SNP: {}".format(len(result)))
print("Genle örtüşen: {}".format(counts["overlap"]))
print("En yakın gen atanan: {}".format(counts["nearest"]))
print("Excel doğrulandı:\n{}".format(output.resolve()))
