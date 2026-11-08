# TÜSEB — 404 ve SLD tüm otozom PC5 GWAS kodları

Bu klasör 6.003.495 varyantlık 404–SLD analizinin kayıtlı GWAS kodlarını içerir. Ana üç Bash betiği, önceki sohbetten tam metin olarak geri alınmış ve değiştirilmeden saklanmıştır. NAS'taki dosyaların sonradan değiştirilip değiştirilmediği hash ile henüz karşılaştırılmamıştır.

## Dosyalar ve sıra

| Sıra | Dosya | İşlem |
|---|---|---|
| 1 | `01_prepare_allchr.sh` | Otozom bazında ortak varyantların bulunması, alel hizalama ve projeye özgü SLD doz yönü düzeltmesi. |
| 2 | `02_merge_and_gwas_allchr.sh` | SLD ve EA dozajlarının VCF üzerinden birleştirilmesi ve PC1–PC5 ile GWAS. |
| 3 | `03_combine_allchr_results.sh` | Geçerli ADD testlerini toplar; tüm genomda Bonferroni, BH-FDR, sıra ve özet hesaplar. |
| Yardımcı | `05_gwas_chr_ozet.py` | Mevcut tüm-genom sonuçlarını kromozom bazında özetler; GWAS veya düzeltme hesaplarını tekrarlamaz. |
| Tarihsel | `historical/04_annotate_gwas_rf_excel.py` | Eski GWAS Top 1000 ve ortak 22 SNP Excel'leri. Güncel ortak 25 veya RF eşik üstü listesinin kodu değildir. |

Proje taslağı için [METHODS.md](METHODS.md) kodla doğrulanabilen yöntemi, sayıları ve eksik kaynakları açıklar. Güncel PC'siz RF ayrı klasördedir: [tuseb-rf](../tuseb-rf/README.md).

**GWAS PC1–PC5 kullanır; güncel RF PC kullanmaz.** Ek MAF %5/%10 denemeleri bu ana GWAS yöntemine eklenmemiştir.

## Girdiler ve bağımlılıklar

Kayıtlı betiklerde ana yol `/volume1/homes/tuseb/new_data_to_be_used/GWAS_404_vs_SLD_allchr` olarak sabittir. Klasörü değiştirmek, betikleri Kasım'a kopyalamakla otomatik gerçekleşmez; yollar çalıştırmadan önce kontrol edilmelidir.

Hazırlık ayrıca şu kaynakları kullanır:

```text
/volume1/homes/tuseb/new_data_to_be_used/QC/EA_c2_qc_dedup.pgen/.pvar/.psam
/volume1/homes/tuseb/new_data_to_be_used/imputation_for_SLD_data/output/plink_final/SLD_chrN_qc.pgen/.pvar/.psam
GWAS_404_vs_SLD_allchr/EA_c2_1580.keep
GWAS_404_vs_SLD_allchr/SLD_EA_c2.pheno.tsv
GWAS_404_vs_SLD_allchr/SLD_EA_c2_PC5.covar.tsv
```

İlk betik tek bir tüm-genom SLD final dosyası yerine **kromozom bazındaki `SLD_chrN_qc`** dosyalarını bekler. Fenotip 1=kontrol, 2=olgu; toplam 2.351 birey (771 SLD kontrol, 1.580 EA olgu). Kovaryat başlıkları PC1–PC5'tir.

PLINK 2, bcftools, Bash, GNU sort ve Python >=3.8 gerekir. Python paketleri `requirements.txt` içindedir. bcftools yolu ikinci betikte sabittir. Orijinal GWAS betikleri 4 PLINK iş parçacığı, bcftools merge için 2 iş parçacığı kullanır; 2 GB bellek garantisi vermez.

## Çalıştırma

Bu komutlar analizleri yeniden başlatabilir. Proje taslağı yazmak veya mevcut sonuçları okumak için çalıştırmak gerekmez. Kaynak ve çıktı yollarını doğruladıktan sonra:

```bash
bash 01_prepare_allchr.sh 1 22
bash 02_merge_and_gwas_allchr.sh
bash 03_combine_allchr_results.sh
python3 05_gwas_chr_ozet.py
```

Hazırlık, daha önce projede tespit edilen SLD alel/doz yönü sorununa özel olarak **tüm seçilen SLD varyantlarında** yön düzeltmesi uygular. Bu adım genel amaçlı bir harmonizasyon reçetesi değildir; önceden düzeltilmiş yeni verilere tekrar uygulanmamalıdır. Betikler arşivlenen yöntemi belgelemek için özgün bırakılmıştır.

İkinci betik `.done`, sonuç ve kromozom özeti varsa ilgili kromozomu atlar. Bu dosyalar tek başına girdilerin aynı sürümde olduğunu garanti etmez. Üçüncü betik bu çalışma için 6.003.487 geçerli test bekler; farklı veri seti için kontrolü bilinçli olarak uyarlamak gerekir.

## Çıktılar

```text
harmonization/allchr_prepare_summary.tsv
gwas_results/GWAS_chrN_corrected_PC5.PHENO.glm.logistic.hybrid
gwas_results/allchr_gwas_summary.tsv
combined_results/GWAS_allchr_corrected_PC5_all_results.tsv.gz
combined_results/GWAS_allchr_corrected_PC5_top1000_basic.tsv
combined_results/GWAS_allchr_corrected_PC5_summary.tsv
combined_results/GWAS_chr_summary_for_excel.tsv
```

Raporlanan çalışmada 6.003.495 ADD testi, 6.003.487 geçerli p-değeri ve 8 geçersiz test vardır. Genom genelindeki düzeltmeler yalnızca geçerli testleri kullanır. FDR<1 veya Bonferroni<1, 0,05 anlamlılık sınırı değildir.

## Kaynak ve doğrulama

- İlk üç betik: “TÜSEB 404 vs SLD chr17 GWASı” sohbetindeki tam kodlar (`01a11cfb-52f4-7341-89f8-bda2245b7cd0`).
- Tarihsel anotasyon: 8 Ekim 2026 yerel çıktı dosyası.
- Kromozom özeti: 18 Ekim 2026 yerel çıktı dosyası.
- `SOURCE_SHA256SUMS`: buradaki beş kod dosyasının özetleri; NAS sürümleriyle karşılaştırılabilir.

Bash ve Python sözdizimi kontrol edilmiştir. Tüm-genom toplama ve Bonferroni/BH hesapları kontrollü 22 kromozomluk sentetik sonuçlarla sınanır; bu test gerçek genotiplerle GWAS'ın yeniden çalıştırılması değildir. NAS araç sürümleri henüz sabitlenmemiştir.

Veri, birey listeleri, fenotip ve kovaryat değerleri, referans genom dosyaları ve sonuçlar Git'e yüklenmez. Girdi QC, fenotip oluşturma ve PC5 hesaplama betiklerinin tam kaynakları henüz bu pakette değildir; bunlara ait parametreler tahmin edilmemiştir.
