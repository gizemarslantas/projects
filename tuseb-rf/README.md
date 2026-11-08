# TÜSEB — NAS'ta kullanılan güncel RF kodları

8 Kasım 2026 tarihli PC'siz RF, eşik üstü Excel ve hg19 gen anotasyonu akışı.

Proje taslağı için: [RF yöntem metni](METHODS.md) ve [parametrelerin açıklamaları](PARAMETERS.md). Bu belgeler güncel kodu, seçim gerekçelerini ve yorum sınırlarını açıklar.

| Sıra | Dosya | İşlem |
|---|---|---|
| 1 | `run_rf_10k_nopc_1thread.sh` | 22 otozomu sırayla işler ve genom çapındaki son modeli kurar. |
| 1a | `rf_rank_chr_nopc.py` | Blok taraması ve kromozom RF modeli. |
| 1b | `rf_combine_allchr_nopc.py` | Adayları birleştirip son RF sıralamasını oluşturur. |
| 2 | `export_rf_mean_excel.py` | Ortalama importance'ın kesin olarak üzerindeki SNP'leri Excel'e çıkarır. |
| 3 | `annotate_rf_hg19.py` | Gen, örtüşme/yakınlık, mesafe, biyotip ve Ensembl kimliği ekler. |

Eski `run_rf_nopc.sh` komutu da aynı güncel tek iş parçacıklı başlatıcıya yönlenir. Önceki dört iş parçacıklı başlatıcı Git geçmişinde kalır.

## Analiz ayarları

- PC okunmaz, modele eklenmez ve genotiplerden PC etkisi çıkarılmaz. Kalan PC değişkenleri sıfır sütunlu matrislerdir.
- PHENO=1 kontrol, PHENO=2 olgu; modelde sırasıyla 0 ve 1.
- Blok: 10.000 SNP, 50 ağaç, en fazla 50 aday.
- Kromozom modeli: 1.000 ağaç. Genom modeli: 2.000 ağaç. Minimum yaprak: 5.
- Başlatıcı `--threads 1` kullanır; sayısal kütüphanelerin iş parçacıkları da 1 ile sınırlandırılır.
- Kromozom seed'i 40417, genom modeli seed'i 4042201.

Python RF dosyaları değiştirilmemiştir; doğrudan çağrıldıklarında eski `--threads 4` varsayılanları vardır. Güncel başlatıcıyı kullanın veya doğrudan çağrıda `--threads 1` verin.

## NAS'ta çalıştırma

Kod klasöründe:

```bash
python3 -m pip install --user -r requirements.txt
ANALYSIS_BASE="/volume1/homes/tuseb/new_data_to_be_used/Kasım/yardimci_dosya_adaylari/GWAS_404_vs_SLD_allchr"
```

Bu klasörde `SLD_EA_c2.pheno.tsv` ve 1–22 kromozomları için `per_chr/chrN/SLD_EA_chrN_merged.pgen/.pvar/.psam` bulunmalıdır. Fenotip tablosunda IID ve PHENO gerekir. IID benzersiz, kromozomların birey sırası aynı olmalıdır.

### 1. Yalnızca yeni RF analizi başlatılacaksa

```bash
RF_LOG="rf_launcher_$(date +%Y%m%d_%H%M%S)_$$.log"
nohup bash run_rf_10k_nopc_1thread.sh "$ANALYSIS_BASE" > "$RF_LOG" 2>&1 &
echo "PID: $! | Log: $RF_LOG"
```

Girdi klasöründe `RF_noPC_10k_<tarih>_<saat>_<PID>` oluşturulur. Kullanılan iki Python dosyası da sonuç klasörüne kopyalanır.

### 2. Tamamlanmış sonuçtan Excel

Aşağıdaki yol tamamlanmış mevcut çalışmaya aittir. Yeni analizde `RUN_DIR` değerini değiştirin. Excel oluşturmak için RF'yi yeniden başlatmayın.

```bash
RUN_DIR="$ANALYSIS_BASE/RF_noPC_10k_20261108_181230_4597"
python3 export_rf_mean_excel.py --combined "$RUN_DIR/combined"
```

Exporter bu çalışmadaki 30.550 adayı ve importance toplamının 1 olduğunu doğrular. Filtre: `GENOME_RF_IMPORTANCE > ortalama` (yaklaşık 0,000032733224). Farklı aday sayılı yeni analizde sayı kontrolü bilinçli olarak güncellenmelidir.

### 3. Gen anotasyonu

`--input` olarak ikinci adımın yazdırdığı Excel yolunu kullanın. Örnek, son kullanılan dosyadır:

```bash
python3 annotate_rf_hg19.py \
  --base "$ANALYSIS_BASE" \
  --input "$RUN_DIR/combined/RF_ortalama_ustu_SNPler_20261108_221824_581056.xlsx" \
  --reference "/volume1/homes/tuseb/new_data_to_be_used/GWAS_404_vs_SLD_allchr/reference/Homo_sapiens.GRCh37.87.gtf.gz"
```

Girdi Excel'inde `SNPler` sayfası gerekir. Yeni, zaman damgalı Excel girdi dosyasının yanına kaydedilir. RF yeniden çalıştırılmaz. `--reference` verilmezse bilinen NAS yolları aranır; dosya yoksa Ensembl'den indirilir.

- `overlap`: SNP gen gövdesinin içinde, intronlar dahil; mesafe 0.
- `nearest`: gen dışında; en yakın gen sınırına baz çifti uzaklığı.
- Tüm biyotipler kullanılır; çoklu örtüşmeler/eşit uzaklıklar birlikte yazılır.
- Gen adları Ensembl GRCh37 release 87 (hg19) sürümündendir. Yakınlık işlevsel etki kanıtı değildir.

## Sürüm ve doğrulama

Sohbetteki son kodlar dosyalara ayrılmış; giriş yolları komut satırı seçenekleriyle belirtilebilir hale getirilmiştir. RF algoritması ve model ayarları değiştirilmemiştir. `SOURCE_SHA256SUMS` altı güncel analiz/çalıştırma dosyasının özetlerini içerir.

Python 3.8 sözdizimi ve Bash sözdizimi kontrolleri yapılır; anotasyonda `str.removeprefix()` kullanılmaz. Excel eşiği kontrollü sentetik veride, gen anotasyonu sağlanan RF listesi ve GRCh37.87 referansında yerelde doğrulanır. NAS'ta RF yeniden çalıştırılmamıştır; 2 GB tepe bellek garantisi yoktur. Bağımlılık sürümleri NAS ortamından henüz sabitlenmemiştir.

RF importance p-değeri değildir. Son model tüm SNP'leri aynı anda değil, blok finalistlerini kullanır. Aday seçimi tüm bireyleri kullandığından kromozom OOB metrikleri bağımsız test başarısı sayılmaz.

Genotip, birey kimliği, fenotip tabloları, referanslar ve sonuçlar Git'e eklenmez. GWAS kodları ve yöntem açıklaması ayrı klasördedir: [tuseb-gwas](../tuseb-gwas/README.md). GWAS PC1–PC5 kullanır; güncel RF PC kullanmaz.
