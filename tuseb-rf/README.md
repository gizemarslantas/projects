# TÜSEB — PC kullanmayan Random Forest

Mevcut analiz kodunun sürüm kontrolü için başlangıç kopyasıdır. Üç kaynak betik değiştirilmeden alınmıştır. Genotip, birey kimliği, fenotip tablosu, kovaryat veya analiz sonucu içermez.

## Analiz

- Otozomlar (1–22) sırayla işlenir.
- PHENO=1 kontrol; PHENO=2 olgu. Modelde 0/1 olarak kodlanır.
- PC okunmaz, modele eklenmez ve genotiplerden PC etkisi çıkarılmaz. Kodda kalan PC değişkenleri sıfır sütunlu matrislerdir.
- Blok boyutu 10.000 SNP; blok başına 50 ağaç ve en fazla 50 aday SNP.
- Kromozom son modeli 1.000 ağaç; genom çapında aday modeli 2.000 ağaç.
- Sonuç tüm SNP'leri aynı anda kullanan tek bir RF değil, bloklardan seçilen adaylarla kurulan modellerdir.

## Girdiler (NAS üzerinde kalır)

`--base` klasörü altında:

```text
SLD_EA_c2.pheno.tsv
per_chr/chr1/SLD_EA_chr1_merged.pgen
per_chr/chr1/SLD_EA_chr1_merged.pvar
per_chr/chr1/SLD_EA_chr1_merged.psam
... chr22 için aynı yapı ...
```

Fenotip tablosunda IID ve PHENO sütunları gerekir. IID benzersiz olmalı; kromozomların birey sırası aynı olmalıdır.

## Kurulum ve çalıştırma

```bash
python3 -m pip install -r requirements.txt
bash run_rf_nopc.sh /path/to/analysis_base
```

Başlatıcıya Kasım içindeki analiz klasörünü açıkça verin; mevcut varsayılan yol eski analiz klasörüdür. Başlatıcı yeni, zaman damgalı sonuç klasörü oluşturur.

## Kaydedilen sürümün sınırları

- Başlatıcı ve Python varsayılanları 4 iş parçacığı kullanır; bu kopya 2 GB RAM için doğrulanmamıştır. Genom çapında aday matrisi bellekte tutulur.
- Bağımlılık sürümleri bilinmediği için requirements.txt sabitlenmemiştir. NAS ortamındaki gerçek sürümler ayrıca kaydedilmelidir.
- OOB metrikleri aday seçiminden bağımsız bir test değerlendirmesi değildir; aday seçimi aynı bireylerden yapılır.
- Yerelde yalnızca Python sözdizimi ve Bash sözdizimi kontrol edilmiştir; gerçek veriyle model çalıştırılmamıştır.
- GWAS betiklerinin tam kaynakları bu pakette yoktur; ayrıca eklenmelidir.

`SOURCE_SHA256SUMS`, değiştirilmeden alınan üç betiğin SHA-256 özetlerini içerir. Bu değerler NAS kopyalarıyla karşılaştırılarak aynı sürüm oldukları doğrulanabilir.
