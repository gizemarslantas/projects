# Proje taslağı için Random Forest yöntemi

Bu belge `run_rf_10k_nopc_1thread.sh`, `rf_rank_chr_nopc.py` ve `rf_combine_allchr_nopc.py` dosyalarının mevcut davranışını açıklar. Sayılar 8 Kasım 2026 tarihli tamamlanan çalışmaya aittir. Parametrelerin tek tek açıklaması: [PARAMETERS.md](PARAMETERS.md).

## Taslakta kullanılabilecek yöntem metni

Random Forest (RF), 404/EA olguları ile SLD kontrollerini ayırt etmede kullanılan SNP'lerin önceliklendirilmesi amacıyla uygulanmıştır. Analiz 1.580 olgu ve 771 kontrol olmak üzere 2.351 bireyi ve otozomlar üzerinde toplam 6.003.495 SNP'yi kapsamıştır. Birleşik ve harmonize PGEN dosyalarındaki alel dozajları `pgenlib` ile okunmuştur. Fenotip tablosundaki 1=kontrol ve 2=olgu etiketleri model içinde sırasıyla 0 ve 1 olarak kodlanmıştır. Bireyler IID üzerinden eşleştirilmiş; tekrarlanan kimlikler, eksik fenotipler ve dosyalar arası örnek/varyant sayısı uyumsuzlukları kontrol edilmiştir. Genom çapındaki birleştirmede kromozomların birey sıraları ayrıca chr1 ile karşılaştırılmıştır.

Bellek kullanımını yönetmek için 22 otozom sırayla işlenmiş, her kromozom dosyadaki varyant sırasına göre ardışık 10.000 SNP'lik bloklara ayrılmıştır. Son blok daha küçük olabilmektedir. Her blok için 50 ağaçlı bir RF sınıflandırıcısı eğitilmiş ve Gini safsızlık azalmasına dayalı önem değeri en yüksek en fazla 50 SNP sonraki aşamaya alınmıştır. Böylece GWAS p-değerine dayalı ön seçim yapılmadan tüm SNP'ler ilk taramaya dahil edilmiştir. Bloklar biyolojik lokus veya LD bloğu olarak tanımlanmamıştır.

Her kromozomda bloklardan gelen adayların birleşimiyle 1.000 ağaçlı ikinci bir RF modeli kurulmuş ve kromozom içi önem sıralaması elde edilmiştir. Ardından 22 kromozomun blok finalistlerinin tamamı birleştirilerek 30.550 aday SNP üzerinde 2.000 ağaçlı genom çapında RF modeli eğitilmiştir. Genom modeline yalnızca kromozom Top 1000 listeleri alınmamıştır; her kromozomun tam aday sıralama dosyası kullanılmıştır. Son SNP sıralaması genom modelinin normalize edilmiş Gini önem değerlerine göre oluşturulmuştur.

Tüm RF modellerinde Gini ölçütü, her bölünmede toplam özellik sayısının karekökü kadar rastgele özellik değerlendirilmesi, en az 5 örnekli yapraklar, üst sınır belirtilmemiş ağaç derinliği, yerine koyarak örnekleme ve ağaç başına örnek sayısının %80'i kadar bootstrap çekimi kullanılmıştır. Sınıf dengesizliği `class_weight="balanced_subsample"` ile her ağacın bootstrap örnekleminde sınıf frekansına göre ağırlıklandırılmıştır. Güncel başlatıcı modelleri tek iş parçacığıyla çalıştırmaktadır.

RF modellerine PC kovaryatları eklenmemiş ve RF ön işlemesinde genotiplerden PC etkisi çıkarılmamıştır. Kodda bulunan PC değişkenleri sıfır sütunlu matrislerdir. Bu RF, PC1–PC5 içeren GWAS modelinden bu yönüyle ayrılır.

Kodun eksik olarak tanımladığı negatif dozajlar, ilgili SNP'nin o okumadaki gözlenen bireyler üzerindeki ortalama dozajıyla doldurulmaktadır. Bir SNP için bütün dozajlar eksikse kod ortalamayı 0 kabul eder. Dozaj ve model girdi matrisleri `float32` biçimindedir. Tamamlanan genom modelinde raporlanan `missing_dosages_imputed=0`, son modele alınan aday dozajlarında doldurma gerekmediğini gösterir; ilk taramadaki bütün SNP'lerde eksik veri olmadığını tek başına göstermez.

Kromozom düzeyindeki modellerde out-of-bag (OOB) tahminleri üzerinden ROC-AUC ve balanced accuracy hesaplanmıştır; sınıflama eşiği 0,5'tir. Blok modellerinde OOB kapalıdır. Genom çapındaki son modelde OOB etkinleştirilmemiş, son model için bağımsız test veya çapraz doğrulama performansı hesaplanmamıştır.

## RF importance eşiği ve raporlama

Kullanılan importance, eğitim verisindeki ağırlıklı Gini azalmasının ağaçlar boyunca birleştirilip normalize edilmesidir. Son modeldeki 30.550 SNP için değerlerin toplamı 1'dir. Genom importance'ı, blok veya kromozom importance'ıyla sayısal olarak aynı ölçekte bir karşılaştırma gibi yorumlanmamalıdır; farklı modellerde normalize edilmişlerdir.

Excel raporu için son modelin adayları arasındaki ortalama importance eşik olarak alınmıştır:

```text
eşik = importance toplamı / son modeldeki SNP sayısı
     = 1 / 30.550 ≈ 0,000032733224

seçim: GENOME_RF_IMPORTANCE > eşik
```

Bu filtre 768 SNP'yi korumuştur. Ortalama, 6.003.495 taranan SNP üzerinden değil, son RF'ye alınan 30.550 aday üzerinden hesaplanmıştır. Bu eşik aday önceliklendirme tercihidir; p-değeri, FDR, Bonferroni sınırı veya evrensel RF anlamlılık eşiği değildir. Koddaki Top 1000 çıktısı ayrı bir sabit sıra raporudur; ortalama üstü 768 SNP listesiyle karıştırılmamalıdır.

Gen anotasyonu Ensembl GRCh37 release 87 (hg19) gen gövdelerine göre yapılmıştır. Gen sınırları arasında kalan SNP'ler `overlap`; diğerleri en yakın gen sınırına uzaklıkla `nearest` olarak raporlanır. Tüm gen biyotipleri, birden çok örtüşme ve eşit uzaklıklar korunur. Yakınlık işlevsel etki veya nedensellik kanıtı değildir.

## Yorum sınırları ve tamamlanacak doğrulamalar

- Bu üç aşamalı aday seçimi, tüm 6 milyon SNP'nin aynı anda yer aldığı tek bir RF ile eşdeğer değildir. İlk blok aşamasında elenen SNP'ler geri gelmez; farklı bloklardaki etkileşimler ilk taramada kaçabilir.
- Blok başına 50 ağaç/50 aday, kromozomda 1.000 ağaç, genomda 2.000 ağaç ve minimum yaprak 5 mevcut çalıştırma ayarlarıdır. Bu betiklerde optimizasyon, birden fazla seed konsensüsü veya ayar kararlılığı analizi yoktur; bu çalışma için kanıtlanmış optimum değerler olarak sunulmamalıdır.
- SNP seçimi bütün bireylerin etiketleriyle yapılmıştır. Bu nedenle kromozom OOB metrikleri tüm seçim sürecinin bağımsız testi değildir. Dış örneklem veya seçim adımlarını eğitim bölümü içinde tekrarlayan doğrulama yapılmadan genellenebilir performans iddiası kurulamaz.
- Olgu ve kontrol etiketleri farklı veri kohortlarıyla çakışır. PC'siz RF popülasyon yapısı veya ölçüm/imputasyon/batch farklarını ayırıcı sinyal olarak kullanabilir. Yüksek importance biyolojik ilişkiyi tek başına doğrulamaz.
- Sabit seed tekrarlanabilirliği destekler; aynı sonuç için birey ve varyant sırası, girdi dosyaları, bağımlılık sürümleri ve model ayarları da kaydedilmelidir. NAS ortamının tam bağımlılık sürümleri henüz sabitlenmemiştir.
- 10.000 SNP'lik bir blok için yalnızca `2351 × 10000 × 4` baytlık float32 matris yaklaşık 89,7 MiB'dir; son aday matrisi yaklaşık 274,0 MiB'dir. Kopyalar, eksiklik maskeleri, ağaçlar ve metadata ek bellek kullanır. Bu değerler 2 GB tepe bellek garantisi değildir.
- Raporlanan 3,05 dakika yalnızca genom çapındaki son model betiğinin süresidir; tüm blok/kromozom analizinin toplam süresi değildir.

## Kod ve teknik kaynaklar

- [Güncel başlatıcı](run_rf_10k_nopc_1thread.sh)
- [Blok/kromozom RF kodu](rf_rank_chr_nopc.py)
- [Genom RF kodu](rf_combine_allchr_nopc.py)
- [Ortalama üstü Excel](export_rf_mean_excel.py)
- [hg19 anotasyonu](annotate_rf_hg19.py)
- [Scikit-learn RandomForestClassifier](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)
- [Scikit-learn veri sızıntısı ve özellik seçimi](https://scikit-learn.org/stable/common_pitfalls.html#data-leakage)

Sayısal sonuçlar kullanıcının NAS çalışma özeti ve sağladığı RF Excel'inden alınmıştır; bu belgeyi hazırlamak için model yeniden eğitilmemiştir.
