# RF parametreleri — değer, anlam ve seçim gerekçesi

Bu değerler güncel NAS başlatıcısının verdiği ayarlardır. Başlatıcıdaki değerler, Python dosyalarının varsayılanlarından önce gelir. Ayrıntılı çalışma yöntemi: [METHODS.md](METHODS.md).

## 1. Analiz akışını belirleyen ayarlar

| Ayar | Kullanılan değer | Ne yapar? | Gerekçe ve sınır |
|---|---:|---|---|
| Kromozom aralığı | 1–22 | Otozomlar sırayla işlenir. | X, Y ve mitokondri bu çalışmaya dahil değildir. |
| `--block-size` | 10.000 SNP | Varyant dosyasındaki ardışık SNP sütunlarını böler. Bireyler parçalanmaz. | Bellek ve çalışma süresini yönetir. LD/lokus tabanlı biyolojik bloklama değildir. |
| `--block-trees` | 50 ağaç | Her blokta bağımsız RF kurar. | Hızlı tarama tercihi; yeterlilik/kararlılık bu betikte optimize edilmez. |
| `--top-per-block` | en fazla 50 SNP | Her bloğun en yüksek önem değerli SNP'lerini seçer; küçük blokta `min(50, blok SNP sayısı)`. | Aday havuzunu küçültür. Sabit kontenjan olduğundan zayıf bloktan da aday seçilebilir; önemler sıfırsa da kontenjan uygulanır. |
| `--final-trees` | 1.000 ağaç | Bir kromozomdaki bütün blok adaylarıyla model kurar. | Kromozom önem sıralaması sağlar. Genom aşamasına sadece bu modelin Top 1000'i değil, tüm adayları gider. |
| Genom `--trees` | 2.000 ağaç | 22 kromozomun adaylarıyla son model kurar. | Tek ve ortak önem sıralaması sağlar; erken elenen SNP'leri geri getirmez. |
| PC sütunları | 0 | Genotip matrisine PC eklenmez. | Genotipler PC'ye göre residualize de edilmez. Popülasyon yapısı bu yolla kontrol edilmiş değildir. |
| Fenotip | 1→0 kontrol; 2→1 olgu | Hedef sınıfı oluşturur. | IID eşleştirmesiyle PGEN birey sırası korunur. |
| Veri tipi | `float32` | Dozaj ve X matrisleri 32 bit kayan noktalıdır. | Matris belleğini azaltır. Betikler tüm-genom final matrisini RAM'de tutar; otomatik disk matrisi kullanmaz. |
| Eksik dozaj | `g < 0` | SNP başına gözlenen ortalama dozajla doldurulur; tamamen eksik SNP'de 0. | Bu betikteki eksiklik tanımıdır; NaN/sonsuz değerler için ayrı bir genel temizlik uygulanmaz. |

## 2. RandomForestClassifier ayarları

Aşağıdaki ortak ayarlar blok, kromozom ve genom modellerinde kullanılır. Kavramsal tanımlar [scikit-learn belgeleri](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html) ile uyumludur.

| Parametre | Değer | Anlamı | Kullanım gerekçesi ve sınır |
|---|---|---|---|
| `n_estimators` | 50 / 1.000 / 2.000 | Kurulan karar ağacı sayısı. | Aşama ilerledikçe daha fazla ağaç kullanılır. Bu değerlerin optimum olduğu bu kodla gösterilmez. |
| `criterion` | `gini` | Bölünmenin sınıf karışıklığını ne kadar azalttığını ölçer. | İkili sınıflandırma için kullanılan ölçüt; importance de bu safsızlık azalmasına dayanır. |
| `max_features` | `sqrt` | Her düğümde yaklaşık toplam özellik sayısının karekökü kadar rastgele özellik adayına bakılır. | Ağaçları çeşitlendirir; 10.000 SNP'lik blokta yaklaşık 100 aday özellik. |
| `max_depth` | `None` | Açık ağaç derinliği sınırı yoktur. | Bölünmeler yaprak ve diğer durdurma koşullarıyla sonlanır; bellek veya aşırı uyumu tek başına sınırlamaz. |
| `min_samples_leaf` | 5 | Bölünme sonrası her yaprakta en az 5 eğitim örneği gerekir. | Çok küçük gruplara özel dalları azaltır; küçük alt gruplardaki sinyalleri de zayıflatabilir. Örneğin 12 örnek 8+4 bölünemez, 7+5 bölünebilir. |
| `class_weight` | `balanced_subsample` | Sınıf ağırlıkları her ağacın bootstrap örneklemindeki frekanslardan hesaplanır. | 771 kontrol/1.580 olgu dengesizliğinin etkisini azaltmayı hedefler. Yeni birey üretmez, SMOTE yapmaz veya sınıf sayılarını zorunlu eşitlemez. |
| `bootstrap` | `True` | Her ağaç için yerine koyarak birey örneklenir. | Ağaçlar farklı örneklemler görür; aynı birey bir ağaç için birden fazla kez çekilebilir. |
| `max_samples` | `0.80` | Ağaç başına yaklaşık N×0,8 bootstrap çekimi. | 2.351 birey için yaklaşık 1.881 çekim; benzersiz birey sayısı bu kadar olmak zorunda değildir. %80 eğitim/%20 sabit test ayrımı anlamına gelmez. |
| `n_jobs` | başlatıcıda 1 | Ağaçlar için paralel iş sayısı. | NAS belleğini yönetmek için tek iş parçacığı. Python dosyalarının doğrudan CLI varsayılanı 4'tür; başlatıcı bunu 1 ile değiştirir. |
| `random_state` | aşağıdaki seed düzeni | Rastgele bootstrap ve özellik seçimlerini belirler. | Tekrar üretilebilirliği destekler; bağımlılık ve veri sırası değişikliklerinde birebir sonuç garantisi değildir. |

## 3. Seed ve OOB düzeni

| Aşama | Gerçek seed hesabı | OOB |
|---|---|---|
| Blok | `40417 + block_no`; blok numarası her kromozomda yeniden başlar. | Kapalı (`oob_score=False`). |
| Kromozom modeli | `40417 + 100000 + chrom`; chr1 için 140418. | Açık (`oob_score=True`). |
| Genom modeli | 4042201 | Etkinleştirilmemiş; sklearn varsayılanı kapalı. |

Kromozom OOB olasılıklarından ROC-AUC hesaplanır; olasılık ≥0,5 için olgu tahmini yapılıp balanced accuracy hesaplanır. `oob_score_` özelliğinin kendi accuracy değeri ayrı bir çıktı olarak yazılmaz. Genom modelinin özetinde bağımsız AUC/accuracy bulunmaz. Önceki aday seçimi tüm bireyleri kullandığından OOB, tüm akışın bağımsız değerlendirmesi değildir.

## 4. Sıralama ve raporlama

| Alan / ayar | İşlem ve yorum |
|---|---|
| `BLOCK_RF_IMPORTANCE` | Yalnızca ilgili blok modelindeki normalize önem. Eşit önemlerde küçük varyant indeksi öne alınır. |
| `RF_IMPORTANCE` / `CHR_RF_IMPORTANCE` | Kromozom adaylarıyla kurulan RF'nin önemi. Kromozom sıralamasında eşitlikler blok önemi ve pozisyonla çözülür. |
| `GENOME_RF_IMPORTANCE` | 30.550 adayı birlikte değerlendiren son modelin önemi; nihai liste bu alana göre azalan sıralanır. |
| `GENOME_RF_RANK` | Son modeldeki özgün 1 tabanlı sıra. Eşitlikler kromozom RF önemi, blok önemi, CHR ve POS ile çözülür. |
| `VARIANT_INDEX` | Kromozom PGEN/PVAR dosyasındaki 0 tabanlı indeks; importance sırası değildir. |
| Top 1000 | Sıralamanın ilk en fazla 1.000 satırının ayrıca yazılması; istatistiksel eşik değildir. |
| Ortalama üstü eşik | `GENOME_RF_IMPORTANCE > importance.mean()`. Mevcut çalışmada 1/30.550 ≈ 0,000032733224 ve 768 SNP. Yuvarlanmış sayı yerine dosyadaki gerçek ortalama kullanılır. |
| `--min-leaf` / `--threads` | Son genom modeli için de 5 ve 1 olarak açıkça verilir. |
| `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS` | Başlatıcı her birini 1 yapar; ek sayısal kütüphane paralelliğini sınırlar. |
| `--rf-root` | Zaman damgalı sonuç klasörü. Modelin biyolojik ayarı değildir; farklı çalışmaları ayırır. |

## 5. Bu ayarları nasıl gerekçelendirmeliyiz?

10.000 SNP'lik blok ve tek iş parçacığı bellek yönetimi tercihleridir. Blok başına 50 aday son model boyutunu sınırlar. Ağaç ve yaprak ayarları mevcut başlangıç yapılandırmasıdır. Bu betiklerde grid search, nested CV veya çoklu seed kararlılık analizi olmadığı için “en iyi parametreler bulundu” denmemelidir. Daha önce farklı veri veya farklı modelde yapılan denemeler bu 6 milyon SNP'lik akışın optimizasyonu olarak sunulmamalıdır.

Değerleri değiştirmek aday havuzunu ve importance dağılımını değiştirebilir. Bu yüzden yeni modelde eski sayısal importance eşiği otomatik taşınmaz; eşik ve aday sayısı o model için yeniden belirlenir.
