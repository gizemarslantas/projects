# Proje taslağı için GWAS yöntemi

## Kodla doğrulanabilen yöntem metni

404/EA ve SLD kohortlarının genotip dozajları, otozomlar (1–22) ayrı ayrı işlenerek analiz edilmiştir. Çalışma 1.580 EA olgusu ve 771 SLD kontrolü olmak üzere 2.351 bireyi içermektedir. Fenotip dosyasında kontrol=1 ve olgu=2 kodlanmıştır. Gen anotasyonu için kullanılan koordinat sistemi GRCh37/hg19'dur.

PLINK 2 ile tek bazlı A/C/G/T SNP'leri ve en fazla iki alelli kayıtlar seçilmiş; varyant kimlikleri kromozom, pozisyon, REF ve ALT biçiminde düzenlenmiştir. Yinelenen kimlikler `--rm-dup exclude-all` ile dışlanmıştır. İki kohort arasında aynı alel sırasındaki ve REF/ALT sırası ters kayıtlar belirlenmiş, ortak test edilebilir varyantlar çıkarılarak alel düzenleri hizalanmıştır. Bu işlem indel normalizasyonu veya genom sürümü dönüşümü içermez. Hazırlık kodu, proje sırasında tespit edilmiş SLD doz yönü sorununa özel bir düzeltme uygular; bu adımın yeni veya zaten düzeltilmiş bir veri setine genellenmesi uygun değildir.

Varyant sayısı ve kimlik/REF/ALT düzeni karşılaştırılmış; düzeltme sonrası kohort alel frekansları arasında korelasyon >0,90 ve medyan mutlak fark <0,10 koşulları aranmıştır. Bu toplu kontroller, her varyantın doğruluğunu veya kohortlar arasındaki teknik farkların tümüyle giderildiğini tek başına kanıtlamaz.

Dozajları korumak amacıyla PLINK 2 çıktıları DS alanı içeren sıkıştırılmış VCF biçiminde dışa aktarılmış, bcftools ile bireyler birleştirilmiş ve DS değerleri yeniden PGEN biçimine alınmıştır. Her kromozomda birleşik birey sayısı ve beklenen ortak varyant sayısı doğrulanmıştır.

İkili fenotip için PLINK 2 `--glm firth-fallback hide-covar` kullanılmıştır. Model, additive alel dozajı etkisini ve popülasyon yapısını temsil eden PC1–PC5 kovaryatlarını içerir. `firth-fallback`, gerektiğinde Firth regresyonuna geçiş sağlar; bütün testlerin zorunlu olarak Firth ile yapıldığı anlamına gelmez. Bu kayıtlı GWAS komutunda yaş ve cinsiyet kovaryatları belirtilmemiştir. RF'nin güncel PC'siz sürümü ayrı bir analizdir.

Kromozom sonuçları birleştirilirken yalnızca `TEST=ADD` satırları ve sonlu, 0–1 aralığındaki p-değerleri korunmuştur. Toplam 6.003.495 ADD testinden 6.003.487'si geçerli p-değeri üretmiş; 8 test düzeltme hesabının dışında kalmıştır. Geçersiz testlerin nedenleri bu birleştirme kodundan belirlenemez; ilgili hata kodları ve loglar incelenmelidir.

Çoklu test düzeltmesi kromozom başına ayrı ayrı değil, 6.003.487 geçerli testin tamamı üzerinden uygulanmıştır. Bonferroni değeri `min(1, P × m)` olarak; Benjamini–Hochberg FDR değeri sıralı p-değerlerinde `P × m / sıra` ve sondan monoton kümülatif minimum kullanılarak hesaplanmıştır. Bonferroni için aile düzeyinde alfa=0,05; FDR için q<0,05 ve ayrıca P<5×10^-8 eşikleri raporlanmıştır. Tüm sonuçlar p-değerine göre artan sırada kaydedilmiş, ilk 1.000 varyant ayrıca dışa aktarılmıştır.

Genomik kontrol katsayısı, p-değerlerinden elde edilen 1 serbestlik dereceli ki-kare istatistiklerinin medyanının 0,4549364231195727'ye bölünmesiyle hesaplanmıştır. Kod bu katsayıyı raporlar; p-değerlerine genomik kontrol düzeltmesi uygulamaz.

## Taslağa eklemeden önce tamamlanacak bilgiler

- Girdi genotip QC ve imputasyon kalite eşikleri, panel ve yazılım sürümleri ana üç betikte tanımlı değildir. Ek MAF %5/%10 denemeleri ana GWAS filtresi olarak yazılmamalıdır.
- PC1–PC5 dosyasının üretildiği birey/varyant kümesi, LD pruning ayarları ve PCA komutu bu paket içinde yoktur.
- `EA_c2_1580.keep` ve fenotip tablosunun üretimindeki dahil etme/dışlama ölçütleri ek kaynak gerektirir.
- Olgu/kontrol etiketleri aynı zamanda farklı kohortları ayırır. Popülasyon, ölçüm ve imputasyon/batch farkları sonuçları etkileyebilir; PC5 bunların tamamını giderdiğinin kanıtı değildir.
- Kaynak kodlar geçmiş sohbetten alınmıştır. Son NAS betikleriyle hash karşılaştırması ve kullanılan PLINK/bcftools/Python sürümlerinin kaydı henüz yapılmamıştır.

Bu metin mevcut kodu ve gözlenen analiz sayılarını açıklar; eksik protokol ayrıntılarını veya nedensel biyolojik sonuçları varsaymaz.
