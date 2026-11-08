#!/usr/bin/env python3
"""Eski PC5 GWAS sonuçlarını kromozom bazında özetler; GWAS çalıştırmaz."""
import csv,gzip,math
from pathlib import Path
base=Path('/volume1/homes/tuseb/new_data_to_be_used/GWAS_404_vs_SLD_allchr')
src=base/'combined_results/GWAS_allchr_corrected_PC5_all_results.tsv.gz'
out=base/'combined_results/GWAS_chr_summary_for_excel.tsv'
rows={str(c):dict(CHR=c,INPUT_VARIANTS=0,RESULT_ROWS=0,VALID_P=0,INVALID_P=0,P_LT_005=0,P_LT_5E8=0,BONF_LT_005=0,FDR_LT_005=0,BONF_LT_1=0,FDR_LT_1=0,MIN_P=None) for c in range(1,23)}
for chrom,r in rows.items():
    p=base/f'per_chr/chr{chrom}/SLD_EA_chr{chrom}_merged.pvar'
    with p.open() as f: r['INPUT_VARIANTS']=sum(bool(l.strip()) and not l.startswith('#') for l in f)
def numeric(v):
    try: x=float(v)
    except (ValueError,TypeError): return None
    return x if math.isfinite(x) and 0<=x<=1 else None
with gzip.open(src,'rt') as f:
    rd=csv.DictReader(f,delimiter='\t')
    names=rd.fieldnames or []
    chrcol=next((c for c in ['CHROM','#CHROM','CHR'] if c in names),None)
    assert chrcol and all(c in names for c in ['P','BONFERRONI','FDR_BH']),f'Beklenmeyen sütunlar: {names}'
    for x in rd:
        if 'TEST' in x and x['TEST']!='ADD': continue
        chrom=x[chrcol]
        if chrom.startswith('chr'): chrom=chrom[3:]
        assert chrom in rows,f'Beklenmeyen kromozom: {chrom}'
        r=rows[chrom]; r['RESULT_ROWS']+=1
        p=numeric(x['P'])
        if p is None:
            r['INVALID_P']+=1; continue
        r['VALID_P']+=1; r['P_LT_005']+=p<0.05; r['P_LT_5E8']+=p<5e-8
        r['MIN_P']=p if r['MIN_P'] is None else min(p,r['MIN_P'])
        for column,field in [('BONFERRONI','BONF_LT_005'),('FDR_BH','FDR_LT_005')]:
            q=numeric(x[column])
            assert q is not None,f'Geçerli P için eksik düzeltme: chr{chrom}'
            r[field]+=q<0.05
            r['BONF_LT_1' if column=='BONFERRONI' else 'FDR_LT_1']+=q<1
assert all(r['RESULT_ROWS']>0 for r in rows.values()),'Bir veya daha fazla kromozom sonucu eksik'
with out.open('w') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows['1']),delimiter='\t'); wr.writeheader(); wr.writerows(rows.values())
print('Kaynak:',src)
print('Düzeltmeler: mevcut tüm-genom BONFERRONI ve FDR_BH sütunları; chr içinde yeniden hesaplanmadı.')
for col in ['INPUT_VARIANTS','RESULT_ROWS','VALID_P','P_LT_005','P_LT_5E8','BONF_LT_005','FDR_LT_005','BONF_LT_1','FDR_LT_1']:
    print(col, sum(r[col] for r in rows.values()))
print('Excel için gönderilecek dosya:',out)
