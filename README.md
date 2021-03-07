# BOVERI-598
### Annotating black list indels
##### Cedric Chauve, March 5, 2021

The indel caller performances on clinical and commercial dilution series is
heavily impacted by the use of a blacklist of 7 indels. It is important to have
more information about these indels to see if they should indeed be blacklisted.

The blacklist is composed of the following 7 indels, for which I could find some
references linked to diseases for five of them:

```
id                    reference
chr3:178936116:GT:C   COSMIC:2253661,dbSNP:rs587777795,ClinVar:446487
chr19:4117563:GCTT:G  nan
chr15:90631917:TC:T   COSMIC:28065520
chr12:25398297:C:CCG  nan
chr12:25380252:TC:T   COSMIC:21681080
chr1:156846307:GC:G   COMIC:74827609
chr10:89717769:TA:T   COSMIC:35243223,ClinVar:92828
```

Notes.
For `chr3:178936116:GT:C`, the provided COSMIC reference is only for the GT
deletion.  
For `chr15:90631917:TC:T`, `chr12:25380252:TC:T`, `chr1:156846307:GC:G` and
`chr10:89717769:TA:T` the COSMIC entries assume deletion of the last base of the
homopolymer while my encoding is based on deletion of the first base.

### SnpEff annotation

```
id                    annotation
chr3:178936116:GT:C   C|frameshift_variant&missense_variant|HIGH|PIK3CA|PIK3CA|transcript|NM_006218.3|protein_coding|10/21|c.1658_1659delGTinsC|p.Ser553fs|1815/9093|1658/3207|553/1068||
chr19:4117563:GCTT:G  G|conservative_inframe_deletion|MODERATE|MAP2K2|MAP2K2|transcript|NM_030662.3|protein_coding|2/11|c.154_156delAAG|p.Lys52del|410/1733|154/1203|52/400||
chr15:90631917:TC:T   T|frameshift_variant|HIGH|IDH2|IDH2|transcript|NM_002168.3|protein_coding|4/11|c.435delG|p.Thr146fs|599/1810|435/1359|145/452||,T|frameshift_variant|HIGH|IDH2|IDH2|transcript|NM_001289910.1|protein_coding|4/11|c.279delG|p.Thr94fs|366/1577|279/1203|93/400||,T|frameshift_variant|HIGH|IDH2|IDH2|transcript|NM_001290114.1|protein_coding|2/9|c.45delG|p.Thr16fs|341/1552|45/969|15/322||
chr12:25398297:C:CCG  CCG|frameshift_variant|HIGH|KRAS|KRAS|transcript|NM_033360.3|protein_coding|2/6|c.21_22insCG|p.Val8fs|213/5889|21/570|7/189||,CCG|frameshift_variant|HIGH|KRAS|KRAS|transcript|NM_004985.4|protein_coding|2/5|c.21_22insCG|p.Val8fs|213/5765|21/567|7/188||
chr12:25380252:TC:T   T|frameshift_variant|HIGH|KRAS|KRAS|transcript|NM_033360.3|protein_coding|3/6|c.205delG|p.Asp69fs|397/5889|205/570|69/189||,T|frameshift_variant|HIGH|KRAS|KRAS|transcript|NM_004985.4|protein_coding|3/5|c.205delG|p.Asp69fs|397/5765|205/567|69/188||
chr1:156846307:GC:G   G|frameshift_variant|HIGH|NTRK1|NTRK1|transcript|NM_002529.3|protein_coding|14/17|c.1753delC|p.Leu585fs|1809/2655|1753/2391|585/796||INFO_REALIGN_3_PRIME,G|frameshift_variant|HIGH|NTRK1|NTRK1|transcript|NM_001007792.1|protein_coding|14/17|c.1645delC|p.Leu549fs|1725/2571|1645/2283|549/760||INFO_REALIGN_3_PRIME,G|frameshift_variant|HIGH|NTRK1|NTRK1|transcript|NM_001012331.1|protein_coding|13/16|c.1735delC|p.Leu579fs|1791/2637|1735/2373|579/790||INFO_REALIGN_3_PRIME
chr10:89717769:TA:T   T|frameshift_variant&splice_region_variant|HIGH|PTEN|PTEN|transcript|NM_001304717.2|protein_coding|8/10|c.1319delA|p.Lys440fs|1831/8701|1319/1731|440/576||INFO_REALIGN_3_PRIME,T|frameshift_variant&splice_region_variant|HIGH|PTEN|PTEN|transcript|NM_000314.6|protein_coding|7/9|c.800delA|p.Lys267fs|1832/8702|800/1212|267/403||INFO_REALIGN_3_PRIME,T|frameshift_variant&splice_region_variant|HIGH|PTEN|PTEN|transcript|NM_001304718.1|protein_coding|7/9|c.209delA|p.Lys70fs|1946/8816|209/621|70/206||INFO_REALIGN_3_PRIME
```

### Position in amplicons

```
id                    amplicons coordinates                                                         comments
chr3:178936116:GT:C   CG001v5.0.59:chr3:178936049:178936141,CG001v4.0.9:chr3:178936036:178936133    within reverse primer in CG001v4.0.9, so not called in this amplicon
chr19:4117563:GCTT:G  CG001v3.4.62:chr19:4117483:4117588                                            close to reverse primer
chr15:90631917:TC:T   CG001v3.4.43:chr15:90631812:90631959                                          --
chr12:25398297:C:CCG  CG001v3.4.56:chr12:25398214:25398331,CG001v5.0.57:chr12:25398250:25398330     --
chr12:25380252:TC:T   CG001v3.4.59:chr12:25380228:25380356,CG001v5.0.73:chr12:25380218:25380317     close to forward primer in CG001v3.4.59
chr1:156846307:GC:G   CG001v5.0.127b:chr1:156846258:156846363,CG001v5.0.11:chr1:156846258:156846348 --
chr10:89717769:TA:T   CG001v5.0.43:chr10:89717712:89717833                                          --
```

### Flanking sequences

Five out of the 7 indels are linked to low-complexity regions:
- `chr19:4117563:GCTT:G` is part of a repeat [GCTT]CTTCTGCTG  
- `chr15:90631917:TC:T` is part of a homopolymer [TC]CCCCCC
- `chr12:25380252:TC:T` is part of a small repeat-like region [TC]CCTC  
- `chr1:156846307:GC:G` is part of a GC-rich homopolymer region GGGCC[GC]CCCCTGC
- `chr10:89717769:TA:T` is part of a homopolymer [TA]AAAAA

### Occurrences in NextSeq v51 runs

The 461 non-Blank samples of the following runs were analyzed to see if these
indels occur in them:
```
- CG001Qv51Next001,200831_NB551381_0081_AHFY5KAFX2
- CG001Qv51Next002,201002_NB551381_0087_AHFW2FAFX2
- CG001Qv51Next003,201007_NB551381_0088_AHFYYHAFX2
- CG001Qv51Next004,201009_NB551381_0089_AHFWCMAFX2
- CG001Qv51Next005,201013_NB551381_0090_AHFW2GAFX2
- CG001Qv51Next006,201016_NB551381_0091_AHFW22AFX2
- CG001Qv51Next007,201023_NB551381_0093_AHFWLJAFX2
- CG001Qv51Next008,201026_NB551381_0094_AHHH73AFX2
- CG001Qv51Next009,201027_NB551381_0095_AHHGYHAFX2
- CG001Qv51Next010,201029_NB551381_0096_AHHKM7AFX2
- CG001Qv51Next011,201102_NB551381_0097_AHHKHWAFX2
- CG001Qv51Next012,201104_NB551381_0098_AHHKK5AFX2
- CG001Qv51Next013,201112_NB551381_0099_AHJM32AFX2
- CG001Qv51Next014,201117_NB551381_0100_AHJKK3AFX2
- CG001Qv51Next015,201120_NB551381_0101_AHJLKMAFX2
- CG001Qv51Next016,201123_NB551381_0102_AHJLHYAFX2
- CG001Qv51Next017,201125_NB551381_0103_AHJLHWAFX2
- CG001Qv51Next018,201127_NB551381_0104_AHJLGGAFX2
- CG001Qv51Next019,201130_NB551381_0105_AHJLL7AFX2
- CG001Qv51Next020,201202_NB551381_0106_AHJLHLAFX2
- CG001Qv51Next021,201207_NB551381_0107_AHJLJYAFX2
- CG001Qv51Next024,210121_NB551381_0110_AHJM5HAFX2
- CG001Qv51Next025,210125_NB551381_0111_AHJLY3AFX2
- CG001Qv51Next026,210129_NB551381_0114_AHKFMWAFX2
- CG001Qv51Next032,210205_NB551381_0118_AHMHT3AFX2
- CG001Qv51Next034,210217_NB551381_0120_AHMHG7AFX2
```

The table below shows that these indels are observed (which does not imply they
are called) in most samples, with a high mean VAF, although the min. VAF can be
low, and that they generally occur in the normal female with a relatively
similar VAF (column mean_ctrl) but for the indel in chr19.

```
indel                   nb_occ  min_vaf max_vaf mean_vaf  min_ctrl  max_ctrl  mean_ctrl
chr3:178936116:GT:C     461     17.12   57.03   43.82     0.74      1.0       0.94
chr19:4117563:GCTT:G    456     0.02    4.01    1.27      0.02      1.0       0.15
chr15:90631917:TC:T     461     0.64    32.05   3.77      0.03      1.0       0.82
chr12:25398297:C:CCG    450     0.0     10.58   4.35      0.34      1.0       0.78
chr12:25380252:TC:T     460     0.0     11.99   1.74      0.02      1.0       0.75
chr1:156846307:GC:G     461     0.09    34.42   6.47      0.0       1.0       0.63
chr10:89717769:TA:T     456     0.32    3.05    0.51      0.14      1.0       0.84
```

Last I looked for indels that appear in at least 450 samples with a mean VAF at
least 0.51, i.e. have similar features than the blacklisted indels. The list is
below:

```
indel                   nb_occ  min_vaf max_vaf mean_vaf  min_ctrl  max_ctrl  mean_ctrl
chr12:25398298:CA:C     456     0.0     10.61   4.35      0.34      1.0       0.78
chr15:90631917:T:TC     456     0.32    0.96    0.58      0.54      1.0       0.9
chr17:7572962:GT:G      456     0.43    2.58    0.77      0.29      1.0       0.91
chr19:3118933:TG:T      456     0.36    0.94    0.54      0.53      1.0       0.92
chr2:29443611:TC:T      456     0.37    1.09    0.65      0.54      1.0       0.92
chr2:29445297:GGGA:G    456     0.99    2.23    1.38      0.64      1.0       0.92
chr4:1806180:G:GC       456     0.34    1.12    0.69      0.41      1.0       0.85
chr4:1806180:GC:G       456     0.46    1.18    0.78      0.54      1.0       0.9
chrX:27130922:TA:T      456     0.21    50.0    3.97      0.01      1.0       0.35
chrX:27130923:AA:T      456     0.0     34.55   1.51      0.08      1.0       0.97
```

Note that there are many groups of widespread indels. For example there are
100 different indels that are observed in all 461 samples, although most of the
time with a very low average VAF which explains why they do not appear in this
list.
