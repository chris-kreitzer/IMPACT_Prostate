# IMPACT_Prostate

This repo contains scripts used to analyse MSK_IMPACT Prostate samples;   
Data fetched from https://github.mskcc.org/knowledgesystems/dmp-2020/tree/master/mskimpact (controlled access)   
date of data retrieval: 09/01/2020
***
#### data: MSKCC github enterprise
- data_mutations_extended.txt
- data_CNA.txt
- data_fusions.txt
- pathway annotation: OncoPath12; Bastien

#### sample selection (top plot): ONE sample per patient   
- if just one sample --> take it
- if multiple samples, priority is given to primary sample
- if multiple samples and no primary available: i) max tumor purity and ii) coverage of sample
- all panels are included (ACCESS129 IMPACT341 IMPACT410 IMPACT468)
- n = 2,588 remaining

#### mutation classification:
- truncating classes: 'Nonsense_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins'
- InFrame_Mutations: 'In_Frame_Ins', 'In_Frame_Del'
- VUS: '3\'Flank', '5\'UTR', 'Intron', 'Splice_Region', 'Nonstop_Mutation', 'Translation_Start_Site', '5\'Flank'


