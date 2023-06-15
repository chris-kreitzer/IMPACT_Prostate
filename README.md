# IMPACT_Prostate

This repo contains scripts used to analyse MSK_IMPACT Prostate samples;   
Data fetched from https://github.mskcc.org/knowledgesystems/dmp-2020/tree/master/mskimpact (controlled access)   
date of data retrieval: 09/01/2020
***
#### data: MSKCC github enterprise
- data_mutations_extended.txt
- data_CNA.txt (unique values: "0", " 0", "-2", "2", " 2", " 0.0", " 2.0", "-2.0", "-1.5")
- data_fusions.txt
- IMPACT 468 is considered (working with those genes)
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

#### Tumor mutational burden:
TMB is defined as the number of somatic, coding, base substitution (including INDEL mutations) per megabase of genome examined). Variants outside coding regions are neglected. Synonymous mutations are counted. Their presence is a signal of mutational processes that will aslo have resulted in nonsynonymous mutations and neoantigens elsewhere in the genome. [PMID:28420421](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395719/)

For MSKCC samples, the mutation count was divided by   
**0.896665**, **1.016478**, and **1.139322** Mb for the 341-, 410- and 468-gene panels, respectively. [PMID: 31832578](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6907021/)

#### INDEL load
INDEL == 'Ins' & 'Del'   
rel.INDEL.load = ('INS' + 'DEL') / ('INS' + 'DEL' + 'SNP' + 'DNP' + 'ONP')

'DNP', 'ONP' are not considered as INDELS

#### OncoPrint:
- data fetched from GitHub Enterprise MSKCC dmp2020 (see above)
- only genes with alteration frequency > 2% are shown
- TMB/INDEL load (top annotation) are calculated like points mentioned above
- SOMATIC alterations are excluded from the OncoPrint

Harmonization of TMB!
