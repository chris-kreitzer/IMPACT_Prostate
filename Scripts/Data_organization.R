## MSK IMPACT: Data files and sample organization:
## Specifically working with IMPACT-PRAD samples: 


## raw data:
## obtained from https://github.mskcc.org/knowledgesystems/dmp-2020/tree/master/mskimpact
## retrieval date: 09/01/2020
## data stored locally


PRAD_meta.data = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/IMPACT_PRAD.tsv', sep = '\t')

## mutations
IMPACT_mutations = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_mutations_extended.txt', sep = '\t')
PRAD.IMPACT.mut = IMPACT_mutations[IMPACT_mutations$Tumor_Sample_Barcode %in% PRAD.samples$sample.ID, ]
PRAD.IMPACT.mut.somatic = PRAD.IMPACT.mut[which(PRAD.IMPACT.mut$Mutation_Status == 'SOMATIC'), ]

## cna
IMPACT_cna = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_CNA.txt', sep = '\t')
data.cna = as.data.frame(IMPACT_cna)
data.cna = t(data.cna)
colnames(data.cna) = data.cna[1, ]
data.cna = as.data.frame(data.cna)
data.cna = data.cna[-1, ]
row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
PRAD.IMPACT.cna = data.cna[row.names(data.cna) %in% PRAD.samples$sample.ID, ]

## fusion
IMPACT_fusion = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_fusions.txt', sep = '\t')
PRAD.IMPACT.fusion = IMPACT_fusion[IMPACT_fusion$Tumor_Sample_Barcode %in% PRAD.samples$sample.ID, ]

## structural variants:
Impact_SV = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_Structural_variant.txt', sep = '\t')
PRAD.IMPACT.SV = Impact_SV[Impact_SV$SampleId %in% PRAD.samples$sample.ID, ]



## function which selects ONE sample per patient
## if just one sample per patient --> use it (regardless of primary or metastasis)
## if multiple samples, priority is given to i) primaries, ii) tumor purity and iii) coverage

sample_selection = function(data){
  data.analysis = data.frame()
  data.in = as.data.frame(data)
  
  for(i in unique(data.in$Patient.ID)){
    # check how many samples per patient
    data.sub = data.in[which(data.in$Patient.ID == i), ]
    
    # one sample; go for it
    if(nrow(data.sub) == 1){
      chosen.sample = data.sub$Sample.ID
      } 
    
    # if multiple samples and primary available, choose primary
    else if (nrow(data.sub) != 1 && 'Primary' %in% data.sub$Sample.Type){
      chosen.sample = data.sub$Sample.ID[which(data.sub$Sample.Type == 'Primary' && data.sub$Sample.Class == 'Tumor')]
      chosen.sample = ifelse(length(chosen.sample) != 1, 
                             data.sub$Sample.ID[which.max(data.sub$Tumor.Purity)],
                             chosen.sample)
    }
    
    # if multiple samples and no primary, choose sample with max purity
    else if (nrow(data.sub) != 1 && 'Primary' != data.sub$Sample.Type){
      chosen.sample = data.sub$Sample.ID[which.max(data.sub$Tumor.Purity)]
    }
    
    data.out = data.frame(patient.ID = rep(i, length(chosen.sample)),
                          sample.ID = chosen.sample,
                          panel = data.sub$Gene.Panel[which(data.sub$Sample.ID == chosen.sample)])
    
    data.analysis = rbind(data.analysis, data.out)
  }
  
  # clean up the mess
  rm(data.sub)
  rm(chosen.sample)
  rm(data.out)
  
  return(data.analysis)
}

PRAD.samples = sample_selection(data = PRAD_meta.data)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## save a comprehensive RData file for all components needed
save(PRAD.samples, 
     PRAD.IMPACT.mut, 
     PRAD.IMPACT.mut.somatic, 
     PRAD.IMPACT.cna, 
     PRAD.IMPACT.fusion, 
     PRAD.IMPACT.SV,
     PRAD_IMPACT_TMB,
     PRAD.IMPACT.indel,
     PRAD_meta.data,
     file = '~/Documents/MSKCC/03_Prostate/tmp_data/data.files.RData')

