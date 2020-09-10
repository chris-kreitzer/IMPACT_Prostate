## Calculation of Tumor mutational burden for MSKCC IMPACT and TCGA
## based on PMID: 31832578
## PRE-DEFINE mutation data and sample list before running function

## Input data
## MSK.data = MAF format (sequencing) containing @least
## 'Tumor_Sample_Barcode
## 'Variant_Classification
## loaded = whether data needs to be loaded, or already present in WD
## meta.data.MSK containing @least
## 'Sample.ID
## 'gene.panel (needed for mut count normalisation)

TMB = function(TCGA = NULL,
               MSK.data = NULL,
               loaded = NULL,
               meta.data.MSK){
  
  message('provide an absolute path for every input data when loaded == NULL')
  message('make sure, data input (maf) and sample list are pre-defined')
  message('make sure germline variants are excluded from data input')
  
  # define non-synonymous mutations
  nonsynonymous = c("Missense_Mutation",
                    "Silent",
                    "In_Frame_Del",
                    "Nonsense_Mutation",
                    "Splice_Site",
                    "In_Frame_Ins",
                    "Frame_Shift_Del",
                    "Nonstop_Mutation",
                    "Frame_Shift_Ins",
                    "Translation_Start_Site",
                    "Splice_Region")
  
  
  # import mutation data (TCGA)
  if(!is.null(TCGA)){
    if(is.null(loaded)){
      mutations = read.csv(TCGA, sep = '\t')
    } else {
      
      mutations = TCGA
      mutations$Sample = substr(mutations$Tumor_Sample_Barcode,
                                start = 1, stop = 15)
      
      # calculate TMB
      TMB.mutations = mutations
      
      # mutations per Mbs (TCGA);
      # excluding germline variants
      # excluding non-coding variants
      TMB_samples = as.character(unique(TMB.mutations$Sample))
      TMB_out = matrix(ncol = 3, nrow = length(TMB_samples))
      colnames(TMB_out) = c('Sample', 'TMB', 'TMB_rounded')
      
      for(i in 1:length(TMB_samples)){
        print(paste0('Currently analysing: ', TMB_samples[i]))
        TMB_calculation_sub = TMB.mutations[which(TMB.mutations$Sample == TMB_samples[i]), ]
        TMB_calculation_out = TMB_calculation_sub[TMB_calculation_sub$Variant_Classification %in% nonsynonymous, ]

        # calculate TMB for TCGA
        TMB_out = ifelse(nrow(TMB_calculation_out) == 0, 0,
                         nrow(TMB_calculation_out) / 38) # 38 Mb as suggested for TCGA (exome)
        
        TMB_out[i, 'Sample'] = TMB_samples[i]
        TMB_out[i, 'TMB'] = TMB_out
        TMB_out[i, 'TMB_rounded'] = round(TMB_out)
        
        }
      
      TMB_out = as.data.frame(TMB_out)
      TMB_out$TMB.scaled = TMB_out$TMB.scaled = (TMB_out$TMB - mean(TMB_out$TMB)) / sd(TMB_out$TMB)
      
      }
    }
    
  # ~~~~~~~~~~~~~~~
  # MSK impact data (PRAD)
  if(!is.null(MSK.data)){
    if(is.null(loaded)){
      MSK.mut = read.csv(MSK.data, sep = '\t')
      MSK.meta = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/IMPACT_PRAD.tsv', sep = '\t')
    
      } else {
      
      MSK.mut = MSK.data
      MSK.meta = meta.data.MSK
      
      ## calculate TMB MSK
      MSK_samples = as.character(unique(MSK.mut$Tumor_Sample_Barcode))
      TMB_MSK = matrix(ncol = 3, nrow = length(MSK_samples))
      colnames(TMB_MSK) = c('Sample', 'TMB', 'TMB_rounded')
      
      for(i in 1:length(MSK_samples)){
        print(paste0('Currently analysing: ', MSK_samples[i]))
        MSK_calculation = MSK.mut[MSK.mut$Tumor_Sample_Barcode == MSK_samples[i], ]
        MSK_calculation_out = MSK_calculation[MSK_calculation$Variant_Classification %in% nonsynonymous, ]
        MSK_calculation_out$gene.panel = rep(MSK.meta[MSK.meta$Sample.ID == MSK_samples[i], 'Gene.Panel'], nrow(MSK_calculation_out))

        MSK_out = ifelse(nrow(MSK_calculation_out) == 0, 0,
                         ifelse(nrow(MSK_calculation_out) != 0 && MSK_calculation_out$gene.panel == 'IMPACT468', 
                                nrow(MSK_calculation_out) / 1.139322,
                                ifelse(nrow(MSK_calculation_out) != 0 && MSK_calculation_out$gene.panel == 'IMPACT410',
                                       nrow(MSK_calculation_out) / 1.016478,
                                            ifelse(nrow(MSK_calculation_out) != 0 && MSK_calculation_out$gene.panel == 'IMPACT341',
                                                   nrow(MSK_calculation_out) / 0.896665, NA)))) # PMID: 31832578
        
        TMB_MSK[i, 'Sample'] = MSK_samples[i]
        TMB_MSK[i, 'TMB'] = MSK_out
        TMB_MSK[i, 'TMB_rounded'] = round(MSK_out)
        
        }
  
      TMB_MSK = as.data.frame(TMB_MSK)
      TMB_MSK$TMB = as.numeric(as.character(TMB_MSK$TMB))
      
      # scale TMB
      TMB_MSK$TMB.scaled = (TMB_MSK$TMB - mean(TMB_MSK$TMB)) / sd(TMB_MSK$TMB)
      
      }
  }
  
  if(!is.null(TCGA)){
    return(TMB_TCGA = TMB_out)
  } else {
    return(TMB_MSK = TMB_MSK)
  }
}

## test run
# PRAD_IMPACT_TMB = TMB(MSK.data = PRAD.IMPACT.mut.somatic, meta.data.MSK = PRAD_meta.data, loaded = 'YES')

