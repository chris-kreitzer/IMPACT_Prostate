## INDEL LOAD

INDEL_Annotation = function(mut.data){
  
  data.mut = mut.data
  input.samples = unique(data.mut$Tumor_Sample_Barcode)
  
  # INDEL LOAD
  # INDEL == INS && DEL [!= DNP | ONP]
  PRAD.IMPACT.indel = matrix(nrow = length(input.samples), ncol = 4)
  colnames(PRAD.IMPACT.indel) = c('Sample.ID', 'Ins', 'Del', 'Snps')
  
  for(i in 1:length(input.samples)){
    data.sub = data.mut[data.mut$Tumor_Sample_Barcode == input.samples[i], ]
    
    if(nrow(data.sub) != 0){
      # make frequency table
      freq.table = table(data.sub$Variant_Type)
      Ins = freq.table['INS'][[1]]
      Del = freq.table['DEL'][[1]]
      Snps = ifelse(length(freq.table[!dimnames(freq.table)[[1]] %in% c('INS', 'DEL')]) == 0, 0,
                    freq.table[!dimnames(freq.table)[[1]] %in% c('INS', 'DEL')][[1]])
    } else {
      Ins = 0
      Del = 0
      Snps = 0
    }

    PRAD.IMPACT.indel[i, 'Sample.ID'] = input.samples[i]
    PRAD.IMPACT.indel[i, 'Ins'] = Ins
    PRAD.IMPACT.indel[i, 'Del'] = Del
    PRAD.IMPACT.indel[i, 'Snps'] = Snps
  }
  
  PRAD.IMPACT.indel = as.data.frame(PRAD.IMPACT.indel)
  PRAD.IMPACT.indel$Ins = ifelse(is.na(PRAD.IMPACT.indel$Ins), 0, PRAD.IMPACT.indel$Ins)
  PRAD.IMPACT.indel$Del = ifelse(is.na(PRAD.IMPACT.indel$Del), 0, PRAD.IMPACT.indel$Del)
  PRAD.IMPACT.indel$Ins = as.numeric(as.character(PRAD.IMPACT.indel$Ins))
  PRAD.IMPACT.indel$Del = as.numeric(as.character(PRAD.IMPACT.indel$Del))
  PRAD.IMPACT.indel$Snps = as.numeric(as.character(PRAD.IMPACT.indel$Snps))
  
  # convert to log10 (easier to handel)
  PRAD.IMPACT.indel$Ins.log = log10(PRAD.IMPACT.indel$Ins + 1)
  PRAD.IMPACT.indel$Del.log = log10(PRAD.IMPACT.indel$Del + 1)
  
  # relative INDEL proporation: https://doi.org/10.1016/S1470-2045(17)30516-8
  PRAD.IMPACT.indel$rel.indel.burden = (PRAD.IMPACT.indel$Ins + PRAD.IMPACT.indel$Del) / 
    (PRAD.IMPACT.indel$Ins + PRAD.IMPACT.indel$Del + PRAD.IMPACT.indel$Snps)
  
  return(PRAD.IMPACT.indel)
  
}

# sample run
# PRAD.IMPACT.indel = INDEL_Annotation(mut.data = PRAD.IMPACT.mut)


