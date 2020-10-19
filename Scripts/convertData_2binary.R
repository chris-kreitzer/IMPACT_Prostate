## modify fetched data from (github) [mut, cna, fusion] and create binary matrix;
## m (samples) x n (genes) format is needed for oncoprint (complexheatmap)
## data.mutation, data.cna and data.fusion are getting converted and integrated

## input data must be in the same format as obtained from github

## values found in data.cna (not sure if there is any difference)
## "0", " 0", "-2", "2", " 2", " 0.0", " 2.0", "-2.0", "-1.5"
     
mutation.matrix = function(mut, cna){
  gene.panel = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/impact468_gene_panel.txt', sep = '\t')
  sample = data.frame(genes = 'Sample.ID')
  gene.panel = rbind(sample, gene.panel)
  
  # data somatic mutations
  data.mut = as.data.frame(mut)
  
  # data somatic copy number alterations
  data.cna = as.data.frame(cna)
  # data.cna = t(data.cna)
  # colnames(data.cna) = data.cna[1, ]
  # data.cna = as.data.frame(data.cna)
  # data.cna = data.cna[-1, ]
  # row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
  # row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
  # row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
  
  # create empty alteration matrix
  alteration_matrix = setNames(data.frame(matrix(ncol = length(gene.panel$genes), nrow = 0)), gene.panel$genes)
  
  counter = 1
  
  for(i in unique(data.mut$Tumor_Sample_Barcode)){
    data.mut.sub = data.mut[data.mut$Tumor_Sample_Barcode == i, ]
    print(i)
    
    if(nrow(data.mut.sub) != 0){
      
      for(j in unique(data.mut.sub$Hugo_Symbol)){
        if(j %in% colnames(alteration_matrix)){
          alteration_matrix[counter, j] = data.mut.sub$Variant_Classification[which(data.mut.sub$Hugo_Symbol == j)][1]
          
        }
      }
    }
    
    alteration_matrix[counter, 'Sample.ID'] = i
    counter = counter + 1
  }
  
  return(list(alteration_matrix = alteration_matrix,
              data.cna = data.cna))
}


## integrate cna and fusion data and add to mut.matrix (mutation.matrix)
## row.names == samples; colnames == genes (IMPACT 468)

onco.matrix = function(mut, cna, fusion){
  out = mutation.matrix(mut = mut, cna = cna)
  
  for(i in unique(row.names(out$data.cna))){
    data.cna.sub = out$data.cna[row.names(out$data.cna) == i, ]
    
    for(k in 1:ncol(data.cna.sub)){
      data.cna.sub[, k] = ifelse(as.numeric(noquote(data.cna.sub[, k])) < 0, -2, 
                                 ifelse(as.numeric(noquote(data.cna.sub[, k])) > 0, 2, data.cna.sub[, k]))
    }
    
    for(j in unique(colnames(data.cna.sub))){
      if(j %in% colnames(out$alteration_matrix)){
        out$alteration_matrix[out$alteration_matrix$Sample.ID == row.names(data.cna.sub), j] = ifelse(as.numeric(noquote(data.cna.sub[, j])) == 0,
                                                                                                      out$alteration_matrix[out$alteration_matrix$Sample.ID == row.names(data.cna.sub), j],
                                                                                                      ifelse(as.numeric(noquote(data.cna.sub[, j])) == -2,
                                                                                                             paste(out$alteration_matrix[out$alteration_matrix$Sample.ID == row.names(data.cna.sub), j], 'deep_deletion', sep = ';'),
                                                                                                             ifelse(as.numeric(noquote(data.cna.sub[, j])) == 2,
                                                                                                                    paste(out$alteration_matrix[out$alteration_matrix$Sample.ID == row.names(data.cna.sub), j], 'AMP', sep = ';'),
                                                                                                                    out$alteration_matrix[out$alteration_matrix$Sample.ID == row.names(data.cna.sub), j])))
      }
    }
  }
  
  # rm NA values
  # remove NA values
  k = out$alteration_matrix
  for(i in 1:ncol(k)){
    k[ , i] = ifelse(grepl('NA', k[, i]), sub(pattern = 'NA;', replacement = '', k[,i]), k[,i])
  }
  
  
  # change mutation classes
  truncating.classes = c('Nonsense_Mutation','Splice_Site','Frame_Shift_Del','Frame_Shift_Ins')
  
  for(i in 1:ncol(k)){
    k[ , i] = ifelse(grepl('Nonsense_Mutation', k[, i]), sub(pattern = 'Nonsense_Mutation', replacement = 'Truncating_Mutation', k[,i]),
                     ifelse(grepl('Splice_Site', k[, i]), sub(pattern = 'Splice_Site', replacement = 'Truncating_Mutation', k[,i]),
                            ifelse(grepl('Frame_Shift_Del', k[, i]), sub(pattern = 'Frame_Shift_Del', replacement = 'Truncating_Mutation', k[,i]),
                                   ifelse(grepl('Frame_Shift_Ins', k[, i]), sub(pattern = 'Frame_Shift_Ins', replacement = 'Truncating_Mutation', k[,i]),
                                          ifelse(grepl('In_Frame_Ins', k[, i]), sub(pattern = 'In_Frame_Ins', replacement = 'Inframe_Mutation', k[,i]),
                                                 ifelse(grepl('In_Frame_Del', k[, i]), sub(pattern = 'In_Frame_Del', replacement = 'Inframe_Mutation', k[,i]),
                                                        ifelse(grepl('3\'Flank', k[, i]), sub(pattern = '3\'Flank', replacement = 'VUS', k[,i]),
                                                               ifelse(grepl('5\'UTR', k[, i]), sub(pattern = '5\'UTR', replacement = 'VUS', k[,i]),
                                                                      ifelse(grepl('Intron', k[, i]), sub(pattern = 'Intron', replacement = 'VUS', k[,i]),
                                                                             ifelse(grepl('Splice_Region', k[, i]), sub(pattern = 'Splice_Region', replacement = 'VUS', k[,i]),
                                                                                    ifelse(grepl('Nonstop_Mutation', k[, i]), sub(pattern = 'Nonstop_Mutation', replacement = 'VUS', k[,i]),
                                                                                           ifelse(grepl('Translation_Start_Site', k[, i]), sub(pattern = 'Translation_Start_Site', replacement = 'VUS', k[,i]),
                                                                                                  ifelse(grepl('5\'Flank', k[, i]), sub(pattern = '5\'Flank', replacement = 'VUS', k[,i]), k[,i])))))))))))))
  }
  
  Somatic_alteration_matrix = k
  
  # look 
  for(i in 2:ncol(Somatic_alteration_matrix)){
    
    for(j in 1:nrow(Somatic_alteration_matrix)){
      
      Somatic_alteration_matrix[j, i] = ifelse(is.na(Somatic_alteration_matrix[j, i]), "", Somatic_alteration_matrix[j, i])
      
    }
  }
  
  # add fusion data
  data.fusion = fusion
  
  for(i in unique(data.fusion$Tumor_Sample_Barcode)){
    fusion_sub = data.fusion[data.fusion$Tumor_Sample_Barcode == i, ]
    
    if(i %in% Somatic_alteration_matrix$Sample.ID){
      for(j in unique(fusion_sub$Hugo_Symbol)){
        if(j %in% colnames(Somatic_alteration_matrix)[-1]){
          Somatic_alteration_matrix[Somatic_alteration_matrix$Sample.ID == i, j] = paste(Somatic_alteration_matrix[Somatic_alteration_matrix$Sample.ID == i, j], 'Fusion', sep = ';')
        }
      }
    }
  }
  
  return(Somatic_alteration_matrix)
}

# som = onco.matrix(mut = PRAD.IMPACT.mut.somatic, cna = PRAD.IMPACT.cna, fusion = PRAD.IMPACT.fusion)
# write.table(som, file = '~/Documents/MSKCC/03_Prostate/tmp_data/somatic_alteration_matrix.txt', sep = '\t')  



