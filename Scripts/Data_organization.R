# Data files and sample organization:
# IMPACT-PRAD samples: starting with exploratory data analysis

# raw data
Impact_PRAD = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/IMPACT_PRAD.tsv', sep = '\t')
Impact_mut = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_mutations_extended.txt', sep = '\t')
data.cna = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_CNA.txt', sep = '\t')
data.fusion = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_fusions.txt', sep = '\t')

# selection hierachy for selecting samples considered for further analysis
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

PRAD.samples = sample_selection(data = Impact_PRAD)
# write.table(PRAD.samples, file = '~/Documents/MSKCC/03_Prostate/tmp_data/PRAD_IMPACT_Samples.txt', sep = '\t', quote = F, row.names = F)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create data for oncoprint/plot
# working exclusively with somatic mutations first
Impact_mut = Impact_mut[Impact_mut$Tumor_Sample_Barcode %in% Impact_Prad_Samples$sample.ID, ]
data.mut = Impact_mut[which(Impact_mut$Mutation_Status == 'SOMATIC'), ]

# convert mutations into binary matrix
mutation.matrix = function(mut, cna){
  gene.panel = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/impact468_gene_panel.txt', sep = '\t')
  sample = data.frame(genes = 'Sample.ID')
  gene.panel = rbind(sample, gene.panel)
  
  # data somatic mutations
  data.mut = as.data.frame(mut)
  
  # data somatic copy number alterations
  data.cna = as.data.frame(cna)
  data.cna = t(data.cna)
  colnames(data.cna) = data.cna[1, ]
  data.cna = as.data.frame(data.cna)
  data.cna = data.cna[-1, ]
  row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
  row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
  row.names(data.cna) = sub(pattern = '\\.', replacement = '-', x = row.names(data.cna))
  
  # data.fusion = as.data.frame(fusion)
  
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add SCNA data to SNVs
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

# run function
matrix.out = onco.matrix(mut = data.mut, cna = data.cna, fusion = data.fusion)



# write.table(matrix.out, file = '~/Documents/MSKCC/03_Prostate/tmp_data/somatic_alteration_matrix.txt', sep = '\t', quote = F)

# # only work with IMPACT-PRAD data.cna
# colnames(data.cna) = sub(pattern = '\\.', replacement = '-', colnames(data.cna))
# sample.index = which(colnames(data.cna) %in% Impact_Prad_Samples$sample.ID)
# k = data.cna[, c(1, sample.index)]
# data.cna = k
# colnames(data.cna) = sub(pattern = '-', replacement = '\\.', x = colnames(data.cna))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ComplexHeatmap: oncoplot of IMPACT samples
## modify data:
row.names(Somatic_alteration_matrix)  = Somatic_alteration_matrix[, 1]
Somatic_alteration_matrix = Somatic_alteration_matrix[, -1]
som.mat = t(as.matrix(Somatic_alteration_matrix))


library(ComplexHeatmap)
# define color-space for alterations:
col = c("deep_deletion" = "#2f4e9d", 
        "AMP" = "#a20414",
        'Fusion' = '#960b79',
        "Missense_Mutation" = "#008000",
        'Truncating_Mutation' = '#b88e38',
        'Inframe_Mutation' = '#000000',
        'VUS' = '#88d332')

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = "#e6e6e6", col = NA))
    },
  
  # big blue
  deep_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["deep_deletion"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Fusion"], col = NA))
  },
  
  # small green
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  
  # small black
  Inframe_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Inframe_Mutation"], col = NA))
  },
  
  # small brown 
  Truncating_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Truncating_Mutation"], col = NA))
  },
  
  VUS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["VUS"], col = NA))
  }
  
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make meta-data frame for heatmap annotation
# select pathways
pathways = load(file = '~/Documents/MSKCC/00_Data/OncoPath12.Rdata')
pathways = get(pathways)

meta.rows = data.frame(genes = row.names(som.mat))
for(i in 1:nrow(meta.rows)){
  meta.rows$pathway[i] = ifelse(sum(grepl(meta.rows$genes[i], pathways)) == 1, unique(names(pathways)[grepl(meta.rows$genes[i], pathways)]),
                                'not.available')
}

row.names(meta.rows) = meta.rows[,'genes']

# add specific category to certain genes; TP53, AR, SPOP and FOXA1
meta.rows[meta.rows$genes == 'AR', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'SPOP', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'TP53', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'FOXA1', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'RB1', 'pathway'] = 'first'

# exclude those not available
meta.rows = meta.rows[meta.rows$pathway != 'not.available', ]
meta.rows$pathway = factor(meta.rows$pathway, levels = c('first',
                                                         'RTK_RAS',
                                                         'Epigenetic',
                                                         'DDR',
                                                         'PI3K',
                                                         'Cell_Cycle',
                                                         'NOTCH',
                                                         'WNT',
                                                         'TGF-Beta',
                                                         'HIPPO',
                                                         'MYC',
                                                         'NRF2',
                                                         'TP52'))

## exclude genes with frequency < 1 %
count.matrix = as.data.frame(rowSums(som.mat != ""))
colnames(count.matrix) = 'Frequency'
count.matrix$alteration_frequency = count.matrix$Frequency / ncol(som.mat)
count.matrix = count.matrix[count.matrix$alteration_frequency >= 0.01, ]
genes.remaining = as.character(row.names(count.matrix))

meta.rows = meta.rows[row.names(meta.rows) %in% genes.remaining, ]
meta.rows$genes = NULL

# ~~
# column annotation;
meta.column = data.frame(sample.id = colnames(som.mat))
meta.column = merge(meta.column, Impact_PRAD[, c('Sample.ID', 'Race.Category')], 
                    by.x = 'sample.id', by.y = 'Sample.ID', all.x = T)
row.names(meta.column) = meta.column$sample.id
meta.column$Race.Category = ifelse(is.na(meta.column$Race.Category), 'OTHER',
                                   ifelse(meta.column$Race.Category == 'PT REFUSED TO ANSWER', 'REFUSED TO ANSWER',
                                          ifelse(meta.column$Race.Category == 'NO VALUE ENTERED', 'OTHER',
                                                 ifelse(meta.column$Race.Category == 'UNKNOWN', 'OTHER',
                                                        ifelse(meta.column$Race.Category == 'NATIVE HAWAIIAN OR PACIFIC ISL', 'OTHER',
                                                               meta.column$Race.Category)))))

meta.column$Race.Category = factor(meta.column$Race.Category, levels = c('WHITE',
                                                                         'BLACK OR AFRICAN AMERICAN',
                                                                         'ASIAN-FAR EAST/INDIAN SUBCONT',
                                                                         'REFUSED TO ANSWER',
                                                                         'OTHER'))
# add info if Sample.id has ERG Fusion or not
meta.column$ERG_status = NA
for(i in colnames(som.mat)){
  if(i %in% meta.column$sample.id){
    if(grepl('Fusion', som.mat[row.names(som.mat) == 'ERG', i])){
      meta.column[meta.column$sample.id == i, 'ERG_status'] = 'Fusion'
    }
    else {
      meta.column[meta.column$sample.id == i, 'ERG_status'] = 'no_Fusion'
    }
  } 
}

meta.column$sample.id = NULL

## subset original matrix to meta data
som.mat.print = som.mat[row.names(som.mat) %in% row.names(meta.rows), colnames(som.mat) %in% row.names(meta.column)]

ComplexHeatmap::oncoPrint(som.mat.print,
                          name = 'race',
                          row_names_gp = gpar(fontsize = 8),
                          row_names_side = 'left',
                          pct_side = 'right',
                          alter_fun = alter_fun, 
                          column_gap = unit(2.5, 'mm'),
                          col = col,
                          top_annotation = top_annox,
                          right_annotation = NULL,
                          remove_empty_columns = T,
                          remove_empty_rows = T, 
                          row_split = meta.rows$pathway,
                          column_split = meta.column$Race.Category,
                          row_gap = unit(2.5, "mm"),
                          row_title_rot = 0,
                          row_title_gp = gpar(col = 'black'),
                          use_raster = F,
                          column_title = '',
                          width = unit(25, 'cm'),
                          height = unit(25, 'cm'))



top_annox <- HeatmapAnnotation(df = meta.column,
                               col = list(Race.Category = c("WHITE" = "#b84050",
                                                            'REFUSED TO ANSWER' = '#87abce',
                                                            "BLACK OR AFRICAN AMERICAN" = '#e9ac66',
                                                            "ASIAN-FAR EAST/INDIAN SUBCONT" = '#b8d8e8',
                                                            "OTHER" = '#d76c47'),
                                          ERG_status = c('Fusion' = '#2f4e9d',
                                                         'no_Fusion' = '#b8d8e8')),
                               gap = unit(c(2.5, 5, 10), "mm"),
                               simple_anno_size = unit(2.5, "mm"),
                               annotation_label = c("Race", "ERG Status"),
                               border = c(Race.Category = TRUE))









# mutations Mbs per 

ComplexHeatmap::oncoPrint(white.mat[,1:1758],
                          name = 'race',
                          row_names_gp = gpar(fontsize = 6),
                          row_names_side = 'left',
                          pct_side = 'right',
                          alter_fun = alter_fun, 
                          top_annotation = top_annotation,
                          
                          
                          col = col,
                          remove_empty_columns = T,
                          remove_empty_rows = T, 
                          row_split = df2$category,
                          row_gap = unit(5, "mm"),
                          row_title_rot = 0,
                          row_title_gp = gpar(col = c("red", "blue"), font = 1:2))
                          


top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                   df = df, col = list(Race.Category = c("WHITE" = "#6C9E09")),
                                   simple_anno_size = unit(2.5, "mm"))










# playing around: groups for race
white = test[,colnames(test) %in% anno$Sample.ID[which(anno$Race.Category == 'WHITE')]]
black = test[,colnames(test) %in% anno$Sample.ID[which(anno$Race.Category == 'BLACK OR AFRICAN AMERICAN')]]

p.white = ComplexHeatmap::oncoPrint(white,
                          row_names_gp = gpar(fontsize = 8),
                          row_names_side = 'left',
                          pct_side = 'right',
                          alter_fun = alter_fun,
                          col = col,
                          remove_empty_columns = F,
                          remove_empty_rows = F,
                          top_annotation = rowAnnotation(Race = rep(1, 42),
                                                         col = list(Race = u)),
                          right_annotation = NULL)



p.white

p.black = ComplexHeatmap::oncoPrint(black,
                                    row_names_gp = gpar(fontsize = 8),
                                    row_names_side = 'left',
                                    pct_side = 'right',
                                    alter_fun = alter_fun,
                                    col = col,
                                    remove_empty_columns = F,
                                    remove_empty_rows = F,
                                    top_annotation = ha)



df = data.frame(type = c(rep("Black", 2)))
ha = HeatmapAnnotation(df = df)



p.white + p.black


















truncating_classes <- c('Nonsense_Mutation','Splice_Site','Frame_Shift_Del','Frame_Shift_Ins')

important_classes <- c('Missense_Mutation','Nonsense_Mutation','Splice_Site','In_Frame_Del','In_Frame_Ins',
                       'Frame_Shift_Del','Frame_Shift_Ins','TERT promoter','Translation_Start_Site')


# IMPACT - PROSTATE:


Impact_SV = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/data_Structural_variant.txt', sep = '\t')


# Sequencing over time
library(readxl)
sequencing = read_excel(path = '~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/Sequencing_Dates_MSKIMPACT.xlsx')



tcga = read.csv('~/Documents/MSKCC/CPNA_analysis/TCGA/RawData/prad_tcga_pan_can_atlas_2018_clinical_data.tsv', sep = '\t')
