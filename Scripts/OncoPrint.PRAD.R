## OncoPrint of IMPACT_PRAD

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

# function to draw geometric shapes for oncoprint
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

## change global options (padding between top annotation and heatmap)
ht_opt$COLUMN_ANNO_PADDING = unit(6, "mm")
ComplexHeatmap::oncoPrint(som.mat.print,
                          name = 'race',
                          row_names_gp = gpar(fontsize = 7),
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
                          column_title = NULL,
                          width = unit(25, 'cm'),
                          height = unit(25, 'cm'))


# top annotation
top_annox <- HeatmapAnnotation(df = meta.column,
                               col = list(Race.Category = c("WHITE" = "#b84050",
                                                            'REFUSED TO ANSWER' = '#87abce',
                                                            "BLACK OR AFRICAN AMERICAN" = '#e9ac66',
                                                            "ASIAN-FAR EAST/INDIAN SUBCONT" = '#b8d8e8',
                                                            "OTHER" = '#d76c47'),
                                          ERG_status = c('Fusion' = '#2f4e9d',
                                                         'no_Fusion' = '#b8d8e8')),
                               gap = unit(c(2), "mm"),
                               simple_anno_size = unit(2.5, "mm"),
                               annotation_label = c("Race", "ERG Status"),
                               border = c(Race.Category = TRUE))

