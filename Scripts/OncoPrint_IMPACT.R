## creating meta data for oncoPrint (complexheatmap)

# data below depicts alteration (binary matrix); created with function convertData2binary
som.mat = Somatic_alteration_matrix
row.names(som.mat) = Somatic_alteration_matrix$Sample.ID
som.mat$Sample.ID = NULL
som.mat = t(som.mat)

# to do: add primary/metastasis and split top annotation (like tracerx lung)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# column annotation
# Race annotation;
meta.column = data.frame(sample.id = colnames(som.mat))
meta.column = merge(meta.column, PRAD_meta.data[, c('Sample.ID', 'Race.Category')], 
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

# add info about ERG* fusion; Fusion postive or negative
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


# TMB annotation
meta.column = merge(meta.column, PRAD_IMPACT_TMB[, c('Sample', 'TMB')],
                    by.x = 'sample.id', by.y = 'Sample', all.x = T)

# INDEL annotation
meta.column = merge(meta.column, PRAD.IMPACT.indel[,c('Sample.ID', 'rel.indel.burden')],
                    by.x = 'sample.id', 'Sample.ID', all.x = T)

# Sequencing coverage info
meta.column = merge(meta.column, PRAD_meta.data[,c('Sample.ID', 'Sample.coverage')],
                    by.x = 'sample.id', 'Sample.ID',
                    all.x = T)

# Purity info
meta.column = merge(meta.column, PRAD_meta.data[,c('Sample.ID', 'Tumor.Purity')],
                    by.x = 'sample.id', 'Sample.ID',
                    all.x = T)

# some downstream modies
row.names(meta.column) = meta.column$sample.id
meta.column$Tumor.Purity = as.numeric(as.character(meta.column$Tumor.Purity))
meta.column$sample.id = NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make top annotation
library(ComplexHeatmap)
coverage.color = circlize::colorRamp2(seq(min(meta.column$Sample.coverage), max(meta.column$Sample.coverage), length = 3),  
                                          c("white", "grey", "black"))

purity.color = circlize::colorRamp2(seq(10, 95, length = 3),
                                    c('#7c97ad', '#3e4245', '#292929'))


top = HeatmapAnnotation(TMB = anno_barplot(meta.column$TMB,
                                           which = 'column', 
                                           border = T, 
                                           height = unit(2, "cm"), 
                                           baseline = 0, 
                                           axis = T,
                                           bar_width = 2,
                                           gp = gpar(fontsize = 8,
                                                     side = 'left',
                                                     border = NA, 
                                                     fill = "red", 
                                                     lty = "blank")),
                                           
                        INDEL = anno_barplot(meta.column$rel.indel.burden,
                                             which = 'column', 
                                             border = T, 
                                             height = unit(2, "cm"), 
                                             baseline = 0, 
                                             axis = T,
                                             gp = gpar(fontsize = 8,
                                                       side = 'left')),
                        
                        ERG = anno_simple(meta.column$ERG_status,
                                          col = c('Fusion' = '#2f4e9d',
                                                  'no_Fusion' = '#b8d8e8')),
                        
                        Coverage = anno_simple(meta.column$Sample.coverage, col = coverage.color),
                        
                        Purity = anno_simple(meta.column$Tumor.Purity, col = purity.color),
                        
                        Race = anno_simple(meta.column$Race.Category, col = c("WHITE" = "#b84050",
                                                                              'REFUSED TO ANSWER' = '#87abce',
                                                                              "BLACK OR AFRICAN AMERICAN" = '#e9ac66',
                                                                              "ASIAN-FAR EAST/INDIAN SUBCONT" = '#b8d8e8',
                                                                              "OTHER" = '#d76c47')),
                        
                        gap = unit(2, "mm"),
                        
                        annotation_name_side = "left",
                        
                        annotation_legend_param = list())
                        
## look at NEJM (TRACERx lung)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# row.labels: split into pathways:
# select pathways
pathways = load(file = '~/Documents/MSKCC/00_Data/OncoPath12.Rdata')
pathways = get(pathways)

meta.rows = data.frame(genes = row.names(som.mat))
for(i in 1:nrow(meta.rows)){
  meta.rows$pathway[i] = ifelse(sum(grepl(meta.rows$genes[i], pathways)) == 1, unique(names(pathways)[grepl(meta.rows$genes[i], pathways)]),
                                'not.available')
}

row.names(meta.rows) = meta.rows[, 'genes']

# add specific category to certain genes; TP53, AR, SPOP and FOXA1, RB1, IDH
meta.rows[meta.rows$genes == 'AR', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'SPOP', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'TP53', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'FOXA1', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'RB1', 'pathway'] = 'first'
meta.rows[meta.rows$genes == 'IDH', 'pathway'] = 'first'

# set genes with no pathway annotation to NA
meta.rows[meta.rows$pathway == 'not.available', 'pathway'] = 'NA'
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
                                                         'TP53',
                                                         'NA'))

## Exclude genes with frequency < 2 %
count.matrix = as.data.frame(rowSums(som.mat != ""))
colnames(count.matrix) = 'Frequency'
count.matrix$alteration_frequency = count.matrix$Frequency / ncol(som.mat)
count.matrix = count.matrix[count.matrix$alteration_frequency >= 0.02, ]
genes.remaining = as.character(row.names(count.matrix))

meta.rows = meta.rows[row.names(meta.rows) %in% genes.remaining, ]

# exclude TMPRSS2 and ERG from meta.rows
meta.rows = meta.rows[!meta.rows$genes %in% c('TMPRSS2', 'ERG'), ]

meta.rows$genes = NULL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heatmap Body: OncoPrint
som.mat = Somatic_alteration_matrix
row.names(som.mat) = Somatic_alteration_matrix$Sample.ID
som.mat$Sample.ID = NULL
som.mat = t(som.mat)

# subset data frame to those observatons where meta data is available
som.mat.print = som.mat[row.names(som.mat) %in% row.names(meta.rows), colnames(som.mat) %in% row.names(meta.column)]
# exclude TMP and ERG (as those display a own group)
som.mat.print = som.mat.print[!row.names(som.mat.print) %in% c('TMPRSS2', 'ERG'), ]

## change global options (padding between top annotation and heatmap)
ht_opt$COLUMN_ANNO_PADDING = unit(6, "mm")

x.first = ComplexHeatmap::oncoPrint(som.mat.print,
                          name = 'race',
                          row_names_gp = gpar(fontsize = 7),
                          row_names_side = 'left',
                          pct_side = FALSE,
                          alter_fun = alter_fun, 
                          column_gap = unit(2.5, 'mm'),
                          col = col,
                          top_annotation = top,
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
                          width = unit(18, 'cm'),
                          height = unit(18, 'cm'))



# with draw we can make additional modification on the output

draw(x.first, 
     padding = unit(c(1, 1, 1, 1), "cm"),
     heatmap_legend_list = list(
  Legend(labels = c('Fusion', 'no_Fusion'), legend_gp = gpar(fill = c('Fusion' = '#2f4e9d',
                                                                      'no_Fusion' = '#b8d8e8')), 
         title = "ERG_Status"),
  Legend(col_fun = coverage.color, title = "Coverage", at = c(0, 1000, 2000, 2500)),
  Legend(col_fun = purity.color, title = "Purity", at = c(10, 50, 95)),
  Legend(labels = c('WHITE', 'REFUSED TO ANSWER', 'BLACK OR AFRICAN AMERICAN' , 'ASIAN-FAR EAST/INDIAN SUBCONT', 'OTHER'),
         legend_gp = gpar(fill = c("WHITE" = "#b84050",
                 'REFUSED TO ANSWER' = '#87abce',
                 "BLACK OR AFRICAN AMERICAN" = '#e9ac66',
                  "ASIAN-FAR EAST/INDIAN SUBCONT" = '#b8d8e8',
                  "OTHER" = '#d76c47')))
))


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

