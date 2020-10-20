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


som.mat.df = as.data.frame(som.mat)
som.mat.df.out = as.data.frame(colSums(som.mat.df != ""))
colnames(som.mat.df.out) = 'freq'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make top annotation
library(ComplexHeatmap)
coverage.color = circlize::colorRamp2(seq(min(meta.column$Sample.coverage), max(meta.column$Sample.coverage), length = 3),  
                                          c("white", "grey", "black"))

purity.color = circlize::colorRamp2(seq(10, 95, length = 3),
                                    c('#7c97ad', '#3e4245', '#292929'))



top = HeatmapAnnotation(
                        mutations = anno_barplot(
                          som.mat.df.out$freq,
                          bar_width = 1,
                          gp = gpar(col = 'black'),
                          border = T,
                          height = unit(2, 'cm')),
                        
                        # ERG Status
                        ERG = anno_simple(meta.column$ERG_status,
                                          col = c('Fusion' = '#2f4e9d',
                                                  'no_Fusion' = '#b8d8e8'),
                                          height = unit(3, 'mm')),
                        # Coverage
                        Coverage = anno_simple(meta.column$Sample.coverage, col = coverage.color,
                                               height = unit(3, 'mm')),
                        
                        # Purity
                        Purity = anno_simple(meta.column$Tumor.Purity, col = purity.color,
                                             height = unit(3, 'mm')),
                        
                        # Race annotation
                        Race = anno_simple(meta.column$Race.Category, col = c("WHITE" = "#b84050",
                                                                              'REFUSED TO ANSWER' = '#87abce',
                                                                              "BLACK OR AFRICAN AMERICAN" = '#e9ac66',
                                                                              "ASIAN-FAR EAST/INDIAN SUBCONT" = '#b8d8e8',
                                                                              "OTHER" = '#d76c47'),
                                           height = unit(3, "mm")),
                        
                        gap = unit(2, "mm"),
                        
                        annotation_name_side = "left",
                        
                        annotation_legend_param = list())
                        


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
meta.rows = data.frame(genes = row.names(som.mat))
meta.rows$pathway = NA

# add specific category to certain genes; TP53, AR, SPOP and FOXA1, RB1, IDH
meta.rows[meta.rows$genes == 'AR', 'pathway'] = ''
meta.rows[meta.rows$genes == 'PTEN', 'pathway'] = ''
meta.rows[meta.rows$genes == 'TP53', 'pathway'] = ''
meta.rows[meta.rows$genes == 'RB1', 'pathway'] = ''
meta.rows[meta.rows$genes == 'IDH1', 'pathway'] = ''
meta.rows[meta.rows$genes == 'SPOP', 'pathway'] = ''
meta.rows[meta.rows$genes == 'FOXA1', 'pathway'] = ''

meta.rows = meta.rows[!is.na(meta.rows$pathway),, drop = F]

# DNA repair
DNA.repair = data.frame(genes = c('BRCA2',
               'BRCA1',
               'ATM',
               'FANCA',
               'CDK12',
               'MSH2',
               'MLH1'),
               pathway = 'DNA repair')
               
# PI3K pathway
PI3K = data.frame(genes = c('PIK3CA',
                            'PIK3R1',
                            'AKT1',
                            'AKT3'),
                  pathway = 'PI3K pathway')


# Wnt pathway
Wnt = data.frame(genes = c('APC',
                           'CTNNB1',
                           'RNF43'),
                 pathway = 'Wnt pathway')



# MAPK pathway 
MAPK = data.frame(genes = c('BRAF',
                            'HRAS',
                            'KRAS'),
                  pathway = 'MAPK pathway')


# chromatin remodeling
Chromatin = data.frame(genes = c('KMT2A',
                                 'KMT2C',
                                 'KMT2D',
                                 'KDM6A'),
                       pathway = 'Chromatin Remodeling')


meta.rows = rbind(meta.rows,
                  DNA.repair,
                  PI3K,
                  Wnt,
                  MAPK,
                  Chromatin)

row.names(meta.rows) = meta.rows$genes
meta.rows$pathway = factor(meta.rows$pathway, levels = c('', 'DNA repair', 'PI3K pathway', 'Wnt pathway', 
                                                         'MAPK pathway', 'Chromatin Remodeling'))
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

som.mat.print = som.mat.print[order(match(row.names(som.mat.print), row.names(meta.rows))), ]

## change global options (padding between top annotation and heatmap)
ht_opt$COLUMN_ANNO_PADDING = unit(2.5, "mm")

x.first = ComplexHeatmap::oncoPrint(som.mat.print,
                          row_names_gp = gpar(fontsize = 10),
                          row_names_side = 'left',
                          alter_fun = alter_fun, 
                          pct_side = 'right',
                          show_pct = T,
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
col = c("deep_deletion" = "#0C02FF", 
        "AMP" = "#FF0101",
        'Fusion' = '#8B02C9',
        "Missense_Mutation" = "#008004",
        'Truncating_Mutation' = '#000000',
        'Inframe_Mutation' = '#A78029',
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

