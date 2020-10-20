## Racial comparisons:

## mutation occurrence distribution
alteration.occurrence = Somatic_alteration_matrix[, c(2:ncol(Somatic_alteration_matrix))]
row.names(alteration.occurrence) = Somatic_alteration_matrix$Sample.ID
gene_alterations = rowSums(alteration.occurrence != "")
gene_alterations = as.data.frame(gene_alterations)
gene_alterations$Sample.ID = row.names(gene_alterations)
gene_alterations = merge(gene_alterations, PRAD_meta.data[,c('Sample.ID', 'Race.Category', 'Sample.Type')],
                          by.x = 'Sample.ID', by.y = 'Sample.ID', all.x = T)

gene_alterations$Race.Category =  ifelse(is.na(gene_alterations$Race.Category), 'OTHER',
                                          ifelse(gene_alterations$Race.Category == 'PT REFUSED TO ANSWER', 'REFUSED TO ANSWER',
                                                 ifelse(gene_alterations$Race.Category == 'NO VALUE ENTERED', 'OTHER',
                                                        ifelse(gene_alterations$Race.Category == 'UNKNOWN', 'OTHER',
                                                               ifelse(gene_alterations$Race.Category == 'NATIVE HAWAIIAN OR PACIFIC ISL', 'OTHER',
                                                                      gene_alterations$Race.Category)))))

gene_alterations$grouping = ifelse(gene_alterations$gene_alterations == 1, 'group1',
                                   ifelse(gene_alterations$gene_alterations == 2, 'group2',
                                          ifelse(gene_alterations$gene_alterations > 2 & gene_alterations$gene_alterations < 6, 'group3',
                                                 ifelse(gene_alterations$gene_alterations > 5 & gene_alterations$gene_alterations < 11, 'group4',
                                                        ifelse(gene_alterations$gene_alterations > 10 & gene_alterations$gene_alterations < 21, 'group5', 'group6')))))


gene_alterations = gene_alterations[gene_alterations$Race.Category %in% c('ASIAN-FAR EAST/INDIAN SUBCONT', 
                                                                          'BLACK OR AFRICAN AMERICAN', 'WHITE'), ]

## make grouped barchart
library(ggplot2)
primary_samples = gene_alterations[gene_alterations$Sample.Type == 'Primary', ]

alter.matrix = as.data.frame(xtabs(~primary_samples$Race.Category + primary_samples$grouping))
alter.matrix$Freq = (alter.matrix$Freq / c(43, 129, 1165)) * 100

colnames(alter.matrix) = c('Race', 'group', 'freq')

primary.distribution = ggplot(alter.matrix, aes(fill = Race, x = group, y = freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  
  scale_fill_manual(values = c('#FF6237', '#FEFF9D', '#00A488'),
                    labels = c("Asian", "Black", "White")) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50), breaks = c(0, 10, 20, 30, 40, 50)) +
  
  theme(axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        panel.background = element_blank(),
        aspect.ratio = 0.5,
        legend.position = 'top',
        axis.text = element_text(size = 10, colour = 'black'),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  
  scale_x_discrete(labels=c("group1" = "1", "group2" = "2",
                            "group3" = "3-5", 'group4' = '6-10',
                            'group5' = '11-20', 'group6' = '>20')) +
  
  labs(y = 'Percentage of Patients', x = 'No. of mutations')
  

primary.distribution


## look at metastatic samples
metastatic_samples = gene_alterations[gene_alterations$Sample.Type == 'Metastasis', ]

alter.matrix = as.data.frame(xtabs(~metastatic_samples$Race.Category + metastatic_samples$grouping))
alter.matrix$Freq = (alter.matrix$Freq / c(45, 84, 860)) * 100

colnames(alter.matrix) = c('Race', 'group', 'freq')

metastatic.distribution = ggplot(alter.matrix, aes(fill = Race, x = group, y = freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  
  scale_fill_manual(values = c('#FF6237', '#FEFF9D', '#00A488'),
                    labels = c("Asian", "Black", "White")) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50), breaks = c(0, 10, 20, 30, 40, 50)) +
  
  theme(axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        panel.background = element_blank(),
        aspect.ratio = 0.5,
        legend.position = 'top',
        axis.text = element_text(size = 10, colour = 'black'),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  
  scale_x_discrete(labels=c("group1" = "1", "group2" = "2",
                            "group3" = "3-5", 'group4' = '6-10',
                            'group5' = '11-20', 'group6' = '>20')) +
  
  labs(y = 'Percentage of Patients', x = 'No. of mutations')


metastatic.distribution



## looking at specific genes:
gene.mutation = Somatic_alteration_matrix
genes = c('AR',
          'TP53',
          'FAT1',
          'TMPRSS2',
          'ATR',
          'ATM',
          'ERG',
          'DOT1L',
          'ZFHX3',
          'PTEN',
          'APC',
          'BRAF',
          'SPOP',
          'FOXA1',
          'RB1',
          'IDH1',
          'BRCA2',
          'BRCA1',
          'CDK12',
          'KMT2C',
          'KMT2D',
          'KDM6A',
          'KMT2A',
          'CTNNB1',
          'RNF43',
          'PIK3CA',
          'PIK3R1',
          'AKT1',
          'AKT3')

# some modifications for data display
gene.mutation = Somatic_alteration_matrix[, colnames(Somatic_alteration_matrix) %in% genes]
row.names(gene.mutation) = Somatic_alteration_matrix$Sample.ID
gene.mutation$Sample.ID = row.names(gene.mutation)
gene.mutation = merge(gene.mutation, PRAD_meta.data[,c('Sample.ID', 'Race.Category', 'Sample.Type')],
                      by.x = 'Sample.ID', by.y = 'Sample.ID', all.x = T)

gene.mutation = gene.mutation[gene.mutation$Sample.Type == 'Primary', ]
gene.mutation = gene.mutation[gene.mutation$Race.Category %in% c('ASIAN-FAR EAST/INDIAN SUBCONT', 
                                                                                          'BLACK OR AFRICAN AMERICAN', 'WHITE'), ]
gene.mutation$Sample.ID = NULL
gene.mutation$Sample.Type = NULL


all.out = data.frame()
for(i in unique(gene.mutation$Race.Category)){
  print(i)
  data.sub = gene.mutation[gene.mutation$Race.Category == i,, drop = F]
  data.process = as.data.frame(colSums(data.sub != ""))
  colnames(data.process) = 'freq'
  sum.group = data.process[nrow(data.process), ]
  data.process$freq = (data.process$freq / sum.group) * 100
  data.process$group = rep(i, nrow(data.process))
  data.process$gene = row.names(data.process)
  
  all.out = rbind(all.out, data.process)

}


all.out = all.out[!all.out$gene == 'Race.Category', ]


## make gene alteration plot
gene.primary.distribution = ggplot(all.out, aes(fill = group, x = gene, y = freq)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.7) +
  
  scale_fill_manual(name = '', values = c('#FF6237', '#FEFF9D', '#00A488'),
                    labels = c("Asian", "Black", "White")) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60), breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  
  theme(axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, colour = 'grey'),
        aspect.ratio = 1.75,
        legend.position = 'top',
        axis.text = element_text(size = 10, colour = 'black'),
        axis.title = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  
  
  labs(y = 'Patients with Mutations (%)', x = '') +
  
  coord_flip()


gene.primary.distribution


##
## look into metastatic samples:
gene.mutation = Somatic_alteration_matrix[, colnames(Somatic_alteration_matrix) %in% genes]
row.names(gene.mutation) = Somatic_alteration_matrix$Sample.ID
gene.mutation$Sample.ID = row.names(gene.mutation)
gene.mutation = merge(gene.mutation, PRAD_meta.data[,c('Sample.ID', 'Race.Category', 'Sample.Type')],
                      by.x = 'Sample.ID', by.y = 'Sample.ID', all.x = T)

gene.mutation = gene.mutation[gene.mutation$Sample.Type == 'Metastasis', ]
gene.mutation = gene.mutation[gene.mutation$Race.Category %in% c('ASIAN-FAR EAST/INDIAN SUBCONT', 
                                                                 'BLACK OR AFRICAN AMERICAN', 'WHITE'), ]
gene.mutation$Sample.ID = NULL
gene.mutation$Sample.Type = NULL


all.out = data.frame()
for(i in unique(gene.mutation$Race.Category)){
  print(i)
  data.sub = gene.mutation[gene.mutation$Race.Category == i,, drop = F]
  data.process = as.data.frame(colSums(data.sub != ""))
  colnames(data.process) = 'freq'
  sum.group = data.process[nrow(data.process), ]
  data.process$freq = (data.process$freq / sum.group) * 100
  data.process$group = rep(i, nrow(data.process))
  data.process$gene = row.names(data.process)
  
  all.out = rbind(all.out, data.process)
  
}


all.out = all.out[!all.out$gene == 'Race.Category', ]


## make gene alteration plot
gene.metastasis.distribution = ggplot(all.out, aes(fill = group, x = gene, y = freq)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.7) +
  
  scale_fill_manual(name = '', values = c('#FF6237', '#FEFF9D', '#00A488'),
                    labels = c("Asian", "Black", "White")) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60), breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  
  theme(axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(size = 0.3, colour = 'grey'),
        aspect.ratio = 1.75,
        legend.position = 'top',
        axis.text = element_text(size = 10, colour = 'black'),
        axis.title = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  
  
  labs(y = 'Patients with Mutations (%)', x = '') +
  
  coord_flip()


gene.metastasis.distribution




# further categories
# Actionable genes:
# Actionable = c('ABL1', 'EGFR', 'ERBB2', 'BRAF', 'BRCA1',
#                'BRAC2', 'FGFR2', 'FGFR3', 'KIT', 'NTRK1',
#                'NTRK2', 'NTRK3', 'PDGFRA', 'RET', 'ROS1', 'ALK', 'PIK3CA')
# 
# DNA repair = c('ERCC5', 'MRE11', 'TP53BP1', 'POLE', 'RAD21', 
#                'MSH2', 'MSH6', 'BRCA1', 'BRCA2', 'ATR', 'ATM')
