## ICGC Data Portal data analysis

sample.sheet = read.csv('~/Documents/MSKCC/00_Data/ICGC_PCAWG_data/pcawg_sample_sheet.txt', sep = '\t')

## prostate codes:

# n = 1,419 donors; 500 TCGA PRAD
# 'EOPC-DE'
EOPC = sample.sheet[grep('EOPC', sample.sheet$donor_unique_id), ]
EOPC = EOPC[grep('tumour', EOPC$dcc_specimen_type), ]
EPOC = EOPC[EOPC$donor_wgs_exclusion_white_gray == 'Whitelist', ]

# 'PRAD-CA'
PRAD.CA = sample.sheet[grep('PRAD-CA', sample.sheet$donor_unique_id), ]
PRAD.CA = PRAD.CA[grep('tumour', PRAD.CA$dcc_specimen_type), ]
PRAD.CA = PRAD.CA[PRAD.CA$donor_wgs_exclusion_white_gray == 'Whitelist', ]

# PRAD-UK
PRAD.UK = sample.sheet[grep('PRAD-UK', sample.sheet$donor_unique_id), ]
PRAD.UK = PRAD.UK[grep('tumour', PRAD.UK$dcc_specimen_type), ]
PRAD.UK = PRAD.UK[PRAD.UK$donor_wgs_exclusion_white_gray == 'Whitelist', ]
PRAD.UK = PRAD.UK[grep('Primary', PRAD.UK$dcc_specimen_type), ]

# combine all non TCGA samples
Prad_combined = rbind(EOPC, PRAD.CA, PRAD.UK)

# maf annotation file:
maf = vroom::vroom(file = '~/Documents/MSKCC/00_Data/ICGC_PCAWG_data/final_consensus_passonly.snv_mnv_indel.icgc.public.maf', delim = '\t')
maf.subset = maf[maf$Tumor_Sample_Barcode %in% Prad_combined$aliquot_id, ]
# maf.subset = maf.subset[maf.subset$Hugo_Symbol %in% genes.remaining, ]



## mutation matrix
mutation.matrix = function(mut){
  gene.panel = read.csv('~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/impact468_gene_panel.txt', sep = '\t')
  sample = data.frame(genes = 'Sample.ID')
  gene.panel = rbind(sample, gene.panel)
  
  # data somatic mutations
  data.mut = as.data.frame(mut)
  
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
  
  row.names(alteration_matrix) = alteration_matrix$Sample.ID
  alteration_matrix$Sample.ID = NULL
  alteration_matrix[is.na(alteration_matrix)] = ''
  
  return(list(alteration_matrix = alteration_matrix,
              gene.panel = gene.panel))
}


x = mutation.matrix(mut = maf.subset)


#~~~~~~~~~``

onco.matrix = function(mut){

  mut.matrix.out = mutation.matrix(mut = mut)
  alteration_matrix = mut.matrix.out$alteration_matrix
  gene.panel = mut.matrix.out$gene.panel
  gene.panel = as.character(gene.panel$genes)
  gene.panel = gene.panel[-1]
  
  # load genomic coordiantes
  load('~/Documents/MSKCC/00_Data/gene_coordinates.rda')
  genes_hg19 = as.data.frame(genes_hg19)
  gene.panel = intersect(gene.panel, genes_hg19$gene)
  
  ## how does SCNA data look like:
  message('preparing SCNA data')
  icg.cna.part = list.files('~/Documents/MSKCC/00_Data/ICGC_PCAWG_data/SCNA_PCAWG/', full.names = F)
  icg.cna = gsub('\\..*', '', icg.cna.part)
  icg.cna = intersect(mut$Tumor_Sample_Barcode, icg.cna)
  icg.cna.full = list.files('~/Documents/MSKCC/00_Data/ICGC_PCAWG_data/SCNA_PCAWG/', full.names = T)
  
  
  SCNA.paths.out = c()
  for(i in 1:length(icg.cna)){
    print(i)
    SCNA.paths = grep(icg.cna[i], icg.cna.full, value = T)
    SCNA.paths.out = c(SCNA.paths.out, SCNA.paths)
  }
  
  
  # make transformation
  out = data.frame()
  out = data.frame(matrix(ncol = length(gene.panel), nrow = length(SCNA.paths.out)))
  
  colnames(out) = gene.panel
  row.names(out) = SCNA.paths.out
  
  
  for(i in 1:length(SCNA.paths.out)){
    
    data.cna = read.csv(file = SCNA.paths.out[i], sep = '\t')
    print(i)
    
    for(j in 1:length(gene.panel)){
      gene.selected = genes_hg19[genes_hg19$gene == gene.panel[j], ]
      Segmentation_gene = data.cna[data.cna$chromosome == gene.selected$chrom, ]
      Segmentation_gene = Segmentation_gene[which(Segmentation_gene$start <= gene.selected$start &
                                                    Segmentation_gene$end >= gene.selected$end |
                                                    Segmentation_gene$start >= gene.selected$start &
                                                    Segmentation_gene$start < gene.selected$end &
                                                    Segmentation_gene$end >= gene.selected$end |
                                                    Segmentation_gene$start <= gene.selected$start &
                                                    Segmentation_gene$end <= gene.selected$end &
                                                    Segmentation_gene$end > gene.selected$start), ]
      if(nrow(Segmentation_gene) == 1){
        status = ifelse(Segmentation_gene$total_cn == 0, 'deep_deletion',
                        ifelse(Segmentation_gene$total_cn > 3, 'amp', 'neutral'))
        status = unique(status)
        gene = j
        
      } else {
        status = 'NA'
      }
      out[row.names(out) == SCNA.paths.out[i], gene.panel[j]] = status
      
    }
  }
  
  # modify matrix
  out[out == 'neutral'] = ''
  row.names(out) = sub('.*//', replacement = '', row.names(out))
  row.names(out) = gsub('\\..*', '', row.names(out))
  
  x = alteration_matrix[, colnames(alteration_matrix) %in% colnames(out)]
  
  # merge
  L <- Map(function(...) paste(..., sep = ";"), x, out)
  k = replace(x, TRUE, L)
  
  k[k == ';NA'] = ''
  k[k == ';'] = ''
  
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
                                                                                                  ifelse(grepl('Silent', k[, i]), sub(pattern = 'Silent', replacement = 'VUS', k[,i]),
                                                                                                         ifelse(grepl('IGR', k[, i]), sub(pattern = 'IGR', replacement = 'VUS', k[,i]),
                                                                                                                ifelse(grepl('3\'UTR', k[, i]), sub(pattern = '3\'UTR', replacement = 'VUS', k[,i]),
                                                                                                                       ifelse(grepl('5\'Flank', k[, i]), sub(pattern = '5\'Flank', replacement = 'VUS', k[,i]), k[,i]))))))))))))))))
  }
  
  
  
  return(k)

}
  
x = onco.matrix(mut = maf.subset)









