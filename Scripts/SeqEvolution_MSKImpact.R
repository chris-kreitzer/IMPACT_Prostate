## sequencing evolution at MSKCC

seq = read_excel(path = '~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/Sequencing_Dates_MSKIMPACT.xlsx')
seq = seq[!is.na(seq$Year), ]
months = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

for(i in 1:nrow(seq)){
  seq$quartal[i] = ifelse(seq$Month[i] %in% months[1:6], 'first', 'second')
}

seq$SAMPLE_ID = NULL


# patient.ID
clini = read.table(file = '~/Documents/MSKCC/00_Data/IMPACT_DATA_2020.08/mskimpact.clinical_data_092020.tsv', 
                      sep = "\t", 
                      header = T, 
                      comment.char = "#",
                      na.strings = ".", 
                      stringsAsFactors = FALSE,
                      quote = "", fill = FALSE)

clini = clini[!duplicated(clini$Patient.ID), ]

meta.impact = merge(seq, clini[, c('Patient.ID', 'Cancer.Type', 'Race.Category')],
                    by.x = 'PATIENT_ID',
                    by.y = 'Patient.ID', 
                    all.x = T)


library(ggplot2)
all.samples = data.frame(with(meta.impact, table(Year)))
white = data.frame(with(subset(meta.impact, Race.Category == 'WHITE'), table(Year)))
asian = data.frame(with(subset(meta.impact, Race.Category == 'ASIAN-FAR EAST/INDIAN SUBCONT'), table(Year)))
black = data.frame(with(subset(meta.impact, Race.Category == 'BLACK OR AFRICAN AMERICAN'), table(Year)))
refused = data.frame(with(subset(meta.impact, Race.Category == 'PT REFUSED TO ANSWER'), table(Year)))
other = subset(meta.impact, Race.Category %in% c('NO VALUE ENTERED',
                                                 'UNKNOWN',
                                                 'NATIVE AMERICAN-AM IND/ALASKA',
                                                 'NATIVE HAWAIIAN OR PACIFIC ISL',
                                                 'OTHER'))
other = data.frame(with(other, table(Year)))

# make plot
Seq.Impact = ggplot() +
  geom_line(data = all.samples, aes(x = Year, y = cumsum(Freq) / 1000),
            group = 0,
            color = '#2d4753',
            size = 2) + 
  geom_point(data = all.samples, aes(x = Year, y = cumsum(Freq) / 1000)) +
  
  geom_line(data = white, aes(x = Year, y = cumsum(Freq) / 1000),
            group = 0,
            color = '#5a9f8e',
            size = 1) + 
  geom_point(data = white, aes(x = Year, y = cumsum(Freq) / 1000)) +
  
  geom_line(data = asian, aes(x = Year, y = cumsum(Freq) / 1000),
            group = 0,
            color = '#f5de7a',
            size = 1) + 
  geom_point(data = asian, aes(x = Year, y = cumsum(Freq) / 1000)) +
  
  geom_line(data = black, aes(x = Year, y = cumsum(Freq) / 1000),
            group = 0,
            color = '#e48654',
            size = 1) + 
  geom_point(data = black, aes(x = Year, y = cumsum(Freq) / 1000)) +
  
  geom_line(data = refused, aes(x = Year, y = cumsum(Freq) / 1000),
            group = 0,
            color = '#eebebd',
            size = 1) + 
  geom_point(data = refused, aes(x = Year, y = cumsum(Freq) / 1000)) +
  
  geom_line(data = other, aes(x = Year, y = cumsum(Freq) / 1000),
            group = 0,
            color = 'grey55',
            size = 1) + 
  geom_point(data = other, aes(x = Year, y = cumsum(Freq) / 1000)) + 
  
  scale_y_continuous(expand = c(0.05, 0.05), 
                     labels = c('0K', '10K', '20K', '30K', '40K', '50K')) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  
  labs(x = 'Year', 
       y = "Cumulative sum of patients",
       title = 'The evolution of IMPACT sequencing',
       subtitle = "Patients enrolled in MSK cancer care program where at least one tumor specimen\nwas already sequenced. Colors depict different ethnicities",
       caption = "considered are patients") + 
  
  theme(aspect.ratio = 2) +
  theme_ipsum(base_family = "Arial Narrow",
              base_size = 14,
              plot_title_family = "Arial Narrow",
              plot_title_size = 18,
              plot_title_face = "bold",
              plot_title_margin = 10,
              subtitle_family = "Arial Narrow",
              subtitle_size = 12,
              subtitle_face = "plain",
              subtitle_margin = 15,
              strip_text_family = "Arial Narrow",
              strip_text_size = 12,
              strip_text_face = "plain",
              caption_family = "Arial Narrow",
              caption_size = 9,
              caption_face = "italic",
              caption_margin = 10,
              axis_text_size = 12,
              axis_title_family = "Arial Narrow",
              axis_title_size = 12,
              axis_title_face = "plain",
              axis_title_just = "rt",
              plot_margin = margin(30, 30, 30, 30),
              grid_col = "#cccccc",
              grid = TRUE,
              axis_col = "#cccccc",
              axis = T,
              ticks = F)
  




##
plot_list = list()

for(i in unique(meta.impact$Cancer.Type)){
  
  data = meta.impact[meta.impact$Cancer.Type == i, ]
  
  all.samples = data.frame(with(data, table(Year)))
  white = data.frame(with(subset(data, Race.Category == 'WHITE'), table(Year)))
  asian = data.frame(with(subset(data, Race.Category == 'ASIAN-FAR EAST/INDIAN SUBCONT'), table(Year)))
  black = data.frame(with(subset(data, Race.Category == 'BLACK OR AFRICAN AMERICAN'), table(Year)))
  refused = data.frame(with(subset(data, Race.Category == 'PT REFUSED TO ANSWER'), table(Year)))
  other = subset(data, Race.Category %in% c('NO VALUE ENTERED',
                                                   'UNKNOWN',
                                                   'NATIVE AMERICAN-AM IND/ALASKA',
                                                   'NATIVE HAWAIIAN OR PACIFIC ISL',
                                                   'OTHER'))
  other = data.frame(with(other, table(Year)))
  
  # define y cutoff
  y.all = round(sum(all.samples$Freq)) + 100

  # make plot
  plot = ggplot() +
    # all
    geom_line(data = all.samples, aes(x = Year, y = cumsum(Freq)),
              group = 0,
              color = '#2d4753',
              size = 2) + 
    geom_point(data = all.samples, aes(x = Year, y = cumsum(Freq))) +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0, y.all), 
                       breaks = c(0, round(y.all * 0.25), round(y.all * 0.5), round(y.all * 0.75), y.all)) +
    
    # white
    geom_line(data = white, aes(x = Year, y = cumsum(Freq)),
              group = 0,
              color = '#5a9f8e',
              size = 1) + 
    geom_point(data = white, aes(x = Year, y = cumsum(Freq))) +

    # asian
    geom_line(data = asian, aes(x = Year, y = cumsum(Freq)),
              group = 0,
              color = '#f5de7a',
              size = 1) + 
    geom_point(data = asian, aes(x = Year, y = cumsum(Freq))) +

    # black
    geom_line(data = black, aes(x = Year, y = cumsum(Freq)),
              group = 0,
              color = '#e48654',
              size = 1) + 
    geom_point(data = black, aes(x = Year, y = cumsum(Freq))) +

    # refused
    geom_line(data = refused, aes(x = Year, y = cumsum(Freq)),
              group = 0,
              color = '#eebebd',
              size = 1) + 
    geom_point(data = refused, aes(x = Year, y = cumsum(Freq))) +

    # other
    geom_line(data = other, aes(x = Year, y = cumsum(Freq)),
              group = 0,
              color = 'grey55',
              size = 1) + 
    geom_point(data = other, aes(x = Year, y = cumsum(Freq))) + 

    labs(title = i,
         x = NULL,
         y = NULL) +

    theme(aspect.ratio = 1) +
    theme_ipsum(base_family = "Arial Narrow",
                base_size = 14,
                plot_title_family = "Arial Narrow",
                plot_title_size = 18,
                plot_title_face = "bold",
                plot_title_margin = 10,
                subtitle_family = "Arial Narrow",
                subtitle_size = 12,
                subtitle_face = "plain",
                subtitle_margin = 15,
                strip_text_family = "Arial Narrow",
                strip_text_size = 12,
                strip_text_face = "plain",
                caption_family = "Arial Narrow",
                caption_size = 9,
                caption_face = "italic",
                caption_margin = 10,
                axis_text_size = 12,
                axis_title_family = "Arial Narrow",
                axis_title_size = 12,
                axis_title_face = "plain",
                axis_title_just = "rt",
                plot_margin = margin(30, 30, 30, 30),
                grid_col = "#cccccc",
                grid = F,
                axis_col = "#cccccc",
                axis = T,
                ticks = F) +
    theme(axis.text.y = element_blank()) +
    scale_x_discrete(expand = c(0.01, 0.01), drop = FALSE)
    
  
  plot_list[[i]] = plot
  
  rm(data)
  rm(plot)
  
}

ggarrange(plot_list[[1]], 
          plot_list[[2]],
          plot_list[[3]],
          plot_list[[4]],
          plot_list[[5]],
          plot_list[[6]],
          plot_list[[8]],
          plot_list[[9]],
          plot_list[[10]],
          plot_list[[13]],
          plot_list[[14]],
          plot_list[[15]],
          plot_list[[16]],
          plot_list[[17]],
          plot_list[[18]],
          plot_list[[22]],
          nrow = 4)












library(hrbrthemes)

ggplot() +
  geom_line(data = all.samples, aes(x = Year, y = cumsum(Freq), color = group),
            group = 0,
            color = '#2d4753',
            size = 2) + 
  geom_point(data = all.samples, aes(x = Year, y = cumsum(Freq), color = group)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y), breaks = c(0, round(y * 0.25), round(y * 0.5), round(y * 0.75), y)) +
  theme_ipsum(base_family = "Arial Narrow",
              base_size = 14,
              plot_title_family = "Arial Narrow",
              plot_title_size = 18,
              plot_title_face = "bold",
              plot_title_margin = 10,
              subtitle_family = "Arial Narrow",
              subtitle_size = 12,
              subtitle_face = "plain",
              subtitle_margin = 15,
              strip_text_family = "Arial Narrow",
              strip_text_size = 12,
              strip_text_face = "plain",
              caption_family = "Arial Narrow",
              caption_size = 9,
              caption_face = "italic",
              caption_margin = 10,
              axis_text_size = 12,
              axis_title_family = "Arial Narrow",
              axis_title_size = 12,
              axis_title_face = "plain",
              axis_title_just = "rt",
              plot_margin = margin(30, 30, 30, 30),
              grid_col = "#cccccc",
              grid = TRUE,
              axis_col = "#cccccc",
              axis = T,
              ticks = F) +
  theme(legend.position="bottom")
  





