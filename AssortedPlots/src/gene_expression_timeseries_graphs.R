library(tidyverse)
library(ggforce)

setwd('~/')

#import data, format columns 
data <- read_csv('Thesis/RiceRNASeq/AssortedPlots/data/salmon_counts.csv')
data <- data %>% pivot_longer( !Name ,names_to = c('group','replicate'), values_to = 'count', names_sep='_') %>% 
  mutate(genotype=str_split(group, 'e', simplify = TRUE)[,1], time=as.integer(str_split(group, 'e', simplify = TRUE)[,2]))
data <- rename(data, Gene=Name)
data <- rename(data, Genotype=genotype)

#create id col to group by
data$id <- paste(data$Gene, data$group)

#group by id to merge replicates, average counts and get se
standard_error <- function(x) sd(x) / sqrt(length(x))
data <- data %>% group_by(id) %>% summarise(avg_count = mean(count), across(!c(replicate))) %>% summarize(se_count = standard_error(count), across()) %>% summarize(sd_count = sd(count), across())

#select only certain genes
genes <- c( "Os01t0100466-00", "Os01t0100200-01", "Os01t0100300-00", "Os01t0100650-00", "Os01t0100700-01", "Os01t0101175-00")

#color genes accordingly
colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#06d6a0','#118ab2', '#073b4c')

#set ymax to 20 counts higher than max
data <- data %>% mutate(orig_avg_count = avg_count, orig_se_count = se_count)
data <- data %>% mutate(avg_count = avg_count / 100, se_count = se_count / 100)
ymax <- round(max(filter(data, Gene %in% genes)$avg_count), 1)*1.15

#format data for excel output
subset <- data %>% filter(Gene %in% genes) %>% ungroup() %>%  select(c(time, Gene, avg_count, se_count, sd_count, Genotype)) %>% 
  unique() %>% pivot_wider(names_from=time, values_from=c(avg_count, se_count, sd_count)) %>% 
  select(Gene,  Genotype, avg_count_0, avg_count_4, avg_count_8, avg_count_12, avg_count_24, avg_count_36, 
         se_count_0, se_count_4, se_count_8, se_count_12, se_count_24, se_count_36,
         sd_count_0, sd_count_4, sd_count_8, sd_count_12, sd_count_24, sd_count_36) %>% 
  pivot_longer(cols=contains('_count_'), names_to=c('statistic', 'trash', 'time'), names_sep='_',values_to='count') %>%
  select(c(Gene,Genotype, statistic, time, count)) %>% pivot_wider(names_from=time, values_from=count) %>% 
  arrange(statistic, Gene, desc(Genotype)) %>% mutate(Gene_Genotype = paste(Gene, Genotype, sep='-'), .before = statistic) %>% 
  select(!c(Gene, Genotype)) %>% relocate(statistic, .after = last_col())


s1 <- subset %>% filter(statistic=='avg')
s2 <- subset %>% filter(statistic=='se') %>% select(!c(Gene_Genotype))
s3 <- subset %>% filter(statistic=='sd') %>% select(!c(Gene_Genotype))

write_delim(cbind(s1,s2,s3), 'Downloads/expression_timeseries_stats.csv', delim=',',quote='all')


#plot that graph
gg <- ggplot(filter(data, Gene %in% genes), mapping=aes(x=time, y=avg_count, color=Gene))+
  geom_line(mapping=aes(linetype=Genotype),alpha=0.7)+
  geom_point(alpha=0.7,shape='square')+
  geom_errorbar(mapping=aes(ymax = avg_count+se_count, ymin=avg_count-se_count), alpha=0.7, width=1.5)+
  scale_linetype_discrete(c('w','m'))+
  xlab('Time')+
  ylab('Count / 100')+
  scale_color_manual(values=colors)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,ymax), breaks=c(0, 0.01, 0.025 ,0.5,0.75,1.15))+
  scale_x_continuous(expand = c(0, 0),breaks=c(0,4,8,12,24,36))+
  theme_classic()+
  theme(axis.title.x = element_text(family='Arial', size=12),
        axis.title.y = element_text(family='Arial', size=12),
        axis.text.x = element_text(family='Arial', size=10),
        axis.text.y = element_text(family='Arial', size=10),
        axis.line = element_line(size=1))

#show the graph
gg
