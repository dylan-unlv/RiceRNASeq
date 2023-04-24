library(tidyverse)
library(ggalluvial)

#this script creates a flowchart showing which modules genes are grouped between two runs of WGCNA

filtered <- read_delim('Downloads/gene_modules/gene_modules.txt', delim='\t')
nofilter <- read_delim('Downloads/gene_modules/gene_modules.nofilter.txt', delim='\t')

merged <- merge(filtered, nofilter, by='gene_id')
merged <- merged %>%  select(gene_id, fcolors, nfcolors)

links <- merged %>% group_by(fcolors, nfcolors) %>% count()
links$fcolors <- sprintf('%s.f', links$fcolors)
links$nfcolors <- sprintf('%s.nf', links$nfcolors)
links$fcolors <- as.factor(links$fcolors)
links$nfcolors <- as.factor(links$nfcolors)

ggplot(links, aes(axis1 = fcolors, axis2 = nfcolors, y=n))+
  geom_alluvium()+
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("fcolors", "nfcolors"),
                   expand = c(.1, .1)) +
  labs(title = "WGCNA filtered vs non filtered",
       subtitle = "genes grouped by module",
       y = "Frequency") +
  scale_fill_viridis_b()+
  theme_minimal()
  
  