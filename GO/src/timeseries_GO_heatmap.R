#library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
library(extrafont)
library(tidyverse)
loadfonts(device = "postscript")

goterms <- read_delim('Thesis/RiceRNASeq/GO/data/timecourse.up.aggregated.slim.tsv', delim='\t')
#goterms %>% arrange(FDR) 
#fix label on MF 
s <- 'oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor'
ns <- 'oxidoreductase activity,\n acting on the CH-OH group of donors, NAD or NADP as acceptor'
goterms$label <- gsub(s,ns,goterms$label)
goterms$module <- factor(goterms$module, levels=c('36hr', '24hr','12hr','8hr','4hr','0hr'))
goterms <- filter(goterms, annotation=='MF')
top_go <- unique(arrange(goterms, -`fold-enrichment`)$Go_Term)[1:15]
ggplot(filter(goterms, Go_Term %in% top_go))+
  geom_tile(mapping = aes(x=module, y=paste(Go_Term, label), fill=`fold-enrichment`))+
  theme_classic(base_size = 8, base_family='Arial Black')+
  scale_fill_gradient(low='white',high='#dc3b3b',name='Upregulated\nfold-enrichment')+
  theme(legend.position = 'left')+
  labs(x='', y='')
ggsave('Thesis/RiceRNASeq/GO/data/DEG/figs/GO_DEG_UP_MF.png', dpi=1300, width=7.00, height=6.00, units='in')

goterms <- read_delim('Thesis/RiceRNASeq/GO/data/timecourse.down.aggregated.slim.tsv', delim='\t')
#goterms %>% arrange(FDR) 
#fix label on MF 
s <- 'oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen'
ns <- 'oxidoreductase activity, acting on single donors with incorporation of molecular oxygen,\n incorporation of two atoms of oxygen'
goterms$label <- gsub(s,ns,goterms$label)

goterms$module <- factor(goterms$module, levels=c('0hr', '4hr', '8hr', '12hr', '24hr', '36hr'))
goterms <- filter(goterms, annotation=='MF')
top_go <- unique(arrange(goterms, -`fold-enrichment`)$Go_Term)[1:15]
ggplot(filter(goterms, Go_Term %in% top_go))+
  geom_tile(mapping = aes(x=module, y=paste(label,Go_Term), fill=`fold-enrichment`))+
  scale_fill_gradient(low='white',high='#186aa1',name='Downregulated\nfold-enrichment')+
  theme_classic(base_size = 8, base_family='Arial Black')+
  theme(legend.position='right')+
  scale_y_discrete(position="right")+
  labs(x='', y='')
ggsave('Thesis/RiceRNASeq/GO/data/DEG/figs/GO_DEG_DOWN_MF.png', dpi=1300, width=7.50, height=6.00, units='in')
