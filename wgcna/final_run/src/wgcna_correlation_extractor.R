library(tidyverse)
library(magrittr)
library(WGCNA)
library(DESeq2)
library(genefilter)

#set wd to find precalculated TOM 
setwd('Thesis/RiceRNASeq/wgcna/final_run')

#import data
data <- readr::read_delim('raw_data/fC_all_MSU_v2.txt', delim='\t')

#clean data, prepare for input to wgcna
data <- select(data, -c('Chr','Start','End','Strand','Length'))

#only keep replicates in meta file (same as their DESeq analysis)
replicates <- readr::read_delim('raw_data/Meta_v2.txt', delim='\t')
keepers <- as.list(replicates[1])[[1]]
data <- select(data, c(Geneid, all_of(keepers)))

#remove file extension
for ( col in 1:ncol(data)){
  colnames(data)[col] <-  sub(".sam", "", colnames(data)[col])
}

#########################################
### normalize count matrix with DESeq ###
#########################################

de_input <- as.matrix(data[, -1])
rownames(de_input) <- data$Geneid 

#split by condition (each wt/mutant + time combination is its own condition)
meta_df <- data.frame( Sample = names(data[-1])) %>%
  mutate(
    Type = gsub("_.*","", Sample)
  )

dds <- DESeqDataSetFromMatrix(round(de_input),
                              meta_df,
                              design = ~Type)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
#q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
#expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
expr_normalized <- wpn_vsd[ rv_wpn > q75_wpn, ]


expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)


#########################################
###       prepare data for WGCNA      ###
#########################################

#multithread
allowWGCNAThreads(8)

#transform data
input_mat = t(expr_normalized)

#initialize WGCNA
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;



#########################################
###             Run WGCNA             ###
#########################################


#choose a power near the upper limit of the sigmoidal curve in the scale independence plot
#I've chosen 18, but we can experiment with other powers
picked_power = 18

#Run WGCNA
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Load archived TOM if exists
                          saveTOMFileBase = "ER",
                          loadTOM=TRUE,
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

#plot results
mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


#once you have network, get TOM similarity for a particular gene
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors))
rownames(module_df) <- module_df$gene_id

target_gene <- 'LOC_Os02g08440'

#find the module which this gene lies in
module_of_interest <- module_df %>% filter(gene_id==target_gene) %>% pull(colors) %>% as.character()

#extract genes from module
genes_of_interest = module_df %>%
  filter(colors == module_of_interest)

#subset the expression data and rerun TOM to get network similarity scores within module
expr_of_interest <- expr_normalized[as.character(genes_of_interest$gene_id),]
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

#extract similarity scores as edges of a network (can import to cytoscape)
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, TOM_similarity = value) %>%
  unique() %>%
  subset(!(gene1==gene2))

#write_delim(edge_list,
#            file = "data/correlated_WRKY71_edgelist.tsv",
#            delim = "\t")

#filter edge list by gene of interest, it's symmetrical so only need one side
#gene_edges <- edge_list %>% filter((gene1==target_gene | gene2==target_gene))
gene_edges <- edge_list %>% filter(gene1==target_gene)

#pick a percentile since TOM correlation doesn't exactly map to pearson correlation
threshold <- 0.995 #percentile from 0 - 1
quant <- quantile( gene_edges$TOM_similarity, threshold)
filtered_edges <- gene_edges %>% filter(TOM_similarity>quant)

genes <- c(target_gene, filtered_edges$gene2)
print('Genes extracted from TOM:')
print(genes)
#now load in expression data to a graph with the labeled target gene
gdata <- expr_normalized_df %>% 
  mutate(group=str_split(name, '_', simplify = TRUE)[,1], replicate=str_split(name, '_', simplify = TRUE)[,2]) %>% 
  mutate(Genotype=str_split(group, 'e', simplify = TRUE)[,1], time=as.integer(str_split(group, 'e', simplify = TRUE)[,2])) %>% 
  filter(Gene_id %in% genes)

#id col needs to be gene and group specific to merge replicates
gdata$id <- paste0(gdata$Gene_id, gdata$group) 

#group by id to merge replicates, average counts and get se
standard_error <- function(x) sd(x) / sqrt(length(x))
gdata <- gdata %>% group_by(id) %>% summarise(avg_count = mean(value), across(!c(replicate))) %>% summarize(se_count = standard_error(value), across())

#set ymax to 20 counts higher than max
ymax <- round(max(filter(gdata, Gene_id %in% genes)$avg_count), 1)+5
gdata$GeneType <- 'Correlated Gene'
gdata[gdata$Gene_id==target_gene,]$GeneType <- target_gene

#label the facets
flabels <- c('Mutant','Wild Type')
names(flabels) <- c('m','w')
gdata$Genotype <- factor(gdata$Genotype, levels=c('w','m'))

#plot that graph
ggplot(gdata, mapping=aes(x=time, y=avg_count, group=interaction(Gene_id, Genotype), color=GeneType))+
  geom_line(mapping=aes(linetype=Genotype),alpha=0.6)+
  geom_line(data=filter(gdata, Gene_id==target_gene), mapping=aes(linetype=Genotype),alpha=0.7, color='red')+
  #geom_errorbar(mapping=aes(ymax = avg_count+se_count, ymin=avg_count-se_count), alpha=0.7, width=1.5)+
  xlab('Time')+
  ylab('Normalized Count')+
  scale_color_manual(values=c('grey', 'red'))+
  scale_linetype(guide='none')+
  scale_y_continuous(expand = c(0, 0), limits=c(1e-5,ymax))+
  scale_x_continuous(expand = c(0, 0),breaks=c(0,4,8,12,24,36), limits=c(0,37))+
  facet_wrap(~Genotype, labeller = labeller(Genotype=flabels))+
  theme_classic()+
  theme(axis.title.x = element_text(family='Arial', size=12),
        axis.title.y = element_text(family='Arial', size=12),
        axis.text.x = element_text(family='Arial', size=10),
        axis.text.y = element_text(family='Arial', size=10),
        axis.line = element_line(size=1))

###
# export raw data for santi
###

#add TOM similarity score
gdata$TOM_Similarity = 0
for (i in unique(gdata$Gene_id)){
  if (i==target_gene){
    gdata[gdata$Gene_id==i,'TOM_Similarity'] = 1
  }
  else{
    gdata[gdata$Gene_id==i,'TOM_Similarity'] = gene_edges[gene_edges$gene2==i,]$TOM_similarity
  }
}

#export data

write_delim(gdata, paste0('data/',target_gene,'_correlated_genes.csv'),delim=',')

########################################################
###     TOM similarity for anti correlated gene      ###
########################################################

#connect modules (colors) to treatment groups
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order),
    group=paste0(str_split(treatment, '_', simplify = TRUE)[,1], '_',name),
    replicate=str_split(treatment, '_', simplify = TRUE)[,2],
    treatment=str_split(treatment, '_', simplify = TRUE)[,1],
    genotype=str_split(treatment, 'e', simplify = TRUE)[,1],
    time=str_split(treatment, 'e', simplify = TRUE)[,2],
  )

#create eigenspace
edf <- mME %>% group_by(group) %>% mutate(Mean=mean(value, na.rm=TRUE)) %>% select(-c(replicate, value)) %>% unique() %>% 
  ungroup() %>% select(-c(group)) %>% pivot_wider(names_from=name, values_from=Mean)

#find module that's the farthest away in eigenspace by calculating distance between eigenvectors
colors<-unique(as.character(module_df$colors))
colors <- colors[colors != module_of_interest]
max_dist <- 0
opposite_module <- ''
for (i in colors){
  tmp_dist <- sum(abs(edf[,module_of_interest] - edf[,i]))
  if (tmp_dist > max_dist){
    max_dist <- tmp_dist
    opposite_module<-i
  }
  }

#run TOM similarity with target gene and opposite module
genes_of_interest = module_df %>%
  filter(colors == opposite_module)
expr_of_interest <- expr_normalized[c(as.character(genes_of_interest$gene_id), target_gene),]
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

#filter, graph alongside original expression data
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, TOM_similarity = value) %>%
  unique() %>%
  subset(!(gene1==gene2))

#write_delim(edge_list,
#            file = "data/anticorrelated_WRKY71_edgelist.tsv",
#            delim = "\t")

#filter edge list by gene of interest, it's symmetrical so only need one side
#gene_edges <- edge_list %>% filter((gene1==target_gene | gene2==target_gene))
gene_edges <- edge_list %>% filter(gene1==target_gene)

#pick a percentile since TOM correlation doesn't exactly map to pearson correlation
threshold <- 0.995 #percentile from 0 - 1
quant <- quantile( gene_edges$TOM_similarity, threshold)
filtered_edges <- gene_edges %>% filter(TOM_similarity>quant)

genes <- c(target_gene, filtered_edges$gene2)
print('Genes extracted from TOM:')
print(genes[genes!=target_gene])
#now load in expression data to a graph with the labeled target gene
gdata <- expr_normalized_df %>% 
  mutate(group=str_split(name, '_', simplify = TRUE)[,1], replicate=str_split(name, '_', simplify = TRUE)[,2]) %>% 
  mutate(Genotype=str_split(group, 'e', simplify = TRUE)[,1], time=as.integer(str_split(group, 'e', simplify = TRUE)[,2])) %>% 
  filter(Gene_id %in% genes)

#id col needs to be gene and group specific to merge replicates
gdata$id <- paste0(gdata$Gene_id, gdata$group) 
gdata <- gdata %>% group_by(id) %>% summarise(avg_count = mean(value), across(!c(replicate))) %>% summarize(se_count = standard_error(value), across())

#set ymax to 20 counts higher than max
ymax <- round(max(filter(gdata, Gene_id %in% genes)$avg_count), 1)+5
gdata$GeneType <- 'Anti-correlated Gene'
gdata[gdata$Gene_id==target_gene,]$GeneType <- target_gene

#label the facets
flabels <- c('Mutant','Wild Type')
names(flabels) <- c('m','w')
gdata$Genotype <- factor(gdata$Genotype, levels=c('w','m'))

#plot that graph
ggplot(gdata, mapping=aes(x=time, y=avg_count, group=interaction(Gene_id, Genotype), color=GeneType))+
  geom_line(mapping=aes(linetype=Genotype),alpha=0.6)+
  geom_line(data=filter(gdata, Gene_id==target_gene), mapping=aes(linetype=Genotype),alpha=0.7, color='red')+
  #geom_errorbar(mapping=aes(ymax = avg_count+se_count, ymin=avg_count-se_count), alpha=0.7, width=1.5)+
  xlab('Time')+
  ylab('Normalized Count')+
  scale_color_manual(values=c('grey', 'red'))+
  scale_linetype(guide='none')+
  scale_y_continuous(expand = c(0, 0), limits=c(1e-5,ymax))+
  scale_x_continuous(expand = c(0, 0),breaks=c(0,4,8,12,24,36), limits=c(0,37))+
  facet_wrap(~Genotype, labeller = labeller(Genotype=flabels))+
  theme_classic()+
  theme(axis.title.x = element_text(family='Arial', size=12),
        axis.title.y = element_text(family='Arial', size=12),
        axis.text.x = element_text(family='Arial', size=10),
        axis.text.y = element_text(family='Arial', size=10),
        axis.line = element_line(size=1))

###
# export raw data for santi
###

#add TOM similarity score
gdata$TOM_Similarity = 0
for (i in unique(gdata$Gene_id)){
  if (i==target_gene){
    gdata[gdata$Gene_id==i,'TOM_Similarity'] = 1
  }
  else{
    gdata[gdata$Gene_id==i,'TOM_Similarity'] = gene_edges[gene_edges$gene2==i,]$TOM_similarity
  }
}

#export data

write_delim(gdata, paste0('data/',target_gene,'_anticorrelated_genes.csv'),delim=',')

########################################################
###          heatmap to choose module                ###
########################################################
mME$id <- paste0(mME$treatment,'_', mME$replicate)
level_order <- c('me0_1','me0_2','me0_3','me4_1','me4_2','me4_3','me8_1','me8_2','me8_3','me12_1','me12_2','me12_3','me24_1','me24_2','me24_4r','me36_1','me36_2r','me36_3r',
                 'we0_1','we0_2','we0_3','we4_1','we4_2','we4_3','we8_1','we8_2','we8_3','we12_1','we12_2','we12_3','we24_1','we24_2','we24_3','we36_1r','we36_2r','we36_3')
mod_treat_heatmap <- mME %>% ggplot(., aes(x=factor(id, level = level_order), y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", x='', fill="Corr")

#tan and purple modules show correlation with mutant genotype
#this means they may contain genes downregulated by WRKY71
#salmon module shows anti correlation with mutant genotype
#this means it may contain genes upregulated by WRKY71

########################################################
###       TOM similarity for specific module         ###
########################################################

#run TOM similarity with target gene and chosen module
module_of_interest <- 'purple'
genes_of_interest = module_df %>%
  filter(colors == module_of_interest)
expr_of_interest <- expr_normalized[c(as.character(genes_of_interest$gene_id), target_gene),]
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

#filter, graph alongside original expression data
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, TOM_similarity = value) %>%
  unique() %>%
  subset(!(gene1==gene2))

#write_delim(edge_list,
#            file = "data/tan_WRKY71_edgelist.tsv",
#            delim = "\t")

#filter edge list by gene of interest, it's symmetrical so only need one side
#gene_edges <- edge_list %>% filter((gene1==target_gene | gene2==target_gene))
gene_edges <- edge_list %>% filter(gene1==target_gene)

#pick a percentile since TOM correlation doesn't exactly map to pearson correlation
threshold <- 0.85 #percentile from 0 - 1
quant <- quantile( gene_edges$TOM_similarity, threshold)
filtered_edges <- gene_edges %>% filter(TOM_similarity>quant)

genes <- c(target_gene, filtered_edges$gene2)
print('Genes extracted from TOM:')
print(genes[genes!=target_gene])
#now load in expression data to a graph with the labeled target gene
gdata <- expr_normalized_df %>% 
  mutate(group=str_split(name, '_', simplify = TRUE)[,1], replicate=str_split(name, '_', simplify = TRUE)[,2]) %>% 
  mutate(Genotype=str_split(group, 'e', simplify = TRUE)[,1], time=as.integer(str_split(group, 'e', simplify = TRUE)[,2])) %>% 
  filter(Gene_id %in% genes)

#id col needs to be gene and group specific to merge replicates
gdata$id <- paste0(gdata$Gene_id, gdata$group) 
gdata <- gdata %>% group_by(id) %>% summarise(avg_count = mean(value), across(!c(replicate))) %>% summarize(se_count = standard_error(value), across())

#set ymax to 20 counts higher than max
ymax <- round(max(filter(gdata, Gene_id %in% genes)$avg_count), 1)+5
gdata$GeneType <- str_to_title(paste0(module_of_interest,'-Correlated Gene'))
gdata[gdata$Gene_id==target_gene,]$GeneType <- target_gene

#label the facets
flabels <- c('Mutant','Wild Type')
names(flabels) <- c('m','w')
gdata$Genotype <- factor(gdata$Genotype, levels=c('w','m'))6

#plot that graph
ggplot(gdata, mapping=aes(x=time, y=avg_count, group=interaction(Gene_id, Genotype), color=GeneType))+
  geom_line(mapping=aes(linetype=Genotype),alpha=0.6)+
  geom_line(data=filter(gdata, Gene_id==target_gene), mapping=aes(linetype=Genotype),alpha=0.7, color='red')+
  #geom_errorbar(mapping=aes(ymax = avg_count+se_count, ymin=avg_count-se_count), alpha=0.7, width=1.5)+
  xlab('Time')+
  ylab('Normalized Count')+
  scale_color_manual(values=c('red', 'grey'))+
  scale_linetype(guide='none')+
  scale_y_continuous(expand = c(0, 0), limits=c(1e-5,ymax))+
  scale_x_continuous(expand = c(0, 0),breaks=c(0,4,8,12,24,36), limits=c(0,37))+
  facet_wrap(~Genotype, labeller = labeller(Genotype=flabels))+
  theme_classic()+
  theme(axis.title.x = element_text(family='Arial', size=12),
        axis.title.y = element_text(family='Arial', size=12),
        axis.text.x = element_text(family='Arial', size=10),
        axis.text.y = element_text(family='Arial', size=10),
        axis.line = element_line(size=1))
