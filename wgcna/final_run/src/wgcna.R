library(tidyverse)
library(magrittr)
library(WGCNA)
library(DESeq2)
library(genefilter)

#import data
data <- readr::read_delim('raw_data/fC_all_MSU_v2.txt', delim='\t')

#clean data, prepare for input to wgcna
data <- select(data, -c('Chr','Start','End','Strand','Length'))

#only keep replicates in meta file (same as their DESeq analysis)
replicates <- readr::read_delim('raw_data/Meta_v2.txt', delim='\t')
keepers <- as.list(replicates[1])$X1
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


norm_expr_plot <- expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

ggsave(plot=norm_expr_plot, filename='figs/norm_expr_plot.nofilter.png', device='png')


#########################################
###       prepare data for WGCNA      ###
#########################################

#multithread
allowWGCNAThreads(40)

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

#plot scale independence
png(file='figs/scale_independence.nofilter.png',width=600,height=600)
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
dev.off()

#plot connectivity
png(file='figs/mean_connectivity.nofilter.png',width=600, height=600)
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
dev.off()



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

                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",

                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

#plot results
mergedColors = labels2colors(netwk$colors)
png(file='figs/dendrogram.nofilter.png', width=1000, height=600)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
dev.off()


#write file connecting genes to their modules (colors)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

write_delim(module_df,
            path = "data/gene_modules.nofilter.txt",
            delim = "\t")


#########################################
###     post-analysis processing      ###
#########################################

#connect modules (colors) to treatment groups
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )


#str(mME)
#Do a linear modeling test to see which modules are significant for time and genotype
mME <- mutate(mME, genotype = factor(substr(treatment, 1, 1)))
mME$time <- as.integer(str_replace(mME$treatment, pattern = ".*e([^_]*)_.*", replacement = "\\1") )

#summary(mME$genotype)
#summary(mME$time)

sink(file='data/glm_module_output.txt')
for (module in unique(mME$name)){
	data <- filter(mME, name==module)
	model <- glm(formula= value ~ genotype + time + genotype*time, data=data)
	print(module)
	print(summary(model))
}
sink()

#plot heatmap tying expression modules to treatment groups
level_order <- c('me0_1','me0_2','me0_3','me4_1','me4_2','me4_3','me8_1','me8_2','me8_3','me12_1','me12_2','me12_3','me24_1','me24_2','me24_4r','me36_1','me36_2r','me36_3r',
	'we0_1','we0_2','we0_3','we4_1','we4_2','we4_3','we8_1','we8_2','we8_3','we12_1','we12_2','we12_3','we24_1','we24_2','we24_3','we36_1r','we36_2r','we36_3')
mod_treat_heatmap <- mME %>% ggplot(., aes(x=factor(treatment, level = level_order), y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

ggsave(plot=mod_treat_heatmap, file='figs/treatment_heatmap.nofilter.png', device='png')

# pick out a few modules of interest here using the output from linear modeling
modules_of_interest = c("pink", "blue", "salmon", "magenta","purple")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

#dissect name vble into time and genotype
submod_df <- mutate(submod_df, genotype = factor(substr(name, 1, 1)))
submod_df$time <- as.integer(str_replace(submod_df$name, pattern = ".*e([^_]*)_.*", replacement = "\\1") )

#plot expression levels amongst treatments for specific modules
mod_expr_plt <- submod_df %>% ggplot(., aes(x=time, y=value, group=gene_id)) +
  geom_line(aes(color = genotype),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

ggsave(plot=mod_expr_plt, file='figs/module_norm_expression.nofilter.png', device='png')




