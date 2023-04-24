#############################################################
# Draw a Scatter Plot to show the expression patters of selected genes
# Written by Dylan Barth
# Edited by Jeff Shen
# Modified based on Dylan's emails on to export data for drawing graphs in Excel and changing Y scales
# 02/09/23 Further modified to make graphs based on CGSNL, MSU.ID or RAP.ID
#############################################################
if (!require("ggbreak", quietly = TRUE))
  install.packages("ggbreak")
if (!require("patchwork", quietly = TRUE))
  install.packages("patchwork")


library(tidyverse)
library(extrafont)
library(ggplot2)
library(ggbreak)
library(patchwork)


#import data, format columns 
data <- read_csv("Thesis/RiceRNASeq/Amato/Normalized_2KO.csv")
#data <- read.delim2(file="salmon_RAP112421 embryos three replicates each 090822.txt", head=T, sep="/t")

########################################################
#Step 1 - Propare a new table by deleting unneeded columns
# Which name/id to use? Five choices ((enter abbriviation, e.g. "c", in the following code line
#CGSNL.Gene.Symbol (c)   Gene.symbol.synonym(s)	Shen.Lab (sl)   CGSNL.Gene.Name (g)   	RAP.ID (r)     MSU.ID (m)
# Shen.Lab names = CGSNL names if there is one. If not, use the Gene.symbol.synonym names. If Gene.symbol.synonym contains several names, pick one only for drawing figures.
#gene_name <-  "sl"


#dylan code review: 
#here you're selecting the col with the gene name column by index
#while this may have worked on the first dataset, once the order / num of cols changes
#this strategy fails. Let's implement a solution that's flexible and readable 
#for any students that may read this code in the future
#this way, if there's more/less gene id cols, you can edit gene_cols here

gene_cols <- c("CGSNL.Gene.symbol","Gene.symbol","Shen.Lab","CGSNL.Gene.Name","RAP.ID","MSU.ID")
gene_type <- 'Shen.Lab'
gene_cols <- gene_cols[gene_cols!=gene_type]
data_by_name <- data %>% select(-all_of(gene_cols))
data_by_name <- data_by_name %>% mutate(Gene = unlist(data[,`gene_type`]), .before=1) %>% select(-all_of(gene_type))
data_by_name <- data_by_name %>% separate(Gene, "Gene",extra = "drop", fill = "right", sep=',')
#above command drops everything in Gene col after first id in list 


########################################################
#Step 2- Producing group data for plot; 

####################
#dylan code review:
#this command is specific to the experimental setup of Santi's data
#so we will have to update it to get your data in the correct format

####################
#first remove na
data_by_name <- data_by_name[which(!is.na(data_by_name$MU.DT.rep1_raw)),]

####################
#split data into counts and stats comparing groups
data <- data_by_name %>% select(-matches('log2|pvalue|padj'))
stats <- data_by_name %>% select(matches('log2|pvalue|padj'))
stats$Gene = data$Gene

####################
#now we need to restructure the data because columns contain many variables
data <- data %>% pivot_longer(!Gene, 
                              names_to = c('Genotype','Treatment','Replicate', 'Process'), 
                              values_to = 'count', 
                              names_pattern="(.*)\\.(.*)\\.(.*)_(.*)") #this is a regular expression written especially for victoria's dataset

####################
#filter the data by processing
data <- data %>% filter(Process!='raw')


#########################################################
#Step 3 - Group replicates, get summary statistics

####################
#create id col to group by
data$id <- paste(data$Gene, data$Genotype, data$Treatment, sep='_')

####################
#group by id to merge replicates, average counts and get se, sd
standard_error <- function(x) sd(x) / sqrt(length(x))
data <- data %>% group_by(id) %>% 
  summarise(avg_count = mean(count), across(!c(Replicate))) %>% 
  summarize(se_count = standard_error(count), across()) %>% 
  summarize(sd_count = sd(count), across()) %>% 
  select(-count) %>% 
  unique()

####################
#Optional -- check a given gene name or locus ID is included in the dataframe
subset(data, Gene=="ARP")

###############################################################################
# Step 4 - Filter genes for the summary figure
# genes <- c("AMY1A", "AMY1C", "AMY3A", "AMY5A", "AMY3E", "DOF3")
#genes <- c("WSI18", "EM", "LEA3-1", "HVA22-like", "RAB16D", "RAB16C", "RAB16B", "RAB16A", "HVA22-Ortho")
#genes <- c("GID1", "GID2", "SLR1", "SLRL1", "GAMYB", "GAMYBL1", "GAMYBL2")
#genes <- c("OsRHP1", "OsRING-1", "OsRING151", "OsRING158", "OsRING159", "OsRING277", "OsRING327", "OsRING331", "OsRING389", "OsRING443", "OsRINGC2-2", "OsRING252", "OsRING73")
#genes <- c('OsRHP1', 'OsRING73', 'OsRING443')
#genes <- c("HUB1")
genes <- c('ARP','OsAH1','OsMT3a')
# To assign defined colors to genes, use one of the two options below, Dylan 021623
#data$Gene <- factor(data$Gene, levels= c('gene5','gene2','gene3','gene4','gene1') )
#replace the list of genes with the order you would like them to be in. If the order is already stored in the 'genes' variable, you can simplify the code and just write the following
#If the order is already stored in the 'genes' variable, you can simplify the code and just write the following
  #one extra piece is needed now, filter data by genes
data_sub <- data %>% filter(Gene %in% genes)
data_sub$Gene <- factor(data_sub$Gene, levels=genes)
data_sub$group <- paste(data_sub$Genotype,data_sub$Treatment,sep='_')
data_sub$group <- factor(data_sub$group, levels=c('WT_WW','MU_WW','WT_DT','MU_DT'))
  
###############################################################################
#Step 5 - Set colors; The number of colors must be the same as that of factors in group
#best website for hex color code
#https://www.w3schools.com/colors/colors_hexadecimal.asp
#colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#06d6a0','#118ab2', '#073b4c', '#ff0000', '#c53229', '#fffe29', '#7b2cbf', '#ff476f', '#ffe766', '#000000')
#colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#06d6a0','#118ab2', '#073b4c', '#ff0000', '#c53229', '#fffe29')
#colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#06d6a0','#118ab2', '#073b4c', '#ff0000')
#colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#06d6a0','#118ab2', '#073b4c', '#ff0000')
#colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#06d6a0','#118ab2', '#ff0000')
#colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#06d6a0','#118ab2')
colors <- c('#7b2cbf', '#ef476f', '#ffd166', '#118ab2')
#colors <- c('#7b2cbf', '#ef476f', '#ffd166')
#colors <- c('#ef476f')

###############################################################################
#step 6 - Reduce numbers of zeros in the values on the Y-axis, by dividing the DEseq normalized read values with a number
# Adjust as needed

#reduction <-  1000
# If change this, remember to change the value in "ylab(bquote('Normalized Count ' (x10^3)))+" in the gg block below

#data <- data %>% mutate(orig_avg_count = avg_count, orig_se_count = se_count)
#data_reduced <- data %>% mutate(avg_count = avg_count / reduction, se_count = se_count /reduction)
ymax <- round(max(data_sub$avg_count), 1)*1.05
#May need to adjust the value 1)*1.5 based on the expression levels of selected genes

#Munaully adjust ymax, old way
#set ymax to a number (x) higher than max after "1)+" at the next line. Or the error bar won't be displayed correctly.
#ymax <- round(max(filter(data, Gene %in% genes)$avg_count), 1)+10000  #use thie value for highly expressed genes such as alpha amylases
#ymax <- round(max(filter(data, Gene %in% genes)$avg_count), 1)+5000  #use thie value for genes such as WRKY71 
#ymax <- round(max(filter(data, Gene %in% genes)$avg_count), 1)+3000  #use thie value for lowly expressed gene such as VP1 

###############################################################################
# Step 7 - Plotting
# Show WT data in solid lines and mutant data in dashed lines.
#data_sub$Genotype <- factor(data_sub$Genotype,levels=c('WT','MU'))

# Important!!!
#make sure that the value for in the line of step 8, "ylab('Normalized Count (x1000000)')+" or "ylab(bquote('Normalized Count ' (x10^5)))+"
# is the same as the reduction value in step 6.

# Also adjust the parameters if necessary
# geom_point size value for size of the data points in the graph
# geom_errorbar size value for thickness of the error bars
# theme element_text size for the size of texts - 4 values total in this script
# theme axis.line size for the size of x and y axes
# geom_errorbar "width=" value to adjust the width of error bars and size for the thickness of the error bars
# theme(legend.key.size for the length of the legend bars
yticks <- c(0,5,10,15,20)
gg <- ggplot(filter(data_sub, Gene %in% genes), mapping=aes(x=Gene, y=avg_count, fill=group))+
  geom_bar(stat='identity', pos='dodge')+
  scale_fill_manual(values=colors)+
  geom_errorbar(mapping=aes(ymax = avg_count+se_count, ymin=avg_count-se_count), 
                position = position_dodge(0.9),
                width=0.7, alpha=0.7)+
  geom_point(alpha=0.9, size=0.7, position=position_dodge(0.9), show.legend = FALSE)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)), limits=c(0,ymax), breaks=yticks)+
  theme_classic()+
  theme(axis.title.x = element_text(family='Arial', size=24),
        axis.title.y = element_text(family='Arial', size=24),
        axis.text.x = element_text(family='Arial', size=24),
        axis.text.y = element_text(family='Arial', size=24),
        axis.line = element_line(size=2),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  xlab('Gene')+
  ylab('Normalized Count')+
  guides(fill=guide_legend(title="Treatment"))






#show the graph. To save a plot as a PNG file, use the Export function at the bottom right panel
#or just use ggsave after running gg
gg
ggsave('Downloads/test_graph.png')


##################################################
# if using an axis break, use this plot
#set the range you wish to compress and tick marks
brange <- c(12,19)
yticks <- c(0,5,10,15,20)

gg <- ggplot(filter(data_sub, Gene %in% genes), mapping=aes(x=Gene, y=avg_count, fill=group))+
  geom_bar(stat='identity', pos='dodge')+
  scale_fill_manual(values=colors)+
  geom_errorbar(mapping=aes(ymax = avg_count+se_count, ymin=avg_count-se_count), 
                position = position_dodge(0.9),
                width=0.7, alpha=0.7)+
  geom_point(alpha=0.9, size=0.7, position=position_dodge(0.9), show.legend = FALSE)+
  scale_y_break(brange, scale=.25, space=0.05)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)), limits=c(0,ymax), breaks=yticks)+
  theme_classic()+
  theme(axis.title.x = element_text(family='Arial', size=24),
        axis.title.y = element_text(family='Arial', size=24),
        axis.text.x = element_text(family='Arial', size=24),
        axis.text.y = element_text(family='Arial', size=24),
        axis.line = element_line(size=2),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  xlab('Gene')+
  ylab('Normalized Count')+
  guides(fill=guide_legend(title="Treatment"))
gg
