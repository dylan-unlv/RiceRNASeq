library(tidyverse)
library(riceidconverter)

#translate DEG from MSU 2 RAP gene id format

data <- read_delim('DEG/DEG_0hr.csv', delim='\t')
updata <- filter(data, log2FoldChange>0)
downdata <- filter(data, log2FoldChange<0)
uptrans <- RiceIDConvert(updata$gene, 'MSU', 'RAP')[,2] %>% na.omit 
downtrans <- RiceIDConvert(downdata$gene, 'MSU', 'RAP')[,2] %>% na.omit
write_lines(uptrans,'DEG/upregulated/DEG_0hr.MSU2RAP.txt')
write_lines(downtrans,'DEG/downregulated/DEG_0hr.MSU2RAP.txt')

data <- read_delim('DEG/DEG_4hr.csv', delim='\t')
updata <- filter(data, log2FoldChange>0)
downdata <- filter(data, log2FoldChange<0)
uptrans <- RiceIDConvert(updata$gene, 'MSU', 'RAP')[,2] %>% na.omit 
downtrans <- RiceIDConvert(downdata$gene, 'MSU', 'RAP')[,2] %>% na.omit
write_lines(uptrans,'DEG/upregulated/DEG_4hr.MSU2RAP.txt')
write_lines(downtrans,'DEG/downregulated/DEG_4hr.MSU2RAP.txt')

data <- read_delim('DEG/DEG_8hr.csv', delim='\t')
updata <- filter(data, log2FoldChange>0)
downdata <- filter(data, log2FoldChange<0)
uptrans <- RiceIDConvert(updata$gene, 'MSU', 'RAP')[,2] %>% na.omit 
downtrans <- RiceIDConvert(downdata$gene, 'MSU', 'RAP')[,2] %>% na.omit
write_lines(uptrans,'DEG/upregulated/DEG_8hr.MSU2RAP.txt')
write_lines(downtrans,'DEG/downregulated/DEG_8hr.MSU2RAP.txt')

data <- read_delim('DEG/DEG_12hr.csv', delim='\t')
updata <- filter(data, log2FoldChange>0)
downdata <- filter(data, log2FoldChange<0)
uptrans <- RiceIDConvert(updata$gene, 'MSU', 'RAP')[,2] %>% na.omit 
downtrans <- RiceIDConvert(downdata$gene, 'MSU', 'RAP')[,2] %>% na.omit
write_lines(uptrans,'DEG/upregulated/DEG_12hr.MSU2RAP.txt')
write_lines(downtrans,'DEG/downregulated/DEG_12hr.MSU2RAP.txt')

data <- read_delim('DEG/DEG_24hr.csv', delim='\t')
updata <- filter(data, log2FoldChange>0)
downdata <- filter(data, log2FoldChange<0)
uptrans <- RiceIDConvert(updata$gene, 'MSU', 'RAP')[,2] %>% na.omit 
downtrans <- RiceIDConvert(downdata$gene, 'MSU', 'RAP')[,2] %>% na.omit
write_lines(uptrans,'DEG/upregulated/DEG_24hr.MSU2RAP.txt')
write_lines(downtrans,'DEG/downregulated/DEG_24hr.MSU2RAP.txt')

data <- read_delim('DEG/DEG_36hr.csv', delim='\t')
updata <- filter(data, log2FoldChange>0)
downdata <- filter(data, log2FoldChange<0)
uptrans <- RiceIDConvert(updata$gene, 'MSU', 'RAP')[,2] %>% na.omit 
downtrans <- RiceIDConvert(downdata$gene, 'MSU', 'RAP')[,2] %>% na.omit
write_lines(uptrans,'DEG/upregulated/DEG_36hr.MSU2RAP.txt')
write_lines(downtrans,'DEG/downregulated/DEG_36hr.MSU2RAP.txt')