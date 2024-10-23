#need R/3.6.2 NOT 4

library(tidyr)
library(data.table)
library(qqman)

setwd("pkd/paper/")
trait="all_cystic"

#get gene coordinates
gene <- fread("genes_positions.txt", header=T, data.table=F)
gene <- unique(gene)
gene <- gene[!duplicated(gene$`Gene name`),]
gene<- gene[,c(4,1,2,3)]
names(gene) <- c("CHR", "START", "END", "Gene")
gene$START <- as.numeric(gene$START)
gene$END <- as.numeric(gene$END)
gene$BP <- (gene$START+gene$END)/2

#nvd lof
x <- fread("saige_LoF_pheno_nmd_v15.txt_for_export.txt")
x <- separate(x, col="Gene", into=c("Gene", "ENS"), sep="\\_")
b <- x[,c(1,3)]
c <- merge(b, gene, by="Gene", all.x=T)
c <- c[,c(1:3,6)]
names(c) <- c("SNP", "P", "CHR", "BP")
c$CHR=gsub("chrX","chr23",c$CHR)
c$CHR=as.numeric(gsub("chr","",c$CHR))
c <- c[!is.na(c$P), ]
c$BP <- as.numeric(c$BP)
c$P <- as.numeric(c$P)
d <- na.omit(c)
nvd_lof <- d

#nvd missense
x <- fread("for_supplementary/missense+_unsolved_no_sv_maf0.001_cadd20.SAIGE.gene.txt_for_export.txt")
x <- separate(x, col="Gene", into=c("Gene", "ENS"), sep="\\_")
b <- x[,c(1,3)]
c <- merge(b, gene, by="Gene", all.x=T)
c <- c[,c(1:3,6)]
names(c) <- c("SNP", "P", "CHR", "BP")
c$CHR=gsub("chrX","chr23",c$CHR)
c$CHR=as.numeric(gsub("chr","",c$CHR))
c <- c[!is.na(c$P), ]
c$BP <- as.numeric(c$BP)
c$P <- as.numeric(c$P)
d <- na.omit(c)
nvd_missense <- d

#all lof
x <- fread("saige_LoF_pheno_v15_all_cystic.txt_for_export.txt")
x <- separate(x, col="Gene", into=c("Gene", "ENS"), sep="\\_")
b <- x[,c(1,3)]
c <- merge(b, gene, by="Gene", all.x=T)
c <- c[,c(1:3,6)]
names(c) <- c("SNP", "P", "CHR", "BP")
c$CHR=gsub("chrX","chr23",c$CHR)
c$CHR=as.numeric(gsub("chr","",c$CHR))
c <- c[!is.na(c$P), ]
c$BP <- as.numeric(c$BP)
c$P <- as.numeric(c$P)
d <- na.omit(c)
all_lof <- d
p=length(unique(d$SNP))

#all missense
x <- fread("missense+_all_cystic_maf0.001_cadd20.SAIGE.gene.txt_for_export.txt")
x <- separate(x, col="Gene", into=c("Gene", "ENS"), sep="\\_")
b <- x[,c(1,3)]
c <- merge(b, gene, by="Gene", all.x=T)
c <- c[,c(1:3,6)]
names(c) <- c("SNP", "P", "CHR", "BP")
c$CHR=gsub("chrX","chr23",c$CHR)
c$CHR=as.numeric(gsub("chr","",c$CHR))
c <- c[!is.na(c$P), ]
c$BP <- as.numeric(c$BP)
c$P <- as.numeric(c$P)
d <- na.omit(c)
all_missense <- d

#splice
splice <- fread("spliceAI_splie_summary_stats_all_combined.txt")
splice[13830,1]<- "PKD1 Splice"
splice[13839,1]<- "PKD2 Splice"
splice <- splice[-c(13831,13832),]

#sv
sv <- fread("pkdmaf<0.001_filtered_exome_wide_SV_summary_stats.txt")
sv <- sv[,c(1,7)]
c <- merge(sv, gene, by.x="GENE",by.y="Gene", all.x=T)
c <- c[,c(1,2,3,6)]
names(c) <- c("SNP", "P", "CHR", "BP")
c$CHR=gsub("chrX","chr23",c$CHR)
c$CHR=as.numeric(gsub("chr","",c$CHR))
c <- c[!is.na(c$P), ]
c$BP <- as.numeric(c$BP)
c$P <- as.numeric(c$P)
d <- na.omit(c)
sv <- d
sv <- sv[order(sv$P),]
sv[3,1] <- "17q12 CNV"
sv[1,1]<- "PKD1 SV"
sv[2,1]<- "PKD2 SV"
sv[c(4:22),1] <- ""
#rename all the variables
all_lof$SNP <- paste(all_lof$SNP,"LoF")
all_missense$SNP <- paste(all_missense$SNP,"Missense+")
nvd_lof$SNP <- paste(nvd_lof$SNP,"NVD LoF")
nvd_missense$SNP <- paste(nvd_missense$SNP,"NVD Missense+")
#sv and splice already done 
all <- rbind(all_lof, all_missense, nvd_lof, nvd_missense, sv, splice)
all<- all[order(all$P),]
all$SNP <- ifelse(all$P>1.257990e-06, "", all$SNP)

# Plot Manhattan
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("ggmanh")

library(ggmanh)
library(ggrepel)
png("for_abstract/all_associations_pkd_for_visual_abstract.png", type="cairo", units="in", res=600, width=7, height=5)
manhattan_plot(x = all, pval.colname = "P", 
               chr.colname = "CHR", 
               pos.colname = "BP", 
               plot.title = "", 
               label.colname="SNP",
               rescale.ratio.threshold = 0.9,
               force=32,
               max.overlaps=100,
               rescale=T, 
               signif=1.257990e-06,
               point.size = 0.9,
               label.font.size = 3)
dev.off()

#bigger
manhattan_plot(x = all, pval.colname = "P", 
               chr.colname = "CHR", 
               pos.colname = "BP", 
               plot.title = "", 
               y.label="Association with Cystic Kidney Disease (-log10(P))",
               label.colname="SNP",
               rescale.ratio.threshold = 0.9,
               force=45,
               max.overlaps=100,
               rescale=T, 
               signif=1.257990e-06,
               point.size = 1.5,
               label.font.size = 3.5)


#qqman
manhattan(all, ylim=c(0,120), cex.axis=0.6, col = c("darkorange", "deepskyblue3"), annotateTop = FALSE,annotatePval = 0.00005,genomewideline =5, suggestiveline = FALSE)

#full plot
png("hes_stones_missense_saigegene_manhattan_maf0.005.png", type="cairo", units="in", res=300, width=12, height=6)
par(omi=c(0,0.5,0,0))
#manhattan(d, ylim=c(0,10), cex.axis=0.6, col = c("darkorange", "deepskyblue3"), annotateTop = FALSE,annotatePval = 0.00005,genomewideline =-log10(0.05/p), suggestiveline = FALSE)
manhattan(d, ylim=c(0,12), cex.axis=0.6, col = c("darkorange", "deepskyblue3"), annotateTop = FALSE,annotatePval = 0.00005,genomewideline =5, suggestiveline = FALSE)
dev.off()


###ggplot
library(ggplot2)
library(ggmanh)
library(tidyverse)
library(data.table)
library(ggtext)
library(ggbreak)
library(GWASTools)

#prep
data_cum <- d %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHR, bp_add)

d <- d %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add )
#set plotting axes
axis_set <- d%>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- d%>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)

#0.05/length(summary stats)
sig <- 0.05/p

#plotting
manhplot <- ggplot(d, aes(x = BP, y = -log10(d$P), 
                             color = as_factor(CHR))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(P-value)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

manhplot 
t <- manhplot + scale_y_cut(breaks=c(0,10)) 
t
y <- t + scale_y_cut(breaks=c(149,150)) 
y
o<- y + scale_y_cut(breaks=c(300,310), space = 0.2)
o

#qq
#Genomic inflation factor function
lambda<-function(pvalues){
  chisq <- qchisq(1-pvalues,1)
  lambda=median(chisq,na.rm=TRUE)/qchisq(0.5,1)
  return(lambda)
}
png("metal_metanalysis_qqplot.png", type = "cairo", units="in", res=300, width=6, height=6)
par(omi=c(0,0.5,0,0))
qqPlot(gwas$`P-value`,main=paste0("lambda=",round(lambda(gwas$`P-value`),4)), cex=0.6)
dev.off()

qq(gwas$`P-value`)



