library(dplyr)
library(tidyr)
library(seqinr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)

my.cols <- brewer.pal(5, "Blues")

data1 = read.delim("/Volumes/data/AbX/experiments/170323/IMGT_download/AK170_2_S2_IgG1_ECBC21_comb/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
data2 = read.delim("/Volumes/data/AbX/experiments/170323/IMGT_download/AK170_2_S2_IgG2_ECBC21/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
data3 = read.delim("/Volumes/data/AbX/experiments/170323/IMGT_download/AK170_2_S2_IgG3_ECBC21/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
data4 = read.delim("/Volumes/data/AbX/experiments/170323/IMGT_download/AK170_2_S2_IgG4_ECBC21/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
dataL = read.delim("/Volumes/data/AbX/experiments/170323/IMGT_download/AK170_2_S2_l_ECBC21_comb/1_Summary.txt", header=T, sep = "\t", fill = TRUE)

#Analysis for IgG1 sample, last line has to be commented out if all reads are to be considered
HC1 = data1 %>%
  mutate(V_gene=gsub("Homsap ", "", V.GENE.and.allele)) %>%
  select(Sequence.ID, Functionality, V_gene) %>%
  separate(V_gene, c("var1", "var2"), sep = ",") %>%
  separate(var1, c("gene", "allele"), sep = "\\*") %>%
  mutate(family=gene) %>%
  separate(family, c("family", "rest"), sep = "\\-") %>%
  mutate(gene=gsub("IGHV", "", gene)) %>%
  select(family, gene, Functionality) 

HC1_fam = HC1 %>%
  group_by(family) %>%
  summarise(sum_fam=n()) %>%
  mutate(subtype="IgG1") %>%
  mutate(sum_reads=sum(sum_fam))

HC1_gene = HC1 %>%
  group_by(gene) %>%
  summarise(sum_gene=n()) %>%
  mutate(subtype="IgG1") %>%
  mutate(sum_reads=sum(sum_gene))

#Figures for IgG1 sample
fig_HC1_fam = ggplot(HC1_fam, aes(x=family, y=sum_fam)) +
  geom_bar(stat="identity", color="steelblue", fill = "steelblue") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  scale_y_log10() +
  scale_x_discrete(labels=c("no results", "IGHV1", "IGHV3", "IGLV1")) +
  xlab("\nVH Family") +
  ylab("Number of Reads\n") +
  ggtitle("Reads per Family (tp2, IgG1)") 

fig_HC1_fam

fig_HC_gene = ggplot(HC1_gene, aes(x=gene, y=sum_gene)) +
  geom_bar(stat="identity", color="steelblue", fill = "steelblue") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  scale_y_log10() +
  xlab("\nVH Gene") +
  ylab("Number of Reads\n") +
  ggtitle("Reads per Gene (tp2, IgG1)") 

fig_HC1_gene

#Analysis for other IgG subtypes
HC2 = data2 %>%
  mutate(V_gene=gsub("Homsap ", "", V.GENE.and.allele)) %>%
  select(Sequence.ID, Functionality, V_gene) %>%
  separate(V_gene, c("var1", "var2"), sep = ",") %>%
  separate(var1, c("gene", "allele"), sep = "\\*") %>%
  mutate(family=gene) %>%
  separate(family, c("family", "rest"), sep = "\\-") %>%
  mutate(gene=gsub("IGHV", "", gene)) %>%
  select(family, gene, Functionality) 

HC2_fam = HC2 %>%
  group_by(family) %>%
  summarise(sum_fam=n()) %>%
  mutate(subtype="IgG2") %>%
  mutate(sum_reads=sum(sum_fam))

HC2_gene = HC2 %>%
  group_by(gene) %>%
  summarise(sum_gene=n()) %>%
  mutate(subtype="IgG2") %>%
  mutate(sum_reads=sum(sum_gene))

HC3 = data3 %>%
  mutate(V_gene=gsub("Homsap ", "", V.GENE.and.allele)) %>%
  select(Sequence.ID, Functionality, V_gene) %>%
  separate(V_gene, c("var1", "var2"), sep = ",") %>%
  separate(var1, c("gene", "allele"), sep = "\\*") %>%
  mutate(family=gene) %>%
  separate(family, c("family", "rest"), sep = "\\-") %>%
  mutate(gene=gsub("IGHV", "", gene)) %>%
  select(family, gene, Functionality) 

HC3_fam = HC3 %>%
  group_by(family) %>%
  summarise(sum_fam=n()) %>%
  mutate(subtype="IgG3") %>%
  mutate(sum_reads=sum(sum_fam))

HC3_gene = HC3 %>%
  group_by(gene) %>%
  summarise(sum_gene=n()) %>%
  mutate(subtype="IgG3") %>%
  mutate(sum_reads=sum(sum_gene))

HC4 = data4 %>%
  mutate(V_gene=gsub("Homsap ", "", V.GENE.and.allele)) %>%
  select(Sequence.ID, Functionality, V_gene) %>%
  separate(V_gene, c("var1", "var2"), sep = ",") %>%
  separate(var1, c("gene", "allele"), sep = "\\*") %>%
  mutate(family=gene) %>%
  separate(family, c("family", "rest"), sep = "\\-") %>%
  mutate(gene=gsub("IGHV", "", gene)) %>%
  select(family, gene, Functionality) 

HC4_fam = HC4 %>%
  group_by(family) %>%
  summarise(sum_fam=n()) %>%
  mutate(subtype="IgG4") %>%
  mutate(sum_reads=sum(sum_fam))

HC4_gene = HC4 %>%
  group_by(gene) %>%
  summarise(sum_gene=n()) %>%
  mutate(subtype="IgG4") %>%
  mutate(sum_reads=sum(sum_gene))

#Combined figure
HC12_fam = full_join(HC1_fam, HC2_fam)
HC123_fam = full_join(HC12_fam, HC3_fam)

HC_all_fam = full_join(HC123_fam, HC4_fam) %>%
  mutate(family=sub("^$", "No Results", family)) %>%
  spread(family, sum_fam) %>%
  replace(is.na(.), 0) %>%
  gather(family, sum_fam, 3:6) %>%
  group_by(subtype) %>%
  arrange(subtype) %>%
  mutate(freq=sum_fam/sum_reads) 

HC12_gene = full_join(HC1_gene, HC2_gene)
HC123_gene = full_join(HC12_gene, HC3_gene) 

HC_all_gene = full_join(HC123_gene, HC4_gene) %>%
  mutate(gene=sub("^$", "No Results", gene)) %>%
  spread(gene, sum_gene) %>%
  replace(is.na(.), 0) %>%
  gather(gene, sum_gene, 3:43) %>%
  group_by(subtype) %>%
  arrange(subtype) %>%
  mutate(freq=sum_gene/sum_reads)

fig_HC_all_fam = ggplot(HC_all_fam, aes(x=family, y=sum_fam, fill=subtype)) +
  geom_bar(stat="identity", position="dodge") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  #scale_y_log10(breaks =c(1,100,10000,1000000), labels=c("1","100","10`000","1`000`000")) +
  #scale_x_discrete(labels=c("no results", "IGHV1", "IGHV3", "IGLV1")) +
  xlab("VH Family") +
  ylab("Number of Reads") +
  ggtitle("Reads per Family (wk89)") 

fig_HC_all_fam

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/wk89_Reads_per_Family.pdf", plot=fig_HC_all_fam, width = 30/2.54, height = 21/2.54)


fig_HC_all_fam_freq = ggplot(HC_all_fam, aes(x=family, y=freq, fill=subtype)) +
  geom_bar(stat="identity", position="dodge") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  #scale_x_discrete(labels=c("no results", "IGHV1", "IGHV3", "IGLV1")) +
  xlab("VH Family") +
  ylab("Number of Reads") +
  ggtitle("Reads per Family (tp2, IgG1)") 

fig_HC_all_fam_freq


fig_HC_all_gene = ggplot(HC_all_gene, aes(x=gene, y=sum_gene, fill=subtype)) +
  geom_bar(stat="identity", position="dodge") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  #scale_y_log10(breaks =c(1,100,10000,1000000), labels=c("1","100","10`000","1`000`000")) +
  #scale_x_discrete(labels=c("no results", "IGHV1", "IGHV3", "IGLV1")) +
  xlab("VH Gene (IGHV)") +
  ylab("Number of Reads") +
  ggtitle("Reads per Gene (wk89)") 

fig_HC_all_gene

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/wk89_Reads_per_Gene.pdf", plot=fig_HC_all_gene, width = 30/2.54, height = 21/2.54)

fig_HC_all_gene_freq = ggplot(HC_all_gene, aes(x=gene, y=freq, fill=subtype)) +
  geom_bar(stat="identity", position="dodge") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust=1, size = 20, vjust = 0.5), 
        text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  #scale_x_discrete(labels=c("no results", "IGHV1", "IGHV3", "IGLV1")) +
  xlab("VH Gene (IGHV)") +
  ylab("Relative Number of Reads per Gene\n") +
  ggtitle("Reads per Gene (tp2)") 

fig_HC_all_gene_freq

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/tp2_Reads_per_Gene_freq.pdf", plot=fig_HC_all_gene_freq, width = 30/2.54, height = 21/2.54)

#Analysis for IgL sample
LC = dataL %>%
  mutate(V_gene=gsub("Homsap ", "", V.GENE.and.allele)) %>%
  select(Sequence.ID, Functionality, V_gene) %>%
  separate(V_gene, c("var1", "var2"), sep = ",") %>%
  separate(var1, c("gene", "allele"), sep = "\\*") %>%
  mutate(family=gene) %>%
  separate(family, c("family", "rest"), sep = "\\-") %>%
  mutate(gene=gsub("IGLV", "", gene)) %>%
  select(family, gene, Functionality) %>%
  mutate(gene=sub("^$", "No Results", gene)) %>%
  mutate(family=sub("^$", "No Results", family))

LC_fam = LC %>%
  group_by(family) %>%
  summarise(sum_fam=n()) %>%
  mutate(subtype="IgL") %>%
  mutate(sum_reads=sum(sum_fam))

LC_gene = LC %>%
  group_by(gene) %>%
  summarise(sum_gene=n()) %>%
  mutate(subtype="IgL") %>%
  mutate(sum_reads=sum(sum_gene)) %>%
  mutate(freq=sum_gene/sum_reads)

#Figures for IgL sample
fig_LC_fam = ggplot(LC_fam, aes(x=family, y=sum_fam)) +
  geom_bar(stat="identity", color="steelblue", fill = "steelblue") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  #scale_y_log10(breaks =c(1,100,10000,1000000), labels=c("1","100","10`000","1`000`000")) +
  #scale_x_discrete(labels=c("no results", "IGHV1", "IGHV3", "IGLV1")) +
  xlab("\nVL Family") +
  ylab("Number of Reads\n") +
  ggtitle("Reads per Family (wk89, IgL)") 

fig_LC_fam

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/wk89_Reads_per_Fam_LC.pdf", plot=fig_LC_fam, width = 30/2.54, height = 21/2.54)

fig_LC_gene_freq = ggplot(LC_gene, aes(x=gene, y=freq)) +
  geom_bar(stat="identity", color="steelblue", fill = "steelblue") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  xlab("\nVL Gene") +
  ylab("Number of Reads\n") +
  ggtitle("Reads per Gene (tp2, IgL)") 

fig_LC_gene_freq

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/tp2_Reads_per_Gene_LC.pdf", plot=fig_LC_gene, width = 30/2.54, height = 21/2.54)

# fig LC frequencies
fig_LC_gene = ggplot(LC_gene, aes(x=gene, y=sum_gene)) +
  geom_bar(stat="identity", color="steelblue", fill = "steelblue") +
  theme(plot.title=element_text(size=24), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16, vjust = 0.5), text = element_text(size=24), 
        axis.text.x=element_text(size=50), 
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70"),
        legend.text = element_text(size = 18)) +
  #scale_y_log10(breaks =c(1,100,10000,1000000), labels=c("1","100","10`000","1`000`000")) +
  xlab("VL Gene") +
  ylab("Number of Reads") +
  ggtitle("Reads per Gene (wk89, IgL)") 

fig_LC_gene
ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/wk89_Reads_per_Gene_LC.pdf", plot=fig_LC_gene, width = 30/2.54, height = 21/2.54)


plot_grid(fig_HC_all_fam, fig_HC_all_gene, fig_LC_fam, fig_LC_gene, ncol=2, nrow=2 )

 #############
group1 = read.delim("/Volumes/data/AbX/ECBC/ECBC21/IMGT_grouped/2_IgG1/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
group2 = read.delim("/Volumes/data/AbX/ECBC/ECBC21/IMGT_grouped/2_IgG2/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
group3 = read.delim("/Volumes/data/AbX/ECBC/ECBC21/IMGT_grouped/2_IgG3/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
groupL = read.delim("/Volumes/data/AbX/ECBC/ECBC21/IMGT_grouped/2_IgL/1_Summary.txt", header=T, sep = "\t", fill = TRUE)

gr1 = group1 %>%
  select(Sequence.ID) %>%
  separate(Sequence.ID, c("reads", "ID"), sep= "_reads-") %>%
  mutate(reads=as.numeric(reads))

gr2 = group2 %>%
  select(Sequence.ID) %>%
  separate(Sequence.ID, c("reads", "ID"), sep= "_reads-") %>%
  mutate(reads=as.numeric(reads))

gr3 = group3 %>%
  select(Sequence.ID) %>%
  separate(Sequence.ID, c("reads", "ID"), sep= "_reads-") %>%
  mutate(reads=as.numeric(reads))

grL = groupL %>%
  select(Sequence.ID) %>%
  separate(Sequence.ID, c("reads", "ID"), sep= "_reads-") %>%
  mutate(reads=as.numeric(reads))
  
#select IGHV3-23 and IGHJ6 for HCs, IGLV2-23 and IGLJ3 or IGLJ5
raw1 = data1 %>%
  select(Sequence.ID, Functionality, V.GENE.and.allele, J.GENE.and.allele) %>%
  filter(grepl("^productive", Functionality)) %>%
  filter(grepl("IGHJ6", J.GENE.and.allele)) %>%
  filter(grepl("IGHV3-23", V.GENE.and.allele))

nn1=sum(nrow(raw1))

raw2 = data2 %>%
  select(Sequence.ID, Functionality, V.GENE.and.allele, J.GENE.and.allele) %>%
  filter(grepl("^productive", Functionality)) %>%
  filter(grepl("IGHJ6", J.GENE.and.allele)) %>%
  filter(grepl("IGHV3-23", V.GENE.and.allele))

nn2=sum(nrow(raw2))

raw3 = data3 %>%
  select(Sequence.ID, Functionality, V.GENE.and.allele, J.GENE.and.allele) %>%
  filter(grepl("^productive", Functionality)) %>%
  filter(grepl("IGHJ6", J.GENE.and.allele)) %>%
  filter(grepl("IGHV3-23", V.GENE.and.allele))

nn3=sum(nrow(raw3))

rawL = dataL %>%
  select(Sequence.ID, Functionality, V.GENE.and.allele, J.GENE.and.allele) %>%
  filter(grepl("^productive", Functionality)) %>%
  filter(grepl("IGLJ3|IGLJ5", J.GENE.and.allele)) %>%
  filter(grepl("IGLV2-23", V.GENE.and.allele))

nnL=sum(nrow(rawL))
  
#these are numbers for dataframe to avoid reading large files again
n1=93524
n2=13061
n3=15460
nL=107531

#create dataframe for figure
Analysis=c("ECBC", "all", "ECBC", "all", "ECBC", "all", "ECBC", "all")
subtype=c("IgG1", "IgG1", "IgG2", "IgG2", "IgG3", "IgG3", "IgL", "IgL") 
no_reads=c(sum(gr1$reads), n1, sum(gr2$reads), n2, sum(gr3$reads), n3, sum(grL$reads), nL)

df = data.frame(Analysis, subtype, no_reads)

fig_a = ggplot(df, (aes(y=no_reads, x=subtype, fill=Analysis))) +
  geom_bar(stat="identity", position="dodge") +
  theme(plot.title=element_text(size=18), 
        axis.text=element_text(size=18), 
        axis.title=element_text(size=24),
        legend.text = element_text(size = 18),
        legend.title=element_text(size=18),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  xlab("\nSubtype") +
  ylab("Number of Reads\n") +
  ggtitle("Reads included in ECBC consensus (tp2)") 

fig_a

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/tp2_Reads_in_ECBCs.pdf", plot=fig_a, width = 30/2.54, height = 21/2.54)

#data for relative figure
subtype_rel=c("IgG1", "IgG2", "IgG3", "IgL") 
no_reads_rel=c(sum(gr1$reads)/n1, sum(gr2$reads)/n2, sum(gr3$reads)/n3, sum(grL$reads)/nL)

df_rel=data.frame(subtype_rel, no_reads_rel)

fig_rel = ggplot(df_rel, (aes(x=subtype_rel, y=no_reads_rel))) +
  geom_bar(stat="identity", fill = "steelblue", width=0.5) +
  theme(plot.title=element_text(size=16), 
        axis.text=element_text(size=14), 
        axis.title=element_text(size=14),
        legend.text = element_text(size = 14),
        legend.title=element_text(size=14),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  xlab("\nSubtype") +
  ylab("relative #Reads in ECBCS\n") +
  ggtitle("Reads included in ECBC consensus (tp2)") 
  
fig_rel

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/tp2_rel_Reads_in_ECBCs.pdf", plot=fig_rel, width = 30/2.54, height = 21/2.54)

#read stats from excel file (numbers extracted from fasta and 1_Summary files by grep, awk etc)
summary_file = read.table("/Volumes/data/AbX/ECBC/ECBC21/read_stats.txt", header=T, sep = "\t", fill = TRUE)

data_clean = summary_file %>%
  filter(tp>0) 

data_clean$tp <- factor(data_clean$tp, levels=c("46", "89", "115", "179"))

reads_tp = ggplot(data_clean, aes(x=subtype, y=reads_panda, fill=tp)) +
  geom_bar(stat="identity", position= "dodge") +
  theme(plot.title=element_text(size=14), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) +
  xlab("") +
  ylab("Reads after panda") +
  ggtitle("Reads after panda") 

reads_tp

reads_tp_prod = ggplot(data_clean, aes(x=subtype, y=reads_productive, fill=tp)) +
  geom_bar(stat="identity", position= "dodge") +
  theme(plot.title=element_text(size=14), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) +
  ylab("Productive Reads") +
  xlab("") +
  ggtitle("Productive Ab Reads") 

reads_tp_prod

reads_tp_prod_rel = ggplot(data_clean, aes(x=subtype, y=rel_reads_prod, fill=tp)) +
  geom_bar(stat="identity", position= "dodge") +
  theme(plot.title=element_text(size=14), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  scale_y_continuous(limits=c(0, 100)) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) +
  ylab("% Productive Reads") +
  xlab("") +
  ggtitle("Productive Ab Reads") 

reads_tp_prod_rel

reads_tp_subtype = ggplot(data_clean, aes(x=subtype, y=rel_reads_subtype, fill=tp)) +
  geom_bar(stat="identity", position= "dodge") +
  theme(plot.title=element_text(size=14), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) +
  scale_x_discrete(limits=c("IgG1", "IgG2", "IgG3", "IgG4")) +
  ylab("% Reads IgG subtype") +
  xlab("") +
  ggtitle("IgG Subtype distribution") 

reads_tp_subtype

reads_tp_VJ_freq = ggplot(data_clean, aes(x=subtype, y=rel_reads_V_J, fill=tp)) +
  geom_bar(stat="identity", position= "dodge") +
  theme(plot.title=element_text(size=14), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) +
  scale_y_continuous(limits=c(0, 16)) +
  ylab("% Reads with correct V-, J-Gene\nof productive Reads") +
  xlab("") +
  ggtitle("Productive Ab Reads, filtered for V- and J-Gene") 

reads_tp_VJ_freq

plot_grid(reads_tp, reads_tp_subtype, reads_tp_prod, reads_tp_prod_rel, ncol=2, nrow=2 )
plot_grid(reads_tp_VJ_freq, reads_tp_VJ_CDRlength, reads_tp_VJ_CDRlength_freq, ncol=2, nrow=2 )

reads_tp_VJ_CDRlength = ggplot(data_clean, aes(x=subtype, y=reads_prod_V_J_CDR3length, fill=tp)) +
  geom_bar(stat="identity", position= "dodge") +
  theme(plot.title=element_text(size=14), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) +
  #scale_y_continuous(limits=c(0, 16)) +
  ylab("No. of Reads with \ncorrect V-, J-Gene, CDR3 length") +
  xlab("") +
  ggtitle("Productive Ab Reads, filtered for V-, J-Gene and CDR3 length") 

reads_tp_VJ_CDRlength

reads_tp_VJ_CDRlength_freq = ggplot(data_clean, aes(x=subtype, y=reads_prod_V_J_CDR3length_rel, fill=tp)) +
  geom_bar(stat="identity", position= "dodge") +
  theme(plot.title=element_text(size=14), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color="gray70")) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) +
  scale_y_continuous(limits=c(0, 16)) +
  ylab("% Reads with correct V-, J-Gene, \nCDR3 length of productive Reads") +
  xlab("") +
  ggtitle("Productive Ab Reads, \nfiltered for V-, J-Gene and CDR3 length") 

reads_tp_VJ_CDRlength_freq

#not needed any more
mutate(Sequence.ID=gsub("or_IGLJ._", "", Sequence.ID)) %>%
  mutate(Sequence.ID=gsub("or_IGLV.-.._", "", Sequence.ID)) %>%
  separate(Sequence.ID, c("V_gene", "J_gene", "CDR3_length", "ECBC", "x"), sep = "_") %>%
  separate(ECBC, c("ECBC", "reads"), sep = "-") %>%
  filter(reads>2) %>%
  mutate(reads=as.numeric(reads))


df <- data.frame(b=1:10, c=c("z", "a"))
ggplot(df, aes(x=1, y=b, fill=c)) + 
  geom_bar(stat="identity", position="dodge")

df$c <- factor(df$c, levels=c("z", "a"))
ggplot(df, aes(x=1, y=b, fill=c)) + 
  geom_bar(stat="identity", position="dodge")
