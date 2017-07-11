library(dplyr)
library(tidyr)
library(seqinr)
library(ggplot2)

# this writes fasta files of CDR3 and full-length VDJ seqeunces (of IMGT summary files, uploaded after error correction)
data = read.delim("/Volumes/data/AbX/ECBC/results/script170529/IMGT_grouped/2_IgG1/3_Nt-sequences.txt", header=T, sep = "\t", fill = TRUE)
dataG = read.delim("/Volumes/data/AbX/ECBC/IMGT_grouped/S2_IgG1/1_Summary.txt", header=T, sep = "\t", fill = TRUE)
#dataA = read.delim("/Volumes/data/AbX/ECBC/IMGT_grouped/S2_IgA/1_Summary.txt", header=T, sep = "\t", fill = TRUE)


write.fasta(as.list(data$CDR3.IMGT), data$Sequence.ID, "/Volumes/data/AbX/ECBC/results/script170529/2_IgG1_ntCDR3.fasta")

write.fasta(as.list(data$V.D.J.REGION), data$Sequence.ID, "/Volumes/data/AbX/ECBC/IMGT_grouped/S5_IgG3/S5_IgG3_nt.fasta")

# this writes fasta files of mapped sequences after bwa sam
mapped_names = read.table("/Volumes/data/AbX/ECBC/IMGT_grouped/S2_IgA/2_IgA_CDR3_mapped85.txt", header=F, col.names="Sequence.ID")
#mapped_namesG = read.table("/Volumes/data/AbX/ECBC/IMGT_grouped/S2_IgG1/2_IgG1_CDR3_mapped85.txt", header=F, col.names="Sequence.ID")

mapped_seq = left_join(mapped_namesG, dataG) %>%
  mutate(Sequence.ID=gsub("Homsap-", "", Sequence.ID)) %>%
  mutate(subtype="S2_IgG1") %>%
  unite(Sequence_ID, subtype, Sequence.ID)

mapped_IMGT_G = left_join(mapped_namesA, dataA) %>%
  mutate(Sequence.ID=gsub("Homsap-", "", Sequence.ID)) %>%
  mutate(subtype="S2_IgA") %>%
  mutate(subtype2="S2_IgA") %>%
  unite(Sequence_ID, subtype, Sequence.ID) %>%
  arrange(V.REGION.identity..) %>%
  filter(V.REGION.identity.. <86.1) %>%
  select(Sequence_ID)

mapped_IMGT_A = left_join(mapped_namesA, dataA) %>%
  mutate(Sequence.ID=gsub("Homsap-", "", Sequence.ID)) %>%
  mutate(subtype="S2_IgA") %>%
  mutate(subtype2="S2_IgA") %>%
  unite(Sequence_ID, subtype, Sequence.ID) %>%
  arrange(V.REGION.identity..) 

mapped_IMGT = rbind(mapped_IMGT_A, mapped_IMGT_G)
  
write.fasta(as.list(mapped_seq$V.D.J.REGION), mapped_seq$Sequence_ID, "/Volumes/data/AbX/ECBC/results/S2_IgG1_CDR3_mapped85_full.fasta")
write.table(mapped_IMGT_G, "/Volumes/data/AbX/ECBC/results/2_IgA_CDR3_mapped85_IMGT_mut.txt", row.names = F, col.names=F, quote = F)

figure = ggplot(mapped_IMGT_G, aes(x=V.REGION.identity.., fill=subtype2, color=subtype2)) +
  geom_histogram(alpha=0.7) +
  scale_fill_manual(values=c("steelblue")) +
  scale_color_manual(values=c("steelblue")) +
  geom_vline(xintercept = 86.11, color = "red") +
  xlab(" \nV Region identity (%)") +
  ylab("Count\n") +
  ggtitle("V Region identity distribution") +
  scale_x_continuous(limits=c(82,94), breaks=seq(84, 92, 2)) +
  scale_y_continuous(limits=c(0,30), expand = c(0, 0)) +
  theme(plot.title=element_text(size=24), text = element_text(size=24), legend.title=element_blank(),
        axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey"), legend.text = element_text(size = 18)) 

figure

ggsave(filename="/Volumes/data/AbX/ECBC/results/figures/S5_VRegionIdent_CDR385.pdf", plot=figure, width = 30/2.54, height = 21/2.54)

# this gets fastas of y-motiv mapped sequences
data = read.delim("/Volumes/data/AbX/ECBC/IMGT_grouped/5_IgG1_Ymotivs.txt", header=F, sep = "\t", fill = TRUE)

data1 = data %>%
  select(V2, V29) %>%
  mutate(Sequence_ID=gsub("Homsap-", "", V2)) %>%
  separate(Sequence_ID, c("HV", "HJ", "CDR3", "ECBC", "junk"), sep="_") %>%
  separate(ECBC, c("ECBC", "reads"), sep="-") %>%
  group_by(V29) %>%
  mutate(reads=as.numeric(reads)) %>%
  summarise(n_reads=sum(reads), n_ECBC=n()) %>%
  mutate(subtype="5G1") %>%
  unite(Sequence_ID, subtype, n_ECBC, n_reads)

write.fasta(as.list(data1$V29), data1$Sequence_ID, "/Volumes/data/AbX/ECBC/results/S5_IgG1_Ymotivs_comb.fasta")
