library(dplyr)
library(tidyr)
library(seqinr)
library(ggplot2)

# read files (Nt sequences and list of names of mapped sequences)
data = read.delim("/Volumes/data/AbX/ECBC/ECBC21/IMGT_grouped/5_IgL/3_Nt-sequences.txt", header=T, sep = "\t", fill = TRUE)

mapped_names = read.table("/Volumes/data/AbX/ECBC/ECBC21/results/5_IgL_CDR3_mapped85.txt", header=F, col.names = "Sequence.ID")

# prepare Sequence_ID: select V.D.J.REGION for HCs, V.J.REGION for LCs
mapped_seq = left_join(mapped_names, data) %>%
  select(Sequence.ID, V.J.REGION) %>%
  separate(Sequence.ID, c("reads", "HV", "HJ", "CDR3", "ECBC", "junk"), sep="_") %>%
  mutate(reads=as.numeric(reads))

# combine identical sequences, keep total numbers of reads and numbers of ECBC per sequence in Sequence_ID, group V.D.J.REGION for HCs, V.J.REGION for LCs
comb_seq = mapped_seq %>%
  group_by(V.J.REGION) %>%
  mutate(reads=as.numeric(reads)) %>%
  summarise(n_reads=sum(reads), n_ECBC=n()) %>%
  mutate(subtype="5L") %>%
  unite(Sequence_ID, subtype, n_ECBC, n_reads)

for_fig = comb_seq %>%
  group_by(n_reads, n_ECBC) %>%
  summarise(Count = n())

write.fasta(as.list(comb_seq$V.J.REGION), comb_seq$Sequence_ID, "/Volumes/data/AbX/ECBC/ECBC21/results/5_IgL_CDR3_mapped85_comb.fasta")

#gsub IGHJ6 only necessary for data sets with "or" in J gene, does nothing if no "or" option: 
data_all = data %>%
  mutate(Sequence.ID=gsub("_or_IGHJ6", "", Sequence.ID)) %>%
  separate(Sequence.ID, c("HV", "HJ", "CDR3", "ECBC", "junk"), sep="_") %>%
  separate(ECBC, c("ECBC", "reads"), sep="-") %>%
  mutate(reads=as.numeric(reads))

#this is version for LCs
data_all = data %>%
  mutate(Sequence.ID=gsub("_or_IGLJ3|_or_IGLJ6|_or_IGLJ7", "", Sequence.ID)) %>%
  separate(Sequence.ID, c("LV", "LJ", "CDR3", "ECBC", "junk"), sep="_") %>%
  separate(ECBC, c("ECBC", "reads"), sep="-") %>%
  filter(reads!="NA") %>%
  mutate(reads=as.numeric(reads))

# figure: number of ECBCs vs number of reads per combined seq
figure = ggplot(for_fig, aes(x=n_ECBC, y=n_reads)) +
  geom_point(aes(size=Count), alpha=0.5, color = "darkblue") +
  scale_size(range = c(3, 7)) +
  theme(
    axis.text = element_text(size = 18, color="black"),
    panel.grid.major = element_line(colour = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white", color="gray70"),
    legend.key = element_rect(fill = "white"),
    text = element_text(size=20)
  ) +
  xlab("\nNumber of Barcodes per\nidentical consensus sequence") +
  ylab("Number of identical reads per\nidentical consensus sequence\n") +
  #scale_x_continuous(limits=c(0,13)) +
  #scale_y_continuous(limits=c(0,56)) +
  ggtitle("Combined Sequences (tp2, IgG1)") 

figure

ggsave(filename="/Volumes/data/AbX/ECBC/results/figures/comb_Seqs.pdf", plot=figure, width = 30/2.54, height = 21/2.54)

# figure: histogram of reads per unique ECBCs
d = ggplot(data_all, aes(x=reads)) +
  geom_histogram(binwidth=0.5, origin=-0.25, fill = "darkblue", alpha=0.7) +
  theme_bw() +
  xlab("\n Reads per Barcode") +
  ylab("Count") +
  scale_x_continuous(limits=c(2,8), breaks=c(3,10,20,30,40,50), labels=c(3,10,20,30,40,50), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,350)) +
  ggtitle("Reads per unique Barcode (tp2, IgL)") 

d

ggsave(filename="/Volumes/data/AbX/ECBC/ECBC21/results/figures/2_IgL_Hist_reads_per_ECBC.pdf", plot=d, width = 30/2.54, height = 21/2.54)


# read files (Nt sequences and list of names of mapped sequences)
mapped_names = read.table("/Volumes/data/AbX/ECBC/IMGT_grouped/1_IgA_Ymotivs.txt", header=F, fill = TRUE, sep = "\t")

mapped_seq = mapped_names %>%
  mutate(Sequence.ID=gsub("Homsap-", "", V2)) %>%
  mutate(V.D.J.REGION=V29) %>%
  select(Sequence.ID, V.D.J.REGION) %>%
  separate(Sequence.ID, c("HV", "HJ", "CDR3", "ECBC", "junk"), sep="_") %>%
  separate(ECBC, c("ECBC", "reads"), sep="-") %>%
  mutate(reads=as.numeric(reads))

comb_seq = mapped_seq %>%
  group_by(V.D.J.REGION) %>%
  summarise(n_reads=sum(reads), n_ECBC=n()) %>%
  mutate(subtype="1A") %>%
  unite(Sequence_ID, subtype, n_ECBC, n_reads)

#this contains V.D.J.REGION of sequences with >1 ECBC
extr_ECBC = comb_seq %>%
  separate(Sequence_ID, c("subtype", "ECBC", "reads"), sep="_") %>%
  filter(ECBC>1) %>%
  select(V.D.J.REGION)

dd = left_join(extr_ECBC, mapped_seq) %>%
  group_by(V.D.J.REGION)

#combine CDR3 mapped and Y motif filtered sequences
data_CDR3 = read.delim("/Volumes/data/AbX/ECBC/ECBC21/results/all_mapped85_comb.fasta", header=F, sep = "\t", fill = TRUE)

CDR3_names = data_CDR3 %>%
  filter(grepl(">", V1)) %>%
  rename(ID = V1) %>%
  mutate(ID=gsub(">", "", ID))

CDR3_seq = data_CDR3 %>%
  filter(!grepl(">", V1)) %>%
  rename(nt = V1) 

CDR3 = bind_cols(CDR3_names, CDR3_seq) %>%
  mutate(selection="CDR")

data_Y = read.delim("/Volumes/data/AbX/ECBC/ECBC21/results/all_Ymotifs_comb.fasta", header=F, sep = "\t", fill = TRUE)

Y_names = data_Y %>%
  filter(grepl(">", V1)) %>%
  rename(ID = V1) %>%
  mutate(ID=gsub(">", "", ID))

Y_seq = data_Y %>%
  filter(!grepl(">", V1)) %>%
  rename(nt = V1)

Y = bind_cols(Y_names, Y_seq) %>%
  mutate(selection2="Y")


data_comb = full_join(CDR3, Y, by=c("nt", "ID")) %>%
  unite(ID_comb, ID, selection, selection2) %>%
  mutate(ID_comb=gsub("NA", "", ID_comb))

write.fasta(as.list(data_comb$nt), data_comb$ID_comb, "/Volumes/data/AbX/ECBC/ECBC21/results/all_mapped.fasta")



mutate(Sequence.ID=gsub("Homsap-", "", Sequence.ID)) %>%
