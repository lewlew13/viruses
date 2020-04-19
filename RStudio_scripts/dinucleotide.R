setwd("~/Downloads/NCBI_files/dat/homopoly/dinucleotide/redone")
colnames(detailed) <- c("genome", "glen", "start", "end", "motif", "orientation")
detailed$mlen <- nchar(as.character(detailed$motif))
library(tidyverse)

detailed <- detailed %>%
  mutate(abbrev = str_extract(genome, ".*\\.[0-9]"), mlen = nchar(as.character(detailed$motif))) %>%
  select(abbrev, everything())

#dinucleotide only
detailed <- detailed %>%
  arrange(abbrev, start) %>%
  mutate(lagged = start - lag(start, n=1)) %>%
  filter(lagged > 1 | is.na(lagged))

View(table(detailed$mlen))

pdf("histogram.pdf", width = 60, height = 60)
ggplot(detailed, aes(x=mlen)) +
  geom_histogram(binwidth = 1, fill="pink", boundary=0, colour="black") +
  scale_x_continuous( breaks= 3:10, limits=c(3,10)) +
  facet_wrap(abbrev~., scales="free_y")
dev.off()


ccrep <- table(detailed$abbrev, detailed$mlen)
ccrep <- as.data.frame(ccrep)
colnames(ccrep) <- c("Accession", "MotLen", "Freq")
View(ccrep)
View(filter(ccrep, Freq>0))

genome_data <- detailed %>% select(abbrev, glen) %>% unique()
View(genome_data)

ccrep <- left_join(ccrep, genome_data, by = c("Accession" = "abbrev")) %>%
  mutate(rep_density = 1000 * Freq/ (glen/2))
View(ccrep)

ccrep$Sqrt_Rep_Density <- sqrt(ccrep$rep_density)

pdf("heatmapnew.pdf", width = 50, height = 20)
ggplot(ccrep, aes(x=Accession, y=MotLen, fill = Sqrt_Rep_Density)) +
  geom_tile() +
  xlab(label = "Accession") +
  scale_fill_gradient(name = "Sqrt Rep Density", low = "#DEFEFE", high = "#D000C7") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y = element_blank())
dev.off()
