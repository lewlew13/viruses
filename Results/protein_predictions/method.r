library(tidyverse)
library(readr)
seq144 <- read_csv("144.fasta_capsids.scores.fixed.csv", col_types="cdf")
View(seq144)
ggplot(seq144, aes(x=Repeats, y=Score)) + geom_boxplot()
ggplot(seq144, aes(x=Score)) + geom_histogram(binwidth = 0.5, fill="lightblue", colour="black")

t.test(Score ~ Repeats,seq144_t, alternative = "less")

p.adjust(p = c(0.2354,0.003594,0.8143,0.5224,0.5269,0.1451,0.6613,0.79994,0.037,0.0003078,0.912,9.94E-05), method = "fdr")