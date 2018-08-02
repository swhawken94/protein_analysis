# 07.31.2018
# Creator: Susana Wilson Hawken

#-------------------------------------------------------------------#

#      analysis of cancer mutations in protein activation domain    #

#-------------------------------------------------------------------#

#download packages
library(dplyr)
library(data.table)
library(protr)
library(factoextra)
library(cluster) 
library(gridExtra)
library(plyr)
library(readr)
library(gridExtra)
library(grid)

#-------------------------------------------------------------------#

#                        user-defined variables                     #

#-------------------------------------------------------------------#

# enter peptide ID
peptide.id <- "ENSP00000367207"
uniprot.id <- "P01106"

# transactivation domain coordinates
TAD_start <- 1
TAD_end <- 143

# DBD coordinates
DBD_start <- 265
DBD_end <- 318

# IDR coordinates
IDR1 <- c(49,57)
IDR2 <- c(107,108)
IDR3 <- c(167,185)
IDR4 <- c(214, 215)
IDR5 <- c(218,313)
IDR6 <- c(318,340)
IDR7 <- c(342,343)
IDR8 <- c(350,384)
IDR9 <- c(395,396)


#-------------------------------------------------------------------#

#                        unchanged variables                        #

#-------------------------------------------------------------------#

# create url for protein of interest from peptide id and uniprot url
url.firsthalf <- "https://www.uniprot.org/uniprot/"
url.secondhalf <- ".fasta"
url <- paste(url.firsthalf, uniprot.id, url.secondhalf, sep = "", collapse = "")

# columns for data frame
numcol.sampleid = 1
numcol.wildtype = 2
numcol.position = 3
numcol.altered = 4

#-------------------------------------------------------------------#

#                        define functions                           #

#-------------------------------------------------------------------#

get_wildtype <- function(y) {
  x <- c()
  for (i in 1:length(y)){
    start <- substr(as.character(y[i]),0,1)
    x <- c(x,start)
  }
  return(x)
}

get_altered <- function(y) {
  x <- c()
  for(i in 1:length(y)){
    end <- substr(as.character(y[i]),length(strsplit(as.character(y[i]), "")[[1]]), 
                  length(strsplit(as.character(y[i]), "")[[1]]))
    x <- c(x, end)
  }
  return(x)
}

get_position <- function(y) {
  x <- c()
  for(i in 1:length(y)){
    end <- substr(as.character(y[i]),2, 
                  length(strsplit(as.character(y[i]), "")[[1]])-1)
    x <- c(x, end)
  }
  return(x)
}

#-------------------------------------------------------------------#

#                        load and clean files                       #

#-------------------------------------------------------------------#

# set path and read the file
setwd("~/Desktop/Young_Lab_2018/myc_Analysis/")
args = commandArgs(trailingOnly=TRUE)

# read command line arguments
if (length(args)==0) {
     stop("At least two argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
     stop("At least two arguments must be supplied (input file).n", call.=FALSE)
} else {
  args[3] = "out.txt"
}

# read amino acid statistics
# TK change back "AA_stats.csv" to arg[2]
aa.stats <- read.csv("AA_stats.csv", sep = "\t", col.names = c("peptide.id", "length",
                                                                       "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                                                                       "perc.A", "perc.C", "perc.D", "perc.E", "perc.F", "perc.G", "perc.H", "perc.I", "perc.K", "perc.L",
                                                                       "perc.M", "perc.N", "perc.P", "perc.Q", "perc.R", "perc.S", "perc.T", "perc.V", "perc.W", "perc.Y",
                                                                       "aliphatic", "perc.aliphatic", "aromatic", "perc.aromatic", "acidic", "perc.acidic", "basic", 
                                                                       "perc.basic", "hydrophobic", "perc.hydrophobic", "polar", "perc.polar", "nonpolar", "perc.nonpolar", 
                                                                       "essential", "perc.essential"))

# read in the cBioportal mutation data
# TK change back "myc_mutations.tsv to args[1]
mut.df = read.csv("myc_mutations.tsv", sep = "\t", col.names = c("study",	"sample.id", "cancer.type",	"protein.change",	"annotation",	"functional.impact",	"mutation.type",	"copy.num",
                                                 "cosmic",	"MS",	"VS",	"center",	"chromosome",	"start.pos",	"end.pos",	"ref", "var",	"allele.freq.T",	"allele.freq.N",	
                                                 "variant.reads",	"ref.reads",	"variant.reads.N",	"ref.reads.N",	"num.mut.sample"))


# capture amino acid stats for protein of interest
protein.stats <- aa.stats[which(as.character(aa.stats$peptide.id) == peptide.id),2:ncol(aa.stats)]
protein.stats <- protein.stats %>%
  as.character() %>%
  as.numeric() %>%
  unlist()

# protein length in first column
protein.len <- protein.stats[1]

# amino acid counts in next 20 columns
whole.counts <- protein.stats[2:21]

# amino acid frequencies (calculate percentages)
whole.freq <- 100 *(whole.counts/protein.len)

# amino acid letters
whole.names <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", 
                 "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

whole.df <- data.frame(whole.names, whole.freq)

# cleaning data frame
# capture only the amino acid change and position
mut.df <- filter(mut.df, mutation.type == "Missense_Mutation")

mutants <- mut.df$protein.change

# get wildtype aa
wildtype <- get_wildtype(mutants)
altered <- get_altered(mutants)
position <- get_position(mutants)

# create dataframe with wiltype aa, position of aa, and what it's altered to
sample.id <- mut.df$sample.id
mut.dt <- data.frame(sample.id,wildtype, position, altered)

#-------------------------------------------------------------------#

#    Plot the amino acid frequencies in the protein TAD, DBD,	      #
#          	      and whole sequence			                          #

#-------------------------------------------------------------------#

# downlaod the protein sequence from uniprot and delete the 
# information line (very first line) in the file

download.file(url, destfile = "sequence.fasta", method = "curl")
file_name <- "sequence.fasta"
sequence <- read.csv(file_name, skip=1)
write.table(sequence, "sequence.fasta", quote= FALSE, row.names = FALSE)

# break up TAD sequence into amino acid frequencies
# eliminate spaces and new line characters
sequence <- readChar(file_name, file.info(file_name)$size)
sequence <- sequence %>% 
  gsub("[\r\n]", "",.) %>% 
  gsub("\\s", "", .)

strsplit_cleaned <- strsplit(sequence, "")[[1]]


# frequencies of amino acids in TADs and DBDs 
sequence_unlist <- unlist(strsplit_cleaned)

# test length of protein is correct
print("length of protein")
print(protein.len)
print("length from uniprot")
print(length(sequence_unlist))

# TAD sequence
TAD <- sequence %>% substr(TAD_start, TAD_end)

# composition of TADs and DBDs
TAD.freq <- 100*extractAAC(TAD)
TAD.names <- names(TAD.freq)

# dataframe for TAD
TAD.df <- data.frame(TAD.names, TAD.freq)
TAD.df <- TAD.df[order(TAD.df$TAD.names),]


# max frequency
whole.max <- max(whole.freq)
TAD.max <- max(TAD.freq)
max <- max(whole.max, TAD.max)

# ggplot representation
whole.plot <- ggplot(data = whole.df, aes(x = whole.names, y = whole.freq)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  xlab("Amino Acids") + 
  ylab("Amino acid frequency (%)") + 
  ylim(0,max)

TAD.plot <- ggplot(data = TAD.df, aes(x = TAD.names, y = TAD.freq)) + 
  geom_bar(stat = "identity", fill = "green") + xlab("Amino Acids") +
  ylab("Amino acid frequency (%)")  +
  ylim(0,max) 

# plot just the TAD and whole sequence
#grid.arrange(whole.plot, TAD.plot,nrow = 2)


#-------------------------------------------------------------------#

#      analysis of amino acid relative fold change                  #

#-------------------------------------------------------------------#

TAD.df <- TAD.df %>% mutate(.,foldchange = ifelse(TAD.freq >whole.df$whole.freq, 
                                                  TAD.freq/whole.df$whole.freq, 
                                                  ifelse((TAD.freq < whole.df$whole.freq) & (TAD.freq != 0), 
                                                         -(whole.df$whole.freq/TAD.freq), 
                                                         0)))

# ggplot 
aa.foldchange <- ggplot(data = TAD.df, aes(x = TAD.names, y = foldchange)) + 
  geom_bar(stat = "identity", fill = "purple") + 
  xlab("Amino Acids") + 
  ylab("Relative fold change")

#grid.arrange(whole.plot, TAD.plot, aa.foldchange, nrow = 3)


#-------------------------------------------------------------------#

#		analysis of mutation frequency	                                #

#-------------------------------------------------------------------#

whole.mutants <- paste(as.character(mut.dt$wildtype), sep = "", collapse = "")

TAD.mutants <- mut.dt[which(as.numeric(as.character(mut.dt$position)) < TAD_end),]
TAD.mutants <- paste(TAD.mutants[which(as.numeric(as.character(TAD.mutants$position)) > TAD_start),
                                 numcol.wildtype], 
                     sep  = "", 
                     collapse = "")

print(TAD.mutants)

whole.mut.freq <- 100*extractAAC(whole.mutants)
TAD.mut.freq <- 100*extractAAC(TAD.mutants)

TAD.mut.names<- names(TAD.mut.freq)
whole.mut.names <- names(whole.mut.freq)


TAD.mut.df <- data.frame(TAD.mut.freq, TAD.mut.names)
whole.mut.df <- data.frame(whole.mut.freq, whole.mut.names)

TAD.mut.df <- TAD.mut.df[order(TAD.mut.df$TAD.mut.names),]
whole.mut.df <- whole.mut.df[order(whole.mut.df$whole.mut.names),]

TAD.df$mut.freq <- TAD.mut.df$TAD.mut.freq
whole.df$mut.freq <- whole.mut.df$whole.mut.freq

whole.max <- max(whole.df$mut.freq)
TAD.max <- max(TAD.df$mut.freq)
max <- max(whole.max, TAD.max)

whole.mut.plot <- ggplot(data = whole.df, aes(x = whole.names, y = mut.freq)) + 
  geom_bar(stat = "identity", fill = "blue")  + 
  xlab("Amino Acids") +
  ylab("Mutation frequency (%)") + 
  ylim(0,max)

TAD.mut.plot <- ggplot(data = TAD.df, aes(x = TAD.names, y = mut.freq)) + 
  geom_bar(stat = "identity", fill = "green") + xlab("Amino Acids") + 
  ylab("Mutation frequency (%)") +
  ylim(0,max)

#grid.arrange(whole.mut.plot,  TAD.mut.plot, nrow = 2)

#-------------------------------------------------------------------#

#               analysis of mutation relative fold change                      #

#-------------------------------------------------------------------#

TAD.df <- TAD.df %>% mutate(.,foldchange.mut = ifelse(mut.freq > whole.df$mut.freq, 
                                                      mut.freq/whole.df$mut.freq, 
                                                      ifelse((mut.freq < whole.df$mut.freq) & (mut.freq != 0), 
                                                             -(whole.df$mut.freq/mut.freq), 
                                                             0)))

mut.foldchange <- ggplot(data = TAD.df, aes(x = TAD.names, y = foldchange.mut)) + 
  geom_bar(stat = "identity", fill = "purple") + 
  xlab("Amino Acids") +
  ylab("Relative fold change") 

#grid.arrange(whole.mut.plot, TAD.mut.plot, mut.foldchange, nrow = 2)

#-------------------------------------------------------------------#

#                specific amino acid analysis                       #

#-------------------------------------------------------------------#

# number of annotated mutants in TAD
TAD.muts = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% nrow()


#-------------------------------------------------------------------#

#                proline amino acid analysis                        #

#-------------------------------------------------------------------#

# number of annotated proline mutants in TAD
TAD.pro.muts = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "P") %>%
  nrow(.)

#-------------------------------------------------------------------#

#               acidic amino acid analysis                          #

#-------------------------------------------------------------------#
# number of annotated glutamic acid mutants in TAD
TAD.glu.muts = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "E") %>%
  nrow(.)

# number of annotated aspartic acid mutants in TAD
TAD.asp.muts = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "D") %>%
  nrow(.)

TAD.acid.muts = TAD.glu.muts + TAD.asp.muts

#-------------------------------------------------------------------#

#                aromatic amino acid analysis                       #

#-------------------------------------------------------------------#

TAD.phe.muts = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "F") %>%
  nrow()

TAD.tyr.muts = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "Y") %>%
  nrow()

TAD.trp.muts = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "W") %>%
  nrow()

TAD.aro.muts = TAD.phe.muts + TAD.tyr.muts + TAD.trp.muts

#-------------------------------------------------------------------#

#               most frequently mutated proline (TAD)               #

#-------------------------------------------------------------------#

g1 <- mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "P") %>%
  ddply(., .(wildtype,position, altered), nrow) %>%
  .[order(.$V1, decreasing = TRUE),] %>%
  .[1:5,] %>% 
  dplyr::rename(., Number.annotated.mutations = V1)
rownames(g1) <- c("1","2","3","4", "5")

#-------------------------------------------------------------------#

#               how many prolines mutated in (TAD)                  #

#-------------------------------------------------------------------#
num.pro.TAD = which(sequence_unlist[TAD_start:TAD_end] == "P") %>% 
  length()


##
num.pro.mut.TAD = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "P") %>%
  ddply(., .(position), nrow) %>%
  nrow()

# make dataframe
p1 <- data.frame(num.pro.TAD, num.pro.mut.TAD)
colnames(p1) <- c("Number prolines in TAD", "Number proline mutated in TAD")
#-------------------------------------------------------------------#

#               most frequently mutated acidic (TAD)                #

#-------------------------------------------------------------------#
cat("Most frequently mutated acidics (TAD)\n")
g2 <- mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "E") %>%
  ddply(., .(wildtype,position, altered), nrow) %>%
  .[order(.$V1, decreasing = TRUE),] %>%
  .[1:5,] %>%
  dplyr::rename(., Number.annotated.mutations = V1)
rownames(g2) <- c("1","2","3","4", "5")
g3 <- mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "D") %>%
  ddply(., .(wildtype,position, altered), nrow) %>%
  .[order(.$V1, decreasing = TRUE),] %>%
  .[1:5,] %>%
  dplyr::rename(., Number.annotated.mutations = V1)
rownames(g3) <- c("1","2","3","4", "5")
#-------------------------------------------------------------------#

#               how many acidics mutated in (TAD)                   #

#-------------------------------------------------------------------#
num.glu.TAD = which(sequence_unlist[TAD_start:TAD_end] == "E") %>% 
  length()

num.asp.TAD = which(sequence_unlist[TAD_start:TAD_end] == "D") %>% 
  length()
cat("Number acidics in TAD\n")
num.acid.TAD <- (num.glu.TAD + num.asp.TAD)

##
num.glu.mut.TAD = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "E") %>%
  ddply(., .(position), nrow) %>%
  nrow()

num.asp.mut.TAD = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "D") %>%
  ddply(., .(position), nrow) %>%
  nrow()
cat("Number acidics mutated in TAD\n")
num.acid.mut.TAD <- (num.glu.mut.TAD + num.asp.mut.TAD)

# make dataframe
col1 <- c(num.acid.TAD)
col2 <- c(num.acid.mut.TAD)
p2 <- data.frame(col1,col2)
colnames(p2) <- c("Number acidics in TAD", "Number acidics mutated in TAD")
#-------------------------------------------------------------------#

#               most frequently mutated aromatic (TAD)              #

#-------------------------------------------------------------------#
cat("Most frequently mutated aromatics (TAD)\n")
g4 <- mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "F") %>%
  ddply(., .(wildtype,position, altered), nrow) %>%
  .[order(.$V1, decreasing = TRUE),] %>%
  .[1:5,] %>%
  dplyr::rename(., Number.annotated.mutations = V1)

rownames(g4) <- c("1","2","3","4", "5")

g5 <- mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "Y") %>%
  ddply(., .(wildtype,position, altered), nrow) %>%
  .[order(.$V1, decreasing = TRUE),] %>%
  .[1:5,] %>%
  dplyr::rename(., Number.annotated.mutations = V1)

rownames(g5) <- c("1","2","3","4", "5")

mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "W") %>%
  ddply(., .(wildtype,position, altered), nrow) %>%
  .[order(.$V1, decreasing = TRUE),] %>%
  .[1:5,]

#-------------------------------------------------------------------#

#               how many aromatics mutated in (TAD)                #

#-------------------------------------------------------------------#
num.trp.TAD = which(sequence_unlist[TAD_start:TAD_end] == "W") %>% 
  length()

num.tyr.TAD = which(sequence_unlist[TAD_start:TAD_end] == "Y") %>% 
  length()

num.phe.TAD = which(sequence_unlist[TAD_start:TAD_end] == "F") %>% 
  length()
cat("Number aromatics in TAD\n")
num.aro.TAD <- (num.trp.TAD + num.tyr.TAD + num.phe.TAD)

##
num.trp.mut.TAD = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "W") %>%
  ddply(., .(position), nrow) %>%
  nrow()

num.tyr.mut.TAD = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "Y") %>%
  ddply(., .(position), nrow) %>%
  nrow()

num.phe.mut.TAD = mut.dt %>% subset(position %>% as.character() %>% as.numeric() < TAD_end) %>% 
  subset(wildtype %>% as.character() == "F") %>%
  ddply(., .(position), nrow) %>%
  nrow()
cat("Number aromatics mutated in TAD\n")
num.aro.mut.TAD <- (num.trp.mut.TAD + num.tyr.mut.TAD + num.phe.mut.TAD)

col1 <- c(num.aro.TAD)
col2 <- c(num.aro.mut.TAD)
p3 <- data.frame(col1,col2)
colnames(p3) <- c("Number aromatics in TAD", "Number aromatics mutated in TAD")
#-------------------------------------------------------------------#

#                putting together the data (TAD)                    #

#-------------------------------------------------------------------#

TAD.pro <- TAD.df %>% subset(.,TAD.names == "P") %>% .$TAD.freq
TAD.glu <- TAD.df %>% subset(.,TAD.names == "E") %>% .$TAD.freq
TAD.asp <- TAD.df %>% subset(.,TAD.names == "D") %>% .$TAD.freq
TAD.phe <- TAD.df %>% subset(.,TAD.names == "F") %>% .$TAD.freq
TAD.trp <- TAD.df %>% subset(.,TAD.names == "W") %>% .$TAD.freq
TAD.tyr <- TAD.df %>% subset(.,TAD.names == "Y") %>% .$TAD.freq

perc.pro = as.integer(TAD.pro) 
perc.acid = as.integer(TAD.glu + TAD.asp)
perc.aro = as.integer(TAD.phe + TAD.tyr + TAD.trp)

perc.mut.pro <- as.integer(TAD.pro.muts/TAD.muts*100)
perc.mut.acid <- as.integer(TAD.acid.muts/TAD.muts*100)
perc.mut.aro <- as.integer(TAD.aro.muts/TAD.muts*100)

acidics <- c(perc.acid, perc.mut.acid)
prolines <- c(perc.pro, perc.mut.pro)
aromatics <- c(perc.aro, perc.mut.aro)

df <- data.frame()
df <- rbind(df, acidics)
df <- rbind(df, prolines)
df <- rbind(df, aromatics)

colnames(df) <- c("percentage AA in TAD", "percentage mutation in TAD")
rownames(df) <- c("acidic", "proline", "aromatic")

g6 <- df

title1 <- textGrob("Most frequently mutated prolines in TAD")

pdf("TAD.mutation.stats.pdf")
grid.arrange(tableGrob(g6), ncol =1 , nrow = 1, top = textGrob("Percentage amino acids and mutations in TAD"))
  grid.arrange(tableGrob(g1),tableGrob(p1), ncol =1 , nrow = 2, top = textGrob("Most frequently mutated prolines in TAD"))
  grid.arrange(tableGrob(g2), tableGrob(p2), ncol =1 , nrow = 2, top = textGrob("Most frequently mutated glutamic acids in TAD"))
  grid.arrange(tableGrob(g3), ncol =1 , nrow = 1, top = textGrob("Most frequently mutated aspartic acids in TAD"))
  grid.arrange(tableGrob(g4), ncol =1 , nrow = 1, top = textGrob("Most frequently mutated phenylalanines in TAD"))
  grid.arrange(tableGrob(g5), tableGrob(p3), ncol =1 , nrow = 2, top = textGrob("Most frequently mutated tyrosines in TAD"))
dev.off()

#-------------------------------------------------------------------#

#                putting together the data (whole)                    #

#-------------------------------------------------------------------#

# all number of mutants
num.muts = mut.dt  %>% nrow()

whole.pro <- whole.df %>% subset(.,whole.names == "P") %>% .$whole.freq
perc.whole.pro <- whole.df %>% subset(., whole.names == "P") %>% .$mut.freq

whole.acid <- (whole.df %>% subset(.,whole.names == "E") %>% .$whole.freq) +
  (whole.df %>% subset(.,whole.names == "D") %>% .$whole.freq)

perc.whole.acid <- (whole.df %>% subset(., whole.names == "E") %>% .$mut.freq) + 
  (whole.df %>% subset(., whole.names == "D") %>% .$mut.freq)

whole.aro <- (whole.df %>% subset(.,whole.names == "W") %>% .$whole.freq) + 
  (whole.df %>% subset(.,whole.names == "F") %>% .$whole.freq) +
  (whole.df %>% subset(.,whole.names == "Y") %>% .$whole.freq) 

perc.whole.aro <- (whole.df %>% subset(., whole.names == "W") %>% .$mut.freq) + 
  (whole.df %>% subset(., whole.names == "F") %>% .$mut.freq) +
  (whole.df %>% subset(., whole.names == "Y") %>% .$mut.freq) 

col1 <- c( ceiling(whole.acid), ceiling(whole.pro),ceiling(whole.aro))
col2 <- c( ceiling(perc.whole.acid), ceiling(perc.whole.pro),ceiling(perc.whole.aro))

g7 <- data.frame(col1, col2)
colnames(g7) = c("Percentage AA", "Percentage mutations")
rownames(g7) = c("acidic", "proline", "aromatic")

pdf("TAD.mutation.stats.pdf")
  grid.arrange(tableGrob(g7), ncol =1 , nrow = 2, top = textGrob("Percentage amino acids and mutations in whole sequence"))
  grid.arrange(tableGrob(g6), ncol =1 , nrow = 2, top = textGrob("Percentage amino acids and mutations in TAD"))
  grid.arrange(tableGrob(g1),tableGrob(p1), ncol =1 , nrow = 2, top = textGrob("Most frequent proline mutations in TAD"))
  grid.arrange(tableGrob(g2), tableGrob(p2), ncol =1 , nrow = 2, top = textGrob("Most frequent glutamic acid mutations in TAD"))
  grid.arrange(tableGrob(g3), ncol =1 , nrow = 2, top = textGrob("Most frequent aspartic acid mutations in TAD"))
  grid.arrange(tableGrob(g4), tableGrob(p3),ncol =1 , nrow = 2, top = textGrob("Most frequent phenylalanine mutations in TAD"))
  grid.arrange(tableGrob(g5),  ncol =1 , nrow = 2, top = textGrob("Most frequent tyrosine mutations in TAD"))
dev.off()
