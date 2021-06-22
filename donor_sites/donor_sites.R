### This file will blast sequences for donor site segments (for use with internal repeat logic), 
### to be fed into find_ds_combinations.py.

list.of.packages <- c("ggplot2","stringr","Biostrings","dplyr","phylotools","data.table","cowplot")
suppressMessages(invisible(lapply(list.of.packages,library,character.only=T)))

setwd("/Users/uwvirongs/Documents/Michelle/longitudinal_tprk/redo_longitudinal_tprk_with_doubleprep/validating_dna_reads")

# Read in all of the validated over5*.csv files to get the nucleotide sequences of the variable regions
filenames <- Sys.glob("*validated_over5*")

# Convert csv tables into fasta format
for (i in c(1:length(filenames))){
  seq <- NULL
  sample<-as.character((strsplit(filenames[i],"_")[[1]][1]))
  print(sample)
  mycsv <- read.table(filenames[i],sep=',',header = TRUE)
  colnames(mycsv) <- c("Region","Read","RelativeFreq","Count")
  newnames <- paste(sample,mycsv$Region,rownames(mycsv),sep='_')
  mycsv <- mutate(mycsv,newnames=paste(sample,mycsv$Region,as.character(Count),rownames(mycsv),sep='_'))
  mycsv <- mycsv[mycsv$Count >= 0,]
  dna <- DNAStringSet(mycsv$Read,use.names=TRUE)
  names(dna) <- mycsv$newnames
  filename <- paste0(sample,"_variable.fasta")
  writeXStringSet(dna, filename)
}

syscommand <- "cat *_variable.fasta > all-variable.fasta"
system(syscommand)

# Filter variable regions based on prior knowledge of conserved sequences
fastaFile <- readDNAStringSet("all-variable.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
allvariable <- data.frame(seq_name, sequence)
colnames(allvariable) <- c("seq.name", "seq.text")
allvariable$region <- str_split_fixed(allvariable$seq.name,"_",4)[,2]
allvariable[allvariable$region=="V3",]$seq.text <- str_sub(allvariable[allvariable$region=="V3",]$seq.text, end=-24)
allvariable[allvariable$region=="V2",]$seq.text <- str_sub(allvariable[allvariable$region=="V2",]$seq.text, end=-15)
allvariable[allvariable$region=="V5",]$seq.text <- str_sub(allvariable[allvariable$region=="V5",]$seq.text, end=-14)
allvariable[grep("N",allvariable$seq.text),]$seq.text <- NA
allvariable <- allvariable[complete.cases(allvariable), ]
allvariable <- mutate(allvariable,seq.length = nchar(allvariable$seq.text))
dnaString <- BStringSet(allvariable$seq.text)
names(dnaString) = paste0(allvariable$seq.name)
writeXStringSet(dnaString, "all-variable_edited.fasta",append=FALSE,format="fasta")

# Blast the edited allvariable fasta against our 17.2kb locus
syscommand <- "blastn -db Nichols_tprD_new.fasta -dust no -soft_masking false -query all-variable_edited.fasta -evalue 1e1 -word_size 5 -ungapped -perc_identity 100 -penalty -10 -outfmt 6 -out donorsite_search.blast"
system(syscommand)

# Read in relevant reference files
tprd <- read_file("nichols_tprD_locus.txt")
tprd_comp <- read_file("nichols_tprD_locus_comp.txt")
internal_repeats <- read_tsv("acl_internal_repeats.txt")
#internal_repeats <- read_tsv("new_repeats_list.tsv")

hits <- read_tsv("donorsite_search.blast",col_names=c("qseqid","sseqid","pident","length","mismatch",
                                                  "gapopen","qstart","qend","sstart","send","evalue","bitscore"))


# Now we make it go the right direction relative to tprD locus numbering
rel_start <- ifelse(hits$sstart > hits$send, hits$send, hits$sstart)
rel_end <- ifelse(hits$sstart < hits$send, hits$send, hits$sstart)
hits$rel_start <- rel_start
hits$rel_end <- rel_end
x <- as.character(NULL)
for (i in 1:nrow(hits)) {
  seq <- substr(tprd, hits$rel_start[i], hits$rel_end[i])
  x <- c(x, seq)
}
hits$fwd <- x
y <- as.character(NULL)
for (i in 1:nrow(hits)) {
  seq <- substr(tprd_comp, hits$rel_start[i], hits$rel_end[i])
  seq <- reverse(seq)
  y <- c(y, seq)
}
hits$rev <- y

sequence <- ifelse(hits$sstart > hits $send, hits$rev, hits$fwd)
hits$sequence <- sequence
hits <- subset(hits, select = -c(fwd, rev)) #avoid confusion later.....
hits <- rownames_to_column(hits)
colnames(hits)[1] <- "IDnum"

# Add region column
region <- str_split_fixed(hits$qseqid,"_",4)[,2]
hits <- cbind(hits,region)


annotations <- read.table("Nichols_tprD_newannotations_edited.csv",sep=",",header=TRUE)

# Add donor site information to hits
hits2 <- as.data.frame(NULL)
for(i in 1:nrow(hits)) {
  sstart <- (hits[i,10])[[1]]
  send <- (hits[i,11])[[1]]
  min <- min(sstart,send)
  max <- max(sstart,send)
  vregion <- (hits[i,16])[[1]]
  
  annotations_subset <- filter(annotations, grepl(vregion,ds_name))
  for(ds in 1:nrow(annotations)) {
    full_sstart <- annotations[ds,1]
    full_send <- annotations[ds,2]
    full_min <- min(full_sstart, full_send)
    full_max <- max(full_sstart, full_send)
    
    ds_name <- annotations[ds,3]
    if(min >= full_min && max <= full_max) {
      add_row <- hits[i,]
      add_row$ds_name <- ds_name
      hits2 <- rbind(hits2,add_row)
    }
  }
}


##### Find largest donor site hits, not needed for internal repeat logic
### Set up the blasthits file
blasthits <- read.table("Step3_search.blast",
                        col.names = c("qseqid","sseqid","pident","length","mismatch",
                                      "gapopen","qstart","qend","sstart","send","evalue","bitscore"))

count <- as.numeric(str_split_fixed(blasthits$qseqid,"_",4)[,3])
sample <- str_split_fixed(blasthits$qseqid,"_",4)[,1]
region <- str_split_fixed(blasthits$qseqid,"_",4)[,2]
interval<-paste(blasthits$sstart,blasthits$send,sep="-")
blasthits <- dplyr::mutate(blasthits,count=count,sample=sample,interval=interval,region=region)
blasthits <- blasthits %>% dplyr::group_by(sample, region)  %>% dplyr::mutate(regionsamplecounts = sum(count))
blasthits <- ungroup(blasthits)
blasthits <- dplyr::mutate(blasthits, percentage = round(count / regionsamplecounts * 100,3))


##Filter blasthits above an arbitrary count, in this case 50 reads and within-sample percentage of 0.2%
blasthits_100 <- blasthits[blasthits$count > 50 & blasthits$percentage > 0.1,]

##Put the longest blasthits at the top to ease the sort below.
blasthits_100 <- blasthits_100[order(-blasthits_100$length),]

# Function to find the largest locus for each donor site
UniqueDonorSites <- function(df) {
  variableReg <- data.frame(df$sstart[1],df$send[1],df$region[1],df$count[1],df$sample[1],df$percentage[1])
  variableReg$morethanone <- 1
  colnames(variableReg)<- c("start","end","region","count","sample","percentage","morethanone")
  for (i in c(2:length(df$sstart)))
  {
    indf <- FALSE
    for (j in c(1:length(variableReg$start)))
    {
      if (df$region[i] != variableReg$region[j]) next()
      if (inrange(df$sstart[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE) & 
          inrange(df$send[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE)){
        indf <- TRUE
        variableReg$count[j] <- variableReg$count[j] + df$count[i]
        variableReg$percentage[j] <- variableReg$percentage[j] + df$percentage[i]
        if (df$sample[i]!=variableReg$sample[j]) variableReg$morethanone[j] <- 2
      } else if (inrange(df$sstart[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE) & 
                 inrange(df$send[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE)){
        indf <- TRUE
        variableReg$count[j] <- variableReg$count[j] + df$count[i]
        variableReg$percentage[j] <- variableReg$percentage[j] + df$percentage[i]
        if (df$sample[i]!=variableReg$sample[j]) variableReg$morethanone[j] <- 2
      } 
    }
    if (!indf) {
      counter <- counter+1
      print(c(counter,i))
      addition <- data.frame(df$sstart[i],df$send[i],df$region[i],df$count[i],df$sample[1],df$percentage[i])
      addition <- mutate(addition,morethanone=1)
      if ((inrange(df$sstart[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE) | 
           inrange(df$send[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE)) | 
          (inrange(df$sstart[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE) | 
           inrange(df$send[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE))) addition$morethanone <- 2
      variableReg <- rbind(variableReg,setNames(addition,names(variableReg)))
      indf <- FALSE
    }
  }
  return(variableReg)
}

# Get unique donor sites (can take a minute)
variableRegions <- UniqueDonorSites(blasthits_100)

### Read in the fasta locus sequence to store in variableRegions dataframe
sequence <- phylotools::read.fasta("tprDlocus.fasta")
tprDlocusSeq <- sequence$seq.text
seqname <- sequence$seq.name
variableRegions$region <- as.character(variableRegions$region)
variableRegions$region <- factor(variableRegions$region, levels = c("V1","V2","V3","V4","V5","V6","V7"))
variableRegions <- dplyr::mutate(variableRegions,length=abs(start-end)+1,
                                 seqname=seqname,source=".",feature="variable",score=".",frame=0,
                                 strand=ifelse(start-end>0,"-","+"))

##Function to get a list of sequences of the variable regions by parsing the original fasta file
getVariableSequences <- function(variableRegions){
  loci <- NULL
  for (i in c(1:length(variableRegions$start))){
    if (variableRegions$start[i] > variableRegions$end[i]) { 
      start <- (variableRegions$end[i])
      end <- (variableRegions$start[i])
      loci <- getSequence(rev(comp(getSequence(substr(tprDlocusSeq,start,end)
                                               ,split=""),forceToLower = FALSE)),as.string = TRUE)[[1]] } 
    else {
      start <- variableRegions$start[i]
      end <- variableRegions$end[i]
      loci <-substring(tprDlocusSeq,start,end) }
    variableRegions$sequence[i] <- loci }
  return(variableRegions)
}

variableRegions <- getVariableSequences(variableRegions)

## Get a count of number of times a putative donor site sequence is represented in the variable region dataset
for (i in c(1:length(variableRegions$start))){
  variableRegions$sequence_count[i] <- sum(str_count(variableRegions[variableRegions$region==variableRegions$region[i],]$sequence,
                                                     variableRegions$sequence[i])) }

## Keep just the variable regions that are unique
variableRegions <- variableRegions[variableRegions$sequence_count==1,]
variableRegions$sequence_count <- NULL

### Expand donor sites that are overlapping to get full-length of potential donor site
### And count how many samples they are in
variableReg <- variableRegions
variableReg <- variableReg[order(-variableReg$length),]
for (i in c(1:(length(variableReg$feature)-1))){
  for (j in c(i:length(variableReg$feature))){
    if (variableReg$region[i] != variableReg$region[j]) next()
    if (i==j) next()
    if ((inrange(variableReg$start[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE)) | 
        (inrange(variableReg$start[i], variableReg$end[j], variableReg$start[j], incbounds=TRUE))){
      variableReg$start[i] <- variableReg$start[j]
      variableReg$count[i] <- variableReg$count[j] + variableReg$count[i]
      variableReg$morethanone[i] <- 2
      variableReg$count[j] <- 0
      variableReg$start[j] <- 0
      variableReg$end[j] <- 0
    } else if ((inrange(variableReg$end[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE)) | 
               (inrange(variableReg$end[i], variableReg$end[j], variableReg$start[j], incbounds=TRUE))){
      variableReg$end[i] <- variableReg$end[j]
      variableReg$count[i] <- variableReg$count[j] + variableReg$count[i]
      variableReg$morethanone[i] <- 2
      variableReg$start[j] <- 0
      variableReg$end[j] <- 0
      variableReg$count[j] <- 0
    }
  }
}

### Get rid of donor site sequences that are now redundant or within another site
variableRegions <- variableReg[!variableReg$start==0,]

## For specificity sake, we are only calling variable regions if they are present in more than one sample.
## Cull all variable regions not present in more than one sample
variableRegions <- variableRegions[variableRegions$morethanone==2,]
variableRegions <- dplyr::mutate(variableRegions,length=abs(start-end)+1)


## Get a count of number of times a sequence is represented in the variable region dataset
## This in case two variable regions became one but had not yet been culled.
variableRegions <- getVariableSequences(variableRegions)
for (i in c(1:length(variableRegions$start))){
  variableRegions$sequence_count[i] <- sum(str_count(variableRegions[variableRegions$region==variableRegions$region[i],]$sequence,
                                                     variableRegions$sequence[i])) }

variableRegions <- variableRegions[variableRegions$sequence_count==1,]

ds_names <- unique.data.frame(df3 %>% select(bitscore,interval,donorsite,color))
ds_seqs <- variableRegions %>% select(name,sequence,start,end)
colnames(ds_seqs) <- c("donorsite","full_ds_seq","start","end")
df <- culture_hits %>% mutate(ds_seq = substr(sequence,qstart,qend)) %>% filter(region=="V6") %>% filter(week==1 | week==10)

# grab donor site names
df1 <- merge(ds_names, df,by=c("bitscore","interval"))
df1 <- merge(df1,ds_seqs,by=c("donorsite"))
write.csv(df1,"donorsites.csv",quote=FALSE,row.names=FALSE) # this file can be used for plotting in plot_donor_sites.py.


######
## Back to internal repeats

write.table(hits2,"hits2.csv",sep=",",row.names=FALSE,quote=FALSE)
write.table(allvariable_edited,"allvariable_edited.csv",sep=",",row.names=FALSE,quote=FALSE)

allvariable_edited <- allvariable
colnames(allvariable) <- c("qseqid","qseqseq")

hits3 <- merge(hits2, allvariable, by=c("qseqid")) 
hits3$count <- as.numeric(str_split_fixed(hits3$qseqid,"_",4)[,3])
hits3 <- hits3 %>% filter(count>=10)
write.table(hits3,"hits3.csv",sep=",",row.names=FALSE,quote=FALSE)

### Run python script find_ds_combinations.py here
### Create figures based on output of above python script.

big_ds <- read.table("acl_donorsite_combinations.csv",sep=",",header=TRUE, fill=TRUE) 
big_ds <- big_ds[-c(6:31)]

# Get number of repeats used
num_repeats <- rowSums(big_ds_edited == "repeat")
big_ds_tophits <- big_ds %>% group_by(qseqid) %>% top_n(n=1, wt=percent_covered) %>% top_n(n=-1, wt=num_segments)
big_ds_tophits <- unique(big_ds_tophits)

# Get count
big_ds_tophits$count <- as.numeric(str_split_fixed(big_ds_tophits$qseqid,"_",4)[,3])

# Get some statistics
all_covered <- (big_ds_tophits[big_ds_tophits$percent_covered<0.95,]) %>% group_by(region) %>% summarise(sum(count))
nrow(all_covered)
sum(all_covered$count)
nrow(big_ds_tophits)
sum(big_ds_tophits$count)

# Grab only "full" coverage
only100 <- big_ds_edited %>% filter(percent_covered >=0.95)
only100$num_segments[only100$percent_covered!=1] <- only100$num_segments[only100$percent_covered!=1] + 1

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

## Look at %GC content, not significant
only100_immuno <- only100 %>% filter(!grepl("W",qseqid)) %>% separate(qseqid, c("SampleName","Other"),sep="_") %>% separate(SampleName, c("Sample","Week"),sep="-")
only100_immuno$Week[only100_immuno$Week=="inoculum"] <- 0
only100_culture <- only100 %>% filter(grepl("W",qseqid)) %>% mutate(Week = substr(qseqid, 2,3))
only100_immuno$qseqlength <- nchar(only100_immuno$qseqseq)
only100_immuno$qseq_g <- str_count(only100_immuno$qseqseq, "G")
only100_immuno$qseq_c <- str_count(only100_immuno$qseqseq, "C")
only100_immuno <- only100_immuno %>% mutate(qseq_gc = qseq_g + qseq_c, qseq_gc_content = qseq_gc / qseqlength * 100)
only100_culture$qseqlength <- nchar(only100_culture$qseqseq)
only100_culture$qseq_g <- str_count(only100_culture$qseqseq, "G")
only100_culture$qseq_c <- str_count(only100_culture$qseqseq, "C")
only100_culture <- only100_culture %>% mutate(qseq_gc = qseq_g + qseq_c, qseq_gc_content = qseq_gc / qseqlength * 100)
only100_culture_V6 <- only100_culture %>% filter(region == "V6")
only100_culture_V4 <- only100_culture %>% filter(region=="V1")
t.test(only100_culture_V6$qseq_gc_content,only100_culture_V4$qseq_gc_content, var.equal=TRUE,
       alternative="greater")
allvar_immuno <- allvar %>% filter(!grepl("W",seq_name)) %>% separate(seq_name, c("SampleName","Other"),sep="_") %>% separate(SampleName, c("Sample","Week"),sep="-")
allvar_immuno$Week[allvar_immuno$Week=="inoculum"] <- 0
allvar_culture <- allvar %>% filter(grepl("W",seq_name)) %>% mutate(Week = substr(seq_name, 2,3))
allvar_immuno$qseqlength <- nchar(allvar_immuno$sequence)
allvar_immuno$qseq_g <- str_count(allvar_immuno$sequence, "G")
allvar_immuno$qseq_c <- str_count(allvar_immuno$sequence, "C")
allvar_immuno <- allvar_immuno %>% mutate(qseq_gc = qseq_g + qseq_c, qseq_gc_content = qseq_gc / qseqlength * 100)
allvar_culture$qseqlength <- nchar(allvar_culture$sequence)
allvar_culture$qseq_g <- str_count(allvar_culture$sequence, "G")
allvar_culture$qseq_c <- str_count(allvar_culture$sequence, "C")
allvar_culture <- allvar_culture %>% mutate(qseq_gc = qseq_g + qseq_c, qseq_gc_content = qseq_gc / qseqlength * 100)

ggplot(allvar_immuno,aes(x=Week, y=qseq_gc_content)) + geom_boxplot()
ggplot(allvar_culture, aes(x = Week, y = qseq_gc_content)) + 
  stat_summary(aes(y=qseq_gc_content), fun.y=mean, geom = "line")

# Sup Fig 6: Plot out # segments used per variable region
ggplot(only100, aes(x=region, y=num_segments)) + geom_boxplot() + 
  theme_clean() + 
  scale_y_continuous(limits=c(0,12), breaks=seq(0,12,3)) + labs(y="Number of Donor Site Segments", x = "Variable Region")
ggsave(filename="donorsite_usage.pdf",height = 3, width=5)
