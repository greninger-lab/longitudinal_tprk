# Verifies that reads are contained in both Illumina preps at a mean rf of >0.1%.
# Also makes Illumina_squared files for each individual sample.

list.of.packages <- c("dplyr", "tidyr", "tibble", "foreach","iterators","doParallel")
lapply(list.of.packages,library,character.only=T)

setwd("/Users/uwvirongs/Documents/Michelle/longitudinal_tprk/redo_longitudinal_tprk_with_doubleprep/validating_dna_reads")

metadata <- read.table("metadata.csv", sep=',', header=TRUE)
sample_names = as.character(metadata$SampleName)
illumina1_table = as.character(metadata$Illumina1)
illumina2_table = as.character(metadata$Illumina2)

# Relative frequency threshold
rf_cutoff = 0.1

## Comparing Illuminas ##
Illumina_freq_files = list.files('./',pattern="*_final_dna_data.csv")
Illumina_freq_files_fullpath <- Illumina_freq_files
compare_Illumina_df <- data.frame(Region=character(),Read=character())

sample_names <- gsub("-","\\.",sample_names)

i=1

for (i in 1:length(sample_names)) {
  IlluminaFreqtitle <- paste("Ill_",sample_names[i],"_RelativeFreq",sep='')
  IlluminaCounttitle <- paste("Ill_",sample_names[i],"_Count",sep='')
  IlluminaFreq2title <- paste("Ill2_",sample_names[i],"_RelativeFreq",sep='')
  IlluminaCount2title <- paste("Ill2_",sample_names[i],"_Count",sep='')
  print(IlluminaCounttitle)
  Illuminadf <- read.csv(illumina1_table[i],col.names = c("Region","Read",IlluminaFreqtitle,IlluminaCounttitle),check.names = FALSE)
  Illuminadf <- Illuminadf[order(Illuminadf$Region,-Illuminadf[[IlluminaFreqtitle]]),][]
  
  Illumina2df <- read.csv(illumina2_table[i],col.names = c("Region","Read",IlluminaFreq2title,IlluminaCount2title),check.names = FALSE)
  Illumina2df <- Illumina2df[order(Illumina2df$Region,-Illumina2df[[IlluminaFreq2title]]),][]
  Illumina_vs_Illumina_df <- merge(Illuminadf,Illumina2df,all=TRUE)
  
  write.csv(Illumina_vs_Illumina_df,file=paste(sample_names[i],"_prefiltered.csv",sep=""),row.names=FALSE,quote=FALSE)
  

  # Only take reads that are in both Illumina dataframes
  commondfIllumina <- filter(Illumina_vs_Illumina_df,!(is.na(Illumina_vs_Illumina_df[[IlluminaFreq2title]]) | is.na((Illumina_vs_Illumina_df[[IlluminaFreqtitle]]))))
  Illumina_vs_Illumina_name <- paste(sample_names[i],"_Illumina_squared.csv",sep="")
  write.csv(commondfIllumina, file=Illumina_vs_Illumina_name,row.names=FALSE,quote=FALSE)
  
  system(paste("Rscript /Users/uwvirongs/Documents/Michelle/longitudinal_tprk/redo_longitudinal_tprk_with_doubleprep/longitudinal_tprk_code/recalculate_frequency.R -m metadata.csv -f ",
               Illumina_vs_Illumina_name))
  
  commondfIllumina <- read.csv(Illumina_vs_Illumina_name,sep=",", check.names=FALSE)

  #common_mean_Illumina <- commondfIllumina %>% transmute(Read,Region,
  #                                                       rf = rowMeans(select(.,IlluminaFreqtitle,IlluminaFreq2title)),
  #                                                       count = rowMeans(select(.,IlluminaCounttitle,IlluminaCount2title)))
  
  # Use our second Illumina set of sample for counts and mean the rfs
  common_mean_Illumina <- commondfIllumina %>% transmute(Read,Region,
                                                         rf = rowMeans(select(.,IlluminaFreqtitle)),
                                                         count = rowMeans(select(.,IlluminaCounttitle)))
                                                         #rf = rowMeans(select(.,IlluminaFreqtitle,IlluminaFreq2title)),
                                                         #count = rowMeans(select(.,IlluminaCounttitle,IlluminaCount2title)))

  common_mean_Illumina <- common_mean_Illumina %>% filter(rf >= rf_cutoff)
  names(common_mean_Illumina)[names(common_mean_Illumina) == "rf"] <- IlluminaFreqtitle
  names(common_mean_Illumina)[names(common_mean_Illumina) == "count"] <- IlluminaCounttitle
  common_mean_Illumina <- common_mean_Illumina[c(2,1,3,4)]
  common_mean_Illumina_name <- paste(sample_names[i],"_validated_over5_count.csv",sep="")
  write.csv(common_mean_Illumina,file=common_mean_Illumina_name,row.names=FALSE,quote=FALSE)
  system(paste("Rscript /Users/uwvirongs/Documents/Michelle/longitudinal_tprk/redo_longitudinal_tprk_with_doubleprep/longitudinal_tprk_code/recalculate_frequency.R -m metadata.csv -f ",
               common_mean_Illumina_name))
  
  compare_Illumina_df <- merge(compare_Illumina_df,common_mean_Illumina,all=TRUE)
}

write.csv(compare_Illumina_df,file="Illumina_validated_Illumina_allreads.csv",row.names=FALSE,quote=FALSE)
system(paste("Rscript /Users/uwvirongs/Documents/Michelle/longitudinal_tprk/redo_longitudinal_tprk_with_doubleprep/longitudinal_tprk_code/recalculate_frequency.R -m metadata.csv -f Illumina_validated_Illumina_allreads.csv"))

# Also now convert Illumina_validated_Illumina_allreads.csv to regular allreads.csv with illumina_to_allreads.py. 