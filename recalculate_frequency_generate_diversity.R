# Recalculates relative frequency given a pre-filtered file (i.e. allreads_filtered.csv).
# Also makes diversity csvs. 

require("optparse")

filename <- "allreads_filtered.csv"
allreads_filtered <- read.table(filename, sep=',', header=TRUE)

# Grabs the actual number of samples.
numsamples <- (length(colnames(allreads_filtered)) - 2) / 2
metadata <- read.table(opt$metadata, sep=',', header=TRUE)
sample_names <- c(as.character(metadata$SampleName))

allreads_filtered1 <- allreads_filtered

# Loops through samples and recalculates frequency for each sample (column) for each region.
for (sample in c(1:(numsamples))){
  freqcol <- (sample * 2) + 1
  countcol <- (sample * 2) + 2
  #for (region in unique(allreads_filtered1$Region)){
  #  allreads_filtered1[which(allreads_filtered1[,1]==region),freqcol] <- 
  #    allreads_filtered1[which(allreads_filtered1[,1]==region),countcol] / sum(allreads_filtered1[which(allreads_filtered1[,1]==region),countcol],na.rm=TRUE) * 100
  #}
  
  # Change to recalculate based on rf instead of count
  for (region in unique(allreads_filtered1$Region)){
    allreads_filtered1[which(allreads_filtered1[,1]==region),freqcol] <- 
      allreads_filtered1[which(allreads_filtered1[,1]==region),freqcol] / sum(allreads_filtered1[which(allreads_filtered1[,1]==region),freqcol],na.rm=TRUE) * 100
  }
  allreads_filtered1[,freqcol] <- trimws(format(allreads_filtered1[,freqcol], digits = 4, nsmall = 4))
}


# Rewrites the file with the newly calculated relative frequencies.
write.csv(allreads_filtered1,file=filename,row.names=FALSE,quote=FALSE)

# Also makes diversity csvs
if (grepl("allreads_filtered.csv",filename)) {
  amin_csv <- t(allreads_filtered1)
  variable_region_lists <- split(allreads_filtered1, f = allreads_filtered1$Region)
  for (i in 1:length(variable_region_lists)) {
    amin_csv <- t(variable_region_lists[[i]])
    colnames(amin_csv) <- NULL
    to_remove <- "Region"
    for (row in 1:nrow(amin_csv)) {
      row_name <- rownames(amin_csv)[row]
      if (grepl("_RelativeFreq", row_name)) {
        to_remove <- c(to_remove, row_name)
      }
    }
    removedamin_csv = amin_csv[ !(row.names(amin_csv) %in% to_remove), ]
    path <- gsub("[\r\n]", "", path)
    dir.create(file.path(path, "Diversity"), showWarnings=FALSE)
    write.table(removedamin_csv,file=paste(path,"/Diversity/V",i,".csv",sep=""), sep=",",quote=FALSE, col.names=FALSE)
  }
}