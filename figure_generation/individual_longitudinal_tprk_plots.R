# Generates dot-line plots for comparing the variable regions between longitudinal samples. (Sup Figs 3-4)

list.of.packages <- c("ggplot2", "grid", "dplyr", "scales", "gridExtra", "RColorBrewer", "optparse","randomcoloR", "cowplot",
                      "tidyr", "tibble", "foreach","iterators","doParallel","ggtext")
lapply(list.of.packages,library,character.only=T)

setwd("/Users/uwvirongs/Documents/Michelle/longitudinal_tprk/redo_longitudinal_tprk_with_doubleprep/validating_dna_reads")
path<-getwd()

allreads <- paste(path,"/allreads.csv", sep="")
allreads_filtered <- paste(path,"/allreads_filtered.csv", sep="")

alldata <- read.csv(allreads_filtered,header=TRUE,sep=",",stringsAsFactors = FALSE)

myColors <- distinctColorPalette(length(alldata$Read))
names(myColors) <- levels(alldata$Read)
colScale <- scale_colour_manual(name = NULL, guide = FALSE, values = myColors)

alldata_old <- alldata
colors <- data.frame(color=myColors)
alldata <- cbind(alldata,colors)

sample_names <- c("W01glystock","W02","W03a","W04a","W05a","W06a","W07a","W08a","W09a","W10a")
sample_names <- c("W01glystock","W02","W03b","W04b","W05b","W06b","W07b","W08b","W09b","W10b")
sample_names <- c("W01glystock","W05site3","W06site2","W07site4","W08site5","W08site7","W10site2","W10site4","W10site5")

sample_names <- c("7386.inoculum","7397.1","7397.2","7397.3")#,"7397.4")#,"7397.7")
sample_names <- c("7386.inoculum","7400.1","7400.2","7400.3","7400.4","7400.5")#,"7400.6","7400.7")
sample_names <- c("7386.inoculum","7401.1","7401.2","7401.3","7401.4","7401.5")#,"7401.6","7401.7")
sample_names <- c("7386.inoculum","7403.1","7403.2","7403.3","7403.4","7403.5")#,"7403.6","7403.7")
sample_names <- c("7386.inoculum","7406.1","7406.2","7406.3","7406.4","7406.5")
sample_names <- c("7386.inoculum","7391.1","7391.2","7391.3")
sample_names <- c("7386.inoculum","7396.1","7396.2","7396.3","7396.4")#,"7396.7")
sample_names <- c("7386.inoculum","7399.1","7399.2","7399.3","7399.4")
sample_names <- c("7386.inoculum","7411.1","7411.2","7411.3","7411.5")#,"7411.6","7411.7")
sample_names <- c("7386.inoculum","7413.1","7413.2","7413.3","7413.5")#,"7413.6","7413.7")

samp_name = "7391"

rfcol <- paste("Ill_",sample_names[1],"_RelativeFreq",sep = "")
rfcol2 <- paste("Ill_",sample_names[2],"_RelativeFreq",sep = "")
rfcol3 <- paste("Ill_",sample_names[3],"_RelativeFreq",sep = "")
rfcol4 <- paste("Ill_",sample_names[4],"_RelativeFreq",sep = "")
rfcol5 <- paste("Ill_",sample_names[5],"_RelativeFreq",sep = "")
rfcol6 <- paste("Ill_",sample_names[6],"_RelativeFreq",sep = "")
rfcol7 <- paste("Ill_",sample_names[7],"_RelativeFreq",sep = "")
rfcol8 <- paste("Ill_",sample_names[8],"_RelativeFreq",sep = "")
rfcol9 <- paste("Ill_",sample_names[9],"_RelativeFreq",sep = "")
rfcol10 <- paste("Ill_",sample_names[10],"_RelativeFreq",sep = "")


commondfIllumina <- dplyr::select(alldata,Region,Read,color,rfcol,rfcol2,rfcol3,rfcol4)#,rfcol5)#,rfcol6)#,rfcol7,rfcol8,rfcol9)#,rfcol10)
commondfIllumina <- filter(commondfIllumina,!((commondfIllumina[[rfcol]] == 0) & (commondfIllumina[[rfcol2]] == 0) & 
                                                (commondfIllumina[[rfcol3]] == 0) & (commondfIllumina[[rfcol4]] == 0)))# & 
                                               # (commondfIllumina[[rfcol5]] == 0)))# & 
                                                #(commondfIllumina[[rfcol6]] == 0)))# & 
                                                #(commondfIllumina[[rfcol7]] == 0) & 
                                                #(commondfIllumina[[rfcol8]] == 0) &
                                                #(commondfIllumina[[rfcol9]] == 0)))# & 
                                                #(commondfIllumina[[rfcol10]] == 0)))
                                              
sortedIllumina <- commondfIllumina[order(commondfIllumina$Region,-commondfIllumina[[rfcol]],-commondfIllumina[[rfcol2]]),]
sortedIllumina <- gather(sortedIllumina,rfcol,rfcol2,rfcol3,rfcol4,key="Sample",value="Frequency")

sortedIlluminaold <- sortedIllumina

# Only keep reads that have multiple samples
sortedIllumina <- sortedIllumina[complete.cases(sortedIllumina), ]
sortedIllumina <- sortedIllumina %>% group_by(Region, Read) %>% dplyr::filter(n()>1) %>% ungroup()

sortedIllumina$Sample[sortedIllumina$Sample == rfcol] <- sample_names[1]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol2] <- sample_names[2]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol3] <- sample_names[3]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol4] <- sample_names[4]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol5] <- sample_names[5]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol6] <- sample_names[6]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol7] <- sample_names[7]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol8] <- sample_names[8]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol9] <- sample_names[9]
sortedIllumina$Sample[sortedIllumina$Sample == rfcol10] <- sample_names[10]

#ct_values <- read_csv("copynumbers_culture.csv") %>% filter(SampleName %in% c(sample_names[1],sample_names[2],sample_names[3],sample_names[4],sample_names[5],
#                                                                              sample_names[6],sample_names[7],sample_names[8],sample_names[9],sample_names[10]))

ct_values <- read_tsv("copynumbers.tsv") %>% filter(SampleName %in% c(sample_names[1],sample_names[2],sample_names[3],sample_names[4],sample_names[5],
                                                                             sample_names[6]))#,sample_names[7]))#,sample_names[8]))#,sample_names[9]))#,sample_names[10]))


#sortedIllumina$Sample <- factor(sortedIllumina$Sample, levels = unique(sortedIllumina$Sample), ordered = TRUE)
colScale <- scale_colour_manual(name = NULL, guide = FALSE, values = myColors)

# Time to select which reads have legends. Select any reads that have >=10% at any time point.
legend_reads <- sortedIllumina
legend_reads <- legend_reads %>% group_by(Read) %>% filter(any(Frequency>=10))

sortedIllumina$Sample <- gsub("\\.",'-', sortedIllumina$Sample)
#sortedIllumina_df <- extractDaysFromImmuno(sortedIllumina)

colnames(ct_values) <- c("Sample","TP47","Sequenced","CFTR","Copy_Number","logTP47","logCFTR","logcopy")
ct_values$Sample <- gsub("\\.",'-', ct_values$Sample)
ct_values <- extractDaysFromImmuno(ct_values)

extractDaysFromImmuno <- function(df) {
  new_df <- data.frame()
  for (samplename in unique(df$Sample)) {
    subset <- df %>% filter(Sample==samplename)
    if(grepl("-inoculum",samplename)) {
      subset$day = 0
    }
    if(grepl("-1",samplename)) {
      subset$day = 7
    }
    if(grepl("-2",samplename)) {
      subset$day = 14
    }
    if(grepl("-3",samplename)) {
      subset$day = 21
    }
    if(grepl("-4",samplename)) {
      subset$day = 28
    }
    if(grepl("-5",samplename)) {
      subset$day = 35
    }
    if(grepl("-6",samplename)) {
      subset$day = 42
    }
    if(grepl("-7",samplename)) {
      subset$day = 49
    }
    new_df <- rbind(new_df,subset)
  }
  return(new_df)
}
# Extract days

#sortedIllumina <- sortedIllumina_df
#sortedIllumina$Sample <- sortedIllumina$day

#sortedIllumina <- sortedIllumina%>% mutate(day = substr(Sample, 2,3))
#ct_values <- ct_values %>% mutate(day = substr(Sample,2,3))
#sortedIllumina$Week[sortedIllumina$Week=="inoculum"] <- 0


#sortedIllumina$day <- as.numeric(as.character(sortedIllumina$day))
#ct_values$day <- as.numeric(as.character(ct_values$day))

#sortedIllumina$Sample <- sortedIllumina$day

out <- by(data=sortedIllumina, INDICES = sortedIllumina$Region, FUN = function(m) {
  m <- droplevels(m)
  region_num <- m$Region[1]
  # V1 keeps y-axis scale and title, everything else does not have one
  if (region_num == "V1") {
    m <- ggplot(m) + geom_point(aes(y = Frequency, x = Sample, color=color)) + geom_line(aes(y = Frequency, x = Sample, group=Read, color=color)) +  
      theme_clean() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),plot.background = element_blank()) +
      #scale_x_discrete(name = "") + 
      #scale_x_continuous(name = "",breaks=c(0,7,14,21,28,35),limits = c(0,35)) +
      theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) + labs(x="") + 
      scale_color_identity("",guide = "legend", labels=legend_reads$Read, breaks=legend_reads$color) +
      theme(legend.position = c(0.37, 0.95), legend.background = element_blank(), legend.text = element_text(size = 5), legend.key.size = unit(0.5, 'lines')) +  
      facet_grid(~Region) + theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm")) + ylim(0,119) + 
      theme(axis.text=element_text(size=8)) #+ theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
  } else {
    m <- ggplot(m) + geom_point(aes(y = Frequency, x = Sample, color=color)) + geom_line(aes(y = Frequency, x = Sample, group=Read, color=color)) +  
      theme_clean() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),plot.background = element_blank()) +
      #scale_x_continuous(name = "",breaks=c(0,7,14,21,28,35),limits = c(0,35)) +
      theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +  labs(x="") + 
      theme(axis.title.y = element_blank()) + theme(axis.text.y = element_blank()) + theme(axis.ticks.y = element_blank()) + 
      scale_color_identity("",guide = "legend", labels=legend_reads$Read, breaks=legend_reads$color) +
      facet_grid(~Region) + theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm")) + ylim(0,119) + 
      theme(axis.text=element_text(size=8))# + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) + 
    
    # Now we adjust legend positions manually because widths are variably adjusted later to account for extra y-axis in V1
    if (region_num == "V2" | region_num == "V3" | region_num == "V4") {
      m <- m + theme(legend.position = c(0.4, 0.95), legend.background = element_blank(), 
                     legend.text = element_text(size = 5), legend.key.size = unit(0.5, 'lines'))
    } else if (region_num == "V5") {
      m <- m + theme(legend.position = c(0.45, 0.95), legend.background = element_blank(), 
                     legend.text = element_text(size = 5), legend.key.size = unit(0.5, 'lines'))
    } else if (region_num == "V6" | region_num == "V7") {
      m <- m + theme(legend.position = c(0.51, 0.95), legend.background = element_blank(), 
                     legend.text = element_text(size = 5), legend.key.size = unit(0.5, 'lines'))       
    }
    # 
    # if (region_num == "V4") {
    #   m <- m + scale_x_discrete(name = "")
    # } else {
    #   m <- m + scale_x_discrete(name = "")
    # }
  }
})

a <- grid.arrange(grobs = out, nrow=1, top = textGrob(paste(samp_name, " - Immunocompetent",sep=""), gp=gpar(cex = 1.15)), 
                  widths = c(1.32,1,1,1,1,1,1))

ct_values$logcopy <- log10(ct_values$Copy_Number)
ct_values$group <- "group"

b <- ggplot(ct_values, aes(x=Sample, y=logcopy, group=group)) + geom_line() + geom_point() + theme_clean() + 
  scale_y_continuous(limits = c(-4,4)) + labs(y=(log[10]~frac(italic(tp47)~copies,CFTR~copies))) + theme(plot.background=element_blank())
  #scale_x_continuous(name = "Day", limits=c(0,35),breaks=c(0,7,14,21,28,35))
  #scale_x_continuous(name = "",breaks=c(0,2,4,6,8,10),limits = c(0,10))

plot_grid(a,NULL,b, ncol=1, rel_heights = c(1,0.02,0.3))

suppressMessages(ggsave(paste(samp_name,"_Longitudinal_DotLine_Filtered.pdf",sep=""),width=11,height=7,units="in"))