# Compares two biological cultural replicates and generates figures 5-6. 
Illumina_freq_files = list.files('./',pattern="*_validated_over5_count.csv")
Illumina_freq_files_fullpath <- Illumina_freq_files
compare_Illumina_df <- data.frame(Region=character(),Read=character())

sample_names <- c("W01glystock","W02","W03","W04","W05","W06","W07","W08","W09","W10")

allreads <- read.table("allreads_filtered.csv",sep=",", header=TRUE)
culture_reads <- allreads[-c(3:88)]
culture_reorganized <- data.frame()

i=1
for (i in 1:length(sample_names)) {
  week = substr(sample_names[i],2,3)
  IlluminaFreqtitle <- paste("Ill_",sample_names[i],"a_RelativeFreq",sep='')
  IlluminaCounttitle <- paste("Ill_",sample_names[i],"a_Count",sep='')
  IlluminaFreq2title <- paste("Ill_",sample_names[i],"b_RelativeFreq",sep='')
  IlluminaCount2title <- paste("Ill_",sample_names[i],"b_Count",sep='')
  print(IlluminaFreqtitle)
  
  a <- cbind(culture_reads['Region'],culture_reads['Read'],culture_reads[IlluminaFreqtitle],culture_reads[IlluminaFreq2title])
  a <- a[rowSums(is.na(a[3:4])) < 2L,]
  a[is.na(a)] <- 0
  
  a$Week <- week
  colnames(a) <- c("Region","Read","Culture_A","Culture_B","Week")
  culture_reorganized <- rbind(culture_reorganized,a)
}

# Find r2 for culture replicate relative frequency reproducibility.
lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 10));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}
lm_eqn(lm(Culture_B ~ Culture_A, culture_reorganized))

# Reorganize weeks 
culture_reorganized$Week <- as.numeric(culture_reorganized$Week) - 3
culture_reorganized$Week <- as.character(culture_reorganized$Week)

# Reformat for plots
culture_comparisons <- culture_reorganized %>% filter(!(Read %in% inoculum_reads))
culture_comparisons <- culture_comparisons %>% unite("title", Region, Read, sep=": ", remove = FALSE)

### Figure generation
# Fig 5A
full_plot <- ggplot(culture_reorganized, aes(x=Culture_A, y=Culture_B, color=Region)) + geom_point() +
  theme_clean() + 
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  #geom_smooth(method="lm",se=FALSE) +
  #facet_wrap(~ Region, nrow=2) + 
  scale_color_brewer(palette="Dark2") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Culture A RF (%)", y = "Culture B RF (%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank())

# Fig 5B - novel sequences
full_plot_new <- ggplot(culture_comparisons, aes(x=Culture_A, y=Culture_B, color=Region)) + geom_point() +
  theme_clean() + 
  scale_x_continuous(limits=c(0,0.8)) + scale_y_continuous(limits=c(0,0.8)) + 
  scale_color_brewer(palette="Dark2") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Culture A RF (%)", y = "Culture B RF (%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank())

# Looking at sequences that are only in one culture replicate
onlyone_cult <- culture_comparisons %>% filter(Culture_A == 0 | Culture_B == 0)

# Looking at low frequency reads
culture_reorganized_low <- culture_reorganized_filtered %>% group_by(Read) %>% 
        filter(all(Culture_A <= 4)) %>% filter(all(Culture_A >=0.5))
culture_reorganized_low <- culture_reorganized_filtered %>% group_by(Read) %>% 
  filter(any(Culture_B >=0.5)) #%>% filter(all(Culture_B <= 1))
culture_v1 <- culture_reorganized_filtered %>% filter(Culture_A <= 4) %>% filter(Region == "V1") %>% filter(Read == "GIASEKNGGAQPLKH")
culture_reorganized_low <- rbind(culture_reorganized_low, culture_v1)
culture_reorganized_low <- culture_reorganized_low %>% unite("title", Region, Read, sep=": ", remove = FALSE)

culture_comparisons2 <- culture_comparisons %>% group_by(Read) %>% filter(any(Culture_B >= 0.4))

# Fig 6: new alleles above 0.4 rf
lowplot <- ggplot(culture_comparisons2[order(culture_comparisons2$Week),], 
       aes(x=Culture_A, y=Culture_B, color=Week, group = Read)) + geom_point() + geom_path() + 
  theme_clean() + 
  facet_wrap(~ title, nrow = 5) + 
  scale_x_continuous(limits=c(0,1.5)) + scale_y_continuous(limits=c(0,1.5)) + 
  scale_color_brewer(palette="Spectral") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Culture A Relative Frequency (%)", y = "Culture B Relative Frequency (%)", color="") + 
  theme(legend.background = element_blank(), plot.background = element_blank(), panel.spacing=unit(1,"lines"),
        strip.text = element_text(size = 6), legend.position = "None")
ggsave("Culture_replicates_by_read.pdf", width=6, height = 7, units="in")

# Fig 5C: comparing culture replicates at high rfs
culture_highs <- culture_reorganized %>% filter(Culture_A >= 60)
highplot <- ggplot(culture_reorganized[order(culture_reorganized$Week),], 
       aes(x=Culture_A, y=Culture_B, color=Week, group=Read)) + geom_point() + geom_path() + 
  theme_clean() + 
  facet_wrap(~ Region, nrow=1) + 
  scale_x_continuous(limits=c(70,100), breaks=seq(70,100,10)) + scale_y_continuous(limits=c(70,100), breaks=seq(70,100,10)) + 
  scale_color_brewer(palette="Spectral") + 
   geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Culture A RF (%)", y = "Culture B RF(%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank(), panel.spacing=unit(1,"lines"))

# Fig 5D: regression lines by week
weekplot <- ggplot(culture_reorganized[order(culture_reorganized$Week),], 
                   aes(x=Culture_A, y=Culture_B, color=Week)) + geom_point() + 
  geom_smooth(method = "lm", fill = "NA") + 
  theme_clean() + 
  #facet_wrap(~ Region, nrow=2) + 
  scale_x_continuous(limits=c(0.5,4)) + scale_y_continuous(limits=c(0.5,4)) + 
  scale_color_brewer(palette="Spectral") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Culture A RF (%)", y = "Culture B RF (%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank()) 

# Combine sub-panels in Fig 5
top_row <- plot_grid(full_plot,full_plot_new,nrow=1, labels=c("A","B"))
middle_row <- plot_grid(highplot, weekplot, nrow=2, labels=c("C","D"), rel_heights=c(1,2.5))
plot_grid(top_row, middle_row, labels=c("",""), rel_heights=c(1,2),nrow=2)
ggsave("Fig5_Culture_Replicates_Rf_Plot.pdf", width=8.5, height = 9, units="in")



# Sup Table 2: looking at weekly rate of decrease of clonal sequence
rate <- culture_highs %>% group_by(Region) %>% filter(Week == 0 | Week == 7)
rate <- rate %>% mutate(mean_rf = (Culture_A + Culture_B)/2)

rate_df <- data.frame(matrix(ncol=2,nrow=0))
colnames(rate_df) <- c("Vregion","Rate")
for(variable_region in unique(rate$Region)) {
  rate_region <- rate %>% filter(Region == variable_region)
  diff <- (rate_region$mean_rf[rate_region$Week==7] - rate_region$mean_rf[rate_region$Week==0])/7
  to_add <- c(variable_region, diff)
  rate_df <- rbind(rate_df, to_add)
}

colnames(rate_df) <- c("Vregion","Rate")
write_csv(rate_df, file="culture_rate.csv")

clonal_seqs <- allreads %>% filter(Ill_7386.inoculum_RelativeFreq >= 70)
clonal_seqs <- clonal_seqs %>% select("Region","Ill_7386.inoculum_RelativeFreq",
                                      "Ill_7400.1_RelativeFreq","Ill_7400.2_RelativeFreq","Ill_7400.3_RelativeFreq","Ill_7400.4_RelativeFreq","Ill_7400.5_RelativeFreq",
                                      "Ill_7401.1_RelativeFreq","Ill_7401.2_RelativeFreq","Ill_7401.3_RelativeFreq","Ill_7401.4_RelativeFreq","Ill_7401.5_RelativeFreq",
                                      "Ill_7403.1_RelativeFreq","Ill_7403.2_RelativeFreq","Ill_7403.3_RelativeFreq","Ill_7403.4_RelativeFreq","Ill_7403.5_RelativeFreq",
                                      "Ill_7406.1_RelativeFreq","Ill_7406.2_RelativeFreq","Ill_7406.3_RelativeFreq","Ill_7406.4_RelativeFreq","Ill_7406.5_RelativeFreq")
clonal_seqs <- melt(clonal_seqs)
clonal_seqs$Week <- str_split_fixed(str_split_fixed(clonal_seqs$variable,"\\.",2)[,2],"_",2)[,1]

clonal_seqs_mean <- clonal_seqs %>% filter(Week != "inoculum") %>% group_by(Region,Week) %>% summarise(value = mean(value))
clonal_seqs_inoculum <- clonal_seqs %>% filter(Week == "inoculum")
clonal_seqs_mean <- rbind(clonal_seqs_mean,clonal_seqs_inoculum)

rate_df <- data.frame(matrix(ncol=2,nrow=0))
colnames(rate_df) <- c("Vregion","Rate")
for(variable_region in unique(clonal_seqs_mean$Region)) {
  rate_region <- clonal_seqs_mean %>% filter(Region == variable_region)
  diff <- (rate_region$value[rate_region$Week==5] - rate_region$value[rate_region$Week==1])/5
  to_add <- c(variable_region, diff)
  rate_df <- rbind(rate_df, to_add)
}

colnames(rate_df) <- c("Vregion","Rate")
write_csv(rate_df, file="immunosuppressed_rate.csv")
