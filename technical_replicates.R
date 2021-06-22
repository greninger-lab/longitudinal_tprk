# Compares cultured technical replicates.
# Generates S2 Fig.  

setwd("/Users/uwvirongs/Documents/Michelle/longitudinal_tprk/redo_longitudinal_tprk_with_doubleprep/validating_dna_reads")
metadata <- read.table("metadata.csv", sep=',', header=TRUE)
sample_names = as.character(metadata$SampleName)
Illumina_freq_files = list.files('./',pattern="*_prefiltered.csv")

immuno_tech_comparisons <- data.frame()
culture_tech_comparisons <- data.frame()

immunosup_snames <- c("7397","7400","7401","7403","7406")
immunocomp_snames <- c("7391","7396","7399","7411","7413")
culture_snames <- cultures_ct$SampleName
culture_rab_snames <- rabbit_ct$SampleName

for (i in 1:length(Illumina_freq_files)) {
  print(Illumina_freq_files[i])
  individual_comp <- read.csv(Illumina_freq_files[i],check.names = FALSE)[-c(4,6)]
  sampname = str_split_fixed(colnames(individual_comp[3]),"_",3)[,2]
  colnames(individual_comp) <- c("Region","Read","Prep1_rf","Prep2_rf")
  individual_comp$SampName <- sampname
  if(grepl("W",sampname)) {
    individual_comp <- individual_comp %>% dplyr::mutate(Week = substr(SampName, 2,3))
    if(sampname %in% culture_rab_snames) {
      individual_comp$Status <- "Rabbit"
    } else {
      individual_comp$Status <- "Culture"
    }
    culture_tech_comparisons <- rbind(culture_tech_comparisons, individual_comp)
  }
  else {
    individual_comp <- individual_comp %>% separate(SampName,c("Sample","Week"),remove=FALSE,sep="\\.")
    individual_comp$Week[individual_comp$Week=="inoculum"] <- 0
    
    if(unique(individual_comp$Sample) %in% immunosup_snames) {
      individual_comp$Status <- "Immunosuppressed"
    } else {
      individual_comp$Status <- "Immunocompetent"
    }
    individual_comp$Status[individual_comp$Week==0] <- "Inoculum"
    immuno_tech_comparisons <- rbind(immuno_tech_comparisons, individual_comp)
  }
}
immuno_tech_comparisons <- immuno_tech_comparisons %>% filter(SampName != "7413.5") %>% filter(!grepl(".6",SampName)) %>% filter(!grepl(".7",SampName))

culture_tech_comparisons$Week <- as.character(as.numeric(culture_tech_comparisons$Week) - 3)

immuno_tech_comparisons[is.na(immuno_tech_comparisons)] <- 0
culture_tech_comparisons[is.na(culture_tech_comparisons)] <- 0

tech_plot_immunocomp <- ggplot(subset(immuno_tech_comparisons, Status %in% "Immunocompetent"), aes(x=Prep1_rf, y=Prep2_rf, color=Week)) + geom_point() +
  theme_clean() + 
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  #geom_smooth(method="lm",se=FALSE) +
  #facet_wrap(~ Region, nrow=2) + 
  scale_color_brewer(palette="Spectral") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Library Prep 1 RF (%)", y = "Library Prep 2 RF (%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank())

tech_plot_immunosup <- ggplot(subset(immuno_tech_comparisons, Status %in% "Immunosuppressed"), aes(x=Prep1_rf, y=Prep2_rf, color=Week)) + geom_point() +
  theme_clean() + 
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  #geom_smooth(method="lm",se=FALSE) +
  #facet_wrap(~ Region, nrow=2) + 
  scale_color_brewer(palette="Spectral") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Library Prep 1 RF (%)", y = "Library Prep 2 RF (%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank())

tech_plot_cultrabbit <- ggplot(subset(culture_tech_comparisons, Status %in% "Rabbit"), aes(x=Prep1_rf, y=Prep2_rf, color=Week)) + geom_point() +
  theme_clean() + 
  scale_color_brewer(palette="Spectral") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Library Prep 1 RF (%)", y = "Library Prep 2 RF (%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank())

culture_tech_comparisons <- culture_tech_comparisons %>% filter(Week >=0)
tech_plot_cult <- ggplot(subset(culture_tech_comparisons, Status %in% "Culture"), aes(x=Prep1_rf, y=Prep2_rf, color=Week)) + geom_point() +
  theme_clean() + 
  scale_color_brewer(palette="Spectral") + 
  geom_abline(slope = 1, linetype="dashed", color="gray") + 
  labs(x = "Library Prep 1 RF (%)", y = "Library Prep 2 RF (%)") + 
  theme(legend.background = element_blank(), plot.background = element_blank())

lm_eqn(lm(Prep2_rf ~ Prep1_rf, all_tech_comparisons))

plot_grid(tech_plot_immunocomp, tech_plot_immunosup, tech_plot_cultrabbit, tech_plot_cult, labels=c("A","B","C","D"))
ggsave(filename="SupFig_technical_replicates.pdf",height = 5, width=7, units="in")
  