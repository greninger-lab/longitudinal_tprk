# Calculates diversity statistics from csvs in diversity folder (generated from recalculate_frequency_amin.R.)
# Also generates figures 2-4.

library(vegan)
library(dplyr)

#V1
V1=read.table("Diversity/V1.csv", header=TRUE, sep=",") #%>% filter(Read != "Ill_7413.5_Count")
i <- V1[1]
isolates <- data.frame(i)
V1[1] <- NULL
V1[is.na(as.matrix(V1))] = 0
V1[is.nan(as.matrix(V1))] = 0
V1_variants <- specnumber(V1)
V1_diversity <- diversity(V1)
V1_evenness <- V1_diversity/log(specnumber(V1))
V1_vregion <- V1_diversity
V1_vregion[1:length(V1_vregion)] <- "V1"

#V2
V2 =  read.table("Diversity/V2.csv", header=TRUE, sep=",") #%>% filter(Read != "Ill_7413.5_Count")
V2[1] <- NULL
V2[is.na(as.matrix(V2))] = 0
V2[is.nan(as.matrix(V2))] = 0
V2_variants <- specnumber(V2)
V2_diversity <- diversity(V2)
V2_evenness <- V2_diversity/log(specnumber(V2))
V2_vregion <- V2_diversity
V2_vregion[1:length(V2_vregion)] <- "V2"


#V3
V3 =  read.table("Diversity/V3.csv", header=TRUE, sep=",") #%>% filter(Read != "Ill_7413.5_Count")
V3[1] <- NULL
V3[is.na(as.matrix(V3))] = 0
V3[is.nan(as.matrix(V3))] = 0
V3_variants <- specnumber(V3)
V3_diversity <- diversity(V3)
V3_evenness <- V3_diversity/log(specnumber(V3))
V3_vregion <- V3_diversity
V3_vregion[1:length(V3_vregion)] <- "V3"

#V4
V4 =  read.table("Diversity/V4.csv", header=TRUE, sep=",") #%>% filter(Read != "Ill_7413.5_Count")
V4[1] <- NULL
V4[is.na(as.matrix(V4))] = 0
V4[is.nan(as.matrix(V4))] = 0
V4_variants <- specnumber(V4)
V4_diversity <- diversity(V4)
V4_evenness <- V4_diversity/log(specnumber(V4))
V4_vregion <- V4_diversity
V4_vregion[1:length(V4_vregion)] <- "V4"

#V5
V5 = read.table("Diversity/V5.csv", header=TRUE, sep=",")# %>% filter(Read != "Ill_7413.5_Count")
V5[1] <- NULL
V5[is.na(as.matrix(V5))] = 0
V5[is.nan(as.matrix(V5))] = 0
V5_variants <- specnumber(V5)
V5_diversity <- diversity(V5)
V5_evenness <- V5_diversity/log(specnumber(V5))
V5_vregion <- V5_diversity
V5_vregion[1:length(V5_vregion)] <- "V5"

#V6
V6 =  read.table("Diversity/V6.csv", header=TRUE, sep=",") # %>% filter(Read != "Ill_7413.5_Count")
V6[1] <- NULL
V6[is.na(as.matrix(V6))] = 0
V6[is.nan(as.matrix(V6))] = 0
V6_variants <- specnumber(V6)
V6_diversity <- diversity(V6)
V6_evenness <- V6_diversity/log(specnumber(V6))
V6_vregion <- V6_diversity
V6_vregion[1:length(V6_vregion)] <- "V6"

#V7
V7 =  read.table("Diversity/V7.csv", header=TRUE, sep=",") #%>% filter(Read != "Ill_7413.5_Count")
V7[1] <- NULL
V7[is.na(as.matrix(V7))] = 0
V7[is.nan(as.matrix(V7))] = 0
V7_variants <- specnumber(V7)
V7_diversity <- diversity(V7)
V7_evenness <- V7_diversity/log(specnumber(V7))
V7_vregion <- V7_diversity
V7_vregion[1:length(V7_vregion)] <- "V7"

# Combine all variable regions
all_diversity <- c(V1_diversity, V2_diversity, V3_diversity, V4_diversity, V5_diversity, V6_diversity,V7_diversity)
all_vregion <- c(V1_vregion, V2_vregion, V3_vregion, V4_vregion, V5_vregion, V6_vregion, V7_vregion)
all_variants <- c(V1_variants, V2_variants, V3_variants, V4_variants, V5_variants,V6_variants,V7_variants)
all_evenness <- c(V1_evenness, V2_evenness, V3_evenness, V4_evenness, V5_evenness, V6_evenness, V7_evenness)

df <-cbind(isolates,all_vregion,all_variants,all_diversity,all_evenness)

# join diversity metrics together + remove NaN if no. of variant = 1
#df <- cbind(isolates, V1_variants, V1_evenness, V1_diversity, V2_variants, V2_evenness, V2_diversity, V3_variants, V3_evenness, V3_diversity, V4_variants, V4_evenness, V4_diversity, V5_variants, V5_evenness, V5_diversity, V6_variants, V6_evenness, V6_diversity, V7_variants, V7_evenness, V7_diversity) #Full_Length_variants, #Full_Length_evenness, Full_Length_diversity)
df <- df %>% filter(all_variants > 0)
df[is.na(as.matrix(df))] <- 0
df[is.nan(as.matrix(df))] <- 0
df <- df %>% dplyr::rowwise() %>% dplyr::mutate(Read = strsplit(Read,"_")[[1]][2])
write.table(df,"Table_S1.csv",sep=",",index_col=FALSE)


## Divide into different groups
immuno_set <- df %>% filter(!grepl("W",Read)) %>% separate(Read, c("Sample","Week"),sep="\\.")
culture_set <- df %>% filter(grepl("W",Read)) %>% mutate(Week = substr(Read, 2,3))

immuno_set$Week[immuno_set$Week=="inoculum"] <- 0
immuno_set$Week <- as.numeric(as.character(immuno_set$Week))
culture_set$Week <- as.numeric(as.character(culture_set$Week))

immunosuppressed <- immuno_set %>% filter(grepl('7397',Sample) | grepl('7400',Sample) | grepl('7401',Sample) |
                                            grepl('7403',Sample) | grepl('7406',Sample ))
immunocompetent <- immuno_set %>% filter(grepl('7391',Sample) | grepl('7396',Sample) | grepl('7399',Sample) |
                                           grepl('7411',Sample) | grepl('7413',Sample ))
inoculum <- immuno_set %>% filter(grepl('7386',Sample))

cultures <- culture_set %>% filter(Week != 1 & Week!=2) %>% filter(!grepl("site",Read))
rabbit <- culture_set %>% filter(Week != 1 & Week!=2) %>% filter(grepl("site",Read)) %>% filter(!grepl("W04",Read))
inoculum2 <- culture_set %>% filter(Week==1 | Week==2)

immunosuppressed <- immunosuppressed %>% group_by(Week,all_vregion) %>% summarise_at(vars(all_variants,all_diversity,all_evenness), funs(mean, sd))
immunocompetent <- immunocompetent %>% group_by(Week,all_vregion) %>% summarise_at(vars(all_variants,all_diversity,all_evenness), funs(mean, sd))
inoculum <- inoculum %>% group_by(Week,all_vregion)%>% summarise_at(vars(all_variants,all_diversity,all_evenness), funs(mean, sd))

cultures <- cultures %>% ungroup() %>% group_by(Week,all_vregion)%>% summarise_at(vars(all_variants,all_diversity,all_evenness), funs(mean, sd))
rabbit <- rabbit %>% ungroup() %>% group_by(Week,all_vregion)%>% summarise_at(vars(all_variants,all_diversity,all_evenness), funs(mean, sd))
inoculum2 <- inoculum2 %>% ungroup() %>% group_by(Week,all_vregion)%>% summarise_at(vars(all_variants,all_diversity,all_evenness), funs(mean, sd))

immunosuppressed$Status <- "Immunosuppressed"
immunocompetent$Status <- "Immunocompetent"
inoculum$Status <- "Inoculum"

cultures$Status <- "Cultured"
rabbit$Status <- "Rabbit"
inoculum2$Status <- "Inoculum"

# Grab copy numbers
immuno_ct_values <- read_tsv("copynumbers.tsv") %>% separate(SampleName, c("Sample","Week"),sep="\\.")
immuno_ct_values$Week[immuno_ct_values$Week=="inoculum"] <- 0
immunosuppressed_ct <- immuno_ct_values %>% filter(grepl('7397',Sample) | grepl('7400',Sample) | grepl('7401',Sample) |
                                                     grepl('7403',Sample) | grepl('7406',Sample ))
immunocompetent_ct <- immuno_ct_values %>% filter(grepl('7391',Sample) | grepl('7396',Sample) | grepl('7399',Sample) |
                                                    grepl('7411',Sample) | grepl('7413',Sample ))
inoculum_ct <- immuno_ct_values %>% filter(grepl('7386',Sample))
immunosuppressed_ct$Status <- "Immunosuppressed"
immunocompetent_ct$Status <- "Immunocompetent"
inoculum_ct$Status <- "Inoculum"
immuno_ct_values <- rbind(immunosuppressed_ct,immunocompetent_ct,inoculum_ct)
immuno_ct_sum <- immuno_ct_values %>%group_by(Status,Week) %>% summarise(mean(Copy_Number),sd(Copy_Number))
culture_ct_values <- read_csv("copynumbers_culture.csv")%>% mutate(Week = substr(SampleName, 2,3))
culture_ct_values$Week <- as.numeric(as.character(culture_ct_values$Week))
cultures_ct <- culture_ct_values %>% filter(Week != 1 & Week!=2) %>% filter(!grepl("site",SampleName))
rabbit_ct <- culture_ct_values %>% filter(Week != 1 & Week!=2) %>% filter(grepl("site",SampleName)) #%>% filter(!grepl("W04",SampleName))
inoculum2_ct <- culture_ct_values %>% filter(Week==1 | Week==2)
cultures_ct$Status <- "Cultured"
rabbit_ct$Status <- "Rabbit"
inoculum2_ct$Status <- "Inoculum"
culture_ct_values <- rbind(cultures_ct, rabbit_ct, inoculum2_ct)



df1 <- rbind(inoculum, immunosuppressed, immunocompetent)
df2 <- rbind(inoculum2, cultures, rabbit)
df2$Status <- factor(df2$Status, levels = c("Rabbit","Cultured","Inoculum"))
df3 <- rbind(df1, df2)

## Fig 3: Pielou's evenness plots 
immuno_plots <- ggplot(df1, aes(x = Week, y = all_evenness_mean, group=Status,color=Status)) + 
  geom_pointrange(aes(ymin=all_evenness_mean-all_evenness_sd, ymax=all_evenness_mean+all_evenness_sd)) + geom_line() +
  #scale_shape_manual(values=c(1,16)) + 
  facet_wrap(~all_vregion, nrow=1) + theme_bw() + scale_color_brewer(palette="Dark2") + 
  scale_x_continuous(breaks=seq(1,10, by=1),limits=c(1,5)) + theme_clean() + theme(plot.background = element_blank()) + 
  scale_y_continuous(limits=c(0,1)) + 
  theme(legend.position = "right", legend.background = element_blank(), legend.text = element_text(size=8)) + 
  labs(y ="Mean Pielou's Evenness Score", color = "", x = "Week Post-Infection")
culture_plots <- ggplot(df2, aes(x = Week - 3, y = all_evenness_mean, group=Status,color=Status)) + 
  geom_pointrange(aes(ymin=all_evenness_mean-all_evenness_sd, ymax=all_evenness_mean+all_evenness_sd)) + geom_line() +
  #scale_shape_manual(values=c(1,16)) + 
  facet_wrap(~all_vregion, nrow=1) + theme_bw() + scale_color_brewer(palette="Dark2") + 
  scale_x_continuous(breaks=seq(1,7, by=1),limits=c(1,7)) + theme_clean() + theme(plot.background = element_blank()) + 
  scale_y_continuous(limits=c(0,1)) + 
  theme(legend.position = "right", legend.background = element_blank(), legend.text = element_text(size=8)) + 
  labs(y ="Mean Pielou's Evenness Score", color = "", x = "Week Post-Infection")

plot_grid(immuno_plots,culture_plots,labels=c("A","B"), nrow=2, align = "v")
ggsave("Fig3_edit.pdf",width = 8.5, height = 6, units = "in")


library(lme4)

# Look at some stats (csvs combined into Table 1)
fit <- (lmer(all_variants_mean~Week*all_vregion + (1+Week|all_vregion),data=df2_culture))
coeffs <- summary(fit, ddf = "Kenward-Roger")
write.table(coeffs, "cultured_stats.csv", sep = ",")
df1_immunosuppressed <- df1 %>% filter(Status == "Immunosuppressed") %>% filter(Week >0)
fit <- (lmer(all_variants_mean~Week*all_vregion + (1+Week|all_vregion),data=df1_immunosuppressed))
coeffs <- summary(fit, ddf = "Kenward-Roger")$coefficients
write.table(coeffs, "immunosuppressed_stats.csv", sep = ",")


# Fig 4: # unique variants in immmunosuppressed and cultured
df4 <- df3 %>% filter(Status == "Immunosuppressed" | Status == "Cultured")
df4$Week[df4$Status == "Cultured"] <- df4$Week[df4$Status == "Cultured"] - 3
variants_plots <- ggplot(df4, aes(x = Week, y = all_variants_mean, group=Status,color=Status)) + 
  geom_pointrange(aes(ymin=all_variants_mean-all_variants_sd, ymax=all_variants_mean+all_variants_sd)) + geom_line() +
  #scale_shape_manual(values=c(1,16)) + 
  facet_wrap(~all_vregion, nrow=1) + theme_bw() + scale_color_brewer(palette="Accent") + 
  scale_x_continuous(breaks=seq(1,7, by=1),limits=c(1,7)) + 
  theme_clean() + theme(plot.background = element_blank()) + 
  #scale_y_continuous(limits=c(0,1)) + 
  theme(legend.position = "right", legend.background = element_blank(), legend.text = element_text(size=8)) + 
  labs(y ="Number of Unique Variants", color = "", x = "Week Post-Infection")
ggsave("Fig4.pdf",width = 7, height = 4, units = "in")

# Fig 1: copy number values
options(scipen = 999)
immuno_ct_values$Week <- as.numeric(as.character(immuno_ct_values$Week))
immuno_ct_values$logcopy <- log10(immuno_ct_values$Copy_Number)
culture_ct_values$logcopy <- log10(culture_ct_values$Copy_Number)

immuno_ct_plot <- ggplot(immuno_ct_values, aes(x=Week, y = logcopy, group=Sample,color=Status, shape=Sequenced)) + 
  geom_point(alpha=0.7) + scale_y_continuous(limits=c(-4,4),breaks=(seq(-4,4,by=1))) + 
  #geom_point(alpha=0.7) + scale_y_continuous(limits=c(0,1.25),breaks=(seq(0,1.25,by=0.25))) + 
  scale_shape_manual(values=c(4,16),guide=FALSE) + 
  theme_clean() + 
  stat_summary(aes(y=logcopy,group=Status), fun.y=mean, geom = "line") +
  scale_color_brewer(palette="Dark2") + labs(y=expression(log[10]~(italic(tp47)~copies/CFTR~copies)),color="") + 
  scale_x_continuous(limits=c(0,5)) + 
  theme(legend.background = element_blank(), plot.background=element_blank(),
                        legend.text = element_text(size = 8), 
                      legend.direction="vertical", legend.margin = margin(0,0,0,0, "cm"),
                      legend.spacing.y = unit(-0.7,"cm"), legend.key.height = unit(0.1,"cm"),
                        legend.position = "bottom") 

culture_ct_values$Status <- factor(culture_ct_values$Status, levels = c("Rabbit","Cultured","Inoculum"))

culture_ct_values$Week <- culture_ct_values$Week - 3
culture_ct_values <- culture_ct_values %>% filter(Status != "Rabbit")
culture_ct_plot <- ggplot(culture_ct_values, aes(x=Week, y = logcopy, color=Status, shape=Sequenced)) + 
  geom_point(alpha=0.7) + scale_y_continuous(limits=c(-4,4),breaks=(seq(-4,4,by=1))) + 
  scale_shape_manual(values=c(4,16), guide=FALSE) + 
  stat_summary(aes(y=logcopy,group=Status), fun.y=mean, geom = "line") + 
  scale_color_brewer(palette="Dark2") + labs(y=expression(log[10]~(italic(tp47)~copies/CFTR~copies))) + 
  scale_x_continuous(limits=c(-2,7),breaks=c(-2,-1,0,1,2,3,4,5,6,7)) + 
  theme_clean() +theme(legend.background = element_blank(), plot.background=element_blank(),
                       legend.text = element_text(size = 8), 
                       legend.direction="vertical", legend.margin = margin(0,0,0,0, "cm"),
                       legend.spacing.y = unit(-0.7,"cm"), legend.key.height = unit(0.1,"cm"),
                       legend.position = "bottom")+ labs(color="")

plot_grid(immuno_ct_plot,culture_ct_plot, labels=c("A","B"),align = "h", rel_widths = c(1,1))
ggsave("Copynumber_plots_denominator_adjusted.pdf", width=7,height=4,units="in")

## Sup Fig 1: non denominator-adjusted copy number plots
immuno_ct_plot <- ggplot(immuno_ct_values, aes(x=Week, y =logTP47, group=Sample,color=Status, shape=Sequenced)) + 
  geom_point(alpha=0.7) + #scale_y_continuous(limits=c(-4,4),breaks=(seq(-4,4,by=1))) + 
  #geom_point(alpha=0.7) + scale_y_continuous(limits=c(0,1.25),breaks=(seq(0,1.25,by=0.25))) + 
  scale_shape_manual(values=c(4,16),guide=FALSE) + 
  theme_clean() + 
  stat_summary(aes(y=logTP47,group=Status), fun.y=mean, geom = "line") +
  scale_color_brewer(palette="Dark2") + labs(y=expression(log[10]~(italic(tp47)~"copies/µl"))) + 
  scale_x_continuous(limits=c(0,6)) + 
  theme(legend.background = element_blank(), plot.background=element_blank(),
        legend.text = element_text(size = 8), 
        legend.direction="vertical", legend.margin = margin(0,0,0,0, "cm"),
        legend.spacing.y = unit(-0.7,"cm"), legend.key.height = unit(0.1,"cm"),
        legend.position = "bottom") 

culture_ct_plot <- ggplot(culture_ct_values, aes(x=Week, y =logTP47, color=Status, shape=Sequenced)) + 
  geom_point(alpha=0.7) + 
  scale_y_continuous(limits=c(0,6),breaks=(seq(0,6,by=1))) + 
  scale_shape_manual(values=c(4,16), guide=FALSE) + 
  stat_summary(aes(y=logTP47,group=Status), fun.y=mean, geom = "line") + 
  scale_color_brewer(palette="Dark2") + labs(y=expression(log[10]~(italic(tp47)~"copies/µl"))) + 
  scale_x_continuous(limits=c(-2,7),breaks=seq(-2,7,by=1)) + 
  theme_clean() +theme(legend.background = element_blank(), plot.background=element_blank(),
                       legend.text = element_text(size = 8), 
                       legend.direction="vertical", legend.margin = margin(0,0,0,0, "cm"),
                       legend.spacing.y = unit(-0.7,"cm"), legend.key.height = unit(0.1,"cm"),
                       legend.position = "bottom")+ labs(color="")

plot_grid(immuno_ct_plot,culture_ct_plot, labels=c("A","B"),align = "h", rel_widths = c(1,1))
ggsave("sup_copynumber_plots_denominator_adjusted.pdf", width=7,height=4,units="in")
