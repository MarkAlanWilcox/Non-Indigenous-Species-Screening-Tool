## NISST ##

#Load packages
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(ggpubr)

####Generating scores and error using assessor generated distributions-----
#Read in master data.This will be the csv file used during data entry and provided on the GitHub page

data <- read.csv("Master_data.csv")

#Parse into separate files for impacts and invasion. 
#This is necessary due to the different scoring rubrics (3 and 5 point scale) and thus 
# must be separated to facilitate application of Monte Carlo simulations.

data_inv <- data[(data$Category == "Invasion"),]
data_imp <- data[(data$Category != "Invasion"),]


#Apply Monte Carlo techniques

# Invasion module

#Parse data to separate Modifier questions which have NA for distribution values and thus can not
#be incorporated in the for loop

INV_1 <-  data_inv[is.na(data_inv$O1),]
INV_2 <-  data_inv %>% drop_na(O1)

#Apply for loop to generate scores and error 
###NOTE that the total for the generated distribution across all outcomes MUST equal 100

for (i in 1:nrow(INV_2)){
  v <- sample(x=c(1,2,3,4,5), 
              size = 1000, 
              replace = TRUE, 
              prob = c(INV_2$O1[i], 
                       INV_2$O2[i], 
                       INV_2$O3[i], 
                       INV_2$O4[i], 
                       INV_2$O5[i]))  
  
  INV_2$Mean[i] <- mean(v)  
  INV_2$Sd[i]  <- sd(v)
  
}

#Recombine the Modifier questions and new scores where the climate modifier score has been 
#applied to the mean column

INV_1$Mean <- INV_1$Modifier.Score
INV_1$Sd <- NA
Inv_mod <- rbind(INV_2,INV_1)

#Square standard deviations for error propagation

Inv_mod$Sd <- (Inv_mod$Sd)^2

## Impact modules

#Parse data to separate Modifier questions which have NA for distribution values and thus can not
#be incorporated in the for loop

IMP_1 <-  data_imp[is.na(data_imp$O1),]
IMP_2 <-  data_imp %>% drop_na(O1)

#Apply for loop to generate scores and error
###NOTE that the total for the generated distribution across all outcomes MUST equal 100

for (i in 1:nrow(IMP_2)){
  v <- sample(x=c(1,2,3), 
              size = 1000, 
              replace = TRUE, 
              prob = c(IMP_2$O1[i], 
                       IMP_2$O2[i], 
                       IMP_2$O3[i]))  
  
  IMP_2$Mean[i] <- mean(v)  
  IMP_2$Sd[i]  <- sd(v)
  
}

#Recombine the Modifier questions and new scores where the climate modifier score has been 
# applied to the mean column

IMP_1$Mean <- IMP_1$Modifier.Score
IMP_1$Sd <- NA
IMP_mod <- rbind(IMP_2,IMP_1)

#Square standard deviations for error propagation

IMP_mod$Sd <- (IMP_mod$Sd)^2

#Separate impact modules

Soc_imp_mod <- IMP_mod[IMP_mod$Category == "Socio-Econ",]
Eco_imp_mod <- IMP_mod[IMP_mod$Category == "Ecological",]


#####For each species, for each assessor, for each ecoregion, calculate the total scores for each 
#module and error

#Invasion module

## Reframe the mean values and error such that each question occupies a column to 
#facilitate calculation of total scores. Taxa and P/A are removed and will be add back later
Inv_Output <- Inv_mod %>% group_by(Species, 
                                   Assessor, 
                                   Ecoregion.Assessed, 
                                   Category)%>% reframe("S1.1" = Mean[Question.No==1.1], 
                                                        "S1.2" = Mean[Question.No==1.2],
                                                        "S1.3" = Mean[Question.No==1.3],
                                                        "S2.1" = Mean[Question.No==2.1],
                                                        "S2.2" = Mean[Question.No==2.2],
                                                        "S3.1" = Mean[Question.No==3.1],
                                                        "S3.2" = Mean[Question.No==3.2],
                                                        "S3.3" = Mean[Question.No==3.3],
                                                        "S4.1" = Mean[Question.No==4.1],
                                                        "S4.2" = Mean[Question.No==4.2],
                                                        "S4.3" = Mean[Question.No==4.3],
                                                        Module_error = sum(Sd,na.rm = TRUE))

# Create new columns to calculate module modified scores
Inv_Output$S1 <- NA
Inv_Output$S1C <- NA
Inv_Output$S2 <- NA
Inv_Output$S2C <- NA
Inv_Output$S3 <- NA
Inv_Output$S3C <- NA
Inv_Output$S4 <- NA
Inv_Output$S4C <- NA
Inv_Output$Module_score <- NA
Inv_Output$Module_climate_score <- NA

#Apply modifiers

for (i in 1:nrow(Inv_Output)){
  
  Inv_Output$S1[i] <- if(sum(Inv_Output$S1.1[i],Inv_Output$S1.2[i])>5){5}else{sum(Inv_Output$S1.1[i],Inv_Output$S1.2[i])}
  Inv_Output$S1C[i] <- if(sum(Inv_Output$S1[i],Inv_Output$S1.3[i])>5){5}else{sum(Inv_Output$S1[i],Inv_Output$S1.3[i])}
  Inv_Output$S2[i] <- Inv_Output$S2.1[i]
  Inv_Output$S2C[i] <- if(sum(Inv_Output$S2[i],Inv_Output$S2.2[i])>5){5}else{sum(Inv_Output$S2[i],Inv_Output$S2.2[i])}
  Inv_Output$S3[i] <- min(Inv_Output$S3.1[i],Inv_Output$S3.2[i])
  Inv_Output$S3C[i] <- if(sum(Inv_Output$S3[i],Inv_Output$S3.3[i])>5){5}else{sum(Inv_Output$S3[i],Inv_Output$S3.3[i])}
  Inv_Output$S4[i] <- if(sum(Inv_Output$S4.1[i],Inv_Output$S4.2[i])>5){5}else{sum(Inv_Output$S4.1[i],Inv_Output$S4.2[i])}
  Inv_Output$S4C[i] <- if(sum(Inv_Output$S4[i],Inv_Output$S4.3[i])>5){5}else{sum(Inv_Output$S4[i],Inv_Output$S4.3[i])}
  
}

#Calculate module scores

for (i in 1:nrow(Inv_Output)){  
  Inv_Output$Module_score[i] <- sum(Inv_Output$S1[i],Inv_Output$S2[i],Inv_Output$S4[i],Inv_Output$S3[i])
  Inv_Output$Module_climate_score[i] <- sum(Inv_Output$S1C[i],Inv_Output$S2C[i],Inv_Output$S4C[i],Inv_Output$S3C[i])
}

#Scale (average) scores. If modification to the total number of questions, 
#this value of 4 must be reflective of the no of questions

Inv_Output$Module_score <- Inv_Output$Module_score/4
Inv_Output$Module_climate_score <- Inv_Output$Module_climate_score/4

#Calculate error by taking the sqrt and dividing by the number of questions (4)
Inv_Output$Module_error <- sqrt(Inv_Output$Module_error)/4

#Ecological module

## Reframe the mean values and error such that each question occupies a column to 
#facilitate calculation of total scores. Taxa and P/A are removed and will be add back later
Ecol_Output <- Eco_imp_mod %>% group_by(Species, 
                                        Assessor, 
                                        Ecoregion.Assessed,
                                        Category) %>% reframe("S1.1" = Mean[Question.No==1.1],
                                                              "S1.2" = Mean[Question.No==1.2],
                                                              "S1.3" = Mean[Question.No==1.3],
                                                              "S1.4" = Mean[Question.No==1.4],
                                                              "S2.1" = Mean[Question.No==2.1],
                                                              "S2.2" = Mean[Question.No==2.2],
                                                              "S2.3" = Mean[Question.No==2.3],
                                                              "S2.4" = Mean[Question.No==2.4],
                                                              "S3.1" = Mean[Question.No==3.1],
                                                              "S3.2" = Mean[Question.No==3.2],
                                                              "S3.3" = Mean[Question.No==3.3],
                                                              "S3.4" = Mean[Question.No==3.4],
                                                              "S4.1" = Mean[Question.No==4.1],
                                                              "S4.2" = Mean[Question.No==4.2],
                                                              "S4.3" = Mean[Question.No==4.3],
                                                              Module_error = sum(Sd,na.rm = TRUE))

# Create new columns to calculate module modified scores
Ecol_Output$S1 <- NA
Ecol_Output$S1C <- NA
Ecol_Output$S2 <- NA
Ecol_Output$S2C <- NA
Ecol_Output$S3 <- NA
Ecol_Output$S3C <- NA
Ecol_Output$S4 <- NA
Ecol_Output$S4C <- NA
Ecol_Output$Module_score <- NA
Ecol_Output$Module_climate_score <- NA

#Apply modifiers

for (i in 1:nrow(Ecol_Output)){
  
  Ecol_Output$S1[i] <- sum(Ecol_Output$S1.1[i],Ecol_Output$S1.2[i],Ecol_Output$S1.3[i])
  Ecol_Output$S1C[i] <- if(sum(Ecol_Output$S1[i],Ecol_Output$S1.4[i])>9){9}else{sum(Ecol_Output$S1[i],Ecol_Output$S1.4[i])}
  Ecol_Output$S2[i] <- sum(Ecol_Output$S2.1[i],Ecol_Output$S2.2[i],Ecol_Output$S2.3[i])
  Ecol_Output$S2C[i] <- if(sum(Ecol_Output$S2[i],Ecol_Output$S2.4[i])>9){9}else{sum(Ecol_Output$S2[i],Ecol_Output$S2.4[i])}
  Ecol_Output$S3[i] <- sum(Ecol_Output$S3.1[i],Ecol_Output$S3.2[i],Ecol_Output$S3.3[i])
  Ecol_Output$S3C[i] <- if(sum(Ecol_Output$S3[i],Ecol_Output$S3.4[i])>9){9}else{sum(Ecol_Output$S3[i],Ecol_Output$S3.4[i])}
  Ecol_Output$S4[i] <- sum(Ecol_Output$S4.1[i],Ecol_Output$S4.2[i])
  Ecol_Output$S4C[i] <- if(sum(Ecol_Output$S4[i],Ecol_Output$S4.3[i])>6){6}else{sum(Ecol_Output$S4[i],Ecol_Output$S4.3[i])}
  
}

#Calculate module scores

for (i in 1:nrow(Ecol_Output)){  
  Ecol_Output$Module_score[i] <- sum(Ecol_Output$S1[i],Ecol_Output$S2[i],Ecol_Output$S4[i],Ecol_Output$S3[i])
  Ecol_Output$Module_climate_score[i] <- sum(Ecol_Output$S1C[i],Ecol_Output$S2C[i],Ecol_Output$S4C[i],Ecol_Output$S3C[i])
}

#Scale (average) scores. If modification to the total number of questions, 
#this value of 11 must be reflective of the no of questions

Ecol_Output$Module_score <- Ecol_Output$Module_score/11
Ecol_Output$Module_climate_score <- Ecol_Output$Module_climate_score/11

#Calculate error by sqrt and dividing by number of questions (11)
Ecol_Output$Module_error <- sqrt(Ecol_Output$Module_error)/11

#Socioeconomic module

## Reframe the mean values and error such that each question occupies a column to 
#facilitate calculation of total scores. Taxa and P/A are removed and will be add back later

Soc_Output <- Soc_imp_mod %>% group_by(Species, 
                                       Assessor, 
                                       Ecoregion.Assessed, 
                                       Category) %>% reframe("S1.1" = Mean[Question.No==1.1],
                                                             "S1.2" = Mean[Question.No==1.2],
                                                             "S1.3" = Mean[Question.No==1.3],
                                                             "S1.4" = Mean[Question.No==1.4],
                                                             "S2.1" = Mean[Question.No==2.1],
                                                             "S2.2" = Mean[Question.No==2.2],
                                                             "S3.1" = Mean[Question.No==3.1],
                                                             "S3.2" = Mean[Question.No==3.2],
                                                             "S3.3" = Mean[Question.No==3.3],
                                                             "S4.1" = Mean[Question.No==4.1],
                                                             "S4.2" = Mean[Question.No==4.2],
                                                             "S4.3" = Mean[Question.No==4.3],
                                                             Module_error = sum(Sd,na.rm = TRUE))
                                                                                                    
# Create new columns to calculate module modified scores
Soc_Output$S1 <- NA
Soc_Output$S1C <- NA
Soc_Output$S2 <- NA
Soc_Output$S2C <- NA
Soc_Output$S3 <- NA
Soc_Output$S3C <- NA
Soc_Output$S4 <- NA
Soc_Output$S4C <- NA
Soc_Output$Module_score <- NA
Soc_Output$Module_climate_score <- NA

#Apply modifiers

for (i in 1:nrow(Soc_Output)){
  
  Soc_Output$S1[i] <- sum(Soc_Output$S1.1[i],Soc_Output$S1.2[i],Soc_Output$S1.3[i])
  Soc_Output$S1C[i] <- if(sum(Soc_Output$S1[i],Soc_Output$S1.4[i])>9){9}else{sum(Soc_Output$S1[i],Soc_Output$S1.4[i])}
  Soc_Output$S2[i] <- Soc_Output$S2.1[i]
  Soc_Output$S2C[i] <- if(sum(Soc_Output$S2[i],Soc_Output$S2.2[i])>3){3}else{sum(Soc_Output$S2[i],Soc_Output$S2.2[i])}
  Soc_Output$S3[i] <- sum(Soc_Output$S3.1[i],Soc_Output$S3.2[i])
  Soc_Output$S3C[i] <- if(sum(Soc_Output$S3[i],Soc_Output$S3.3[i])>6){6}else{sum(Soc_Output$S3[i],Soc_Output$S3.3[i])}
  Soc_Output$S4[i] <- sum(Soc_Output$S4.1[i],Soc_Output$S4.2[i])
  Soc_Output$S4C[i] <- if(sum(Soc_Output$S4[i],Soc_Output$S4.3[i])>6){6}else{sum(Soc_Output$S4[i],Soc_Output$S4.3[i])}
  
}

#Calculate module scores

for (i in 1:nrow(Soc_Output)){  
  Soc_Output$Module_score[i] <- sum(Soc_Output$S1[i],Soc_Output$S2[i],Soc_Output$S4[i],Soc_Output$S3[i])
  Soc_Output$Module_climate_score[i] <- sum(Soc_Output$S1C[i],Soc_Output$S2C[i],Soc_Output$S4C[i],Soc_Output$S3C[i])
}

#Scale (average) scores. If modification to the total number of questions, 
#this value of 8 must be reflective of the no of questions

Soc_Output$Module_score <- Soc_Output$Module_score/8
Soc_Output$Module_climate_score <- Soc_Output$Module_climate_score/8

#Calculate error by sqrt and dividing by number of questions (8)
Soc_Output$Module_error <- sqrt(Soc_Output$Module_error)/8

#Remove unwanted columns and bind all

Inv_Output <- Inv_Output[c(1:4,16,25,26)]
Ecol_Output <- Ecol_Output[c(1:4,20,29,30)]
Soc_Output <- Soc_Output[c(1:4,17,26,27)]
Total_assessment <- rbind(Inv_Output, Soc_Output, Ecol_Output)

#Create new file for the averaged scores across all assessors

Final_df <- Total_assessment

#Square the errors to facilitate the sum of squares

Final_df$Module_error <- (Final_df$Module_error)^2

#Summarise results by averaging across assessors

Final_output <- Final_df %>% group_by(Species,
                                      Ecoregion.Assessed, 
                                      Category)%>% summarise(Mod_error = (sqrt(sum(Module_error))/n()),
                                                             Mod_score = mean(Module_score),
                                                             Mod_score_clim = mean(Module_climate_score))

# Rearrange data frame such that module scores are specific columns organized by species

# Create common column for reintegration

Final_output$id <- paste(Final_output$Species,Final_output$Ecoregion.Assessed)

Final_Ecol <- Final_output[Final_output$Category == "Ecological",]
Final_Soc <- Final_output[Final_output$Category == "Socio-Econ",]
Final_Inv <- Final_output[Final_output$Category == "Invasion",]

#Rename columns
Final_Inv <- Final_Inv %>% rename("Inv_Module_Score" = Mod_score, 
                                  "Inv_Module_Error" = Mod_error,
                                  "Inv_Module_Climate_Score" = Mod_score_clim)
Final_Inv <- Final_Inv[-c(1:3)]

Final_Ecol <- Final_Ecol %>% rename("Ecol_Module_Score" = Mod_score, 
                                    "Ecol_Module_Error" = Mod_error,
                                    "Ecol_Module_Climate_Score" = Mod_score_clim)
Final_Ecol <- Final_Ecol[-c(1:3)]

Final_Soc <- Final_Soc %>% rename("Soc_Module_Score" = Mod_score, 
                                  "Soc_Module_Error" = Mod_error,
                                  "Soc_Module_Climate_Score" = Mod_score_clim)
Final_Soc <- Final_Soc[-c(1:3)]

Final_plotting <- merge(Final_Inv,Final_Ecol,by = c("id"))
Final_plotting <- merge(Final_plotting,Final_Soc,by = c("id"))

# Add presence/absence and taxonomic data back into the plot
# Create an data frame with desired information and same id used for merging "Final_plotting"

test <- data_imp
test$id <- paste(test$Species,test$Ecoregion.Assessed)
test <- test[c(2,3,5,6,16)]
Final_plotting <- merge(Final_plotting, test, by = "id")
Final_plotting <- Final_plotting[!duplicated(Final_plotting$id),]

#Create Total Scores
#Create new dataframe for total scores

Totalscore.df <- Final_plotting

#Create Total Scores by combining modules
Totalscore.df$Inv_Ecol_Error <- sqrt((Totalscore.df$Inv_Module_Error^2)+(Totalscore.df$Ecol_Module_Error^2))
Totalscore.df$Inv_Ecol_Score <- (Totalscore.df$Ecol_Module_Score+Totalscore.df$Inv_Module_Score)
Totalscore.df$Inv_Ecol_Climate_Score <- (Totalscore.df$Ecol_Module_Climate_Score+Totalscore.df$Inv_Module_Climate_Score)

Totalscore.df$Ecol_Soc_Score <- (Totalscore.df$Ecol_Module_Score+Totalscore.df$Soc_Module_Score)
Totalscore.df$Ecol_Soc_Error <- sqrt(((Totalscore.df$Soc_Module_Error)^2)+((Totalscore.df$Ecol_Module_Error)^2))
Totalscore.df$Ecol_Soc_Climate_Score <- (Totalscore.df$Ecol_Module_Climate_Score+Totalscore.df$Soc_Module_Climate_Score)

Totalscore.df$TotalxScore <- (Totalscore.df$Ecol_Soc_Score*Totalscore.df$Inv_Module_Score)
Totalscore.df$TotalxError <- sqrt((Totalscore.df$Ecol_Soc_Error/Totalscore.df$Ecol_Soc_Score)^2+(Totalscore.df$Inv_Module_Error/Totalscore.df$Inv_Module_Score)^2)*Totalscore.df$TotalxScore
Totalscore.df$Total_ClimatexScore <- (Totalscore.df$Ecol_Soc_Climate_Score*Totalscore.df$Inv_Module_Climate_Score)

Totalscore.df$Species <- as.factor(Totalscore.df$Species)
Totalscore.df$Taxa <- as.factor(Totalscore.df$Taxa)
Totalscore.df$Presence.Absense <- as.factor(Totalscore.df$Presence.Absense)

####Ploting-----

library(tidytext) # allows for reordering within facets
library(ggarchery)

#Plot of total scores with Climate modified scores represented using geom_segment. If multiple 
#ecoregions are assessed, use subset(Totalscore.df, Ecoregion.Assessed %in% "Insert Category Here") 
#for data in ggplot

TS_plot <- ggplot(data = Totalscore.df, 
                  aes(y = fct_reorder(Species, TotalxScore),
                      x = TotalxScore,
                      colour = Presence.Absense))+
  geom_errorbar(aes(xmin = TotalxScore-TotalxError, 
                    xmax = TotalxScore+TotalxError), colour = "grey35", alpha = 0.75)+
  geom_point(aes(shape = Taxa), size = 3)+
  geom_segment(aes(y = Species,
                   x = TotalxScore,
                   yend = Species,
                   xend = Total_ClimatexScore),
               arrow = arrow(type = "open", length = unit(0.2, "cm")), 
               colour = "grey35")+
  labs(colour = "Dectected in Area") +
  scale_shape_manual(values= c(15, 16, 17, 18, 4, 5)) +
  scale_colour_manual(labels = c("Not Detected", "Detected"), values = c("#D55E00", "#56B4E9")) +
  xlab("Total Score") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.justification = "left",
        legend.box = "vertical",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8, face = "italic"),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"));TS_plot

ggsave("TotalScores.png", width = 6.5, height = 9)

#Biplots of module scores

## Separate biplots are generated for each impact component plotted against Invasion Risk
#Ecological Impacts vs Invasion Risk (Module scores and separate climate change module scores 
#depicted as a vector plot where the tip of the are represents the climate modified score)
# and Socioeconomic impacts vs Invasion Risk. If multiple ecoregions are assessed, 
# use subset(Totalscore.df, Ecoregion.Assessed %in% "Insert Category Here") 
#for data in ggplot
#To view multiple ecoregions for the same plot, do not use the subset for data but instead use
#facet wrap

ECOvINV_Module <- ggplot(Totalscore.df, 
                         aes(y = Ecol_Module_Score,
                             x = Inv_Module_Score, 
                             colour = Presence.Absense))+
  geom_errorbar(aes(ymin = Ecol_Module_Score-Ecol_Module_Error, 
                    ymax = Ecol_Module_Score+Ecol_Module_Error))+
  geom_errorbarh(aes(xmin = Inv_Module_Score-Inv_Module_Error,
                     xmax = Inv_Module_Score+Inv_Module_Error))+
  geom_point(aes(shape = Taxa), size = 2)+
  labs(colour = "Dectected in Area") +
  scale_colour_manual(labels = c("Not Detected", "Detected"), values = c("#D55E00", "#56B4E9")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
    #    legend.justification = "left",
        legend.box = "vertical",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_x_continuous("Invasion Potential Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3,4,5), 
                     limits = c(0.9,5.1))+
  scale_y_continuous("Ecological Impacts Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1)); ECOvINV_Module

ECOvINV_Climate <- ggplot(Totalscore.df)+
  geom_segment(aes(y = Ecol_Module_Score,
                   x = Inv_Module_Score,
                   yend = Ecol_Module_Climate_Score,
                   xend = Inv_Module_Climate_Score,
                   colour = Presence.Absense),
               arrow = arrow(type = "open", length = unit(0.2, "cm")))+
  labs(colour = "Dectected in Area") +
  scale_colour_manual(labels = c("Not Detected", "Detected"), values = c("#D55E00", "#56B4E9")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
    #    legend.justification = "left",
        legend.box = "vertical",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_x_continuous("Invasion Potential Climate Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3,4,5), 
                     limits = c(0.9,5.1))+
  scale_y_continuous("Ecological Impacts Climate Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1));ECOvINV_Climate

SOCvINV_Module <- ggplot(Totalscore.df, 
                         aes(y = Soc_Module_Score,
                             x = Inv_Module_Score, 
                             colour = Presence.Absense))+
  geom_errorbar(aes(ymin = Soc_Module_Score-Soc_Module_Error, 
                    ymax = Soc_Module_Score+Soc_Module_Error))+
  geom_errorbarh(aes(xmin = Inv_Module_Score-Inv_Module_Error,
                     xmax = Inv_Module_Score+Inv_Module_Error))+
  geom_point(aes(shape = Taxa), size = 2)+
  labs(colour = "Dectected in Area") +
  scale_colour_manual(labels = c("Not Detected", "Detected"), values = c("#D55E00", "#56B4E9")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
 #       legend.justification = "left",
        legend.box = "vertical",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_x_continuous("Invasion Potential Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3,4,5), 
                     limits = c(0.9,5.1))+
  scale_y_continuous("Socio-Economic Impacts Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1)); SOCvINV_Module

SOCvINV_Climate <- ggplot(Totalscore.df)+
  geom_segment(aes(y = Soc_Module_Score,
                   x = Inv_Module_Score,
                   yend = Soc_Module_Climate_Score,
                   xend = Inv_Module_Climate_Score,
                   colour = Presence.Absense),
               arrow = arrow(type = "open", length = unit(0.2, "cm")))+
  labs(colour = "Dectected in Area") +
  scale_colour_manual(labels = c("Not Detected", "Detected"), values = c("#D55E00", "#56B4E9")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
#        legend.justification = "left",
        legend.box = "vertical",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_x_continuous("Invasion Potential Climate Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3,4,5), 
                     limits = c(0.9,5.1))+
  scale_y_continuous("Socio-Economic Impacts Climate Score", 
                     expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1));SOCvINV_Climate

ggarrange(ECOvINV_Module, ECOvINV_Climate, SOCvINV_Module, SOCvINV_Climate, common.legend = TRUE, legend = "bottom")
ggsave("Biplotpanel.png", height = 6.5, width = 6.5, bg="white")

#### Calculate inter assessor scores------

#By Question and Module

##Develop csv files of differences in scores between assessors as a means of QA/QC to be performed
#by an independent reviewer. Highly deviant answers to questions, particularly in the invasion 
#module, likely indicates a difference in available information to assess species, 
#or potentially a data inputting error and should be addressed before finalizing the lists.

##Invasion module Inter-assessor
#Make new working file
IAI.df <- Inv_mod

#Creating unique identifiers to be used to combine data frames
IAI.df$id <- paste(IAI.df$Species,IAI.df$Ecoregion.Assessed,IAI.df$Question.No)

#Separate data frames
IAIAssess.B <- IAI.df[IAI.df$Assessor=="Assess.B",]
IAIAssess.A <- IAI.df[IAI.df$Assessor=="Assess.A",]
IAIAssess.C <- IAI.df[IAI.df$Assessor=="Assess.C",]

#Remove unnecessary columns
IAIAssess.B <- IAIAssess.B[-c(1,10:15,17)]
IAIAssess.A <- IAIAssess.A[-c(1,10:15,17)]
IAIAssess.C <- IAIAssess.C[-c(1,10:15,17)]

#Rename columns
colnames(IAIAssess.B)[9] <- "Assess.B.Mean"
colnames(IAIAssess.A)[9] <- "Assess.A.Mean"
colnames(IAIAssess.C)[9] <- "Assess.C.Mean"

#Remove unnecessary columns
IAIAssess.A <- IAIAssess.A[-c(1:8)]
IAIAssess.C <- IAIAssess.C[-c(1:8)]

#Merge the data files together
IAImaster <- merge(IAIAssess.B, IAIAssess.A, by = "id", all.x=TRUE)
IAImaster1 <- merge(IAImaster, IAIAssess.C, by = "id", all.x=TRUE)

#Remove modifier questions
IAImaster1 <- IAImaster1[!(IAImaster1$Question.No==1.1),]
IAImaster1 <- IAImaster1[!(IAImaster1$Question.No==1.3),]
IAImaster1 <- IAImaster1[!(IAImaster1$Question.No==2.2),]
IAImaster1 <- IAImaster1[!(IAImaster1$Question.No==3.3),]
IAImaster1 <- IAImaster1[!(IAImaster1$Question.No==4.2),]
IAImaster1 <- IAImaster1[!(IAImaster1$Question.No==4.3),]

#Produce columns of differences in question scores within the dataframe
IAImaster1$Assess.B.Assess.C <- IAImaster1$Assess.B.Mean-IAImaster1$Assess.C.Mean
IAImaster1$Assess.B.Assess.A <- IAImaster1$Assess.B.Mean-IAImaster1$Assess.A.Mean
IAImaster1$Assess.C.Assess.A <- IAImaster1$Assess.C.Mean-IAImaster1$Assess.A.Mean

## write a new file for checking in excel
write.csv(IAImaster1, "IAI_check.csv")

##Ecological Impacts module Inter-assessor
#Make new working file
IAEco.df <- Eco_imp_mod

#Creating unique identifiers to be used to combine data frames
IAEco.df$id <- paste(IAEco.df$Species,IAEco.df$Ecoregion.Assessed,IAEco.df$Question.No)

#Separate data frames
IAEcoAssess.B <- IAEco.df[IAEco.df$Assessor=="Assess.B",]
IAEcoAssess.A <- IAEco.df[IAEco.df$Assessor=="Assess.A",]
IAEcoAssess.C <- IAEco.df[IAEco.df$Assessor=="Assess.C",]

#Remove unnecessary columns
IAEcoAssess.B <- IAEcoAssess.B[-c(1,10:15,17)]
IAEcoAssess.A <- IAEcoAssess.A[-c(1,10:15,17)]
IAEcoAssess.C <- IAEcoAssess.C[-c(1,10:15,17)]

#Rename columns
colnames(IAEcoAssess.B)[9] <- "Assess.B.Mean"
colnames(IAEcoAssess.A)[9] <- "Assess.A.Mean"
colnames(IAEcoAssess.C)[9] <- "Assess.C.Mean"

#Remove unnecessary columns
IAEcoAssess.A <- IAEcoAssess.A[-c(1:8)]
IAEcoAssess.C <- IAEcoAssess.C[-c(1:8)]

IAEcomaster <- merge(IAEcoAssess.B, IAEcoAssess.A, by = "id", all.x=TRUE)
IAEcomaster1 <- merge(IAEcomaster, IAEcoAssess.C, by = "id", all.x=TRUE)

#Remove modifier questions
IAEcomaster1 <- IAEcomaster1[!(IAEcomaster1$Question.No==1.4),]
IAEcomaster1 <- IAEcomaster1[!(IAEcomaster1$Question.No==2.4),]
IAEcomaster1 <- IAEcomaster1[!(IAEcomaster1$Question.No==3.4),]
IAEcomaster1 <- IAEcomaster1[!(IAEcomaster1$Question.No==4.3),]

#Produce columns of differences in question scores within the dataframe
IAEcomaster1$Assess.B.Assess.C <- IAEcomaster1$Assess.B.Mean-IAEcomaster1$Assess.C.Mean
IAEcomaster1$Assess.B.Assess.A <- IAEcomaster1$Assess.B.Mean-IAEcomaster1$Assess.A.Mean
IAEcomaster1$Assess.C.Assess.A <- IAEcomaster1$Assess.C.Mean-IAEcomaster1$Assess.A.Mean

## write a new file for checking in excel
write.csv(IAEcomaster1, "IAEco_check.csv")

##SCio-economic Impacts module Inter-assessor
#Make new working file
IASC.df <- Soc_imp_mod

#Creating unique identifiers to be used to combine data frames
IASC.df$id <- paste(IASC.df$Species,IASC.df$Ecoregion.Assessed,IASC.df$Question.No)

#Separate data frames
IASCAssess.B <- IASC.df[IASC.df$Assessor=="Assess.B",]
IASCAssess.A <- IASC.df[IASC.df$Assessor=="Assess.A",]
IASCAssess.C <- IASC.df[IASC.df$Assessor=="Assess.C",]

#Remove unnecessary columns
IASCAssess.B <- IASCAssess.B[-c(1,10:15,17)]
IASCAssess.A <- IASCAssess.A[-c(1,10:15,17)]
IASCAssess.C <- IASCAssess.C[-c(1,10:15,17)]

#Rename columns
colnames(IASCAssess.B)[9] <- "Assess.B.Mean"
colnames(IASCAssess.A)[9] <- "Assess.A.Mean"
colnames(IASCAssess.C)[9] <- "Assess.C.Mean"

#Remove unnecessary columns
IASCAssess.A <- IASCAssess.A[-c(1:8)]
IASCAssess.C <- IASCAssess.C[-c(1:8)]

IASCmaster <- merge(IASCAssess.B, IASCAssess.A, by = "id", all.x=TRUE)
IASCmaster1 <- merge(IASCmaster, IASCAssess.C, by = "id", all.x=TRUE)

#Remove modifier questions
IASCmaster1 <- IASCmaster1[!(IASCmaster1$Question.No==1.4),]
IASCmaster1 <- IASCmaster1[!(IASCmaster1$Question.No==2.2),]
IASCmaster1 <- IASCmaster1[!(IASCmaster1$Question.No==3.3),]
IASCmaster1 <- IASCmaster1[!(IASCmaster1$Question.No==4.3),]

#Produce columns of differences in question scores within the dataframe
IASCmaster1$Assess.B.Assess.C <- IASCmaster1$Assess.B.Mean-IASCmaster1$Assess.C.Mean
IASCmaster1$Assess.B.Assess.A <- IASCmaster1$Assess.B.Mean-IASCmaster1$Assess.A.Mean
IASCmaster1$Assess.C.Assess.A <- IASCmaster1$Assess.C.Mean-IASCmaster1$Assess.A.Mean

## write a new file for checking in excel
write.csv(IASCmaster1, "IASC_check.csv")


####Create interassessor score comparisons for module scores
###NOTE Current issue where using an assessor who has not scored a particular species
###Leads to NA for categorical variables. If possible any assessor who has assessed all species
###should be Assessor A



#Separate by category
MS <- Final_df
MS$id <- paste(MS$Species,MS$Ecoregion.Assessed)
MSI <- MS[MS$Category=="Invasion",]
MSE <- MS[MS$Category=="Ecological",]
MSS <- MS[MS$Category=="Socio-Econ",]

##Invasion
#Separate data frames
MSIAssess.B <- MSI[MSI$Assessor=="Assess.B",]
MSIAssess.A <- MSI[MSI$Assessor=="Assess.A",]
MSIAssess.C <- MSI[MSI$Assessor=="Assess.C",]

#Rename columns
colnames(MSIAssess.B)[6] <- "Assess.B.MSI"
colnames(MSIAssess.A)[6] <- "Assess.A.MSI"
colnames(MSIAssess.C)[6] <- "Assess.C.MSI"

#Remove unnecessary columns
MSIAssess.B <- MSIAssess.B[-c(1:5,7)]
MSIAssess.C <- MSIAssess.C[-c(1:5,7)]
MSIAssess.A <- MSIAssess.A[-c(2,4,5,7)]

#Merge dataframes
MSImaster <- merge(MSIAssess.A, MSIAssess.B, by = "id", all.x=TRUE)
MSImaster1 <- merge(MSImaster, MSIAssess.C, by = "id", all.x=TRUE)

##Ecological
#Separate data frames
MSEAssess.B <- MSE[MSE$Assessor=="Assess.B",]
MSEAssess.A <- MSE[MSE$Assessor=="Assess.A",]
MSEAssess.C <- MSE[MSE$Assessor=="Assess.C",]

#Rename columns
colnames(MSEAssess.B)[6] <- "Assess.B.MSE"
colnames(MSEAssess.A)[6] <- "Assess.A.MSE"
colnames(MSEAssess.C)[6] <- "Assess.C.MSE"

#Remove unnecessary columns
MSEAssess.B <- MSEAssess.B[-c(1:5,7)]
MSEAssess.C <- MSEAssess.C[-c(1:5,7)]
MSEAssess.A <- MSEAssess.A[-c(2,4,5,7)]

#Merge dataframes
MSEmaster <- merge(MSEAssess.A, MSEAssess.B, by = "id", all.x=TRUE)
MSEmaster1 <- merge(MSEmaster, MSEAssess.C, by = "id", all.x=TRUE)
MSEmaster1 <- MSEmaster1[!duplicated(MSEmaster1$id),]

#Socioeconomic 
#Separate data frames
MSSAssess.B <- MSS[MSS$Assessor=="Assess.B",]
MSSAssess.A <- MSS[MSS$Assessor=="Assess.A",]
MSSAssess.C <- MSS[MSS$Assessor=="Assess.C",]

#Rename columns
colnames(MSSAssess.B)[6] <- "Assess.B.MSS"
colnames(MSSAssess.A)[6] <- "Assess.A.MSS"
colnames(MSSAssess.C)[6] <- "Assess.C.MSS"

#Remove unnecessary columns
MSSAssess.B <- MSSAssess.B[-c(1:5,7)]
MSSAssess.C <- MSSAssess.C[-c(1:5,7)]
MSSAssess.A <- MSSAssess.A[-c(2,4,5,7)]

#Merge dataframes
MSSmaster <- merge(MSSAssess.A, MSSAssess.B, by = "id", all.x=TRUE)
MSSmaster1 <- merge(MSSmaster, MSSAssess.C, by = "id", all.x=TRUE)

##Total Scores

#Merge module dataframes
TSmaster <- merge(MSEmaster1, MSImaster1, by = "id", all.x=TRUE)
TSmaster <- merge(TSmaster, MSSmaster1, by = "id", all.x=TRUE)



TSmaster$Assess.ATS <- TSmaster$Assess.A.MSI * (TSmaster$Assess.A.MSE +TSmaster$Assess.A.MSS)  
TSmaster$Assess.BTS <- TSmaster$Assess.B.MSI * (TSmaster$Assess.B.MSE +TSmaster$Assess.B.MSS)  
TSmaster$Assess.CTS <- TSmaster$Assess.C.MSI * (TSmaster$Assess.C.MSE +TSmaster$Assess.C.MSS)  



###Plotting of interassessor scores. This section of code provides four plots; one for each module
#comparing question scores and one comparing the total scores for the assessment. If multiple 
#ecoregions are assessed, use subset([dataframe], Ecoregion.Assessed %in% "Insert Category Here") 
#for data in ggplot. Use data from MSSmaster1, MSEmaster1, and MSImaster1 to 
##examine the module scores rather than the question scores.

pINV <- ggplot(IAImaster1)+
  geom_point(aes(x = Assess.B.Mean,
                 y = Assess.C.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.B.Mean,
                  y = Assess.C.Mean,
                  colour = "A-B"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.B.Mean,
                 y = Assess.A.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.B.Mean,
                  y = Assess.A.Mean,
                  colour = "A-C"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.C.Mean,
                 y = Assess.A.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.C.Mean,
                  y = Assess.A.Mean,
                  colour = "B-C"),
              method = lm, size = 0.75)+
  geom_abline()+
  ylab("") +
  xlab("") +
  labs(colour = "Assessor combination (1 vs 2)") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0.5, size = 10),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"));pINV

pECO <- ggplot(IAEcomaster1)+
  geom_point(aes(x = Assess.B.Mean,
                 y = Assess.C.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.B.Mean,
                  y = Assess.C.Mean,
                  colour = "A-B"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.B.Mean,
                 y = Assess.A.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.B.Mean,
                  y = Assess.A.Mean,
                  colour = "A-C"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.C.Mean,
                 y = Assess.A.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.C.Mean,
                  y = Assess.A.Mean,
                  colour = "B-C"),
              method = lm, size = 0.75)+
  geom_abline()+
  xlab("") +
  ylab("") +
  labs(colour = "Assessor combination (1 vs 2)") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0.5, size = 10),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_x_continuous(expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1))+
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1));pECO

pSoc <- ggplot(IASCmaster1)+
  geom_point(aes(x = Assess.B.Mean,
                 y = Assess.C.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.B.Mean,
                  y = Assess.C.Mean,
                  colour = "A-B"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.B.Mean,
                 y = Assess.A.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.B.Mean,
                  y = Assess.A.Mean,
                  colour = "A-C"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.C.Mean,
                 y = Assess.A.Mean),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.C.Mean,
                  y = Assess.A.Mean,
                  colour = "B-C"),
              method = lm, size = 0.75)+
  geom_abline()+
  xlab("") +
  ylab("") +
  labs(colour ="Assessor combination (1 vs 2)") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0.5, size = 10),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_x_continuous(expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1))+
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(1,2,3), 
                     limits = c(0.9,3.1));pSoc

pTotal <- ggplot(data = TSmaster)+
  geom_point(aes(x = Assess.BTS,
                 y = Assess.CTS),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.BTS,
                  y = Assess.CTS,
                  colour = "A-B"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.BTS,
                 y = Assess.ATS),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.BTS,
                  y = Assess.ATS,
                  colour = "A-C"),
              method = lm, size = 0.75)+
  geom_point(aes(x = Assess.CTS,
                 y = Assess.ATS),
             alpha = 0.2)+
  geom_smooth(aes(x = Assess.CTS,
                  y = Assess.ATS,
                  colour = "B-C"),
              method = lm, size = 0.75)+
  geom_abline()+
  xlab("") +
  ylab("") +
  labs(colour ="Assessor combination") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0.5, size = 10),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(colour = "white"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_x_continuous(expand = c(0,0), 
                     breaks = c(5,10,15,20,25), 
                     limits = c(1,26))+
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(5,10,15,20,25), 
                     limits = c(1,26));pTotal


ggarrange(pINV, pECO, pSoc, pTotal, common.legend = T, nrow=2, ncol=2, legend = "bottom")
ggsave("InterAssess.tiff", width = 9, height = 9.5, bg = "white")

#Correlation testing----
#Question scores by module (May need to be separated by ecoregion)
#Invasion
A.B <- lm(Assess.A.Mean ~ Assess.B.Mean, data = IAImaster1)
summary(A.B)
A.C <- lm(Assess.A.Mean ~ Assess.C.Mean, data = IAImaster1)
summary(A.C)
B.C <- lm(Assess.B.Mean ~ Assess.C.Mean, data = IAImaster1)
summary(B.C)

cor.test(IAImaster1$Assess.B.Mean, IAImaster1$Assess.A.Mean, method = "pearson")
cor.test(IAImaster1$Assess.C.Mean, IAImaster1$Assess.A.Mean, method = "pearson")
cor.test(IAImaster1$Assess.C.Mean, IAImaster1$Assess.B.Mean, method = "pearson")

#Ecological
A.B <- lm(Assess.A.Mean ~ Assess.B.Mean, data = IAEcomaster1)
summary(A.B)
A.C <- lm(Assess.A.Mean ~ Assess.C.Mean, data = IAEcomaster1)
summary(A.C)
B.C <- lm(Assess.B.Mean ~ Assess.C.Mean, data = IAEcomaster1)
summary(B.C)

cor.test(IAEcomaster1$Assess.B.Mean, IAEcomaster1$Assess.A.Mean, method = "pearson")
cor.test(IAEcomaster1$Assess.C.Mean, IAEcomaster1$Assess.A.Mean, method = "pearson")
cor.test(IAEcomaster1$Assess.C.Mean, IAEcomaster1$Assess.B.Mean, method = "pearson")

#Socioeconomic
A.B <- lm(Assess.A.Mean ~ Assess.B.Mean, data = IASCmaster1)
summary(A.B)
A.C <- lm(Assess.A.Mean ~ Assess.C.Mean, data = IASCmaster1)
summary(A.C)
B.C <- lm(Assess.B.Mean ~ Assess.C.Mean, data = IASCmaster1)
summary(B.C)

cor.test(IASCmaster1$Assess.B.Mean, IASCmaster1$Assess.A.Mean, method = "pearson")
cor.test(IASCmaster1$Assess.C.Mean, IASCmaster1$Assess.A.Mean, method = "pearson")
cor.test(IASCmaster1$Assess.C.Mean, IASCmaster1$Assess.B.Mean, method = "pearson")

#Correlation testing of Module scores
#Invasion
A.B <- lm(Assess.A.MSI ~ Assess.B.MSI, data = MSImaster1)
summary(A.B)
A.C <- lm(Assess.A.MSI ~ Assess.C.MSI, data = MSImaster1)
summary(A.C)
B.C <- lm(Assess.B.MSI ~ Assess.C.MSI, data = MSImaster1)
summary(B.C)

cor.test(MSImaster1$Assess.B.MSI, MSImaster1$Assess.A.MSI, method = "pearson")
cor.test(MSImaster1$Assess.C.MSI, MSImaster1$Assess.A.MSI, method = "pearson")
cor.test(MSImaster1$Assess.C.MSI, MSImaster1$Assess.B.MSI, method = "pearson")

#Ecological
A.B <- lm(Assess.A.MSE ~ Assess.B.MSE, data = MSEmaster1)
summary(A.B)
A.C <- lm(Assess.A.MSE ~ Assess.C.MSE, data = MSEmaster1)
summary(A.C)
B.C <- lm(Assess.B.MSE ~ Assess.C.MSE, data = MSEmaster1)
summary(B.C)

cor.test(MSEmaster1$Assess.B.MSE, MSEmaster1$Assess.A.MSE, method = "pearson")
cor.test(MSEmaster1$Assess.C.MSE, MSEmaster1$Assess.A.MSE, method = "pearson")
cor.test(MSEmaster1$Assess.C.MSE, MSEmaster1$Assess.B.MSE, method = "pearson")

#Socioeconomic
A.B <- lm(Assess.A.MSS ~ Assess.B.MSS, data = MSSmaster1)
summary(A.B)
A.C <- lm(Assess.A.MSS ~ Assess.C.MSS, data = MSSmaster1)
summary(A.C)
B.C <- lm(Assess.B.MSS ~ Assess.C.MSS, data = MSSmaster1)
summary(B.C)

cor.test(MSSmaster1$Assess.B.MSS, MSSmaster1$Assess.A.MSS, method = "pearson")
cor.test(MSSmaster1$Assess.C.MSS, MSSmaster1$Assess.A.MSS, method = "pearson")
cor.test(MSSmaster1$Assess.C.MSS, MSSmaster1$Assess.B.MSS, method = "pearson")

##Correlation on total score (may need to be separated by ecoregion)
cor.test(TSmaster$Assess.CTS, TSmaster$Assess.ATS, method = "pearson")
cor.test(TSmaster$Assess.BTS, TSmaster$Assess.ATS, method = "pearson")
cor.test(TSmaster$Assess.CTS, TSmaster$Assess.BTS, method = "pearson")




