load(file="Selection3_import0.1.RData")
getwd()

locationsDF <- data.frame(Location = getLocation(c(age0$Mel, age1$Mel), collapse = TRUE),
                          Beekeeper = c(rep("Beekeeper1", nColonies(age0$Mel)),
                                        rep("Beekeeper2", nColonies(age1$Mel))))
ggplot(data = locationsDF, aes(x = Location.1, y = Location.2, colour = Beekeeper)) +  geom_point()
############################################################################################################################

#to plot IBD
library(AlphaSimR)
library(ggplot2)
library(tictoc)
library(R6)
library(nadiv)
library(Matrix)
library(SIMplyBee)
library(dplyr)
library(ggpubr)
#If we have them on separate files
getwd()
setwd("C:/Users/Irene/Desktop/Galway/Code/lugh/tests/stopimports")
resultsdf<- read.csv("results.csv")
resultsdf2<- read.csv("results04.csv")
resultsdf<-rbind(MeanVarMel,MeanVarCar)
nuevamel<- resultsdf %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                              quanIBDl= quantile(MeanIBD,p=0.025),
                                                              quanIBDh= quantile(MeanIBD,p=0.975))
Combined<-resultsdf2 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                 quanIBDl= quantile(MeanIBD,p=0.025),
                                                                 quanIBDh= quantile(MeanIBD,p=0.975))
mel <- Combined %>%
  filter(Population == "Mel")
pel <- nuevamel %>%
  filter(Population == "Mel")
ggplot(nuevamel, aes(x = Year, y = meanIBD, color = Population),linewidth=1.15) +
  geom_line(linewidth=1.15) +
  geom_ribbon(aes(ymin = quanIBDl, 
                  ymax = quanIBDh, 
                  fill = Population), alpha=0.05, color = NA) +
  #geom_line(data = mel, aes(x = Year, y = meanIBD, color = "Mel 20"), linewidth=1.15) +
  #geom_ribbon(data = mel, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel 20"), alpha = 0.05, color = NA) +
  labs(title = "Introgression over time",
       x = "Years",
       y = "pMelHaplo") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c(  "black","#FC4E07", "#286ceb","purple"))

############################################################################################################################
#Para plotear la miel (aplicable a fitness)

resultsdf <- read.csv("results.csv")
resultsdf2 <- read.csv("results04.csv")

# Calculate summary statistics for honey yield
nuevamel <- resultsdf %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
            quanIBDl = quantile(MeanIBD, p = 0.025),
            quanIBDh = quantile(MeanIBD, p = 0.975),
            honey_yield = ifelse(Population %in% c("Mel", "MelnI"), mean(HoneyYieldBrit), mean(HoneyYieldEu)),
          quanhyl = ifelse(Population %in% c("Mel", "MelnI"), quantile(HoneyYieldBrit, p = 0.025), quantile(HoneyYieldEu, p = 0.025)),
          quanhyh = ifelse(Population %in% c("Mel", "MelnI"), quantile(HoneyYieldBrit, p = 0.975), quantile(HoneyYieldEu, p = 0.975)))

Combined <- resultsdf2 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
            quanIBDl = quantile(MeanIBD, p = 0.025),
            quanIBDh = quantile(MeanIBD, p = 0.975),
            honey_yield = ifelse(Population == "Car", mean(HoneyYieldEu), mean(HoneyYieldBrit)),
          quanhyl = ifelse(Population =="Car", quantile(HoneyYieldEu, p = 0.025), quantile(HoneyYieldBrit, p = 0.025)),
          quanhyh = ifelse(Population  =="Car", quantile(HoneyYieldEu, p = 0.975), quantile(HoneyYieldBrit, p = 0.975)))

# Filter honey yield for Mel and MelnI populations
mel <- Combined %>%
  filter(Population %in% "Mel")

# Plotting
ggplot(nuevamel, aes(x = Year, y = honey_yield, color = Population)) +
  geom_line(linewidth = 1.15) +
  geom_ribbon(aes(ymin = quanhyl, ymax = quanhyh, fill = Population), alpha = 0.05, color = NA) +
  geom_line(data = mel, aes(color = "Mel 20"), linewidth = 1.15) +
  geom_ribbon(data = mel, aes(ymin = quanhyl, ymax = quanhyh, fill = "Mel 20"), alpha = 0.05, color = NA) +
  labs(title = "Honeyyield over time",
       x = "Years",
       y = "Honey Yield") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("red", "green", "#286ceb", "purple"))


################################################################################################################################################
#Fitness
resultsdf1 <- read.csv("results.csv")
resultsdf20 <- read.csv("results04.csv")

# Calculate summary statistics for fitness
nuevamelf <- resultsdf %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population %in% c("Mel", "MelnI"), mean(FitnessBrit), mean(FitnessEu)),
          quanfl = ifelse(Population %in% c("Mel", "MelnI"), quantile(FitnessBrit, p = 0.025), quantile(FitnessEu, p = 0.025)),
          quanfh = ifelse(Population %in% c("Mel", "MelnI"), quantile(FitnessBrit, p = 0.975), quantile(FitnessEu, p = 0.975)))

Combinedf <- resultsdf20 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population == "Car", mean(FitnessEu), mean(FitnessBrit)),
          quanfl = ifelse(Population =="Car", quantile(FitnessEu, p = 0.025), quantile(FitnessBrit, p = 0.025)),
          quanfh = ifelse(Population  =="Car", quantile(FitnessEu, p = 0.975), quantile(FitnessBrit, p = 0.975)))

# Filterfitness for Mel and MelnI populations
melf <- Combinedf %>%
  filter(Population %in% "Mel")

# Plotting
ggplot(nuevamelf, aes(x = Year, y = fitness, color = Population)) +
  geom_line(linewidth = 1.15) +
  geom_ribbon(aes(ymin = quanfl, ymax = quanfh, fill = Population), alpha = 0.05, color = NA) +
 # geom_line(data = melf, aes(color = "Mel -0.7"), linewidth = 1.15) +
 # geom_ribbon(data = melf, aes(ymin = quanfl, ymax = quanfh, fill = "Mel -0.7"), alpha = 0.05, color = NA) +
  labs(title = "fitness over time",
       x = "Years",
       y = "Fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("red", "green", "#286ceb", "purple"))

############################################################################################################################
# plotting for spatial

locationsDF <- data.frame(Location = getLocation(c(age0$Mel), collapse = TRUE),
                          Beekeeper = c(rep("Beekeeper1", nColonies(age0$Mel))))
ggplot(data = locationsDF, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
  geom_point()




# getting the number of carnica queens that collapse each winter
resultsdf <- read.csv("resultscor07.csv")
nuevamelf <- resultsdf %>%
  group_by(Year, Population) %>%
  reframe(survivingCar = mean(survivingCar))
melf <- nuevamelf %>%
  filter(Population %in% "Mel")
ggplot(melf, aes(x = as.factor(Year), y = survivingCar)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Year", y = "Number of Surviving Queens", fill = "Replica") +
  ggtitle("Number of Surviving Queens Over the Years") +
  theme_minimal()
melf

#If we have them in just one file
#para tenerlos todos en una lista hay que cambiar la forma de guardarlo en simulation.R para que los csv se guarden
#como ScenarioData.csv y luego cambiar analyze.R para que lea los ScenarioData.csv y los ponga en una tabla unica
#
resultsdf2<- read.csv("results2.csv")
summary_df <- results %>%
  group_by(Scenario,Year, Population) %>%
  summarise(meanIBD = mean(MeanIBD),
            quanIBDl= quantile(MeanIBD,p=0.025),
            quanIBDh= quantile(MeanIBD,p=0.975))

# Filter the summary data for plotting
plot_data <- summary_df %>%
  filter((Scenario == "name of scenario" & Population %in% c("Mel", "Car", "MelnI")) |
           (Scenario == "name of scenario2" & Population == "Mel"))

# Plotting using ggplot2
ggplot(plot_data, aes(x = Year, y = Weight, color = Population, linetype = Scenario)) +
  geom_line() +
  labs(title = "Weight per Year per Population per Scenario",
       x = "Year", y = "Weight") +
  theme_minimal()        
##########################################################################################################################################################
#plots for the correlation
getwd()
setwd("C:/Users/Irene/Desktop/Galway/Results/scenario2/scenario2gv")
resultsdf<- read.csv("resultscor25.csv")
resultsdf2<- read.csv("resultscor07.csv")
resultsdf3<- read.csv("resultscor-07.csv")
nuevamel<- resultsdf %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                quanIBDl= quantile(MeanIBD,p=0.025),
                                                                quanIBDh= quantile(MeanIBD,p=0.975))
Combined<-resultsdf2 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                quanIBDl= quantile(MeanIBD,p=0.025),
                                                                quanIBDh= quantile(MeanIBD,p=0.975))
cal<-resultsdf3 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                quanIBDl= quantile(MeanIBD,p=0.025),
                                                                quanIBDh= quantile(MeanIBD,p=0.975))
melc <- cal %>%
  filter(Population == "Mel")
mel <- Combined %>%
  filter(Population == "Mel")
pel <- nuevamel %>%
  filter(Population == "Mel")
ggplot(nuevamel, aes(x = Year, y = meanIBD, color = Population),linewidth=1.15) +
  geom_line(linewidth=1.15) +
  geom_ribbon(aes(ymin = quanIBDl, 
                  ymax = quanIBDh, 
                  fill = Population), alpha=0.05, color = NA) +
  geom_line(data = mel, aes(x = Year, y = meanIBD, color = "Mel 0.7"), linewidth=1.15) +
  geom_ribbon(data = mel, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel 0.7"), alpha = 0.05, color = NA) +
  geom_line(data = melc, aes(x = Year, y = meanIBD, color = "Mel -0.7"), linewidth=1.15) +
  geom_ribbon(data = melc, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel -0.7"), alpha = 0.05, color = NA) +
  labs(title = "Introgression over time",
       x = "Years",
       y = "pMelHaplo") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c(  "#FC4E07","green", "#286ceb","purple"))
getgv





#plots for scenario1
getwd()
setwd("C:/Users/Irene/Desktop/Galway/Results/quickhaploscen1")
resultsdf<- read.csv("results.csv")
resultsdf2<- read.csv("results04.csv")
resultsdf3<- read.csv("results015.csv")
resultsdf4<- read.csv("results0.csv")

nuevamel<- resultsdf %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                quanIBDl= quantile(MeanIBD,p=0.025),
                                                                quanIBDh= quantile(MeanIBD,p=0.975))
Combined<-resultsdf2 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                quanIBDl= quantile(MeanIBD,p=0.025),
                                                                quanIBDh= quantile(MeanIBD,p=0.975))
cal<-resultsdf3 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                           quanIBDl= quantile(MeanIBD,p=0.025),
                                                           quanIBDh= quantile(MeanIBD,p=0.975))
cero<-resultsdf4 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                           quanIBDl= quantile(MeanIBD,p=0.025),
                                                           quanIBDh= quantile(MeanIBD,p=0.975))

mel <- Combined %>%
  filter(Population == "Mel")
melc <- cal %>%
  filter(Population == "Mel")
pel <- cero %>%
  filter(Population == "Mel")
ggplot(nuevamel, aes(x = Year, y = meanIBD, color = Population),linewidth=1.15) +
  geom_line(linewidth=1.15) +
  geom_ribbon(aes(ymin = quanIBDl, 
                  ymax = quanIBDh, 
                  fill = Population), alpha=0.05, color = NA) +
  geom_line(data = mel, aes(x = Year, y = meanIBD, color = "Mel 10%"), linewidth=1.15) +
  geom_ribbon(data = mel, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel 4%"), alpha = 0.05, color = NA) +
  geom_line(data = melc, aes(x = Year, y = meanIBD, color = "Mel 1.5%"), linewidth=1.15) +
  geom_ribbon(data = melc, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel 1.5%"), alpha = 0.05, color = NA) +
  geom_line(data = pel, aes(x = Year, y = meanIBD, color = "Mel 0%"), linewidth=1.15) +
  geom_ribbon(data = pel, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel 0%"), alpha = 0.05, color = NA) +
  labs(title = "Introgression over time",
       x = "Years",
       y = "pMelHaplo") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c(  "#FC4E07","green", "#286ceb","purple","grey"))

#Fitness
resultsdf<- read.csv("results.csv")
resultsdf2<- read.csv("results05.csv")
resultsdf3<- read.csv("results015.csv")
resultsdf4<- read.csv("results0.csv")

# Calculate summary statistics for fitness
nuevamelf <- resultsdf %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population %in% c("Mel"), mean(FitnessBrit), mean(FitnessEu)),
          quanfl = ifelse(Population %in% c("Mel"), quantile(FitnessBrit, p = 0.025), quantile(FitnessEu, p = 0.025)),
          quanfh = ifelse(Population %in% c("Mel"), quantile(FitnessBrit, p = 0.975), quantile(FitnessEu, p = 0.975)))

Combinedf <- resultsdf2 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population == "Car", mean(FitnessEu), mean(FitnessBrit)),
          quanfl = ifelse(Population =="Car", quantile(FitnessEu, p = 0.025), quantile(FitnessBrit, p = 0.025)),
          quanfh = ifelse(Population  =="Car", quantile(FitnessEu, p = 0.975), quantile(FitnessBrit, p = 0.975)))

Combinedf1 <- resultsdf3 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population == "Car", mean(FitnessEu), mean(FitnessBrit)),
          quanfl = ifelse(Population =="Car", quantile(FitnessEu, p = 0.025), quantile(FitnessBrit, p = 0.025)),
          quanfh = ifelse(Population  =="Car", quantile(FitnessEu, p = 0.975), quantile(FitnessBrit, p = 0.975)))



# Filterfitness for Mel and MelnI populations
melf <- Combinedf %>%
  filter(Population %in% "Mel")
melf2 <- Combinedf1 %>%
  filter(Population %in% "Mel")
# Plotting
ggplot(nuevamelf, aes(x = Year, y = fitness, color = Population)) +
  geom_line(linewidth = 1.15) +
  geom_ribbon(aes(ymin = quanfl, ymax = quanfh, fill = Population), alpha = 0.05, color = NA) +
  geom_line(data = melf, aes(color = "Mel 4%"), linewidth = 1.15) +
  geom_ribbon(data = melf, aes(ymin = quanfl, ymax = quanfh, fill = "Mel 4%"), alpha = 0.05, color = NA) +
  geom_line(data = melf2, aes(color = "Mel 1.5%"), linewidth = 1.15) +
  geom_ribbon(data = melf2, aes(ymin = quanfl, ymax = quanfh, fill = "Mel 1.5%"), alpha = 0.05, color = NA) +
  labs(title = "fitness over time",
       x = "Years",
       y = "Fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("red", "green", "#286ceb", "purple"))

#HY
resultsdf<- read.csv("results.csv")
resultsdf2<- read.csv("results04.csv")
resultsdf3<- read.csv("results015.csv")
resultsdf4<- read.csv("results0.csv")
# Calculate summary statistics for fitness
nuevamelf <- resultsdf %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          Honeyyield = ifelse(Population %in% c("Mel"), mean(HoneyYieldBrit), mean(HoneyYieldEu)),
          quanhl = ifelse(Population %in% c("Mel"), quantile(HoneyYieldBrit, p = 0.025), quantile(HoneyYieldEu, p = 0.025)),
          quanhh = ifelse(Population %in% c("Mel"), quantile(HoneyYieldBrit, p = 0.975), quantile(HoneyYieldEu, p = 0.975)))

Combinedf <- resultsdf2 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          Honeyyield = ifelse(Population == "Car", mean(HoneyYieldEu), mean(HoneyYieldBrit)),
          quanhl = ifelse(Population =="Car", quantile(HoneyYieldEu, p = 0.025), quantile(HoneyYieldBrit, p = 0.025)),
          quanhh = ifelse(Population  =="Car", quantile(HoneyYieldEu, p = 0.975), quantile(HoneyYieldBrit, p = 0.975)))

Combinedf1 <- resultsdf3 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          Honeyyield = ifelse(Population == "Car", mean(HoneyYieldEu), mean(HoneyYieldBrit)),
          quanhl = ifelse(Population =="Car", quantile(HoneyYieldEu, p = 0.025), quantile(HoneyYieldBrit, p = 0.025)),
          quanhh = ifelse(Population  =="Car", quantile(HoneyYieldEu, p = 0.975), quantile(HoneyYieldBrit, p = 0.975)))



# Filterfitness for Mel and MelnI populations
melf <- Combinedf %>%
  filter(Population %in% "Mel")
melf2 <- Combinedf1 %>%
  filter(Population %in% "Mel")
# Plotting
ggplot(nuevamelf, aes(x = Year, y = Honeyyield, color = Population)) +
  geom_line(linewidth = 1.15) +
  geom_ribbon(aes(ymin = quanhl, ymax = quanhh, fill = Population), alpha = 0.05, color = NA) +
  geom_line(data = melf, aes(color = "Mel 04"), linewidth = 1.15) +
  geom_ribbon(data = melf, aes(ymin = quanhl, ymax = quanhh, fill = "Mel 04"), alpha = 0.05, color = NA) +
  geom_line(data = melf2, aes(color = "Mel 015"), linewidth = 1.15) +
  geom_ribbon(data = melf2, aes(ymin = quanhl, ymax = quanhh, fill = "Mel 015"), alpha = 0.05, color = NA) +
  labs(title = "fitness over time",
       x = "Years",
       y = "Fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("red", "green", "#286ceb", "purple"))









#plots for scenario2
getwd()
setwd("C:/Users/Irene/Desktop/Galway/Results/scenario2")
resultsdf<- read.csv("resultscor25.csv")
resultsdf2<- read.csv("results07.csv")
resultsdf3<- read.csv("results-07.csv")
nuevamel<- resultsdf %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                quanIBDl= quantile(MeanIBD,p=0.025),
                                                                quanIBDh= quantile(MeanIBD,p=0.975))
Combined<-resultsdf2 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                                quanIBDl= quantile(MeanIBD,p=0.025),
                                                                quanIBDh= quantile(MeanIBD,p=0.975))
cal<-resultsdf3 %>% group_by(Year,Population)%>% summarise(meanIBD = mean(MeanIBD),
                                                           quanIBDl= quantile(MeanIBD,p=0.025),
                                                           quanIBDh= quantile(MeanIBD,p=0.975))

mel <- Combined %>%
  filter(Population == "Mel")
melc <- cal %>%
  filter(Population == "Mel")

ggplot(nuevamel, aes(x = Year, y = meanIBD, color = Population),linewidth=1.15) +
  geom_line(linewidth=1.15) +
  geom_ribbon(aes(ymin = quanIBDl, 
                  ymax = quanIBDh, 
                  fill = Population), alpha=0.05, color = NA) +
  geom_line(data = mel, aes(x = Year, y = meanIBD, color = "Mel 0.7"), linewidth=1.15) +
  geom_ribbon(data = mel, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel 0.7"), alpha = 0.05, color = NA) +
  geom_line(data = melc, aes(x = Year, y = meanIBD, color = "Mel -0.7"), linewidth=1.15) +
  geom_ribbon(data = melc, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel -0.7"), alpha = 0.05, color = NA) +
  labs(title = "Introgression over time",
       x = "Years",
       y = "pMelHaplo") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c(  "#FC4E07","green", "#286ceb","purple"))

#Fitness
resultsdf<- read.csv("resultscor25.csv")
resultsdf2<- read.csv("results07.csv")
resultsdf3<- read.csv("results-07.csv")

# Calculate summary statistics for fitness
nuevamelf <- resultsdf %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population %in% c("Mel"), mean(FitnessBrit), mean(FitnessEu)),
          quanfl = ifelse(Population %in% c("Mel"), quantile(FitnessBrit, p = 0.025), quantile(FitnessEu, p = 0.025)),
          quanfh = ifelse(Population %in% c("Mel"), quantile(FitnessBrit, p = 0.975), quantile(FitnessEu, p = 0.975)))

Combinedf <- resultsdf2 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population == "Car", mean(FitnessEu), mean(FitnessBrit)),
          quanfl = ifelse(Population =="Car", quantile(FitnessEu, p = 0.025), quantile(FitnessBrit, p = 0.025)),
          quanfh = ifelse(Population  =="Car", quantile(FitnessEu, p = 0.975), quantile(FitnessBrit, p = 0.975)))

Combinedf1 <- resultsdf3 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          fitness = ifelse(Population == "Car", mean(FitnessEu), mean(FitnessBrit)),
          quanfl = ifelse(Population =="Car", quantile(FitnessEu, p = 0.025), quantile(FitnessBrit, p = 0.025)),
          quanfh = ifelse(Population  =="Car", quantile(FitnessEu, p = 0.975), quantile(FitnessBrit, p = 0.975)))



# Filterfitness for Mel and MelnI populations
melf <- Combinedf %>%
  filter(Population %in% "Mel")
melf2 <- Combinedf1 %>%
  filter(Population %in% "Mel")
# Plotting
ggplot(nuevamelf, aes(x = Year, y = fitness, color = Population)) +
  geom_line(linewidth = 1.15) +
  geom_ribbon(aes(ymin = quanfl, ymax = quanfh, fill = Population), alpha = 0.05, color = NA) +
  geom_line(data = melf, aes(color = "Mel 0.7"), linewidth = 1.15) +
  geom_ribbon(data = melf, aes(ymin = quanfl, ymax = quanfh, fill = "Mel 0.7"), alpha = 0.05, color = NA) +
  geom_line(data = melf2, aes(color = "Mel -0.7"), linewidth = 1.15) +
  geom_ribbon(data = melf2, aes(ymin = quanfl, ymax = quanfh, fill = "Mel -0.7"), alpha = 0.05, color = NA) +
  labs(title = "fitness over time",
       x = "Years",
       y = "Fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("red", "green", "#286ceb", "purple"))

#HY
resultsdf<- read.csv("resultscor25.csv")
resultsdf2<- read.csv("results07.csv")
resultsdf3<- read.csv("results-07.csv")

# Calculate summary statistics for fitness
nuevamelf <- resultsdf %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          Honeyyield = ifelse(Population %in% c("Mel"), mean(HoneyYieldBrit), mean(HoneyYieldEu)),
          quanhl = ifelse(Population %in% c("Mel"), quantile(HoneyYieldBrit, p = 0.025), quantile(HoneyYieldEu, p = 0.025)),
          quanhh = ifelse(Population %in% c("Mel"), quantile(HoneyYieldBrit, p = 0.975), quantile(HoneyYieldEu, p = 0.975)))

Combinedf <- resultsdf2 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          Honeyyield = ifelse(Population == "Car", mean(HoneyYieldEu), mean(HoneyYieldBrit)),
          quanhl = ifelse(Population =="Car", quantile(HoneyYieldEu, p = 0.025), quantile(HoneyYieldBrit, p = 0.025)),
          quanhh = ifelse(Population  =="Car", quantile(HoneyYieldEu, p = 0.975), quantile(HoneyYieldBrit, p = 0.975)))

Combinedf1 <- resultsdf3 %>%
  group_by(Year, Population) %>%
  reframe(meanIBD = mean(MeanIBD),
          quanIBDl = quantile(MeanIBD, p = 0.025),
          quanIBDh = quantile(MeanIBD, p = 0.975),
          Honeyyield = ifelse(Population == "Car", mean(HoneyYieldEu), mean(HoneyYieldBrit)),
          quanhl = ifelse(Population =="Car", quantile(HoneyYieldEu, p = 0.025), quantile(HoneyYieldBrit, p = 0.025)),
          quanhh = ifelse(Population  =="Car", quantile(HoneyYieldEu, p = 0.975), quantile(HoneyYieldBrit, p = 0.975)))



# Filterfitness for Mel and MelnI populations
melf <- Combinedf %>%
  filter(Population %in% "Mel")
melf2 <- Combinedf1 %>%
  filter(Population %in% "Mel")
# Plotting
ggplot(nuevamelf, aes(x = Year, y = Honeyyield, color = Population)) +
  geom_line(linewidth = 1.15) +
  geom_ribbon(aes(ymin = quanhl, ymax = quanhh, fill = Population), alpha = 0.05, color = NA) +
  geom_line(data = melf, aes(color = "Mel 0.7"), linewidth = 1.15) +
  geom_ribbon(data = melf, aes(ymin = quanhl, ymax = quanhh, fill = "Mel 0.7"), alpha = 0.05, color = NA) +
  geom_line(data = melf2, aes(color = "Mel -0.7"), linewidth = 1.15) +
  geom_ribbon(data = melf2, aes(ymin = quanhl, ymax = quanhh, fill = "Mel -0.7"), alpha = 0.05, color = NA) +
  labs(title = "fitness over time",
       x = "Years",
       y = "Fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("red", "green", "#286ceb", "purple"))


#############################PLOTING_EUCLIDEAN#######################################################################

setwd("C:/Users/Irene/Desktop/Galway/Code/lugh/tests/dublin")

resultsdf<- read.csv("Eucdist_rep1_import0.1.csv")
resultsdf2<- read.csv("Eucdist_rep2_import0.1.csv")
resultsdf3<- read.csv("Eucdist_rep3_import0.1.csv")
num_rows <- nrow(resultsdf)
resultsdf$Rep <- rep(1, num_rows)

num_rows2 <- nrow(resultsdf2)
resultsdf2$Rep <- rep(2, num_rows)

num_rows3 <- nrow(resultsdf3)
resultsdf3$Rep <- rep(3, num_rows)

#bigdf<-rbind(resultsdf,resultsdf2)

sorted_df <- resultsdf %>%
  arrange(eucyear, EuclAge)
sorted_df2 <- resultsdf2 %>%
  arrange(eucyear, EuclAge)
sorted_df3 <- resultsdf3 %>%
  arrange(eucyear, EuclAge)

sorted_df$Introgression <- 1- sorted_df$IBD
sorted_df2$Introgression <- 1- sorted_df2$IBD
sorted_df3$Introgression <- 1- sorted_df3$IBD

newdf<-rowMeans(abind::abind(sorted_df, sorted_df2, sorted_df3, along = 3), dims=2)
newdf<-as.data.frame(newdf)
au<-newdf[1801:2000,]
aw<-newdf[801:1000,]
av<-newdf[1:200,]
a3<-newdf[601:800,]
a20<-newdf[3801:4000,]
#hay que hacerlo over the years
ggplot(a20, aes(x = EuclAge, y = Introgression, color = Introgression)) +
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  geom_smooth(method = "lm", se = FALSE) +  # Add trend line (linear regression)+
  labs(title = "Introgression vs. Euclidean Distance",
       x = "Euclidean Distance", y = "IBD") +
  theme_minimal()

model_0_years <- lm(Introgression ~ EuclAge, data = av)
model_3_years <- lm(Introgression ~ EuclAge, data = a3)
model_5_years <- lm(Introgression ~ EuclAge, data = aw)
model_10_years <- lm(Introgression ~ EuclAge, data = au)
model_20_years <- lm(Introgression ~ EuclAge, data = a20)

# Extract slope coefficients
slope_0_years <- coef(model_0_years)[2]
slope_3_years <- coef(model_3_years)[2]
slope_5_years <- coef(model_5_years)[2]
slope_10_years <- coef(model_10_years)[2]
slope_20_years <- coef(model_20_years)[2]

Euclideandistance <- ahaha %>%
  group_by(eucyear) %>%
  reframe(EuclAge=EuclAge,
          IBD=IBD,
          IdAge=IdAge,)
ggplot(Euclideandistance, aes(x = EuclAge, y = IBD, color = IBD)) +
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Scatter Plot of IBD vs. Euclidean Distance",
       x = "Euclidean Distance", y = "IBD") +
  theme_minimal()
ggplot(Euclideandistance, aes(x = EuclAge, y = IBD)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Smoothed Curve of IBD vs. Euclidean Distance",
       x = "Euclidean Distance", y = "IBD") +
  theme_minimal()
