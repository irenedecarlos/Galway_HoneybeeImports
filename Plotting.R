load(file="Selection1_import0.1.RData")
getwd()

locationsDF <- data.frame(Location = getLocation(c(age0$Mel, age1$Mel), collapse = TRUE),
                          Beekeeper = c(rep("Beekeeper1", nColonies(age0$Mel)),
                                        rep("Beekeeper2", nColonies(age1$Mel))))
ggplot(data = locationsDF, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
  geom_point()
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
setwd("C:/Users/Irene/Desktop/Galway/Results/test160224")
resultsdf<- read.csv("results.csv")
resultsdf2<- read.csv("results20.csv")
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
  geom_line(data = mel, aes(x = Year, y = meanIBD, color = "Mel 20"), linewidth=1.15) +
  geom_ribbon(data = mel, aes(x = Year, ymin = quanIBDl, ymax = quanIBDh, fill = "Mel 20"), alpha = 0.05, color = NA) +
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
resultsdf2 <- read.csv("results20.csv")

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
resultsdf20 <- read.csv("results20.csv")

# Calculate summary statistics for fitness
nuevamelf <- resultsdf1 %>%
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
  geom_line(data = melf, aes(color = "Mel -0.7"), linewidth = 1.15) +
  geom_ribbon(data = melf, aes(ymin = quanfl, ymax = quanfh, fill = "Mel -0.7"), alpha = 0.05, color = NA) +
  labs(title = "fitness over time",
       x = "Years",
       y = "Fitness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("red", "green", "#286ceb", "purple"))

############################################################################################################################
# plotting for spatial

locationsDF <- data.frame(Location = getLocation(c(age1$Mel), collapse = TRUE),
                          Beekeeper = c(rep("Beekeeper1", nColonies(age1$Mel))))
ggplot(data = locationsDF, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
  geom_point()




# getting the number of carnica queens that collapse each winter
resultsdf <- read.csv("results20.csv")
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