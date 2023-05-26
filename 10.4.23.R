install.packages('usethis')
library(usethis)
library(ggplot2)
use_git_config(user.name="Natlehhh", user.email= "natalieroberts2031@gmail.com")

library(lme4)
library(multcomp)

# hypothesis 1

prey_ab_vol <- read.csv('Dataa/prey_ab_vol.csv')
prey_ab_vol$row_id <- factor(1:nrow(prey_ab_vol))
Summaryprey_ab_vol <- summary(prey_ab_vol)
print(Summaryprey_ab_vol)

sd_N.H..perlevis <- sd(prey_ab_vol$N.H..perlevis)
sd_Volumen.Hym..mL. <- sd(prey_ab_vol$Volumen.Hym..mL.)
SDresults <- paste("Standard Deviation (N.H..perlevis):", sd_N.H..perlevis, "Standard Deviation (Volumen.Hym..mL.):", sd_Volumen.Hym..mL.)
print(SDresults)
  

boxplot(sqrt(N.H..perlevis)~ Date,prey_ab_vol,col="lemonchiffon", main = "Relationship between abundance and month",cex = 0.2)

# mixed model on counts
m1 <- glm(N.H..perlevis ~ Date,prey_ab_vol,family='poisson')
m1 <- glmer(N.H..perlevis ~ Date + (1|row_id),prey_ab_vol,family='poisson')
drop1(m1,test='Chisq')
summary(glht(m1, linfct = mcp(Date = "Tukey")))
(m1)
# abundance and volume correlated (weakly)
cor.test(log(prey_ab_vol$N.H..perlevis),log(prey_ab_vol$Volumen.Hym..mL.))

# volumetric analysis
m1 <- aov(log(Volumen.Hym..mL.) ~ Date,prey_ab_vol)
anova(m1)
TukeyHSD(m1,ordered=TRUE)
boxplot(log(Volumen.Hym..mL.)~Date,prey_ab_vol,main = "Relationship between volume and month",col= 'pink')


library(car)
leveneTest(m1)
leveneTest(log(Volumen.Hym..mL.) ~ factor(Date), data = prey_ab_vol)

residuals <- resid(m1, type = "pearson")
qqnorm(residuals)
qqline(residuals, distribution = qpois, lambda = mean(residuals))


hist(prey_ab_vol$Volumen.Hym..mL., 
     breaks = "FD",
     freq = FALSE,
     main = "Histogram of Dependent Variable",
     xlab = "Dependent Variable",
     ylab = "Density")




# homogeneity of variance assumption upheld


# new hypothesis 2
#abundance vs temp,salinity and pH, ened to average temperature, salinity and pH values

TpHSal <- read.csv('Dataa/TempsalinitypH.csv')
TpHSal$row_id <- factor(1:nrow(TpHSal))
summaryTpHSal <- summary(TpHSal)
print(summaryTpHSal)

sdT <- sd(TpHSal$T)
sdpH <- sd(TpHSal$PH)
Sal <- sd(TpHSal$S)
SDresults2 <- paste("Standard Deviation (T):", sdT, "Standard Deviation (pH):", sdpH, "Standard Deviation (Sal):", Sal)
print(SDresults2)


library(dplyr)


# Cleaning up data
data <- data.frame(
  ID_SITE = c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
  DATE = as.Date(c("10/6/2020", "18/6/2020", "26/6/2020", "5/11/2020", "12/11/2020", "20/11/2020", "27/11/2020", "3/4/2021", "8/4/2021", "16/4/2021", "7/6/2021", "15/6/2021", "21/6/2021", "29/6/2021"), format = "%d/%m/%Y"),
  T = c(15.416, 15.892, 16.824, 13.417, 13.241, 13.048, 12.848, 11.072, 11.23, 10.971, 15.061, 16.295, 15.688, 15.399),
  S = c(34.733, 34.632, 34.704, 34.552, 34.424, 34.56, 34.31, 33.62, 33.8, 34.156, 34.626, 34.762, 34.865, 35.014),
  PH = c(8.138, 7.959, 7.953, 7.799, 8.055, 7.791, 7.794, 7.9, 7.93, 7.91, 7.94, 7.91, 7.93, 7.88)
)

data$DATE <- as.POSIXlt(data$DATE)

# Extract month and year from the DATE column
data$Month <- data$DATE$mon + 1
data$Year <- data$DATE$year + 1900

# Calculate the average values for each month
averagestemp_ph_S <- aggregate(cbind(T, S, PH) ~ Month + Year, data, FUN = mean)

#sometimes have to rerun this code to get june20 and june21 seperated????
june_2020 <- subset(averagestemp_ph_S, Month == 6 & Year == 2020)
june_2021 <- subset(averagestemp_ph_S, Month == 6 & Year == 2021)

# Print the average values for each month.
print(averagestemp_ph_S)

#checked in excel, to make sure averages are right, they are (:
# making averages into values that can be compared

#comparing averages per month to temperature

# Calculate the average of Volumen.Hym..mL. divided by abundance for each date = Volumen.Hym..mL. for each sponge
averages2 <- aggregate(Volumen.Hym..mL. / N.H..perlevis ~ Date, prey_ab_vol, FUN = mean)
print(averages2)


#making dataset for sponge volume the same as the averages for pH,temp & salinity

averagesvolume <- data.frame(
  Month = c(6, 11, 4, 6),
  Year = c(2020, 2020, 2021, 2021),
  Volumen.Hym..mL..N.H..perlevis = c(3.8781924, 12.6688840, 3.0448443, 0.8122513)
)

#now can compare sponge volume vs temp,pH and salinity
# averagesvolume = average volume of H.perlevis for the months
#averagestemp_ph_S = temperature,ph and salinity values for the months

merged_dataset <- merge(averagesvolume, averagestemp_ph_S, by = c('Year', 'Month'))

# correlation analysis
correlation2 <- cor(merged_dataset[, c("Volumen.Hym..mL..N.H..perlevis", "T", "S", "PH")])
correlation2

residualsc2 <- residuals(correlation2)
qqnormm0 <- qqnorm(residualsc2)
qqnormm0 <- qqline(residualsc2)

# Hypotheis 3 #newdataset

#relationship between biomass and abundance per habitat
data2022 <- read.csv('Dataa/data2022.csv')
data2022$row_id <- factor(1:nrow(data2022))
data2022 <- data2022[data2022$Abundance.sponges > 0, ]
data2022$Sponge.biomass <- as.numeric(gsub("[^0-9.]", "", data2022$Sponge.biomass))

summarydata2022 <- summary(data2022)
print(summarydata2022)


data2022$log.n <- log(data2022$Abundance.sponges)
data2022$log.biomass <- log(data2022$Sponge.biomass)


plot(data2022$log.n, data2022$log.biomass, xlab = "log abundance", ylab = "log Biomass", main = "Relationship between Abundance and Biomass",type = "n")

j <- 0
for (i in unique(data2022$Habitat)) {
  points(data2022$log.n[data2022$Habitat == i], data2022$log.biomass[data2022$Habitat == i],
         pch = j + 1, col = c("pink", "blue", "purple","green", "cyan",'grey')[j + 1])
  j <- j + 1
}

legend("bottomright", legend = unique(data2022$Habitat),
       pch = 1:length(unique(data2022$Habitat)),
       col = c("pink", "blue", "purple", "green", "cyan", 'grey'),
       title = "Habitat", horiz = TRUE, cex = 0.35, pt.cex = 1.2, pt.bg = "white")

correlation3 <- cor(data2022$Abundance.sponges, data2022$Sponge.biomass)
correlation3


m10 <- lm(log.biomass ~ Habitat*log.n, data = data2022)
m10b <- lm(log.biomass ~ Habitat + log.n, data = data2022)
anova(m10b,m10)
summary(m10)

sdabundance <- sd(data2022$Abundance.sponges)
sdbiomass <- sd(data2022$Sponge.biomass)

sd2022 <- paste("Standard Deviation (Abundance):", sdabundance, "Standard Deviation (Biomass):", sdbiomass)
print(sd2022)

##ggplot1
ggplot(data2022, aes(x = log.n, y = log.biomass, color = Habitat)) +
  geom_point() +
  labs(x = "log abundance", y = "log Biomass", title = "Relationship between Abundance and Biomass") +
  theme_minimal() +
  theme(legend.position = "bottomright") +
  scale_color_manual(values = c("lightpink1", "lightgreen", "paleturquoise", "orchid1", "mediumpurple1", "palevioletred1"),
                     labels = unique(data2022$Habitat),
                     name = "Habitat")


library(ggplot2)

ggplot(data2022, aes(x = log.n, y = log.biomass, color = Habitat)) +
  geom_point() +
  labs(x = "log abundance", y = "log Biomass") +
  theme_minimal() +
  facet_wrap(~ Habitat, nrow = 2) +
  scale_color_manual(values = c("lightpink1", "lightgreen", "paleturquoise", "orchid1", "mediumpurple1", "palevioletred1"),
                     labels = unique(data2022$Habitat),
                     name = "Habitat")



*Model 2*
  
  
  ```{r, echo=FALSE}
residualsm2 <- residuals(m10b)
qqnormm2 <- qqnorm(residualsm2)
qqnormm2 <- qqline(residualsm2)

````
