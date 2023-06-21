
library(car)
library(lme4)
library(multcomp)
library(ggplot2)
````
Hypothesis 1: The size-correct abundance of H.perlevis is influenced by month
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

prey_ab_vol <- read.csv('Dataa/prey_ab_vol.csv')
prey_ab_vol$n <- prey_ab_vol$N.H..perlevis
prey_ab_vol$row_id <- factor(1:nrow(prey_ab_vol))
prey_ab_vol$Date <- factor(prey_ab_vol$Date)
Summaryprey_ab_vol <- summary(prey_ab_vol)
summary(prey_ab_vol)
prey_ab_vol$ln.biovolume <- log(prey_ab_vol$Volumen.Hym..mL./prey_ab_vol$N.H..perlevis)
m1 <- glm(n ~ ln.biovolume*Date,family=poisson,data=prey_ab_vol)
overdisp_fun(m1) # reject null hypothesis of no overdispersion

summary(m1)

m2 <- glmer(n ~ ln.biovolume*Date +(1|row_id),family=poisson,data=prey_ab_vol)
anova(m2)
drop1(m2,test='Chisq') # homogeneity of slopes test, don't need interaction
summary(m2)

m3 <- glmer(N.H..perlevis ~ ln.biovolume + Date +(1|row_id),family=poisson,data=prey_ab_vol)
summary(m3)
drop1(m3,test='Chisq')
summary(glht(m3,linfct = mcp(Date="Tukey")))
overdisp_fun(m3)
summary(m3)

#plotting

ggplot(prey_ab_vol, aes(x = Date, y = sqrt(N.H..perlevis))) +
  geom_boxplot(fill = "#FFDAB9", color = "black") +
  labs(title = "Relationship between abundance (n) per quadrant and month") +
  xlab("Date") +
  ylab("sqrt(abundance)")


pred <- predict(m3,type='response',re.form=NA)
ggplot(prey_ab_vol,aes(x=ln.biovolume,y=log(N.H..perlevis),color=Date)) + geom_point() +
  geom_line(aes(y=log(pred)))


#assumptions
residuals <- resid(m3, type = "pearson")
qqnorm(residuals)
qqline(residuals)

qqnorm(ranef(m1)$row_id[,1]); qqline(ranef(m1)$row_id[,1])


Hypothesis 2: Does the month have a effect on the volume of H.perlevis?
  

prey_ab_vol$log.v <-log(prey_ab_vol$Volumen.Hym..mL.)
prey_ab_vol <- prey_ab_vol[, colSums(is.na(prey_ab_vol)) == 0]

m7a <- aov(log.v~Date,data=prey_ab_vol)
anova(m7a)
TukeyHSD(m7a,ordered=TRUE)


library(ggplot2)

ggplot(prey_ab_vol, aes(x = Date, y = sqrt(N.H..perlevis))) +
  geom_boxplot(fill = "thistle", color = "black") +
  labs(title = "Relationship between total volume of H.perlevis (mL) and Date ") +
  xlab("Date") +
  ylab("sqrt(abundance)")

  
residualsm5 <- residuals(m7)
qqnormm5 <- qqnorm(residualsm5)
qqnormm5 <- qqline(residualsm5)

leveneTest(m7a)

```


Hypothesis 3: Is there a relationship between average volume of H.perlevis per month and pH, salinity and temperature averages for the months**</div></font>
  

TpHSal <- read.csv('Dataa/TempsalinitypH.csv')
empty_columns <- which(colSums(is.na(TpHSal)) == nrow(TpHSal))
TpHSal <- subset(TpHSal, select = -empty_columns)
summaryTpHSal <- summary(TpHSal)
print(summaryTpHSal)
sdT <- sd(TpHSal$T)
sdpH <- sd(TpHSal$PH)
Sal <- sd(TpHSal$S)
SDresults2 <- paste("Standard Deviation (T):", sdT, "Standard Deviation (pH):", sdpH, "Standard Deviation (Sal):", Sal)
print(SDresults2)
library(dplyr)
data <- data.frame(
  ID_SITE = c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
  DATE = as.Date(c("10/6/2020", "18/6/2020", "26/6/2020", "5/11/2020", "12/11/2020", "20/11/2020", "27/11/2020", "3/4/2021", "8/4/2021", "16/4/2021", "7/6/2021", "15/6/2021", "21/6/2021", "29/6/2021"), format = "%d/%m/%Y"),
  T = c(15.416, 15.892, 16.824, 13.417, 13.241, 13.048, 12.848, 11.072, 11.23, 10.971, 15.061, 16.295, 15.688, 15.399),
  S = c(34.733, 34.632, 34.704, 34.552, 34.424, 34.56, 34.31, 33.62, 33.8, 34.156, 34.626, 34.762, 34.865, 35.014),
  PH = c(8.138, 7.959, 7.953, 7.799, 8.055, 7.791, 7.794, 7.9, 7.93, 7.91, 7.94, 7.91, 7.93, 7.88)
)
data$DATE <- as.POSIXlt(data$DATE)
data$Month <- data$DATE$mon + 1
data$Year <- data$DATE$year + 1900
averagestemp_ph_S <- aggregate(cbind(T, S, PH) ~ Month + Year, data, FUN = mean)
june_2020 <- subset(averagestemp_ph_S, Month == 6 & Year == 2020)
june_2021 <- subset(averagestemp_ph_S, Month == 6 & Year == 2021)
print(averagestemp_ph_S)
averages2 <- aggregate(Volumen.Hym..mL. / N.H..perlevis ~ Date, prey_ab_vol, FUN = mean)
print(averages2)
averagesvolume <- data.frame(
  Month = c(6, 11, 4, 6),
  Year = c(2020, 2020, 2021, 2021),
  Volumen.Hym..mL..N.H..perlevis = c(3.8781924, 12.6688840, 3.0448443, 0.8122513)
)

merged_dataset <- merge(averagesvolume, averagestemp_ph_S, by = c('Year', 'Month'))
correlation2 <- cor(merged_dataset[, c("Volumen.Hym..mL..N.H..perlevis", "T", "S", "PH")])
print(correlation2)

  
Hypothesis 4: What is the relationship between the abundance of H.perlevis and biomass per haibitat (2022 study)**</div></font>
 
library(ggplot2)
data2022 <- read.csv('Dataa/data2022.csv')
data2022$row_id <- factor(1:nrow(data2022))
data2022 <- data2022[data2022$Abundance.sponges > 0, ]
data2022$Sponge.biomass <- as.numeric(gsub("[^0-9.]", "", data2022$Sponge.biomass))
data2022$log.biomass <- log(data2022$Sponge.biomass)
summarydata2022 <- summary(data2022)
print(summarydata2022)
sdabundance <- sd(data2022$Abundance.sponges)
sdbiomass <- sd(data2022$Sponge.biomass)
sd2022 <- paste("Standard Deviation (Abundance):", sdabundance, "Standard Deviation (Biomass):", sdbiomass)
print(sd2022)
data2022$ln.biovolume <- log(data2022$Sponge.biomass/data2022$Abundance.sponges)
data2022$n <- data2022$Abundance.sponges
data2022$log.n <- log(data2022$n)
data2022$row_id <- factor(1:nrow(data2022))
m10 <- glm(n ~ Habitat*ln.biovolume, family=poisson,data = data2022)
overdisp_fun(m10)
m10a <- glm.nb(n ~ Habitat*ln.biovolume,data = data2022)
drop1(m10a,test='Chisq') #homogeneity of slopes, no need for interaction term
m10b <- glm.nb(n ~ Habitat + ln.biovolume,data = data2022)
overdisp_fun(m10b)

pred <- predict(m10b)
ggplot(data2022,aes(x=ln.biovolume,y=log(n),color=Habitat)) + geom_point() +
  geom_line(aes(y=pred))

plot(fitted.values(m10b),residuals(m10b,type='deviance'))
lines(lowess(fitted.values(m10b),residuals(m10b,type='deviance')))

qqnorm(residuals(m10b,type='deviance'))
qqline(residuals(m10b,type='deviance'))


## Hypothesis 5: biomass

m10b <- aov(log.biomass ~ Habitat, data = data2022)
summary(m10b) 
anova(m10b)
drop1(m10b,test='Chisq')
TukeyHSD(m10b)

dog <- lm(log.biomass~log.n + Habitat,data2022)
pred <- predict(dog)
ggplot(data2022,aes(x=log.n,y=log.biomass,color=Habitat)) + geom_point() +
  geom_line(aes(y=pred))

#assumptions

residualsm3 <- residuals(m10b)
qqnormm3 <- qqnorm(residualsm3)
qqnormm3 <- qqline(residualsm3)



