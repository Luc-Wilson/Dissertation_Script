
setwd("c:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data")

Data_BrisGov = read.csv("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/Bristol_water_quality.csv", header = TRUE) 

if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if(!require("Rtools")) install.packages("Rtools")
if(!require("plotrix")) install.packages("plotrix")
library(plotrix) 
install.packages (c("dataRetrieval","EGRET"))
library(EGRET)
#for NetCDF
install.packages (c("raster","ncdf4"))
if(!require("gridExtra")) install.packages("gridExtra")
#================================================================

library(readr)
PiddleBagggsMill <- read_delim("c:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/PiddleBagggsMill.txt ", 
                               delim = "\t", escape_double = FALSE, 
                               col_types = cols(SMP_DateTime = col_date(format = "%Y-%m-%d"), 
                                                SMP_MeasurementResult = col_number()), 
                               trim_ws = TRUE)
view(PiddleBagggsMill)


#======================================================
#simple plot of concentration and time to show noise
library(tidyverse)

ConcentrationTime = 
  read.table("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/Exported CSVs/Solids, Suspended at 105 C.txt",
                  header = TRUE,
                  sep='\t') %>%   
  data.frame(.)


ConcentrationTime$SMP_DateTime = as.Date(ConcentrationTime$SMP_DateTime, format="%Y-%m-%d")

x_expression = expression(paste(Date))
y_expression = expression(paste(Solids, ~ Suspended ~ at ~ 105~degree~C))


plot(x=ConcentrationTime$SMP_DateTime, 
     y=ConcentrationTime$SMP_MeasurementResult,
#     xlim = NULL,
#    ylim = c(0, 0.4),    
     xlab = x_expression, # paste changes format of brackets so redo
     ylab = y_expression) +
  cor(x=ConcentrationTime$SMP_DateTime, 
      y=ConcentrationTime$SMP_MeasurementResult, 
      method = c("pearson"))


#=====================================================

# make SMP_determinands wide data - must be ran all together
subtable_a = PiddleBagggsMill[,c(4, 6, 7, 12)]
a = levels(as.factor(subtable_a$SMP_Determinand))[34] #number for which determinand
subtable_b = subtable_a[subtable_a$SMP_Determinand == a,]

#===============================================================
subtable_c = PiddleBagggsMill[,c(1, 2, 3, 4)]
c = levels(as.factor(subtable_c$SMP_Determinand))[110] #number for which determinand
subtable_d = subtable_c[subtable_c$SMP_Determinand == c,]

#export and save csv
#write.csv(subtable_d, 
# "C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/Exported CSVs\\Iron.csv",
#  row.names=FALSE)

 #==============================================
 
library(EGRET) 


 #reading in discharge data
 filePath <- "C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/"
 fileName <- "PiddleDischarge.txt"
 #specify discharge in cubic metres per second
 Daily <- readUserDaily(filePath, fileName, hasHeader = TRUE, separator = "\t",
                       qUnit = 2, verbose = TRUE)
 summary(Daily)
 #view(Daily)

 #=============================
 # Flow history analysis

 INFO <- readNWISInfo("", "")

 INFO = readNWISInfo(fileName)

#=============================
 
# Check flow history data:
 eList <- as.egret(INFO, Daily, NA, NA) #sometimes needs to be ran twice 
 #1st NA might be the name of the determinand

 library(gridExtra)
 # Creates a plot of a time series of a particular flow statistic and a loess smooth of that flow statistic
 plotFlowSingle(eList, istat=1,qUnit=2)  
 plotFlowSingle(eList, istat=8,qUnit=2)
 plotFlowSingle(eList, istat=2,qUnit=2)
 plotFlowSingle(eList, istat=7,qUnit=2)
 plotFlowSingle(eList, istat=4,qUnit=2)
 plotFlowSingle(eList, istat=5,qUnit=2) # would be nice to display all of these as a grid
 #istat: 1=minday,2=7daymin,3=30daymin,4=mediandaily,
 #5=meandaily,6=30daymax,7=7daymax,8=maxday
 
 plotSDLogQ(eList)
 plotFour(eList, qUnit=2)
 plotFourStats(eList, qUnit=2) # qUnit=2 specifies m^3/s
 
 
 # save work
 savePath<-"C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/" 
 saveResults(savePath, eList)

 #====================================
 #Water quality analysis 
 
 library(EGRET)
 
 #read in data
 filePathqual <- "C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/Exported CSVs/"
 fileNamequal <- "NitriteasN.txt"
 
 Sample <- readUserSample(filePathqual, fileNamequal, hasHeader = TRUE, separator = "\t",
                          verbose = TRUE)
 
 

INFO$shortName <- "River Piddle at Baggs Mill"

 eList <- mergeReport(INFO, Daily, Sample)

#====================================
 # Check sample data:
 boxConcMonth(eList)
 boxQTwice(eList)
 plotConcTime(eList) # not sure why the title says loess smooth, it is just a scatter plot of concentration and time
 plotConcQ(eList) # plots concentration against discharge
 multiPlotDataOverview(eList)
 
 #====================================
 # Run WRTDS model
 eList <- modelEstimation(eList, minNumUncen =50)
 #default minNumuUcensored is 50
 
 #====================================
 #Check model results:
 

 plotConcTimeDaily(eList, randomCensored = TRUE) 
 # random censored puts any less than values in a random point in their less than position 
 # this is so the plot is more readable and is not used in any computations
 plotFluxTimeDaily(eList, randomCensored = TRUE)
 plotConcPred(eList, randomCensored = TRUE, logScale = TRUE) # line is x=y
 plotFluxPred(eList, randomCensored = TRUE)
 plotResidPred(eList, randomCensored = TRUE)
 plotResidQ(eList, randomCensored = TRUE) 
 #stdResid = TRUE makes it standardized residuals instead of log of actual residuals. 
 plotResidTime(eList, randomCensored = TRUE)
 # NB hollow circles are randomized censored values - just for aesthetics
 boxResidMonth(eList) 
 # boxplot widths are proportional to the square root of the sample size
 boxConcThree(eList)
 
 
 plotConcHist(eList, col.pred = "black")
 plotFluxHist(eList, col.pred = "black")
 

 fluxBiasMulti(eList)
 

 #Contour plots:

 maxDiff<-0.8
 yearStart <- 1978
 yearEnd <- 2021
 
 
 clevel = seq(0,0.1,0.01)
 plotContours(eList, yearStart,yearEnd, 
              qBottom = NA,qTop = NA, 
              contourLevels = clevel, #NA means it is autoset, clevel for clevel
              qUnit=2,
              col = colorRampPalette(c("white","lightgray", "skyblue1","blue","orange", "firebrick1")), # not sure this is the best colour scheme
              whatSurface = 3, # 3 for concentration
              lwd = 0.5,  # line width default is 2
              printTitle = FALSE) 
 # different surfaces (1-3) are log concentration, standard error, and concentration


 plotDiffContours(eList, yearStart,yearEnd, # can specify new years with year0 and year1 - this might make the plot cover the two years.
                  qBottom=NA,qTop=NA,
                  maxDiff=NA,
                  qUnit=2,
                  span=60,
                  lwd = 0.5,
                  color.palette = colorRampPalette(c("blue", "skyblue1", "lightgray", "orange", "firebrick1"))) 
 #span is the number of days smoothing for the black lines
 # maxDiff=NA sets it to 5% and 95% bounds BUT
 #if plotPercent=TURE then it takes all results
 

 savePath<-"C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/" 
 setwd(savePath)
 save(file="eList.RSodium", eList) #or saveRDS
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RDataAmmoniacal Nitrogen as N")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RIron")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RLead")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RNitrate as N")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RNitriteasN") # this is nitrate. redo
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RNitrogen, Total Oxidised as N")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RNitrogenTotalN")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.ROrthophosphate, reactive as P")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RpH")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RPotassium")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RSilica, reactive as SiO2")
 #load("C:/Users/lewil/OneDrive/Documents/Geog BSC/Year 3/Dissertation/Data/EGRETS results/eList.RSodium")
 
 plotConcTime(eList)
 
 
 #===================
# Precipitation
library(tidyverse)
 library(dplyr)

prec = read.table("Precipitation_Data_catchment44.txt",
                  header = TRUE,
                  sep='\t') %>%
data.frame(.)


month_storage = matrix(data=NA, nrow=dim(prec)[1], ncol=1)

split_dates = str_split(prec$date, pattern="/")

for(d in 1:dim(prec)[1]){
  
on_date = split_dates[[d]][c(2,3)]  
on_date = paste0(on_date[2], "/", on_date[1])
month_storage[d,1] = on_date
  
}

prec$month = month_storage[,1]

aggregated_precipitation = aggregate(formula=precipitation ~ month, 
                                     data = prec, 
                                     FUN=mean)
view(aggregated_precipitation)


write.csv(aggregated_precipitation, "aggregated_precipitation.csv")
#=======================
# Discharge
Discharge = read.table("PiddleDischarge.txt",
                       header = TRUE,
                       sep='\t') %>%
  data.frame(.)


month_storage = matrix(data=NA, nrow=dim(Discharge)[1], ncol=1)

split_dates = str_split(Discharge$SMP_DateTime, pattern="-")


for(d in 1:dim(Discharge)[1]){
  
  on_date = split_dates[[d]][c(1,2)]  # was 2,3
  on_date = paste0(on_date[1], "/", on_date[2]) #was 2 and 1
  month_storage[d,1] = on_date
  
}

Discharge$month = month_storage[,1] 

aggregated_Discharge = aggregate(formula=Discharge ~ month, data = Discharge, FUN=mean)

view(aggregated_Discharge)

write.csv(aggregated_Discharge, "aggregated_Discharge.csv")

#=====================


precdis <- merge(x = aggregated_precipitation, 
                     y = aggregated_Discharge, 
                     by = "month")

precdis %>% separate(
  month, c("year", "month"), sep="/") %>%
  group_by("year") %>%
  summarise(mean = mean(precipitation)) %>%
  mutate(.after = "month")


#summarise() # was above


view(precdis)


lm_precdis <- lm(precipitation ~ Discharge, data = precdis)

summary(lm_precdis)
# adjusted r sqrd 0.1632
# p-value: < 2.2e-16


R2_expression = expression(paste("Adjusted" ~ R^2 ~ "= 0.1632"))

plot(precdis$precipitation,
     precdis$Discharge[85:624],
     ylab = "Discharge",
     xlab = "Precipitation")
abline(lm_precdis)
text(x = 8.5,
     y = 11.75,
     labels = R2_expression)


plot(precdis$precipitation)

class(prec$date)

precdis$month <- as.zoo(precdis$month) 


library(zoo)
precdis%>%
  as.yearmon(month)

plot(precdis$month, precdis$precipitation)

# monthly
precdis$month <- parse_date_time(precdis$month, "%Y/%m")
plot(precdis$precipitation ~ precdis$month, 
     type = "l",
     ylab = "Precipitation (mm)",
     xlab = "Time")
points(precdis$precipitation ~ precdis$month)

precdis%>%
  as.yearmon(month)


#========================
library(tidyverse)

#redone for annual
# precipitation
prec = readr::read_table("Precipitation_Data_catchment44.txt") %>% 
  # split the date into day month year
  separate(
    date, 
    sep='/', 
    into=c("day", "month", "year"), 
    remove=F
  )

# discharge 
discharge = readr::read_table("PiddleDischarge.txt") %>% 
  separate(
    SMP_DateTime,
    sep='-',
    into=c("year", "month", "day"),
  ) %>% 
  rename(discharge = Discharge)

# merge discharge and precipitation by day-month-year. 
data = left_join(prec, discharge, by=c("day", "month", "year")) %>%
  #  parse  "date" into datetime objects
  # lubridate to make day/month/year numeric
  mutate(date = lubridate::dmy(date),
         day = as.numeric(day), 
         month = as.numeric(month), 
         year = as.numeric(year)
  )

# make annual plots 
annuals = data %>% 
  group_by(year) %>% 
  summarise(
    annual_precip = mean(precipitation), 
    annual_discharge = sum(discharge)
  ) 


plot(annuals$year, annuals$annual_precip,
     type = "l",
     ylab = "Precipitation (mm)",
     xlab = "Time")

lm_annuals = lm(annual_precip ~ annual_discharge, data = annuals)
summary(lm_annuals) # comparing mean to total discharge
#Adjusted R-squared:  0.3683

plot(annuals$annual_precip ~ annuals$annual_discharge)
abline(lm_annuals)



# check the average
data %>% filter(year == 1978) %>% 
  # not for dis so can use ggplot
  ggplot(., aes(x=precipitation)) + 
  # make a density plot
  geom_density() + 
  # vertical line of annuals mean for 1978
  geom_vline(
    xintercept = annuals %>% filter(year == 1978) %>% pull(annual_precip),
    color='red'
  )
 
 