#Code written 2014-015 by Jeff Atkins (jeffatkins@virginia.edu or via Twitter @atkinsjeff)
#This code analyzes data from the MesoWest Weather Station at Bearden Knob (BDKW2) located at the top of the Weimer Run Watershed in David, WV
#This code also analyzes long term trends in precipitation using NCDC data for the Weimer Run Watershed
#Dataset files are uploaded into figshare
#and is in suppor of research outlined in the manuscript Atkins et al. 2014 Biogeosciences (article title final here)

# Station data downloaded from http://mesowest.utah.edu/ for the Bearden Knob Weather Station (BDKW2)
# 
# http://mesowest.utah.edu/html/help/main_index.html#usage
# 
# Column descriptions and units:
# DATETIME - date (mm/dd/yyyy) with time (hh:mm)
# TEMP - air temperature in degrees Celsius
# RELH - Relative humidity (%)
# WIND - hourly wind speed (m/s)
# GUST - peak wind speed gust during hour (m/s)
# DIR - wind speed direction (compass degrees)
# QFLG - Quality of data (OK = good to use; CAUTION = examine for incongruity)
# SOLAR - Incoming solar radiation (w/m^2)
# PRECIPcum - cumulative precipitation (cm)
# PRECIP - precip during past hour (cm)
# PEAK - duplicate column of max wind gust (m/s)
# PEAKDIR - direction of PEAK and GUST columns (compass degrees)
# DWP - dew point (degrees Celsius)
# 


#required packages
require(xts)
require(ggplot2)
require(psych)
library(Hmisc)
library(lmtest)
library(car)
library(gtable)
library(gridExtra)

#Importing MesoWest data that has been imported to figshare
MAIN2010 <-read.csv("http://files.figshare.com/1893326/BDKW_2010.csv",head=TRUE, strip.white=TRUE )
MAIN2011 <-read.csv("http://files.figshare.com/1893327/BDKW_2011.csv",head=TRUE, strip.white=TRUE )
MAIN2012 <-read.csv("http://files.figshare.com/1893328/BDKW_2012.csv",head=TRUE, strip.white=TRUE )

MAIN2010$PRECIP <-as.numeric(as.character(MAIN2010$PRECIP))
MAIN2011$PRECIP <-as.numeric(as.character(MAIN2011$PRECIP))
MAIN2012$PRECIP <-as.numeric(as.character(MAIN2012$PRECIP))
#converting DATETIME to POSIXct class 
MAIN2010$DATETIME <-  as.POSIXct(strptime(as.character(MAIN2010$DATETIME), format = "%m-%d-%Y %H:%M"))
MAIN2011$DATETIME <-  as.POSIXct(strptime(as.character(MAIN2011$DATETIME), format = "%m-%d-%Y %H:%M"))
MAIN2012$DATETIME <-  as.POSIXct(strptime(as.character(MAIN2012$DATETIME), format = "%m-%d-%Y %H:%M"))

#aggregating precip values to the day from hourly data
dp2010 <- aggregate(list(precip = MAIN2010$PRECIP),
                          list(day = cut(MAIN2010$DATETIME,"day")),
                          sum)
dp2011 <- aggregate(list(precip = MAIN2011$PRECIP),
                          list(day = cut(MAIN2011$DATETIME,"day")),
                          sum)
dp2012 <- aggregate(list(precip = MAIN2012$PRECIP),
                          list(day = cut(MAIN2012$DATETIME,"day")),
                          sum)
  
dp2010 <-data.frame(dp2010)
dp2011 <-data.frame(dp2011)
dp2012 <-data.frame(dp2012)


###creating a temp row to account for JD mismatch from leap year
temprow <- matrix(c(rep.int(NA,length(dp2010))),nrow=1,ncol=length(dp2010))

newrow <- data.frame(temprow)
colnames(newrow) <- colnames(dp2010)

# rbind the empty row to data
dp2010 <- rbind(dp2010,newrow)
dp2011 <-rbind(dp2011,newrow)


#Adding in Julian Day to make life easier and because I like using Julian Day
dp2010$JD <- c(1:366)
dp2011$JD <- c(1:366)
dp2012$JD <- c(1:366)

#adding a year column
dp2010$YEAR <- c("2010")
dp2011$YEAR <- c("2011")
dp2012$YEAR <- c("2012")

#combinging all the data back to one df
DAILYx <- rbind(dp2010,dp2011,dp2012)

# convert precip to mm
DAILYx$precip <-DAILYx$precip*10

#adding the month
attach(DAILYx)
DAILYx$MONTH[JD<=31] <-"JAN"
DAILYx$MONTH[JD>=32 & JD < 60] <-"FEB"
DAILYx$MONTH[JD>=60 & JD<= 90] <-"MAR"
DAILYx$MONTH[JD>=91 & JD<= 120] <-"APR"
DAILYx$MONTH[JD>=121 & JD<= 151] <-"MAY"
DAILYx$MONTH[JD>=152 & JD<= 181] <-"JUN"
DAILYx$MONTH[JD>=182 & JD<= 212] <-"JUL"
DAILYx$MONTH[JD>=213 & JD<= 243] <-"AUG"
DAILYx$MONTH[JD>=244 & JD<= 273] <-"SEP"
DAILYx$MONTH[JD>=274 & JD<= 304] <-"OCT"
DAILYx$MONTH[JD>=305 & JD<= 332] <-"NOV"
DAILYx$MONTH[JD>=333 & JD<= 366] <-"DEC"

#removes the NA values created from the lack of leap year in 2010 and 2011
DAILYx <- na.omit(DAILYx)

#labels for each year with precip totals

labels <-data.frame(YEAR = c("2010","2011","2012"), lab = c("2010(1042 mm)", "2011(1739 mm)", "2012(1243 mm)"))

#Facet grid plots of daily precip (mm) by month

#black and white 
pBW <-ggplot(DAILYx, aes(x=JD, y=precip, fill=YEAR))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#999999", "#666666","#333333"))+
  theme_classic()+
  guides(fill=FALSE)+
#   theme(legend.justification=c(1,1), 
#         legend.position=c(1,1),
#         legend.title=element_blank(),
#         legend.text=element_text(color='black', size=20), 
#         legend.background = element_rect( fill="#FFFFFF"), 
#         legend.key=element_rect(color="#FFFFFF", fill="#FFFFFF") )+
  labs(y="P (mm) ")+
  labs(x="Julian Day")+
  theme(axis.title.x=element_text(size=20))+ 
  theme(axis.title.y= element_text(size=20),
        axis.text.x = element_text(size=(18)),
        axis.text.y = element_text(size=(18)))+
  geom_text(x=340,y=50, aes(label=lab), data=labels, size=8)+
  facet_grid(YEAR ~ .)+
  theme(strip.text.y =  element_blank(),
        strip.background = element_blank())

pBW



p.color <-ggplot(DAILYx, aes(x=JD, y=precip, fill=YEAR))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#FF0000", "#00A08A", "#F2AD00"))+
  theme_classic()+
  guides(fill=FALSE)+
  labs(y="P (mm) ")+
  labs(x="Julian Day")+
  theme(axis.title.x=element_text(size=20))+ 
  theme(axis.title.y= element_text(size=20),
        axis.text.x = element_text(size=(18)),
        axis.text.y = element_text(size=(18)))+
  geom_text(x=340,y=50, aes(label=lab), data=labels, size=8)+
  facet_grid(YEAR ~ .)+
  theme(strip.text.y =  element_blank(),
        strip.background = element_blank())

p.color

# ### Weimer Long term precipitation 
# NCDC data from Davis, WV 26260
# Data submitted/requested September 2, 2014
# Station COOP:461393 CANAAN VALLEY WV US
# 
# DSNW - Number days with snow depth > 1 inch.
# HTDD - Heating degree days
# DP10 - Number of days with greater than or equal to 1.0 inch of precipitation
# MXSD - Maximum snow depth (mm)
# EMXT - Extreme maximum daily temperature
# DPNP - Departure from normal monthly precipitation.
# DP01 - Number of days with greater than or equal to 0.1 inch of precipitation
# MMNT - Monthly Mean minimum temperature
# TPCP - Total precipitation (tenths of mm)
# EMXP - Extreme maximum daily precipitation (tenths of mm)
# DT00 - Number days with minimum temperature less than or equal to 0.0 F
# MMXT - Monthly Mean maximum temperature
# DT90 - Number days with maximum temperature greater than or equal 90.0 F
# DT32 - Number days with minimum temperature less than or equal to 32.0 F
# CLDD - Cooling degree days
# EMNT - Extreme minimum daily temperature
# TSNW - Total snow fall (mm)
# MNTM - Monthly mean temperature
# DP05 - Number of days with greater than or equal to 0.5 inch of precipitation
# DX32 - Number days with maximum temperature less than or equal to 32.0 F
# DPNT - Departure from normal monthly temperature.

#importing NCDC data from figshare
NCDC <-read.csv("http://files.figshare.com/1894690/ncdc_davis_1948_2014.csv", head=TRUE, strip.white=TRUE, na.strings= -9999 )

#subsetting data after 1970
NCDC70 <-subset(NCDC, NCDC$DATE >197000)

#adding year row and formatting
NCDC70$YEAR <- substr(NCDC70$DATE, 1, 4)
NCDC70$YEAR <- as.POSIXlt(NCDC70$YEAR,"%Y")
NCDC70$YEAR <- as.numeric(format(NCDC70$YEAR,"%Y"))

#structure of df
str(NCDC70)

#yearly totals
NCDCsums <- aggregate(formula= TPCP~YEAR, data=NCDC70, FUN=sum)

#yearly totals of extreme precip (above 1 inch)
NCDCex <-aggregate(formula= DP10~YEAR, data=NCDC70, FUN=sum )

#convert to mm
NCDCsums$TPCP <- (NCDCsums$TPCP / 10)
NCDCsums <- NCDCsums[-45,]
NCDCex <- NCDCex[-45,]


#scatter plot of annual precipitation totals (mm) since 1970
pNCDC <- ggplot(NCDCsums, aes(x=as.numeric(YEAR), y=TPCP))+
  geom_point(size=10, alpha=0.8, color="#F98400")+
  theme_classic()+
  theme(legend.justification=c(1,1), 
        legend.position=c(1,1),
        legend.title=element_blank(),
        legend.text=element_text(color='black', size=24), 
        legend.background = element_rect( fill="#FFFFFF"), 
        legend.key=element_rect(color="#FFFFFF", fill="#FFFFFF") )+
  scale_y_continuous(limits=c(0,1500), expand = c(0,0))+
  labs(y="P (mm) ")+
  labs(x="")+
  theme(axis.title.x=element_text(size=28))+ 
  theme(axis.title.y= element_text(size=28),
        axis.text.x = element_text(size=(20)),
        axis.text.y = element_text(size=(20)))+

  stat_smooth(method=lm, se=FALSE, size=1.5, color="black")
  
pNCDC

#linear regression modeland correlation coeff for increasing precip since 1970
fit.1970 <- lm(NCDCsums$TPCP~as.numeric(NCDCsums$YEAR))
summary.lm(fit.1970)
cor(NCDCsums$TPCP, as.numeric(NCDCsums$YEAR))

#testing for heteroskedasticity since 1970 using a Breusch-Pagan test
bptest(fit.1970)
coeftest(fit.1970, vcov=hccm(fit.1970))

#looking at extreme precip days
#scatter plot of extreme precipitation days (where precip was >1 in or 25.4 mm) since 1970
pex <- ggplot(NCDCex, aes(x=as.numeric(YEAR), y=DP10))+
  geom_point(size=10, alpha=0.8, color="#5BBCD6")+
  theme_classic()+
  theme(legend.justification=c(1,1), 
        legend.position=c(1,1),
        legend.title=element_blank(),
        legend.text=element_text(color='black', size=24), 
        legend.background = element_rect( fill="#FFFFFF"), 
        legend.key=element_rect(color="#FFFFFF", fill="#FFFFFF") )+
  #scale_x_continuous(limits=c(0,25))+
  scale_y_continuous(limits=c(0,40), expand = c(0,0))+
  labs(y="No. of EPD")+
  labs(x="")+
  theme(axis.title.x=element_text(size=28))+ 
  theme(axis.title.y= element_text(size=28),
        axis.text.x = element_text(size=(20)),
        axis.text.y = element_text(size=(20)))+
  
  stat_smooth(method=lm, se=FALSE, size=1.5, color="black")

pex 

#linear regression of extreme precip days since 1970
fit.ex <- lm(NCDCex$DP10~as.numeric(NCDCex$YEAR))
cor(NCDCex$DP10, as.numeric(NCDCex$YEAR))
summary.lm(fit.ex)

#testing for heteroskedasticity since 1970 using a Breusch-Pagan test
bptest(fit.ex)
coeftest(fit.ex, vcov=hccm(fit.ex))



####plotting all graphs for final, combined graph in manuscript
# Get the widths
gA <- ggplot_gtable(ggplot_build(p.color))
gB <- ggplot_gtable(ggplot_build(pNCDC))
gC <- ggplot_gtable(ggplot_build(pex))


maxWidth = unit.pmax(gA$widths[1:3],
                     gB$widths[2:3], 
                     gC$widths[2:3])
# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth


# Arrange the four charts
grid.arrange(gA, gB, gC, nrow=3)
