#Code written 2012-2015 by Jeff Atkins (jeffatkins@virginia.edu or via Twitter @atkinsjeff)
#This code analyzes data collected in the Weimer Run watershed from 2010-2012 and supports
#the research artice Atkins et al. 2014 Biogeosciences (article title final here)

# This code supports includes the bulk of the analysis conducte in the aforementioned manuscript and consists of several parts
# including a script to convert volumetric water content measured in the field to water-filled pore space values based
# on empirically derived coefficients, coefficients from the manufacturer and from field-measured soil porosity values.

# This script alsou outputs the plots included in the manuscript including bar graphs of flux, WFPS, and soil temperature (Figure 3), 
# all flux/residual plots (Figure 2), and also creates soil/flux regression plots included in supplemental materials
# There is also a temperature threshold section that is discussed in then interactive discussion

# The last section of the script focuses on soil analyses found in the manuscript (Figure 4, Table A3)


#required packages
require(plyr)
require(ggplot2)
require(segmented)
require(gridExtra)
require(gtable)
require(polynom)
#require(wesanderson)  #Colors used later come from the Darjeeling and Darjeeling2 palettes in this package, but version conflicts didn't work, but 
# you should check out the package anyway though (https://github.com/karthik/wesanderson)

#importing soils and flux data from figshare (http://figshare.com/authors/Jeff_Atkins/665669)
soil.data <- read.csv("http://files.figshare.com/1806287/corrected_master_2014_soil_data.csv", header=TRUE)
c.n.data <- read.csv("http://files.figshare.com/1806288/soil_c_and_n_data_weimer_run.csv", header=TRUE)
CVI <- read.csv("http://files.figshare.com/1806266/weimer_run_flux_and_micromet_data_2010_2012.csv", header=TRUE)

#WFPS script
#This creates the ID column and the WFPS (water-filled pore space) values
################

#making the ID column
attach(CVI)
CVI$ID[ELEV=="low" & VEG=="open" & PLOT == 1] <-'LO1'
CVI$ID[ELEV=="low" & VEG=="open" & PLOT == 2] <-"LO2"
CVI$ID[ELEV=="low" & VEG=="open" & PLOT == 3] <-"LO3"
CVI$ID[ELEV=="low" & VEG=="shrub" & PLOT == 1] <-"LS1"
CVI$ID[ELEV=="low" & VEG=="shrub" & PLOT == 2] <-"LS2"
CVI$ID[ELEV=="low" & VEG=="shrub" & PLOT == 3] <-"LS3"
CVI$ID[ELEV=="low" & VEG=="canopy" & PLOT == 1] <-"LC1"
CVI$ID[ELEV=="low" & VEG=="canopy" & PLOT == 2] <-"LC2"
CVI$ID[ELEV=="low" & VEG=="canopy" & PLOT == 3] <-"LC3"

CVI$ID[ELEV=="mid" & VEG=="open" & PLOT == 1] <-"MO1"
CVI$ID[ELEV=="mid" & VEG=="open" & PLOT == 2] <-"MO2"
CVI$ID[ELEV=="mid" & VEG=="open" & PLOT == 3] <-"MO3"
CVI$ID[ELEV=="mid" & VEG=="shrub" & PLOT == 1] <-"MS1"
CVI$ID[ELEV=="mid" & VEG=="shrub" & PLOT == 2] <-"MS2"
CVI$ID[ELEV=="mid" & VEG=="shrub" & PLOT == 3] <-"MS3"
CVI$ID[ELEV=="mid" & VEG=="canopy" & PLOT == 1] <-"MC1"
CVI$ID[ELEV=="mid" & VEG=="canopy" & PLOT == 2] <-"MC2"
CVI$ID[ELEV=="mid" & VEG=="canopy" & PLOT == 3] <-"MC3"

CVI$ID[ELEV=="high" & VEG=="open" & PLOT == 1] <-"HO1"
CVI$ID[ELEV=="high" & VEG=="open" & PLOT == 2] <-"HO2"
CVI$ID[ELEV=="high" & VEG=="open" & PLOT == 3] <-"HO3"
CVI$ID[ELEV=="high" & VEG=="shrub" & PLOT == 1] <-"HS1"
CVI$ID[ELEV=="high" & VEG=="shrub" & PLOT == 2] <-"HS2"
CVI$ID[ELEV=="high" & VEG=="shrub" & PLOT == 3] <-"HS3"
CVI$ID[ELEV=="high" & VEG=="canopy" & PLOT == 1] <-"HC1"
CVI$ID[ELEV=="high" & VEG=="canopy" & PLOT == 2] <-"HC2"
CVI$ID[ELEV=="high" & VEG=="canopy" & PLOT == 3] <-"HC3"

#attaching measured soil porosity values (m^3 m^-3)
attach(CVI)
CVI$PORE[ID=="LO1"] <-0.67886641
CVI$PORE[ID=="LS1"] <-0.773572016
CVI$PORE[ID=="LC1"] <-0.760341411
CVI$PORE[ID=="LO2"] <-0.593540032
CVI$PORE[ID=="LS2"] <-0.739923893
CVI$PORE[ID=="LC2"] <-0.654671743
CVI$PORE[ID=="LO3"] <-NA
CVI$PORE[ID=="LS3"] <-0.95121
CVI$PORE[ID=="LC3"] <-0.893877
CVI$PORE[ID=="MO1"] <-0.644723
CVI$PORE[ID=="MS1"] <-0.77067
CVI$PORE[ID=="MC1"] <-0.6703927
CVI$PORE[ID=="MO2"] <-0.710620
CVI$PORE[ID=="MS2"] <-0.614471
CVI$PORE[ID=="MC2"] <-0.770767
CVI$PORE[ID=="MO3"] <-0.74201
CVI$PORE[ID=="MS3"] <-0.908828
CVI$PORE[ID=="MC3"] <-0.577685
CVI$PORE[ID=="HO1"] <-0.71795
CVI$PORE[ID=="HS1"] <-0.909273
CVI$PORE[ID=="HC1"] <-0.681454
CVI$PORE[ID=="HO2"] <-0.70935
CVI$PORE[ID=="HS2"] <-0.922143
CVI$PORE[ID=="HC2"] <-0.65663
CVI$PORE[ID=="HO3"] <-0.557315
CVI$PORE[ID=="HS3"] <-0.886787
CVI$PORE[ID=="HC3"] <-0.764728



# this is important as it eliminates all of the LO3 plots as explained in manuscript
CVI <-na.omit(CVI)

#formatting year and adding ln(flux)
CVI$YEAR <- as.character(CVI$YEAR)
CVI$LN_FLUX <- log(CVI$FLUX)

##############################################
#####   TEMPERATURE THRESHOLD SECTION    ##### 
##############################################

# Piecewise regression using package(segmented) to establish breakpoints of slop in temperature threshold

# plotting temperature against the natural log of flux
plot(CVI$TEMP, CVI$LN_FLUX)

# This gets rid of the infinite values in the data set. I think there are only two.
CVI.ln <- subset(CVI, LN_FLUX!=-Inf)

# Creating linear model
cvi.lin.mod <- lm(LN_FLUX~TEMP, na.rm =TRUE, data=CVI.ln)
segmented.cvi.lin.mod <- segmented(cvi.lin.mod, seg.Z= ~TEMP, psi=10)

# summary of piecewise regression
summary(segmented.cvi.lin.mod)

#creating plots
plot(segmented.cvi.lin.mod, conf.level=0.95, shade=TRUE)




####plot average tempertaure to determine when the average measured soil temperature exceeds 11 degree C

# julian day vs. soil temp
plot(CVI$JD, CVI$TEMP)
abline(h=11)

# Using a fourth order polynomial to evaluate relationship
fit.jd2 <- lm(CVI$TEMP ~ poly(CVI$JD, 4, raw=TRUE))
summary(fit.jd2)

# adding polynomial function to plot
plot(CVI$JD, CVI$TEMP)
abline(h=11)
points(CVI$JD, predict(fit.jd2), col="red", lwd=2)


fit.jd2.eq <- polynomial(coef(fit.jd2))

# Solving polynomial for julian day
predict(fit.jd2.eq, newdata=127)
predict(fit.jd2.eq, newdata=287)

##############################################
#####          WFPS SCRIPT               ##### 
##############################################

##### Making appropriate calibration of THETA based on lab results (outlined in manuscript)
CVI$THETA = CVI$THETA /100
CVI$PERIOD <- (-0.3385* (CVI$THETA^2))+(0.7971*CVI$THETA)+0.7702
CVI$Ka <- ((CVI$PERIOD-0.79)/ ((1.37-0.79) * (sqrt(80)-1) +1 ))
CVI$THETA_c <- (7.0341*CVI$Ka)+0.0806         #this is the recalculated THETA value based on soil type

#creates the WFPS column (m^3 m^-3)
CVI$WFPS <- CVI$THETA_c * CVI$PORE

#checking data to see if the structure looks correct
head(CVI)

#Now we restrict data to soil temperature above 11 degrees C
CVI_ABOVE <- subset(CVI, TEMP >=11)
CVI_BELOW <- subset(CVI, TEMP< 11)



# renaming column PLOT in soils data to match with the ID column in the CVI data structure as there is a mismatch
# and the columns need to match
soil.data <- rename(soil.data, c("plot" = "ID"))
c.n.data <- rename(c.n.data, c("PLOT" = "ID"))

#checking it 
head(soil.data)
head(c.n.data)


#Then we want to subset the 0 - 12 soil depth profile for analysis 
#This profile is used to be in phase with THETA and TEMP measurements as outlined in the manuscript
soil.12 <- subset(soil.data, depth == 12)
c.n.12 <- subset(c.n.data, DEPTH == "TWELVE")

head(soil.12)
head(c.n.12)

#making means for soil CO2 flux, soil builk density, and soil organic matter (%) using the plyr package

flux <- ddply(CVI_ABOVE, c("ID"), summarise, 
               flux.mean = mean(FLUX))

bulk <- ddply(soil.12, c("ID"), summarise,
               bulk.mean = mean(soil.bulk.density))

som <- ddply(soil.12, c("ID"), summarise,
                      som.mean = mean(som))

carbon <- ddply(c.n.12, c("ID"), summarise,
               carbon.mean = mean(C_per))

nitrogen <- ddply(c.n.12, c("ID"), summarise,
                nitrogen.mean = mean(N_per))




#combining means into one data frame. Admittedly there is a better way to do this I should look into
df.means <- merge(flux, bulk, "ID")
df.means <- merge(df.means, som, "ID")
df.means <- merge(df.means, carbon, "ID")
df.means <- merge(df.means, nitrogen, "ID")

# manually adding PAI values

df.means$pai.mean <- c(2.53, 3.58,  3.1,	2.14,	2.88,	4.48,	3.72,	5.93,	1.27,	2.39,	1.81,	0.59,	1.54,	2.21,	2.3, 1.31,	2.31,	1.01,	0.94,	1.87,	1.67,	4.93,	2.63,	3.48) 

#######################################################
#################    graph section    #################
#######################################################

#labels for graphs
fluxlabel = expression(paste(F[SOIL]~(µmol ~CO[2] ~m^-2 ~s^-1)), parse=TRUE)
lnfluxlabel = expression(paste(ln(F[SOIL]~(µmol ~CO[2] ~m^-2 ~s^-1))), parse=TRUE)
residfluxlabel = expression(paste(Residuals~(µmol ~CO[2] ~m^-2 ~s^-1)), parse=TRUE)
templabel = expression(paste(T[SOIL] ~(degree~C)), parse=TRUE )
wfpslabel = expression(paste(WFPS~(m^3~m^-3)))
bulklabel = expression(paste(rho[S]~(g ~cm^-3)), parse=TRUE)
somlabel = expression(paste(SOM~(cm^3 ~cm^-3 )), parse=TRUE)
pailabel = expression(paste(PAI~(m^3 ~m^-3 )), parse=TRUE)
sbulklabel = expression(paste(rho[s]~( ~g ~cm^-3 )), parse=TRUE)
tbulklabel = expression(paste(rho[t]~( ~g ~cm^-3 )), parse=TRUE)

######### SOIL/FLUX REGRESSION PLOTS

p.bulk <- ggplot(df.means, aes(x=bulk.mean, y=flux.mean)) +
  geom_point(size=10, alpha=0.75, colour="#FF0000") +
  scale_y_continuous(limits=c(0,10))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='white', size=24), legend.background = element_rect( fill="#404040"), 
        legend.key=element_rect(color="#404040", fill="#404040") )+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  ylab(fluxlabel)+
  xlab(bulklabel)+
  geom_smooth(method=lm, se=FALSE, color="black", size=1)

p.som <- ggplot(df.means, aes(x=som.mean, y=flux.mean)) +
  geom_point(size=10, alpha=0.7, colour="#F2AD00") +
  scale_y_continuous(limits=c(0,10))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='white', size=24), legend.background = element_rect( fill="#404040"), 
        legend.key=element_rect(color="#404040", fill="#404040") )+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  ylab(fluxlabel)+
  xlab(somlabel)

p.carbon <- ggplot(df.means, aes(x=carbon.mean, y=flux.mean)) +
  geom_point(size=10, alpha=0.7, colour="#046C9A") +
  scale_y_continuous(limits=c(0,10))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='white', size=24), legend.background = element_rect( fill="#404040"), 
        legend.key=element_rect(color="#404040", fill="#404040") )+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  ylab(fluxlabel)+
  xlab("Total Soil C (%)")

p.nitrogen<- ggplot(df.means, aes(x=nitrogen.mean, y=flux.mean)) +
  geom_point(size=10, alpha=0.7, colour= "#D69C4E") +
   scale_y_continuous(limits=c(0,10))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='white', size=24), legend.background = element_rect( fill="#404040"), 
        legend.key=element_rect(color="#404040", fill="#404040") )+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  ylab(fluxlabel)+
  xlab("Total Soil N (%)")

p.pai<- ggplot(df.means, aes(x=pai.mean, y=flux.mean)) +
  geom_point(size=10, alpha=0.7, colour="#00A08A") +
  scale_y_continuous(limits=c(0,10))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='white', size=24), legend.background = element_rect( fill="#404040"), 
        legend.key=element_rect(color="#404040", fill="#404040") )+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  ylab(fluxlabel)+
  xlab(pailabel)

grid.arrange(p.bulk, p.som, p.carbon, p.nitrogen, p.pai, nrow=5)


################################3
# Soil/flux regression models and summary


bulk.lm <- lm(flux.mean ~ bulk.mean, data=df.means)
som.lm <- lm(flux.mean ~ som.mean, data=df.means)
carbon.lm <- lm(flux.mean ~ carbon.mean, data=df.means)
nitrogen.lm <- lm(flux.mean ~ nitrogen.mean, data=df.means)
pai.lm <- lm(flux.mean ~ pai.mean, data=df.means)

#making a list because the functions I kept writing didn't work. I need to get better at that
# then arranging data into a table using the adjusted R-squared from the linear models

lm.list = list() 

lm.list[[1]] = bulk.lm
lm.list[[2]] = som.lm
lm.list[[3]] = carbon.lm
lm.list[[4]] = nitrogen.lm
lm.list[[5]] = pai.lm

lm.table <- matrix( c(summary(lm.list[[1]])$adj.r.squared,
                    summary(lm.list[[2]])$adj.r.squared,
                    summary(lm.list[[3]])$adj.r.squared,
                    summary(lm.list[[4]])$adj.r.squared,
                    summary(lm.list[[5]])$adj.r.squared                    
                    ), ncol=1, byrow=TRUE)

colnames(lm.table) <- c( "Model R-Sqaured")
rownames(lm.table) <- c("Bulk Density", "SOM", "Carbon", "Nitrogen", "PAI")

lm.table

#######

#######
# making residuals and graphing those. 

bob <- lm(CVI_ABOVE$LN_FLUX ~ CVI_ABOVE$TEMP)  #puts the ln of flux first against Temp . . . i.e. conventional wisdom

summary(bob)

# pearson's r
cor.test(~CVI_ABOVE$LN_FLUX + CVI_ABOVE$TEMP)
cor.test(~CVI_ABOVE$LN_FLUX + CVI_ABOVE$WFPS)

# adding residuals to the df
CVI_ABOVE$RESID <-residuals(bob)


fitresiduals <-lm(CVI_ABOVE$RESID~CVI_ABOVE$WFPS)
summary(fitresiduals)

cor.test(~CVI_ABOVE$RESID + CVI_ABOVE$WFPS)

#Temp regressions for above and below
fitABOVE <- nls(FLUX ~ REGRESS(TEMP,a,b), data = CVI_ABOVE, start = c(a=0.5, b= 0.5), trace=T)
coeffABOVE <- coef(fitABOVE)

fitBELOW <- nls(FLUX ~ REGRESS(TEMP,a,b), data = CVI_BELOW, start = c(a=0.5, b= 0.5), trace=T)
coeffBELOW <- coef(fitBELOW)

#ABOVE and BELOW Temperature curve
sABOVE<-seq(from=min(CVI_ABOVE$TEMP,na.rm=T),to=max(CVI_ABOVE$TEMP,na.rm=T),length=50)
curveABOVE <- data.frame(TEMP=sABOVE,FLUX=REGRESS(sABOVE, a=coeffABOVE[1], b=coeffABOVE[2]))

sBELOW<-seq(from=min(CVI_BELOW$TEMP,na.rm=T),to=max(CVI_BELOW$TEMP,na.rm=T),length=50)
curveBELOW <- data.frame(TEMP=sBELOW,FLUX=REGRESS(sBELOW, a=coeffBELOW[1], b=coeffBELOW[2]))

#R^2 values
R2TEMP <-1-(deviance(fitTEMP)/sum((CVI$FLUX-mean(CVI$FLUX))^2))
R2ABOVE <-1-(deviance(fitABOVE)/sum((CVI_ABOVE$FLUX-mean(CVI_ABOVE$FLUX))^2))
R2BELOW <-1-(deviance(fitBELOW)/sum((CVI_BELOW$FLUX-mean(CVI_BELOW$FLUX))^2))

R2TEMP
R2ABOVE
R2BELOW

#################################################
# Flux, temp, and WFPS graphs

# FLUX against temperature with ABOVE and BELOW (11 degree) values color coded
pTEMP11 <- ggplot(CVI_ABOVE, aes(x=TEMP, y=FLUX)) +
  geom_point(shape=16, alpha=0.5, size=6, color="#FF0000") +
  
  geom_point(aes(x=TEMP,y=FLUX), data=CVI_BELOW, shape=16, size=6, alpha=0.5,color="#00A08A")+
  scale_x_continuous(limits=c(0.001,23), expand=c(0,0))+
  scale_y_continuous(limits=c(0.001,23), expand = c(0,0))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='black', size=24), legend.background = element_rect( fill="white"), 
        legend.key=element_rect(color="white", fill="white"), legend.title=element_blank() )+
  theme(panel.border = element_rect(colour="black", fill=NA))+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  theme(axis.ticks = element_blank())+
  xlab(templabel)+
  ylab( fluxlabel)+
  geom_line(data=curveABOVE, size=1, color="black", linetype= "solid")+
  geom_line(data=curveBELOW, size=1, color="black", linetype= "dotted")+
  geom_vline(aes(xintercept=11))

# ln(FLUX) values against soil temperature
pABOVE <- ggplot(CVI_ABOVE, aes(x=TEMP, y=LN_FLUX)) +
  geom_point(shape=16, alpha=0.5, size=7, color="#F98400") +
  scale_x_continuous(expand = c(0,0), limits=c(9,23), 
                     breaks = c(seq(from = 10, to = 22, by = 2)),
                     labels = c(seq(from = 10, to = 22, by = 2 )))+
  scale_y_continuous( limits=c(-1.5,3.99), expand = c(0,0))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='black', size=24), legend.background = element_rect( fill="white"), 
        legend.key=element_rect(color="white", fill="white"), legend.title=element_blank() )+
  theme(panel.border = element_rect(colour="black", fill=NA))+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  theme(axis.ticks = element_blank())+
  xlab(templabel)+
  ylab( lnfluxlabel)+
  stat_smooth(method=lm, se=FALSE, size=1, color="black")+
  geom_hline(aes(yintercept=0))

# residuals from ln(FLUX)/TEMP function against WFPS values
pRESID <- ggplot(CVI_ABOVE, aes(x=WFPS, y=RESID)) +
  geom_point(shape=16, alpha=0.5, size=7, color="#046C9A") +
  scale_y_continuous( limits=c(-2.5,2.5), expand = c(0,0))+
  theme_classic()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(color='black', size=24), legend.background = element_rect( fill="white"), 
        legend.key=element_rect(color="white", fill="white"), legend.title=element_blank() )+
  theme(panel.border = element_rect(colour="black", fill=NA))+
  theme(axis.title.x= element_text(size=20), 
        axis.text.x = element_text(size=(18)))+
  theme(axis.title.y= element_text(size=20), 
        axis.text.y = element_text(size=(18)))+
  theme(axis.ticks = element_blank())+
  labs(x=wfpslabel,
       y=residfluxlabel)+
  stat_smooth(method=lm, se=FALSE, size=1, color="black")+
  geom_hline(aes(yintercept=0))

pRESID


# making the major graph that is in the paper
grid.arrange(pTEMP11, pABOVE, pRESID)



###########################################################
##########################   BAR GRAPHS ###################

# FLUX by ELEVATION
CVI_p <-ddply(CVI, c("YEAR","ELEV"), summarise, 
                MEAN = mean(FLUX),
                SD = sd(FLUX),
                N = length(FLUX))

CVI_p$SE <- CVI_p$SD / (sqrt(CVI_p$N))

#significance labels
CVI_p$LABEL <- c("b","a","b","ab","ab","ab","ab","ab","ab")
CVI_p$LABEL2 <- c(NA,NA, NA,"B","A","AB",NA,NA,NA)

#Inputting SAS output data
CVI_p$LSMEAN<-c(6.3231,4.6993,6.1305,4.7621,4.7505,4.8257,4.71,4.4598,4.0434)
CVI_p$SE <- c(0.6688,0.6874,0.6918,0.5517,0.5717,0.5613,0.6812,0.7222,0.7027)


dodge=position_dodge(width=0.9)

p.MEANS <- ggplot( CVI_p, aes(x=toupper(ELEV), y = LSMEAN, fill=YEAR)) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("LOW","MID","HIGH"))+
  scale_fill_grey(start=1, end=.4)+
  scale_y_continuous(limits= c(0,8.5), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin=LSMEAN-SE, ymax=LSMEAN+SE), width=.2,position=dodge)+
  labs(y=fluxlabel)+
  labs(x="")+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=(18)))+
  geom_text(aes(y= 0.5, label=tolower(LABEL)),vjust=0, position=dodge)+
  geom_text(aes(y=1,label=LABEL2), vjust=0, position=dodge)

##################

# FLUX by VEGETATION


#making means

CVI_V <-ddply(CVI, c("YEAR","VEG"), summarise, 
              MEAN = mean(FLUX),
              SD = sd(FLUX),
              N = length(FLUX))

#significance labels
CVI_V$LABEL <- c("ab","ab","c","ab","b","ac","ab","ab","abc")
CVI_V$LABEL2 <- c(NA,NA,NA,"B","A","C",NA,NA,NA)

#Inputting SAS output data
CVI_V$LSMEAN<-c(5.1176,4.5466,7.4887,4.6841,4.0209,5.6332,4.3126,3.7752,5.1254)
CVI_V$SE <- c(0.6746,0.6854,0.6741,0.5572,0.5629,0.5591,0.6979,0.6982,0.705)


p.MEANSV <- ggplot( CVI_V, aes(x=toupper(VEG), y = LSMEAN, fill=YEAR)) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("OPEN","CANOPY","SHRUB"))+
  scale_fill_grey(start=1, end=.4)+
  scale_y_continuous(limits= c(0,8.5), expand=c(0,0))+
  theme_classic()+
  theme(legend.justification=c(1,1), 
        legend.position=c(1,1),
        legend.title=element_blank(),
        legend.text=element_text(color='black', size=24), 
        legend.background = element_rect( fill="#FFFFFF"), 
        legend.key=element_rect(color="#FFFFFF", fill="#FFFFFF") )+
  geom_errorbar(aes(ymin=LSMEAN-SE, ymax=LSMEAN+SE), width=.2,position=dodge)+
  labs(y="")+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=(18)))+
  geom_text(aes(y= 0.5, label=tolower(LABEL)),vjust=0, position=dodge)+
  geom_text(aes(y=1,label=LABEL2), vjust=0, position=dodge)


####################################################################################

# WFPS by ELEVATION

CVI_pW <-ddply(CVI, c("YEAR","ELEV"), summarise, 
               MEAN = mean(WFPS),
               SD = sd(WFPS),
               N = length(WFPS))

#Significance labels
CVI_pW$LABEL <- c("c","a","a","b","b","b","ac","ac","ab")
CVI_pW$LABEL2 <- c(NA,NA,NA,"B","A","A",NA,NA,NA)

#SAS least-squares means
CVI_pW$LSMEAN<-c(0.1413,0.1891,0.1844,0.2487,0.2472,0.2504,0.1839,0.1848,0.2060)
CVI_pW$SE <- c(0.01391,0.01431,0.01440,0.01188,0.01231,0.01208,0.01406,0.01482,0.01455)

p.WFPS_e <- ggplot( CVI_pW, aes(x=toupper(ELEV), y = LSMEAN, fill=YEAR)) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("LOW","MID","HIGH"))+
  scale_fill_grey(start=1, end=.4)+
  scale_y_continuous(limits= c(0,.36), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin=LSMEAN-SE, ymax=LSMEAN+SE), width=.2,position=dodge)+
  labs(y=wfpslabel)+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=20),
        axis.text.x = element_text(size=(18)),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=(18)))+
  geom_text(aes(y= 0.03, label=tolower(LABEL)),vjust=0, position=dodge)+
  geom_text(aes(y=0.05,label=LABEL2), vjust=0, position=dodge)



# WFPS by VEGETATION

CVI_VW <-ddply(CVI, c("YEAR","VEG"), summarise, 
               MEAN = mean(WFPS),
               SD = sd(WFPS),
               N = length(WFPS))

#Significance labels
CVI_VW$LABEL <- c("a","a","ab","c","b","c","ab","ab","ab")
CVI_VW$LABEL2 <- c(NA,NA,NA,"B","A","B",NA,NA,NA)

#SAS least-squares means
CVI_VW$LSMEAN<-c(0.1675,0.1614,0.1858,0.251,0.225,0.2703,0.1988,0.1876,0.1883)
CVI_VW$SE <- c(0.01403,0.01427,0.01402,0.012,0.01212,0.01204,0.01441,0.01434,0.01457)

p.WFPS_v<- ggplot( CVI_VW, aes(x=toupper(VEG), y = LSMEAN, fill=YEAR)) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("OPEN","CANOPY","SHRUB"))+
  scale_fill_grey(start=1, end=.4)+
  scale_y_continuous(limits= c(0,.36), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin=LSMEAN-SE, ymax=LSMEAN+SE), width=.2,position=dodge)+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_blank(),
        axis.text.x = element_text(size=(18)),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=(18)))+
  geom_text(aes(y= 0.03, label=tolower(LABEL)),vjust=0, position=dodge)+
  geom_text(aes(y=0.05,label=LABEL2), vjust=0, position=dodge)


#####################################################################

# TEMP by ELEVATION

CVI_pT <-ddply(CVI, c("YEAR","ELEV"), summarise, 
               MEAN = mean(TEMP),
               SD = sd(TEMP),
               N = length(TEMP))
#Significance labels
CVI_pT$LABEL <- c("ab","a","b","ab","a","ab","b","ab","b")
CVI_pT$LABEL2 <- c(NA,NA,NA,"C","A","B",NA,NA,NA)

#Least-squares means from SAS
CVI_pT$LSMEAN<-c(15.3011,16.2927,14.9018,15.5449,16.6144,15.31,13.9804,15.0877,13.9384)
CVI_pT$SE <- c(0.6546,0.6564,0.6568,0.5182,0.5201,0.5192,0.6561,0.659,0.658)

p.TEMP_e <- ggplot( CVI_pT, aes(x=toupper(ELEV), y = LSMEAN, fill=YEAR)) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("LOW","MID","HIGH"))+
  scale_fill_grey(start=1, end=.4)+
  scale_y_continuous(limits= c(0,20.5), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin=LSMEAN-SE, ymax=LSMEAN+SE), width=.2,position=dodge)+
  labs(y=templabel)+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=(18)))+
  geom_text(aes(y=3,label=LABEL2), vjust=0, position=dodge)

# TEMP by VEGETATION

CVI_VT <-ddply(CVI, c("YEAR","VEG"), summarise, 
               MEAN = mean(TEMP),
               SD = sd(TEMP),
               N = length(TEMP))
#Significance labels
CVI_VT$LABEL <- c("ab","ab","ab","ab","a","ab","b","ab","b")
CVI_VT$LABEL2 <- c(NA,NA,NA,"B","A","C",NA,NA,NA)

#Least-squares means from SAS output
CVI_VT$LSMEAN<-c(15.3969,15.6739,15.4248,15.762,16.3177,15.3896,14.1591,14.864,13.9833)
CVI_VT$SE <- c(0.6551,0.6561,0.6551,0.5187,0.5192,0.5189,0.6576,0.6568,0.6583)

p.TEMP_v<- ggplot( CVI_VT, aes(x=toupper(VEG), y = LSMEAN, fill=YEAR)) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("OPEN","CANOPY","SHRUB"))+
  scale_fill_grey(start=1, end=.4)+
  scale_y_continuous(limits= c(0,20.5), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin=LSMEAN-SE, ymax=LSMEAN+SE), width=.2,position=dodge)+
  theme(axis.title.x=element_blank(), 
        axis.title.y= element_blank(),
        axis.text.x = element_text(size=18),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=(18)))+
  geom_text(aes(y=3,label=LABEL2), vjust=0, position=dodge)


# Get the widths correct for arranging grids
# The maxWidth section still doesn't work, but the graph comes out nicely anyway.
gA <- ggplot_gtable(ggplot_build(p.MEANS))
gB <- ggplot_gtable(ggplot_build(p.MEANSV))
gC <- ggplot_gtable(ggplot_build(p.TEMP_e))
gD <- ggplot_gtable(ggplot_build(p.TEMP_v))
gE <- ggplot_gtable(ggplot_build(p.WFPS_e))
gF <- ggplot_gtable(ggplot_build(p.WFPS_v))

maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3], 
                     gC$widths[2:3], gD$widths[2:3]
                     gE$widths[2:3], gF$widths[2:3])

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth

# Arrange the four charts
grid.arrange(gA, gB, gC, gD, gE, gF, nrow=3)




##################################################################
##################    SOILS ANALYSIS AND PLOTS    ################
##################################################################


#Shapiro-Wilk normality test on data
shapiro.test(soil.data$som)
qqnorm(soil.data$som)

#Mann whitney U wilcoxon test

kruskal.test(som ~ elev, data=soil.data)
kruskal.test(som ~ veg, data=soil.data)
kruskal.test(som ~ depth, data=soil.data)

#subsetting by depth
soil5 <-subset(soil.data, soil.data$depth==5)
soil12 <-subset(soil.data, soil.data$depth==12)
soil20 <-subset(soil.data, soil.data$depth==20)

kruskal.test(som ~ elev, data=soil5)
kruskal.test(som ~ veg, data=soil5)
nrow(soil5)

kruskal.test(som ~ elev, data=soil12)
kruskal.test(som ~ veg, data=soil12)

kruskal.test(som ~ elev, data=soil20)
kruskal.test(som ~ veg, data=soil20)

#calculating differences in soil profiles
fivetotwelve <- merge(soil5, soil12, by="ID")

fivetotwelve$som512 <- ( ( (12*(fivetotwelve$soil.bulk.density.y*(fivetotwelve$som.y)/100)) - (5*(fivetotwelve$soil.bulk.density.x*(fivetotwelve$som.x/100))) )  /   ( (fivetotwelve$soil.bulk.density.y*12)-(fivetotwelve$soil.bulk.density.x*5)))*100 

twelvetotwenty <- merge(soil12, soil20, by="ID")
twelvetotwenty$som1220 <- ( ( (20*(twelvetotwenty$soil.bulk.density.y*(twelvetotwenty$som.y)/100)) - (12*(twelvetotwenty$soil.bulk.density.x*(twelvetotwenty$som.x/100))) )  /   ( (twelvetotwenty$soil.bulk.density.y*20)-(twelvetotwenty$soil.bulk.density.x*12)))*100 


#means
SOMmeans <- aggregate(formula= som~elev+veg, data=soil.data, FUN=mean)

SOMmeansveg <- aggregate(formula= som~veg, data=soil.data, FUN=mean)
SOMmeanselev <- aggregate(formula= som~elev, data=soil.data, FUN=mean)
SOMmeansdepth1 <- aggregate(formula= som~depth, data=soil.data, FUN=mean) 

SOMmeansdepth <- aggregate(formula= som~elev+veg+depth, data=soil.data, FUN=mean)
SOMmeansvegdepth <- aggregate(formula= som~veg+depth, data=soil.data, FUN=mean)
SOMmeanselevdepth <- aggregate(formula= som~elev+depth, data=soil.data, FUN=mean)

#shapiro testin' other things
shapiro.test(soils$soil.bulk.density)
qqnorm(soils$soil.bulk.density)
shapiro.test(soils$total.bulk.density)
qqnorm(soils$total.bulk.density)

#using way too much code to generate means of soil builk density by elev * depth and veg * depth

BULKsmeans <- aggregate(formula= soil.bulk.density~elev+veg, data=soil.data, FUN=mean)
BULKsmeansveg <- aggregate(formula= soil.bulk.density~veg, data=soil.data, FUN=mean)
BULKsmeanselev <- aggregate(formula= soil.bulk.density~elev, data=soil.data, FUN=mean)

BULKsmeansdepth <- aggregate(formula= soil.bulk.density~elev+veg+depth, data=soil.data, FUN=mean)
BULKsmeansvegdepth <- aggregate(formula= soil.bulk.density~veg+depth, data=soil.data, FUN=mean)
BULKsmeanselevdepth <- aggregate(formula= soil.bulk.density~elev+depth, data=soil.data, FUN=mean)

BULKsmeansvegdepth_sd <- aggregate(formula= soil.bulk.density~veg+depth, data=soil.data, FUN=sd)
BULKsmeanselevdepth_sd <- aggregate(formula= soil.bulk.density~elev+depth, data=soil.data, FUN=sd)
BULKsmeansvegdepth_sd <- rename(BULKsmeansvegdepth_sd, c("soil.bulk.density" = "sd"))
BULKsmeanselevdepth_sd <- rename(BULKsmeanselevdepth_sd, c("soil.bulk.density" = "sd"))

BULKsmeansvegdepth_n <- aggregate(formula= soil.bulk.density~veg+depth, data=soil.data, FUN=length)
BULKsmeanselevdepth_n <- aggregate(formula= soil.bulk.density~elev+depth, data=soil.data, FUN=length)
BULKsmeansvegdepth_n <- rename(BULKsmeansvegdepth_n, c("soil.bulk.density" = "n"))
BULKsmeanselevdepth_n <- rename(BULKsmeanselevdepth_n, c("soil.bulk.density" = "n"))


BULKsmeansvegdepth <-merge(BULKsmeansvegdepth, BULKsmeansvegdepth_sd, by=c("veg", "depth") )
BULKsmeansvegdepth <-merge(BULKsmeansvegdepth, BULKsmeansvegdepth_n, by=c("veg", "depth") )
BULKsmeanselevdepth <-merge(BULKsmeanselevdepth, BULKsmeanselevdepth_sd, by=c("elev", "depth") )
BULKsmeanselevdepth <-merge(BULKsmeanselevdepth, BULKsmeanselevdepth_n, by=c("elev", "depth") )


#standard error
BULKsmeansvegdepth$se <- (BULKsmeansvegdepth$sd/ (log(BULKsmeansvegdepth$n-1)))
BULKsmeanselevdepth$se <- (BULKsmeanselevdepth$sd/ (log(BULKsmeanselevdepth$n-1)))

#this renames CLOSED to CANOPY as it was mislabelled in the .csv
levels(BULKsmeansvegdepth$veg) <-c(levels(BULKsmeansvegdepth$veg), "CANOPY")
BULKsmeansvegdepth$veg[BULKsmeansvegdepth$veg == 'CLOSED'] <-"CANOPY"

####making bargraphs with standard error and junk
dodge=position_dodge(width=0.9)
p.BULKs.depth.elev <- ggplot(BULKsmeanselevdepth, aes(x=toupper(elev), y = soil.bulk.density, fill=as.factor(depth))) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("LOW","MID","HIGH"))+
  scale_fill_grey(start=1, end=.4, 
                  name="Soil Depth",
                  breaks=c("5", "12", "20"),
                  labels=c("5 cm", "12 cm", "20 cm"))+
  scale_y_continuous(limits= c(0,1.10), expand=c(0,0))+
  theme_classic()+
  theme(legend.position="none")+
  geom_errorbar(aes(ymin=soil.bulk.density-se, ymax=soil.bulk.density+se), width=.2,position=dodge)+
  labs(y=sbulklabel)+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=28),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))


p.BULKs.depth.veg <- ggplot(BULKsmeansvegdepth, aes(x=toupper(veg), y = soil.bulk.density, fill=as.factor(depth))) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("OPEN","CANOPY","SHRUB"))+
  scale_fill_grey(start=1, end=.4, 
                  name="Soil Depth",
                  breaks=c("5", "12", "20"),
                  labels=c("5 cm", "12 cm", "20 cm"))+
  scale_y_continuous(limits= c(0,1.10), expand=c(0,0))+
  theme_classic()+
  theme(legend.justification=c(1,1), 
        legend.position=c(1,1),
        legend.title=element_text(size=24),
        legend.text=element_text(color='black', size=24), 
        legend.background = element_rect( fill="#FFFFFF"), 
        legend.key=element_rect(color="#FFFFFF", fill="#FFFFFF") )+
  geom_errorbar(aes(ymin=soil.bulk.density-se, ymax=soil.bulk.density+se), width=.2,position=dodge)+
  labs(y=sbulklabel)+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=28),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))

#working with too much code to do means of SOM for elev*depth and veg *depth
SOMmeansvegdepth_sd <- aggregate(formula= som~veg+depth, data=soil.data, FUN=sd)
SOMmeanselevdepth_sd <- aggregate(formula= som~elev+depth, data=soil.data, FUN=sd)
SOMmeansvegdepth_sd <- rename(SOMmeansvegdepth_sd, c("som" = "sd"))
SOMmeanselevdepth_sd <- rename(SOMmeanselevdepth_sd, c("som" = "sd"))

SOMmeansvegdepth_n <- aggregate(formula= som~veg+depth, data=soil.data, FUN=length)
SOMmeanselevdepth_n <- aggregate(formula= som~elev+depth, data=soil.data, FUN=length)
SOMmeansvegdepth_n <- rename(SOMmeansvegdepth_n, c("som" = "n"))
SOMmeanselevdepth_n <- rename(SOMmeanselevdepth_n, c("som" = "n"))


SOMmeansvegdepth <-merge(SOMmeansvegdepth, SOMmeansvegdepth_sd, by=c("veg", "depth") )
SOMmeansvegdepth <-merge(SOMmeansvegdepth, SOMmeansvegdepth_n, by=c("veg", "depth") )
SOMmeanselevdepth <-merge(SOMmeanselevdepth, SOMmeanselevdepth_sd, by=c("elev", "depth") )
SOMmeanselevdepth <-merge(SOMmeanselevdepth, SOMmeanselevdepth_n, by=c("elev", "depth") )


#standard error
SOMmeansvegdepth$se <- (SOMmeansvegdepth$sd/ (log(SOMmeansvegdepth$n-1)))
SOMmeanselevdepth$se <- (SOMmeanselevdepth$sd/ (log(SOMmeanselevdepth$n-1)))

#this renames CLOSED to CANOPY as it was mislabelled in the .csv
levels(SOMmeansvegdepth$veg) <-c(levels(SOMmeansvegdepth$veg), "CANOPY")
SOMmeansvegdepth$veg[SOMmeansvegdepth$veg == 'CLOSED'] <-"CANOPY"

####  SOM plots

p.SOM.depth.elev <- ggplot(SOMmeanselevdepth, aes(x=toupper(elev), y = som, fill=as.factor(depth))) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("LOW","MID","HIGH"))+
  scale_fill_grey(start=1, end=.4, 
                  name="Soil Depth",
                  breaks=c("5", "12", "20"),
                  labels=c("5 cm", "12 cm", "20 cm"))+
  scale_y_continuous(limits= c(0,80), expand=c(0,0))+
  theme_classic()+
  theme(legend.position="none")+
  geom_errorbar(aes(ymin=som-se, ymax=som+se), width=.2,position=dodge)+
  labs(y="SOM (%)")+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=28),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))

p.SOM.depth.veg <- ggplot(SOMmeansvegdepth, aes(x=toupper(veg), y = som, fill=as.factor(depth))) + geom_bar(position=dodge, stat="identity", colour="black")+
  scale_x_discrete(limits=c("OPEN","CANOPY","SHRUB"))+
  scale_fill_grey(start=1, end=.4, 
                  name="Soil Depth",
                  breaks=c("5", "12", "20"),
                  labels=c("5 cm", "12 cm", "20 cm"))+
  scale_y_continuous(limits= c(0,80), expand=c(0,0))+
  theme_classic()+
  theme(legend.position="none")+
  geom_errorbar(aes(ymin=som-se, ymax=som+se), width=.2,position=dodge)+
  labs(y="SOM (%)")+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=28),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))

#making panel graph of soils
#trying another approach


# Get the widths
gG <- ggplot_gtable(ggplot_build(p.BULKs.depth.elev))
gH <- ggplot_gtable(ggplot_build(p.BULKs.depth.veg))
gI <- ggplot_gtable(ggplot_build(p.SOM.depth.elev))
gJ <- ggplot_gtable(ggplot_build(p.SOM.depth.veg))


maxWidth = unit.pmax(gH$widths[2:3],
                     gI$widths[2:3], 
                     gJ$widths[2:3],
                     gK$widths[2:3])
# Set the widths
gG$widths[2:3] <- maxWidth
gH$widths[2:3] <- maxWidth
gI$widths[2:3] <- maxWidth
gJ$widths[2:3] <- maxWidth


# Arrange the four charts
grid.arrange(gG, gH, gI, gJ, nrow=2)




###### The following are extra plots not used in the paper and focus on aggregate elevation/veg class means
 
#this section adds an ordering variable so that ggplot is tricked in to forcing the color gradient to come out right
SOMmeanselev$ORDER [SOMmeanselev$elev == "LOW"] <- "a"
SOMmeanselev$ORDER [SOMmeanselev$elev == "MID"] <- "b"
SOMmeanselev$ORDER [SOMmeanselev$elev == "HIGH"] <- "c"

levels(SOMmeansveg$veg) <-c(levels(SOMmeansveg$veg), "CANOPY")
SOMmeansveg$veg[SOMmeansveg$veg == 'CLOSED'] <-"CANOPY"

SOMmeansveg$ORDER [SOMmeansveg$veg == "OPEN"] <- "a"
SOMmeansveg$ORDER [SOMmeansveg$veg == "CANOPY"] <- "b"
SOMmeansveg$ORDER [SOMmeansveg$veg == "SHRUB"] <- "c"


p.SOMelev <- ggplot( SOMmeanselev, aes(x=elev, y =  som, fill=ORDER)) + geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=c("LOW","MID","HIGH"))+
  #scale_fill_manual(values=c("#FFFFFF", "#999999", "#222222"))+
  scale_fill_grey(start=1, end=.2)+
  scale_y_continuous(limits= c(0,50), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+  
  labs(y="SOM (%)")+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=32),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))

p.SOMelev

p.SOMveg<- ggplot( SOMmeansveg, aes(x=veg, y =  som, fill=ORDER)) + geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=c("OPEN","CANOPY","SHRUB"))+
  #scale_fill_manual(values=c("#FFFFFF", "#999999", "#222222"))+
  scale_fill_grey(start=1, end=.2)+
  scale_y_continuous(limits= c(0,50), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+  
  labs(y="SOM (%)")+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=32),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))

p.SOMveg

#soil bulk density section
BULKsmeanselev$ORDER [BULKsmeanselev$elev == "LOW"] <- "a"
BULKsmeanselev$ORDER [BULKsmeanselev$elev == "MID"] <- "b"
BULKsmeanselev$ORDER [BULKsmeanselev$elev == "HIGH"] <- "c"

levels(BULKsmeansveg$veg) <-c(levels(BULKsmeansveg$veg), "CANOPY")
BULKsmeansveg$veg[BULKsmeansveg$veg == 'CLOSED'] <-"CANOPY"

BULKsmeansveg$ORDER [BULKsmeansveg$veg == "OPEN"] <- "a"
BULKsmeansveg$ORDER [BULKsmeansveg$veg == "CANOPY"] <- "b"
BULKsmeansveg$ORDER [BULKsmeansveg$veg == "SHRUB"] <- "c"



p.BULKselev <- ggplot( BULKsmeanselev, aes(x=elev, y =  soil.bulk.density, fill=ORDER)) + geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=c("LOW","MID","HIGH"))+
  #scale_fill_manual(values=c("#FFFFFF", "#999999", "#222222"))+
  scale_fill_grey(start=1, end=.2)+
  scale_y_continuous(limits= c(0,1), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+  
  labs(y=sbulklabel)+
  #labs(y=expression(paste(rho, [s],"(g cm"^"3"*")")))+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=32),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))

p.BULKselev

p.BULKsveg<- ggplot( BULKsmeansveg, aes(x=veg, y =  soil.bulk.density, fill=ORDER)) + geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=c("OPEN","CANOPY","SHRUB"))+
  #scale_fill_manual(values=c("#FFFFFF", "#999999", "#222222"))+
  scale_fill_grey(start=1, end=.2)+
  scale_y_continuous(limits= c(0,1), expand=c(0,0))+
  theme_classic()+
  guides(fill=FALSE)+  
  labs(y=sbulklabel)+
  theme(axis.title.x=element_blank())+ 
  theme(axis.title.y= element_text(size=32),
        axis.text.x = element_text(size=(21)),
        axis.text.y = element_text(size=(21)))

p.BULKsveg




