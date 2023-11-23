# Stream Water Temperature
# Master code

# Author: Vaughn Grey
# Created: 16/05/2023
# Version for initial submission to Water Research: 10/07/2023
# Updated 05/09/2023 responding to reviewer comments
# Final version as per accepted manuscript on 6 October 2023 
# Final page proofs on 23 November 2023

#' Code parts:
#' Part A: Calculation of water temperature trends
#' Part B: Random Forest model: delta WT
#' Part C: Random Forest model: spatial site mean WT
#' Part D: Random Forest model: spatio-temporal WT 
#' Part E: Methods figures

###############
### Options ###
###############

options.saving <- TRUE # if true saves plots

# Part A
options.hm.plot <- TRUE
options.acf <- FALSE
options.gam.plots <- TRUE
options.gam.ind.plots <- TRUE

# Part B
options.date.start <- 1992  # sets the start date for analysis (must be after 1992)
options.date.end <- 2022     # sets the end date for analysis (must be 2018 or earlier)
options.pdp <- TRUE 

# Part C 
options.plot.vi.sp.RHS <- FALSE

# Part D
options.plot.vi.st.RHS <- FALSE

options(scipen = 9999)

#####################
### Load packages ###
#####################

library(nlme)
library(mgcv)
library(lubridate)
library(car)
library(sae)
library(ggplot2)
library(dplyr)
library(sf)
library(ozmaps)
library(ggspatial)
library(cowplot)
library(caret)
library(gtools)
library(dplyr)
library(lubridate)
library(randomForest)
library(Hmisc)
library(corrplot)
library(grid)

##################
### Load data ###
##################

### Load data from "Water_temp_dataset.R" 
setwd("~/OneDrive/R/")
load("Scripts/WaterResearch_Stream_Temp/data.wt.inputs.Q.RData") # load input data

#####################################################################################################################################################################
### PART A: Calculation of water temperature trends ###
#####################################################################################################################################################################

##########################
### Set up ###
##########################

trend.CI <- 0.95

########################################################
### Removing of sites with insufficient data record  ###
########################################################

### Data counts by year
hm_data <- wt
hm_data$year <- year(hm_data$Date)

hm_count <- hm_data %>%
  group_by(Site,year) %>%
  summarise(
    Count = n()
  )

### Removing of sites with insufficient data record to be used in calculations

# Cut sites that hadn't sampled by 1994
wt.n <- hm_data %>%
  group_by(Site,year) %>% 
  summarise(
    n = n()
  )

# no. samples in first 5 years
wt.first <- wt.n[wt.n$year <= 1996,]
wt.first <- wt.first %>% group_by(Site) %>% summarise(n.sam.first = sum(n))

# no. samples last 5  years
wt.last <- wt.n[wt.n$year >= 2017,]
wt.last <- wt.last %>% group_by(Site) %>% summarise(n.sam.last = sum(n))

wt.long <- left_join(hm_data, wt.first, by = "Site")
wt.long <- left_join(wt.long, wt.last, by = "Site")

wt.long <- wt.long[wt.long$n.sam.first > 6,] # Cut sites with less than 6 samples in first 5 years
wt.long <- wt.long[wt.long$n.sam.last >6,] # Cut sites with less than 6 samples in last 5 years

### Data gaps

# ID length of years with missing data
wt.long.gaps <- wt.long %>% group_by(Site) %>% summarise(n = length(unique(year)))

# remove sites with no data for greater than 10 years
wt.long.gaps <- wt.long.gaps[wt.long.gaps$n >= 20,]

wt.sites <- wt.long.gaps$Site

### Create site numbers
#data.wt <- data.input.deltaWT
#site.numbers <- unique(data.wt$Site)
site.numbers <- wt.sites
site.numbers <- site.numbers[order(site.numbers)]
site.numbers <- data.frame(Site = site.numbers)
site.numbers$S.num <- seq(1,length(site.numbers$Site),by=1)

# and text nudges for plotting
site.numbers$xnudge <- 0
site.numbers$ynudge <- 0
site.numbers$ynudge[5] <- 0.01
site.numbers$xnudge[31] <- -0.015
site.numbers$ynudge[31] <- 0.01
site.numbers$xnudge[33] <- 0.001
site.numbers$xnudge[34] <- -0.027
site.numbers$ynudge[41] <- -0.017
site.numbers$xnudge[44] <- -0.001
site.numbers$ynudge[47] <- 0.005
site.numbers$xnudge[48] <- 0.01
site.numbers$xnudge[70] <- -0.01
site.numbers$ynudge[73] <- 0.01


# Plots, data count by year
if (options.hm.plot == TRUE){
  
hm_count <- hm_count[hm_count$Site %in% wt.sites,]
hm_count$Count <- ifelse(hm_count$Count > 10, 10, hm_count$Count)
hm_count <- left_join(hm_count, site.numbers, by = "Site")
hm_count <- as.data.frame(hm_count)
hm_count$S.num <- as.factor(hm_count$S.num)

hm_plot <- ggplot(hm_count, aes(x = year, y = S.num, fill = Count))+
  geom_tile() +
 # scale_fill_gradient(low="white", high="darkMagenta", limits=c(0,10), oob = scales::squish) +
  scale_fill_distiller(palette = "Purples", direction = 1, limits=c(0,10), breaks = c(0,5,10), labels = c("  0"," 5",">10")) +
  
  #guides(fill="none") +
  theme(panel.background=element_blank()) +
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.text=element_text(colour="black", size=10)) + 
  theme(axis.title=element_text(colour="black", size=12)) +
  theme(axis.ticks=element_line(colour="black")) +
  theme(axis.text.x = element_text(hjust = 0, angle = 90)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(1992, 2021,by=1)) +
  scale_y_discrete(limits = rev) +
  xlab("Year") +
  ylab("Site") 

if (options.saving == TRUE) {
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Suppl_Fig_1_heat_map.pdf", hm_plot, width=18, height=18, unit="cm", dpi=1000)
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Suppl_Fig_1_heat_map.jpeg", hm_plot, width=18, height=18, unit="cm", dpi=330)
}
}

# Environment clean up
rm(hm_data, wt.first, wt.last, wt.long, wt.long.gaps, wt.n)

# Subset to only those required for further analysis
wt <- wt[wt$Site %in% wt.sites,]

########################
### WT data summary  ###
########################

# Summary of WT stats
wt.summary <- wt %>% group_by(Site) %>% summarise(
  n.samples = n(), 
  wt.median = median(Value_Numeric), wt.min = min(Value_Numeric), wt.max = max(Value_Numeric))

# Number of samples per year including years without any samples
wt.summary.peryear <- wt 
wt.summary.peryear$year <- year(wt.summary.peryear$Date)
wt.summary.peryear <- wt.summary.peryear %>% group_by(Site, year) %>% summarise(n.samp = n())

all.years <- seq(1992, 2021, by = 1)
all.years <- tidyr::crossing(all.years, unique(wt$Site))
colnames(all.years) <- c("year","Site")
wt.summary.peryear <- left_join(all.years, wt.summary.peryear, by = c("Site","year"))
wt.summary.peryear$n.samp <- ifelse(is.na(wt.summary.peryear$n.samp),0, wt.summary.peryear$n.samp)
wt.summary.peryear <- wt.summary.peryear %>% group_by(Site) %>% summarise(n.median = median(n.samp), n.min = min(n.samp), n.max = max(n.samp))

# Final table
wt.summary <- left_join(wt.summary, wt.summary.peryear, by = "Site")
wt.summary <- left_join(wt.summary, site.numbers[, c("Site","S.num")], by = "Site")

wt.summary[,3:6] <- round(wt.summary[,3:6], digits = 1)
wt.summary$wt.median <- format(wt.summary$wt.median, nsmall = 1)
wt.summary$wt.min <- format(wt.summary$wt.min, nsmall = 1)
wt.summary$wt.max <- format(wt.summary$wt.max, nsmall = 1)
wt.summary$n.median <- format(wt.summary$n.median, nsmall = 1)
wt.summary$n.min <- format(wt.summary$n.min, nsmall = 0)
wt.summary$n.max <- format(wt.summary$n.max, nsmall = 0)

wt.summary$n.sum <- paste0(wt.summary$n.median, " (",wt.summary$n.min,", ",wt.summary$n.max,")")
wt.summary$wt.sum <- paste0(wt.summary$wt.median, " (",wt.summary$wt.min,", ",wt.summary$wt.max,")")
wt.summary$Site <- wt.summary$S.num
wt.summary <- wt.summary[,c("Site","n.samples","n.sum","wt.sum")]
colnames(wt.summary) <- c("Site","Total number of samples","Median number of samples per year (min, max)",
                          "Median stream temperature observations (°C) (min, max)")

library(flextable)
flextable(wt.summary) %>% save_as_docx(path = "Outputs/WaterResearch_Stream_Temp/Table_S1.docx")


############################
### Preparation for GAM  ###
############################

cs.data <- wt[,c("Site","Date","Time","Value_Numeric")]

# Date
cs.data$Date <- as.Date(cs.data$Date)

# find decimal day of year since sampling began (including proportion of year)
cs.data$Year <- year(cs.data$Date) - 1992 +
  ( yday(cs.data$Date) / yday(as.Date(paste0(year(cs.data$Date), "-12-31"), format = "%Y-%m-%d") ) )

## Hour, and remove those with NA 
cs.data$Hour <- as.numeric(substr(cs.data$Time, start = 1, stop =2)) + (as.numeric(substr(cs.data$Time,4,5))/60)
cs.data <- cs.data[!is.na(cs.data$Hour),]

## Load "Site" as a factor
cs.data$Site <- as.factor(cs.data$Site)

## Fraction year
cs.data$fr_yr <- cs.data$Year - floor(cs.data$Year)


##################################
### Model 1: GAM for all sites ###
##################################

### fit GAM with splines, linear trend
c.gam.L <- mgcv::gam(Value_Numeric ~ Year + Site + s(fr_yr, bs="cc") + s(Hour, bs = "cc"), 
                     data = cs.data)

preds.L.30yrs <- c.gam.L$coefficients[[2]] * 30
print(preds.L.30yrs)

# calculate confidence intervals
preds.30yrs.lowerCI <- confint.lm(c.gam.L, level = trend.CI)[2,1] * 30
preds.30yrs.upperCI <- confint.lm(c.gam.L, level = trend.CI)[2,2] * 30

# calculate p-value and adj-r2
summary(c.gam.L)$r.sq
summary(c.gam.L)$s.table[1,4] #p-value
AIC(c.gam.L)
BIC(c.gam.L)

###############################################
### Impact of less samples on overall trend ###
###############################################

trend.test <- left_join(all.years, hm_count[,c("Site","year","Count")], by = c("Site","year"))
trend.test$Count <- ifelse(is.na(trend.test$Count),0, trend.test$Count)

wt.trend.test.store <- vector('list', length = 3)

for (i in 0:2){
  wt.trend.test <- trend.test %>% group_by(Site) %>% summarise(n.times = sum(Count <= 1)) # 1 or less samples
  wt.trend.test <- wt.trend.test[wt.trend.test$n.times <= i,] # occuring i or less times
  wt.trend.test <- unique(wt.trend.test$Site)
  
  # Now cut wt.sites to only those also in wt.sites.reducedsamples
  wt.trend.test <- cs.data[cs.data$Site %in% wt.trend.test,]
  
  ### fit GAM with splines, linear trend
  wt.trend.test.gam.L <- mgcv::gam(Value_Numeric ~ Year + Site + s(fr_yr, bs="cc") + s(Hour, bs = "cc"), 
                       data = wt.trend.test)
  
  wt.trend.test.preds.L.30yrs <- round(wt.trend.test.gam.L$coefficients[[2]] * 30, digits = 2)

  # calculate confidence intervals
  wt.trend.test.preds.30yrs.lowerCI <- round(confint.lm(wt.trend.test.gam.L, level = trend.CI)[2,1] * 30, digits = 2)
  wt.trend.test.preds.30yrs.upperCI <- round(confint.lm(wt.trend.test.gam.L, level = trend.CI)[2,2] * 30, digits = 2)
  
  wt.trend.test.preds <- data.frame(n.site = length(unique(wt.trend.test$Site)), trend = wt.trend.test.preds.L.30yrs, 
                                    t.LCI = wt.trend.test.preds.30yrs.lowerCI,
                                    t.UCI = wt.trend.test.preds.30yrs.upperCI)
  
  wt.trend.test.store[[i+1]] <- wt.trend.test.preds
  
}

wt.trend.test.summary <- do.call(rbind, wt.trend.test.store)

wt.trend.test.summary[4,1] <- length(unique(cs.data$Site))
wt.trend.test.summary[4,2] <- format(round(preds.L.30yrs, digits = 2), nsmall = 2)
wt.trend.test.summary[4,3] <- round(preds.30yrs.lowerCI, digits = 2)
wt.trend.test.summary[4,4] <- round(preds.30yrs.upperCI, digits = 2)

wt.trend.test.summary <- wt.trend.test.summary[c(4,1,2,3),]
wt.trend.test.summary$t.summary <- paste0(wt.trend.test.summary$trend, " [",wt.trend.test.summary$t.LCI,", ",wt.trend.test.summary$t.UCI,"]")

wt.trend.test.summary$Scenario <- c("All 75 sites", "0 years with 0 or 1 samples", "1 year with 0 or 1 samples", "2 years with 0 or 1 samples")
wt.trend.test.summary <- wt.trend.test.summary[,c("Scenario","n.site","t.summary")]
colnames(wt.trend.test.summary) <- c("Scenario","Number of sites","Trend")

flextable(wt.summary) %>% save_as_docx(path = "Outputs/WaterResearch_Stream_Temp/Table_S1.docx")


#######################
### Autocorrelation ###
#######################

if(options.acf == TRUE){

alpha <- 0.95
c.gamL.acf <- acf(c.gam.L$residuals, plot=FALSE)
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(c.gamL.acf$n.used)

acf.mod4.gamL <- acf(c.gam.L$residuals, plot=FALSE)$acf %>% 
  as_tibble() %>% mutate(lags = 1:n()) %>% 
  ggplot(aes(x=lags, y = V1)) + 
  #as_tibble() %>%
  #ggplot(aes(x=lag, y = V1)) + 
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  geom_hline(yintercept=0, col='black') +
  #scale_x_continuous(breaks=seq(0,10,1), limits = c(2,11), labels = seq(-1,9,1)) +
  #scale_y_continuous(limits = c(-0.05,0.1)) +
  labs(y="Autocorrelations", x="Lag", title= "ACF: GAM with linear trend") +
  #theme(axis.title.x=element_blank())+
  geom_segment(aes(xend=lags, yend=0)) +
  #geom_point() +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) 
#acf.mod4.gamL

acf.c <- acf(c.gam.L$residuals, plot=FALSE)$acf %>% 
  as_tibble() %>% mutate(lags = 1:n())

if (options.saving == TRUE) {
ggsave(plot = acf.mod4.gamL, filename = "Outputs/WaterResearch_Stream_Temp/acf_plot_GAM.pdf", width = 12, height = 12, dpi = 1000)
}
}

#####################
### Model 1 plots ###
#####################

if (options.gam.plots == TRUE){
  
  alllabel <- expression("Stream Temperature" ~ (degree*C))
  ylabel <- expression("Observed Stream Temperature" ~ (degree*C))
  xlabel <- expression("Modelled Stream Temperature" ~ (degree*C))
  
  ### Timeseries plot
  cs.data$c.mod.fv <- c.gam.L$fitted.values
  cs.data$pYear <- cs.data$Year + 1992
  
  cs.data$colourP <- "Observations"
  cs.data$colourL <- "GAM"
  cs.ind.colours <- c("Observations" = "orange","GAM" = "grey40")
  
  c.gam.L.plot.ts <- ggplot(data = cs.data) +
    geom_point(aes(x = Year, y = Value_Numeric, colour = colourP)) +
    geom_line(aes(x = Year, y = c.mod.fv, colour = colourL), linewidth = 0.1)+
    geom_abline(slope = c.gam.L$coefficients[[2]],
                intercept = c.gam.L$coefficients[[1]], colour = "grey40", linetype = "dashed", linewidth = 1) +
    scale_colour_manual(values = cs.ind.colours) +
    theme(panel.background = element_blank())+
    theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
    ggtitle(paste0("Stream temperature observations 1992 - 2021, fitted model and linear trend")) +
    xlab("Year")+
    ylab(alllabel) +
    scale_x_continuous(labels = c(1992,2002,2012,2022), breaks = c(0,10,20,30)) +
    labs(colour = element_blank()) +
    theme(legend.position = c(0.08,0.9), legend.direction = "vertical", 
          legend.spacing.y = unit(1, 'mm'))+
    theme(legend.title = element_text(size = 7), 
          legend.text  = element_text(size = 7),
          legend.key.size = unit(1, "lines"))+
    guides(colour = guide_legend(order = 1, reverse = T))
           
  #c.gam.L.plot.ts
  
  if (options.saving == TRUE) {
  ggsave(plot = c.gam.L.plot.ts, filename = paste0("Outputs/WaterResearch_Stream_Temp/Suppl_Fig_3_GAMTrend.pdf"), width = 18, height = 10, units = "cm",dpi = 1000)
  ggsave(plot = c.gam.L.plot.ts, filename = paste0("Outputs/WaterResearch_Stream_Temp/Suppl_Fig_3_GAMTrend.jpeg"), width = 18, height = 10, units = "cm",dpi = 330)
  }
  
  ### Observed vs modelled
  c.gam.L.plot.om <- ggplot() +
    geom_point(data = cs.data, aes(x = c.mod.fv,y = Value_Numeric), colour = "purple")+
    geom_abline(slope = 1,
                intercept = 0, colour = "blue", linetype = "dashed", linewidth = 1) +
    #geom_point(data = cs.mth.data,aes(x = Year, y = Value_Numeric), colour = "orange") +
    theme(panel.background = element_blank())+
    theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
    #  ggtitle(paste0("Observed vs Modelled Water Temperature: \nGAM with linear trend")) +
    ggtitle(paste0("GAM with linear trend")) +
    xlab(xlabel) +
    ylab(ylabel) +
    xlim(0,35)+
    ylim(0,35)
  #c.gam.L.plot.om
  
  #ggsave(plot = c.gam.L.plot.om, filename = paste0("Outputs/WaterResearch_Stream_Temp/WT_om_gam.pdf"), width = 8, height = 8, units = "cm",dpi = 1000)
  
}


#############################################
### Model 2: GAM for each individual site ###
#############################################

ind.sites <- unique(cs.data$Site)
W <- seq(1, length(ind.sites), by = 1)

wt.ind <- lapply(W, function(w){
  
  cs.ind.site <- ind.sites[[w]]
  cs.ind.data <- cs.data[cs.data$Site == cs.ind.site, ]
  
  c.gam.W <- mgcv::gam(Value_Numeric ~ Year + s(fr_yr, bs="cc") + s(Hour, bs = "cc"), 
                       data = cs.ind.data)
  
  # rate of change per decade
  cs.ind.preds <- data.frame(Site = cs.ind.site,CI = trend.CI)
  cs.ind.preds[1,"delta.dc"]  <- c.gam.W$coefficients[[2]] * 10
  cs.ind.preds[1,"wt.p1.delta.30y"] <- cs.ind.preds[1,"delta.dc"]  *  3  
  
  # calculate confidence intervals
  cs.ind.preds[1,"CI.low.delta.dc"] <- confint.lm(c.gam.W, level = trend.CI)[2,1] * 10
  cs.ind.preds[1,"CI.high.delta.dc"] <- confint.lm(c.gam.W, level = trend.CI)[2,2] * 10
  cs.ind.preds[1,"CI.low.delta.30y"] <- confint.lm(c.gam.W, level = trend.CI)[2,1] * 30
  cs.ind.preds[1,"CI.high.delta.30y"] <- confint.lm(c.gam.W, level = trend.CI)[2,2] * 30
  
  # model fit
  cs.ind.preds[1,"adj.r.squared"] <- summary(c.gam.W)$r.sq
  
  # find median WT value, using intercept and trend
  cs.ind.preds[1,"medianWT"] <- c.gam.W$coefficients[[1]] + c.gam.W$coefficients[[2]] * 15
  
  ### Plot
  
  ## Set up timeseries of prediction dates for plotting
  
  a.ts <- data.frame(TS = seq(0,30*365,by = 1/(24)))
  a.ts$Site <- cs.ind.site
  a.ts$Year <- floor(a.ts$TS) / 365
  a.ts$Hour <- (a.ts$TS - floor(a.ts$TS)) * 24
  a.ts$fr_yr <- a.ts$Year - floor(a.ts$Year)
  
  a.ts$pred <- predict(c.gam.W, a.ts)
  #a.ts$pred <- predict(c.gamm$gam, a.ts)
  
  cs.ind.data$c.mod2.fv <- c.gam.W$fitted.values
  
  cs.plot.W <- ggplot() +
    
    geom_line(data = a.ts, aes(x = Year, y = pred), colour = "blue", linewidth = 0.1, alpha = 0.3)+
    
    #geom_line(data = cs.ind.data, aes(x = Year, y = c.mod2.fv), colour = "blue", linewidth = 0.1)+
    geom_abline(slope = c.gam.W$coefficients[[2]],
                intercept = c.gam.W$coefficients[[1]], colour = "blue", linetype = "dashed", linewidth = 1) +
    
    geom_point(data = cs.ind.data,aes(x = Year, y = Value_Numeric), colour = "orange") +
    
    theme(panel.background = element_blank())+
    theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
    ggtitle(paste0("Timeseries 1992 - 2021 and GAM. Site = ", cs.ind.site)) +
    ylab("Stream Temp degC") 
  #xlim(10,13)
  #ylim(5,30)
  
  cs.plot.W
  
  #ggsave(file = paste0("Outputs/WaterResearch_Stream_Temp/GAM_",cs.ind.site,".pdf"), cs.plot.W, width=16, height=16, unit="cm", dpi=1000)
  print(w)
  return(cs.ind.preds)
  
})

wt.ind <- do.call(rbind, wt.ind)

################################
### Individual sites summary ###
################################

## Site mean water temperatures
wt.site.mean <- wt.ind[,c("Site","medianWT")]
colnames(wt.site.mean) <- c("Site","WT.site.mean")
#write.csv(x = wt.site.mean, file = "Outputs/WaterResearch_Stream_Temp/Site_mean_WT.csv", row.names = F)

## Site delta temp
wt.delta <- wt.ind
#write.csv(x = wt.delta, file = "Outputs/WaterResearch_Stream_Temp/delta_WT.csv", row.names = F)


#####################
### Model 2 plots ###
#####################

if (options.gam.ind.plots == TRUE){

#########################
### Plotting on table ###
#########################

    
wt.ind.p <- wt.ind
wt.ind.p <- wt.ind.p[,c("Site","wt.p1.delta.30y","CI.low.delta.30y","CI.high.delta.30y")]


wt.ind.p$pColour <- ifelse((wt.ind.p$CI.low.delta.30y) > 0, "positive",
                           ifelse((wt.ind.p$CI.high.delta.30y) < 0, "negative", "neutral")
)

wt.ind.p <- tidyr::pivot_longer(data = wt.ind.p, cols = c(2:4)) 

wt.ind.p.mean <- wt.ind.p[wt.ind.p$name == "wt.p1.delta.30y",]
wt.ind.p.mean <- wt.ind.p.mean[,c("Site","value")]
colnames(wt.ind.p.mean) <- c("Site","mean")

wt.ind.p <- left_join(wt.ind.p, wt.ind.p.mean, by = "Site")

wt.ind.p <- wt.ind.p[order(wt.ind.p$Site, decreasing = F),]

wt.ind.p <- left_join(wt.ind.p, site.numbers, by = "Site")
wt.ind.p$S.num <- as.factor(wt.ind.p$S.num)

wt.ind.colours <- c("positive" = "red","negative" = "#1976D2","neutral" = "grey")
label.wt <- expression("Stream temp. change " ~ (degree*C))

# now add in data for "all catchment"

wt.ind.p1 <- ggplot(data = wt.ind.p)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_line(aes(x = value, y = S.num, colour = pColour))+
  geom_point(aes(x = mean, y = S.num, colour = pColour)) +
  #xlim(-3.5,3.5) +
  scale_colour_manual(values = wt.ind.colours) +
  scale_y_discrete(limits=rev, position = "right")+
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(legend.position = "none")+
  # ggtitle("River temperature change 1992 - 2021") +
  ylab("Site") +
  xlab(label.wt) +
  theme(axis.text.y=element_text(size=4),axis.text.x=element_text(size=6),
        axis.title=element_text(size=7))

#wt.ind.p1

#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Fig_1_B_WT_ind_sites.pdf", wt.ind.p1, width=10, height=20, unit="cm", dpi=1000)
#write.csv(x = wt.ind.p, file = "Outputs/WaterResearch_Stream_Temp/WT_delta_indSites.csv", )

#######################
### Plotting on map ###
#######################

oz <- ozmap_data("states")
vic <- oz[oz$NAME == "Victoria",]

GIS.wt.catch <- GIS.wt.catch[GIS.wt.catch$SiteLT %in% wt.sites,]

GIS.wt.delta <- left_join(GIS.lt, wt.ind.p[wt.ind.p$name == "wt.p1.delta.30y",], by = "Site")
GIS.wt.delta <- GIS.wt.delta[!is.na(GIS.wt.delta$pColour),]
GIS.wt.delta$psize <- abs(GIS.wt.delta$value)

GIS.wt.sitenum <- left_join(GIS.wt.delta[,c("geometry","Site")], site.numbers, by = "Site")

cs.ind.colours <- c("positive" = "#F44336","negative" = "#1976D2","neutral" = "black")
cs.ind.fill <- c("positive" = "#F44336","negative" = "#1976D2","neutral" = "grey")
cs.ind.shapes <- c('positive'=24, 'negative'=25, 'neutral'=21)

label.size <- expression("Stream temp. change"~(degree*C))

fig.ind.trends <- ggplot() +
  geom_sf(data = vic, fill = NA, colour = "black") +
  #geom_sf(data = GIS.catch, fill = "white") +
  geom_sf(data = GIS.storages[GIS.storages$ASSET_NAME != "THOMSON RESERVOIR",], fill = "grey70", colour = NA ,linewidth = 0.1,alpha = 0.4) +
  geom_sf(data = GIS.streams, colour = "grey70", alpha = 0.4) +
  geom_sf(data = GIS.wt.catch, fill = "NA", colour = "black") +
  #geom_sf(data = GIS.wt, colour = "black") +
  geom_sf(data = GIS.wt.delta, aes(colour = pColour, fill = pColour, shape = pColour, size = psize)) +
  geom_sf_text(data = GIS.wt.sitenum, aes(label = S.num),nudge_x = GIS.wt.sitenum$xnudge, nudge_y = GIS.wt.sitenum$ynudge, 
               colour = "black", size = 1) +
  scale_colour_manual(values = wt.ind.colours) +
  scale_fill_manual(values = wt.ind.colours) +
  scale_shape_manual(values = cs.ind.shapes) +
  scale_size_continuous(range = c(1,4)) +
  guides(colour = guide_legend(order = 1, reverse = T),
         fill = guide_legend(order = 1, reverse = T),
         shape = guide_legend(order = 1, reverse = T),
         size = guide_legend(order = 2, reverse = T))+
  labs(colour = "Stream temp. change \n(direction)")+
  labs(fill = "Stream temp. change \n(direction)")+
  labs(shape = "Stream temp. change \n(direction)")+
  labs(size = label.size)+
#  theme(legend.position = "right", legend.spacing.y = unit(1, 'mm'))+
  theme(legend.position = c(0.898,0.22), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm'))+
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.55, "lines"))+
  theme(panel.background = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank())+
  ggspatial::annotation_scale(location = "bl", bar_cols = c("grey60", "white"), width_hint = 0.10) +
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.5, "cm"), pad_y = unit(0.75, "cm"),
                                    style = north_arrow_fancy_orienteering()) +
#  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.3), ylim = c(-37.1,-38.6), expand = F)
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.5), ylim = c(-37.1,-38.6), expand = F)+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) 
  
#fig.ind.trends

#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Map_ind_sites.pdf", fig.ind.trends, width=15, height=10, unit="cm", dpi=1000)


### Inset map of Australia
fig.oz <- ggplot() +
  geom_sf(data = oz, fill = "white")+
  geom_sf(data = GIS.catch, fill = "grey60", colour = "grey60") +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(plot.background = element_blank()) +
  coord_sf(xlim = c(112, 155), ylim = c(-10, -44), expand = F) +
#  ggspatial::annotation_scale(location = "bl", bar_cols = c("grey60", "white")) +
  geom_rect(aes(xmin = 143.8, xmax = 146.5,
                ymin = -36.9, ymax = -38.8),
            fill = NA, colour = "black", linewidth = 0.5)
#fig.oz

## Put together

fig.all <- cowplot::ggdraw(plot = fig.ind.trends) +
  draw_plot(plot = fig.oz,
            #            x = 0.5692, y = 0.741,
            #x = 0.6725, y = 0.70,
            x = 0.720, y = 0.70,
            width = 0.3, height = 0.3)
#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Fig_1_A_Map_ind_sites.pdf", fig.all, width=15, height=10, unit="cm", dpi=1000)

#fig1 <- cowplot::plot_grid(fig.all, wt.ind.p1,labels = c('A', 'B'), label_size = 12, rel_widths = c(2.6, 1))
#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Fig_1_Site_Trends.pdf", fig1, width=18, height=10, unit="cm", dpi=1000)

blank <- ggplot() + theme_void()

fig1 <- cowplot::ggdraw(plot = blank) +
  draw_plot(plot = fig.all,
            x = 0, y = 0.063, width = 0.72, height = 0.937) +
  draw_plot(plot = wt.ind.p1,
            x = 0.72, y = 0, width = 0.28, height = 1.0) +
  draw_plot_label(label = "(a)",fontface = "bold",x = 0.02, y = 0.98, hjust = -0.5, vjust = 1.5, size = 10) +
  draw_plot_label(label = "(b)",fontface = "bold",x = 0.73, y = 0.98, hjust = -0.5, vjust = 1.5, size = 10)

if (options.saving == TRUE) {
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_2_Site_Trends.pdf", fig1, width=18, height=10.58, unit="cm", dpi=1000)
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_2_Site_Trends.jpeg", fig1, width=18, height=10.58, unit="cm", dpi=330)
}


}



#####################################################################################################################################################################
### PART B: Random Forest model :: delta-WT ###
#####################################################################################################################################################################

################
### Delta-WT ###
################
data.wt.delta <- wt.delta
delta.wt <- data.wt.delta[,c("Site","delta.dc")]
delta.wt$delta.wt <- (delta.wt$delta.dc / 10) * (options.date.end - options.date.start)
delta.wt$delta.dc <- NULL  ## 
colnames(delta.wt) <- c("Site","data")

###############################################
### Select base variables for deltaWT model ###
###############################################

#data.wt <- data.wt.cov.all
#data.wt <- data.wt[,c("Site","Date","Year","precip.mm","tmax","srad","windSpeed","sm.Pc",
#                      "Urban","Ag","Forest")]

data.wt <- data.input.deltaWT
data.wt <- data.wt[data.wt$Site %in% wt.sites,] # subset to sites

###################################
### Calculate landcover metrics ###
###################################

### Calculate landcover in km2
carea <- data.input.spatial[,c("Site","carea_km2")]
data.wt <- left_join(data.wt, carea, by="Site")
data.wt$UrbanKm2 <- data.wt$Urban * data.wt$carea_km2
data.wt$AgKm2 <- data.wt$Ag * data.wt$carea_km2
data.wt$ForestKm2 <- data.wt$Forest * data.wt$carea_km2
data.wt$carea_km2 <- NULL

### Calculate landcover metrics
delta.preds <- data.wt 
delta.preds.L <- tidyr::pivot_longer(data = delta.preds, cols = c(4:length(colnames(delta.preds))))
delta.preds.L <- delta.preds.L[!is.na(delta.preds.L$value),]

delta.preds.L.urban <- delta.preds.L[delta.preds.L$name == "Urban",]
delta.preds.L.urban <- delta.preds.L.urban %>%
  group_by(Site) %>%
  slice_min(Date)
delta.preds.L.urban <- delta.preds.L.urban[,c("Site","value")]
colnames(delta.preds.L.urban) <- c("Site","UrbanPc.1992") 
delta.preds.L.urban <- left_join(delta.preds.L.urban, carea, by = "Site")
delta.preds.L.urban$UrbanKm2.1992 <- delta.preds.L.urban$UrbanPc.1992*delta.preds.L.urban$carea_km2

delta.preds.L.ag <- delta.preds.L[delta.preds.L$name == "Ag",]
delta.preds.L.ag <- delta.preds.L.ag %>%
  group_by(Site) %>%
  slice_min(Date)
delta.preds.L.ag <- delta.preds.L.ag[,c("Site","value")]
colnames(delta.preds.L.ag) <- c("Site","AgPc.1992")
delta.preds.L.ag <- left_join(delta.preds.L.ag, carea, by = "Site")
delta.preds.L.ag$AgKm2.1992 <- delta.preds.L.ag$AgPc.1992*delta.preds.L.ag$carea_km2

delta.preds.L.forest <- delta.preds.L[delta.preds.L$name == "Forest",]
delta.preds.L.forest <- delta.preds.L.forest %>%
  group_by(Site) %>%
  slice_min(Date)
delta.preds.L.forest <- delta.preds.L.forest[,c("Site","value")]
colnames(delta.preds.L.forest) <- c("Site","ForestPc.1992")
delta.preds.L.forest <- left_join(delta.preds.L.forest, carea, by = "Site")
delta.preds.L.forest$ForestKm2.1992 <- delta.preds.L.forest$ForestPc.1992*delta.preds.L.forest$carea_km2

## now pull together
data.preds.1992 <- left_join(delta.preds.L.urban, delta.preds.L.ag, by = c("Site","carea_km2"))
data.preds.1992 <- left_join(data.preds.1992, delta.preds.L.forest, by = c("Site","carea_km2"))

###################################
### Summarise changes over time ###
###################################

delta.preds.L$YearS <- delta.preds.L$Year

# find decimal day of year since sampling began (including proportion of year)
delta.preds.L$Year <- year(delta.preds.L$Date) - 1992 +
  ( yday(delta.preds.L$Date) / yday(as.Date(paste0(year(delta.preds.L$Date), "-12-31"), format = "%Y-%m-%d") ) )

## find fraction year
delta.preds.L$fr_yr <- delta.preds.L$Year - floor(delta.preds.L$Year)

#####################
### Loop 1 = by site
pred.sites <- unique(delta.preds.L$Site)
pred.sites.store <- vector('list',length = length(pred.sites))

for (ps in 1:length(pred.sites)){
  
  pred.site <- pred.sites[[ps]]
  pred.site.data <- delta.preds.L[delta.preds.L$Site == pred.site,]
  
  ### Loop 2 = by variable
  pred.vars <- unique(delta.preds.L$name)
  pred.vars.store <- vector('list', length = length(pred.vars))

  for (pv in 1:length(pred.vars)){
    
    pv.var <- pred.vars[[pv]]
    pv.data <- pred.site.data[pred.site.data$name == pv.var,]
    
    if (pv.var == "precip.mm"){
      
      # calculate monthly precip.mm and trend per year
      pv.data$pMonth <- month(pv.data$Date)
      pv.data$pYear <- year(pv.data$Date)
      pv.data<- pv.data %>% group_by(pYear, pMonth) %>% summarise(value = sum(value))
      pv.data$Date <- as.Date(paste0(pv.data$pYear,"-",pv.data$pMonth,"-15"), format = "%Y-%m-%d")
      pv.data$Year <- year(pv.data$Date) - 1992 +
        ( yday(pv.data$Date) / yday(as.Date(paste0(year(pv.data$Date), "-12-31"), format = "%Y-%m-%d") ) )
    }
    
    # linear models
    pv.lm <- lm(value ~ Year, data = pv.data)
    pv.lin <- pv.lm$coefficients[[2]]
    pred.vars.store[[pv]] <- pv.lin
    
  }
  
  pred.vars.out.d <- do.call(rbind, pred.vars.store)
  pred.vars.out <- data.frame(pred.vars,pred.vars.out.d)
  colnames(pred.vars.out) <- c("name","delta")
  pred.vars.out$Site <- pred.site
  
  pred.sites.store[[ps]] <- pred.vars.out
  
  print(ps)
}

# Processing
pred.sites.out <- do.call(rbind, pred.sites.store)
pred.sites.out$delta <- pred.sites.out$delta * (options.date.end - options.date.start)
pred.sites.out.W <- tidyr::pivot_wider(data =  pred.sites.out, names_from = "name", values_from = "delta")

# Then combine back into data.wt
delta.all <- left_join(delta.wt, pred.sites.out.W, by = "Site")
data.wt <- delta.all


#############################################
### Calculate changes to hour of sampling ###
#############################################
#delta.hour <- data.wt.temporal.noQ[,c("Site","Date","Value_Numeric","Time")]
delta.hour <- data.input.temporal.ante[,c("Site","Date","Value_Numeric","Time")]
delta.hour$Hour <- as.numeric(substr(delta.hour$Time, start = 1, stop = 2)) + (as.numeric(substr(delta.hour$Time, start = 4, stop = 5))/60)
delta.hour <- delta.hour[delta.hour$Hour < 23 & delta.hour$Hour > 5,]
delta.hour$Time <- NULL
delta.hour <- delta.hour[delta.hour$Site %in% wt.sites,]

# find decimal day of year since sampling began (including proportion of year)
delta.hour$Year <- year(delta.hour$Date) - 1992 +
  ( yday(delta.hour$Date) / yday(as.Date(paste0(year(delta.hour$Date), "-12-31"), format = "%Y-%m-%d") ) )

pred.sites <- unique(delta.hour$Site)
pred.sites.store <- vector('list',length = length(pred.sites))

for (ps in 1:length(pred.sites)){
  
  pred.site <- pred.sites[[ps]]
  pred.site.data <- delta.hour[delta.hour$Site == pred.site,]
  
    # linear models
    pv.lm <- lm(Hour ~ Year, data = pred.site.data)
    pv.lin <- pv.lm$coefficients[[2]]
    
  pred.vars.out <- data.frame("Hour",pv.lin)
  colnames(pred.vars.out) <- c("name","delta")
  pred.vars.out$Site <- pred.site
  
  pred.sites.store[[ps]] <- pred.vars.out
  
  print(ps)
}

# Processing
delta.hour <- do.call(rbind, pred.sites.store)
delta.hour$delta <- delta.hour$delta * (options.date.end - options.date.start)
delta.hour$name <- NULL
colnames(delta.hour) <- c("Hour","Site")

data.wt <- left_join(data.wt, delta.hour, by = "Site")

##############################################
### Combine temporal and spatial variables ###
##############################################

data.wt <- left_join(data.wt, data.preds.1992, by = "Site") # Landuse
spatial.pv <- data.input.spatial[,c("Site","DOR","AF","Elev","GradientPC","Aspect","siteannual.sm.Pc","mean.annual.precip.mm","SoilBulkDensity","ClayPC","MAQ_mm_SM")] # Add other spatial variables suitable to delta-WT
data.wt <- left_join(data.wt, spatial.pv, by = "Site")

# rename delta landuse variables
colnames(data.wt)[colnames(data.wt) == "Urban"] <- "deltaUrbanPc"
colnames(data.wt)[colnames(data.wt) == "UrbanKm2"] <- "deltaUrbanKm2"
colnames(data.wt)[colnames(data.wt) == "Ag"] <- "deltaAgPc"
colnames(data.wt)[colnames(data.wt) == "AgKm2"] <- "deltaAgKm2"
colnames(data.wt)[colnames(data.wt) == "Forest"] <- "deltaForestPc"
colnames(data.wt)[colnames(data.wt) == "ForestKm2"] <- "deltaForestKm2"

# add in EC and pH
data.wt <- left_join(data.wt, ec.delta[,c("Site","EC.delta")], by = c("Site"))
data.wt <- left_join(data.wt, ph.delta[,c("Site","ph.delta")], by = c("Site"))

########################
# Set up response and training data

rf.data <- data.wt # input data all 

rf.train <- data.wt # input predictor variables only
rf.train$data <- NULL 
rf.train$Site <- NULL

rf.response <- data.frame(data = data.wt[,c("data")]) # data (y) for training only

#save.image(file="Outputs/WaterResearch_Stream_Temp/RF_Inputs.RData") 

#######################################
### Random Forest model using caret ### 
#######################################

# Use tuneRF to identify mtry with lowest OOB
#mtry <- tuneRF(x = rf.train, y = rf.response$data, ntreeTry = 500, improve = 0.01, stepFactor = 1.5,
#               trace = F, plot = F)

set.seed(923)

mtry.rf <- 6 # manually set from results of above mtry test

# Fit RF model here
myControl <- trainControl(method="repeatedcv", 
                          number=5, 
                          repeats=10)

tunegrid.rf <- expand.grid(.mtry=mtry.rf)

c.mod.RF<-train(rf.train,
                rf.response$data,
                method = "rf",
                trControl=myControl,
                tuneGrid=tunegrid.rf,
                ntree= 500,
                preProc=NULL
                )

##################
### Model plot ###
##################

### Observed vs modelled
rf.plot.data <- rf.response
rf.plot.data$fv <- predict(c.mod.RF)
rf.plot.data$Site <- data.wt$Site

rf.plot.colours <- data.wt.delta
rf.plot.colours$WTtrend <- ifelse(rf.plot.colours$CI.low.delta.30y > 0, "positive",
                                  ifelse(rf.plot.colours$CI.high.delta.30y < 0, "negative", "neutral")
)
rf.plot.colours <- rf.plot.colours[,c("Site","WTtrend")]

rf.plot.data <- left_join(rf.plot.data, rf.plot.colours, by = "Site")
rf.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")

# Plot
rf.plot.om <- ggplot() +
  geom_point(data = rf.plot.data, aes(x = fv, y = data), colour = "purple")+
#  geom_point(data = rf.plot.data, aes(x = fv, y = data, colour = WTtrend))+
#  scale_colour_manual(values = rf.colours) +
  geom_abline(slope = 1,
              intercept = 0, colour = "blue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(legend.position = "none") +
  ggtitle(paste0("(a) Stream temperatures changes")) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), plot.title = element_text(size=15)) +
  xlab("Predicted stream temperature change (°C)") +
  ylab("Observed stream temperature change (°C)") +
  xlim(-2,2.3)+
  ylim(-2,2.3)

#rf.plot.om

#ggsave(plot = rf.plot.om, filename = paste0("Outputs/WaterResearch_Stream_Temp/deltaWT_Obvs_pred.pdf"),
#       width=14, height=14, unit="cm", dpi=1000)

#########################
### Model performance ###
#########################

print(c.mod.RF)

## Performance from caret
rf.rmse <- getTrainPerf(c.mod.RF)$TrainRMSE
rf.r2 <- getTrainPerf(c.mod.RF)$TrainRsquared
rf.mae <- getTrainPerf(c.mod.RF)$TrainMAE
rf.NSE <- hydroGOF::NSE(rf.plot.data$fv, rf.plot.data$data)
rf.RMSE.calc <- sqrt(mean((rf.plot.data$data - rf.plot.data$fv)^2))

rf.per <- data.frame(RMSE = rf.rmse, R2 = rf.r2, MAE = rf.mae, NSE = rf.NSE)

#print(rf.NSE)

###########################
### Variable importance ###
###########################

rf.varimp <- varImp(c.mod.RF, scale = T)
rf.varimp.df <- data.frame(rf.varimp$importance)
rf.varimp.df$Variable <- row.names(rf.varimp.df)
colnames(rf.varimp.df) <- c("Importance","Variable")
rf.varimp.df <- rf.varimp.df[,c("Variable","Importance")]
row.names(rf.varimp.df) <- NULL

rf.varimp.df <- rf.varimp.df[order(rf.varimp.df$Importance, decreasing = F),]

# Assign plot names to variable codes
var.names <- as.data.frame(c("windSpeed","sm.Pc","srad","precip.mm","tmax",
                             "UrbanKm2.1992","ForestKm2.1992","AgKm2.1992",
                             "UrbanKm2Med","ForestKm2Med","AgKm2Med",
                             "deltaUrbanKm2","deltaForestKm2","deltaAgKm2",
                             "UrbanKm2","ForestKm2","AgKm2",
                             "UrbanPc.1992","ForestPc.1992","AgPc.1992",
                             "UrbanPcMed","ForestPcMed","AgPcMed",
                             "deltaUrbanPc","deltaForestPc","deltaAgPc",
                             "UrbanPc","ForestPc","AgPc",
                             "carea_km2", "Elev", "AF","mean.tmax", 
                             "GradientPC","Aspect", "DOR",
                             "siteannual.sm.Pc","ClayPC","SoilBulkDensity",
                             "mean.srad","mean.annual.precip.mm", "siteannualmean.windSpeed","Hour",
                             "tmax.1day","tmax.2day","tmax.3day","tmax.7day","tmax.14day","tmax.30day","tmax.60day",
                             "precip.mm.1day","precip.mm.2day","precip.mm.3day","precip.mm.7day","precip.mm.14day","precip.mm.30day","precip.mm.60day",
                             "srad.1day","srad.2day","srad.3day","srad.7day","srad.14day","srad.30day","srad.60day",
                             "sm.Pc.1day","sm.Pc.2day","sm.Pc.3day","sm.Pc.7day","sm.Pc.14day","sm.Pc.30day","sm.Pc.60day",
                             "windSpeed.1day","windSpeed.2day","windSpeed.3day","windSpeed.7day","windSpeed.14day","windSpeed.30day","windSpeed.60day",
                             "Time",
                             "EC.delta","ph.delta", "EC.site.mean", "ph.site.mean", "EC", "pH",
                             "MAQ_mm","MAQ_mm_SM"
                             ))
colnames(var.names) <- "Variable"

var.all <- rf.varimp.df$Variable
var.missing <- var.all[!(var.all %in% var.names$Variable)]

var.names$Vname <- c("Wind speed (change, m/s)","Soil moisture (change, %)", "Solar radiation (change, MJ/m2)", "Precipitation (change, mm)","Air temperature (change, °C)",
                     "Urban landcover (km2, 1992)","Forest landcover (km2, 1992)",  "Agricultural landcover (km2, 1992)",
                     "Urban landcover (km2, median)","Forest landcover (km2, median)",  "Agricultural landcover (km2, median)",
                     "Urban landcover (change, km2)","Forest landcover (change, km2)",  "Agricultural landcover (change, km2)",
                     "Urban landcover (km2)","Forest landcover (km2)",  "Agricultural landcover (km2)",
                     "Urban landcover (%, 1992)","Forest landcover (%, 1992)",  "Agricultural landcover (%, 1992)",
                     "Urban landcover (%, median)","Forest landcover (%, median)",  "Agricultural landcover (%, median)",
                     "Urban landcover (change, %)","Forest landcover (change, %)",  "Agricultural landcover (change, %)",
                     "Urban landcover (%)","Forest landcover (%)",  "Agricultural landcover (%)",
                     "Catchment area (km2)", "Elevation (m AHD)", "Attenuated forest cover (%)","Air temperature (site mean)",
                     "Gradient (%)", "Aspect (°)", "Degree of regulation",
                     "Soil moisture (site mean, %)","Soil clay content (%)", "Soil bulk density (g/cm3)",
                     "Solar radiation (site mean)","Precipitation (annual site mean, mm)", "Wind speed (site mean)","Time of day",
                     "Air temperature (1 day antecedent)","Air temperature (mean 2 days antecedent)","Air temperature (mean 3 days antecedent)","Air temperature (mean 7 days antecedent)","Air temperature (mean 14 days antecedent)","Air temperature (mean 30 days antecedent)","Air temperature (mean 60 days antecedent)",
                     "Precipitation (1 day antecedent)","Precipitation (mean 2 days antecedent)","Precipitation (mean 3 days antecedent)","Precipitation (mean 7 days antecedent)","Precipitation (mean 14 days antecedent)","Precipitation (mean 30 days antecedent)","Precipitation (mean 60 days antecedent)",
                     "Solar radiation (1 day antecedent)","Solar radiation (mean 2 days antecedent)","Solar radiation (mean 3 days antecedent)","Solar radiation (mean 7 days antecedent)","Solar radiation (mean 14 days antecedent)","Solar radiation (mean 30 days antecedent)","Solar radiation (mean 60 days antecedent)",
                     "Soil moisture (1 day antecedent)","Soil moisture (mean 2 days antecedent)","Soil moisture (mean 3 days antecedent)","Soil moisture (mean 7 days antecedent)","Soil moisture (mean 14 days antecedent)","Soil moisture (mean 30 days antecedent)","Soil moisture (mean 60 days antecedent)",
                     "Wind speed (1 day antecedent)","Wind speed (mean 2 days antecedent)","Wind speed (mean 3 days antecedent)","Wind speed (mean 7 days antecedent)","Wind speed (mean 14 days antecedent)","Wind speed (mean 30 days antecedent)","Wind speed (mean 60 days antecedent)",
                     "Time of day",
                     "Electical conductivity (change)","pH (change)", "Electical conductivity (site mean)","pH (site mean)","Electrical conductivity","pH",
                     "Total annual discharge (change, mm)","Total annual discharge (mm, site mean)"
                     )    


var.names$Influence <- c("Atmosphere","Stream","Atmosphere","Stream","Atmosphere",
                      "Landcover","Landcover","Landcover",
                      "Landcover","Landcover","Landcover",
                      "Landcover","Landcover","Landcover",
                      "Landcover","Landcover","Landcover",
                      "Landcover","Landcover","Landcover",
                      "Landcover","Landcover","Landcover",
                      "Landcover","Landcover","Landcover",
                      "Landcover","Landcover","Landcover",
                      "Topographic","Topographic","Topographic","Atmosphere",
                      "Topographic","Topographic","Topographic",
                      "Stream","Topographic","Topographic",
                      "Atmosphere","Stream","Atmosphere", "Atmosphere",
                      "Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere",
                      "Stream","Stream","Stream","Stream","Stream","Stream","Stream",
                      "Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere",
                      "Stream","Stream","Stream","Stream","Stream","Stream","Stream",
                      "Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere","Atmosphere",
                      "Atmosphere",
                      "Stream","Stream","Stream","Stream","Stream","Stream",
                      "Stream","Stream"
                      )

rf.varimp.df <- left_join(rf.varimp.df, var.names, by = "Variable")

#    colnames(rf.varimp.df)[colnames(rf.varimp.df) == "Variable"] <- c("Vname")

rf.varimp.df$pImp <- ifelse(rf.varimp.df$Importance < 1.0, round(rf.varimp.df$Importance, digits = 1),
                            round(rf.varimp.df$Importance, digits = 0))

#rf.varimp.df$pImp <- round(rf.varimp.df$Importance * 100, digits = 0)

rf.varimp.df$Vname <- factor(rf.varimp.df$Vname, levels = rf.varimp.df$Vname)

rf.fill <- c("Landcover" = "#edae49","Stream" = "#00798c","Atmosphere" = "grey","Topographic" = "#66a182")

rf.varimp.plot <- ggplot(data = rf.varimp.df, aes(x = Importance, y = Vname, fill = Influence))+
  geom_col() +
  scale_fill_manual(values = rf.fill) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(axis.title.y = element_blank()) +
  geom_text(aes(label = pImp), hjust = -0.1, size = 2)+
#  xlab("% Increase MSE") +
#  xlim(0,4) +
   #xlab("Variable relative importance") +
  theme(axis.title.x = element_blank()) +
 # xlim(0,112) +
  scale_x_continuous(breaks = c(0,25,50,75,100), limits = c(0,105)) +
  ggtitle("(a) Changes in stream temperature") +
  theme(legend.position = c(0.82,0.20), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm'), legend.key.size = unit(3,'mm'))+
  theme(axis.text=element_text(size=4), axis.title=element_text(size=5),
        plot.title = element_text(size=7),
        legend.text = element_text(size = 4), legend.title = element_text(size = 5)
        )
  
#rf.varimp.plot

ggsave(plot = rf.varimp.plot, filename = paste0("Outputs/WaterResearch_Stream_Temp/Fig_3_A_deltaWT_VarImpPlot.pdf"),
       width=9, height=7, unit="cm", dpi=1000)




################################
### Partial dependence plots ###
################################

label.pdp <- expression("Marginal effect" ~ (degree*C))


if (options.pdp == TRUE){
  
  ### Extract PDP
  library(iml)
  predictor.rf <- Predictor$new(model = c.mod.RF, data = rf.train, y = rf.response)
  
  rf.varimp.df <- rf.varimp.df[order(rf.varimp.df$Importance, decreasing = T),]
  colnames(rf.varimp.df)[colnames(rf.varimp.df) == "Vname"] <- c("Variable")
  pdp.preds.list <- rf.varimp.df$Variable
  pdp.preds.list <- as.character(pdp.preds.list)
  
  pdp.preds.p.store <- vector('list', length = length(pdp.preds.list))
  #pdp.preds.p.store <- vector('list', length = 2)

  for (pp in 1:length(pdp.preds.list)){
    #  for (pp in 1:3){   
    
    pdp.pred <- pdp.preds.list[[pp]]  
    pdp.pred.long <- var.names[var.names$Variable == pdp.pred,2]
    
    ## PDP plotting
    rf.pdp <- FeatureEffect$new(predictor.rf, pdp.pred, method = "pdp", grid.size = 30)
    
    pdp.test <- rf.pdp$results

    rf.pdp.plot <- rf.pdp$plot() + 
      #ylim(0,1.1) +
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")+
      theme(panel.background = element_blank())+
      theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
      theme(legend.position = "none")+
      #      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1.10)) +
      scale_y_continuous(breaks = c(0.25,0.50,0.75), limits = c(0.2,0.85)) +
      ylab(label.pdp) +
      xlab(pdp.pred.long) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=5),
            plot.title = element_text(size=7),
            legend.text = element_text(size = 5), legend.title = element_text(size = 5)
      )

    rf.pdp.plot
    
    pdp.preds.p.store[[pp]] <- rf.pdp.plot
    
  }
  
  # Save all pdp's
  pdp.preds.p.all <-  do.call(gridExtra::grid.arrange,pdp.preds.p.store)
  if (options.saving == TRUE) {
  ggsave(plot = pdp.preds.p.all, filename = paste0("Outputs/WaterResearch_Stream_Temp/deltaWT_PDP_plots.pdf"), 
         width=24, height=24, unit="cm", dpi=1000)
  
  ggsave(plot = pdp.preds.p.all, filename = paste0("Outputs/WaterResearch_Stream_Temp/deltaWT_PDP_plots.jpeg"), 
         width=24, height=24, unit="cm", dpi=330)
  }
  
  ######################################
  ### Select PDPs of interest (Figure 3)
  pdp.preds.list <- c("UrbanKm2.1992","deltaUrbanKm2","tmax","MAQ_mm")
  pdp.preds.p.store <- vector('list', length = length(pdp.preds.list))
  label.pdp2 <- bquote(paste("Marginal effect (°C)"))
  label.pdp3 <- bquote(paste(" "))
  
  var.names$pdpVname <- var.names$Vname
  var.names$pdpVname <- gsub("(","\n(",var.names$pdpVname, fixed = T)
  
  for (pp in 1:length(pdp.preds.list)){
    #  for (pp in 1:3){   
    if(pp == 1){label.pdp4 <- label.pdp2}
    if(pp > 1){label.pdp4 <- label.pdp3}
    
    pdp.pred <- pdp.preds.list[[pp]]  
    pdp.pred.long <- var.names[var.names$Variable == pdp.pred,4]
    
    ## PDP plotting
    rf.pdp <- FeatureEffect$new(predictor.rf, pdp.pred, method = "pdp", grid.size = 30)
    
    pdp.test <- rf.pdp$results
    
    rf.pdp.plot <- rf.pdp$plot() + 
      #ylim(0,1.1) +
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")+
      theme(panel.background = element_blank())+
      theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
      theme(legend.position = "none")+
      #      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1.10)) +
      scale_y_continuous(breaks = c(0.25,0.50,0.75), limits = c(0.2,0.85)) +
      theme(axis.title.y = element_blank()) +
      xlab(pdp.pred.long) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=5),
            plot.title = element_text(size=7),
            legend.text = element_text(size = 5), legend.title = element_text(size = 5))
    
    rf.pdp.plot
    
    pdp.preds.p.store[[pp]] <- rf.pdp.plot
    
  }
  
  # select most important
  rf.pdp.preds.top <- gridExtra::grid.arrange(grobs = pdp.preds.p.store[], ncol = 4, 
                                              left = textGrob("Changes to stream \ntemperature (°C)", gp=gpar(fontsize=5), rot = 90 ) )
  #ggsave(plot = rf.pdp.preds.top, filename = paste0("Outputs/WaterResearch_Stream_Temp/Fig_3_A_deltaWT_PDP_plots_top.pdf"), 
  #       width=18, height=4.5, unit="cm", dpi=1000)
  
}


# PDP air temp & urban area

#pdp.pred <- c("UrbanKm2.1992","tmax")  
pdp.pred <- c("tmax","UrbanKm2.1992")  

## PDP plotting
rf.pdp <- FeatureEffect$new(predictor.rf, pdp.pred, method = "pdp", grid.size = 30)

pdp.test <- rf.pdp$results

pdp.y.lab <- expression("Urban landcover " ~ (km^2))
pdp.x.lab <- expression("Air temperature change " ~ (degree*C))
pdp.labs.lab <- bquote(paste("Increases \nin stream \ntemperature (°C)"))

pdp.site.labs <- data.wt
pdp.site.labs <- left_join(pdp.site.labs, site.numbers, by = "Site")
pdp.site.labs <- pdp.site.labs[,c("S.num","UrbanKm2.1992","tmax")]
pdp.site.labs <- pdp.site.labs[pdp.site.labs$S.num %in% c(4,5,10,7,74,75,50,18),]
#pdp.site.labs <- pdp.site.labs[pdp.site.labs$S.num %in% c(4,5,10,7,74,75,50,18,31,35,42,12,23),]

pdp.hm <- pdp.test
colnames(pdp.hm) <- c("tmax","UrbanKm2.1992", "PDPValue","Type")

rf.pdp.plot <- rf.pdp$plot() + 
  #ggplot()+
  geom_tile(data = pdp.hm, aes(x = tmax, y = UrbanKm2.1992, fill = PDPValue, colour = PDPValue)) + 
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  labs(fill = pdp.labs.lab)+
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0.3,0.91), labels = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  scale_colour_distiller(palette = "Reds", direction = 1, limits = c(0.3,0.91), labels = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9), guide = "none") +
  xlab(pdp.x.lab) +
  ylab(pdp.y.lab) +
  geom_text(data = pdp.site.labs, aes(y = UrbanKm2.1992, x = tmax, label = S.num), size = 2)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5),
        plot.title = element_text(size=7),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5)
  )

rf.pdp.plot

if (options.saving == TRUE) {
ggsave(plot = rf.pdp.plot, filename = paste0("Outputs/WaterResearch_Stream_Temp/Figure_5_deltaWT_AT_Urban_PDP_plot.pdf"), 
       width=9, height=6.5, unit="cm", dpi=1000)

ggsave(plot = rf.pdp.plot, filename = paste0("Outputs/WaterResearch_Stream_Temp/Figure_5_deltaWT_AT_Urban_PDP_plot.jpeg"), 
       width=9, height=6.5, unit="cm", dpi=1000)
}

## PDP plotting alternate (for interest) = delta Urban and delta AT
pdp.pred <- c("tmax","deltaUrbanKm2")  

rf.pdp <- FeatureEffect$new(predictor.rf, pdp.pred, method = "pdp", grid.size = 30)

pdp.test <- rf.pdp$results

pdp.y.lab <- expression("Urban landcover change " ~ (km^2))
pdp.x.lab <- expression("Air temperature change " ~ (degree*C))
pdp.labs.lab <- bquote(paste("Increases \nin stream \ntemperature (°C)"))

pdp.site.labs <- data.wt
pdp.site.labs <- left_join(pdp.site.labs, site.numbers, by = "Site")
pdp.site.labs <- pdp.site.labs[,c("S.num","deltaUrbanKm2","tmax")]
pdp.site.labs <- pdp.site.labs[pdp.site.labs$S.num %in% c(4,5,10,7,74,75,50,18),]
#pdp.site.labs <- pdp.site.labs[pdp.site.labs$S.num %in% c(4,5,10,7,74,75,50,18,31,35,42),]

rf.pdp.plot <- rf.pdp$plot() + 
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  labs(fill = pdp.labs.lab)+
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = "grey60", limits = c(0.3,0.91), labels = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  xlab(pdp.x.lab) +
  ylab(pdp.y.lab) +
  geom_text(data = pdp.site.labs, aes(y = deltaUrbanKm2, x = tmax, label = S.num))

rf.pdp.plot

#ggsave(plot = rf.pdp.plot, filename = paste0("Outputs/WaterResearch_Stream_Temp/deltaWT_AT_deltaUrban_PDP_plot.pdf"), 
#       width=14, height=10, unit="cm", dpi=1000)



#############################################
### RF output plotting :: air temperature ###
#############################################

delta.tmax <- data.wt
delta.tmax <- delta.tmax[,c("Site","data","tmax")]

GIS.at.delta <- left_join(GIS.lt, delta.tmax, by = "Site")
GIS.at.delta$data <- abs(GIS.at.delta$data)
GIS.at.delta <- GIS.at.delta[!is.na(GIS.at.delta$tmax),]

GIS.at.delta <- left_join(GIS.at.delta, rf.plot.colours, by = "Site")
GIS.at.delta <- left_join(GIS.at.delta, site.numbers, by = "Site")

rf.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")
rf.shapes <- c('positive'=24, 'negative'=25, 'neutral'=21)

GIS.at.delta$WTtrend <- as.factor(GIS.at.delta$WTtrend)
levels(GIS.at.delta$WTtrend) <- c("negative","neutral","positive")

label.size <- expression("Water temp. change " ~ (degree*C))
label.fill <- expression("Air temp. change" ~ (degree*C))

fig.at.trends <- ggplot() +
  geom_sf(data = vic, fill = NA, colour = "black") +
  geom_sf(data = GIS.catch, fill = "white") +
  geom_sf(data = GIS.storages[GIS.storages$ASSET_NAME != "THOMSON RESERVOIR",], fill = "grey70", colour = NA ,linewidth = 0.1,alpha = 0.4) +
  geom_sf(data = GIS.streams, colour = "grey70", alpha = 0.4) +
  geom_sf(data = GIS.at.delta, aes(fill = tmax,shape = WTtrend, size = data), colour = "grey50") +
  geom_sf_text(data = GIS.wt.sitenum, aes(label = S.num),nudge_x = GIS.wt.sitenum$xnudge, nudge_y = GIS.wt.sitenum$ynudge, 
               colour = "black", size = 1) +
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = "grey60", limits = c(1.7,2.31)) +
  scale_colour_distiller(palette = "Reds", direction = 1, na.value = "grey60") +
  scale_shape_manual(values = rf.shapes) +
  scale_size_continuous(range = c(1,4)) +
  labs(fill = label.fill)+
  labs(shape = "Water temp. change \n(direction)")+
  labs(size = label.size)+
  #theme(legend.position = "right", legend.spacing.y = unit(1, 'mm'))+
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.75, "lines"))+
  guides(shape = guide_legend(order = 1, reverse = T),
         size = guide_legend(order = 2, reverse = T)) +
  #  theme(legend.position = "none")+
  theme(panel.background = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())+
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.8), ylim = c(-37.1,-38.6), expand = F)+
  theme(legend.position = c(0.89,0.5), legend.spacing.y = unit(1, 'mm'))+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) 
  
#fig.at.trends

#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Map_airtemp.pdf", fig.at.trends, width=15, height=10, unit="cm", dpi=1000)


###################################
### RF output plotting :: urban ###
###################################

delta.urb <- data.wt
delta.urb <- delta.urb[,c("Site","data","UrbanKm2.1992")]

GIS.at.delta <- left_join(GIS.lt, delta.urb, by = "Site")
GIS.at.delta$data <- abs(GIS.at.delta$data)
GIS.at.delta <- GIS.at.delta[!is.na(GIS.at.delta$UrbanKm2.1992),]
GIS.at.delta <- GIS.at.delta[order(GIS.at.delta$UrbanKm2.1992),]

GIS.at.delta$UrbanKm2.1992 <- ifelse(GIS.at.delta$UrbanKm2.1992 > 30, 30, GIS.at.delta$UrbanKm2.1992)

GIS.at.delta <- left_join(GIS.at.delta, rf.plot.colours, by = "Site")
GIS.at.delta <- left_join(GIS.at.delta, site.numbers, by = "Site")

rf.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")
rf.shapes <- c('positive'=24, 'negative'=25, 'neutral'=21)

GIS.at.delta$WTtrend <- as.factor(GIS.at.delta$WTtrend)
levels(GIS.at.delta$WTtrend) <- c("negative","neutral","positive")

label.size <- expression("Water temp. change " ~ (degree*C))
label.fill <- expression("Urban area" ~ (km^2))

fig.Urban.trends <- ggplot() +
  geom_sf(data = vic, fill = NA, colour = "black") +
  geom_sf(data = GIS.catch, fill = "white") +
  geom_sf(data = GIS.storages[GIS.storages$ASSET_NAME != "THOMSON RESERVOIR",], fill = "grey70", colour = NA ,linewidth = 0.1,alpha = 0.4) +
  geom_sf(data = GIS.streams, colour = "grey70", alpha = 0.4) +
  geom_sf(data = GIS.at.delta, aes(fill = UrbanKm2.1992,shape = WTtrend,  size = data), colour = "grey50") +
  geom_sf_text(data = GIS.wt.sitenum, aes(label = S.num),nudge_x = GIS.wt.sitenum$xnudge, nudge_y = GIS.wt.sitenum$ynudge, 
               colour = "black", size = 1) +
  scale_fill_distiller(palette = "PuRd", direction = 1, na.value = "grey60", labels = c("  0"," 10"," 20",">30")) +
  scale_colour_distiller(palette = "PuRd", direction = 1, na.value = "grey60") +
  scale_shape_manual(values = rf.shapes) +
  scale_size_continuous(range = c(1,4)) +
  labs(fill = label.fill)+
  labs(shape = "Water temp. change \n(direction)")+
  labs(size = label.size)+
#  theme(legend.position = "right", legend.spacing.y = unit(1, 'mm'))+
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.75, "lines"))+
  guides(shape = guide_legend(order = 1, reverse = T),
         size = guide_legend(order = 2, reverse = T)) +
  theme(panel.background = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())+
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.8), ylim = c(-37.1,-38.6), expand = F)+
  theme(legend.position = c(0.89,0.5), legend.spacing.y = unit(1, 'mm'))+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) 

#fig.Urban.trends
#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Map_Urban.pdf", fig.Urban.trends, width=15, height=10, unit="cm", dpi=1000)


###########################################
### RF output plotting :: precipitation ###
###########################################

delta.precip <- data.wt
delta.precip <- delta.precip[,c("Site","data","precip.mm")]

GIS.at.delta <- left_join(GIS.lt, delta.precip, by = "Site")
GIS.at.delta$data <- abs(GIS.at.delta$data)
GIS.at.delta <- GIS.at.delta[!is.na(GIS.at.delta$precip.mm),]
GIS.at.delta <- GIS.at.delta[order(GIS.at.delta$precip.mm),]

#GIS.at.delta$UrbanKm2.1992 <- ifelse(GIS.at.delta$UrbanKm2.1992 > 30, 30, GIS.at.delta$UrbanKm2.1992)

GIS.at.delta <- left_join(GIS.at.delta, rf.plot.colours, by = "Site")
GIS.at.delta <- left_join(GIS.at.delta, site.numbers, by = "Site")

rf.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")
rf.shapes <- c('positive'=24, 'negative'=25, 'neutral'=21)

GIS.at.delta$WTtrend <- as.factor(GIS.at.delta$WTtrend)
levels(GIS.at.delta$WTtrend) <- c("negative","neutral","positive")

label.size <- expression("Water temp. change " ~ (degree*C))
label.fill <- expression("Precipitation change (mm)")

fig.Precip.trends <- ggplot() +
  geom_sf(data = vic, fill = NA, colour = "black") +
  geom_sf(data = GIS.catch, fill = "white") +
  geom_sf(data = GIS.storages[GIS.storages$ASSET_NAME != "THOMSON RESERVOIR",], fill = "grey70", colour = NA ,linewidth = 0.1,alpha = 0.4) +
  geom_sf(data = GIS.streams, colour = "grey70", alpha = 0.4) +
  geom_sf(data = GIS.at.delta, aes(fill = precip.mm,shape = WTtrend,  size = data), colour = "grey50") +
  geom_sf_text(data = GIS.wt.sitenum, aes(label = S.num),nudge_x = GIS.wt.sitenum$xnudge, nudge_y = GIS.wt.sitenum$ynudge, 
               colour = "black", size = 1) +
  #scale_fill_distiller(palette = "Reds", direction = -1, na.value = "grey60", limits = c(-15.01,0.1),labels = c("-15","-10","-5","0")) +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "grey60", limits = c(-15.01,15.01),labels = c("-15","-10","-5","0","5","10","15")) +
  scale_colour_distiller(palette = "Reds", direction = -1, na.value = "grey60") +
  scale_shape_manual(values = rf.shapes) +
  scale_size_continuous(range = c(1,4)) +
  labs(fill = label.fill)+
  labs(shape = "Water temp. change \n(direction)")+
  labs(size = label.size)+
  #  theme(legend.position = "right", legend.spacing.y = unit(1, 'mm'))+
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.75, "lines"))+
  guides(shape = guide_legend(order = 1, reverse = T),
         size = guide_legend(order = 2, reverse = T)) +
  theme(panel.background = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())+
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.8), ylim = c(-37.1,-38.6), expand = F)+
  theme(legend.position = c(0.89,0.5), legend.spacing.y = unit(1, 'mm'))+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) 

#fig.Precip.trends
#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Map_Precip.pdf", fig.Precip.trends, width=15, height=10, unit="cm", dpi=1000)

###########################################
### RF output plotting :: TAQ ###
###########################################

delta.taq <- data.wt
delta.taq <- delta.taq[,c("Site","data","MAQ_mm")]

GIS.at.delta <- left_join(GIS.lt, delta.taq, by = "Site")
GIS.at.delta$data <- abs(GIS.at.delta$data)
GIS.at.delta <- GIS.at.delta[!is.na(GIS.at.delta$MAQ_mm),]
GIS.at.delta <- GIS.at.delta[order(GIS.at.delta$MAQ_mm),]

#GIS.at.delta$UrbanKm2.1992 <- ifelse(GIS.at.delta$UrbanKm2.1992 > 30, 30, GIS.at.delta$UrbanKm2.1992)

GIS.at.delta <- left_join(GIS.at.delta, rf.plot.colours, by = "Site")
GIS.at.delta <- left_join(GIS.at.delta, site.numbers, by = "Site")

rf.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")
rf.shapes <- c('positive'=24, 'negative'=25, 'neutral'=21)

GIS.at.delta$WTtrend <- as.factor(GIS.at.delta$WTtrend)
levels(GIS.at.delta$WTtrend) <- c("negative","neutral","positive")

label.size <- expression("Water temp. change " ~ (degree*C))
label.fill <- expression("Total annual discharge \nchange (mm)")

fig.TAQ.trends <- ggplot() +
  geom_sf(data = vic, fill = NA, colour = "black") +
  geom_sf(data = GIS.catch, fill = "white") +
  geom_sf(data = GIS.storages[GIS.storages$ASSET_NAME != "THOMSON RESERVOIR",], fill = "grey70", colour = NA ,linewidth = 0.1,alpha = 0.4) +
  geom_sf(data = GIS.streams, colour = "grey70", alpha = 0.4) +
  geom_sf(data = GIS.at.delta, aes(fill = MAQ_mm,shape = WTtrend,  size = data), colour = "grey50") +
  geom_sf_text(data = GIS.wt.sitenum, aes(label = S.num),nudge_x = GIS.wt.sitenum$xnudge, nudge_y = GIS.wt.sitenum$ynudge, 
               colour = "black", size = 1) +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "grey60", limits = c(-140,140),breaks = c(-120,-80,-40,0,40,80,120),labels = c("-120","-80","-40","0","40","80","120")) +
  scale_colour_distiller(palette = "RdBu", direction = 1, na.value = "grey60") +
  scale_shape_manual(values = rf.shapes) +
  scale_size_continuous(range = c(1,4)) +
  labs(fill = label.fill)+
  labs(shape = "Water temp. change \n(direction)")+
  labs(size = label.size)+
  #  theme(legend.position = "right", legend.spacing.y = unit(1, 'mm'))+
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.75, "lines"))+
  guides(shape = guide_legend(order = 1, reverse = T),
         size = guide_legend(order = 2, reverse = T)) +
  theme(panel.background = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())+
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.8), ylim = c(-37.1,-38.6), expand = F)+
  theme(legend.position = c(0.89,0.5), legend.spacing.y = unit(1, 'mm'))+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) 

#fig.Urban.trends
#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Map_SoilM.pdf", fig.SM.trends, width=15, height=10, unit="cm", dpi=1000)

###########################################
### RF output plotting :: soil moisture ###
###########################################

delta.sm <- data.wt
delta.sm <- delta.sm[,c("Site","data","sm.Pc")]

GIS.at.delta <- left_join(GIS.lt, delta.sm, by = "Site")
GIS.at.delta$data <- abs(GIS.at.delta$data)
GIS.at.delta <- GIS.at.delta[!is.na(GIS.at.delta$sm.Pc),]
GIS.at.delta <- GIS.at.delta[order(GIS.at.delta$sm.Pc),]

#GIS.at.delta$UrbanKm2.1992 <- ifelse(GIS.at.delta$UrbanKm2.1992 > 30, 30, GIS.at.delta$UrbanKm2.1992)

GIS.at.delta <- left_join(GIS.at.delta, rf.plot.colours, by = "Site")
GIS.at.delta <- left_join(GIS.at.delta, site.numbers, by = "Site")

rf.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")
rf.shapes <- c('positive'=24, 'negative'=25, 'neutral'=21)

GIS.at.delta$WTtrend <- as.factor(GIS.at.delta$WTtrend)
levels(GIS.at.delta$WTtrend) <- c("negative","neutral","positive")

label.size <- expression("Water temp. change " ~ (degree*C))
label.fill <- expression("Soil moisture change (%)")

fig.SM.trends <- ggplot() +
  geom_sf(data = vic, fill = NA, colour = "black") +
  geom_sf(data = GIS.catch, fill = "white") +
  geom_sf(data = GIS.storages[GIS.storages$ASSET_NAME != "THOMSON RESERVOIR",], fill = "grey70", colour = NA ,linewidth = 0.1,alpha = 0.4) +
  geom_sf(data = GIS.streams, colour = "grey70", alpha = 0.4) +
  geom_sf(data = GIS.at.delta, aes(fill = sm.Pc,shape = WTtrend,  size = data), colour = "grey50") +
  geom_sf_text(data = GIS.wt.sitenum, aes(label = S.num),nudge_x = GIS.wt.sitenum$xnudge, nudge_y = GIS.wt.sitenum$ynudge, 
               colour = "black", size = 1) +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "grey60", limits = c(-0.06,0.06),breaks = c(-0.06,-0.03,0,0.03,0.06),labels = c("-0.06","-0.03","0","0.03","0.06")) +
  scale_colour_distiller(palette = "RdBu", direction = 1, na.value = "grey60") +
  scale_shape_manual(values = rf.shapes) +
  scale_size_continuous(range = c(1,4)) +
  labs(fill = label.fill)+
  labs(shape = "Water temp. change \n(direction)")+
  labs(size = label.size)+
  #  theme(legend.position = "right", legend.spacing.y = unit(1, 'mm'))+
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  theme(legend.title = element_text(size = 5), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.75, "lines"))+
  guides(shape = guide_legend(order = 1, reverse = T),
         size = guide_legend(order = 2, reverse = T)) +
  theme(panel.background = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())+
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.8), ylim = c(-37.1,-38.6), expand = F)+
  theme(legend.position = c(0.89,0.5), legend.spacing.y = unit(1, 'mm'))+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) 

#fig.Urban.trends
#ggsave(file = "Outputs/WaterResearch_Stream_Temp/Map_SoilM.pdf", fig.SM.trends, width=15, height=10, unit="cm", dpi=1000)

###################################
### Figures 5 and Supple. Fig 2 ###
###################################


fig5 <- cowplot::ggdraw(plot = blank) +
  draw_plot(plot = fig.at.trends,
            x = 0, y = 0.5, width = 1.0, height = 0.5) +
  draw_plot(plot = fig.Urban.trends,
            x = 0, y = 0.0, width = 1.0, height = 0.5) +
  draw_plot_label(label = "(a)",fontface = "bold",x = 0.04, y = 0.98, hjust = -0.5, vjust = 1.0, size = 9) +
  draw_plot_label(label = "(b)",fontface = "bold",x = 0.04, y = 0.48, hjust = -0.5, vjust = 1.0, size = 9)

if (options.saving == TRUE) {
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_6_AT_Urban.pdf", fig5, width=15, height=20, unit="cm", dpi=1000)
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_6_AT_Urban.jpeg", fig5, width=15, height=20, unit="cm", dpi=330)
}


fig.suppl.2 <- cowplot::ggdraw(plot = blank) +
  draw_plot(plot = fig.Precip.trends,
            x = 0, y = 0.66, width = 1.0, height = 0.33) +
  draw_plot(plot = fig.TAQ.trends,
            x = 0, y = 0.33, width = 1.0, height = 0.33) +
  draw_plot(plot = fig.SM.trends,
            x = 0, y = 0.0, width = 1.0, height = 0.33) +
  draw_plot_label(label = "(a)", fontface = "bold",x = 0.06, y = 0.97, hjust = 0, vjust = 0, size = 9) +
  draw_plot_label(label = "(b)", fontface = "bold",x = 0.06, y = 0.635, hjust = 0, vjust = 0, size = 9) +
  draw_plot_label(label = "(c)", fontface = "bold",x = 0.06, y = 0.305, hjust = 0, vjust = 0, size = 9)

if (options.saving == TRUE) {
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Suppl_Fig_4_Precip_SM.pdf", fig.suppl.2, width=15, height=30, unit="cm", dpi=1000)
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Suppl_Fig_4_Precip_SM.jpeg", fig.suppl.2, width=15, height=30, unit="cm", dpi=330)
}

##########################
### Cross correlations ###
##########################

# Spearman's rank correlation coefficient between delta-WT and predictors

data.xc <- data.wt
data.xc$Site <- NULL
data.xc$data <- NULL

# Option B : extracts p value but doesn't allow plotting
corr.matrix.B <- rcorr(as.matrix(data.xc), type = "spearman")

# Function to flatten matrix (taken from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

corr.B <- flattenCorrMatrix(corr.matrix.B$r, corr.matrix.B$P)

corr.B$abscor <- abs(corr.B$cor)
corr.B <- corr.B[order(corr.B$abscor, decreasing = T),]

### Summary stats
data.wt.summary <- tidyr::pivot_longer(data = data.wt, cols = c(2:30))
data.wt.summary$Site <- NULL
data.wt.summary <- data.wt.summary %>% group_by(name) %>% summarise(median = median(value), min = min(value), max = max(value), mean = mean(value))
colnames(data.wt.summary) <- c("Variable","Median","Min","Max","Mean")
data.wt.summary <- left_join(data.wt.summary, var.names)
data.wt.summary <- data.wt.summary[,c("Vname","Median","Min","Max","Mean")]
data.wt.summary <- data.wt.summary[!is.na(data.wt.summary$Vname),]

#####################################################################################################################################################################
### PART C: Random Forest model :: spatial, site mean WT ###
#####################################################################################################################################################################

############################################
### Calculate medians for landcover data ###
############################################

#data.landcover <- data.wt.temporal.noQ[,c("Site","Date","Time","Year","Urban","Ag","Forest")]
data.landcover <- data.input.deltaWT[,c("Site","Year","Urban","Ag","Forest")]

data.landcover <- data.landcover[data.landcover$Site %in% wt.sites,]
data.landcover <- data.landcover %>% 
  group_by(Site,Year) %>%
  summarise(
    UrbanPcMed = median(Urban),
    AgPcMed = median(Ag),
    ForestPcMed = median(Forest)
  )

data.landcover <- data.landcover[data.landcover$Year == 2008,]
data.landcover$Year <- NULL

data.landcover <- left_join(data.landcover, data.wt[,c("Site","carea_km2")], by = "Site")

data.landcover$UrbanKm2Med <- data.landcover$UrbanPcMed * data.landcover$carea_km2
data.landcover$AgKm2Med <- data.landcover$AgPcMed * data.landcover$carea_km2
data.landcover$ForestKm2Med <- data.landcover$ForestPcMed * data.landcover$carea_km2
data.landcover$carea_km2 <- NULL 

###################################
### Spatial predictor variables ###
###################################

#data.sp <- data.wt[,c("Site","carea_km2","DOR","AF","Elev","GradientPC","Aspect","mean.tmax")]

#data.sp <- data.wt.spatial[,c("Site","carea_km2","DOR","AF","Elev","GradientPC","Aspect","mean.tmax","mean.srad",
#                              "siteannualmean.windSpeed", "siteannual.sm.Pc","mean.annual.precip.mm","SoilBulkDensity",
#                              "ClayPC")]

data.sp <- data.input.spatial

data.sp <- data.sp[data.sp$Site %in% wt.sites,]

###################################
### Combine predictor variables ###
###################################

data.sp <- left_join(data.sp, data.landcover, by = "Site") # add in landcover variables
data.sp <- left_join(wt.site.mean, data.sp, by = "Site") # add in site mean WT, adding now to ensure same order in sp.train / sp.response below
colnames(data.sp)[colnames(data.sp) == "WT.site.mean"] <- "data"

# add in EC and pH
data.sp <- left_join(data.sp, ec.delta[,c("Site","EC.site.mean")], by = c("Site"))
data.sp <- left_join(data.sp, ph.delta[,c("Site","ph.site.mean")], by = c("Site"))


### Summary stats
data.sp.summary <- tidyr::pivot_longer(data = data.sp, cols = c(2:21))
data.sp.summary$Site <- NULL
data.sp.summary <- data.sp.summary %>% group_by(name) %>% summarise(median = median(value), min = min(value), max = max(value), mean = mean(value))
colnames(data.sp.summary) <- c("Variable","Median","Min","Max","Mean")
data.sp.summary <- left_join(data.sp.summary, var.names)
data.sp.summary <- data.sp.summary[,c("Vname","Median","Min","Max","Mean")]
data.sp.summary <- data.sp.summary[!is.na(data.sp.summary$Vname),]

data.pv.summary <- rbind(data.wt.summary, data.sp.summary)
data.pv.summary <- data.pv.summary[!duplicated(data.pv.summary),]
data.pv.summary[,c(2:4)] <- round(data.pv.summary[,c(2:4)], digits = 4)
data.pv.summary <- data.pv.summary[order(data.pv.summary$Vname),]                                  

data.tod.summary <- data.input.temporal.ante[,c("Site","Date","Value_Numeric","Time")]
data.tod.summary$Hour <- as.numeric(substr(data.tod.summary$Time, start = 1, stop = 2)) + (as.numeric(substr(data.tod.summary$Time, start = 4, stop = 5))/60)
data.tod.summary <- data.tod.summary[data.tod.summary$Hour < 23 & data.tod.summary$Hour > 5,]
data.tod.summary$Time <- NULL
data.tod.summary <- data.tod.summary[data.tod.summary$Site %in% wt.sites,]

data.tod.summary <- data.tod.summary %>% group_by(Site) %>% summarise(median = median(Hour))
data.tod.summary.median <- median(data.tod.summary$median)


########################
# Set up response and training data

sp.train <- data.sp # input predictor variables only
sp.train$data <- NULL 
sp.train$Site <- NULL

sp.response <- data.frame(data = data.sp[,c("data")]) # data (y) for training only

#######################################
### Random Forest model using caret ### 
#######################################

# Use tuneRF to identify mtry with lowest OOB
#mtry <- tuneRF(x = sp.train, y = sp.response$data, ntreeTry = 500, improve = 0.01, stepFactor = 1.5,
#               trace = F, plot = F)

mtry.sp <- 13 # manually set from results of above mtry test

# Fit RF model here
myControl <- trainControl(method="repeatedcv", 
                          number=5, 
                          repeats=10)

tunegrid.sp <- expand.grid(.mtry=mtry.sp)

sp.mod.RF<-train(sp.train,
                sp.response$data,
                method = "rf",
                trControl=myControl,
                tuneGrid=tunegrid.sp,
                ntree= 500,
                preProc=NULL
                )


##################
### Model plot ###
##################

### Observed vs modelled
sp.plot.data <- sp.response
sp.plot.data$fv <- predict(sp.mod.RF)
sp.plot.data$Site <- data.wt$Site

sp.plot.colours <- data.wt.delta
sp.plot.colours$WTtrend <- ifelse(sp.plot.colours$CI.low.delta.30y > 0, "positive",
                                  ifelse(sp.plot.colours$CI.high.delta.30y < 0, "negative", "neutral")
)
sp.plot.colours <- sp.plot.colours[,c("Site","WTtrend")]

sp.plot.data <- left_join(sp.plot.data, sp.plot.colours, by = "Site")
sp.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")

# Plot
sp.plot.om <- ggplot() +
  geom_point(data = sp.plot.data, aes(x = fv, y = data), colour = "purple")+
#  geom_point(data = sp.plot.data, aes(x = fv, y = data, colour = WTtrend))+
#  scale_colour_manual(values = rf.colours) +
  geom_abline(slope = 1,
              intercept = 0, colour = "blue", linetype = "dashed", linewidth = 1) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(legend.position = "none") +
  ggtitle(paste0("(c) Site mean stream temperature")) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), plot.title = element_text(size=15)) +
  xlab("Predicted site mean stream temperature (°C)") +
  ylab("Observed site mean stream temperature (°C)") +
  xlim(10,18)+
  ylim(10,18)

#sp.plot.om

#ggsave(plot = sp.plot.om, filename = paste0("Outputs/WaterResearch_Stream_Temp/spWT_Obvs_pred.pdf"),
#       width=14, height=14, unit="cm", dpi=1000)

#########################
### Model performance ###
#########################

print(sp.mod.RF)

## Pespormance from caret
sp.rmse <- getTrainPerf(sp.mod.RF)$TrainRMSE
sp.r2 <- getTrainPerf(sp.mod.RF)$TrainRsquared
sp.mae <- getTrainPerf(sp.mod.RF)$TrainMAE
sp.NSE <- hydroGOF::NSE(sp.plot.data$fv, sp.plot.data$data)
sp.RMSE.calc <- sqrt(mean((sp.plot.data$data - sp.plot.data$fv)^2))

sp.per <- data.frame(RMSE = sp.rmse, R2 = sp.r2, MAE = sp.mae, NSE = sp.NSE)

#print(rf.NSE)

###########################
### Variable importance ###
###########################

sp.varimp <- varImp(sp.mod.RF, scale = T)
sp.varimp.df <- data.frame(sp.varimp$importance)
sp.varimp.df$Variable <- row.names(sp.varimp.df)
colnames(sp.varimp.df) <- c("Importance","Variable")
sp.varimp.df <- sp.varimp.df[,c("Variable","Importance")]
row.names(sp.varimp.df) <- NULL
#sp.varimp.df <- sp.varimp.df[sp.varimp.df$Importance >=1,]

sp.varimp.df <- sp.varimp.df[order(sp.varimp.df$Importance, decreasing = F),]

# Assign plot names to variable codes
sp.varimp.df <- left_join(sp.varimp.df, var.names, by = "Variable") ## See Part 2 for var.names

#    colnames(sp.varimp.df)[colnames(sp.varimp.df) == "Variable"] <- c("Vname")

sp.varimp.df$pImp <- ifelse(sp.varimp.df$Importance < 1.0, round(sp.varimp.df$Importance, digits = 1),
                            round(sp.varimp.df$Importance, digits = 0))

#sp.varimp.df$pImp <- round(sp.varimp.df$Importance * 100, digits = 0)

sp.varimp.df$Vname <- factor(sp.varimp.df$Vname, levels = sp.varimp.df$Vname)

sp.fill <- c("Landcover" = "#edae49","Stream" = "#00798c","Atmosphere" = "grey","Topographic" = "#66a182")

# Option 1: plot LHS
sp.varimp.plot <- ggplot(data = sp.varimp.df, aes(x = Importance, y = Vname, fill = Influence))+
  geom_col() +
  scale_fill_manual(values = sp.fill) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(axis.title.y = element_blank()) +
  geom_text(aes(label = pImp), hjust = -0.1, size = 2)+
  #  xlab("% Increase MSE") +
  #  xlim(0,4) +
  xlab("Variable relative importance") +
  # xlim(0,112) +
  scale_x_continuous(breaks = c(0,25,50,75,100), limits = c(0,105)) +
  scale_y_discrete(position = "left")+
  ggtitle("(c) Site mean stream temperature") +
  theme(legend.position = c(0.82,0.20), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm'), legend.key.size = unit(3,'mm'))+
  theme(axis.text=element_text(size=4), axis.title=element_text(size=5),
        plot.title = element_text(size=7),
        legend.text = element_text(size = 4), legend.title = element_text(size = 5)
  )

sp.varimp.plot

if (options.plot.vi.sp.RHS == TRUE){

# Option 2: plot RHS
sp.varimp.plot <- ggplot(data = sp.varimp.df, aes(y = -Importance, x = Vname, fill = Influence))+
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = sp.fill) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(axis.title.y = element_blank()) +
  geom_text(aes(label = pImp), hjust = 1.1)+
  #  xlab("% Increase MSE") +
  #  xlim(0,4) +
  ylab("Variable relative importance") +
  # xlim(0,112) +
  scale_y_continuous(limits = c(-110,0), breaks = c(0,-25,-50,-75,-100), labels = c(0,25,50,75,100)) +
  scale_x_discrete(position = "top")+
  ggtitle("(c) Site mean stream temperature") +
  theme(legend.position = c(0.2,0.20), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm')) +
  

sp.varimp.plot
}

ggsave(plot = sp.varimp.plot, filename = paste0("Outputs/WaterResearch_Stream_Temp/Fig_3_C_spWT_VarImpPlot.pdf"),
       width=9, height=7, unit="cm", dpi=1000)



################################
### Partial dependence plots ###
################################

label.pdp <- expression("Marginal effect" ~ (degree*C))

if (options.pdp == TRUE){
  
  ### Extract PDP
  library(iml)
  predictor.sp <- Predictor$new(model = sp.mod.RF, data = sp.train, y = sp.response)
  
  sp.varimp.df <- sp.varimp.df[order(sp.varimp.df$Importance, decreasing = T),]
  colnames(sp.varimp.df)[colnames(sp.varimp.df) == "Vname"] <- c("Variable")
  sp.pdp.preds.list <- sp.varimp.df$Variable
  sp.pdp.preds.list <- as.character(sp.pdp.preds.list)
  
  sp.pdp.preds.p.store <- vector('list', length = length(sp.pdp.preds.list))

  for (pp in 1:length(sp.pdp.preds.list)){
    #  for (pp in 1:3){   
    
    sp.pdp.pred <- sp.pdp.preds.list[[pp]]  
    sp.pdp.pred.long <- var.names[var.names$Variable == sp.pdp.pred,2]
    
    ## PDP plotting
    sp.pdp <- FeatureEffect$new(predictor.sp, sp.pdp.pred, method = "pdp", grid.size = 30)
    
    sp.pdp.test <- sp.pdp$results
    
    sp.pdp.plot <- sp.pdp$plot() + 
      #ylim(0,1.1) +
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")+
      theme(panel.background = element_blank())+
      theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
      theme(legend.position = "none")+
      #      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1.10)) +
      #scale_y_continuous(breaks = c(0.25,0.50,0.75), limits = c(0.2,0.85)) +
      scale_y_continuous(breaks = seq(14.0,15.5,by = 0.5), limits = c(14.0,15.5)) +
      ylab(label.pdp) +
      xlab(sp.pdp.pred.long) +
      xlab(pdp.pred.long) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=5),
            plot.title = element_text(size=7),
            legend.text = element_text(size = 5), legend.title = element_text(size = 5)
      )
    
    
    sp.pdp.plot
    
    sp.pdp.preds.p.store[[pp]] <- sp.pdp.plot
    
  }
  
  # Save all pdp's
  sp.pdp.preds.p.all <-  do.call(gridExtra::grid.arrange,sp.pdp.preds.p.store)
  if (options.saving == TRUE) {
  ggsave(plot = sp.pdp.preds.p.all, filename = paste0("Outputs/WaterResearch_Stream_Temp/spWT_PDP_plots.pdf"), 
         width=24, height=24, unit="cm", dpi=1000)
  }
  # select top5 most important
  #sp.pdp.preds.top <- gridExtra::grid.arrange(grobs = sp.pdp.preds.p.store[1:3], ncol = 3)
  #ggsave(plot = sp.pdp.preds.top, filename = paste0("Outputs/WaterResearch_Stream_Temp/spWT_PDP_plots_top5.pdf"), 
  #       width=18, height=4.5, unit="cm", dpi=1000)
  

  
  ######################################
  ### Select PDPs of interest (Figure 6)
  sp.pdp.preds.list <- c("UrbanKm2Med","AF","ForestPcMed")
  sp.pdp.preds.p.store <- vector('list', length = length(sp.pdp.preds.list))
  
  for (pp in 1:length(sp.pdp.preds.list)){
    #  for (pp in 1:3){ 
    sp.pdp.pred <- sp.pdp.preds.list[[pp]]  
    sp.pdp.pred.long <- var.names[var.names$Variable == sp.pdp.pred,4]
    
    ## PDP plotting
    sp.pdp <- FeatureEffect$new(predictor.sp, sp.pdp.pred, method = "pdp", grid.size = 30)
    
    sp.pdp.test <- sp.pdp$results
    
    sp.pdp.plot <- sp.pdp$plot() + 
      #ylim(0,1.1) +
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")+
      theme(panel.background = element_blank())+
      theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
      theme(legend.position = "none")+
      #      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1.10)) +
      #scale_y_continuous(breaks = c(0.25,0.50,0.75), limits = c(0.2,0.85)) +
      scale_y_continuous(breaks = seq(14.0,15.5,by = 0.5), limits = c(14.0,15.5)) +
      theme(axis.title.y = element_blank()) +
      xlab(sp.pdp.pred.long) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=5),
            plot.title = element_text(size=7),
            legend.text = element_text(size = 5), legend.title = element_text(size = 5)
      )
    
    sp.pdp.plot
    
    sp.pdp.preds.p.store[[pp]] <- sp.pdp.plot
    
    
  }
  
  # select top5 most important
  #sp.pdp.preds.top5 <- gridExtra::grid.arrange(grobs = sp.pdp.preds.p.store[], ncol = 2, top = "Marginal effect of important drivers on site mean stream temperature (°C)")
  #sp.pdp.preds.top <- gridExtra::grid.arrange(grobs = sp.pdp.preds.p.store[], ncol = 3, left = "Marginal effect on site mean \nstream temperature (°C)")
  sp.pdp.preds.top <- gridExtra::grid.arrange(grobs = sp.pdp.preds.p.store[], ncol = 3, 
                                              left = textGrob("Site mean \nstream temperature (°C)", gp=gpar(fontsize=5), rot = 90 ) )
  
  #ggsave(plot = sp.pdp.preds.top, filename = paste0("Outputs/WaterResearch_Stream_Temp/Fig_3_B_spWT_PDP_plots_top.pdf"), 
  #       width=13.5, height=4.5, unit="cm", dpi=1000)
  
}

#####################################################################################################################################################################
### PART D: Random Forest model :: spatio-temporal WT ###
#####################################################################################################################################################################

################################
### Temporal variables

data.st <- data.input.temporal.ante

data.st <- data.st[data.st$Site %in% wt.sites,] # Subset to sites
data.st <- dplyr::select(data.st, -contains("60day")) # remove antecedent for 60 days
data.st <- data.st[complete.cases(data.st), ]  # remove NA's, majority are missing time stamp

colnames(data.st)[colnames(data.st) == "Value_Numeric"] <- "data"

################################
### Create time related fields

# Decimal Year 
#data.st$Year <- year(data.st$Date) - 1992 +
#  ( yday(data.st$Date) / yday(as.Date(paste0(year(data.st$Date), "-12-31"), format = "%Y-%m-%d") ) )
#data.st$Date <- NULL ## remove Date as now have 'Year' instead (Date format won't work in RF models)

# fr_yr
#data.st$fr_yr <- data.st$Year - floor(data.st$Year)   # don't create, not needed for RF

# Decimal time
data.st$Time <- as.numeric(substr(data.st$Time,start = 1, stop =2)) + as.numeric(substr(data.st$Time,start = 4, stop =5)) / 60

################################
### Add spatial variables
data.st.sp <- data.sp
data.st.sp$data <- NULL
data.st <- left_join(data.st, data.st.sp, by = "Site")
rm(data.st.sp)

################################
### Add EC and pH

ec$Time <- as.numeric(substr(ec$Time,1,2)) + (as.numeric(substr(ec$Time,4,5))/60)
ec$Time <- round(ec$Time, digits = 2)
ec$Date <- as.Date(ec$Date)

#data.st1 <- data.st[,c("Site","Date","Time","data")]
data.st$Date <- as.Date(as.character(data.st$Date))
data.st$Time <- round(as.numeric(data.st$Time), digits = 2)

data.st <- left_join(data.st, ec, by = c("Site","Date","Time"))

data.st <- data.st[!is.na(data.st$EC),]


ph$Time <- as.numeric(substr(ph$Time,1,2)) + (as.numeric(substr(ph$Time,4,5))/60)
ph$Time <- round(ph$Time, digits = 2)
ph$Date <- as.Date(ph$Date)

data.st$Date <- as.Date(as.character(data.st$Date))
data.st$Time <- round(as.numeric(data.st$Time), digits = 2)

data.st <- left_join(data.st, ph, by = c("Site","Date","Time"))

data.st <- data.st[!is.na(data.st$pH),]

data.st$Year <- NULL
data.st$Date <- NULL


################################
### Rename landcover variables
colnames(data.st)[colnames(data.st) == "Urban"] <- "UrbanPc"
colnames(data.st)[colnames(data.st) == "Ag"] <- "AgPc"
colnames(data.st)[colnames(data.st) == "Forest"] <- "ForestPc"

data.st$UrbanKm2 <- data.st$UrbanPc * data.st$carea_km2
data.st$AgKm2 <- data.st$AgPc * data.st$carea_km2
data.st$ForestKm2 <- data.st$ForestPc * data.st$carea_km2

########################
# Set up response and training data

#data.st.noscale <- data.st

data.st.scale <- data.st
data.st.scale$Site <- NULL
#data.st.scale <- scale(data.st.scale, center = T, scale = T)

st.train <- as.data.frame(data.st.scale) # input predictor variables only
st.train$data <- NULL 
st.train$Site <- NULL

st.response <- data.frame(data = data.st.scale[,c("data")]) # data (y) for training only

#######################################
### Random Forest model using caret ### 
#######################################

# Use tuneRF to identify mtry with lowest OOB
#mtry <- tuneRF(x = st.train, y = st.response$data, ntreeTry = 500, improve = 0.01, stepFactor = 1.5,
#               trace = F, plot = F)

mtry.st <- 28 # manually set from results of above mtry test, 31

# Fit RF model here
myControl <- trainControl(method="repeatedcv", 
                          number=5, 
                          repeats=10)

tunegrid.sp <- expand.grid(.mtry=mtry.st)

st.mod.RF<-train(st.train,
                 st.response$data,
                 method = "rf",
                 trControl=myControl,
                 tuneGrid=tunegrid.sp,
                 ntree= 500,
                 preProc=NULL
                 )


##################
### Model plot ###
##################

### Observed vs modelled
st.plot.data <- st.response
st.plot.data$fv <- predict(st.mod.RF)
st.plot.data$Site <- data.st$Site

st.plot.colours <- data.wt.delta
st.plot.colours$WTtrend <- ifelse(st.plot.colours$CI.low.delta.30y > 0, "positive",
                                  ifelse(st.plot.colours$CI.high.delta.30y < 0, "negative", "neutral")
)
st.plot.colours <- st.plot.colours[,c("Site","WTtrend")]

st.plot.data <- left_join(st.plot.data, st.plot.colours, by = "Site")
st.colours <- c("positive" = "red","negative" = "blue","neutral" = "grey")

# Plot
st.plot.om <- ggplot() +
  geom_point(data = st.plot.data, aes(x = fv, y = data, colour = WTtrend))+
#  scale_colour_manual(values = rf.colours) +
  scale_colour_manual(values = c("purple","purple","purple")) +
  geom_abline(slope = 1,
              intercept = 0, colour = "blue", linetype = "dashed", linewidth = 1) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), plot.title = element_text(size=15)) +
  ggtitle(paste0("(b) Spatio-temporal stream temperature")) +
  xlab("Predicted stream temperature (°C)") +
  ylab("Observed stream temperature (°C)") 
 # xlim(10,18)+
#  ylim(10,18)

#st.plot.om

#ggsave(plot = st.plot.om, filename = paste0("Outputs/WaterResearch_Stream_Temp/stWT_Obvs_pred.pdf"),
#       width=14, height=14, unit="cm", dpi=1000)

# grid arrange all three OM plots

plot.om.all <- gridExtra::grid.arrange(rf.plot.om, st.plot.om, sp.plot.om, ncol = 3)
if (options.saving == TRUE) {
ggsave(plot = plot.om.all, filename = paste0("Outputs/WaterResearch_Stream_Temp/Suppl_Fig_2_Obvs_pred_all.pdf"),
       width=42, height=14, unit="cm", dpi=1000)

ggsave(plot = plot.om.all, filename = paste0("Outputs/WaterResearch_Stream_Temp/Suppl_Fig_2_Obvs_pred_all.jpeg"),
       width=42, height=14, unit="cm", dpi=330)
}
#########################
### Model performance ###
#########################

print(st.mod.RF)

## Pespormance from caret
st.rmse <- getTrainPerf(st.mod.RF)$TrainRMSE
st.r2 <- getTrainPerf(st.mod.RF)$TrainRsquared
st.mae <- getTrainPerf(st.mod.RF)$TrainMAE
st.NSE <- hydroGOF::NSE(st.plot.data$fv, st.plot.data$data)
st.RMSE.calc <- sqrt(mean((st.plot.data$data - st.plot.data$fv)^2))

st.per <- data.frame(RMSE = st.rmse, R2 = st.r2, MAE = st.mae, NSE = st.NSE)


###########################
### Variable importance ###
###########################

st.varimp <- varImp(st.mod.RF, scale = T)
st.varimp.df <- data.frame(st.varimp$importance)
st.varimp.df$Variable <- row.names(st.varimp.df)
colnames(st.varimp.df) <- c("Importance","Variable")
st.varimp.df <- st.varimp.df[,c("Variable","Importance")]
row.names(st.varimp.df) <- NULL

st.varimp.df <- st.varimp.df[order(st.varimp.df$Importance, decreasing = F),]
st.varimp.df <- st.varimp.df[st.varimp.df$Importance >=1,] # cut those variables with importance less than 1

# Assign plot names to variable codes
st.varimp.df <- left_join(st.varimp.df, var.names, by = "Variable") ## See Part 2 for var.names
st.varimp.df$Vname <- ifelse(st.varimp.df$Vname =="Air temperature (change)", "Air temperature (day of sampling)", st.varimp.df$Vname)

#    colnames(st.varimp.df)[colnames(st.varimp.df) == "Variable"] <- c("Vname")

st.varimp.df$pImp <- ifelse(st.varimp.df$Importance < 1.0, round(st.varimp.df$Importance, digits = 1),
                            round(st.varimp.df$Importance, digits = 0))

#st.varimp.df$pImp <- round(st.varimp.df$Importance * 100, digits = 0)

st.varimp.df$Vname <- factor(st.varimp.df$Vname, levels = st.varimp.df$Vname)

st.fill <- c("Landcover" = "#edae49","Stream" = "#00798c","Atmosphere" = "grey","Topographic" = "#66a182")

# Option 1: plot LHS
st.varimp.plot <- ggplot(data = st.varimp.df, aes(x = Importance, y = Vname, fill = Influence))+
  geom_col() +
  scale_fill_manual(values = st.fill) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(axis.title.y = element_blank()) +
  geom_text(aes(label = pImp), hjust = -0.1, size = 2)+
  #  xlab("% Increase MSE") +
  #  xlim(0,4) +
  #xlab("Variable relative importance") +
  theme(axis.title.x = element_blank()) +
  # xlim(0,112) +
  scale_x_continuous(breaks = c(0,25,50,75,100), limits = c(0,105)) +
  scale_y_discrete(position = "left")+
  ggtitle("Relative importance of drivers of changes in water temperature") +
  ggtitle("(b) Spatio-temporal stream temperature") +
  theme(legend.position = c(0.82,0.20), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm'), legend.key.size = unit(3,'mm'))+
  theme(axis.text=element_text(size=4), axis.title=element_text(size=5),
        plot.title = element_text(size=7),
        legend.text = element_text(size = 4), legend.title = element_text(size = 5)
  )

st.varimp.plot

if (options.plot.vi.st.RHS == TRUE){
# Option 2: plot RHS
st.varimp.plot <- ggplot(data = st.varimp.df, aes(y = -Importance, x = Vname, fill = Influence))+
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = st.fill) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
  theme(axis.title.y = element_blank()) +
  geom_text(aes(label = pImp), hjust = 1.1)+
  #  xlab("% Increase MSE") +
  #  xlim(0,4) +
  ylab("Variable relative importance") +
  # xlim(0,112) +
  scale_y_continuous(limits = c(-110,0), breaks = c(0,-25,-50,-75,-100), labels = c(0,25,50,75,100)) +
  scale_x_discrete(position = "top")+
  ggtitle("(b) Drivers of site mean water temperature") +
  theme(legend.position = c(0.2,0.20), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm'))

st.varimp.plot
}

ggsave(plot = st.varimp.plot, filename = paste0("Outputs/WaterResearch_Stream_Temp/Fig_3_B_stWT_VarImpPlot.pdf"),
       width=9, height=7, unit="cm", dpi=1000)

########################
### Combine VarImp plots 

fig2 <- cowplot::ggdraw(plot = blank) +
  draw_plot(plot = rf.varimp.plot,
            x = 0.05, y = 0.66, width = 0.95, height = 0.33) +
  draw_plot(plot = st.varimp.plot,
            x = 0, y = 0.33, width = 1.0, height = 0.33) +
  draw_plot(plot = sp.varimp.plot,
            x = 0.05, y = 0, width = 0.95, height = 0.33) 
  #draw_plot_label(label = "A",x = 0.02, y = 0.98, hjust = -0.5, vjust = 1.5, size = 10) +
  #draw_plot_label(label = "B",x = 0.73, y = 0.98, hjust = -0.5, vjust = 1.5, size = 10)

if (options.saving == TRUE) {
ggsave(plot = fig2, filename = paste0("Outputs/WaterResearch_Stream_Temp/Figure_3_VarImpPlot.pdf"),
       width=9, height=21, unit="cm", dpi=1000)

ggsave(plot = fig2, filename = paste0("Outputs/WaterResearch_Stream_Temp/Figure_3_VarImpPlot.jpeg"),
       width=9, height=21, unit="cm", dpi=330)
}

################################
### Partial dependence plots ###
################################

if (options.pdp == TRUE){
  
  ### Extract PDP
  library(iml)
  predictor.st <- Predictor$new(model = st.mod.RF, data = st.train, y = st.response)
  
  st.varimp.df <- st.varimp.df[order(st.varimp.df$Importance, decreasing = T),]
  st.varimp.df <- st.varimp.df[st.varimp.df$Importance >=1,] # cut to those with varImp >=1
  colnames(st.varimp.df)[colnames(st.varimp.df) == "Vname"] <- c("Variable")
  st.pdp.preds.list <- st.varimp.df$Variable
  st.pdp.preds.list <- as.character(st.pdp.preds.list)
  
  st.pdp.preds.p.store <- vector('list', length = length(st.pdp.preds.list))
  
  for (pp in 1:length(st.pdp.preds.list)){
    #  for (pp in 1:3){   
    
    st.pdp.pred <- st.pdp.preds.list[[pp]]  
    st.pdp.pred.long <- var.names[var.names$Variable == st.pdp.pred,2]
    
    ## PDP plotting
    st.pdp <- FeatureEffect$new(predictor.st, st.pdp.pred, method = "pdp", grid.size = 30)
    
    st.pdp.test <- st.pdp$results
    
    st.pdp.plot <- st.pdp$plot() + 
      #ylim(0,1.1) +
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")+
      theme(panel.background = element_blank())+
      theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
      theme(legend.position = "none")+
      #      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1.10)) +
      #scale_y_continuous(breaks = c(0.25,0.50,0.75), limits = c(0.2,0.85)) +
      scale_y_continuous(breaks = seq(13.5,16.5,by = 0.5), limits = c(13.5,16.5)) +
      ylab(label.pdp) +
      xlab(st.pdp.pred.long) 
    
    st.pdp.plot
    
    st.pdp.preds.p.store[[pp]] <- st.pdp.plot
    
  }
  
  # Save all pdp's
  st.pdp.preds.p.all <-  do.call(gridExtra::grid.arrange,st.pdp.preds.p.store)
  if (options.saving == TRUE) {
  ggsave(plot = st.pdp.preds.p.all, filename = paste0("Outputs/WaterResearch_Stream_Temp/stWT_PDP_plots.pdf"), 
         width=24, height=24, unit="cm", dpi=1000)
  
  # select top5 most important
  st.pdp.preds.top5 <- gridExtra::grid.arrange(grobs = st.pdp.preds.p.store[1:5], ncol = 5)
  ggsave(plot = st.pdp.preds.top5, filename = paste0("Outputs/WaterResearch_Stream_Temp/stWT_PDP_plots_top5.pdf"), 
         width=35, height=7, unit="cm", dpi=1000)
  }
  
  ######################################
  ### Select PDPs of interest (Figure 6)
  st.pdp.preds.list <- c("tmax.14day","tmax.7day","srad.30day")
  st.pdp.preds.p.store <- vector('list', length = length(st.pdp.preds.list))

  for (pp in 1:length(st.pdp.preds.list)){
    #  for (pp in 1:3){   
    st.pdp.pred <- st.pdp.preds.list[[pp]]  
    st.pdp.pred.long <- var.names[var.names$Variable == st.pdp.pred,4]
    
    ## PDP plotting
    st.pdp <- FeatureEffect$new(predictor.st, st.pdp.pred, method = "pdp", grid.size = 30)
    
    st.pdp.test <- st.pdp$results
    
    st.pdp.plot <- st.pdp$plot() + 
      #ylim(0,1.1) +
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")+
      theme(panel.background = element_blank())+
      theme(panel.border=element_rect(fill=NA, colour="black", size=1)) +
      theme(legend.position = "none")+
      #      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limits = c(0,1.10)) +
      #scale_y_continuous(breaks = c(0.25,0.50,0.75), limits = c(0.2,0.85)) +
      scale_y_continuous(breaks = seq(13.5,16.5,by = 0.5), limits = c(13.5,16.5)) +
      theme(axis.title.y = element_blank()) +
      xlab(st.pdp.pred.long) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=5),
            plot.title = element_text(size=7),
            legend.text = element_text(size = 5), legend.title = element_text(size = 5)
      )
    
    st.pdp.plot
    
    st.pdp.preds.p.store[[pp]] <- st.pdp.plot
    
    
  }
  
  # select top5 most important
  st.pdp.preds.top <- gridExtra::grid.arrange(grobs = st.pdp.preds.p.store[], ncol = 3, 
                                              left = textGrob("Spatio-temporal \nstream temperature (°C)", gp=gpar(fontsize=5), rot = 90 ) )

  #ggsave(plot = st.pdp.preds.top, filename = paste0("Outputs/WaterResearch_Stream_Temp/Fig_3_C_stWT_PDP_plots_top.pdf"), 
  #       width=13.5, height=4.5, unit="cm", dpi=1000)
  

}

###############################
### Figure 3 = All PDP Plots

fig3 <- cowplot::ggdraw(plot = blank) +
  draw_plot(plot = rf.pdp.preds.top,
            x = 0, y = 0.66, width = 1, height = 0.33) +
  draw_plot(plot = st.pdp.preds.top,
            x = 0, y = 0.33, width = 0.75, height = 0.33) +
  draw_plot(plot = sp.pdp.preds.top,
            x = 0, y = 0, width = 0.75, height = 0.33) +
  draw_plot_label(label = "(a)",fontface = "bold",x = 0.01, y = 0.98, hjust = -0.5, vjust = 1.0, size = 9) +
  draw_plot_label(label = "(b)",fontface = "bold",x = 0.01, y = 0.64, hjust = -0.5, vjust = 1.0, size = 9) +
  draw_plot_label(label = "(c)",fontface = "bold",x = 0.01, y = 0.31, hjust = -0.5, vjust = 1.0, size = 9)

if (options.saving == TRUE) {
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_4_pdp.pdf", fig3, width=18, height=15, unit="cm", dpi=1000)
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_4_pdp.jpeg", fig3, width=18, height=15, unit="cm", dpi=330)
}

#####################################################################################################################################################################
### PART E: Methods figures ###
#####################################################################################################################################################################

#########################
### Map of the region ###
#########################

GIS.wt.sitemean <- left_join(GIS.wt.sitenum, wt.site.mean, by = "Site")

site.mean.label <- bquote(paste("Site mean water \ntemperature °C"))

fig.melb <- ggplot() +
  geom_sf(data = vic, fill = NA, colour = "black") +
  # geom_sf(data = GIS.catch, fill = "white") +
  geom_sf(data = GIS.storages[GIS.storages$ASSET_NAME != "THOMSON RESERVOIR",], fill = "grey70", colour = NA ,linewidth = 0.1,alpha = 0.4) +
  geom_sf(data = GIS.streams, colour = "grey70", alpha = 0.4) +
  geom_sf(data = GIS.wt.catch, fill = "NA", colour = "red") +
  #  geom_sf(data = GIS.wt, colour = "black") +
  geom_sf(data = GIS.wt.sitemean, shape = 21, size = 2.5, colour = "black",aes(fill = WT.site.mean)) +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, na.value = "grey60", limits = c(10,18)) +
  #guides(fill = guide_legend(order = 1, reverse = T)) +
  labs(fill = site.mean.label)+
  #  theme(legend.position = "right", legend.spacing.y = unit(1, 'mm'))+
  theme(legend.position = c(0.90,0.135), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm'))+
  theme(legend.title = element_text(size = 6), 
        legend.text  = element_text(size = 6),
        legend.key.size = unit(0.55, "lines"),
        legend.text.align = 0)+
 # theme(legend.position = "none")+
  #  scale_colour_distiller(palette = "RdBu", direction = -1, na.value = "grey60") +
  #  scale_colour_scico(palette = 'batlow') +
  ggspatial::annotation_scale(location = "bl", bar_cols = c("grey60", "white"), width_hint = 0.10) +
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.5, "cm"), pad_y = unit(0.75, "cm"),
                                    style = north_arrow_fancy_orienteering()) +
  #coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.3), ylim = c(-37.1,-38.6), expand = F)
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(144, 146.44), ylim = c(-37.1,-38.77), expand = F)
  

## Put together
fig.site <- cowplot::ggdraw(plot = fig.melb) +
  draw_plot(plot = fig.oz,
            x = 0.747, y = 0.7235,
            width = 0.26, height = 0.26)

#fig.all

#ggsave(plot = fig.site, filename = "Outputs/WaterResearch_Stream_Temp/Map_region.pdf", width = 15, height = 10, unit = "cm", dpi = 1000)



###############################
### Landuse over the region ###
###############################

library(dplyr)
library(sf)
library(terra)
library(sp)
library(raster)

setwd("~/Onedrive/Data/Land_use_GIS/Vic_land_cover_timeseries/")
filename <- "VIC_LANDCOVER_TS/VIC_LANDCOVER_TS_2015_19.tif"
lc90_95 <- brick(filename, varname = "lc_class")

fig2.catch <- st_transform(GIS.catch, st_crs(lc90_95))
fig2.catch <- st_zm(fig2.catch, drop = TRUE)
fig2.catch.crop <- lc90_95 %>% terra::mask(fig2.catch) %>% terra::crop(fig2.catch)

fig2.catch.crop_15_19 <- fig2.catch.crop
fig2.urban.for.ag <- fig2.catch.crop_15_19

## Classes
# Urban = 2
# Ag / Grassland = 4,5, 8, 10, 13
# Forest = 12

# Project to same CRS as for geom_sf plots above
#crs(fig2.urban.for.ag) <- st_crs("EPSG:28355")$proj4string

crs(fig2.urban.for.ag) <- crs("+init=EPSG:4326")
GIS.catch <- st_transform(GIS.catch, crs = 4326)

cat(wkt(fig2.urban.for.ag))

fig2.urban.for.ag.10 <- raster::aggregate(fig2.urban.for.ag, fact=1, FUN=mean)

fig2.urban.for.ag.df <- as.data.frame(fig2.urban.for.ag.10, xy=TRUE)
fig2.urban.for.ag.df2 <- fig2.urban.for.ag.df[fig2.urban.for.ag.df$VIC_LANDCOVER_TS_2015_19 %in% c(2,4,5,8,10,12,13),]

fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 <- ifelse(fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 == 5, 4, fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19)
fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 <- ifelse(fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 == 8, 4, fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19)
fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 <- ifelse(fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 == 10, 4, fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19)
fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 <- ifelse(fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 == 13, 4, fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19)

fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19 <- as.factor(fig2.urban.for.ag.df2$VIC_LANDCOVER_TS_2015_19)

fig7C.box <- data.frame(x1 = min(fig2.urban.for.ag.df2$x), x2 = max(fig2.urban.for.ag.df2$x) + 0.3 * (max(fig2.urban.for.ag.df2$x) - min(fig2.urban.for.ag.df2$x)),
               y1 = min(fig2.urban.for.ag.df2$y), y2 = max(fig2.urban.for.ag.df2$y))



fig7C <- ggplot() +
  # geom_sf(data = GIS.catch, fill = NA)+
  geom_rect(data=fig7C.box, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill = "white",color="white")+
  geom_raster(data = fig2.urban.for.ag.df2, aes(x = x, y = y, fill = `VIC_LANDCOVER_TS_2015_19`))+
  scale_fill_manual(values = c("2" = "grey60", "12" = "forestgreen", "4"="steelblue"),
                    labels = c("2" = "Urban", "12" = "Forest", "4"="Agriculture")) +
  labs(fill = "Landcover")+
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(legend.position = c(0.1,0.1), legend.direction = "hortizontal", 
        legend.spacing.y = unit(1, 'mm'))+
  theme(legend.title = element_text(size = 5, angle = 90), 
        legend.text  = element_text(size = 5, angle = 90),
        legend.key.size = unit(0.7, "lines"),
        legend.text.align = 0,
        legend.title.align=0.5)+
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank()) +
  guides(fill = guide_legend(title.position = "left")) 

  
rasterOptions()
rasterOptions(memfrac=.3)

fig7C_legend <- get_legend(fig7C)

fig7C <- fig7C+theme(legend.position = "none")







#setwd("~/Onedrive/R/")
#ggsave(plot = fig2, filename = "Outputs/WaterResearch_Stream_Temp/Map_info_landuse.pdf", width = 15, height = 13, unit = "cm", dpi = 1000)



#######################################
### Air temperature over the region ###
#######################################


library(ncdf4)
library(ncdf4.helpers)
library(PCICt)

# tmax
setwd("~/Onedrive/R/")
awap <- nc_open("Data/Bureau/agcd_v1_tmax_mean_r005_daily_1970-2020_Melb.nc")
a.var <- "tmax"

print(awap)

lat <- ncvar_get(awap, "lat")
lon <- ncvar_get(awap, "lon")
time <- ncvar_get(awap, "time")
tunits <- ncatt_get(awap, "time", "units") 

at_array <- ncvar_get(awap, "tmax")
fillvalue <- ncatt_get(awap, "tmax", "_FillValue")
at_array[at_array==fillvalue$value] <- NA


#time_obs <- as.POSIXct(time, origin = "1850-01-01", tz="GMT")
time_obs <- as.Date(time, origin = "1850-01-01", tz="GMT")
time_obs <- time_obs[time_obs >= as.Date("1992-01-01", format("%Y-%m-%d"))]
time_obs <- time_obs[time_obs < as.Date("2020-01-01", format("%Y-%m-%d"))]


time_obs2 <- as.Date(time, origin = "1850-01-01", tz="GMT")
time_obs2 <- data.frame(time_obs2)
time_obs3 <- data.frame(cbind(time_obs2, time))
time_obs3$slice <- seq(1, length(time_obs3$time),by =1)

# "1992-01-01" >> 51864.38 >> 8036
# "2020-01-01" >> 62091.38 >> 18263

range(time_obs)

at.array.short <- at_array[,,8036:18262]

# now put array data into dataframe
lonlattime <- as.matrix(expand.grid(lon,lat,time_obs))
at_vec_long <- as.vector(at.array.short)
at_obs <- data.frame(cbind(lonlattime, at_vec_long))

# THIS WAS HARD TO MAKE, DONT OVERWRITE

# Now make summary
at <- at_obs
colnames(at) <- c("lon","lat","Date","tmax")
at$tmax <- as.numeric(at$tmax)
at.summary <- at %>% group_by(lon, lat) %>% summarise(meanTmax = mean(tmax, na.rm = TRUE))

at.sum <- at.summary

# now make into grid
at.sum <- st_as_sf(at.sum, coords = c("lon", "lat"), 
                   crs = 4326)

at.poly <- at.sum %>%
  st_make_grid(cellsize = 0.05) %>% 
  st_as_sf() %>% 
  st_join(at.sum)

#ggplot() +
#  geom_sf(data = at.poly, aes(fill = meanTmax))
sf::sf_use_s2(FALSE)
GIS.catch.all <- st_union(GIS.catch)
at.poly <- st_transform(at.poly, st_crs(GIS.catch.all))

at.poly.within <- at.poly[GIS.catch,] # subset to only within Melbourne area

# make mask
melb_for_mask <- GIS.catch %>%
  st_geometry() %>%
  st_cast('POLYGON') %>%
  st_union()

melb_bbox <- st_bbox(at.poly.within)

xrange <- melb_bbox$xmax - melb_bbox$xmin
yrange <- melb_bbox$ymax - melb_bbox$ymin

melb_bbox[1] <- melb_bbox[1] - (0.005 * xrange)
melb_bbox[3] <- melb_bbox[3] + (0.3 * xrange)
melb_bbox[2] <- melb_bbox[2] - (0.005 * yrange)
melb_bbox[4] <- melb_bbox[4] + (0.005 * yrange)

melb_bbox <- melb_bbox %>%
  st_as_sfc() %>% 
  st_as_sf()

melb_mask <- st_difference(melb_bbox, melb_for_mask)

air.temp.label <- bquote(paste("Mean maximum daily \nair temperature °C"))


fig7B <- ggplot() +
  geom_sf(data = at.poly.within, aes(fill = meanTmax, colour = meanTmax))+
  #  scale_fill_scico(palette = 'bam', direction = -1) +
  scale_fill_distiller(palette = "RdBu", direction = -1, na.value = "grey60") +
  scale_colour_distiller(palette = "RdBu", direction = -1, na.value = "grey60") +
  geom_sf(data = GIS.catch.all, fill = NA)+
  geom_sf(data = melb_mask, fill = "white", colour = "white") +
  theme(panel.background = element_blank())+
  theme(panel.border=element_rect(fill=NA, colour="black", linewidth=1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(plot.background = element_blank()) +
  #theme(legend.position = "none") +
  labs(fill = air.temp.label, colour = air.temp.label)+
  theme(legend.position = c(1.2,0.5), legend.direction = "vertical", 
        legend.spacing.y = unit(1, 'mm'),)+
  theme(legend.title = element_text(size = 5, angle = 90), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.7, "lines"),
        legend.text.align = 1,
        legend.title.align=0.5)+
  guides(colour = guide_colorbar(title.position = "left"), fill = guide_colorbar(title.position = "left")) +
  #guide_legend(title.position = "left") +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank()) 

fig7B_legend <- get_legend(fig7B)

fig7B <- fig7B+theme(legend.position = "none")

#setwd("~/Onedrive/R/")
#ggsave(plot = fig3, filename = "Outputs/WaterResearch_Stream_Temp/Map_info_airtemp.pdf", width = 15, height = 13, unit = "cm", dpi = 1000)

##########################################
### Combine all methods plots into one ###
##########################################

fig7 <- cowplot::ggdraw(plot = blank) +
  draw_plot(plot = fig.site,
            x = 0, y = 0, width = 0.65, height = 1.0) +
  draw_plot(plot = fig7B,
            x = 0.64, y = 0.512, width = 0.36, height = 0.5) +
  draw_plot(plot = fig7B_legend,
            x = 0.83, y = 0.70, width = 0.1, height = 0.1) +
  draw_plot(plot = fig7C,
            x = 0.64, y = 0.055, width = 0.36, height = 0.48) +
  draw_plot(plot = fig7C_legend,
            x = 0.94, y = 0.27, width = 0.1, height = 0.1) +
  draw_plot_label(label = "(a)",fontface = "bold",x = 0.08, y = 0.95, hjust = -0.5, vjust = 1.0, size = 9) +
  draw_plot_label(label = "(b)",fontface = "bold",x = 0.66, y = 0.95, hjust = -0.5, vjust = 1.0, size = 9) +
  draw_plot_label(label = "(c)",fontface = "bold",x = 0.66, y = 0.50, hjust = -0.5, vjust = 1.0, size = 9)

if (options.saving == TRUE) {
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_1_Methods.pdf", fig7, width=18, height=10, unit="cm", dpi=1000)
ggsave(file = "Outputs/WaterResearch_Stream_Temp/Figure_1_Methods.jpeg", fig7, width=18, height=10, unit="cm", dpi=330)
}

#########################
### Save .Rdata image ###
#########################

save.image(file="Outputs/WaterResearch_Stream_Temp/WaterResearch_20230905.RData") 

