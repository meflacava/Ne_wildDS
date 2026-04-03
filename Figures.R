#Figures

##Colors & symbols
# Nc*        = "black",          pch=15 (square)
# Combo Ne   = "darkolivegreen", pch=16 (circle)
# LD Ne      = "darkorange",     pch=17 (triangle)
#quick visualize them all together
#plot(c(1,1,1),col=c("black","darkolivegreen","darkorange"),pch=c(15,16,17))


#### Compile Ne values ####
ld <- read.csv("ByYear/LD_BY1995-2020_minMAF0.02_Dp5-50_acrossChr-RW.csv") #Robin calc CIs (Aug 25, 2025)
combo1 <- read.csv("ComboNe/ComboNe_BY1995-1999.csv")
combo2 <- read.csv("ComboNe/ComboNe_BY2008-2020.csv")
tne1 <- read.csv("ByGroup/tNe_forComboNe_BY1995-1999.csv")
tne2 <- read.csv("ByGroup/tNe_forComboNe_BY2008-2020.csv")
nc <- read.csv("abundance/DeltaSmeltAbundIndices_Jan2025.csv")


## Start df with abundance index estimates
df <- data.frame(nc[nc$CalendarYear %in% 1995:2020 & nc$Type=="MB_winter",
                    c("CalendarYear","Nabund","NabundSE","Lower95CI","Upper95CI")])
names(df) <- c("yr","Nabund","SE_Nabund","lowerCI_Nabund","upperCI_Nabund")

## Add n genetic samples & BestNe values
df$n_genetic_samples <- NA
df$BestNe <- NA
df$lowerCI_BestNe <- NA
df$upperCI_BestNe <- NA
df$type_BestNe <- NA
df$BestNe_over_Nc <- NA
df$nloci_LDNe <- NA
df$LDNe <- NA
df$lowerCI_LDNe <- NA
df$upperCI_LDNe <- NA
df$nloci_tNe <- NA
df$tNe_MaxTempAdjusted <- NA
df$lowerCI_tNe_MaxTempAdjusted <- NA
df$upperCI_tNe_MaxTempAdjusted <- NA
df$LD_weight <- NA
df$temporal_weight <- NA

for (i in 1:nrow(df)){
    
  if (df$yr[i] %in% ld$yr){
    df$n_genetic_samples[i] <- round(ld$n.samples[ld$yr==df$yr[i]])
    df$nloci_LDNe[i] <- ld$n.loci[ld$yr==df$yr[i]]
    df$LDNe[i] <- ld$Ne[ld$yr==df$yr[i]]
    df$lowerCI_LDNe[i] <- ld$lowNe[ld$yr==df$yr[i]] #Waples et al. 2022 method to calc CIs
    df$upperCI_LDNe[i] <- ld$hiNe[ld$yr==df$yr[i]] #Waples et al. 2022 method to calc CIs
  }
  if (df$yr[i] %in% combo1$yr){
    df$BestNe[i] <- round(combo1$ComboNe[combo1$yr==df$yr[i]])
    df$lowerCI_BestNe[i] <- round(combo1$lowerCI[combo1$yr==df$yr[i]])
    df$upperCI_BestNe[i] <- round(combo1$upperCI[combo1$yr==df$yr[i]])
    df$type_BestNe[i] <- "combo_Ne"
    
    df$nloci_tNe[i] <- tne1$L[tne1$yr==df$yr[i]]
    df$tNe_MaxTempAdjusted[i] <- round(tne1$Ne[tne1$yr==df$yr[i]])
    df$lowerCI_tNe_MaxTempAdjusted[i] <- round(tne1$lowerCI[tne1$yr==df$yr[i]])
    df$upperCI_tNe_MaxTempAdjusted[i] <- round(tne1$upperCI[tne1$yr==df$yr[i]])
    
    df$LD_weight[i] <- round(combo1$LD_weight[combo1$yr==df$yr[i]],2)
    df$temporal_weight[i] <- round(combo1$temporal_weight[combo1$yr==df$yr[i]],2)
  }
  else if (df$yr[i] %in% combo2$yr){
    df$BestNe[i] <- round(combo2$ComboNe[combo2$yr==df$yr[i]])
    df$lowerCI_BestNe[i] <- round(combo2$lowerCI[combo2$yr==df$yr[i]])
    df$upperCI_BestNe[i] <- round(combo2$upperCI[combo2$yr==df$yr[i]])
    df$type_BestNe[i] <- "combo_Ne"
    
    df$nloci_tNe[i] <- tne2$L[tne2$yr==df$yr[i]]
    df$tNe_MaxTempAdjusted[i] <- round(tne2$Ne[tne2$yr==df$yr[i]])
    df$lowerCI_tNe_MaxTempAdjusted[i] <- round(tne2$lowerCI[tne2$yr==df$yr[i]])
    df$upperCI_tNe_MaxTempAdjusted[i] <- round(tne2$upperCI[tne2$yr==df$yr[i]])
    
    df$LD_weight[i] <- round(combo2$LD_weight[combo2$yr==df$yr[i]],2)
    df$temporal_weight[i] <- round(combo2$temporal_weight[combo2$yr==df$yr[i]],2)
  }
  else if (df$yr[i] %in% ld$yr){
    df$BestNe[i] <- round(ld$Ne[ld$yr==df$yr[i]])
    df$lowerCI_BestNe[i] <- round(ld$lowNe[ld$yr==df$yr[i]]) #Waples et al. 2022 method to calc CIs
    df$upperCI_BestNe[i] <- round(ld$hiNe[ld$yr==df$yr[i]]) #Waples et al. 2022 method to calc CIs
    df$type_BestNe[i] <- "LD_Ne"
  }
  #calc Ne/Nc
  if (!is.na(df$BestNe[i]) & !is.infinite(df$BestNe[i]) & df$BestNe[i]/df$Nabund[i] > 0){
    df$BestNe_over_Nc[i] <- round(df$BestNe[i]/df$Nabund[i],4)
  }
}
df


## Set neg values to Inf
df[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)], function(x) {
  x[x<0] <- Inf
  return(x)
})


## Write compiled results
options(scipen = 999)
#write.csv(df,"Ne_Nc_DeltaSmelt_2025nov14.csv",row.names=F,quote=F)
options(scipen = 0)









#### Fig 2. Two-panel plot of Ne over time and Nc* over time ####

## Import data
n <- read.csv("Ne_Nc_DeltaSmelt_2025nov14.csv")
max(n$BestNe[n$BestNe!=Inf],na.rm=T)
max(n$upperCI_BestNe[n$upperCI_BestNe!=Inf],na.rm=T)
#set Inf values to 5000
n[sapply(n, is.numeric)] <- lapply(n[sapply(n, is.numeric)], function(x) {
  x[is.infinite(x)] <- 5000
  x[x>5000] <- 4999 #to distinguish for plotting
  return(x)
})


pdf("Figures/Fig2_Ne_Nc_overtime.pdf",
    width=7,height=7)
par(mfrow=c(2,1),mar=c(3.1,4.1,2.1,2.1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)

## Ne over time
plot(BestNe~yr,data=n,col="white",log="y",ylim=c(100,5000),yaxt="n",
     xlab="",ylab="") #set up plotting area
#ylab
par(las = 2)  
mtext(expression(hat(N)[e]), side = 2, line = 2.8) 
par(las = 0) 

#add dashed line for Inf
abline(h=5000,lty=2)

#combo Ne values
points(BestNe~yr,data=n[n$type_BestNe=="combo_Ne",],pch=16,col="darkolivegreen",cex=1.1) #combo Ne values
arrows(n$yr[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe<5000],n$lowerCI_BestNe[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe<5000],
       n$yr[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe<5000],n$upperCI_BestNe[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe<5000],
       length=0,col="darkolivegreen")
arrows(n$yr[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe==5000],n$lowerCI_BestNe[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe==5000],
       n$yr[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe==5000],n$upperCI_BestNe[n$type_BestNe=="combo_Ne" & n$upperCI_BestNe==5000],
       length=0,col="darkolivegreen",lty=5)

#LD Ne
points(BestNe~yr,data=n[n$type_BestNe=="LD_Ne",],pch=17,col="darkorange",cex=0.9) #other values
arrows(n$yr[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe<5000],n$lowerCI_BestNe[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe<5000],
       n$yr[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe<5000],n$upperCI_BestNe[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe<5000],
       angle=90,length=0,col="darkorange")
arrows(n$yr[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe==5000],n$lowerCI_BestNe[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe==5000],
       n$yr[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe==5000],n$upperCI_BestNe[n$type_BestNe=="LD_Ne" & n$upperCI_BestNe==5000],
       angle=90,length=0,col="darkorange",lty=5)
legend(1999,300,legend=c("Combo","LD"),cex=0.8,pch=c(16,17),pt.cex=c(1,0.9),
       col=c("darkolivegreen","darkorange"))

#x-axis ticks
axis(side=1,at=setdiff(1995:2020, seq(1995, 2020, by = 5)),labels=F,tick=T,tcl = -0.3)
#log y-axis
major_ticks <- 10^(2:3)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 2:3)),las=1)
minor_ticks <- unlist(sapply(2:4, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)
axis(2,at=5000,label=expression(infinity),tick=T,las=1,cex.axis=1.5)

## Add panel label
mtext("(a)",font=2,side=2,line=3,adj=1,las=1,padj=-10)

## Nc* over time
par(mar=c(4.5,4.1,1.1,2.1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)
n <- read.csv("Ne_Nc_DeltaSmelt_2025nov14.csv")
plot(Nabund~yr,data=n,log="y",yaxt="n",ylim=c(min(n$lowerCI_Nabund),max(n$upperCI_Nabund)),
     xlab="Year",ylab="",
     pch=15)
arrows(n$yr,n$lowerCI_Nabund,
       n$yr,n$upperCI_Nabund,
       angle=90,length=0)
#ylab
par(las = 2)  
mtext(expression(N[c] * "*"), side = 2, line = 2.8) 
par(las = 0) 
#x-axis ticks
axis(side=1,at=setdiff(1995:2020, seq(1995, 2020, by = 5)),labels=F,tick=T,tcl = -0.3)
#log y-axis (for Nabund & CI range 1063-6,227,231)
major_ticks <- 10^(3:7)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 3:7)),las=1)
minor_ticks <- unlist(sapply(3:7, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)

## Add panel label
mtext("(b)",font=2,side=2,line=3,adj=1,las=1,padj=-10)

dev.off()


##### lin reg Ne ~ yr ####

## Import data
n <- read.csv("Ne_Nc_DeltaSmelt_2025nov14.csv")
n_clean <- n[!is.na(n$BestNe) & is.finite(n$BestNe), ]

cor(log(n_clean$BestNe),n_clean$yr)
#-0.2617383
cor.test(n_clean$BestNe,n_clean$yr)
# Pearson's product-moment correlation
# data:  n_clean$BestNe and n_clean$yr
# t = -0.57164, df = 18, p-value = 0.5746
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5439111  0.3283958
# sample estimates:
#        cor 
# -0.1335298 


#### Fig 3: linear regression Ne ~ Nc ####

## Import data
n <- read.csv("Ne_Nc_DeltaSmelt_2025nov14.csv")
n_clean <- n[!is.na(n$BestNe) & is.finite(n$BestNe), ]

## Pearson's correlation coefficient
cor(log(n_clean$BestNe),log(n_clean$Nabund))
#0.5840851
cor.test(log(n_clean$BestNe), log(n_clean$Nabund))
# Pearson's product-moment correlation
# 
# data:  log(n_clean$BestNe) and log(n_clean$Nabund)
# t = 3.053, df = 18, p-value = 0.006848
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1909083 0.8157572
# sample estimates:
#       cor 
# 0.5840851 


## Linear regression

pdf("NeDeclinePub/Figures/Fig3_Ne_Nc.pdf",
    width=6,height=4)
par(mar=c(4.5,6.5,0.5,1.1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)

#Make the base plot with log scales
plot(BestNe ~ Nabund, data=n_clean, log="xy", pch = 16, col="gray60",ylim=c(100,3000),
     xaxt="n",yaxt="n",xlab=expression(N[c] * "*"),ylab="")
#ylab
par(las = 2)
mtext(expression(hat(N)[e]), side = 2, line = 2.8)
par(las = 0)
#log x-axis
major_ticks <- 10^(3:7)
axis(1, at = major_ticks, labels = parse(text = paste0("10^", 3:7)),las=1)
minor_ticks <- unlist(sapply(3:7, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(1, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)
#log y-axis
major_ticks <- 10^(2:4)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 2:4)),las=1)
minor_ticks <- unlist(sapply(2:4, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)

#Get predicted line on the original scale & add regression line to plot
newx <- seq(min(n_clean$Nabund), max(n_clean$Nabund), length.out=100)
pred <- predict(lm(log(BestNe) ~ log(Nabund), data=n_clean), newdata=data.frame(Nabund=newx))
lines(newx, exp(pred))

dev.off()


#Model summary
summary(lm(log(BestNe) ~ log(Nabund), data=n_clean))
# Call:
#   lm(formula = log(BestNe) ~ log(Nabund), data = n_clean)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.1689 -0.2278  0.2428  0.4122  0.8680 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  3.81338    1.01109   3.772  0.00140 **
#   log(Nabund)  0.24692    0.08088   3.053  0.00685 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5992 on 18 degrees of freedom
# Multiple R-squared:  0.3412,	Adjusted R-squared:  0.3046 
# F-statistic: 9.321 on 1 and 18 DF,  p-value: 0.006848





#### Fig 4. Ne/Nc by year ####

#use clean data
n <- read.csv("Ne_Nc_DeltaSmelt_2025nov14.csv")
#clean n
n_clean <- n[!is.na(n$BestNe) & !is.infinite(n$BestNe),]


pdf("NeDeclinePub/Figures/Fig4_NeOverNc_OverTime_BestNe.pdf",
    width=6,height=4)
par(mar=c(4.5,6.5,0.5,1.1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)

## Plot log(Ne/Nc ratio) by year with regression
plot(BestNe_over_Nc~yr,data=n_clean,log="y",ylim=c(0.0001,1),yaxt="n",
     xlab="Year",ylab="",col="white",
     pch=16,cex=1.2,cex.axis=0.9)
points(BestNe_over_Nc~yr,data=n[n$type_BestNe=="combo_Ne",],pch=16,col="darkolivegreen",cex=1.1) #combo Ne values
points(BestNe_over_Nc~yr,data=n[n$type_BestNe=="LD_Ne",],pch=17,col="darkorange",cex=0.9) #other values

#ylab
par(las = 2)  
mtext(expression(hat(N)[e] * " / " * N[c] * "*"), side = 2, line = 3.2) 
par(las = 0) 
#x-axis
axis(side=1,at=setdiff(1995:2020, seq(1995, 2020, by = 5)),labels=F,tick=T,tcl = -0.3)
#log y-axis
major_ticks <- 10^seq(-3, 0, by = 1) 
axis(2, at = major_ticks, labels = c("0.001", "0.01", "0.1", "1"),las=1,cex.axis=0.9)
minor_ticks <- unlist(sapply(-4:0, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3)
#regression lines
lines(n_clean$yr[!is.na(n_clean$BestNe_over_Nc)], exp(predict(lm(log(BestNe_over_Nc) ~ yr, 
                                                      data = n_clean[!is.na(n_clean$BestNe_over_Nc), ]))))

#legend
legend(1995,1,legend=c(expression("Combo "*hat(N)[e]),expression("LD "*hat(N)[e])),
       cex=0.8,pch=c(16,17),pt.cex=c(1,0.9),col=c("darkolivegreen","darkorange"))

dev.off()


#summary(lm(log(BestNe_over_Nc) ~ yr, data = n_clean[!is.na(n_clean$BestNe_over_Nc),]))
r <- round(summary(lm(log(BestNe_over_Nc) ~ yr, data = n_clean[!is.na(n_clean$BestNe_over_Nc),]))$r.squared,2)
r #0.7
# Call:
#   lm(formula = log(BestNe_over_Nc) ~ yr, data = n_clean[!is.na(n_clean$BestNe_over_Nc), 
#   ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.82850 -0.36744  0.07098  0.46740  1.42753 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -314.31370   47.11501  -6.671 2.94e-06 ***
#   yr             0.15377    0.02346   6.554 3.70e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8128 on 18 degrees of freedom
# Multiple R-squared:  0.7047,	Adjusted R-squared:  0.6883 
# F-statistic: 42.95 on 1 and 18 DF,  p-value: 3.7e-06







#### X-fold variability in values ####

n <- read.csv("Ne_Nc_DeltaSmelt_2025nov14.csv")

#Nc*
max(n$Nabund)/min(n$Nabund) 
#1220x

#Ne (exclude Inf values)
max(n$BestNe[!is.na(n$BestNe) & !is.infinite(n$BestNe)])/min(n$BestNe[!is.na(n$BestNe) & !is.infinite(n$BestNe)])
#13x

#Ne/Nc*
max(n$BestNe_over_Nc,na.rm=T)/min(n$BestNe_over_Nc,na.rm=T)
#562x

#Combo Ne
n_clean <- n[!is.na(n$BestNe) & !is.infinite(n$BestNe),]
max(n_clean$BestNe[n_clean$type_BestNe=="combo_Ne"])/min(n_clean$BestNe[n_clean$type_BestNe=="combo_Ne"])
#13x

#LD Ne
max(n$LDNe,na.rm=T)/min(n$LDNe,na.rm=T)
#6x

#Max-Temp adjusted temporal Ne
n_clean <- n[!is.na(n$tNe_MaxTempAdjusted) & !is.infinite(n$tNe_MaxTempAdjusted),]
max(n_clean$tNe_MaxTempAdjusted)/min(n_clean$tNe_MaxTempAdjusted)
#146x





#### Compare CIs across methods ####

## Ne values
plot(n$yr[n$type_BestNe=="combo_Ne"],n$BestNe[n$type_BestNe=="combo_Ne"],
     type="l",col="darkolivegreen",
     xlab="year",ylab="Ne",ylim=c(100,5000))
lines(n$yr,n$LDNe,col="darkorange")
lines(n$yr[!is.infinite(n$tNe_MaxTempAdjusted)],n$tNe_MaxTempAdjusted[!is.infinite(n$tNe_MaxTempAdjusted)],col="steelblue4")
#some tNe values are massive, so LD Ne and combo Ne are lower and more consistent

## lower Ne CI values
plot(n$yr[n$type_BestNe=="combo_Ne"],n$lowerCI_BestNe[n$type_BestNe=="combo_Ne"],
     type="l",col="darkolivegreen",
     xlab="year",ylab="Ne",ylim=c(100,5000))
lines(n$yr,n$lowerCI_LDNe,col="darkorange")
lines(n$yr,n$lowerCI_tNe_MaxTempAdjusted,col="steelblue4")
#one massive lower CI for tNe, one massive lower CI for combo Ne

## Size of confidence intervals
tightest_df <- data.frame(
  yr = n$yr,
  BestNe_span = ifelse(n$type_BestNe == "combo_Ne",
                       n$upperCI_BestNe - n$lowerCI_BestNe,
                       NA),
  LDNe_span = n$upperCI_LDNe - n$lowerCI_LDNe,
  tNe_span = n$upperCI_tNe_MaxTempAdjusted - n$lowerCI_tNe_MaxTempAdjusted
)

tightest_df$tightest_CI_method <- apply(tightest_df[, 2:4], 1, function(row) {
  valid <- !is.na(row) & is.finite(row)
  if (!any(valid)) return(NA)
  names(row)[which.min(row[valid])]
})
tightest_df
#BestNe has tightest CI span for non-NA, non-Inf CIs



## abundance index CIs
abund <- data.frame(
  yr = n$yr,
  abs_CI_width = n$upperCI_Nabund - n$lowerCI_Nabund,
  rel_CI_width = (n$upperCI_Nabund - n$lowerCI_Nabund) / n$Nabund
)
plot(abund$yr,abund$abs_CI_width)
plot(abund$yr,abund$rel_CI_width)
#relatively wider CIs for lower Nc estimates than for higher Nc estimates


#### Fig S2. LD Ne over time ####

## Import data
ld <- read.csv("ByYear/LD_BY1995-2020_minMAF0.02_Dp5-50_acrossChr-RW.csv")
#set neg values to Inf
ld[sapply(ld, is.numeric)] <- lapply(ld[sapply(ld, is.numeric)], function(x) {
  x[x<0] <- Inf
  return(x)
})
#check values
min(ld$lowNe)
max(ld$Ne)
max(ld$hiNe[!is.infinite(ld$hiNe)])
#set Inf values to 8000
ld[sapply(ld, is.numeric)] <- lapply(ld[sapply(ld, is.numeric)], function(x) {
  x[is.infinite(x)] <- 8000
  x[x>8000] <- 7999 #to distinguish for plotting
  return(x)
})


## Plot
pdf("NeDeclinePub/Figures/FigS2_LDNe_OverTime.pdf",
    width=6,height=3.5)
par(mar=c(4.5,5.1,0.5,2.1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)

plot(Ne~yr,data=ld,log="y",yaxt="n",ylim=c(100,8000),
     xlab="Year",ylab="",
     pch=17,col="darkorange")
arrows(ld$yr[!is.na(ld$Ne) & ld$hiNe<8000],ld$lowNe[!is.na(ld$Ne) & ld$hiNe<8000],
       ld$yr[!is.na(ld$Ne) & ld$hiNe<8000],ld$hiNe[!is.na(ld$Ne) & ld$hiNe<8000],
       angle=90,length=0,col="darkorange")
arrows(ld$yr[!is.na(ld$Ne) & ld$hiNe==8000],ld$lowNe[!is.na(ld$Ne) & ld$hiNe==8000],
       ld$yr[!is.na(ld$Ne) & ld$hiNe==8000],ld$hiNe[!is.na(ld$Ne) & ld$hiNe==8000],
       lty=5,angle=90,length=0,col="darkorange")

#add dashed line for Inf
abline(h=8000,lty=2)

#ylab
par(las = 2)  
mtext(expression("LD " * hat(N)[e]), side = 2, line = 2.8) 
par(las = 0) 

#x-axis ticks
axis(side=1,at=setdiff(1995:2020, seq(1995, 2020, by = 5)),labels=F,tick=T,tcl = -0.3)
#log y-axis
major_ticks <- 10^(2:3)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 2:3)),las=1)
minor_ticks <- unlist(sapply(2:5, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)
axis(2,at=8000,label=expression(infinity),tick=T,las=1,cex.axis=1.5)

dev.off()




#### Fig S3. temporal Ne over time (initial & adjusted) ####


## Import data
n1 <- read.csv("ByGroup/tNe_toPlot_BY1995-1999_JR_0.02.csv")
n1 <- n1[n1$t==1,]
n2 <- read.csv("ByGroup/tNe_toPlot_BY2008-2020.csv")
n2 <- n2[n2$t==1,]
#visually check ranges for data to plot
n1
n2

#set Inf values to 18000
n1[sapply(n1, is.numeric)] <- lapply(n1[sapply(n1, is.numeric)], function(x) {
  x[is.infinite(x)] <- 10000
  x[x>10000] <- 9999 #to distinguish for plotting
  return(x)
})
n2[sapply(n2, is.numeric)] <- lapply(n2[sapply(n2, is.numeric)], function(x) {
  x[is.infinite(x)] <- 10000
  x[x>10000] <- 9999 #to distinguish for plotting
  return(x)
})


## Plot
pdf("NeDeclinePub/Figures/FigS3_init_adj_tNe_OverTime.pdf",
    width=8,height=3.5)
layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(0.4, 1)) 
par(mar=c(4,5.1,1,1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)

## 1995-1999
plot(Ne~yr2,data=n1,log="y",xlim=c(1996,1999.3),ylim=c(20,10000),
     xaxt="n",yaxt="n",
     xlab="Year",ylab=expression("temporal " * hat(N)[e]),
     pch=16,col="darkgray")
arrows(n1$yr2[n1$uCIj<10000],n1$lCIj[n1$uCIj<10000],
       n1$yr2[n1$uCIj<10000],n1$uCIj[n1$uCIj<10000],
       angle=90,length=0,col="darkgray")
arrows(n1$yr2[n1$uCIj==10000],n1$lCIj[n1$uCIj==10000],
       n1$yr2[n1$uCIj==10000],n1$uCIj[n1$uCIj==10000],
       lty=5,angle=90,length=0,col="darkgray")
#add adjusted values, offset a little on x-axis
points(n1$yr2+0.2,n1$Ne_MaxTempAdjusted,
       pch=16,col="steelblue4")
arrows(n1$yr2[n1$upperCI_MaxTempAdjusted<10000]+0.2,n1$lowerCI_MaxTempAdjusted[n1$upperCI_MaxTempAdjusted<10000],
       n1$yr2[n1$upperCI_MaxTempAdjusted<10000]+0.2,n1$upperCI_MaxTempAdjusted[n1$upperCI_MaxTempAdjusted<10000],
       angle=90,length=0,col="steelblue4")
arrows(n1$yr2[n1$upperCI_MaxTempAdjusted==10000]+0.2,n1$lowerCI_MaxTempAdjusted[n1$upperCI_MaxTempAdjusted==10000],
       n1$yr2[n1$upperCI_MaxTempAdjusted==10000]+0.2,n1$upperCI_MaxTempAdjusted[n1$upperCI_MaxTempAdjusted==10000],
       lty=5,angle=90,length=0,col="steelblue4")

#add dashed line for Inf
abline(h=10000,lty=2)

#x-axis
axis(side=1,at=1996:1999,labels=T,tick=T,las=2)
#log y-axis
major_ticks <- 10^(2:3)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 2:3)),las=1)
minor_ticks <- unlist(sapply(2:3, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)
axis(2,at=10000,label=expression(infinity),tick=T,las=1,cex.axis=1.5)

## Add panel label - CHANGE to side 3
mtext("(a)",font=2,side=3,line=0,adj=-0.7)

## 2008-2020
par(mar=c(4,2,1,5.1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)

plot(Ne~yr2,data=n2,log="y",xlim=c(2009,2020.3),ylim=c(20,10000),
     xaxt="n",yaxt="n",
     xlab="Year",ylab=expression("temporal " * hat(N)[e]),
     pch=16,col="darkgray")
arrows(n2$yr2[n2$uCIj<10000],n2$lCIj[n2$uCIj<10000],
       n2$yr2[n2$uCIj<10000],n2$uCIj[n2$uCIj<10000],
       angle=90,length=0,col="darkgray")
arrows(n2$yr2[n2$uCIj==10000],n2$lCIj[n2$uCIj==10000],
       n2$yr2[n2$uCIj==10000],n2$uCIj[n2$uCIj==10000],
       lty=5,angle=90,length=0,col="darkgray")

#add adjusted values, offset a little on x-axis
points(n2$yr2+0.2,n2$Ne_MaxTempAdjusted,
       pch=16,col="steelblue4")
arrows(n2$yr2[n2$upperCI_MaxTempAdjusted<10000]+0.2,n2$lowerCI_MaxTempAdjusted[n2$upperCI_MaxTempAdjusted<10000],
       n2$yr2[n2$upperCI_MaxTempAdjusted<10000]+0.2,n2$upperCI_MaxTempAdjusted[n2$upperCI_MaxTempAdjusted<10000],
       angle=90,length=0,col="steelblue4")
arrows(n2$yr2[n2$upperCI_MaxTempAdjusted==10000]+0.2,n2$lowerCI_MaxTempAdjusted[n2$upperCI_MaxTempAdjusted==10000],
       n2$yr2[n2$upperCI_MaxTempAdjusted==10000]+0.2,n2$upperCI_MaxTempAdjusted[n2$upperCI_MaxTempAdjusted==10000],
       lty=5,angle=90,length=0,col="steelblue4")

#add dashed line for Inf
abline(h=10000,lty=2)

#x-axis
axis(side=1,at=2009:2020,labels=T,tick=T,las=2)
#log y-axis
major_ticks <- 10^(2:3)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 2:3)),las=1)
minor_ticks <- unlist(sapply(2:3, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)
axis(2,at=10000,label=expression(infinity),tick=T,las=1,cex.axis=1.5)
#add legend
par(xpd=T)
legend("right",legend=c("Initial","Adjusted"),
       pch=16,col=c("darkgray","steelblue4"),cex=0.8,inset=c(-0.22,0))

## Add panel label - CHANGE to side 3
mtext("(b)",font=2,side=3,line=0,adj=-0.1)

dev.off()


#### Fig S4. all 3 Ne over time ####

## Import data
n <- read.csv("Ne_Nc_DeltaSmelt_2025nov14.csv")

#Set Inf values for plotting
#Fig 1: Inf=5000
#LD Ne: Inf=8000
#tNe:   Inf=10000 go with this once since it's the largest

n[sapply(n, is.numeric)] <- lapply(n[sapply(n, is.numeric)], function(x) {
  x[is.infinite(x)] <- 10000
  x[x>10000] <- 9999 #to distinguish for plotting
  return(x)
})

#split dataset in 2 year groups
n1 <- n[n$yr %in% 1996:1999,]
n2 <- n[n$yr %in% 2009:2020,]

#plot
pdf("Figures/FigS4_allNe_OverTime_2panel.pdf",
    width=9,height=3.5)
layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(0.4, 1))

## 1996-1999
par(mar=c(4,4,1,1)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)
plot(BestNe~yr,data=n1,col="white",log="y",xlim=c(1995.7,1999.3),ylim=c(50,10000),
     xaxt="n",yaxt="n",xlab="",ylab="") #set up plotting area
#ylab
par(las = 2)
mtext(expression(hat(N)[e]), side = 2, line = 2.8)
par(las = 0)

#add dashed line for Inf
abline(h=10000,lty=2)

#LD Ne
points(n1$yr-0.2,n1$LDNe,pch=17,col="darkorange",cex=0.9)
arrows(n1$yr[n1$upperCI_LDNe<10000]-0.2,n1$lowerCI_LDNe[n1$upperCI_LDNe<10000],
       n1$yr[n1$upperCI_LDNe<10000]-0.2,n1$upperCI_LDNe[n1$upperCI_LDNe<10000],
       angle=90,length=0,col="darkorange")
arrows(n1$yr[n1$upperCI_LDNe==10000]-0.2,n1$lowerCI_LDNe[n1$upperCI_LDNe==10000],
       n1$yr[n1$upperCI_LDNe==10000]-0.2,n1$upperCI_LDNe[n1$upperCI_LDNe==10000],
       angle=90,length=0,col="darkorange",lty=5)

#temporal Ne
points(n1$yr,n1$tNe_MaxTempAdjusted,pch=15,col="steelblue4",cex=0.9)
arrows(n1$yr[n1$upperCI_tNe_MaxTempAdjusted<10000],n1$lowerCI_tNe_MaxTempAdjusted[n1$upperCI_tNe_MaxTempAdjusted<10000],
       n1$yr[n1$upperCI_tNe_MaxTempAdjusted<10000],n1$upperCI_tNe_MaxTempAdjusted[n1$upperCI_tNe_MaxTempAdjusted<10000],
       angle=90,length=0,col="steelblue4")
arrows(n1$yr[n1$upperCI_tNe_MaxTempAdjusted==10000],n1$lowerCI_tNe_MaxTempAdjusted[n1$upperCI_tNe_MaxTempAdjusted==10000],
       n1$yr[n1$upperCI_tNe_MaxTempAdjusted==10000],n1$upperCI_tNe_MaxTempAdjusted[n1$upperCI_tNe_MaxTempAdjusted==10000],
       angle=90,length=0,col="steelblue4",lty=5)

#combo Ne
points(n1$yr[n1$type_BestNe=="combo_Ne"]+0.2,n1$BestNe[n1$type_BestNe=="combo_Ne"],
       pch=16,col="darkolivegreen",cex=1) #combo Ne values
arrows(n1$yr[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe<10000]+0.2,n1$lowerCI_BestNe[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe<10000],
       n1$yr[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe<10000]+0.2,n1$upperCI_BestNe[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe<10000],
       length=0,col="darkolivegreen")
arrows(n1$yr[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe==10000]+0.2,n1$lowerCI_BestNe[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe==10000],
       n1$yr[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe==10000]+0.2,n1$upperCI_BestNe[n1$type_BestNe=="combo_Ne" & n1$upperCI_BestNe==10000],
       length=0,col="darkolivegreen",lty=5)

#x-axis ticks
axis(side=1,at=1996:1999,labels=T,tick=T,las=2)
#log y-axis
major_ticks <- 10^(2:3)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 2:3)),las=1)
minor_ticks <- unlist(sapply(2:4, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)
axis(2,at=10000,label=expression(infinity),tick=T,las=1,cex.axis=1.5)

## Add panel label
mtext("(a)",font=2,side=3,line=0,adj=-0.5)


## 2009-2020
par(mar=c(4,2,1,4.5)) #mar=c(5.1,4.1,4.1,2.1) #(bottom, left, top, right)
plot(BestNe~yr,data=n,col="white",log="y",xlim=c(2008.7,2020.3),ylim=c(50,10000),
     xaxt="n",yaxt="n",xlab="",ylab="") #set up plotting area

#add dashed line for Inf
abline(h=10000,lty=2)

#LD Ne
points(n2$yr-0.2,n2$LDNe,pch=17,col="darkorange",cex=0.9)
arrows(n2$yr[n2$upperCI_LDNe<10000]-0.2,n2$lowerCI_LDNe[n2$upperCI_LDNe<10000],
       n2$yr[n2$upperCI_LDNe<10000]-0.2,n2$upperCI_LDNe[n2$upperCI_LDNe<10000],
       angle=90,length=0,col="darkorange")
arrows(n2$yr[n2$upperCI_LDNe==10000]-0.2,n2$lowerCI_LDNe[n2$upperCI_LDNe==10000],
       n2$yr[n2$upperCI_LDNe==10000]-0.2,n2$upperCI_LDNe[n2$upperCI_LDNe==10000],
       angle=90,length=0,col="darkorange",lty=5)

#temporal Ne
points(n2$yr,n2$tNe_MaxTempAdjusted,pch=15,col="steelblue4",cex=0.9)
arrows(n2$yr[n2$upperCI_tNe_MaxTempAdjusted<10000],n2$lowerCI_tNe_MaxTempAdjusted[n2$upperCI_tNe_MaxTempAdjusted<10000],
       n2$yr[n2$upperCI_tNe_MaxTempAdjusted<10000],n2$upperCI_tNe_MaxTempAdjusted[n2$upperCI_tNe_MaxTempAdjusted<10000],
       angle=90,length=0,col="steelblue4")
arrows(n2$yr[n2$upperCI_tNe_MaxTempAdjusted==10000],n2$lowerCI_tNe_MaxTempAdjusted[n2$upperCI_tNe_MaxTempAdjusted==10000],
       n2$yr[n2$upperCI_tNe_MaxTempAdjusted==10000],n2$upperCI_tNe_MaxTempAdjusted[n2$upperCI_tNe_MaxTempAdjusted==10000],
       angle=90,length=0,col="steelblue4",lty=5)

#combo Ne
points(n2$yr[n2$type_BestNe=="combo_Ne"]+0.2,n2$BestNe[n2$type_BestNe=="combo_Ne"],
       pch=16,col="darkolivegreen",cex=1) #combo Ne values
arrows(n2$yr[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe<10000]+0.2,n2$lowerCI_BestNe[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe<10000],
       n2$yr[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe<10000]+0.2,n2$upperCI_BestNe[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe<10000],
       length=0,col="darkolivegreen")
arrows(n2$yr[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe==10000]+0.2,n2$lowerCI_BestNe[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe==10000],
       n2$yr[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe==10000]+0.2,n2$upperCI_BestNe[n2$type_BestNe=="combo_Ne" & n2$upperCI_BestNe==10000],
       length=0,col="darkolivegreen",lty=5)

#x-axis ticks
axis(side=1,at=2009:2020,labels=T,tick=T,las=2)
#log y-axis
major_ticks <- 10^(2:3)
axis(2, at = major_ticks, labels = parse(text = paste0("10^", 2:3)),las=1)
minor_ticks <- unlist(sapply(2:4, function(i) seq(2 * 10^i, 9 * 10^i, by = 10^i)))
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.3,las=1)
axis(2,at=10000,label=expression(infinity),tick=T,las=1,cex.axis=1.5)

#add legend
par(xpd=T)
legend("right",legend=c("LD","temporal","combo"),cex=0.8,pch=c(17,15,16),pt.cex=c(0.9,0.9,1),
       col=c("darkorange","steelblue4","darkolivegreen"),inset=c(-0.17,0))

## Add panel label
mtext("(b)",font=2,side=3,line=0,adj=-0.1)

dev.off()
