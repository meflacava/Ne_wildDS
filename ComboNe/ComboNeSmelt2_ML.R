## This version iterates a specified number of times
# Modified by MEFL in May 2025

#### Merge LD & MaxTemp-adjusted temporal Ne data for input ####


##### For one set of NeEstimator parameters ####

## Pick group of years you're targeting
#grp <- "BY1995-1999"
#grp <- "BY2008-2020"

## Import
ld <- read.csv("ByYear/LD_BY1995-2020_minMAF0.02_Dp5-50_acrossChr.csv")
tne <- read.csv(paste0("ByGroup/tNe_forComboNe_",grp,".csv"))
outfile <- paste0("ComboNe_",grp,".csv")


## Merge into new df
SmeltData <- tne[,c("yr","Ne","L","Sprime")]
names(SmeltData) <- c("yr","NeTemp","nTemp","STemp")
for (i in 1:nrow(SmeltData)){
  SmeltData$NeLD[i] <- ld$Ne[ld$yr==SmeltData$yr[i]]
  SmeltData$SLD[i] <- ld$w.H.mean[ld$yr==SmeltData$yr[i]]
  SmeltData$LociLD[i] <- ld$n.loci[ld$yr==SmeltData$yr[i]]
}
SmeltData

## Write formatted results for later import
#write.csv(SmeltData,outfile,row.names=F,quote=F)


##### For multiple NeEstimator parameters ####

# ## Pick group of years you're targeting
# #grp <- "BY1995-1999"
# #grp <- "BY2008-2020"
# 
# ## Choose which temporal Ne parameters to subset (if multiple)
# #method <- "Pk"
# #method <- "JR"
# #pcrit <- 0.02
# #pcrit <- 0.05
# outfile <- paste0("ComboNe_",grp,"_",method,"_",pcrit,".csv")
# 
# ## Import
# ld <- read.csv("ByYear/LD_BY1995-2020_minMAF0.02_Dp5-50_acrossChr.csv")
# tne <- read.csv(paste0("ByGroup/tNe_forComboNe_",grp,"_",method,"_",pcrit,".csv"))
# 
# ## Merge into new df
# SmeltData <- tne[,c("yr","Ne","L","Sprime")]
# names(SmeltData) <- c("yr","NeTemp","nTemp","STemp")
# for (i in 1:nrow(SmeltData)){
#   SmeltData$NeLD[i] <- ld$Ne[ld$yr==SmeltData$yr[i]]
#   SmeltData$SLD[i] <- ld$w.H.mean[ld$yr==SmeltData$yr[i]]
#   SmeltData$LociLD[i] <- ld$n.loci[ld$yr==SmeltData$yr[i]]
# }
# SmeltData
# 
# ## Write formatted results for later import
# #write.csv(SmeltData,paste0("ComboNe/ComboNe_",grp,"_in.csv"),row.names=F,quote=F)



##### Or import previously formatted data ####
# setwd("ComboNe")
# 
# ## Pick group of years you're targeting
# #grp <- "1995-1999"
# #grp <- "2008-2020"
# SmeltData = read.csv(paste0("ComboNe_",grp,"_in.csv"))



#### Calc combo Ne with CIs ####

TempLoci  = SmeltData$nTemp
TempS  = SmeltData$STemp
TempNe  = SmeltData$NeTemp
LDLoci = SmeltData$LociLD 
LDS  = SmeltData$SLD
LDNe  = SmeltData$NeLD
# TempLoci  = SmeltData$TempLoci
# TempS  = SmeltData$TempS
# TempNe  = SmeltData$TempNe
# LDLoci = SmeltData$LDLoci
# LDS  = SmeltData$LDS
# LDNe  = SmeltData$LDNe
Chr = 26
t=1
NIterate = 10

###########
GetFnprime <- function(Inputs)  {
  N = Inputs[1]  
  C = Inputs[2] 
  SS= Inputs[3]  
  Loci= Inputs[4] 
  
  p = 1:3
  ######
  ff <- function(x,p){
    x/(1/p[1]^p[3]+(x/p[2])^p[3])^(1/p[3])}
  ######
  
  q1 = c(0.958633285 ,0.013469143 , 0.007135748, -0.001863156) 
  q2 = c(2.08304896, -0.24271278, 0.00923448, -0.57920739, 1.21567379, 0.05929167 , 0.12727216, 0.09241921) 
  q5 = c(-8.8963881, 2.7271474, 1.8244296, 3.8598313, -0.1480143, -0.4126190, -0.3886260 )
  
  p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(C)*log(N)  
  p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(SS) + q2[5]*SS/N + q2[6]*log(C)*log(N) + q2[7]*log(SS)*log(N) + q2[8]*log(SS)*log(C)
  p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N) + q5[7]*log(C)*log(SS)  
  A = ff(log10(Loci),p)
  Enprime = 10^A
  
  #### ensure predicted nprime does not exceed number of loci
  Enprime = min(Enprime, Loci)
  
  return(Enprime)
} # end function  
############

###########
GetLDNprime <- function(Inputs)  {
  N = Inputs[1]  
  C = Inputs[2] 
  SS= Inputs[3]  
  Loci= Inputs[4] 
  
  p = 1:3
  ######
  ff <- function(x,p){
    x/(1/p[1]^p[3]+(x/p[2])^p[3])^(1/p[3])}
  ######
  
  q1 = c(1.53450230,  0.12518487,  0.06692614, -0.06710636,-0.24470293, -0.01474650) 
  q2 = c(0.52949689,0.12559264,0.61411490,0.24493530,0.50730040,-1.38178677,-0.01697811,-0.03484065)
  q5 =  c(1.2317425, 0.5675758,-0.3671743, 1.2211834, 0.3815077, 0.1881963,-0.6756503 )
  p[1] = q1[1] + q1[2]*log(C) + q1[3]*log(N) + q1[4]*log(SS) + q1[5]*SS/N + q1[6]*log(C)*log(N)  
  p[2] = q2[1] + q2[2]*log(C) + q2[3]*log(N) + q2[4]*log(SS) + q2[5]*SS/N + q2[6]/C + q2[7]*log(C)*log(N) + q2[8]*log(SS)*log(N)
  p[3] = q5[1] + q5[2]*log(C) + q5[3]*log(N) + q5[4]*log(SS) + q5[5]*log(C)*log(N) + q5[6]*log(SS)*log(N) + q5[7]*log(C)*log(SS)  
  A = ff(log10(Loci),p)
  Enprime = 10^A
  
  #### ensure predicted nprime does not exceed number of locus pairs
  Enprime = min(Enprime, Loci*(Loci-1)/2)
  
  return(Enprime)
} # end function  
############

testNe = SmeltData$NeTemp
testNe[testNe < 0 | testNe > 2000] <- 2000

lnR =  0.627 - 0.477*log(testNe) + 0.0304*log(testNe)^2 + 0.169*log(SmeltData$STemp) - 0.0236*log(SmeltData$nTemp)
Factor =  exp(lnR)
EndR = 0.2688 - 0.2042*log(testNe) + 0.01079*log(testNe)^2 + 0.09408*log(SmeltData$STemp) - 0.0128*log(SmeltData$nTemp)
Factor[1] = exp(EndR[1])
Factor[nrow(SmeltData)] = exp(EndR[nrow(SmeltData)])
Factor

Results = matrix(NA,nrow(SmeltData),9)
colnames(Results) = c("Year","Ne0","Ne1","Ne2","Ne3","FinalNe","Var","LowerCI","UpperCI")
VLD = rep(NA,nrow(SmeltData))
VTemp = VLD
ComboVar = VLD
Ne0 = VLD

nprimeLD = matrix(NA,nrow(SmeltData),NIterate)
nprimeTemp = nprimeLD
HMNe = nprimeLD
WT = nprimeLD
WR = nprimeLD

InputsLD = rep(NA,4)
InputsTemp = InputsLD

NeLD2 = SmeltData$NeLD
NeLD2[NeLD2 < 0 | NeLD2 > 5000] <- 5000
NeT2 = SmeltData$NeTemp
NeT2[NeT2 < 0 | NeT2 > 5000] <- 5000

##First, get unweighted harmonic mean Ne to represent common parameter both methods are seeking to estimate
HMNe0 = 2/(1/SmeltData$NeLD + 1/SmeltData$NeTemp)

for (j in 1:nrow(SmeltData))  {
  
  k = 1
  
  ##Get initial nprime
  InputsLD[1] = HMNe0[j]
  if(HMNe0[j] < 0 | HMNe0[j] > 5000) { InputsLD[1] = 5000} 
  InputsLD[2] = Chr
  InputsLD[3] = LDS[j]
  InputsLD[4] = LDLoci[j]
  nprimeLD[j,k] = GetLDNprime(InputsLD)
  
  InputsTemp[1] = HMNe0[j]
  if(HMNe0[j] < 0 | HMNe0[j] > 5000) { InputsTemp[1] = 5000} 
  InputsTemp[2] = Chr
  InputsTemp[3] =TempS[j]
  InputsTemp[4] = TempLoci[j]
  nprimeTemp[j,k]  = GetFnprime(InputsTemp)
  
  EF = 1/(2*HMNe0[j]) + 1/SmeltData$STemp[j]
  EVF = (2*EF^2/nprimeTemp[j,k])*Factor[j]^2
  EVNeF = 4*EVF
  ER = 1/(3*HMNe0[j]) + 1/(SmeltData$STemp[j]-1)
  EVR = 2*ER^2/nprimeLD[j,k]
  EVNeR = 9*EVR
  #EVNeF/EVNeR
  
  ## Get initial weights
  WT[j,k] = 1/EVNeF 
  WR[j,k] = 1/EVNeR
  WTot = WT[j,k] + WR[j,k]
  WT[j,k] = WT[j,k]/WTot
  WR[j,k] = WR[j,k]/WTot
  # WT[j,k] + WR[j,k]  # should be 1
  
  ###Get first weighted ComboNe
  HMNe[j,k] = 1/(WR[j,k]/SmeltData$NeLD[j] + WT[j,k]/SmeltData$NeTemp[j])
  
  
  for (k in 2:NIterate)  {
    
    ##Update nprime with new Nehat
    InputsLD[1] = HMNe[j,(k-1)]
    if(HMNe[j,(k-1)] < 0 | HMNe[j,(k-1)] > 5000) { InputsLD[1] = 5000} 
    InputsLD[2] = Chr
    InputsLD[3] = LDS[j]
    InputsLD[4] = LDLoci[j]
    nprimeLD[j,k] = GetLDNprime(InputsLD)
    
    InputsTemp[1] = InputsLD[1]
    InputsTemp[2] = Chr
    InputsTemp[3] =TempS[j]
    InputsTemp[4] = TempLoci[j]
    nprimeTemp[j,k]  = GetFnprime(InputsTemp)
    
    EF = 1/(2*HMNe[j,(k-1)]) + 1/SmeltData$STemp[j]
    EVF = (2*EF^2/nprimeTemp[j,k])*Factor[j]^2
    EVNeF = 4*EVF
    
    ER = 1/(3*HMNe[j,(k-1)]) + 1/(SmeltData$STemp[j]-1)
    EVR = 2*ER^2/nprimeLD[j,k]
    EVNeR = 9*EVR
    
    ## Get updated weights
    WT[j,k] = 1/EVNeF 
    WR[j,k] = 1/EVNeR
    WTot = WT[j,k] + WR[j,k]
    WT[j,k] = WT[j,k]/WTot
    WR[j,k] = WR[j,k]/WTot
    
    ###Update weighted ComboNe
    HMNe[j,k] = 1/(WR[j,k]/SmeltData$NeLD[j] + WT[j,k]/SmeltData$NeTemp[j])
    
  } # end for k
  
  VTemp[j] = 8*(1/(2*HMNe[j,NIterate]) + 1/TempS[j])^2/nprimeTemp[j,NIterate]
  VTemp[j] = VTemp[j]*Factor[j]^2
  VLD[j] = 18*(1/(3*HMNe[j,NIterate]) + 1/(LDS[j]-1))^2/nprimeLD[j,NIterate]
  ComboVar[j] = WT[j,NIterate]^2*VTemp[j] + WR[j,NIterate]^2*VLD[j]
  lowInv = 1/HMNe[j,NIterate] - 1.96*sqrt(ComboVar[j])
  hiInv = 1/HMNe[j,NIterate] + 1.96*sqrt(ComboVar[j])
  lowCI = 1/hiInv
  hiCI = 1/lowInv
  Results[j,] = c(SmeltData$yr[j],HMNe0[j],HMNe[j,1],HMNe[j,2],HMNe[j,3],HMNe[j,NIterate],ComboVar[j],lowCI,hiCI)
  
} # end for j

Results
HMNe # rows are years; columns are iterations
HMNe[,10]/HMNe[,5] # proportional change over last 5 iterations


#### Write results ####
#note - I added the weights for LD vs. temporal to the output
out <- data.frame(Results[,c(1,6,8,9)],WR[,10])
names(out) <- c("yr","ComboNe","lowerCI","upperCI","LD_weight")
out$temporal_weight <- 1 - out$LD_weight
out
#write.csv(out,paste0("ComboNe/",outfile),quote=F,row.names=F)

#range of weights for LD method
range(WR)
#1995-1999: 0.6508264 0.9167398
#2008-2020: 0.6245516 0.9458352
