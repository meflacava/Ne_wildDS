### MaxTemp general script
# July 2024 MEFL annotated script and added plotting code

#### Import data ####

## Raw output from NeEstimator
# - Run temporal Ne with Pollak's estimator, check box to output files in tabular format
#   - pcrit=0.05, delete others and uncheck run without freq restriction
#   - Under Methods->Temporal, Plan II, specify generations in file (e.g., 0,1,2,3,4 if you are
#     running 5 consecutive years (or 0,2,3,4 if you're missing the 2nd year of data))
#   - Options->population range to run = specify if you don't want to run all pops in file
#     (e.g., 2-4 if your file has more generations but you only want to analyze 2-4). NOTE that
#     here the numbers refer to the pops as listed, so if you have 5 pops in your input genepop
#     file, even if they equate temporally to generations 0,1,3,4,5, you need to specify to run
#     2-4 if you want 1,3,4 run
# - Use command line to extract results lines from NeEstimator output
# - Use R code in ByGroup_temporalNe.R to further modify data for MaxTemp
# - Important columns in results:
#   - yr1 & yr2 = actual years in pairwise comparison (e.g., 2013 & 2014)
#   - gen1 & gen2 = generations in set, starting from 1 (e.g., 1 & 2)
#   - t = yr2 - yr1
#   - L = # loci (copy from NeEstimator output)
#   - F = variance in allele frequencies (copy from NeEstimator output)
#   - Fprime = adjust F for sample size (copy from NeEstimator output)
#   - Ne = estimated effective pop size (copy from NeEstimator output)

##### For one set of NeEstimator parameters ####

## Choose group of years you're targeting
#grp <- "BY1995-1999"
#grp <- "BY2008-2020"

Vitals = read.csv(file=paste0("PreMaxTemp_",grp,".csv"))
outfile1 <- paste0("tNe_forComboNe_",grp,".csv") #use this if not choosing temporal Ne parameters
outfile2 <- paste0("tNe_toPlot_",grp,".csv") #use this if not choosing temporal Ne parameters

## Add Sprime (harmonic mean sample size) - if it's not already in data frame
# - Calc from F and F', don't use S output from NeEstimator
#Vitals$Sprime <- 1/(Vitals$F-Vitals$Fprime)




#### Set up ####

## Number of generations in dataset
NSamples = length(unique(Vitals$yr1)) + 1

## Set up F prime matrices
FprimeMatrix  = matrix(NA,NSamples,NSamples)
for(j in 1:nrow(Vitals)){
  FprimeMatrix[Vitals$gen1[j],Vitals$gen2[j]] = Vitals$Fprime[j]  ## Data spanning all generations
}
FprimeMatrix2 = FprimeMatrix
FprimeMatrix3 = FprimeMatrix


#### One generation F ####
OneGen = subset(Vitals,Vitals$t == 1)  ## Data for single generations
InitialF = OneGen$Fprime



#### 1. Adjust one gen Fprime using two gen Fprime ####

## Iterate diagonals (AKA 1st adjustment of one gen F prime)
NewF = rep(NA,length(InitialF))
#first generation (can only use next gen)
j = 1
addF = FprimeMatrix[j,j+2] - InitialF[j+1] #F3_1 - initialF3
NewF[j] = 0.6*InitialF[j] + 0.4*addF
#middle generations (can use gens on each side of target gen)
q = length(InitialF)-1
for (j in 2:q)  {
addF1 = FprimeMatrix[j-1,j+1] - InitialF[j-1]
addF2 = FprimeMatrix[j,j+2] - InitialF[j+1]
NewF[j] = 0.5*InitialF[j] + 0.25*(addF1 + addF2)
}
#last generation (can only use previous gen)
j = length(InitialF)
addF = FprimeMatrix[j-1,j+1] - InitialF[j-1]
NewF[j] = 0.6*InitialF[j] + 0.4*addF

rbind(NewF,InitialF) #view change in F prime

## Update F prime matrix 2 with new one gen F prime values
for (j in 1:length(NewF)){
  FprimeMatrix2[j,j+1] = NewF[j]
}

## Adjust Ne^ after 1 iteration
NewNe = 1/(2*NewF)
#rbind(OneGen$Ne,NewNe) #view change in Ne^



#### 2. Adjust 2+ gen F prime values to then adjust one gen F prime a second time ####

## Iterate off-diagonals (AKA 2+ gen F prime)
for (igen in 2:(length(NewF)-2)){ #adjust two gen F prime, then three gen, ... all 2+ gen F prime values in dataset 
  #first generation (can only use next gen)
  j = 1 #Example (igen=2,j=1): F3_1 = 0.6*F3_1 + 0.4*(F4_1 - adjF4)
  FprimeMatrix2[j,j+igen] = 0.6*FprimeMatrix[j,j+igen] + 0.4*(FprimeMatrix[j,j+igen+1]-NewF[j+igen])
  #print(paste0("F",j,"_",j+igen)) #see which F prime values are being updated
  #middle generations (can use gens on each side of target gen)
  q = length(InitialF)-igen
  for (j in 2:q){ #Example (igen=2,j=2): F4_2 = 0.25*(F4_1 - adjF1) + 0.5*(F4_2) + 0.25*(F5_2 - adjF4)
    FprimeMatrix2[j,j+igen] = 0.25*(FprimeMatrix[j-1,j+igen]-NewF[j-1]) + 0.5*FprimeMatrix[j,j+igen] + 0.25*(FprimeMatrix[j,j+igen+1]-NewF[j+igen])
    #print(paste0("F",j,"_",j+igen)) #see which F prime values are being updated
  }
  #last generation (can only use previous gen)
  j = length(InitialF)-igen+1 #Example (igen=2,j=6): F8_6 = 0.4*(F8_5 - adjF5) + 0.6*(F8_6)
  FprimeMatrix2[j,j+igen] = 0.4*(FprimeMatrix[j-1,j+igen]-NewF[j-1]) + 0.6*FprimeMatrix[j,j+igen]
  #print(paste0("F",j,"_",j+igen))  #see which F prime values are being updated
}
#update the largest F prime possible (for 2013-2020 data, it's six gen, first F7_1 then F8_2)
igen = length(NewF)-1
j = 1
FprimeMatrix2[j,j+igen] = 0.6*FprimeMatrix[j,j+igen] + 0.4*(FprimeMatrix[j,j+igen+1]-NewF[j+igen])
j = length(InitialF)-igen+1
FprimeMatrix2[j,j+igen] = 0.4*(FprimeMatrix[j-1,j+igen]-NewF[j-1]) + 0.6*FprimeMatrix[j,j+igen]

## 2nd adjustment of one gen F prime values using newly adjusted 2+ gen values
NewF2 = rep(NA,length(NewF))
#first generation (can only use next gen)
j = 1
addF = FprimeMatrix2[j,j+2] - NewF[j+1] #F3_1 - adjF2
NewF2[j] = 0.6*NewF[j] + 0.4*addF #2nd adjF1 based on above
#middle generations (can use gens on each side of target gen)
q = length(NewF)-1
for (j in 2:q)  {
  addF1 = FprimeMatrix2[j-1,j+1] - NewF[j-1] #Example: F3_1 - adjF1
  addF2 = FprimeMatrix2[j,j+2] - NewF[j+1] #Example: F4_2 - adjF3
  NewF2[j] = 0.5*NewF[j] + 0.25*(addF1 + addF2) #Example: 2nd adjF2 based on both above
}
#last generation (can only use previous gen)
j = length(NewF)
addF = FprimeMatrix2[j-1,j+1] - NewF[j-1] #Example: F8_6 - adjF6
NewF2[j] = 0.6*NewF[j] + 0.4*addF #Example: 2nd adjF7 based on both above

rbind(InitialF,NewF,NewF2) #view change in F prime

## Update F prime matrix 3 with new one gen F prime values
for (j in 1:length(NewF)){
  FprimeMatrix3[j,j+1] = NewF2[j]
}

## Adjust Ne^ after 2 iterations
NewNe2 = 1/(2*NewF2)
rbind(OneGen$Ne,NewNe,NewNe2) #view change in Ne^



#### 3. Adjust 2+ gen F prime values again  to then adjust one gen F prime a 3rd time ####

## Iterate off-diagonals (AKA 2+ gen F prime)
for (igen in 2:(length(NewF)-2)){ #adjust two gen F prime, then three gen, ... all 2+ gen F prime values in dataset
  #first generation (can only use next gen)
  j = 1
  FprimeMatrix3[j,j+igen] = 0.6*FprimeMatrix2[j,j+igen] + 0.4*(FprimeMatrix2[j,j+igen+1]-NewF2[j+igen])
  #middle generations (can use gens on each side of target gen)
  q = length(InitialF)-igen
  for (j in 2:q){
    FprimeMatrix3[j,j+igen] = 0.25*(FprimeMatrix2[j-1,j+igen]-NewF2[j-1]) + 0.5*FprimeMatrix2[j,j+igen] + 0.25*(FprimeMatrix2[j,j+igen+1]-NewF2[j+igen])
    }
  #last generation (can only use previous gen)
  j = length(InitialF)-igen+1
  FprimeMatrix3[j,j+igen] = 0.4*(FprimeMatrix2[j-1,j+igen]-NewF2[j-1]) + 0.6*FprimeMatrix2[j,j+igen]
  }
#update the largest F prime possible (for 2013-2020 data, it's six gen, first F7_1 then F8_2)
igen = length(NewF)-1
j = 1
FprimeMatrix3[j,j+igen] = 0.6*FprimeMatrix2[j,j+igen] + 0.4*(FprimeMatrix2[j,j+igen+1]-NewF2[j+igen])
j = length(InitialF)-igen+1
FprimeMatrix3[j,j+igen] = 0.4*(FprimeMatrix2[j-1,j+igen]-NewF2[j-1]) + 0.6*FprimeMatrix2[j,j+igen]

## 3rd adjustment of one gen F prime values using newly adjusted 2+ gen values
NewF3 = rep(NA,length(NewF))
#first generation (can only use next gen)
j = 1
addF = FprimeMatrix3[j,j+2] - NewF2[j+1]
NewF3[j] = 0.6*NewF2[j] + 0.4*addF
#middle generations (can use gens on each side of target gen)
q = length(NewF)-1
for (j in 2:q)  {
  addF1 = FprimeMatrix3[j-1,j+1] - NewF2[j-1]
  addF2 = FprimeMatrix3[j,j+2] - NewF2[j+1]
  NewF3[j] = 0.5*NewF2[j] + 0.25*(addF1 + addF2)
  }
#last generation (can only use previous gen)
j = length(NewF)
addF = FprimeMatrix3[j-1,j+1] - NewF2[j-1]
NewF3[j] = 0.6*NewF2[j] + 0.4*addF

#rbind(InitialF,NewF,NewF2,NewF3) #view change in F prime
#NewF/InitialF
#NewF2/NewF
#NewF3/NewF2
#var(NewF/InitialF)
#var(NewF2/NewF)
#var(NewF3/NewF2)

## Adjust Ne^ after 3 iterations
NewNe3 = 1/(2*NewF3)
rbind(OneGen$Ne,NewNe,NewNe2,NewNe3)


#### Calculate CIs ####

## Original CIs, assuming var(F^) = 2F^/n
df = OneGen$L
varFhat = 2*OneGen$F^2/df
sdFhat = sqrt(varFhat)
upperFhat = OneGen$F + 1.96*sdFhat
lowerFhat = OneGen$F - 1.96*sdFhat
upperFprime = upperFhat - 1/OneGen$Sprime
lowerFprime = lowerFhat - 1/OneGen$Sprime
rbind(upperFprime,OneGen$Fprime,lowerFprime)

upperNehat = 1/(2*lowerFprime)
lowerNehat = 1/(2*upperFprime)
rbind(upperNehat,OneGen$Ne, lowerNehat)

#save Ne and CIs as df
old <- data.frame(yr=OneGen$yr2,Ne=OneGen$Ne,upperNehat=upperNehat,lowerNehat=lowerNehat)


##New CIs

#REPLACE negative Ne values or Ne <2000 with 2000 to calculate CIs
# Robin comment 8/9/24: "Our model fits only considered finite Ne values up to 2000.  
# I will need to discuss this with Michele.  But our results indicate the effect of 
# Ne asymptotes near Ne = 1000.  So for time being, if any adjusted Ne is negative or 
# > 2000, replace with 2000.  That should give a realistic value for Factor, which 
# indicates how much sd(F) is reduced.
NewNe3_pos <- NewNe3
NewNe3_pos[NewNe3_pos > 2000 | NewNe3_pos < 0] <- 2000


NewvarF = 2*(NewF3+1/OneGen$Sprime)^2/df

lnR =  0.627 - 0.477*log(NewNe3_pos) + 0.0304*log(NewNe3_pos)^2 + 0.169*log(OneGen$Sprime) - 0.0236*log(OneGen$L)
Factor =  exp(lnR)
EndR = 0.2688 - 0.2042*log(NewNe3_pos) + 0.01079*log(NewNe3_pos)^2 + 0.09408*log(OneGen$Sprime) - 0.0128*log(OneGen$L)
Factor[1] = exp(EndR[1])
Factor[length(InitialF)] = exp(EndR[length(InitialF)])

NewsdF = Factor*sqrt(NewvarF)
NewupperFprime = NewF3 + 1.96*NewsdF
NewlowerFprime = NewF3 - 1.96*NewsdF
rbind(NewupperFprime,NewF3,NewlowerFprime)

NewupperNehat = 1/(2*NewlowerFprime)
NewlowerNehat = 1/(2*NewupperFprime)
rbind(NewupperNehat,NewNe3, NewlowerNehat)




#### Export results ####

## Export for combo Ne
#save Ne and CIs as df
new <- data.frame(yr=OneGen$yr2,Ne=NewNe3,lowerCI=NewlowerNehat,upperCI=NewupperNehat)
new$L <- Vitals$L[Vitals$t==1]
new$Sprime <- Vitals$Sprime[Vitals$t==1]
write.csv(new,outfile1,quote=F,row.names=F)

## Export for Fig S3
all <- Vitals
all$Ne_MaxTempAdjusted <- NA
all$lowerCI_MaxTempAdjusted <- NA
all$upperCI_MaxTempAdjusted <- NA

for (i in 1:nrow(all)){
  if (all$t[i]==1){
    all$Ne_MaxTempAdjusted[i] <- new$Ne[new$yr==all$yr2[i]]
    all$lowerCI_MaxTempAdjusted[i] <- new$lowerCI[new$yr==all$yr2[i]]
    all$upperCI_MaxTempAdjusted[i] <- new$upperCI[new$yr==all$yr2[i]]
  }
}

# Set neg values to Inf
all[sapply(all, is.numeric)] <- lapply(all[sapply(all, is.numeric)], function(x) {
  x[x<0] <- Inf
  return(x)
})

# Export for reporting & plotting
write.csv(all,outfile2,quote=F,row.names=F)


## Export as pre and post values
# old$dataset <- "pre-MaxTemp"
# new$dataset <- "post-MaxTemp"
# write.csv(rbind(old,new),paste0("PostMaxTemp_BY",grp,".csv"),quote=F,row.names=F)
