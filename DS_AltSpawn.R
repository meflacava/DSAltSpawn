#Delta smelt alternative spawning strategy experiments
# R script by ML


#### Import data ####
tank <- read.csv("https://raw.githubusercontent.com/meflacava/DSAltSpawn/main/FinalGenotypedOffspring.csv")
parents <- read.csv("https://raw.githubusercontent.com/meflacava/DSAltSpawn/main/ParentTanks.csv")
pc <- read.csv("https://raw.githubusercontent.com/meflacava/DSAltSpawn/main/PairCrosses.csv")
lar <- read.csv("https://raw.githubusercontent.com/meflacava/DSAltSpawn/main/EstLarvaePerTank.csv")



#### Calculate contributions of PCs and parents ####

## offspring per PC
pc$prog <- NA
pc$tank.prop <- NA
for (i in 1:nrow(pc)){ 
  x <- tank[tank$expt==pc$expt[i] & tank$tank==pc$tank[i],]
  pc$prog[i] <- nrow(x[x$Candidate.mother.ID==pc$dam[i] & x$Candidate.father.ID==pc$sire[i],])
  pc$tank.prop[i] <- pc$prog[i]/nrow(x)
}
pc
sum(pc$prog,na.rm=T) == nrow(tank) #check all offspring in tank were assigned to PCs


#add index of 10 replicates within each spawning type
pc$rep <- rep(rep(1:10,each=9),2)


## offspring per parent
parents$prog <- NA
for (i in 1:nrow(parents)){ 
  x <- tank[tank$expt==parents$expt[i] & tank$tank==parents$tank[i],]
  parents$prog[i] <- sum(x$Candidate.mother.ID==parents$ParentID[i]) + sum(x$Candidate.father.ID==parents$ParentID[i])
  parents$tank.prop[i] <- parents$prog[i]/nrow(x)
  }
parents
sum(parents$prog) == 2*nrow(tank) #NOTE: this is 2x the # of genotyped offspring because each offspring is counted for both mother and father

#add index of 10 replicates within each spawning type
parents$rep <- rep(rep(1:10,each=6),2)




#### Summary ####

#Number of offspring genotyped per tank
table(tank$SpawnType,tank$rep)
#           1  2  3  4  5  6  7  8  9 10
#pooled    91 91 91 90 91 90 90 61 87 91
#factorial 90 91 89 89 91 90 89 88 90 90

#Number of failed PCs per strategy
nrow(pc[pc$SpawnType=="pooled" & pc$tank.prop==0,]) #25 of 90 pooled PCs (28%)
nrow(pc[pc$SpawnType=="factorial" & pc$tank.prop==0,]) #2 of 90 factorial PCs (2%)

#Number of failed parents per strategy
nrow(parents[parents$SpawnType=="pooled" & parents$sex=="F" & parents$tank.prop==0,]) #4 of 30 pooled dams (13%)
nrow(parents[parents$SpawnType=="pooled" & parents$sex=="M" & parents$tank.prop==0,]) #2 of 30 pooled sires (7%)
nrow(parents[parents$SpawnType=="factorial" & parents$sex=="F" & parents$tank.prop==0,]) #0 of 30 factorial dams
nrow(parents[parents$SpawnType=="factorial" & parents$sex=="M" & parents$tank.prop==0,]) #0 of 30 factorial sires

#Range of contributions by PCs
range(pc$tank.prop[pc$SpawnType=="pooled"]) #0-0.80
range(pc$tank.prop[pc$SpawnType=="factorial"]) #0-0.27

#Range of contributions by parents
range(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="F"]) #0-1
range(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="M"]) #0-0.98
range(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="F"]) #0.01-0.59
range(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="M"]) #0.19-0.50

#median contributions by PCs
median(pc$tank.prop[pc$SpawnType=="pooled"]) #0.05
median(pc$tank.prop[pc$SpawnType=="factorial"]) #0.10

#Median contributions by parents
median(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="F"]) #0.29
median(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="M"]) #0.32
median(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="F"]) #0.33
median(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="M"]) #0.32


#### Test for differences between strategies ####

#Wilcoxon rank-sum test
# - paired=F causes the Wilcoxon test to be a rank-sum test (AKA Mann-Whitney test) to work with unpaired data

#PC
wilcox.test(pc$tank.prop[pc$SpawnType=="pooled"],pc$tank.prop[pc$SpawnType=="factorial"],
            paired=F,alternative="two.sided")
#    Wilcoxon rank sum test with continuity correction
#
#data:  pc$tank.prop[pc$SpawnType == "pooled"] and pc$tank.prop[pc$SpawnType == "factorial"]
#W = 3127.5, p-value = 0.008227
#alternative hypothesis: true location shift is not equal to 0

#dams
wilcox.test(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="F"],
            parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="F"],
            paired=F,alternative="less")
# Wilcoxon rank sum test with continuity correction
#
# data:  parents$tank.prop[parents$SpawnType == "pooled" & parents$sex == "F"] and parents$tank.prop[parents$SpawnType == "factorial" & parents$sex == "F"]
# W = 404, p-value = 0.5011
# alternative hypothesis: true location shift is not equal to 0


#sires
wilcox.test(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="M"],
            parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="M"],
            paired=F,alternative="two.sided")
# Wilcoxon rank sum test with continuity correction
#
# data:  parents$tank.prop[parents$SpawnType == "pooled" & parents$sex == "M"] and parents$tank.prop[parents$SpawnType == "factorial" & parents$sex == "M"]
# W = 429.5, p-value = 0.7674
# alternative hypothesis: true location shift is not equal to 0




#### Chi-square test for deviation from expected even contribution ####
#Chi-square test - compare observed to expected frequencies 
# Expected/null = equal contribution of PCs to offspring, so leave blank because that's the default

## Contributions of PCs
chi.pc <- data.frame(SpawnType=c(rep("factorial",10),rep("pooled",10)),rep=rep(1:10,2))
for (i in 1:nrow(chi.pc)){
  alt <- pc$prog[pc$SpawnType==chi.pc$SpawnType[i] & pc$rep==chi.pc$rep[i]]
  chi.pc$pvalue[i] <- format(chisq.test(alt,y=NULL, correct=T)$p.value,scientific=F)
}
chi.pc
# p value < 0.05 indicates deviation from chi-square distribution, so deviation from expected even proportions

#proportion of tanks where 9 PCs had significantly different contributions
sum(chi.pc$pvalue[chi.pc$SpawnType=="factorial"]<0.05)/nrow(chi.pc[chi.pc$SpawnType=="factorial",])
#0.7 #AKA 70% of factorial tanks deviated from expected equal contributions of PCs
sum(chi.pc$pvalue[chi.pc$SpawnType=="pooled"]<0.05)/nrow(chi.pc[chi.pc$SpawnType=="pooled",])
#1 #AKA all pooled tanks deviated from expected equal contributions of PCs


## Contributions of individual dams or sires

#dams
chi.dam <- data.frame(SpawnType=c(rep("factorial",10),rep("pooled",10)),rep=rep(1:10,2)) 
for (i in 1:nrow(chi.dam)){
  temp <- parents[parents$sex=="F" & parents$SpawnType==chi.dam$SpawnType[i] & parents$rep==chi.dam$rep[i],]
  alt <- temp$prog
  chi.dam$pvalue[i] <- format(chisq.test(alt,y=NULL,correct=T)$p.value,scientific=F)
}
chi.dam
# p value < 0.05 indicates deviation from chi-square distribution, so deviation from expected even proportions

#proportion of tanks where 3 dams had significantly different contributions
sum(chi.dam$pvalue[chi.dam$SpawnType=="factorial"]<0.05)/nrow(chi.dam[chi.dam$SpawnType=="factorial",])
#0.6 #AKA 60% of factorial tanks deviated from expected contributions of dams
sum(chi.dam$pvalue[chi.dam$SpawnType=="pooled"]<0.05)/nrow(chi.dam[chi.dam$SpawnType=="pooled",])
#0.9 #AKA 90% of pooled tanks deviated from expected contributions of dams

#sires
chi.sire <- data.frame(SpawnType=c(rep("factorial",10),rep("pooled",10)),rep=rep(1:10,2)) 
for (i in 1:nrow(chi.sire)){
  temp <- parents[parents$sex=="M" & parents$SpawnType==chi.sire$SpawnType[i] & parents$rep==chi.sire$rep[i],]
  alt <- temp$prog
  chi.sire$pvalue[i] <- format(chisq.test(alt,y=NULL,correct=T)$p.value,scientific=F)
}
chi.sire
# p value < 0.05 indicates deviation from chi-square distribution, so deviation from expected even proportions

#proportion of tanks where 3 sires had significantly different contributions
sum(chi.sire$pvalue[chi.sire$SpawnType=="factorial"]<0.05)/nrow(chi.sire[chi.sire$SpawnType=="factorial",])
#0.2 #AKA 20% of factorial tanks deviated from expected contributions of sires
sum(chi.sire$pvalue[chi.sire$SpawnType=="pooled"]<0.05)/nrow(chi.sire[chi.sire$SpawnType=="pooled",])
#0.8 #AKA 80% of pooled tanks deviated from expected contributions of sires




#### Ne ####

##### Expected Ne (Crow and Kimura 1973) ####
# Within each tank
Nd <- 3 #number of dams
Ns <- 3 #number of sires
NeC <- (4 * Nd * Ns)/(Nd + Ns)
NeC #6


##### Observed Ne (equalize family sizes) ####
#Use to detect reduction in Ne due to failed  (Crow and Kimura 1973)
Ne.eq <- data.frame(SpawnType=c(rep("factorial",10),rep("pooled",10)),rep=rep(1:10,2))
for (i in 1:nrow(Ne.eq)){
  Nd <- nrow(parents[parents$SpawnType==Ne.eq$SpawnType[i] 
                        & parents$rep==Ne.eq$rep[i] 
                        & parents$sex=="F"
                        & parents$prog>0,])
  Ns <- nrow(parents[parents$SpawnType==Ne.eq$SpawnType[i] 
                        & parents$rep==Ne.eq$rep[i] 
                        & parents$sex=="M"
                        & parents$prog>0,])
  Ne.eq$Nd[i] <- Nd
  Ne.eq$Ns[i] <- Ns
  Ne.eq$NeObs[i] <- (4 * Nd * Ns)/(Nd + Ns)
}
Ne.eq

## Reduction in Ne due to failed spawners
mean(Ne.eq$NeObs[Ne.eq$SpawnType=="pooled"]) #5.22
sd(Ne.eq$NeObs[Ne.eq$SpawnType=="pooled"]) #0.981835
(6-mean(Ne.eq$NeObs[Ne.eq$SpawnType=="pooled"]))/6 #13% average decline in Ne due to failed spawners
mean(Ne.eq$NeObs[Ne.eq$SpawnType=="factorial"]) #6
(6-mean(Ne.eq$NeObs[Ne.eq$SpawnType=="factorial"]))/6 #0% decline in Ne due to failed spawners
#Nd
mean(Ne.eq$Nd[Ne.eq$SpawnType=="pooled"]) #2.6
(3-mean(Ne.eq$Nd[Ne.eq$SpawnType=="pooled"]))/3 #13% average decline in Ned due to failed spawners
mean(Ne.eq$Nd[Ne.eq$SpawnType=="factorial"]) #3
(3-mean(Ne.eq$Nd[Ne.eq$SpawnType=="factorial"]))/3 #0% average decline in Ned due to failed spawners
#Ns
mean(Ne.eq$Ns[Ne.eq$SpawnType=="pooled"]) #2.8
(3-mean(Ne.eq$Ns[Ne.eq$SpawnType=="pooled"]))/3 #7% average decline in Ned due to failed spawners
mean(Ne.eq$Ns[Ne.eq$SpawnType=="factorial"]) #3
(3-mean(Ne.eq$Ns[Ne.eq$SpawnType=="factorial"]))/3 #0% average decline in Ned due to failed spawners


##### Observed Ne (account for family sizes) ####
#Use to detect reduction in Ne due to failed spawners plus variation in family sizes (Lacy 1989)

Ne.df <- data.frame(SpawnType=c(rep("factorial",10),rep("pooled",10)),rep=rep(1:10,2))
for (i in 1:nrow(Ne.df)){
  Ned <- 1/sum(parents$tank.prop[parents$sex=="F" & parents$SpawnType==Ne.df$SpawnType[i] & parents$rep==Ne.df$rep[i]]^2)
  Nes <- 1/sum(parents$tank.prop[parents$sex=="M" & parents$SpawnType==Ne.df$SpawnType[i] & parents$rep==Ne.df$rep[i]]^2)
  Ne.df$Ned[i] <- round(Ned,3)
  Ne.df$Nes[i] <- round(Nes,3)
  Ne.df$NeObs[i] <- round((4 * Ned * Nes)/(Ned + Nes),3)
}
Ne.df


## Summary stats for manuscript table 3
mean(Ne.df$NeObs[Ne.df$SpawnType=="factorial"]) #5.514
sd(Ne.df$NeObs[Ne.df$SpawnType=="factorial"])
mean(Ne.df$NeObs[Ne.df$SpawnType=="pooled"]) #3.8959
sd(Ne.df$NeObs[Ne.df$SpawnType=="pooled"])
mean(Ne.df$Ned[Ne.df$SpawnType=="factorial"]) #2.655
sd(Ne.df$Ned[Ne.df$SpawnType=="factorial"])
mean(Ne.df$Ned[Ne.df$SpawnType=="pooled"]) #2.0342
sd(Ne.df$Ned[Ne.df$SpawnType=="pooled"])
mean(Ne.df$Nes[Ne.df$SpawnType=="factorial"]) #2.8909
sd(Ne.df$Nes[Ne.df$SpawnType=="factorial"])
mean(Ne.df$Nes[Ne.df$SpawnType=="pooled"]) #2.0883
sd(Ne.df$Nes[Ne.df$SpawnType=="pooled"])

## Reduction in Ne due to variation in family size
(6-mean(Ne.df$NeObs[Ne.df$SpawnType=="pooled"]))/6 #35% average decline in Ne due to failed spawners & var in family size
#family size var alone accounted for 35-13 = 22% decline (remove effect of failed spawners)
(6-mean(Ne.df$NeObs[Ne.df$SpawnType=="factorial"]))/6 #8% average decline in Ne due to failed spawners & var in family size
#family size var alone account for 8-0 = 8% decline 
(3-mean(Ne.df$Ned[Ne.df$SpawnType=="pooled"]))/3 #32% average decline in Ned due to failed spawners & var in family size
#family size var alone accounted for 32-13 = 19% decline (remove effect of failed spawners)
(3-mean(Ne.df$Ned[Ne.df$SpawnType=="factorial"]))/3 #12% average decline in Ned due to failed spawners & var in family size
#family size var alone account for 12-0 = 12% decline 
(3-mean(Ne.df$Nes[Ne.df$SpawnType=="pooled"]))/3 #30% average decline in Nes due to failed spawners & var in family size
#family size var alone accounted for 30-7 = 23% decline (remove effect of failed spawners)
(3-mean(Ne.df$Nes[Ne.df$SpawnType=="factorial"]))/3 #4% average decline in Nes due to failed spawners & var in family size
#family size var alone account for 4-0 = 4% decline 



##### T-test for differences in Ne between methods ####

#Ne
t.test(Ne.df$NeObs[Ne.df$SpawnType=="pooled"],Ne.df$NeObs[Ne.df$SpawnType=="factorial"],alternative="two.sided")
# Welch Two Sample t-test
#
# data:  Ne.df$NeObs[Ne.df$SpawnType == "pooled"] and Ne.df$NeObs[Ne.df$SpawnType == "factorial"]
# t = -4.7442, df = 11.512, p-value = 0.0005341
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.3647244 -0.8714756
# sample estimates:
#   mean of x mean of y
# 3.8959    5.5140

#Ned
t.test(Ne.df$Ned[Ne.df$SpawnType=="pooled"],Ne.df$Ned[Ne.df$SpawnType=="factorial"],alternative="two.sided")
# Welch Two Sample t-test
#
# data:  Ne.df$Ned[Ne.df$SpawnType == "pooled"] and Ne.df$Ned[Ne.df$SpawnType == "factorial"]
# t = -2.5191, df = 12.431, p-value = 0.02636
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.15568305 -0.08591695
# sample estimates:
#   mean of x mean of y
# 2.0342    2.6550

#Nes
t.test(Ne.df$Nes[Ne.df$SpawnType=="pooled"],Ne.df$Nes[Ne.df$SpawnType=="factorial"],alternative="two.sided")
# Welch Two Sample t-test
#
# data:  Ne.df$Nes[Ne.df$SpawnType == "pooled"] and Ne.df$Nes[Ne.df$SpawnType == "factorial"]
# t = -3.8857, df = 9.613, p-value = 0.003264
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.265346 -0.339854
# sample estimates:
#   mean of x mean of y
# 2.0883    2.8909





#### Ne of release ####

##Add the estimated number of larvae to Ne.df
Ne.df <- cbind(Ne.df,Larvae=lar$Larvae) 
Ne.df


## Factorial design

#Expected Ne (Crow and Kimura 1973)
# Across all tanks
Nd <- 30
Ns <- 30
NeC <- (4 * Nd * Ns)/(Nd + Ns)
NeC #60

# NeR assuming equal contributions of each tank (Ryman and Laikre 1991)
#NeR <- 1/sum(xn^2/Nen) where 
# - xn is the proportion of offspring contributed by tank n to the release population 
#   (assuming equal contribution for now)
# - Nen is Ne of tank n
NeR.eq <- 1/sum((1/10)^2/Ne.df$NeObs[Ne.df$SpawnType=="factorial"]) #assuming equal contribution (1/10) of each tank to release pop
NeR.eq
#54.89054
NeR.eq/NeC
#0.9148423

# NeR based on estimated number of larvae per tank
NeR <- 1/sum((Ne.df$Larvae[Ne.df$SpawnType=="factorial"]/sum(Ne.df$Larvae[Ne.df$SpawnType=="factorial"]))^2/Ne.df$NeObs[Ne.df$SpawnType=="factorial"])
NeR
#52.81104
NeR/NeC
#0.8801839



## Pooled design

#Expected Ne (Crow and Kimura 1973)
# Across all tanks
Nd <- 30
Ns <- 30
NeC <- (4 * Nd * Ns)/(Nd + Ns)
NeC #60

# NeR assuming equal contributions of each tank (Ryman and Laikre 1991)
#NeR <- 1/sum(xn^2/Nen) where 
# - xn is the proportion of offspring contributed by tank n to the release population 
#   (assuming equal contribution for now)
# - Nen is Ne of tank n
NeR.eq <- 1/sum((1/10)^2/Ne.df$NeObs[Ne.df$SpawnType=="pooled"]) #assuming equal contribution (1/10) of each tank to release pop
NeR.eq
#36.46216
NeR.eq/NeC
#0.6077027


# NeR based on estimated number of larvae per tank
NeR <- 1/sum((Ne.df$Larvae[Ne.df$SpawnType=="pooled"]/sum(Ne.df$Larvae[Ne.df$SpawnType=="pooled"]))^2/Ne.df$NeObs[Ne.df$SpawnType=="pooled"])
NeR
#29.24533

NeR/NeC
#0.4874221


### Plots ####

#make spawning method a factor for easy plotting
tank$SpawnType <- factor(tank$SpawnType,levels=c("pooled","factorial"))
pc$SpawnType <- factor(pc$SpawnType,levels=c("pooled","factorial"))
parents$SpawnType <- factor(parents$SpawnType,levels=c("pooled","factorial"))
Ne.df$SpawnType <- factor(Ne.df$SpawnType,levels=c("pooled","factorial"))


##### Contributions of PCs, dams, sires to offspring ####

## Summary stats
median(pc$tank.prop[pc$SpawnType=="pooled"]) #0.04969475
sd(pc$tank.prop[pc$SpawnType=="pooled"]) #0.1472336
median(pc$tank.prop[pc$SpawnType=="factorial"]) #0.1011236
sd(pc$tank.prop[pc$SpawnType=="factorial"]) #0.05801825

median(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="F"]) #0.2872405
sd(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="F"]) #0.2826456
median(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="F"]) #0.3277565
sd(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="F"]) #0.1298704
median(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="M"]) #0.3186813
sd(parents$tank.prop[parents$SpawnType=="pooled" & parents$sex=="M"]) #0.2617509
median(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="M"]) #0.3240325
sd(parents$tank.prop[parents$SpawnType=="factorial" & parents$sex=="M"]) #0.06728139


## Boxplot - counts per PC by expt  (manuscript figure 2a)
boxplot(tank.prop ~ SpawnType,data=pc,ylim=c(0,1),xaxt="n",yaxt="n",
        border=rgb(120/255,120/255,120/255),
        main="",xlab="",ylab="")
axis(1,at=1,labels="pooled",tick=F,line=-0.5)
axis(1,at=2,labels="partial-\nfactorial",tick=F,line=0)
axis(2,at=seq(0,1,by=0.2),las=1,cex.axis=0.8)
title(ylab="Proportion of offspring per PC",line=2.5)
abline(h=1/9,lty=2)
legend("topright",legend="expected",lty=2,cex=0.75)

## Boxplot - counts per parent by expt (manuscript figure 2b)
#red = #DC3220, rgb(220/255,50/255,31/255,alpha=0.75)
#blue = #005AB5, rgb(0/255,90/255,181/255,alpha=0.75)

boxplot(tank.prop ~ sex + SpawnType,data=parents,xaxt="n",yaxt="n",
        xlab="",ylab="",
        col=c(rgb(220/255,50/255,31/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75),
              rgb(220/255,50/255,31/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75)),
        border=c("#DC3220","#005AB5","#DC3220","#005AB5"))
axis(1,at=c(1.5,3.5),labels=c("pooled","partial-factorial"),tick=F,line=-0.5)
axis(2,at=seq(0,1,by=0.2),las=1,cex.axis=0.8)
title(ylab="Proportion of offspring per parent",line=2.5)
abline(h=1/3,lty=2)
legend("topright",legend=c("dams","sires","expected"),pch=c(15,15,NA),lty=c(NA,NA,2),
       col=c(rgb(220/255,50/255,31/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75),"black"),
       cex=0.75)



## manuscript figure S1: prop offspring for each PC in pooled strategy 
par(mfrow=c(4,3))
for (i in 1:10){
  x <- pc[pc$SpawnType=="pooled" & pc$rep==i,]
  x <- x[order(x$tank.prop,decreasing=T),]
  barplot(x$tank.prop,ylim=c(0,1),border=F,
          main=paste("Pooled",i),ylab="Prop. offspring per PC",
          names.arg=c("PC 1","PC 2","PC 3","PC 4","PC 5","PC 6","PC 7","PC 8","PC 9"),las=2)
  title(xlab="Pair crosses",line=3.5)
  abline(h=1/9,lty=2)
}
par(mfrow=c(1,1))



## manuscript figure S2: prop offspring for each PC in factorial strategy 
par(mfrow=c(4,3))
for (i in 1:10){
  x <- pc[pc$SpawnType=="factorial" & pc$rep==i,]
  x <- x[order(x$tank.prop,decreasing=T),]
  barplot(x$tank.prop,ylim=c(0,1),border=F,
          main=paste("Factorial",i),ylab="Prop. offspring per PC",
          names.arg=c("PC 1","PC 2","PC 3","PC 4","PC 5","PC 6","PC 7","PC 8","PC 9"),las=2)
  title(xlab="Pair crosses",line=3.5)
  abline(h=1/9,lty=2)
}
par(mfrow=c(1,1))


#red = #DC3220, rgb(220/255,50/255,31/255,alpha=0.75)
#blue = #005AB5, rgb(0/255,90/255,181/255,alpha=0.75)

## manuscript figure S3: prop offspring for each parent in pooled strategy
par(mfrow=c(3,4))
for (i in 1:10){
  x <- parents[parents$SpawnType=="pooled" & parents$rep==i,]
  x <- x[order(x$tank.prop,decreasing=T),]
  x <- x[order(x$sex),]
  x$ParentID <- 
  barplot(x$tank.prop,border=F,ylim=c(0,1),ylab="Prop. offspring assigned",
          col=c(rgb(220/255,50/255,31/255,alpha=0.75),rgb(220/255,50/255,31/255,alpha=0.75),rgb(220/255,50/255,31/255,alpha=0.75),
                rgb(0/255,90/255,181/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75)),
          names.arg=substr(x$ParentID,1,4),las=2)
  title(main=paste("Pooled",i),adj=0)
  title(xlab="Parents",line=3.5)
  abline(h=1/3,lty=2)
}
par(mfrow=c(1,1))


## manuscript figure S4: prop offspring for each parent in factorial strategy
par(mfrow=c(3,4))
for (i in 1:10){
  x <- parents[parents$SpawnType=="factorial" & parents$rep==i,]
  x <- x[order(x$tank.prop,decreasing=T),]
  x <- x[order(x$sex),]
  barplot(x$tank.prop,ylim=c(0,1),border=F,ylab="Prop. offspring assigned",
          col=c(rgb(220/255,50/255,31/255,alpha=0.75),rgb(220/255,50/255,31/255,alpha=0.75),rgb(220/255,50/255,31/255,alpha=0.75),
                rgb(0/255,90/255,181/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75)),
          names.arg=substr(x$ParentID,1,4),las=2)
  title(main=paste("Factorial",i),adj=0)
  title(xlab="Parents",line=3.5)
  abline(h=1/3,lty=2)
}
par(mfrow=c(1,1))




##### Ne ####

#boxplot - overall Ne by expt (manuscript figure 3a)
boxplot(NeObs ~ SpawnType,data=Ne.df,xaxt="n",yaxt="n",ylim=c(0.22,6),
        border=rgb(120/255,120/255,120/255),
        main="",xlab="",ylab="")
        #main="",xlab="",ylab=expression("N"*italic(""[e])~"in replicates"))
axis(1,at=1,labels="pooled",tick=F,line=-0.5)
axis(1,at=2,labels="partial-\nfactorial",tick=F,line=0)
axis(2,at=0:6,labels=NA)
axis(2,at=c(0,2,4,6),las=1,tick=F)
title(ylab=expression(italic("N"[e])~"in multi-family crosses"),line=2.5)
abline(h=6,lty=2)
legend("bottomright",legend="expected",lty=2,cex=0.75)


##boxplot - Ned and Nes by expt (manuscript figure 3b)

#red = #DC3220, rgb(220/255,50/255,31/255,alpha=0.75)
#blue = #005AB5, rgb(0/255,90/255,181/255,alpha=0.75)

library(tidyr)
Ne.long <- gather(data=Ne.df,key=calc,value=value,Ned:NeObs)
Ne.long$sex[Ne.long$calc=="Ned"] <- "F"
Ne.long$sex[Ne.long$calc=="Nes"] <- "M"
#boxplot(value ~ sex + SpawnType,data=Ne.long)

boxplot(value ~ sex + SpawnType,data=Ne.long,xaxt="n",yaxt="n",ylim=c(0.11,3),
        xlab="",ylab="",
        col=c(rgb(220/255,50/255,31/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75),
              rgb(220/255,50/255,31/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75)),
        border=c("#DC3220","#005AB5","#DC3220","#005AB5"))
axis(1,at=c(1.5,3.5),labels=c("pooled","partial-factorial"),tick=F,line=-0.5)
axis(2,at=0:3,las=1)
title(ylab=expression(italic("N"[ed])*" and "*italic("N"[es])*~"in multi-family crosses   "),line=2.5)
abline(h=3,lty=2)
legend("bottomright",cex=0.75,legend=c(expression(italic("N"[ed])),expression(italic("N"[es])),"expected"),
       pch=c(15,15,NA),lty=c(NA,NA,2),
       col=c(rgb(220/255,50/255,31/255,alpha=0.75),rgb(0/255,90/255,181/255,alpha=0.75),"black"))




#### Compare 2021 vs 2022 pooled replicates ####

#Contributions of PCs to offspring
boxplot(tank.prop ~ expt,data=pc[pc$expt %in% 2:3,])
t.test(pc$tank.prop[pc$expt==2],pc$tank.prop[pc$expt==3],alternative="two.sided")
#NO SIG DIFF between pooled replicates in 2021 and 2022
#	Welch Two Sample t-test
#
#data:  pc$tank.prop[pc$expt == 2] and pc$tank.prop[pc$expt == 3]
#t = 2.552e-15, df = 68.919, p-value = 1
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.06509219  0.06509219
#sample estimates:
#  mean of x mean of y 
#0.1111111 0.1111111 

#Contributions of dams to offspring
boxplot(tank.prop~expt,data=parents[parents$sex=="F" & parents$expt %in% 2:3,])
t.test(parents$tank.prop[parents$sex=="F" & parents$expt==2],
       parents$tank.prop[parents$sex=="F" & parents$expt==3],alternative="two.sided")
#NO SIG DIFF in contributions of dams to offspring in 2021 vs 2022 pooled replicates
# Welch Two Sample t-test
# 
# data:  parents$tank.prop[parents$sex == "F" & parents$expt == 2] and parents$tank.prop[parents$sex == "F" & parents$expt == 3]
# t = 0, df = 21.495, p-value = 1
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2287808  0.2287808
# sample estimates:
#   mean of x mean of y 
# 0.3333333 0.3333333 


#Contributions of sires to offspring
boxplot(tank.prop~expt,data=parents[parents$sex=="M" & parents$expt %in% 2:3,])
t.test(parents$tank.prop[parents$sex=="M" & parents$expt==2],
       parents$tank.prop[parents$sex=="M" & parents$expt==3],alternative="two.sided")
#NO SIG DIFF in contributions of sires to offspring in 2021 vs 2022 pooled replicates
# Welch Two Sample t-test
# 
# data:  parents$tank.prop[parents$sex == "M" & parents$expt == 2] and parents$tank.prop[parents$sex == "M" & parents$expt == 3]
# t = 1.1066e-15, df = 22.867, p-value = 1
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2076169  0.2076169
# sample estimates:
#   mean of x mean of y 
# 0.3333333 0.3333333 


#Ne
boxplot(NeObs~expt,data=Ne.df[Ne.df$expt %in% 2:3,])
t.test(Ne.df$NeObs[Ne.df$expt==2],Ne.df$NeObs[Ne.df$expt==3],alternative="two.sided")
# NO SIG DIFF in Ne between 2021 and 2022 pooled replicates
# Welch Two Sample t-test
# 
# data:  Ne.df$NeObs[Ne.df$expt == 2] and Ne.df$NeObs[Ne.df$expt == 3]
# t = -0.082821, df = 4.934, p-value = 0.9372
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.991433  1.867600
# sample estimates:
#   mean of x mean of y 
# 3.858750  3.920667 



#Ned
boxplot(Ned~expt,data=Ne.df[Ne.df$expt %in% 2:3,])
t.test(Ne.df$Ned[Ne.df$expt==2],Ne.df$Ned[Ne.df$expt==3],alternative="two.sided")
# NO SIG DIFF in Ned between 2021 and 2022 pooled replicates
# Welch Two Sample t-test
# 
# data:  Ne.df$Ned[Ne.df$expt == 2] and Ne.df$Ned[Ne.df$expt == 3]
# t = -0.30055, df = 5.8798, p-value = 0.7741
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.380361  1.079694
# sample estimates:
#   mean of x mean of y 
# 1.944000  2.094333 


#Nes
boxplot(Nes~expt,data=Ne.df[Ne.df$expt %in% 2:3,])
t.test(Ne.df$Nes[Ne.df$expt==2],Ne.df$Nes[Ne.df$expt==3],alternative="two.sided")
# NO SIG DIFF in Nes between 2021 and 2022 pooled replicates
# Welch Two Sample t-test
# 
# data:  Ne.df$Nes[Ne.df$expt == 2] and Ne.df$Nes[Ne.df$expt == 3]
# t = 0.24645, df = 5.1044, p-value = 0.8149
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.08419  1.31569
# sample estimates:
#   mean of x mean of y 
# 2.15775   2.04200 
