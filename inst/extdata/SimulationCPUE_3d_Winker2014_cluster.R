#CODE CREATED BY Henning Winker
# July 2014
# Operating model after Winker et al. (2014)
# email: henning.winker@fao.org


#************************************************************************************
# Performance evaluation of four multispecies CPUE standardization procedures
#************************************************************************************

# R code for a simple operating model run to illustrate the performance of four alternative 
# multi-species CPUE standardization procedures: 
#(i) objective subsetting based on the P/A of non-target species after Stephens & McCall (2004) 
#(ii) clustering of fishing tactics based on catch composition information (e.g. Carvalho et al. 2010)
#(iii)?Direct Principal Component? with a GAM framework (DPC-GAM) 
# (iv) a simplified DPC-GLM application, for which the PC predictors are converted into binary vectors. 



#********************************
# 4 Habitat 10 Species Scenario #
#********************************

rm(list=ls())

require("mgcv")
require("nFactors")
require("scatterplot3d")

###### Habitats #####################################################

# 1 DeepReef
# 2 ShallowReef
# 3 DeepSediment
# 4 ShallowSediment

######################################################################
# Models
# DPC = model1; DPC.GLM =  model2; Stephens & McCall = model3 
# Clustering Fishing Tactics (Carvalho et al. 2010) = model4; 
# Stephens & McCall = model4 
###################################3###################################

# for model4, set number of clusters
# here for cluster given four habitats
nclus = 4


# Select between:
# control scenario = 0 (constant effort allocation) 
# test scenario = 1 (constant effort allocation)

scenario = 1

# Scenario set up 
# number of species
nspec <- 10
# number of obsersation per year
reps <- 500
# number of years
years <- 20
year <- c(1:years)

# set tweedie model parameters
p=1.3
disp=10

# determine the probability for a habitat switch during trip t

switch = 0.25
pSwitch <- 1-switch

# Set up data generation 

# Species order
# KOB, GLBK, HAKE, PANG,CRPN,SNTR,ROMN,DRGD,RSTM,SHRK
sp.names = c("KOB","GLBK","HAKE","PANG","CRPN","SNTR","ROMN","DRGD","RSTM","SHRK")

#random rate of population increase/decrease  (r)
r <- runif(nspec)*0.1*(-1)^round((runif(nspec)*4),0)

#random biomass trends for the start year
B1 <- c(200*exp(rnorm(3,0,0.5)),50*exp(rnorm(1,0,0.5)),200*exp(rnorm(1,0,0.5)),50*exp(rnorm(5,0,0.5)))

### get 'true' biomass trajectories for 20 years
B <- mat.or.vec(years,nspec)
for (i in 1:nspec)
{
B[,i] <- B1[i] * exp((year-1) * r[i])
}

# Biomass data
biomass_dat <- data.frame(year,B)

ylim = c(0,max(biomass_dat*1.05))

x11(width=7, height=10)
par(mfrow=c(2,1))
matplot(year, B,type="b",ylim=ylim, ylab="Biomass",pch=1:4,main="Biomass dynamics")


# Catchabilties for species i and habitat h
# KOB, GLBK, HAKE, PANG,CRPN,SNTR,ROMN,DRGD,RSTM,SHRK
q1 <- c(1.0,0.9,0.05,0.1,0.01,0.01,0.01,0.01,0.01,0.5)
q2 <- c(0.02,0.1,1.0,0.8,0.01,0.01,0.01,0.01,0.01,0.05)
q3 <- c(0.01,0.01,0.01,0.01,1.0,0.05,0.2,0.2,0.3,0.05)
q4 <- c(0.01,0.01,0.01,0.01,0.1,1.0,0.8,0.7,0.9,0.3)

# Effort error (lognormal)
error <- 0.2

#dataset structure
sim_year <- rep(c(1:20),each = reps)

KOB <- rep(B[,1], each = reps)
GLBK <- rep(B[,2], each = reps)
HAKE <- rep(B[,3], each = reps)
PANG <- rep(B[,4], each = reps)
CRPN <- rep(B[,5], each = reps)
SNTR <- rep(B[,6], each = reps)
ROMN <- rep(B[,7], each = reps)
DGRD <- rep(B[,8], each = reps)
RSTM <- rep(B[,9], each = reps)
SHRK <- rep(B[,10], each = reps)

n <- length(sim_year)
half <- c(1:10)
habitats <- c(1:4)

# create effort trends in proportion allocation per habitat and year

if(scenario==0)
{
vecP1 <- c(rep(0.2,20))*exp(rnorm(years,0,error))
vecP2 <- c(rep(0.25,20))*exp(rnorm(years,0,error))
vecP3 <- c(rep(0.3,20))*exp(rnorm(years,0,error))
vecP4 <- 1 - (vecP1+vecP2+vecP3)
}

if(scenario==1)
{
vecP1 <- c(0.45+0*(c(1:10)-1),0.11+0*(c(1:10)-1))*exp(rnorm(years,0,error))
vecP2 <- c(0.1+0*(c(1:10)-1),0.45+0*(c(1:10)-1))*exp(rnorm(years,0,error))
vecP3 <- (0.3-0.008*(year-1))*exp(rnorm(years,0,error))
vecP4 <- 1 - (vecP1+vecP2+vecP3)
}


#effort trends pi(h,y)
effort<- (rbind(vecP1,vecP2,vecP3,vecP4))
rand.effort <- effort[sample(1:4,replace=F),]

x11(height=4,width=8)
par(mfrow=c(1,1),mex=0.75, cex=1, cex.lab=1.1,mar=c(5,5,1.5,1))
plot(year,vecP1,type="b",pch=0, xlab="Year",ylab=expression(pi[h][","][y]),main="Effort dynamics",ylim=c(0,1))
lines(year,vecP2,type="b",pch=1,lty=2)
lines(year,vecP3,type="b",pch=1,lty=3)
lines(year,vecP4,type="b",pch=1,lty=4)


# assigne probabilities for habitat to trip 

prob1 <- rep((rand.effort[1,]), each=reps)
prob2 <- rep(rand.effort[2,]/(1-rand.effort[1,]), each=reps)
prob3 <- rep(rand.effort[3,]/(1-(rand.effort[1,]+rand.effort[2,])), each=reps)

# assign up to 2 habitats to each fishing trip
habitatA <- ifelse(runif(n)<prob1,1,ifelse(runif(n)<prob2,2,ifelse(runif(n)<prob3,3,4)))
habitatB <- ifelse(runif(n)<prob1,1,ifelse(runif(n)<prob2,2,ifelse(runif(n)<prob3,3,4)))

######################## Generate dataset ########################################################################

#Determine the fraction of effort (f(h,t))
propA <- 1/(1+exp(-(runif(n)-0.5)/0.5))
propB <- 1-propA

sim_dat0 <- data.frame(sim_year, habitatA, habitatB,propA,propB,KOB,GLBK,HAKE,PANG,CRPN,SNTR,ROMN,DGRD,RSTM,SHRK)

######## generate vectors with habitats 1 - 2

habA <- as.vector(sim_dat0$habitatA)
habB <- as.vector(sim_dat0$habitatB)

#######Generate CPUE records as a function of f(q,B,f)

Biomass <- sim_dat0[,6:(5+nspec)]
qA <- (if(habA[1] == 1) (rTweedie((q1*t(Biomass[1,])),p=p,disp)) else if(habA[1] == 2)(rTweedie((q2*t(Biomass[1,])),p=p,disp)) else if(habA[1] == 3)(rTweedie((q3*t(Biomass[1,])),p=p,disp)) else (rTweedie((q4*t(Biomass[1,])),p=p,disp)))
qB <- (if(runif(1)<pSwitch) (qA) else if(habB[1] == 1) (rTweedie((q1*t(Biomass[1,])),p=p,disp)) else if(habB[1] == 2)(rTweedie((q2*t(Biomass[1,])),p=p,disp)) else if(habB[1] == 3)(rTweedie((q3*t(Biomass[1,])),p=p,disp)) else (rTweedie((q4*t(Biomass[1,])),p=p,disp)))

for(i in 2:n)
{

rqA <- (if(habA[i] == 1) (rTweedie((q1*t(Biomass[i,])),p=p,disp)) else if(habA[i] == 2)(rTweedie((q2*t(Biomass[i,])),p=1.2,disp)) else if(habA[i] == 3)(rTweedie((q3*t(Biomass[i,])),p=1.2,disp)) else (rTweedie((q4*t(Biomass[i,])),p=p,disp)))
qA <- rbind(qA,rqA)
qB <- rbind(qB,if(runif(1)<pSwitch) (rqA) else if(habB[i] == 1) (rTweedie((q1*t(Biomass[i,])),p=p,disp)) else if(habB[i] == 2)(rTweedie((q2*t(Biomass[i,])),p=p,disp)) else if(habB[i] == 3)(rTweedie((q3*t(Biomass[i,])),p=p,disp)) else (rTweedie((q4*t(Biomass[i,])),p=p,disp)))
}

Bname <- data.frame(ifelse(Biomass>0,1,1))
cpueA <- propA*qA*Bname
cpueB <- propB*qB*Bname
cpue <- (cpueA+cpueB)



### create simulationed dataset 
dat <- data.frame(year=sim_year,habitatA,habitatB,cpue)

ncols <- 3+nspec
dat1 <- subset(dat,apply(dat[,4:ncols], 1,sum)>0)


##############################################################################
###########################  PCA  ############################################
##############################################################################

#extract catch composition
spec_dat <- dat1[,4:ncols]
spec_dat$sum <- apply(spec_dat, 1,sum)
mmult_dat <- subset(spec_dat,sum>0)

## prepare species composition for PCA
pca_dat <- (sqrt(sqrt(mmult_dat[,1:nspec]/mmult_dat$sum)))

#run PCA
pca<-prcomp(pca_dat, scale=T) 

# Predict PCA Loadings
pr_pca <- predict(pca,pca_dat)

pcs <- data.frame(pr_pca[,1:nspec])





# get Eigenvalue
eig = summary(pca)$sdev^2

# run Optimal Coordinates test
OCtest = nScree(eig)

x11()
plotnScree(OCtest)
# retain number of PCs in combination with Eigenvalue > 1

nPC = min(OCtest$Components$nkaiser,OCtest$Components$noc)


#### GLM-DPC variable as binary PCs

pcsBin = ifelse(pcs>0,1,0)


############################################################################
#################### Clustering fishing tactics ############################
##########################################################################

clus_dat <- mmult_dat[,1:nspec]/mmult_dat$sum

FT = kmeans(clus_dat,nclus)$cluster






###### Create datasets ######################################################


dataset <- cbind(dat1,pcs,FT,h1=qA[,1],h2=qB[,1])

dataset1 <- cbind(dat1,pcsBin,h1=qA[,1],h2=qB[,1]) # for DPC-GLM 


################## 3D plot of tactics ################################################################

hab.colors =c("blue","red","green4","cyan2","black")

hab.col = ifelse(dataset$h1==dataset$h2,hab.colors[habA],hab.colors[5])

tac.col = hab.colors[FT]

x11(width=8,height=4)
par(mfrow=c(1,2),mar=c(5,5,1,1))


scatterplot3d(pcs[,1],pcs[,2],pcs[,3],hab.col,pch=16,xlab="PC1", ylab="PC2", zlab="PC3")

scatterplot3d(pcs[,1],pcs[,2],pcs[,3],tac.col,pch=16,xlab="PC1", ylab="PC2", zlab="PC3")

hab.pch = ifelse(dataset$h1==dataset$h2,16,1)
 
x11(width=13,height=6)
par(mfrow=c(1,2),mar=c(5,3,1,1),mex=0.8)
scatterplot3d(pcs[,1],pcs[,2],pcs[,3],pch=21,bg="gray",xlab="PC1", ylab="PC2", zlab="PC3",angle=35,cex.symbols=0.8)
scatterplot3d(pcs[,1],pcs[,2],pcs[,3],hab.col,pch=hab.pch,xlab="PC1", ylab="PC2", zlab="PC3",angle=35,cex.symbols=0.8)


###########################################################################################
#################### Run standardization models for each species ##########################                
###########################################################################################



#  make CPUE histograms
br=50

x11(width=10, height=12)
par(mfrow=c(5,2),mex=0.65, cex=0.6, cex.lab=1,mar=c(5,5,1,1))
for (i in 1:nspec)
{
xlim <- c(0,0.5*max(cpue[,i]))
br <- round((max(cpue[,i])/log(max(cpue[,i])/20)),0)
hist(cpue[,i],breaks=br,col="gray",xlab="CPUE",main="",xlim=xlim)
legend('topright',paste(names(cpue[i])),bty="n") 
}




############################################################################
################ Stephens McCall method #####################################
###########################################################################

critV = 0

for(i in 1:nspec)

{


dataset$y <- dataset[,i+3]

# create binary response
p.sp = ifelse(dataset$y>0,1,0)   

# P/A of non-target species
PA.mat = ifelse(mmult_dat[,1:nspec]>0,1,0)

PA.mat = PA.mat[,-i]

### Dataset for Binomial model
Pdat = data.frame(obs = p.sp,pred= fitted.values(glm(p.sp ~ PA.mat,family="binomial")))



pred = Pdat$pred
obs = Pdat$obs
thresh=seq(0.01,1,by=0.01)
thresh.effect = thresh
thresh.count = thresh 



# search for p's between 0.1 and 0.8 for target species
for(s in 1:length(thresh))
{

  thresh.effect[s] = abs(sum(obs) - sum(pred > thresh[s]))
  thresh.count[s] = sum(pred > thresh[s]) 
}

thresh.dat =cbind(thresh,thresh.effect,thresh.count)
best = min(thresh.effect)
critV[i] = thresh[thresh.effect == best]

if(i ==1)
{
x11(width=9,height=6)
par(mex=0.8,cex=0.8)
if(i ==1)plot(thresh, thresh.effect, main='Subsetting (Stephens & McCall, 2004)', ylab = 'Difference between actual and predicted trips', xlab = 'Probablility of capture', type="l")
legend('bottomright',col=(1:10),paste(sp.names),lty=1,cex=1, bty="n")
}
if(i > 1) lines(thresh, thresh.effect,col=i) 
}



#***********************************************
#  RUN MULTISPECIES STANDARDIZATION PROCEFURES #
#***********************************************

#SPECIES LIST: KOB, GLBK, HAKE, PANG, CRPN, SNTR, ROMN, DRGD, RSTM, SHRK


# max number of note points for GAM splines
kn = 6
 
# run standardization models for each species

for(i in 1:nspec)

{

# true biomass for species i
B <- biomass_dat[,i+1]
# response vectors for species i
dataset$y <- dataset[,i+3]
dataset1$y <- dataset1[,i+3]

############################################################
# create objective subset of data (Stephens & McCall 2004)
############################################################


dataset$Prob = pred
StMc.dat = subset(dataset, Prob>critV[i])
# to avoid failures require a minimum of N = 500
if(nrow(StMc.dat)<500) StMc.dat = dataset

############################################################################

# fitting year-effect (Nominal CPUE)
m_base <- glm(y~as.factor(year),data=dataset,family=Tweedie(p=p))




# Original DPC-GAM
model1 <- gam(y~as.factor(year)+s(PC1,bs="cr",k=kn)+s(PC2,bs="cr",k=kn)+s(PC3,bs="cr",k=kn),data = dataset,family=Tweedie(p=p))

# DPC-GLM (wwith binary PCs)
model2  <- glm(y~as.factor(year)+PC1+PC2+PC3,data = dataset1,family=Tweedie(p=p))

# Clustering fishing tactics (Carvalho et al. 2010; 2014) 
model3  <- glm(y~as.factor(year)+as.factor(FT),data = dataset,family=Tweedie(p=p))

# Subsetting (Stephens & MacCall, 2004)
model4  <- glm(y~as.factor(year),data = StMc.dat,family=Tweedie(p=p))
  
  
# calc. normalized trend for known (true) biomass
trueB <- B/(mean(B))

## create prediction dataset 
sel_pcs<- apply(dataset[,(4+nspec):(7+nspec)],2,FUN=mean)

tn <-1
tactic <- tn
for( j in 2:20)
{
sel_pcs <- rbind(apply(dataset[,(4+nspec):(7+nspec)],2,FUN=mean))
}

pred_dat <- data.frame(year=year,sel_pcs,FT=1) 

# predict nominal CPUE
bc <- (predict(m_base,pred_dat,type="response"))
bc_mu <- as.matrix(bc)
norm.bc_mu <- bc_mu/mean(bc_mu)

# predict DPC-GAM CPUE
pm1 <- (predict(model1,pred_dat,type="response"))
pc_mu1 <- as.matrix(pm1)
norm.pc_mu1 <- pc_mu1/mean(pc_mu1)

# predict DPC-GLM CPUE
pm2 <- (predict(model2,pred_dat,type="response"))
pc_mu2 <- as.matrix(pm2)
norm.pc_mu2 <- pc_mu2/mean(pc_mu2)

# predict clustering FT CPUE
pm3 <- (predict(model3,pred_dat,type="response"))
pc_mu3 <- as.matrix(pm3)
norm.pc_mu3 <- pc_mu3/mean(pc_mu3)

# predict  objective subsetting
pm4 <- (predict(model4,pred_dat,type="response"))
pc_mu4 <- as.matrix(pm4)
norm.pc_mu4 <- pc_mu4/mean(pc_mu4)

# fit log-linear trends (rhat)to stardized CPUE
lm.bc <- lm(log(norm.bc_mu)~year)

lm.pc1 <- lm(log(norm.pc_mu1)~year)

lm.pc2 <- lm(log(norm.pc_mu2)~year)

lm.pc3 <- lm(log(norm.pc_mu3)~year)

lm.pc4 <- lm(log(norm.pc_mu4)~year)



# predict log-linear trends (slopes)
slope.bc <- exp(predict(lm.bc,data.frame(year)))
slope.pc1 <- exp(predict(lm.pc1,data.frame(year)))
slope.pc2 <- exp(predict(lm.pc2,data.frame(year)))
slope.pc3 <- exp(predict(lm.pc3,data.frame(year)))
slope.pc4 <- exp(predict(lm.pc4,data.frame(year)))


# Plot trends
if(i==1) 
{
x11(width=16., height=9)
par(mfrow=c(5,5),mex=0.65, cex=0.6, cex.lab=1,mar=c(5,5,1,1))
}


if(i==6) 
{
x11(width=16., height=9)
par(mfrow=c(5,5),mex=0.65, cex=0.6, cex.lab=1,mar=c(5,5,1,1))
}

ylim = c(0,3.4)
plot(year,trueB,ylim=ylim,col=c(2),lty=c(1),lwd=c(2),type=("l"), ylab="rel. Abundance",xlab="Year",main="")
points(year,norm.bc_mu,cex=1.1)
lines(year,slope.bc,lty=5,lwd=2)
text(2,2.8,paste(sp.names[i]))
text(18,2.8,"Nominal")
legend('topleft',c("True", "Pred. Slope","Year-effect"),lty = c(1,5,-1),pch=c(-1,-1,1),col=c(2,1,1),bty="n",pt.cex=1.2,cex=1,bg="white",horiz=T,lwd=c(2,2,1),x.intersp=0.5)



plot(year,trueB,ylim=ylim,col=c(2),lty=c(1),lwd=c(2),type=("l"), ylab="rel. Abundance",xlab="Year",main="")
points(year,norm.pc_mu1,cex=1.1)
lines(year,slope.pc1,lty=5,lwd=2)
text(2,2.8,paste(sp.names[i]))
text(18,2.8,paste("DPC GAM",sep=""))

plot(year,trueB,ylim=ylim,col=c(2),lty=c(1),lwd=c(2),type=("l"), ylab="rel. Abundance",xlab="Year",main="")
points(year,norm.pc_mu2,cex=1.1)
lines(year,slope.pc2,lty=5,lwd=2)
text(2,2.8,paste(sp.names[i]))
text(18,2.8,paste("DPC GLM",sep=""))

plot(year,trueB,ylim=ylim,col=c(2),lty=c(1),lwd=c(2),type=("l"), ylab="rel. Abundance",xlab="Year",main="")
points(year,norm.pc_mu3,cex=1.1)
lines(year,slope.pc3,lty=5,lwd=2)
text(2,2.8,paste(sp.names[i]))
text(18,2.8,paste("Cluster FT",sep=""))

plot(year,trueB,ylim=ylim,col=c(2),lty=c(1),lwd=c(2),type=("l"), ylab="rel. Abundance",xlab="Year",main="")
points(year,norm.pc_mu4,cex=1.1)
lines(year,slope.pc4,lty=5,lwd=2)
text(2,2.8,paste(sp.names[i]))
text(15,2.8,paste("Stephens & McCall",sep=""))

}








#***********************************************
#  RUN  3d visualisations #
#***********************************************

#SPECIES LIST: KOB, GLBK, HAKE, PANG, CRPN, SNTR, ROMN, DRGD, RSTM, SHRK


# max number of note points for GAM splines
kn = 6
 
# run standardization models for each species

for(i in 1:nspec)

{

# true biomass for species i
B <- biomass_dat[,i+1]
# response vectors for species i
dataset$y <- dataset[,i+3]
dataset1$y <- dataset1[,i+3]

############################################################
# create objective subset of data (Stephens & McCall 2004)
############################################################


dataset$Prob = pred
StMc.dat = subset(dataset, Prob>critV[i])
# to avoid failures require a minimum of N = 500
if(nrow(StMc.dat)<500) StMc.dat = dataset

############################################################################

# fitting year-effect (Nominal CPUE)
m_base <- glm(y~as.factor(year),data=dataset,family=Tweedie(p=p))




# Original DPC-GAM
model1 <- gam(y~as.factor(year)+s(PC1,bs="cr",k=kn)+s(PC2,bs="cr",k=kn)+s(PC3,bs="cr",k=kn),data = dataset,family=Tweedie(p=p))

# DPC-GLM (wwith binary PCs)
model2  <- glm(y~as.factor(year)+PC1+PC2+PC3,data = dataset1,family=Tweedie(p=p))

# Clustering fishing tactics (Carvalho et al. 2010; 2014) 
model3  <- glm(y~as.factor(year)+as.factor(FT),data = dataset,family=Tweedie(p=p))

################## 3D

par(mfrow=c(1,1),mar=c(5,5,1,1))

obs = dataset$y
DPC.fit = fitted.values(model1)
nom.fit = fitted.values(m_base) 
FT.fit = fitted.values(model3)
colors <- (topo.colors(100))
#colors <- rev(rainbow(100))

DPC.res = DPC.fit-obs
nom.res = nom.fit-obs

factor=300
obs.col <- colors[(obs - min(obs))/diff(range(obs))*factor + 1]
fit.col <- colors[(DPC.fit - min(DPC.fit))/diff(range(DPC.fit))*factor + 1]
nom.col <- colors[(nom.fit - min(nom.fit))/diff(range(nom.fit))*factor + 1]
ft.col <- colors[(FT.fit - min(FT.fit))/diff(range(FT.fit))*factor + 1]




x11(width=20,height=4)
par(mfrow=c(1,5),mar=c(5,3,1,1),mex=0.8)
scatterplot3d(pcs[,1],pcs[,2],pcs[,3],hab.col,pch=hab.pch,xlab="PC1", ylab="PC2", zlab="PC3",angle=35,cex.symbols=0.8,main="Habitats")
legend("topleft",paste(sp.names[i]),bty="n",cex=1.5,x.intersp=0.1)
scatterplot3d(pcs[,1],pcs[,2],pcs[,3],obs.col,pch=15,xlab="PC1", ylab="PC2", zlab="PC3",angle=35,cex.symbols=0.8,main="Observed")
scatterplot3d(pcs[,1],pcs[,2],pcs[,3],nom.col,pch=15,xlab="PC1", ylab="PC2", zlab="PC3",angle=35,cex.symbols=0.8,main="Nominal")
scatterplot3d(pcs[,1],pcs[,2],pcs[,3],fit.col,pch=15,xlab="PC1", ylab="PC2", zlab="PC3",angle=35,cex.symbols=0.8,main="DPC-predicted")
scatterplot3d(pcs[,1],pcs[,2],pcs[,3],ft.col,pch=15,xlab="PC1", ylab="PC2", zlab="PC3",angle=35,cex.symbols=0.8,main="Clustering-FT")

}



