############################################################################################################
###
### Operating model for the spatial DFA analysis
###
###
############################################################################################################

rm(list=ls())
gc()
cat("\014") # clear console# set working directory
library(ggplot2)
pacman::p_load(MASS, RandomFields, fields, geoR, kml, fpc, gtools, tweedie)
library("ggfortify")
library(ggplot2)
library(reshape2)
library(ggpubr)
# Packages
require("mgcv")
require("nFactors")

#### sourcing codesf
	#setwd("C:/Users/Kotou/Dropbox/Postdoc_projects/side_projects/spatial_DFA")

source("Generate_scenario_data22.r")

#### Scenario 2: targetting changes
### year 1:10 	= base price
### year 11:20	= target petrale, avoid rockfish and sablefish, same for dover
### year 21:30	= back to "normal"


	Sim_Settings2 <- list(
	"n_years" = 30,

	### Species parameters
		"n_species"  = 4,														# petrale, dover, darkblotched, sablefish
		"price_fish" = matrix(c(rep(c(1.5,1,1,0.5),10),								# important because controls vessel dynamics
								rep(c(0.5,1.5,1,1),10),
								rep(c(0.5,1,0.5,2),10)), ncol=4, byrow=T),

		## Fish habitat preference/movement control paramters
			"func_mvt_dist" = c("Exp", "Exp", "Exp", "Exp"),					# shape of mvt p_function
			"func_mvt_depth" = c("Lognorm", "Lognorm", "Lognorm", "Lognorm"),	# shape of depth preference p_function (choice of Normal, Exp, Lognormal, Uniform)
			"func_mvt_lat" = c("Unif", "Unif", "Unif", "Unif"),					# shape of range preference p_function (choice of Normal, Exp, Lognormal, Uniform)

			"Fish_dist_par1" = c(rep(0.4, 3), 1),								# mvt distance
			"Fish_dist_par2" = rep(0, 4), 										# mvt distance
			"Fish_depth_par1" = c(220, 550, 250, 1000),         				# depth preference
			"Fish_depth_par2" = c(0.5, 0.5, 0.3, 0.5),        					# depth preference
			"Fish_range_par1" = rep(0, 4),                    					# X-axis range
			"Fish_range_par2" = rep(50,4),                    					# X axis range

		## Species Schaefer pop dyn parameters
			"B0" =  c(320000,420000,56000,540000),								# petrale, dover, darkblotched, sablefish (mt)
			"r" = c(0.3, 0.17, 0.12, 0.08),
			"no.mvt" = FALSE,				# whether animals redistribute based on habitat preference

	### Parameters controlling the vessel dynamics
		"Nvessels" = 30,					# nb vessels
		"Tot_effort" = 500,					# end year nb of effort
		"CV_effort" = 0.2, 					# CV of effort around the yearly mean effort
		"qq_original" = 0.1,				# the average catchability coef for ALL vessels
		"CV_vessel" = 0.2,					# the between vessel variability in mean catchability
		"do.tweedie" = TRUE,				# include obsrevation error in catchability
		"xi" = 1.2,							# controls observation error in vessel catchability
		"phi" = 0.1,						# controls observation error in vessel catchability
		"Preference" = 1, 					# controls how effort concentrate into areas
		"Changing_preference" = FALSE, 		# whether effort concentration changes over time

	### Habitat (depth) parameters
		"Range_X" = 1:50,
		"Range_Y" = 1:50,
		"SD_O" = 150,						# SD of spatial observation
		"SpatialScale" = 4,					# scale of spatial variation
		"year.depth.restriction" = 1000,	# If MPA, then years when MPA implemented (default > n_years)
		"depth.restriction" = c(0,1),	# location of MPA

	### Other OM control features
		"plotting" = FALSE,
		"Interactive" = FALSE
	)										# End of simulation setting




#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Start simulation
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# Number of simulations
nsim = 1
nyrs = 30

#### Generate data
Data <- Generate_scenario_data(Sim_Settings2)
#save(Data,file="CPUEsim.rdata")

cat("\014") # clear console# set working directory

Data$Plots[[1]]
Data$Plots[[2]]
Data$Plots[[3]]
Data$Plots[[4]]+facet_wrap(~Sp,scale="free")

# Prepare data for 2-stage methods
# Flip Species Matrix
nsp = length(unique(Data$Data$Sp))
d =  data.frame(ID = rep(1:(nrow(Data$Data)/nsp),nsp),Data$Data)
d= reshape2::dcast(d,ID+year+vessel+X+Y+Xobs+Yobs+depth~Sp,value.var = "CPUE",fun.aggregate = mean)


sc <- scale_x_continuous(breaks=seq(min(d$year)+1, max(d$year), 2)) #, labels=as.character(vp)
th <- theme(axis.text.x  = element_text(angle=90, vjust=0.5), panel.grid.minor = element_blank())
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))




# Create dataset with species coulums
spd = d[,paste0("Sp",1:nsp)]
spd$sums = apply(spd,1,sum)
d = d[spd$sums>0,]
spd = spd[spd$sums>0,1:nsp]
# create data for PCA input
pca_dat <- sqrt(sqrt(spd/apply(spd,1,sum)))
#run PCA
names(pca_dat)
pca<-prcomp(pca_dat, scale=T)

pr_pca <- predict(pca,pca_dat)
pcs <- data.frame(pr_pca[,1:(ncol(pr_pca))])

# get Eigenvalue
eig = summary(pca)$sdev^2
# run Optimal Coordinates test
OCtest = nScree(eig)
npc = as.numeric(OCtest$Components[1])
ncl = npc+1

# clus_dat <- pca_dat
cl_sel = clara(pca_dat, ncl,sample=100,sampsize=250,pamLike = TRUE)

pcs$FT = as.matrix(cl_sel$cluster)[,1]
pca_dat$FT = factor(pcs$FT)
ppca= autoplot(pca,data = pca_dat,x=1,y=2,colour="FT",
             loadings = TRUE, loadings.colour = 'blue',
             loadings.label = TRUE,  loadings.label.size = 3)+theme_bw()
ppca

d = data.frame(d,pcs)

# ><> CHECK cluster composition here
sp.comp = aggregate(as.matrix(spd)~FT,d,sum)
sp.prop = round(sp.comp[,2:ncol(sp.comp)]/apply(sp.comp[,2:ncol(sp.comp)],1,sum)*100,1)
# Group composition
sp.prop
nspec = ncol(sp.prop)
tactics = 1:nrow(sp.prop)
col.sp = NULL
FT.sp = NULL
# Find the dominant species per cluster
d$Tactic = paste0("Tactic",d$FT)
ggplot(d,aes(year,fill=Tactic))+geom_bar(position="fill")+theme_bw()+sc+th




xyf  =reshape2::dcast(aggregate(FT~X+Y+Tactic,d,length),X+Y~Tactic,fun.aggregate = sum)
xyf$FT = factor(apply(xyf[,-c(1:2)],1,which.max))

ggplot() +
  geom_raster(data = xyf, aes(x = X, y = Y, fill = FT)) +
  scale_fill_discrete() +
  coord_quickmap(xlim = c(1,NA), ylim=c(1,NA),expand = FALSE)+ggtitle("Cluster")

xypc = reshape2::melt(aggregate(cbind(PC1,PC2)~X+Y,d,mean), id.vars = c("X", "Y"))
ggplot() +
  geom_raster(data = xypc, aes(x = X, y = Y, fill =value)) + facet_wrap(~variable)+
  scale_fill_gradientn(colors = (jet.colors(10))) +
  coord_quickmap(xlim = c(1,NA), ylim=c(1,NA),expand = FALSE)+ggtitle("Principle Components")






models = c("Base","Spatial","DPC","Cluster")
#nominal
mods = list(
  f1 = formula(y~as.factor(year)),
  f2 = formula(y~as.factor(year)+s(Xobs,Yobs,k=8)),
  # DPC GAM
  f3 = formula(y~as.factor(year)+s(Xobs,Yobs,k=8)+s(PC1,k=3)+s(PC2,k=3)),
  # Tactic
  f4 = formula(y~as.factor(year)+s(Xobs,Yobs,k=8)+as.factor(FT))
)


# Estimate disperstion first with GAMs
# fit tweedie GAMMs with random vessel effect

gams = popt = gamms =  list()
for(i in 1:nsp){
d$y = d[,paste0("Sp",i)]
  gams[[i]] = lapply(mods,function(x){
  gam(x,data=d,family=tw())
})
names(gams[[i]]) = models
popt[[i]] = do.call(c,lapply(gams[[i]],function(x){as.numeric(gsub(")","",strsplit(x$family[[1]], "=")[[1]][2]))}))
}
names(gams) = paste0("Sp",1:4)

# ><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Make prediction data for year, month
# ><>><>><>><>><>><>><>><>><>><>><>><>><>><>
pdat = data.frame(year=sort(unique(d$year)),Xobs=25,Yobs=25,FT=1,PC1=0,PC2=0)

idx = NULL

for(i in 1:nsp){
pr = do.call(rbind,Map(function(x,y){
  p <- predict(x,pdat,type="link",se=T)
  df = data.frame(year=pdat$year,Model=y,mu=exp(p$fit),se=p$se.fit)
  df$mu= df$mu/mean(df$mu) # normalize
  df$lci = exp(log(df$mu)-1.96*df$se)
  df$uci = exp(log(df$mu)+1.96*df$se)
  df$lci = df$lci/mean(df$mu)
  df$uci = df$uci/mean(df$mu)
  df
},x=gams[[i]],y=models))
idx = rbind(idx,data.frame(pr,Sp=paste0("Sp",i)))
}
rownames(idx) = 1:nrow(idx)

dt = Data$Biomass
dt$Model = "True"
ggplot(idx,aes(year,mu,color=Model))+facet_wrap(~Sp,scales = "free")+
  geom_ribbon(aes(ymin = lci, ymax = uci,fill=Model),alpha=0.2,col=0)+
  geom_line()+geom_point()+geom_hline(yintercept = 1,linetype=2)+sc+th+
  ylab("Normalized CPUE")+xlab("Year")+theme_bw()+geom_line(data=dt,aes(year,True),size=1)


# Select Best Model to include random vessel effect: Here DPC
smod = 3
# fit tweedie GAMMs with random vessel effect
gamms= NULL
for(i in 1:nsp){
  d$y = d[,paste0("Sp",i)]
  gamms[[i]] = Map(function(x,y){
    gamm(x,data=d,family=Tweedie(p=y),random=list(vessel=~1),method="REML")
  },x=mods[smod],y=popt[[i]][smod]) # change selmod
  names(gamms[[i]]) = models[smod]
}
names(gamms) = paste0("Sp",1:4)

# Final indices over all spatial cells


spdat = d
spdat$year=max(spdat$year)
spdat$xy = paste0(spdat$X,".",spdat$X)
spdat$Xobs=spdat$X
spdat$Yobs=spdat$Y
spdat = aggregate(cbind(PC1,PC2)~year+Xobs+Yobs+X+Y, spdat,mean)
spidx = NULL
for(i in 1:nsp){
  pr = do.call(rbind,Map(function(x,y){
    p <- predict(x$gam,spdat,type="link",se=T)
    df = data.frame(spdat,Model=y,CPUE=exp(p$fit)/mean(exp(p$fit)),se=p$se.fit)
    df
  },x=gamms[[i]],y=models[smod]))
  spidx = rbind(spidx,data.frame(pr,Sp=paste0("Sp",i)))
}
rownames(spidx) = 1:nrow(spidx)

# Spatial only
spoidx = NULL
for(i in 1:nsp){
  pr = do.call(rbind,Map(function(x,y){
    p <- predict(x,spdat,type="link",se=T)
    df = data.frame(spdat,Model=y,CPUE=exp(p$fit)/mean(exp(p$fit)),se=p$se.fit)
  },x=gams[[i]][2],y=models[2]))
  spoidx = rbind(spoidx,data.frame(pr,Sp=paste0("Sp",i)))
}
rownames(spoidx) = 1:nrow(spoidx)

prspo =  ggplot() +
  geom_raster(data = spoidx, aes(x = X, y = Y, fill = CPUE)) + facet_wrap(~Sp)+
  scale_fill_gradientn(colors = (jet.colors(10))) +
  coord_quickmap(xlim = c(0.5,NA), ylim=c(0.5,NA),expand = FALSE)+ggtitle("Predicted Spatial Distribution: Spatial Only")

prsp =  ggplot() +
  geom_raster(data = spidx, aes(x = X, y = Y, fill = log(CPUE+1))) + facet_wrap(~Sp)+
  scale_fill_gradientn(colors = (jet.colors(10))) +
  coord_quickmap(xlim = c(0.5,NA), ylim=c(0.5,NA),expand = FALSE)+ggtitle("Predicted Spatial Distribution: Spatial + DPC")


# Partial influence
yrs= unique(d$year)
pid = d
pid$year = max(pid$year)
pid$vessel = pid$vessel[1]
pid.xy = pid
pid.xy$PC1=0
pid.xy$PC2=0
pid.pc1 = pid
pid.pc1$PC2 = 0
pid.pc1$Xobs = median(pid$Xobs)
pid.pc1$Yobs = median(pid$Yobs)
pid.pc2 = pid
pid.pc2$PC1 = 0
pid.pc2$Xobs = median(pid$Xobs)
pid.pc2$Yobs = median(pid$Yobs)

deltas = NULL
for(i in 1:nsp){
  prid.xy = predict(gamms[[i]][[1]]$gam,pid.xy)
  prid.pc1 = predict(gamms[[i]][[1]]$gam,pid.pc1)
  prid.pc2 = predict(gamms[[i]][[1]]$gam,pid.pc2)
  p.xy = mean(prid.xy)
  p.pc1 = mean(prid.pc1)
  p.pc2 = mean(prid.pc2)

  delta =do.call(rbind,lapply(yrs,function(x){
  rbind(
  data.frame(year = x,Species=paste0("Sp",i),variable="Spatial",
  delta  = exp(sum(prid.xy[d$year==x]-p.xy)/length(prid.xy[d$year==x]))),
  data.frame(year = x,Species=paste0("Sp",i),variable="PC1",
             delta  = exp(sum(prid.pc1[d$year==x]-p.pc1)/length(prid.pc2[d$year==x]))),
  data.frame(year = x,Species=paste0("Sp",i),variable="PC2",
             delta  = exp(sum(prid.pc2[d$year==x]-p.pc2)/length(prid.pc2[d$year==x])))
  )
}))
deltas = rbind(deltas,delta)
}
deltas$variable = factor(deltas$variable,levels=unique(deltas$variable))

ggplot(deltas)+
  geom_line(aes(year,delta,color=variable),size=0.9)+
  facet_wrap(~Species)+geom_hline(yintercept = 1)+
  theme_bw()+sc+th+ylab("Influence")+xlab("Year")


# Extract Random Effect for DPC model
RE = NULL
for(i in 1:nsp){
re = as.numeric(summary(gamms[[i]][[1]]$lme)$coefficient$random$vessel)
re = data.frame(Vessel=1:length(re),Effect=re)
nB = nrow(re)
SD = data.frame(lci=c(-1*(sd(re$Effect))),uci=c(1*(sd(re$Effect))))
CI = data.frame(lsd=c(-1*(sd(re$Effect)*1.96)),usd=c(1*(sd(re$Effect)*1.96)))
re = data.frame(re,SD,CI,Sp=paste0("Sp",i))
RE = rbind(RE,re)
}

ggplot(RE,aes(Vessel,Effect))+facet_wrap(~Sp)+theme_bw()+geom_point()+ylab("Random Effect")+
  geom_ribbon(aes(ymin = lci, ymax = uci),alpha=0.4,col="lightgrey")+
  geom_ribbon(aes(ymin = lsd, ymax = usd),alpha=0.4,col="grey")+
  scale_x_continuous(breaks=c(0,seq(1,nB, 2)))+geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-0.2,0.2),linetype=2,color="red")

# Final standardized Index
stidx = NULL
for(i in 1:nsp){
  pr = do.call(rbind,Map(function(x,y){
    p <- predict(x$gam,pdat,type="link",se=T)
    df = data.frame(year=pdat$year,Model=y,mu=exp(p$fit),se=p$se.fit)
    df$mu= df$mu/mean(df$mu) # normalize
    df$lci = exp(log(df$mu)-1.96*df$se)
    df$uci = exp(log(df$mu)+1.96*df$se)
    df$lci = df$lci/mean(df$mu)
    df$uci = df$uci/mean(df$mu)
    df
  },x=gamms[[i]],y=models[smod]))
  stidx = rbind(stidx,data.frame(pr,Sp=paste0("Sp",i)))
}
rownames(idx) = 1:nrow(idx)

dt = Data$Biomass
dt$Model = "True"
ggplot(stidx,aes(year,mu,color=Model))+facet_wrap(~Sp,scales = "free")+
  geom_ribbon(aes(ymin = lci, ymax = uci,fill=Model),alpha=0.2,col=0)+
  geom_line()+geom_point()+geom_hline(yintercept = 1,linetype=2)+sc+th+
  ylab("Normalized CPUE")+xlab("Year")+theme_bw()+geom_line(data=dt,aes(year,True),size=1)

