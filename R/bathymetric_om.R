	                                                     
## Main function to generate data (the OM)
generate_bathymetric_scenario_data <- function(Sim_Settings){
		
	# Required packages are declared in DESCRIPTION.
	# fields supplies a circulant-embedding Gaussian-field simulator.
	required <- c("fields", "ggplot2", "reshape2", "tweedie")
	missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
	if(length(missing)) stop("Install these required packages first: ",
	                         paste(missing, collapse=", "), call.=FALSE)

  ## generate mvt matrix of an animal from a cell to the other 
	Mvt_matrix <- function(Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, x,y,data.bathym){
		Mvt_mat <- matrix(0,(length(x)*length(y)),(length(x)*length(y)))
		for (loc in 1:(length(x)*length(y))){ 
			Mvt_mat[loc,] <- Mvt_upd1(loc, Pop_param, Depth_eff=Depth_eff, Dist_eff=Dist_eff, Lat_eff=Lat_eff, Depth_func=Depth_func, Dist_func=Dist_func, Lat_func=Lat_func,x,y, data.bathym)
		}	
		return(Mvt_mat)
	}
		
  ## for updating the population distribution according to different different covariates (depth and/or dist), using joint prob
	Mvt_upd1 <- function(loc, Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, x,y, data.bathym){
		Di1 <- Pop_param[1]  	# 1st value for dist function
		Di2 <- Pop_param[2]  	# 2nd value for dist function
		De1 <- Pop_param[3]  	# 1st value for depth function
		De2 <- Pop_param[4]	 	# 2nd value for depth function
		Lat1 <- Pop_param[5]  	# 1st value for lat/long function
		Lat2 <- Pop_param[6]	# 2nd value for lat/long function
		
		a1 <- which(abs(data.bathym$x-data.bathym$x[loc])<=length(x)/2)
		a2 <- which(abs(data.bathym$x-data.bathym$x[loc])>length(x)/2)
		a3 <- which(abs(data.bathym$y-data.bathym$y[loc])<=length(y)/2)
		a4 <- which(abs(data.bathym$y-data.bathym$y[loc])>length(y)/2)
		Data.bathym <- list()
		Data.bathym$x[a1] <- data.bathym$x[a1]
		Data.bathym$x[a2] <- data.bathym$x[a2]#-length(x)
		Data.bathym$y[a3] <- data.bathym$y[a3]
		Data.bathym$y[a4] <- data.bathym$y[a4]#-length(y)
		Dist <- sqrt((Data.bathym$x-Data.bathym$x[loc])^2+(Data.bathym$y-Data.bathym$y[loc])^2)
		
		Prob1 <- 1
		if(Dist_eff=="T"){
			Prob1 <- dexp(Dist, rate=1/Di1)
		}
		
		Prob2 <- 1
		if(Depth_eff=="T"){
			if(Depth_func=="Norm") Prob2 <- dnorm(data.bathym$z, De1, De2*De1)
			if(Depth_func=="Exp") Prob2 <- dexp(data.bathym$z, rate=1/De1)
			if(Depth_func=="Lognorm") { SD <- sqrt(log(De2^2+1)); Prob2 <- dlnorm(data.bathym$z, (log(De1)-SD^2/2), SD) }
			if(Depth_func=="Unif") Prob2 <- dunif(data.bathym$z, De1, De2)
		}
		
		Prob3 <- 1
		if(Lat_eff=="T"){
			if(Lat_func=="Norm") Prob3 <- dnorm(data.bathym$x, Lat1, Lat2)
			if(Lat_func=="Exp") Prob3 <- dexp(data.bathym$x, rate=1/Lat1)
			if(Lat_func=="Lognorm") Prob3 <- dlnorm(data.bathym$x, (log(Lat1)-Lat2^2/2), Lat2)
			if(Lat_func=="Unif") Prob3 <- dunif(data.bathym$x, Lat1, Lat2)
		}
		
		Prob <- Prob1*Prob2*Prob3/sum(Prob1*Prob2*Prob3)
		return(Prob)
	}

  ## Calculate the initial population distribution
	Stable_pop_dist <- function(Mvt_mat_adult, B0){
		Pop1 <- rep(B0/(length(Range_X)*length(Range_Y)), (length(Range_X)*length(Range_Y)))
		Mvt_mat_adult1 <- t(Mvt_mat_adult)
		for (Time in seq_len(1000)){
			Pop <- Mvt_mat_adult1%*%Pop1
			change <- max(abs(Pop-Pop1))/max(1,max(abs(Pop1)))
			Pop1 <- Pop
			if(change < 1e-8) break
		}
		Initpop_adult_mat <- matrix(Pop1, nrow=length(Range_X), ncol=length(Range_Y))
		return(Initpop_adult_mat)
	}

	attach(Sim_Settings, name="sim_settings")
	    on.exit(detach("sim_settings"), add=TRUE)

	#### Set the bathymetry field
	# Same purpose as the former RandomFields::RMgauss/RFsimulate calls.
	# `fields::circulantEmbedding` is maintained and avoids the archived
	# RandomFields dependency.
	grid <- list(x=Range_X, y=Range_Y)
	cov_obj <- fields::circulantEmbeddingSetup(
	  grid, cov.args=list(Covariance="Exponential", aRange=SpatialScale, p=2))
	Bathym <- SD_O * fields::circulantEmbedding(cov_obj)
	data.bathym <- data.frame(ID = length(Range_X)*length(Range_Y), x = rep(Range_X, length(Range_Y))+0.5, y = rep(Range_Y, each=length(Range_X))+0.5, z=as.vector(Bathym))
	data.bathym$z = data.bathym$z+abs(min(data.bathym$z))	# to avoid getting negative depth here depth takes positive real number
	#data.bathym$Depth = -data.bathym$z
	jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	
	
	df = data.bathym
	df$Depth = -df$z
	
	pdepth = ggplot() +
	  geom_raster(data = df, aes(x = x, y = y, fill = Depth)) + 
	  scale_fill_gradientn(colors = (jet.colors(10))) + 
	  coord_quickmap(xlim = c(1,NA), ylim=c(1,NA),expand = FALSE)+ggtitle("Bathymetry")
	
	#### Set the pop mvt param
	Par_mvt_adult <- rbind(Fish_dist_par1, Fish_dist_par2, Fish_depth_par1, Fish_depth_par2, Fish_range_par1, Fish_range_par2);
	Mvt_mat_adult <- lapply(1:n_species, function(xxx) Mvt_matrix(Par_mvt_adult[,xxx], Depth_eff="T", Dist_eff="T", Lat_eff="T", Dist_func=func_mvt_dist[xxx], Depth_func=func_mvt_depth[xxx], Lat_func=func_mvt_lat[xxx], Range_X, Range_Y, data.bathym))

	#### Biomass dynamics: similarly to the work from Caruthers, we assume that there is a regional carrying capacity
	# Initial population distribution 
	B1 = B0*0.5*rlnorm(4,0,0.2)
	Pop_adult <- sapply(1:n_species, function(x) Stable_pop_dist(Mvt_mat_adult[[x]], B1[x]))				
	

	
	Bxy = data.frame(Pop_adult)
	colnames(Bxy) = paste0("Species",1:n_species)
	Bxy = cbind(df,Bxy)[,-c(1,4)]
	names(Bxy)
	
	Bxy = reshape2::melt(Bxy, id.vars = c("x", "y","Depth"))
	colnames(Bxy) = c(names(Bxy)[-5],"Density")
	#reshape2::dcast(Bxy1, x+ y+Depth ~ variable)               
	
	pbxy =  ggplot() +
	  geom_raster(data = Bxy, aes(x = x, y = y, fill = log(Density+1))) + facet_wrap(~variable)+
	  scale_fill_gradientn(colors = (jet.colors(10))) + 
	  coord_quickmap(xlim = c(0.5,NA), ylim=c(0.5,NA),expand = FALSE)+ggtitle("True spatial densisities")
	
	Biomass <- list()
	Biomass[[1]] <- Pop_adult
	Biomass[[1]] <- apply(Biomass[[1]], 2,function(x) replace(x, which(x<=0), 0))
		
	#### Vessel dynamics

	# generate total effort by year (total effort is randomly sampled with log normal error around a logistic type function)
	Sig_effort <- sqrt(log(CV_effort^2+1))
	Effort <- c(); Mean_Effort <- c()		
	Effort_area_year <- list()
	Sig_vessel <- sqrt(log(CV_vessel^2+1))
	qq_vessel <- rlnorm(Nvessels, log(qq_original)-Sig_vessel^2/2, Sig_vessel)		
		
	# use of tweedie distribution to generate catch data for individual boat
	I <- list()	
	xi <- xi	 # that controls for the value of power such that the variance is var[Y] = phi*mu^xi
	phi <- phi	  	
	Catch_area_year <- list(); Catch_area_year_ind <- c()
	set_id_counter <- 0L

	#### Updating the population and generate catch data according to the above specifications
	for (iyear in 1:n_years){
		
		### do vessel change their fishing preference over the year (this ONLY controls vessel concentration factor) 
		if(Changing_preference==TRUE) Preference <- 2-exp(-0.1*iyear)
		
		### Assigning total effort in a year		
		Mean_Effort[iyear] <- (Tot_effort/(1+exp(-0.1*iyear)))		# mean effort by year
		Effort[iyear] <- trunc(rlnorm(1, log(Mean_Effort[iyear])-Sig_effort^2/2, Sig_effort))		# effort by year
			
		### Determining the predicted area specific CPUE
		map.size <- length(Range_X)*length(Range_Y)
			
		### changing the catchability "qq" every year depending on whether a specific area is closed or not
		if (iyear>=year.depth.restriction){
			catchability <- rep(1,map.size)
			if(length(depth.restriction)==2) catchability <- replace(catchability, which(data.bathym$z >= depth.restriction[1] & data.bathym$z <= depth.restriction[2]), 0)
			if(is.character(depth.restriction) && length(depth.restriction)==1 && depth.restriction=="random") catchability <- replace(catchability, sample(1:map.size, size=min(20,map.size)), 0)
			qq <- catchability
		} else { qq <- qq_original }
		
		### Expected catch per cell uses biomass, catchability and species prices.
		catchability_cell <- if(length(qq)==1) rep(qq,map.size) else qq
		qmat <- matrix(1-exp(-catchability_cell),nrow=map.size,ncol=n_species)
		pmat <- matrix(price_fish[iyear,],nrow=map.size,ncol=n_species,byrow=TRUE)
		CPUE_tot <- rowSums(Biomass[[iyear]] * qmat * pmat)
		CPUE_tot[!is.finite(CPUE_tot) | CPUE_tot<0] <- 0
		if(sum(CPUE_tot^Preference)<=0) stop("No positive expected catch available for effort allocation.")
		
		### Distribute the total effort according to the above predicted CPUE
		p <- CPUE_tot^Preference/sum(CPUE_tot^Preference);min(p);max(p)		
		Effort_area_year[[iyear]] <- rmultinom(1, size=Effort[iyear], p=p)
			
		### Then calculate the total catch in each area  
		aa <- c(); AA <- c()			
		Catch_year_area_mat <- matrix(0, nrow=(length(Range_X)*length(Range_Y)), ncol=n_species)
		
		for (xyz in which(Effort_area_year[[iyear]]>0)){
			which.vessel.fished <- sample(1:Nvessels, Effort_area_year[[iyear]][xyz], replace=TRUE)
				
			nsets <- Effort_area_year[[iyear]][xyz]
			set_ids <- set_id_counter + seq_len(nsets)
			set_id_counter <- set_id_counter + nsets
			
			# This is the code using the continuous function (as suggested by Andre: more realistic)
			if(do.tweedie==TRUE) rand.vessel <- matrix(sapply(1:n_species, function(x) tweedie::rtweedie(nsets, xi=xi, mu=qq_vessel[which.vessel.fished], phi=phi)),ncol=n_species)
			if(do.tweedie==FALSE) rand.vessel <- matrix(sapply(1:n_species, function(x) qq_vessel[which.vessel.fished]),ncol=n_species)
			effective_q <- rand.vessel
			catch_area_vessel <- sapply(1:n_species, function(Y) { hum <- Biomass[[iyear]][xyz,Y]*(1-exp(-sum(effective_q[,Y]))); ifelse(hum==0, rep(0,nsets), hum*effective_q[,Y]/sum(effective_q[,Y]))})
			catch_area_vessel <- data.frame(
			  year=rep(iyear,nsets), SetID=set_ids,
			  X=rep(data.bathym$x[xyz],nsets), Y=rep(data.bathym$y[xyz],nsets),
			  depth=rep(data.bathym$z[xyz],nsets), vessel=which.vessel.fished,
			  as.data.frame(catch_area_vessel), check.names=FALSE)
			colnames(catch_area_vessel) <- c("year", "SetID", "X", "Y", "depth",
			  "vessel", paste0("Sp",1:n_species))
			Catch_year_area_mat[xyz,] <- colSums(catch_area_vessel[,paste0("Sp",1:n_species),drop=FALSE])
			
			AA <- rbind(AA, catch_area_vessel)
		}
		
		Catch_area_year_ind <- rbind(Catch_area_year_ind, AA)
		Catch_area_year[[iyear]] <- Catch_year_area_mat
		
   	### Then update the population
		Biomass[[iyear+1]] <- Biomass[[iyear]]+r*Biomass[[iyear]]*(1-Biomass[[iyear]]/Biomass[[1]])-Catch_area_year[[iyear]]
		Biomass[[iyear+1]] <- apply(Biomass[[iyear+1]],2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
		
		### population mvt at the end of year			
		if (no.mvt==FALSE)	Biomass[[iyear+1]] <- sapply(1:n_species, function(x) t(Mvt_mat_adult[[x]])%*%Biomass[[iyear+1]][,x])
		if (no.mvt==TRUE)	Biomass[[iyear+1]] <- Biomass[[iyear+1]]

		Biomass[[iyear+1]] <- apply(Biomass[[iyear+1]],2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))

	}
		
	    dhist = reshape2::melt(data.frame(Catch_area_year_ind[,paste0("Sp",1:n_species)]))
    ph = ggplot(data=dhist)+geom_histogram(aes(x=value),bins=80,fill="grey",col=1)+facet_wrap(~variable, scales = "free")+theme_bw()+xlab("CPUE")+ylab("Frequency")
    
		
	#### Return data                                                                # v.names="CPUE", 
	#Raw_data <- reshape(Catch_area_year_ind, varying=paste0("Sp",1:n_species), timevar="Species", times=1:n_species, direction="long") 
	Obs <- NULL
	for(i in 1:n_species) Obs = rbind(Obs, cbind(Catch_area_year_ind[,c("year","SetID","X","Y","depth","vessel")], "Species"=i, "CPUE"=Catch_area_year_ind[,paste0("Sp",i)]))
	Bio <- t(sapply(1:n_years, function(x) colSums(Biomass[[x]])))
	
	
	Obs = data.frame(Obs)
	Obs$Xobs = round(Obs$X/5)*5+5/2
	Obs$Yobs = round(Obs$Y/5)*5+5/2
	Obs$Sp = paste0("Sp",Obs$Species)	
	Obs = Obs[,c("year","SetID","X","Y","Xobs","Yobs","Species","Sp","CPUE","depth","vessel")]
	
	bio = t(t(Bio)/(apply(Bio,2,mean)))
	idx = reshape2::dcast(aggregate(CPUE~year+Sp,Obs,mean),year~Sp)
	t(t(idx[-1])/(apply(idx[-1],2,mean)))
	bio = data.frame(reshape2::melt(bio))
	names(bio) = c("year","Sp","True")
  
	idx = reshape2::dcast(aggregate(CPUE~year+Sp,Obs,mean),year~Sp)
  idx  = t(t(idx[-1])/(apply(idx[-1],2,mean)))
  idx = data.frame(reshape2::melt(idx))
  names(idx) = c("year","Sp","Nominal")
  idx = data.frame(idx,True=bio$True)
  
  pidx = ggplot(idx)+
    geom_line(aes(year,Nominal))+geom_point(aes(year,Nominal),pch=21,fill="grey")+
    theme_bw()+facet_wrap(~Sp)+geom_line(aes(year,True),linewidth=1.,col="darkblue")+
    geom_hline(yintercept = 1,linetype=2)
  
  plots = list()
  plots$Depth = pdepth
  plots$SpDens = pbxy
  plots$Hist = ph
  plots$CPUE = pidx+facet_wrap(~Sp,scale="free")  
  
  Return = list("Data"= Obs, "Biomass"=idx,"Biomass.XY"=Bxy,"Plots"=plots)
	return(Return)
}		
