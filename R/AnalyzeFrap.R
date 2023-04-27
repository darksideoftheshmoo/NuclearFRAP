
#' Calculate tau and nucFrap from kI and kEVfrac
#'
#' @param x a frap, frapPair or experimentSummary object
#' @return input object with modified param and PL.db elements
#' @export
transformTauNucFrac <- function(x){

	if(!all("param" %in% names(x))) stop("element param required")

	x$param <- transform(x$param, tau=1/((10^log10_kI) + (10^log10_kEVfrac)))
	x$param <- transform(x$param, nucFrac= (10^log10_kI) * tau)

	if("PL.db" %in% names(x)){
		if(!is.null(x$PL.db)){
			if(!"theta"%in%names(x$PL.db)) stop("PL.db should be from contour, not grid method")
			x$PL.db <- transform(x$PL.db, tau=1/((10^log10_kI) + (10^log10_kEVfrac)))
			x$PL.db <- transform(x$PL.db, nucFrac= (10^log10_kI) * tau)

			idVars <- intersect(c("pos","type","fit.index"),names(x$PL.db))
			CI.db<-
			ddply(x$PL.db[,c(idVars,"tau","nucFrac")],idVars,function(df){
				data.frame(tau.Upper=max(df$tau),tau.Lower=min(df$tau),nucFrac.Upper=max(df$nucFrac),nucFrac.Lower=min(df$nucFrac))
			})

			if("pos" %in% names(x$param)) x$param$pos <- as.factor(x$param$pos)
			if("type" %in% names(x$param)) x$param$type <- as.factor(x$param$type)
			x$param <- join(x$param,CI.db,by=idVars)
		}
	}

	return(x)
}


#' get simulated data from a experiment
#'
#' @param experiment a list containing 'frapPair' or 'frap' objects, as returned by \code{\link{analyzeFrapPair}} or \code{\link{analyzeFrap}}
#' @param bestFitOnly boolean indicating if only the parameters of the best fit for each cell should be used.
#' @param reRun boolean indicating if the simulated data should be re run. If FALSE the data from experiment is used.
#' @param param list with fixed values for model parameters. Used to overwrite the values of param for each cell.  
#' @return a data.frame
#' @export
experimentSimData <- function(experiment,bestFitOnly=TRUE,reRun=!is.null(param),param=NULL){
	#experiment <- frapExperiment
	if(reRun){
		if("type"%in%names(experiment[[1]]$db)){
			sim.db <- data.frame()
			for(ipos in names(experiment)){
				frap <- experiment[[ipos]]
				frapParam <- frap$param
				if(bestFitOnly) frapParam <- frap$param1
				if(!is.null(param)) frapParam <- do.call(transform,c(list(frapParam),param))
				sim.db <- rbind(sim.db, data.frame(pos=as.numeric(ipos),getSimDataPair(frap$db,frapParam,p=frap$p)))
			}
		} else {
			browser()
		}
	} else { #loading saved sim data in experiment object
		sim.rdb <- ldply(experiment,function(frap) frap$sim.rdb)
		if(bestFitOnly) sim.rdb <- subset(sim.rdb, fit.index==1)
		if("type"%in%names(sim.rdb)){	#Frap pairs
			sim.db<-cast(subset(sim.rdb,select=-sd),pos+fit.index+FRAP+time+type~variable)
		} else {
			browser()
			sim.db<-cast(subset(sim.rdb,select=-sd),pos+fit.index+FRAP+time~variable)
		}
	}
	return(sim.db)
}

#' get summary tables from a experiment
#'
#' @param experiment a list containing 'frapPair' or 'frap' objects, as returned by \code{\link{analyzeFrapPair}} or \code{\link{analyzeFrap}}
#' @return a object of class 'experimentSummary' with elements PL, param and param1
#' @export
experimentSummary <- function(experiment){

	PL.all<-ldply(experiment,function(l) l$PL.db)
	PL.all<-rename(PL.all,c(".id"="pos"))
	PL.all$pos<-factor(PL.all$pos)
	if(nrow(PL.all)==0) PL.all <- NULL

	param.all<-ldply(experiment,function(l) l$param)
	param.all<-rename(param.all,c(".id"="pos"))

	param1.all<-ldply(experiment,function(l) l$param1)
	param1.all<-rename(param1.all,c(".id"="pos"))
	if("type" %in% names(param1.all)){
		param1.all<-transform(param1.all,posType=paste0(pos,substr(type,1,1)))
		param.all<-transform(param.all,posType=paste0(pos,substr(type,1,1)))
	} else {
		param1.all<-transform(param1.all,posType=paste0(pos,"c"))
		param.all<-transform(param.all,posType=paste0(pos,"c"))
	}
	param1.all<-transform(param1.all,label=paste0(posType,".",fit.index))
	param.all<-transform(param.all,label=paste0(posType,".",fit.index))
	param1.all$pos<-factor(param1.all$pos)

	output <- list(PL.db=PL.all,param=param.all,param1=param1.all,p=experiment[[1]]$p)
	class(output) <- c("experimentSummary","list")
	return(output)
}

#' print method for experimentSummary
#'
#' @param x the experimentSummary object
#' @export
print.experimentSummary <- function(x,...) print(x$param)

#' create a scatter plot of kI and kEVfrac for all pairs
#'
#' @param frapPairList a list containing 'frapPair' objects, as returned by \code{\link{analyzeFrapPair}}
#' @export
plotExperimentSummary <- function(experiment,main=""){

	es <- experimentSummary(experiment)

	p5<-
	plotPL2D(es$PL.db,x="log10_kI",y="log10_kEVfrac",signifCost=qchisq(p=c(0.05), df=sum(!experiment$p$FIXED.PARAMETERS)),param=es$param1,color="pos",alpha=0.2) +
		labs(title=main)
	return(p5)
}

#' Analyze all positions of a Nuclear FRAP pair experiment
#'
#' @details this high level function calls \code{\link{analyzeFrapPair}} for all the pairs in the experiemnt
#'			the output is saved to a new directery named as the current experiment
#' @param d a cell.data object with the relevant experiment
#' @param profileLikelihood boolean indicanting if the profile likelihood for kI and kEVfrac will be run
#' @param allFitsResiduals boolean indicating if residuals plots should be made
#' @param montage boolean indicating if image montages should be made
#' @param PL.max.fit.index Profile Likelihood maximal fit index. Profile Likelihood will be run on fit.index <= PL.max.fit.index
#' @export
analyzeFrapPairExperiment<-function(d,profileLikelihood=TRUE,allFitsResiduals=TRUE,montage=TRUE,PL.max.fit.index=1){
	dir.create(d$experiment)
	#frapPair<-analyzeFrapPair(subset(d$data,pos==2),p=d$parameters,saveBasename=paste0(d$experiment,"/Analysis-"),profileLikelihood=TRUE)
	frapExperiment <- dlply(d$data,.(pos),analyzeFrapPair,p=d$parameters,saveBasename=paste0(d$experiment,"/Analysis-")
														 ,profileLikelihood=profileLikelihood,PL.max.fit.index=PL.max.fit.index
														 ,.progress = "text")
	save(frapExperiment,d,file=paste0(d$experiment,"/Analysis-",d$experiment,".RData"))

	#frapPair<-frapExperiment[[1]]
	l_ply(frapExperiment,function(frapPair){
		png(paste0(d$experiment,"/","Analysis-Pos",frapPair$pos,".png"),width=17,height=8,units="in",res=300)
		plot(frapPair)
		dev.off()

		if(allFitsResiduals){
			png(paste0(d$experiment,"/","AllFitsResiduals-Pos",frapPair$pos,".png"),width=16,height=4*nrow(frapPair$param),units="in",res=300)
			plotResidualsAllFits(frapPair)
			dev.off()
		}
	},.progress = "text")

	if(montage){
		for(i in unique(d$data$pos)){
			img<-try(imageCellPair(d,pos=i,display=FALSE))
			if(!(is.null(img)|(class(img)=="try-error"))) 
				writeImage(img,paste0(d$experiment,"/","Montage-Pos-",i,".jpg"))
		}
	}

	png(paste0(d$experiment,"/","kI-vs-kEVfrac-allCells.png"),width=8,height=7,units="in",res=300)
	print(plotExperimentSummary(frapExperiment,main=d$experiment))
	dev.off()

	return(invisible(frapExperiment))
}


#' Analyze all positions of a Nuclear FRAP experiment
#'
#' @details this high level function calls \code{\link{analyzeFrap}} for all the cells in the experiemnt
#'			the output is saved to a new directery named as the current experiment
#' @param d a cell.data object with the relevant experiment
#' @param profileLikelihood boolean indicanting if the profile likelihood for kI and kEVfrac will be run
#' @param allFitsResiduals boolean indicating if residuals plots should be made
#' @param montage boolean indicating if image montages should be made
#' @export
analyzeFrapExperiment<-function(d,profileLikelihood=TRUE,allFitsResiduals=TRUE,montage=TRUE){

	dir.create(d$experiment)
	frapExperiment <- dlply(d$data,.(pos),analyzeFrap,p=d$parameters,saveBasename=paste0(d$experiment,"/Analysis-"),profileLikelihood=profileLikelihood,.progress = "text")
	save(frapExperiment,d,file=paste0(d$experiment,"/Analysis-",d$experiment,".RData"))

	l_ply(frapExperiment,function(frap){
		while(!is.null(dev.list())) dev.off()

		png(paste0(d$experiment,"/","Analysis-Pos",frap$pos,".png"),width=17,height=8,units="in",res=300)
		plot(frap)
		dev.off()

		if(allFitsResiduals){
			png(paste0(d$experiment,"/","AllFitsResiduals-Pos",frap$pos,".png"),width=16,height=4*nrow(frap$param),units="in",res=300)
			plotResidualsAllFits(frap)
			dev.off()
		}
	},.progress = "text")

	if(montage){
		for(i in unique(d$data$pos)){
			img<-try(imageCell(d,pos=i,display=FALSE))
			if(!(is.null(img)|(class(img)=="try-error"))) 
				writeImage(img,paste0(d$experiment,"/","Montage-Pos-",i,".jpg"))
		}
	}

	png(paste0(d$experiment,"/","kI-vs-kEVfrac-allCells.png"),width=8,height=7,units="in",res=300)
	print(plotExperimentSummary(frapExperiment,main=d$experiment))
	dev.off()

	while(!is.null(dev.list())) dev.off()

	return(invisible(frapExperiment))
}

#' run a basic analysis of a Nuclear FRAP for a single cell
#'
#' @details this funcion runs a fit of the nuclear FRAP model to the data, a profile likelihood analysis of parameters
#'			kI and kEVfrac, and calculates the residuals of the model, for a single cell. 
#' @inheritParams analyzeFrapPair
#' @return a object of type 'frap' containing all the results of the analysis
#' @export
analyzeFrap<-function(db,p,nStarts=5,theta=seq(0,(2-1/8)*pi,pi/8),profileLikelihood=TRUE,saveBasename=NULL,overwrite=FALSE){

	pos<-unique(db$pos)
	if(length(pos)!=1) stop("only one pos can be passed to analyzeFrapPair")

	if(!is.null(saveBasename)){
		outputFilename <- paste0(saveBasename,"pos",pos,".RData")
		if(file.exists(outputFilename) & !overwrite){
			load(outputFilename)
			if("output" %in% ls()) return(output)
		}
	}

	param<-runModelFit(db,p=p)
	param<-param[!is.na(param$fit.index),]
	param1<-subset(param,fit.index==1)
	sim.rdb<-restructureSimData(getSimData(db,param,p=p))
	PL.db<-NULL
	if(profileLikelihood) 
		PL.db<-runProfileLikelihoodContour(db,param1,p=p,scan=c("log10_kEVfrac","log10_kI"),theta=theta,nStarts=nStarts
											,rescale=FALSE,boundaryCost=qchisq(p=c(0.05), df=sum(!p$FIXED.PARAMETERS)))

	#dealing with residuals
	res.db<-getResiduals(db,param=param1,p=p,melt=T)
	resPred.db<-ddply(res.db,.(variable),function(df){
		r<-range(df$value)
		v<-seq(r[1],r[2],length.out=100)
		data.frame(value=v,pred=length(df$value)*dnorm(v,mean=mean(df$value),sd=sd(df$value)))
	})
	resPred2.db<-subset(getResiduals(db,param=param1,p=p,melt=F),select=c("sd.f.nuc","sd.f.cyt"))
	resPred2.db<-rename(resPred2.db,c("sd.f.nuc"="res.nuc","sd.f.cyt"="res.cyt"))
	resPred2.db<-melt(resPred2.db,id.vars=NULL)
	predSD<-aggregateBy(resPred2.db,.(variable),FUN=mean)
	resPred2.db<-ddply(resPred2.db,.(variable),function(df){
		r1<-range(res.db[res.db$variable==unique(df$variable),"value"])
		v<-seq(r1[1],r1[2],length.out=100)
		data.frame(value=v,pred=length(df$value)*dnorm(v,mean=0,sd=mean(df$value)))
	})
	res2.db<-melt(getResiduals(db,param=param1,p=p,melt=F),measure.vars=c("res.nuc","res.cyt"))
	res2.db<-transform(res2.db,typeVar=variable)

	output<-list(db=db,param=param,param1=param1,sim.rdb=sim.rdb,PL.db=PL.db
				,res.db=res.db,res2.db=res2.db,resPred.db=resPred.db,resPred2.db=resPred2.db,predSD=predSD,pos=pos,p=p)
	class(output)<-"frap"
	if(!is.null(saveBasename)) save(output,file=outputFilename)
	return(output)
}



#' run a basic analysis of a Nuclear FRAP pair
#'
#' @details this funcion runs a fit of the nuclear FRAP model to the data, a profile likelihood analysis of parameters
#'			kI and kEVfrac, and calculates the residuals of the model 
#' @param db a data.frame with the data of a pair (required if expData=NULL). 'time', 'FRAP', 'type' and nuclear and cytoplasmic fluorescence variables required. 
#' @param p list of options for the FRAP experiment
#' @param profileLikelihood boolean indicanting if the profile likelihood for kI and kEVfrac will be run
#' @param saveBasename if not NULL, a charcater with the basename of the RData file where the results of the analysis will be saved. 
#' @param overwrite boolean indicating what to do if a .RData with the same name as the one to be created is found. 
#'		            If overwrite=FALSE, the .RData is loaded and its content returned by the function, therefore avoiding to re-run the analysis.
#'					If overwrite=TRUE, the analysis is redone and the file overwritten
#' @param nStarts number of random start points use by matlab's MultiStart function in order to get the global minimum
#' @param quadrantBisection number of times each quadrant is bisected to define the direction of the boundary points
#' @param PL.max.fit.index Profile Likelihood maximal fit index. Profile Likelihood will be run on fit.index <= PL.max.fit.index
#' @return a object of type 'frapPair' containing all the results of the analysis
#' @export
analyzeFrapPair<-function(db,p,profileLikelihood=TRUE,saveBasename=NULL,overwrite=FALSE,nStarts=5,quadrantBisections=1,PL.max.fit.index=1){
	pos<-unique(db$pos)
	if(length(pos)!=1) stop("only one pos can be passed to analyzeFrapPair")

	if(!is.null(saveBasename)){
		outputFilename <- paste0(saveBasename,"pos",pos,".RData")
		if(file.exists(outputFilename) & !overwrite){
			load(outputFilename)
			if("output" %in% ls()) return(output)
		}
	}

	param<-runModelFitPair(db,p=p)
	param<-param[!is.na(param$fit.index),]
	param1<-subset(param,fit.index<=PL.max.fit.index)
	sim.rdb<-restructureSimDataPair(getSimDataPair(db,param,p=p))
	PL.db<-NULL
	if(profileLikelihood) 
		PL.db<-runProfileLikelihoodContourPair(db,param1,p=p,scan=c("log10_kEVfrac","log10_kI"),quadrantBisections=quadrantBisections,nStarts=nStarts,precision=0.0001
				, boundaryCost=qchisq(p=c(0.05), df=sum(!p$FIXED.PARAMETERS)))
	#PL.db<-PL.db[PL.db$cost>0,]

	if(all(c("mother","daughter")%in%param1$type)){
		swapPoints<-param1
		swapPoints[swapPoints$type=="daughter",c("log10_kEVfrac","log10_kI")]<-param1[param1$type=="mother",c("log10_kEVfrac","log10_kI")]
		swapPoints[swapPoints$type=="mother",c("log10_kEVfrac","log10_kI")]<-param1[param1$type=="daughter",c("log10_kEVfrac","log10_kI")]
		swapPoints<-runProfileLikelihoodCostPair(db,swapPoints,p=p,nStart=10,fixed=c("log10_kEVfrac","log10_kI",names(which(p$FIXED.PARAMETERS==1))))
		swapPoints$fit.index<-1
		swapPointsSim<-restructureSimDataPair(getSimDataPair(db,swapPoints,p=p))
	} else {
		swapPoints<-NULL
		swapPointsSim<-NULL
	}

	#dealing with residuals
	res.db<-getResidualsPair(db,param=param1,p=p,melt=T)
	resPred.db<-ddply(res.db,.(type,variable),function(df){
		r<-range(df$value)
		v<-seq(r[1],r[2],length.out=100)
		data.frame(value=v,pred=length(df$value)*dnorm(v,mean=mean(df$value),sd=sd(df$value)))
	})
	resPred2.db<-subset(getResidualsPair(db,param=param1,p=p,melt=F),select=c("type","sd.f.nuc","sd.f.cyt"))
	resPred2.db<-rename(resPred2.db,c("sd.f.nuc"="res.nuc","sd.f.cyt"="res.cyt"))
	resPred2.db<-melt(resPred2.db,id.vars="type")
	predSD<-aggregateBy(resPred2.db,.(type,variable),FUN=mean)
	resPred2.db<-ddply(resPred2.db,.(type,variable),function(df){
		r1<-range(res.db[res.db$variable==unique(df$variable),"value"])
		v<-seq(r1[1],r1[2],length.out=100)
		data.frame(value=v,pred=length(df$value)*dnorm(v,mean=0,sd=mean(df$value)))
	})
	res2.db<-melt(getResidualsPair(db,param=param1,p=p,melt=F),measure.vars=c("res.nuc","res.cyt"))
	res2.db<-transform(res2.db,typeVar=interaction(type,variable,sep = " - "))

	output<-list(db=db,param=param,param1=param1,sim.rdb=sim.rdb,PL.db=PL.db,swapPoints=swapPoints,swapPointsSim=swapPointsSim
				,res.db=res.db,res2.db=res2.db,resPred.db=resPred.db,resPred2.db=resPred2.db,predSD=predSD,pos=pos,p=p)
	class(output)<-"frapPair"
	if(!is.null(saveBasename)) save(output,file=outputFilename)
	return(output)
}


#' plot a FRAP pair experiment
#'
#' @param frapPair a frapPair object, as returned by \code{\link{analyzeFrapPair}}
#' @param swapParam boolean indicating if swap fits parameters should be included in the parameters table
#' @export
plot.frapPair<-function(frapPair,swapParam=TRUE){
	p1 <- plotPL2D(frapPair$PL.db,x="log10_kI",y="log10_kEVfrac",signifCost=qchisq(p=c(0.05), df=sum(!frapPair$p$FIXED.PARAMETERS)),param=frapPair$param,tauBoundary=TRUE)+
		rcell2::zoom(x=c(-3,1),y=c(-3,1),expand.x=c(0,0.25),expand.y=c(0,0.25))+
		labs(title=paste0("pos ",frapPair$pos))

	p2<-
	cplot(frapPair$db,c(f.cyt.d,f.cyt.m,f.nuc.d,f.nuc.m)~time,geom=c("point"),group=interaction(FRAP,variable)
			,size=1,alpha=1,facets=~type, QC.filter = FALSE)+
	cplot(frapPair$db,c(f.cyt.d,f.cyt.m,f.nuc.d,f.nuc.m)~time,geom=c("line"),group=interaction(FRAP,variable)
			,size=0.5,alpha=0.25,facets=~type,layer=T, QC.filter = FALSE)+
	cplot(frapPair$sim.rdb,value~time,ymax=value+sd,ymin=value-sd,geom="smooth",stat="identity",group=interaction(FRAP,variable,fit.index)
			,layer=T,na.rm=FALSE,size=1,linetype=factor(fit.index),facets=~type,color=variable,alpha=0.2)+
	rcell2::yzoom(c(0,max(c(frapPair$db$f.nuc.d,frapPair$db$f.cyt.d,frapPair$db$f.nuc.m,frapPair$db$f.cyt.m))))
	if(!is.null(frapPair$swapPointsSim)){
		p2<-p2+
		cplot(frapPair$swapPointsSim,value~time,geom="line",group=interaction(FRAP,variable,fit.index)
				,layer=T,na.rm=FALSE,size=0.5,facets=~type,color="black",linetype=2)
	}
	p2 <- p2 + facet_wrap(~type,scales="free") + theme(legend.position="none") 
	#p2

	p3<-
	cplot(frapPair$res.db,~value,facets=~type+variable,binwidth=3,alpha=0.3) +
	cplot(frapPair$resPred.db,3*pred~value,geom="line",facets=~type+variable,layer=T,color=I(gg_color_hue(4)[2]),size=1,alpha=1)+
	cplot(frapPair$resPred2.db,3*pred~value,geom="line",facets=~type+variable,layer=T,color=I(gg_color_hue(4)[4]),size=1.5,alpha=1)+
	facet_wrap(~type+variable,nrow=1,scale="free") + theme(legend.position="none")
	#p3

	p4<-
	cplot(frapPair$res2.db,value~time,facets=~type+variable,color=typeVar)+
	cplot(frapPair$res2.db,value~time,geom="line",linetype=2,size=0.5,alpha=0.5
		,layer=T,facets=~type+variable,group=FRAP,color=typeVar)+
	geom_hline(aes(yintercept=value),linetype=3,data=frapPair$predSD)+
	geom_hline(aes(yintercept=-value),linetype=3,data=frapPair$predSD)+
	facet_wrap(~type+variable,nrow=1,scale="free_x") + 
	theme(legend.position="none") + scale_color_manual(values=gg_color_hue(4)[c(3,1,4,2)])
	#p4

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(6, 15)))

	print(p2, vp = vplayout(1:4, 1:8))
	print(p1, vp = vplayout(1:4, 9:12))
	print(p3, vp = vplayout(5:6, 1:8))
	print(p4, vp = vplayout(5:6, 9:15))

	paramTable <- transform(subset(frapPair$param,fit.index<=2),label=paste0(substr(type,1,1),fit.index))
	if(swapParam & !is.null(frapPair$swapPoints)) 
		paramTable <- rbind(paramTable, 
			rcell2:::conform(
				transform(
					subset(frapPair$swapPoints,index==1)
				,label=paste0("s",substr(type,1,1)))
			,to=paramTable)
		)
	paramTable <- transform(paramTable, tau=1/((10^log10_kI) + (10^log10_kEVfrac)))
	paramTable <- transform(paramTable, nucFrac= (10^log10_kI) * tau)
	paramTable <- transform(paramTable, intensityNucFrac= nucFrac * 2^log2_geomVfrac)
	if(all(c("log2_geomVfrac","log2_geomVfrac_mean","log2_geomVfrac_sd") %in% names(paramTable))) {
		#paramTable <- transform(paramTable, log2_geomVfrac_norm=(log2_geomVfrac-log2_geomVfrac_mean)/log2_geomVfrac_sd) 
		paramTable <- transform(paramTable, geomVfrac_sd=2^log2_geomVfrac_mean * log2_geomVfrac_sd) 
		paramTable <- transform(paramTable, geomVfrac_mean=2^log2_geomVfrac_mean, geomVfrac = 2^log2_geomVfrac) 
	}
	if(all(c("log10_kI","log10_kEVfrac") %in% names(paramTable))) paramTable <- transform(paramTable,kI = 10^log10_kI, kEVfrac = 10^log10_kEVfrac) 
	if(all(c("tot0","tot0_mean","tot0_sd") %in% names(paramTable))) paramTable <- transform(paramTable, tot0_norm=(tot0-tot0_mean)/tot0_sd)
	paramTable <- paramTable[,"log"!=substr(names(paramTable),1,3)]
	paramTable <- recast(subset(paramTable,select=c(-type,-fit.index)),variable~label,id.var="label")
	rownames(paramTable) <- paramTable$variable
	paramTable$variable <- NULL
	paramTable <- signif(paramTable,3)
	rmParam <- intersect(names(which(frapPair$p$FIXED.PARAMETERS==1)),names(which(apply(as.data.frame(paramTable),1,function(x)length(unique(x)))==1)))
	paramTable <-subset(paramTable,!rownames(paramTable)%in%rmParam)
	grid.table(as.data.frame(paramTable)
			,gpar.coretext = gpar(col = "black", cex = 0.6)
			,gpar.coltext = gpar(col = "black", cex = 0.8, fontface = "bold")
			,gpar.rowtext = gpar(col = "black", cex = 0.6, fontface = "italic")
			,vp=vplayout(1:4, 13:15))
	
}


#' plot a FRAP time course
#'
#' @param frap a frap object, as returned by \code{\link{analyzeFrap}}
#' @return a ggplot2 plot
#' @export
plotFrapTimeCourse <- function(frap){
	cplot(frap$db,c(f.cyt,f.nuc)~time,geom=c("point"),group=interaction(FRAP,variable),size=1,alpha=1, QC.filter = FALSE)+
	cplot(frap$db,c(f.cyt,f.nuc)~time,geom=c("line"),group=interaction(FRAP,variable),size=0.5,alpha=0.25,layer=T, QC.filter = FALSE)+
	cplot(frap$sim.rdb,value~time,ymax=value+sd,ymin=value-sd,geom="smooth",stat="identity",group=interaction(FRAP,variable,fit.index)
		,layer=T,na.rm=FALSE,size=1,linetype=factor(fit.index),color=variable,alpha=0.2)+
	rcell2::yzoom(c(0,max(c(frap$db$f.nuc,frap$db$f.cyt))))
}

#' get the parameters table for a FRAP fit
#'
#' @param frap a frap object, as returned by \code{\link{analyzeFrap}}
#' @return a data.frame with the parameters of each fit
#' @export
getParamTable <- function(frap){
	paramTable <- transform(subset(frap$param,fit.index<=2),label=paste0("c",fit.index))
	paramTable <- transform(paramTable, tau=1/((10^log10_kI) + (10^log10_kEVfrac)))
	paramTable <- transform(paramTable, nucFrac= (10^log10_kI) * tau)
	paramTable <- transform(paramTable, intensityNucFrac= nucFrac * 2^log2_geomVfrac)
	if(all(c("log2_geomVfrac","log2_geomVfrac_mean","log2_geomVfrac_sd") %in% names(paramTable))) {
		#paramTable <- transform(paramTable, log2_geomVfrac_norm=(log2_geomVfrac-log2_geomVfrac_mean)/log2_geomVfrac_sd) 
		paramTable <- transform(paramTable, geomVfrac_sd=2^log2_geomVfrac_mean * log2_geomVfrac_sd) 
		paramTable <- transform(paramTable, geomVfrac_mean=2^log2_geomVfrac_mean, geomVfrac = 2^log2_geomVfrac) 
	}
	if(all(c("log10_kI","log10_kEVfrac") %in% names(paramTable))) paramTable <- transform(paramTable,kI = 10^log10_kI, kEVfrac = 10^log10_kEVfrac) 
	if(all(c("tot0","tot0_mean","tot0_sd") %in% names(paramTable))) paramTable <- transform(paramTable, tot0_norm=(tot0-tot0_mean)/tot0_sd)
	paramTable <- paramTable[,"log"!=substr(names(paramTable),1,3)]
	paramTable <- recast(subset(paramTable,select=c(-fit.index)),variable~label,id.var="label")
	rownames(paramTable) <- paramTable$variable
	paramTable$variable <- NULL
	paramTable <- signif(paramTable,3)
	rmParam <- intersect(names(which(frap$p$FIXED.PARAMETERS==1)),names(which(!is.na(frap$p$VALUE.PARAMETERS))))
	#paramTable <-subset(paramTable,!rownames(paramTable)%in%rmParam)
	return(as.data.frame(paramTable))
}

#' plot a FRAP experiment
#'
#' @param frap a frap object, as returned by \code{\link{analyzeFrap}}
#' @export
plot.frap<-function(frap){

	p1 <- plotPL2D(frap$PL.db,x="log10_kI",y="log10_kEVfrac",signifCost=qchisq(p=c(0.05), df=sum(!frap$p$FIXED.PARAMETERS)),param=frap$param,color="fit.index")+
		rcell2::zoom(x=c(-3,1),y=c(-3,1),expand.x=c(0,0.25),expand.y=c(0,0.25))+
		labs(title=paste0("pos ",frap$pos))

	p2 <- plotFrapTimeCourse(frap)
	p2 <- p2 + theme(legend.position="none") 
	#p2

	p3<-
	cplot(frap$res.db,~value,facets=~variable,binwidth=5,alpha=0.3) +
	cplot(frap$resPred.db,5*pred~value,geom="line",facets=~variable,layer=T,color=I(gg_color_hue(4)[2]),size=1,alpha=1)+
	cplot(frap$resPred2.db,5*pred~value,geom="line",facets=~variable,layer=T,color=I(gg_color_hue(4)[4]),size=1.5,alpha=1)+
	facet_wrap(~variable,nrow=1,scale="free") + theme(legend.position="none")
	#p3

	p4<-
	cplot(frap$res2.db,value~time,facets=~variable,color=typeVar)+
	cplot(frap$res2.db,value~time,geom="line",linetype=2,size=0.5,alpha=0.5
		,layer=T,facets=~variable,group=FRAP,color=typeVar)+
	geom_hline(aes(yintercept=value),linetype=3,data=frap$predSD)+
	geom_hline(aes(yintercept=-value),linetype=3,data=frap$predSD)+
	facet_wrap(~variable,nrow=1,scale="free_x") + 
	theme(legend.position="none") + scale_color_manual(values=gg_color_hue(4)[c(1,3)])
	#p4

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(6, 15)))

	print(p2, vp = vplayout(1:4, 1:6))
	print(p1, vp = vplayout(1:4, 7:12))
	print(p3, vp = vplayout(5:6, 1:6))
	print(p4, vp = vplayout(5:6, 7:12))

	paramTable <- transform(subset(frap$param,fit.index<=2),label=paste0("c",fit.index))
	paramTable <- transform(paramTable, tau=1/((10^log10_kI) + (10^log10_kEVfrac)))
	paramTable <- transform(paramTable, nucFrac= (10^log10_kI) * tau)
	paramTable <- transform(paramTable, intensityNucFrac= nucFrac * 2^log2_geomVfrac)
	if(all(c("log2_geomVfrac","log2_geomVfrac_mean","log2_geomVfrac_sd") %in% names(paramTable))) {
		#paramTable <- transform(paramTable, log2_geomVfrac_norm=(log2_geomVfrac-log2_geomVfrac_mean)/log2_geomVfrac_sd) 
		paramTable <- transform(paramTable, geomVfrac_sd=2^log2_geomVfrac_mean * log2_geomVfrac_sd) 
		paramTable <- transform(paramTable, geomVfrac_mean=2^log2_geomVfrac_mean, geomVfrac = 2^log2_geomVfrac) 
	}
	if(all(c("log10_kI","log10_kEVfrac") %in% names(paramTable))) paramTable <- transform(paramTable,kI = 10^log10_kI, kEVfrac = 10^log10_kEVfrac) 
	if(all(c("tot0","tot0_mean","tot0_sd") %in% names(paramTable))) paramTable <- transform(paramTable, tot0_norm=(tot0-tot0_mean)/tot0_sd)
	paramTable <- paramTable[,"log"!=substr(names(paramTable),1,3)]
	paramTable <- recast(subset(paramTable,select=c(-fit.index)),variable~label,id.var="label")
	rownames(paramTable) <- paramTable$variable
	paramTable$variable <- NULL
	paramTable <- signif(paramTable,3)
	rmParam <- intersect(names(which(frap$p$FIXED.PARAMETERS==1)),names(which(!is.na(frap$p$VALUE.PARAMETERS))))
	#paramTable <-subset(paramTable,!rownames(paramTable)%in%rmParam)
	grid.table(as.data.frame(paramTable)
			,gpar.coretext = gpar(col = "black", cex = 0.6)
			,gpar.coltext = gpar(col = "black", cex = 0.8, fontface = "bold")
			,gpar.rowtext = gpar(col = "black", cex = 0.6, fontface = "italic")
			,vp=vplayout(1:6, 13:15))
	
}


#' plot a FRAP experiment, with best lower and upper estimations for a parameter
#'
#' @param frap a frap object, as returned by \code{\link{analyzeFrap}}, with best estimation of parameters
#' @param upper a frap object, as returned by \code{\link{analyzeFrap}}, with upper boundary for a parameter
#' @param lower a frap object, as returned by \code{\link{analyzeFrap}}, with lower boundary for a parameter
#' @export
plotBestUpperLowerFrap<-function(frap,upper,lower){

	frap$PL.db$estimateType <- "best" 
	upper$PL.db$estimateType <- "upper" 
	lower$PL.db$estimateType <- "lower"
	PL.db <- rbind(frap$PL.db, upper$PL.db, lower$PL.db)

	frap$param$estimateType <- "best" 
	upper$param$estimateType <- "upper" 
	lower$param$estimateType <- "lower"
	param <- rbind(frap$param, upper$param, lower$param)

	frap$sim.rdb$estimateType <- "best" 
	upper$sim.rdb$estimateType <- "upper" 
	lower$sim.rdb$estimateType <- "lower"
	sim.rdb <- rbind(frap$sim.rdb, upper$sim.rdb, lower$sim.rdb)

	p1 <- plotPL2D(PL.db,x="log10_kI",y="log10_kEVfrac",signifCost=qchisq(p=c(0.05), df=sum(!frap$p$FIXED.PARAMETERS)),param=param,color="estimateType")+
		rcell2::zoom(x=c(-3,1),y=c(-3,1),expand.x=c(0,0.25),expand.y=c(0,0.25))+
		labs(title=paste0("pos ",frap$pos))

	p2<-
	cplot(frap$db,c(f.cyt,f.nuc)~time,geom=c("point"),group=interaction(FRAP,variable)
			,size=1,alpha=1, QC.filter = FALSE, color="black")+
	cplot(frap$db,c(f.cyt,f.nuc)~time,geom=c("line"),group=interaction(FRAP,variable)
			,size=0.5,alpha=0.25,layer=T, QC.filter = FALSE, color="black")+
	cplot(sim.rdb
		,value~time,ymax=value+sd,ymin=value-sd,geom="smooth",stat="identity",group=interaction(FRAP,variable,fit.index,estimateType)
		,layer=T,na.rm=FALSE,size=1,linetype=factor(fit.index),color=estimateType,alpha=0.2)+
	rcell2::yzoom(c(0,max(c(frap$db$f.nuc,frap$db$f.cyt))))
	p2 <- p2 + theme(legend.position="none") 

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(1, 3)))

	print(p2, vp = vplayout(1, 1))
	print(p1, vp = vplayout(1, 2))

	paramTable <- transform(subset(param,fit.index<=2),label=paste0(substr(estimateType,1,1),fit.index))
	paramTable <- recast(subset(paramTable,select=c(-fit.index,-estimateType)),variable~label,id.var="label")
	rownames(paramTable) <- paramTable$variable
	paramTable$variable <- NULL
	paramTable <- signif(paramTable,3)
	grid.table(as.data.frame(paramTable)
			,gpar.coretext = gpar(col = "black", cex = 0.6)
			,gpar.coltext = gpar(col = "black", cex = 0.8, fontface = "bold")
			,gpar.rowtext = gpar(col = "black", cex = 0.6, fontface = "italic")
			,vp=vplayout(1, 3))
	
}


#' plot a FRAP pair experiment, with best lower and upper estimations for a parameter
#'
#' @param frap a frapPair object, as returned by \code{\link{analyzeFrapPair}}, with best estimation of parameters
#' @param upper a frapPair object, as returned by \code{\link{analyzeFrapPair}}, with upper boundary for a parameter
#' @param lower a frapPair object, as returned by \code{\link{analyzeFrapPair}}, with lower boundary for a parameter
#' @export
plotBestUpperLowerFrapPair <- function(frapPair,upper,lower){

	frapPair$PL.db$estimateType <- "best" 
	upper$PL.db$estimateType <- "upper" 
	lower$PL.db$estimateType <- "lower"
	PL.db <- rbind(frapPair$PL.db, upper$PL.db, lower$PL.db)

	frapPair$param$estimateType <- "best" 
	upper$param$estimateType <- "upper" 
	lower$param$estimateType <- "lower"
	param <- rbind(frapPair$param, upper$param, lower$param)

	frapPair$sim.rdb$estimateType <- "best" 
	upper$sim.rdb$estimateType <- "upper" 
	lower$sim.rdb$estimateType <- "lower"
	sim.rdb <- rbind(frapPair$sim.rdb, upper$sim.rdb, lower$sim.rdb)


	p1 <- plotPL2D(PL.db,x="log10_kI",y="log10_kEVfrac",signifCost=qchisq(p=c(0.05), df=sum(!frapPair$p$FIXED.PARAMETERS)),param=param,tauBoundary=TRUE)+
		rcell2::zoom(x=c(-3,1),y=c(-3,1),expand.x=c(0,0.25),expand.y=c(0,0.25))+
		labs(title=paste0("pos ",frapPair$pos))

	p2<-
	cplot(frapPair$db,c(f.cyt.d,f.cyt.m,f.nuc.d,f.nuc.m)~time,geom=c("point"),group=interaction(FRAP,variable)
			,size=1,alpha=1,facets=~type, QC.filter = FALSE)+
	cplot(frapPair$db,c(f.cyt.d,f.cyt.m,f.nuc.d,f.nuc.m)~time,geom=c("line"),group=interaction(FRAP,variable)
			,size=0.5,alpha=0.25,facets=~type,layer=T, QC.filter = FALSE)+
	cplot(sim.rdb,value~time,fit.index==1,ymax=value+sd,ymin=value-sd,geom="smooth",stat="identity"
			,group=interaction(FRAP,variable,fit.index,estimateType)
			,layer=T,na.rm=FALSE,size=1,linetype=factor(estimateType),facets=~type,color=variable,alpha=0.2)+
	rcell2::yzoom(c(0,max(c(frapPair$db$f.nuc.d,frapPair$db$f.cyt.d,frapPair$db$f.nuc.m,frapPair$db$f.cyt.m))))
	p2 <- p2 + facet_wrap(~type,scales="free") + theme(legend.position="none") 
	#p2

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(1, 4)))

	print(p2, vp = vplayout(1, 1:2))
	print(p1, vp = vplayout(1, 3))

	paramTable <- transform(subset(param,fit.index<=2),label=paste0(substr(type,1,1),fit.index,substr(estimateType,1,1)))
	paramTable <- recast(subset(paramTable,select=c(-type,-fit.index,-estimateType)),variable~label,id.var="label")
	rownames(paramTable) <- paramTable$variable
	paramTable$variable <- NULL
	paramTable <- signif(paramTable,2)
	grid.table(paramTable
			,gpar.coretext = gpar(col = "black", cex = 0.6)
			,gpar.coltext = gpar(col = "black", cex = 0.8, fontface = "bold")
			,gpar.rowtext = gpar(col = "black", cex = 0.6, fontface = "italic")
			,vp=vplayout(1, 4))
	
}


