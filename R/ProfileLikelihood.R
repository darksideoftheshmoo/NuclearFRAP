
#' merge profileLikelihoodFrap objects
#'
#' @details combines two profileLikelihoodFrap objects into one
#' @param x a profileLikelihoodFrap object
#' @param y a profileLikelihoodFrap object
#' @return a combined profileLikelihoodFrap object 
#' @export
merge.profileLikelihoodFrap <- function(x,y){
	for(i in intersect(names(x$PL1D),names(y$PL1D))){
		maxFi <- max(x$PL1D[[i]]$fit.index)
		x$PL1D[[i]] <- rbind(x$PL1D[[i]],rcell2:::conform(transform(y$PL1D[[i]],fit.index=fit.index+maxFi),to=x$PL1D[[i]]))
	}

	for(i in intersect(names(x$PL2D),names(y$PL2D))){
		maxFi <- max(x$PL2D[[i]]$fit.index)
		x$PL2D[[i]] <- rbind(x$PL2D[[i]],rcell2:::conform(transform(y$PL2D[[i]],fit.index=fit.index+maxFi),to=x$PL2D[[i]]))
	}
	
	maxFi <- max(x$param$fit.index)
	x$param <- rbind(x$param,rcell2:::conform(transform(y$param,fit.index=fit.index+maxFi),to=x$param))
	
	return(x)
}


#' run profile likelihood analysis for all cells in a experiment
#'
#' @details runs \code{\link{profileLikelihoodFrap}} in pairwise mode for all cells in the experiment
#' @note this function can take very long to run. In my Intel i5 quadcore a 10 cells experiment runs overnight (~12 hours)
#'		 Saves results to disk in .RData files, in the directory named as d$experiment
#' @param d a cell.data object with the experiment to analyze
#' @param vars a character vector with the names of the model parameters that are going to be used for the profile likelihood analysis
#' @inheritParams profileLikelihoodFrap
#' @return a list with profileLikelihoodFrap objects
#' @export
profileLikelihoodFrapExperiment <- function(d,vars=names(which(!d$parameters$FIXED.PARAMETERS)),precision=NULL,quadrantBisections=1){
	dir.create(d$experiment)
	profLikeExp <- dlply(d$data,.(pos),profileLikelihoodFrap,p=d$parameters,vars=vars,single=FALSE,pairwise=TRUE
			,saveBasename=paste0(d$experiment,"/ProfileLikelihoodMatrix-"),precision=precision,quadrantBisections=quadrantBisections,.progress = "text")
	save(profLikeExp,d,file=paste0(d$experiment,"/ProfileLikelihoodMatrix-",d$experiment,".RData"))

	l_ply(profLikeExp,function(pl){
		png(paste0(d$experiment,"/","ProfileLikelihoodMatrix-Pos",pl$pos,".png"),width=2*length(vars),2*length(vars),units="in",res=300)
		plot(pl)
		dev.off()
	},.progress = "text")

	return(invisible(profLikeExp))
}

#' run profile likelihood analysis for all pairs in a experiment
#'
#' @details runs \code{\link{profileLikelihoodFrapPair}} in pairwise mode for all pairs in the experiment
#' @note this function can take very long to run. In my Intel i5 quadcore a 10 cells experiment runs overnight (~12 hours)
#'		 Saves results to disk in .RData files, in the directory named as d$experiment
#' @param d a cell.data object with the experiment to analyze
#' @param vars a character vector with the names of the model parameters that are going to be used for the profile likelihood analysis
#' @inheritParams profileLikelihoodFrap
#' @return a list with profileLikelihoodFrap objects
#' @export
profileLikelihoodFrapPairExperiment <- function(d,vars=names(which(!d$parameters$FIXED.PARAMETERS)),precision=NULL,quadrantBisections=1){
	dir.create(d$experiment)
	profLikeExp <- dlply(d$data,.(pos),profileLikelihoodFrapPair,p=d$parameters,vars=vars,pairwise=TRUE,saveBasename=paste0(d$experiment,"/ProfileLikelihoodMatrix-")
		,precision=precision,quadrantBisections=quadrantBisections,.progress = "text")
	save(profLikeExp,d,file=paste0(d$experiment,"/ProfileLikelihoodMatrix-",d$experiment,".RData"))

	l_ply(profLikeExp,function(pl){
		png(paste0(d$experiment,"/","ProfileLikelihoodMatrix-Pos",pl$pos,".png"),width=2*length(vars),2*length(vars),units="in",res=300)
		plot(pl)
		dev.off()
	},.progress = "text")

	return(invisible(profLikeExp))
}

#' plot method for profileLikelihoodFrapPair objects
#'
#' @details plots all 1D profile likelihood analysis and, if present, all 2D analysis
#' @param object 'profileLikelihoodFrapPair' object as returned by \code{\link{profileLikelihoodFrapPair}}
#' @param ncol number of columns for facets, when only 1D analysis are persent
#' @param nrow number of rows for facets, when only 1D analysis are persent
#' @param addParam additional parameter points to be plotted
#' @param bounds list containing matrices UB and LB, with the upper and lower bounds of each parameter. (as returned by calcInitBound) 
#' @param color string indicating wich variable map to the aesthetic color
#' @export
plot.profileLikelihoodFrapPair <- function(object,ncol=NULL,nrow=NULL,addParam=NULL,bounds=NULL,color="type"){

	if(length(object$PL1D)>0){
		nVars<-length(object$PL1D)
		nameVars<-names(object$PL1D)
	} else if(length(object$PL2D)>0){
		nameVars<-unique(do.call("c",strsplit(names(object$PL2D),"-")))
		nVars<-length(nameVars)
	}

	boundaryCost <- qchisq(p=c(0.05), df=sum(!object$p$FIXED.PARAMETERS)) 
	if(!is.null(addParam)) object$param <- rbind(object$param,addParam)
	object$param1<-transform(object$param1,ncb=1+boundaryCost/cost) #norm cost boundary 
	is.pair <- "type" %in% names(object$param)
	if(is.pair){
		object$param<-transform(object$param,label=paste0(substr(type,1,1),fit.index,"\n",round(cost)))
	} else {
		object$param<-transform(object$param,label=paste0("c",fit.index,"\n",round(cost)))
	}
	object$param$fit.index<-factor(object$param$fit.index)
	object$param1$fit.index<-factor(object$param1$fit.index)
	pairwise <- length(object$PL2D)>0

	if(pairwise){
		ncol <- nVars
		nrow <- nVars
	} else {
		if(is.null(ncol)&is.null(nrow)) ncol <- ceiling(sqrt(nVars))
		if(is.null(nrow)) nrow <- ceiling(nVars/ncol)
		if(is.null(ncol)) ncol <- ceiling(nVars/nrow)
	}

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow, ncol)))

	if(length(object$PL1D)>0){
		for(i in seq_len(nVars)){
			if(pairwise){
				iRow <- i
				iCol <- i
			} else {
				iRow <- 1+(i-1)%/%ncol
				iCol <- 1+(i-1)%%ncol
			}
			if(is.pair){
				myAes <- aes_string(x=nameVars[i],y="norm.cost",color="type")
			} else {
				myAes <- aes_string(x=nameVars[i],y="norm.cost",color="factor(fit.index)")
			}
			try(print(
				ggplot(myAes,data=subset(object$PL1D[[i]],cost>0)) + rcell2::zoom(y=c(1,max(c(3,2*object$param1$ncb))),expand.y=c(-0.2,0.2)) +
					geom_line(linetype=1) +
					#geom_line(linetype=1,data=subset(object$PL1D[[i]],cost>0 & ((norm.cost<=ncb.m & type=="mother")|(norm.cost<=ncb.d & type=="daughter")))) +
					geom_point(data=object$param,y=1) + 
					geom_errorbar(aes(y=ncb,ymax=ncb,ymin=ncb),linetype=2,data=object$param1) +
					geom_text(aes(label=label),y=1,vjust=0.5,size=4,data=object$param) + 
					theme(legend.position="none") + 
					scale_color_manual(values=gg_color_hue(4)[c(3,1)])
			, vp = vplayout(iRow,iCol)))
		}
	}

	if(pairwise){
		for(i in 1:(nVars-1)){
			for (j in (i+1):nVars){
				xlim <- (ylim <- NULL)
				if(!is.null(bounds)){
					xlim <- c(bounds$LB[,nameVars[j]],bounds$UB[,nameVars[j]])
					ylim <- c(bounds$LB[,nameVars[i]],bounds$UB[,nameVars[i]])
				}
				PL2D.db <- object$PL2D[[paste0(nameVars[i],"-",nameVars[j])]]
				print(plotPL2D(PL2D.db,x=nameVars[j],y=nameVars[i],param=object$param, signifCost=qchisq(p=c(0.05), df=sum(!p$FIXED.PARAMETERS))
							   ,tauBoundary=FALSE,xlim=xlim,ylim=ylim,color=color)
						,vp = vplayout(i,j))
			}
		}
	}
}

#' plot method for profileLikelihoodFrap objects
#'
#' @details plots all 1D profile likelihood analysis and, if present, all 2D analysis
#' @inheritParams plot.profileLikelihoodFrapPair 
#' @export
plot.profileLikelihoodFrap <- plot.profileLikelihoodFrapPair


#' run a profile likelihood analysis for a FRAP pair
#'
#' @param db a data.frame with the data of a pair (required if expData=NULL). 'time', 'FRAP', 'type' and nuclear and cytoplasmic fluorescence variables required. 
#' @inheritParams profileLikelihoodFrap
#' @export
profileLikelihoodFrapPair <- function(db,p,vars=c("log10_kEVfrac","log10_kI"),param=NULL,single=TRUE,pairwise=FALSE,nGrid=10,nStarts=5
										,saveBasename=NULL,overwrite=FALSE,precision=NULL,quadrantBisections=1){


	pos<-unique(db$pos)
	if(length(pos)!=1) stop("only one pos can be passed to analyzeFrapPair")

	if(!is.null(saveBasename)){
		outputFilename <- paste0(saveBasename,"pos",pos,".RData")
		if(file.exists(outputFilename) & !overwrite){
			load(outputFilename)
			if("output" %in% ls()) return(output)
		}
	}

	if(is.null(param)) param<-runModelFitPair(db,p=p)
	param1<-subset(param,fit.index==1)
	nVars<-length(vars)

	PL1D<-list()
	if(single){
		cat("1D Profile Likelihood: ")
		for(v in vars){
			cat(v," ")
			PL.db<-runProfileLikelihoodGridPair(db,param1,p=p,scan=v,nGrid=nGrid,nStarts=nStarts)
			PL.db<-ddply(PL.db,.(type),function(df){
				myCost<-param1[param1$type==unique(df$type),"cost"]
				transform(df,norm.cost=cost/myCost,p005=cost<myCost)
			})
			PL1D[[v]]<-PL.db
		}
		cat("\n")
	}

	PL2D<-list()
	if(pairwise){
		cat("2D Profile Likelihood: \n")
		for(i in 1:(nVars-1)){
			for (j in (i+1):nVars){
				pairwiseName <- paste0(vars[i],"-",vars[j])
				cat(pairwiseName," ")
				PL2D[[pairwiseName]] <- runProfileLikelihoodContourPair(db,param,p=p,scan=vars[c(i,j)],nStarts=nStarts
																		  ,precision=precision,quadrantBisections=quadrantBisections)
			}
			cat("\n")
		}
		cat("\n")
	}

	output<-list(PL1D=PL1D,PL2D=PL2D,param=param,param1=param1,pos=pos,p=p)
	class(output) <- c("profileLikelihoodFrapPair","list")
	if(!is.null(saveBasename)) save(output,file=outputFilename)
	return(output)
}


#' run a profile likelihood analysis for FRAP data 
#'
#' @param db a data.frame with the data of a pair (required if expData=NULL). 'time', 'FRAP' and nuclear and cytoplasmic fluorescence variables required. 
#' @param p list of options for the FRAP experiment
#' @param param a data.frame with the best fits parameters
#' @param single boolean indicating if 1D profile likelihoods should be calculated
#' @param pairwise boolean indicating if pairwise profile likelihood contours should be done
#' @param nStarts number of random start points use by matlab's MultiStart function in order to get the global minimum
#' @param nGrid number of points for the 1D analysis for each var 
#' @param saveBasename if not NULL, a charcater with the basename of the RData file where the results of the profile likelihood analysis will be saved.
#' @param overwrite boolean indicating what to do if a .RData with the same name as the one to be created is found. 
#'		            If overwrite=FALSE, the .RData is loaded and its content returned by the function, therefore avoiding to re-run the analysis.
#'					If overwrite=TRUE, the analysis is redone and the file overwritten
#' @inheritParams runProfileLikelihoodContour
#' @export
profileLikelihoodFrap <- function(db,p,vars=c("log10_kEVfrac","log10_kI"),param=NULL,single=TRUE,pairwise=FALSE,nGrid=10
									,nStarts=5,saveBasename=NULL,overwrite=FALSE,precision=NULL,quadrantBisections=1){

	pos<-unique(db$pos)
	if(length(pos)!=1) stop("only one pos can be passed to analyzeFrapPair")

	if(!is.null(saveBasename)){
		outputFilename <- paste0(saveBasename,"pos",pos,".RData")
		if(file.exists(outputFilename) & !overwrite){
			load(outputFilename)
			if("output" %in% ls()) return(output)
		}
	}

	if(is.null(param)) param<-runModelFit(db,p=p)
	param1<-subset(param,fit.index==1)
	nVars<-length(vars)

	PL1D<-list()
	if(single){
		cat("1D Profile Likelihood: ")
		for(v in vars){
			cat(v," ")
			PL.db <- runProfileLikelihoodGrid(db,param1,p=p,scan=v,nGrid=nGrid,nStarts=nStarts)
			PL.db <- transform(PL.db,norm.cost=cost/param1$cost,p005=cost<param1$cost)
			PL1D[[v]]<-PL.db
		}
		cat("\n")
	}

	PL2D<-list()
	if(pairwise){
		cat("2D Profile Likelihood: \n")
		for(i in 1:(nVars-1)){
			for (j in (i+1):nVars){
				pairwiseName <- paste0(vars[i],"-",vars[j])
				cat(pairwiseName," ")
				PL2D[[pairwiseName]] <- runProfileLikelihoodContour(db,param,p=p,scan=vars[c(i,j)],nStarts=nStarts
																	  ,precision=precision,quadrantBisections=quadrantBisections)
			}
			cat("\n")
		}
		cat("\n")
	}

	output<-list(PL1D=PL1D,PL2D=PL2D,param=param,param1=param1,pos=pos,p=p)
	class(output) <- c("profileLikelihoodFrap","list")
	if(!is.null(saveBasename)) save(output,file=outputFilename)
	return(output)
}



#' create a 2D plot of a profile likelihood boundary
#'
#' @param PL2D.db data.frame as returned by \code{\link{runProfileLikelihoodGrid}} or \code{\link{runProfileLikelihoodContour}}
#' @param x string with name of model parameter to use in the x axis
#' @param y string with name of model parameter to use in the y axis
#' @param param a data.frame with the best fits parameters
#' @param signifCost chi square cost significant difference, defines boundary for \code{\link{runProfileLikelihoodGrid}} datasets
#' @param color string indicating wich variable map to the aesthetic color
#' @param shape string indicating wich variable map to the aesthetic shape. Use "fit.index" or "type".
#' @param alpha	transparency level of the contours
#' @param tauBoundary boolean indicating if tau boundaries should be drawn
#' @param costLabel boolean indicating if cost should be indicated in the label
#' @param xlim numeric vector of length with boundaries of x axes. If NULL boundaries are determined from the data. 
#' @param ylim numeric vector of length with boundaries of y axes. If NULL boundaries are determined from the data.
#' @return a ggplot object with the required plot.  
#' @export
plotPL2D<- function(PL2D.db,x,y,param,signifCost=NULL,color="type",shape="fit.index",alpha=1,tauBoundary=TRUE,costLabel=FALSE,xlim=NULL,ylim=NULL){
	if(!is.null(PL2D.db)){ #if there is profile likelihood data
		PL2D.db$fit.index<-factor(PL2D.db$fit.index)
		if(all(c('theta','lambda')%in%names(PL2D.db))){ #if there is only contour data
			contour.pts<-rep(TRUE,times=nrow(PL2D.db))
		} else { #if there is grid, and countour data
			contour.pts<-apply(is.na(PL2D.db),1,any) 
		}
		PL2D.contourPts.db<-PL2D.db[contour.pts,]
		PL2D.db<-PL2D.db[!contour.pts,]
	}
	param$fit.index<-factor(param$fit.index)
	
	if(!"type" %in% names(param) & !is.null(PL2D.db)){
		param <- transform(param,type=factor("cell"))
		if(nrow(PL2D.db)>0) 
			PL2D.db <- transform(PL2D.db,type=factor("cell"))
		if(nrow(PL2D.contourPts.db)>0) 
			PL2D.contourPts.db <- transform(PL2D.contourPts.db,type=factor("cell"))
	}

	if(!"type" %in% names(param)) param$type <- factor("cell")		
	if("pos"%in%names(param)){
		param$pos<-factor(param$pos)
		param<-transform(param,label=paste0(pos,substr(type,1,1),".",fit.index)
						  ,labelCost=paste0(pos,substr(type,1,1),".",fit.index,"\n",round(cost)))
	} else {
		param<-transform(param,label=paste0(substr(type,1,1),".",fit.index)
						  ,labelCost=paste0(substr(type,1,1),".",fit.index,"\n",round(cost)))
	}

	p1<-ggplot(PL2D.db,aes_string(x=x, y=y, z = "cost", color=color, shape=shape, linetype="type")) + 
		geom_point(data=param,size=2) 
	if(!costLabel){
		p1 <- p1 + geom_text(aes(label=label),data=param, vjust=-0.25, hjust=1.1, size=4)
	} else {
		p1 <- p1 + geom_text(aes(label=labelCost),data=param, vjust=0.5, size=4)
	}

	if(!is.null(PL2D.db)){
		for(k in 1:nrow(param)){ 
			if(nrow(PL2D.db)>0){
				tmp <- subset(PL2D.db,pos==param$pos[k]&type==param$type[k]&fit.index==as.character(param$fit.index[k]))
				b <- param$cost[k] + signifCost
				if(!is.na(tmp$cost[1]))
					if(any(tmp$cost<b) & any(tmp$cost>b))	
						p1 <- p1 + stat_contour(breaks=b,data=tmp,linetype=2+(param$type[k]=="daughter"),alpha=alpha)
			}
		}
		#tmpC <- subset(PL2D.contourPts.db,pos==param$pos[k]&type==param$type[k]&fit.index==as.character(param$fit.index[k]))
		if("pos"%in%names(PL2D.contourPts.db)) {
			tmpC <- ddply(PL2D.contourPts.db,.(pos,type,fit.index),function(df){rbind(df,df[1,])})
		} else {
			tmpC <- ddply(PL2D.contourPts.db,.(type,fit.index),function(df){rbind(df,df[1,])})
		}
		if(!is.na(tmpC$cost[1]))
			p1 <- p1 + geom_path(aes(linetype=type),data=subset(tmpC,abs(cost)>sqrt(.Machine$double.eps)),fill=NA,alpha=alpha)
				#if("lambda" %in% names(tmpC)) p1 <- p1 + geom_point(data=subset(tmpC,lambda > 0),shape=3)
	}
	#color scale
	p1 <- p1 + theme(legend.position="none") 
	if(color=="type"){
		p1 <- p1 + scale_color_manual(values=gg_color_hue(4)[c(1,3)])
	} else if(color=="fit.index") {
		p1 <- p1 + scale_colour_brewer(palette="Set2")
	}
	#x and y scale based on xlim and ylim
	if(!is.null(xlim)&!is.null(ylim)){
		p1 <- p1 + rcell2::zoom(x=xlim,y=ylim)
	} else if(!is.null(xlim)) {
		p1 <- p1 + rcell2::zoom(x=xlim)
	} else if(!is.null(ylim)) {
		p1 <- p1 + rcell2::zoom(y=ylim)
	}

	if(isTRUE(tauBoundary)){
		#parameter boundaries
		bound.db<-data.frame(kEVfrac=seq(-3,1,0.001))
		bound.db<-rbind(bound.db,data.frame(kEVfrac=log10((1/c(0.2,0.1))-1e-3)))
		bound.db<-transform(bound.db,l1=log10(1/0.2-10^kEVfrac),l1b=log10(1/(0.2/2)-10^kEVfrac),l2=log10(1/15-10^kEVfrac),l2b=log10(1/(2*15)-10^kEVfrac))
		p1 <- p1 +
		cplot(bound.db,l1~kEVfrac,color="orange",geom="line",layer=T,linetype=2) + 
		cplot(bound.db,l2~kEVfrac,color="orange",geom="line",layer=T,linetype=2) + 
		cplot(bound.db,l1b~kEVfrac,color="orange",geom="line",layer=T,linetype=3) + 
		cplot(bound.db,l2b~kEVfrac,color="orange",geom="line",layer=T,linetype=3) 
	}

	return(p1)
}



#' run the profile likelihhod analysis on 'fix directions' from the fitted point 
#'
#' @inheritParams runProfileLikelihoodGrid
#' @param scan a character vector of length 2 with the names of the parameters that define the plane for the profile likelihood contour.
#' @param boundaryCost chi square cost increment between best fit and boundary. 
#' @param quadrantBisection number of times each quadrant is bisected to define the direction of the boundary points
#' @param precision precision of the boundary. If NULL min(x0,y0)/100 is used
#' @return a data.frame with the points along the directions defined by theta where cost = min.cost + boundaryCost
#' @export
runProfileLikelihoodContour<-function(db=NULL,param=NULL,expData=NULL,fixed=p$FIXED.PARAMETERS,scan=c("kEVfrac","kI")
										,boundaryCost=qchisq(p=c(0.05), df=sum(!p$FIXED.PARAMETERS)),p=param.FRAP.pairs,nStarts=5
										,precision=NULL,quadrantBisections=1){

	if(is.null(param)) stop("param required")
	if(is.null(expData)&is.null(db)) return(NULL)	
	if(is.null(expData)) expData<-makeEXPDATA(db,p=p)
	ib <- calcInitBound(expData,p)
	scan <- p$MODEL.PARAMETERS[p$MODEL.PARAMETERS%in%scan] #reordering scan variables
	maxLambda <- sqrt(sum((ib$UB[ib$names%in%scan]-ib$LB[ib$names%in%scan])^2))
	if(!all(scan%in%p$MODEL.PARAMETERS)) stop("some scan params not in models parameters")
	scanVector <- as.numeric(p$MODEL.PARAMETERS%in%scan)
	
	x0<-abs(param[param$fit.index==1,scan[1]])
	y0<-abs(param[param$fit.index==1,scan[2]])
	if(abs(x0) < sqrt(.Machine$double.eps)) x0 <- 1
	if(abs(y0) < sqrt(.Machine$double.eps)) y0 <- 1
	if(is.null(precision)) precision=min(x0/3,y0/3,maxLambda/1000)

	#df <- param
	output<-
	ddply(param,.(fit.index),function(df){	
		theta <- seq(-pi/2,pi,pi/2)
		df.x0<-as.matrix(df[,p$MODEL.PARAMETERS])
		pl.db<-data.frame()
		for(tb in seq_len(quadrantBisections+1)){

			tmp_new <- profileLikelihoodContour( EXPDATA=expData, theta=t(as.matrix(theta))
							, x0=df.x0, LB=ib$LB, UB=ib$UB, scan=scanVector, fixed=fixed
							, nStarts=nStarts, cost=df$cost + boundaryCost, maxLambda=maxLambda, precisionLambda=precision)
							
			tmp <- tmp_new
							
			if(is.null(tmp)) {
				pl.db <- rbind(pl.db,data.frame(t(as.matrix(rep(NA,length(scan) + 3)))))
			} else {
				pl.db <- rbind(pl.db,data.frame(tmp))
			}
			pl.db<-pl.db[order(pl.db[,3]),]
			theta<-numeric()
			for(i in seq_len(nrow(pl.db))){
				x1<-mean(c(pl.db[i,1],pl.db[i%%nrow(pl.db)+1,1]))
				y1<-mean(c(pl.db[i,2],pl.db[i%%nrow(pl.db)+1,2]))
				dx <- x1-df[1,scan[1]]
				dy <- y1-df[1,scan[2]]
				if(abs(dx) > x1/1000 & abs(dy) > y1/1000){
					theta<-c(theta,atan2(dy,dx))
				} else {
					theta<-c(theta,pl.db[i,3]+pi/(4*tb))
				}
			}
		}
		names(pl.db)<-c(scan,"theta","lambda","cost")
		return(pl.db)
	})

	# #df <- param
	# output<-
	# ddply(param,.(fit.index),function(df){	
		# #gridParam <- as.matrix(df[,scan])  
		# df.x0<-as.matrix(df[,p$MODEL.PARAMETERS])
		# tmp <- callMatlabFun("profileLikelihoodContour", EXPDATA=expData, theta=t(as.matrix(theta))
						# , x0=df.x0, LB=ib$LB, UB=ib$UB, scan=scanVector, fixed=fixed
						# , nStarts=nStarts, cost=df$cost + boundaryCost, maxLambda=maxLambda, precisionLambda=abs(min(x0,y0)/100))
		# if(is.null(tmp)) {
			# return(data.frame(t(as.matrix(rep(NA,length(scan) + 3)))))
		# } else {
			# return(data.frame(tmp))
		# }
	# })

	names(output)<-c("fit.index",scan,"theta","lambda","cost")
	return(output)
}

#' run the profile likelihhod analysis on 'fix directions' from the fitted point, for a pair of cells
#'
#' @inheritParams runProfileLikelihoodContour
#' @return a data.frame with the points along the directions defined by theta where cost = min.cost + boundaryCost
#' @export
runProfileLikelihoodContourPair<-function(...) pairFun(runProfileLikelihoodContour,...)


#' run the profile likelihhod analysis on a regular grid of points in the parameter space
#'
#' @details wrapper over matlabs function 'profileLikelihoodGrid.m'
#' @param db a data.frame with the data of a single cell (required if expData=NULL). 'time', 'FRAP' and nuclear and cytoplasmic fluorescence variables required. 
#' @param param a data.frame with the parameters for which the profile likelihood cost is wanted. 
#' @param expData a EXPDATA list as returned by \code{makeEXPDATA}, which contains the data for a cell (required if db=NULL). 
#' @param p list of options for the FRAP experiment
#' @param fixed a named interger vector specifying wich parameters are fixed (1) and which are free (0)
#' @param scan either a character vector with a variables names to be scaned, or a named list with the values of the parameters with which to make a regular grid
#' @param absolute boolean indicating if the values passed in scan should be treated as 'absolute' values, or increments with respect to the optimum fit passed in \code{param}
#' @param nStarts number of random start points use by matlab's MultiStart function in order to get the global minimum
#' @param nGrid number of points in the regular grid, when scan is a character vector. 
#' @param boundaryCost chi square cost increment between best fit and boundary. Used to estimate scan region when scan is a character vector
#' @param theta numeric vector with angles (in radians) that define the directions for \code{runProfileLikelihoodContour}
#' @param rescale boolean indicating if the coordinates should be rescale before setting the direction defined by theta. See \code{runProfileLikelihoodContour}.
#' @return a data.frame with the points in the parameter space and the profile likelihood cost
#' @export
runProfileLikelihoodGrid<-function(db=NULL,param=NULL,expData=NULL,p=param.FRAP.pairs,fixed=p$FIXED.PARAMETERS,scan=c("kEVfrac","kI"),absolute=FALSE
									,nStarts=5,nGrid=10,boundaryCost=qchisq(p=c(0.05), df=sum(!p$FIXED.PARAMETERS)),theta=seq(pi/4,2*pi,pi/4),rescale=TRUE){

	if(is.null(param)) stop("param required")
	if(is.null(expData)&is.null(db)) return(NULL)	
	if(is.null(expData)) expData<-makeEXPDATA(db,p=p)
	ib<-calcInitBound(expData,p)

	#df<-param[1,]
	output<-
	ddply(param,.(fit.index),function(df){	
		#estimating size of parameter scan
		plContour.db <- NA
		if(is.character(scan)) {
			if(length(scan)==2 & !any(is.na(theta))) 
				plContour.db<-runProfileLikelihoodContour(db=db,param=df,expData=expData,scan=scan,fixed=fixed,boundaryCost=boundaryCost,p=p,nStarts=nStarts,theta=theta,rescale=rescale)
			fixedList<-list()
			if(any(is.na(plContour.db))){
				for(fp in scan){
					B.index<-which(ib$names==fp)
					fixedList[[fp]] <- seq(ib$LB[B.index],ib$UB[B.index],length.out=nGrid) 
				}
			} else {
				for(fp in scan){
					fixedList[[fp]] <- makeSeq(range(plContour.db[,fp]),nGrid) 
				}
			}
			absolute <- TRUE
			scan <- fixedList
		}

		gridParam <- as.matrix(do.call(expand.grid,scan))
		scanParam<-names(scan)
		if(!all(scanParam%in%p$MODEL.PARAMETERS)) stop("some scan params not in models parameters")
		scanVector <- as.numeric(p$MODEL.PARAMETERS%in%scanParam)
		if(!absolute) gridParam <- gridParam + matrix(as.numeric(df[1,scanParam]),nrow=dim(gridParam)[1],ncol=length(scanParam),byrow = TRUE)

		df.x0<-as.matrix(df[,p$MODEL.PARAMETERS])

		profL_new <- profileLikelihoodGrid( EXPDATA=expData
						,x0=df.x0, LB=ib$LB, UB=ib$UB, fixed=t(fixed)
						,gridParam=gridParam, scan=t(scanVector), nStarts=nStarts )

		profL <- profL_new

		profL<-data.frame(profL)
		if(length(profL)==0) profL <- data.frame(t(rep(NA,length(ib$names)+1)))
		names(profL)<-c(ib$names,"cost")

		#adding contour points to profL 
		if(!any(is.na(plContour.db))){
			plContour.db<-transform(plContour.db,cost=cost+boundaryCost)
			profL<-rbind(profL,rcell2:::conform(plContour.db,to=profL))
		}

		return(data.frame(profL))
	})

	names(output)<-c("fit.index",ib$names,"cost")
	return(output)
}

#' run the profile likelihhod analysis on a regular grid of points in the parameter space, for a cell pair
#'
#' @details pair version of \code{\link{runProfileLikelihoodGrid}}. db, param, and expData should contain info about both cells
#' @inheritParams runProfileLikelihoodGrid
#' @return a data.frame with the points in the parameter space and the profile likelihood cost, for the cell pair
#' @export
runProfileLikelihoodGridPair<-function(...) pairFun(runProfileLikelihoodGrid,...)

#' get the profile likelihood cost for a point in the parameter space for a single cell
#'
#' @details wrapper over matlabs function 'profileLikelihoodCost.m'
#' @param db a data.frame with the data of a single cell (required if expData=NULL). 'time', 'FRAP' and nuclear and cytoplasmic fluorescence variables required. 
#' @param param a data.frame with the parameters for which the profile likelihood cost is wanted. 
#' @param expData a EXPDATA list as returned by \code{makeEXPDATA}, which contains the data for a cell (required if db=NULL). 
#' @param fixed a character vector with the names of the model variables that are to be kepet fixed during the profile likelihood cost calculation
#' @param p list of options for the FRAP experiment
#' @param nStarts number of random start points use by matlab's MultiStart function in order to get the global minimum
#' @return a data.frame with the point in the parameter space and the profile likelihood cost
#' @export
runProfileLikelihoodCost<-function(db=NULL,param=NULL,expData=NULL,fixed=c("kEVfrac","kI"),p=param.FRAP.pairs,nStarts=5){

	if(is.null(param)) stop("param required")
	if(is.null(expData)&is.null(db)) return(NULL)	
	if(is.null(expData)) expData<-makeEXPDATA(db,p=p)

	ib<-calcInitBound(expData,p)
	param$index<-seq_len(dim(param)[1])
	
	if(!all(fixed%in%p$MODEL.PARAMETERS)) stop("some fixed params not in models parameters")
	fixedVector <- as.numeric(p$MODEL.PARAMETERS%in%fixed)

	output<-
	ddply(param,.(index),function(df){	
		df.x0<-as.matrix(df[,p$MODEL.PARAMETERS])

		plCost_new<-profileLikelihoodCost(EXPDATA=expData
						, x0=df.x0, LB=ib$LB, UB=ib$UB, fixed=fixedVector
						, nStarts=nStarts)

		plCost <- plCost_new

		return(data.frame(plCost))
	})
	names(output)<-c("index",ib$names,"cost")
	return(output)
}

#' get the profile likelihood cost for a point in the parameter space for a cell pair
#'
#' @details pair version of \code{runProfileLikelihoodCost}
#' @param ... pair arguments for runProfileLikelihoodCost
#' @return a data.frame with the point in the parameter space and the profile likelihood cost, and a 'type' column
#' @export
runProfileLikelihoodCostPair<-function(...) pairFun(runProfileLikelihoodCost,...)


#' get a sample point of a profile likelihood analysis
#'
#' @param PL.db data.frame as returned by \code{\link{runProfileLikelihoodGrid}}
#' @param cost target chi square cost. The returned point will be the one in PL.db, closest in cost to this parameter.
#' @param point target point in paramter space. The returned point will be the one in PL.db closest to this point.  
#' @return a data.frame containing an element of PL.db
#' @export
getProfileLikelihoodPoint<-function(PL.db,param,cost=NULL, point=NULL){
	output<-data.frame()
	if("type"%in%names(param) & "type"%in%names(PL.db)){
		for(i in seq_along(cost)){
			output<-rbind(output,
				data.frame(ddply(param,.(fit.index,type),function(df){
					tmp<-PL.db[PL.db$type==df$type & PL.db$fit.index==df$fit.index,]
					tmp[which.min(abs(tmp$cost-(df$cost + cost[i]))),]
				}),cost.index=i))
		}
		if(!is.null(point)){
			output<-rbind(output,
				data.frame(ddply(param,.(fit.index,type),function(df){
					tmp<-PL.db[PL.db$type==df$type & PL.db$fit.index==df$fit.index,]
					distance <- numeric()
					for(i in seq_along(point))
						distance <- distance + (tmp[,names(point)[i]]-point[i])^2
					distance <- sqrt(distance)		
					tmp[which.min(distance),]
				}),cost.index=NA))
		}
	} else {
		stop("not implemented yet")
	}

	return(output)
}
