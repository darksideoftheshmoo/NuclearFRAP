#' plot all the fits and respective residuals
#'
#' @param frap a frap or frapPair object, as returned by \code{\link{analyzeFrap}} or \code{\link{analyzeFrapPair}}
#' @export
plotResidualsAllFits <- function(frap){

	if(class(frap)=="frapPair"){
		isPair = TRUE
		idVars = c("FRAP","time","type")
		measureVars = c("f.nuc.d","f.cyt.d","f.nuc.m","f.cyt.m")
	} else if (class(frap)=="frap") {
		isPair = FALSE
		idVars = c("FRAP","time")
		measureVars = c("f.nuc","f.cyt")
	} else {
		stop("unknown class for frap")
	}

	mdb <- melt(subset(frap$db,select=c(idVars,measureVars)),id.var=idVars)
	frap$sim.rdb <- rename(frap$sim.rdb,c("value"="pred"))
	dataSim.db <- join(frap$sim.rdb,mdb,by=c(idVars,"variable"))
	dataSim.db <- transform(dataSim.db, res = value - pred, norm.res = (value - pred)/sd )
	dataSim.db$fit.index <- factor(dataSim.db$fit.index)
	if(isPair) dataSim.db$fit.index <- interaction(dataSim.db$fit.index,dataSim.db$type,drop=TRUE)
	nFits <- length(levels(dataSim.db$fit.index))
	nFraps <- length(unique(dataSim.db$FRAP))

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nFits, 6 + (nFraps>1))))
	for( i in seq_len(nFits)){
		pdb <- subset(dataSim.db, fit.index==levels(dataSim.db$fit.index)[i])
		acf.db<-ddply(pdb,.(variable),function(df){
			a <- acf(df$res,lag.max=10,plot=FALSE)
			data.frame(autoCorr=a$acf,lag=a$lag)
		})

		crossCorrelationBetweenFraps <- function(db){
			fs<-unique(db$FRAP)
			output <- data.frame()
			for(i in 1:(length(fs)-1)){
				for(j in (i+1):length(fs)){
					output <- rbind(output, data.frame(
						desc=paste0(fs[i]," vs ",fs[j])
						,crossCorr=ccf(db$res[db$FRAP==fs[i]],db$res[db$FRAP==fs[j]],lag.max=0,plot=FALSE)$acf[1]
					))
				}
			}
			return(output)
		}
		ccf.db <- NULL
		if(length(unique(pdb$FRAP))>1) ccf.db <- ddply(pdb,.(variable),crossCorrelationBetweenFraps)

		p1<-
		cplot(pdb,value~time,geom=c("point"),group=interaction(FRAP,variable),color=variable
				,size=1,alpha=1, QC.filter = FALSE, main=levels(dataSim.db$fit.index)[i])+
		cplot(pdb,value~time,geom=c("line"),group=interaction(FRAP,variable),color=variable
				,size=0.5,alpha=0.25,layer=T, QC.filter = FALSE)+
		cplot(pdb,pred~time,ymax=pred+sd,ymin=pred-sd,geom="smooth",stat="identity",group=interaction(FRAP,variable)
			,layer=T,na.rm=FALSE,size=1,color=variable,alpha=0.2)
		p1 <- p1 + theme(legend.position="none") 

		p2<-
		cplot(pdb,norm.res~time,facets=~variable,color=variable)+
		cplot(pdb,norm.res~time,geom="line",linetype=2,size=0.5,alpha=0.5
			,layer=T,facets=~variable,group=FRAP,color=variable)+
		geom_hline(yintercept=1,linetype=3)+
		geom_hline(yintercept=-1,linetype=3)+
		facet_wrap(~variable,ncol=1) + 
		theme(legend.position="none") 

		p3<-
		ggplot(pdb,aes(x=norm.res)) + geom_histogram(aes(y = ..density..),binwidth=0.25,alpha=0.3) +
		stat_function(fun = dnorm, colour = gg_color_hue(4)[4], arg = list(mean = 0, sd=1),size=1.5,alpha=1) +
		facet_wrap(~variable,ncol=1,scale="free") + theme(legend.position="none")

		p4<-
		cplot(acf.db, autoCorr ~ lag, facets=~variable, geom=c("bar"), stat="identity") + 
		facet_wrap(~variable,ncol=1) + theme(legend.position="none") +
		yzoom(c(-1,1),expand.y=c(-0.1,0.1))+
		geom_text(aes(label=autoCorr),x=0,y=-1,hjust=0,vjust=0,size=4
			,data=aggregateBy(acf.db,.(variable),subset=lag>0,FUN=function(x)signif(sum(abs(x)),3)))

		if(!is.null(ccf.db))
			p5<-
			cplot(ccf.db, crossCorr ~ desc, facets=~variable, geom=c("bar"), stat="identity") + 
			theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=0.5),axis.title.x=element_blank())+
			yzoom(c(-1,1),expand.y=c(-0.1,0.1)) + facet_wrap(~variable,ncol=1) + labs(y="cross corrleation among FRAPs")+
			geom_text(aes(label=crossCorr),x=1,y=-1,hjust=0,vjust=0,size=4
				,data=aggregateBy(ccf.db,.(variable),select="crossCorr",FUN=function(x)signif(sum(abs(x)),3)))

		print(p1, vp = vplayout(i, 1:2))
		print(p2, vp = vplayout(i, 3:4))
		print(p3, vp = vplayout(i, 5))
		print(p4, vp = vplayout(i, 6))
		if(!is.null(ccf.db)) print(p5, vp = vplayout(i, 7))
	}

}



#' gets the residuals of the data to the model, for a pair of cells
#' 
#' @param db a data.frame with the data of a pair cell (required if expData=NULL). 'time', 'FRAP', 'type', nuclear and cytoplasmic fluorescence variables for the mother and daughter required. 
#' @param param a data.frame with the parameters of the fit (required). Must have 'type' vatiable and a variable for each free parameter of the model. 
#' @param a list of EXPDATA as returned by \code{makePairEXPDATA}, which contains the data for a pair of cells (required if db=NULL). 
#' @param  p list of options for the FRAP experiment
#' @param melt boolean specifying if the output should be 'melted' by \code{melt}
#' @return a data.frame. If melt=FASE => variables "FRAP","time","res.nuc","res.cyt","type". If melt=TRUE => variables "FRAP","time","type","variable","value"
#' @export
getResidualsPair<-function(...) pairFun(getResiduals,...)

#' gets the residuals of the data to the model
#' 
#' @param db a data.frame with the data of a single cell (required if expData=NULL). 'time', 'FRAP', nuclear and cytoplasmic fluorescence variables required. 
#' @param param a data.frame with the parameters of the fit (required). Must have a variable for each free parameter of the model, and a sigle row.
#' @param expData a EXPDATA list as returned by \code{makeEXPDATA}, which contains the data for a cell (required if db=NULL). 
#' @param  p list of options for the FRAP experiment
#' @param melt boolean specifying if the output should be 'melted' by \code{melt}
#' @return a data.frame. If melt=FASE => variables "FRAP","time","res.nuc","res.cyt". If melt=TRUE => variables "FRAP","time","variable","value"
#' @export
getResiduals<- function(db=NULL,param=NULL,expData=NULL,p=param.FRAP.pairs,melt=FALSE){
	if(is.null(param)) stop("param required")
	if(is.null(expData)&is.null(db)) return(NULL)	
	if(is.null(expData)) expData<-makeEXPDATA(db,p=p)

	output<-rename(getSimData(param=param,expData=expData,p=p),c("f.nuc"="pred.nuc","f.cyt"="pred.cyt"))
	output<-cbind(output,rename(data.frame(do.call(rbind,replicate(dim(param)[1],expData$y,simplify=FALSE))),c("X1"="exp.nuc","X2"="exp.cyt")))
	output<-transform(output,res.nuc=exp.nuc-pred.nuc,res.cyt=exp.cyt-pred.cyt)

	if(isTRUE(melt)) output<-melt(subset(output,select=c("fit.index","FRAP","time","res.nuc","res.cyt")),id.vars=c("fit.index","FRAP","time"))

	return(output)
}

#' gets the residuals of a local estimation model, on the variables were the FRAP is NOT taking place
#' 
#' @details gets the residuals of a \code{loess} model for the mother nuc and cyt when the FRAP is done on the daughter and vice-versa
#' @param db a data.frame with the data of a single pair. 'time', 'FRAP', 'type', nuclear and cytoplasmic fluorescence variables for mother and daughter required
#' @param  p list of options for the FRAP experiment
#' @param melt boolean specifying if the output should be 'melted' by \code{melt}
#' @return a data.frame. If melt=FASE => variables "FRAP","time","res.nuc","res.cyt","type". If melt=TRUE => variables "FRAP","time","type","variable","value"
#' @export
getCrossResidualsPair <- function(db,p=param.FRAP.pairs,melt=FALSE){
	ut<-unique(db$type)

	p$SIGMA.FUN<-I
	m.expData<-makeEXPDATA(db[db$type=="mother",c("time","f.nuc.m","f.cyt.m","FRAP")],p=p)
	d.expData<-makeEXPDATA(db[db$type=="daughter",c("time","f.nuc.d","f.cyt.d","FRAP")],p=p)

	#crossData has the variables of the mother in the daughter FRAPs and vice-versa
	m.crossData<-makeEXPDATA(db[db$type=="daughter",c("time","f.nuc.m","f.cyt.m","FRAP")],p=p)
	d.crossData<-makeEXPDATA(db[db$type=="mother",c("time","f.nuc.d","f.cyt.d","FRAP")],p=p)

	#calculating sd for each variable from the "unused" data 
	#standard deviation from the "loess" model calculated by makeEXPDATA with p$SIGMA.FUN<-I
	m.res<-data.frame(FRAP=m.crossData$FRAP,time=m.crossData$time,m.crossData$y-m.crossData$std,type="mother")
	d.res<-data.frame(FRAP=d.crossData$FRAP,time=d.crossData$time,type="daughter",d.crossData$y-d.crossData$std)

	res<-rbind(m.res,d.res)
	names(res)<-c("FRAP","time","res.nuc","res.cyt","type")

	if(isTRUE(melt)){
		return(melt(res,id.vars=c("FRAP","time","type")))
	}else{
		return(res)
	}
}
