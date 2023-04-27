#' gets the simulated model with added noise for a single cell and a set of parameters
#'
#' @details use this function to get simulated data ready to be re-fitted
#' @param db a data.frame with the data of a single cell (required if expData=NULL). 'time', 'FRAP' and nuclear and cytoplasmic fluorescence variables required. 
#' @param p list of options for the FRAP experiment
#' @param param a data.frame with the parameters of the fit. If missing \code{runModelFit} is called. 
#' @param sigmaFactor multiplying factor to multiplicate the sigma before generating the artificial residuals
#' @return a data.frame, as required by runModelFit and other functions 
#' @export
getSimulatedNoiseData <- function(db,p,param=NULL,sigmaFactor=1){
	if(is.null(param)) param <- subset(runModelFit(db,p=p),fit.index==1)
	simData <- subset(getSimData(db=db,param=param,p=p),select=-fit.index)
	N <- dim(simData)[1]
	simData <- transform(simData, f.nuc=f.nuc + sigmaFactor*sd.f.nuc*rnorm(N)
								, f.cyt=f.cyt + sigmaFactor*sd.f.cyt*rnorm(N)) 
	simDb <- cbind(simData,db[,setdiff(names(db),names(simData))])
	return(simDb)
}


#' restrucutres the simulated data for a single cell
#'
#' @details returns a restructured data.frame more suitable for plotting
#' @param db a data.frame as returned by \code{getSimData}
#' @return a restructured data.frame
#' @export
restructureSimData<-function(db){
	idVars<-c("fit.index","FRAP","time")
	meanVars<-c("f.cyt","f.nuc")
	sdVars<-paste0("sd.",meanVars)
	return(join(
		subset(melt(db[,c(idVars,meanVars)],id.vars=idVars),!is.na(value))
		,rename(transform(subset(melt(db[,c(idVars,sdVars)],id.vars=idVars),!is.na(value)),variable=substr(variable,4,10)),c("value"="sd"))
		,by=c(idVars,"variable")))
}


#' restrucutres the simulated data for a cell pair
#'
#' @details returns a restructured data.frame more suitable for plotting
#' @param db a data.frame as returned by \code{getSimDataPair}
#' @return a restructured data.frame
#' @export
restructureSimDataPair<-function(db){
	idVars<-c("fit.index","FRAP","time","type")
	meanVars<-c("f.cyt.d","f.cyt.m","f.nuc.d","f.nuc.m")
	sdVars<-paste0("sd.",meanVars)
	return(join(
		subset(melt(db[,c(idVars,meanVars)],id.vars=idVars),!is.na(value))
		,rename(transform(subset(melt(db[,c(idVars,sdVars)],id.vars=idVars),!is.na(value)),variable=substr(variable,4,10)),c("value"="sd"))
		,by=c(idVars,"variable")))
}

#' gets the simulated model for a pair and a set of parameters
#'
#' @details wrapper over matlabs function 'nuclearFrapModel.m'
#' @param db a data.frame with the data of a cell pair (required if expData=NULL). 'time', 'FRAP', 'type', and nuclear and cytoplasmic fluorescence variables required. 
#' @param param a data.frame with the parameters of the fit (required). Must have a variable for each free parameter of the model, and a 'type' variable. 
#' @param p list of options for the FRAP experiment
#' @return a data.frame with the predicted values for the nuclear and cytoplasmic fluorescence, and a 'type' variable
#' @export
getSimDataPair<-function(db,param,p=param.FRAP.pairs){
	ut<-unique(db$type)
	output<-data.frame()
	pairExpData<-makePairEXPDATA(db,p=p)

	#FRAPs on mother
	if("mother"%in%ut & "mother"%in%param$type){
		m.param<-param[param$type=="mother",]
		m.output<-data.frame(getSimData(expData=pairExpData$mother,param=m.param,p=p),type="mother")	
		output<-data.frame(rename(m.output,c("f.nuc"="f.nuc.m","f.cyt"="f.cyt.m","sd.f.nuc"="sd.f.nuc.m","sd.f.cyt"="sd.f.cyt.m"))
												 ,f.nuc.d=NA,f.cyt.d=NA,sd.f.nuc.d=NA,sd.f.cyt.d=NA,type="mother")	
	}

	#FRAPs on daughter
	if("daughter"%in%ut & "daughter"%in%param$type){
		d.param<-param[param$type=="daughter",]
		d.output<-data.frame(getSimData(expData=pairExpData$daughter,param=d.param,p=p),type="daughter") 
		d.output<-data.frame(rename(d.output,c("f.nuc"="f.nuc.d","f.cyt"="f.cyt.d","sd.f.nuc"="sd.f.nuc.d","sd.f.cyt"="sd.f.cyt.d"))
												 ,f.nuc.m=NA,f.cyt.m=NA,sd.f.nuc.m=NA,sd.f.cyt.m=NA,type="daughter")
		if("mother"%in%ut & "mother"%in%param$type){
			output<-rbind(output,rcell2:::conform(d.output,to=output))
		}else{
			output<-d.output
		}
	}

	return(output)
}

#' gets the simulated model for a single cell and a set of parameters
#'
#' @details wrapper over matlabs function 'nuclearFrapModel.m'
#' @param db a data.frame with the data of a single cell (required if expData=NULL). 'time', 'FRAP' and nuclear and cytoplasmic fluorescence variables required. 
#' @param param a data.frame with the parameters of the fit (required). Must have a variable for each free parameter of the model. 
#' @param expData a EXPDATA list as returned by \code{makeEXPDATA}, which contains the data for a cell (required if db=NULL). 
#' @param p list of options for the FRAP experiment
#' @return a data.frame with the predicted values for the nuclear and cytoplasmic fluorescence
#' @export
getSimData<-function(db=NULL,param=NULL,expData=NULL,p=param.FRAP.pairs){
	if(is.null(param)) stop("param required")
	if(is.null(expData)&is.null(db)) return(NULL)	
	if(is.null(expData)) expData<-makeEXPDATA(db,p=p)

	param<-param[,c("fit.index",p$MODEL.PARAMETERS)]
	
	output<-
	ddply(param,.(fit.index),function(df){	

		simData<-nuclearFrapModel( EXPDATA=expData,x=as.matrix(subset(df,select=-fit.index)))

		return(as.data.frame(cbind(expData$FRAP,expData$time,simData,expData$std)))
	})
	names(output)<-c("fit.index","FRAP","time","f.nuc","f.cyt","sd.f.nuc","sd.f.cyt")
	return(output)
}
