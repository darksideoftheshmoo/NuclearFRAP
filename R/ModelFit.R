#' runs the model fit to the data, for a cell pair
#'
#' @details wrapper over matlabs function 'modelFit.m'
#' @param db a data.frame with the data of a single pair. 'time', 'FRAP', 'type', nuclear and cytoplasmic fluorescence variables for mother and daughter required
#' @param p list of options for the FRAP experiment
#' @param fixed a named integer vector indicating if a given parameter is fixed (1) or not(0). Default is p$FIXED.PARAMETERS
#' @return a data.frame of fitted parameters, with a 'type' variable
#' @export
runModelFitPair<-function(db,start.kEVfrac.kI=NULL,p=param.FRAP.pairs,fixed=p$FIXED.PARAMETERS){
	ut<-unique(db$type)
	output<-data.frame()
	pairExpData<-makePairEXPDATA(db,p=p)

	output<-ldply(pairExpData,function(ed) data.frame(runModelFit(expData=ed,start.kEVfrac.kI=start.kEVfrac.kI,p=p)))
	output<-rename(output,c(".id"="type"))

	return(output)
}

#' runs the model fit to the data, for a single cell
#'
#' @details wrapper function over matlabs function 'modelFit.m'
#' @param db a data.frame with the data of a single cell (required if expData=NULL). 'time', 'FRAP', nuclear and cytoplasmic fluorescence variables required. 
#' @param start.kEVfrac.kI a matrix of Nx2 containing the initial guesses for kEVfrac and kI
#' @param expData a EXPDATA list as returned by \code{makeEXPDATA}, which contains the data for a cell (required if db=NULL). 
#' @param p list of options for the FRAP experiment
#' @param fixed a named integer vector indicating if a given parameter is fixed (1) or not(0). Default is p$FIXED.PARAMETERS
#' @return a data.frame of fited parameters
#' @export
runModelFit<-function(db=NULL,start.kEVfrac.kI=NULL,expData=NULL,p=param.FRAP.pairs,fixed=p$FIXED.PARAMETERS){
	if(is.null(expData)&is.null(db)) return(NULL)	
	#defining start kEVfrac and kI if not in argument
	if(is.null(start.kEVfrac.kI)) start.kEVfrac.kI<-as.matrix(expand.grid(log10_kEVfrac=seq(-2,1,1),log10_kI=seq(-2,1,1)))
	#creating exp data if not passed as argument
	if(is.null(expData)) expData<-makeEXPDATA(db,p=p)
	ib<-calcInitBound(expData,p)

	output_new <- try( modelFit( EXPDATA=expData
						,x0=ib$x0, LB=ib$LB, UB=ib$UB, fixed=t(fixed)
						,kEVfrac_kI_guess = start.kEVfrac.kI ))

	if(class(output_new)=="try-error"){
		output_new <- data.frame()
	} else {
        output_new<-data.frame(output_new)
	}

	output <- output_new

	if (length(output)==0) output <- data.frame(t(rep(NA,times=length(ib$names)+2)))
	names(output)<-c("fit.index",ib$names,"cost")
	return(output)
}

#' creates a 'Pair' function from a 'cell' function
#'
#' @details transforms a function that operates over a single cell, to a function that operates over a mother-daughter pair
#' @param FUN the single cell function to be transformed to a pair function
#' @param db either a data.frame with a 'type' column, or a list with 'mother' and 'daughter' elements containing the EXPDATA objects. 
#' @param param a data.frame with the parameters of the fit (required). Must have 'type' vatiable and a variable for each free parameter of the model. 
#' @param p list of options for the FRAP experiment
#' @param ... further arguments for FUN
#' @return a data.frame, resulting from combining both calls to FUN, with a 'type' column
#' @export
pairFun<-function(FUN,db,param,p=param.FRAP.pairs,...){
	output<-data.frame()
	if(class(db)=="list"){
		ut<-names(db)
		pairExpData<-db
	}else if(class(db)=="data.frame"){
		ut<-unique(db$type)
		pairExpData<-makePairEXPDATA(db,p=p)
	}else{
		stop("db should be a list (expDataPair) or a data.frame")
	}	

	#FRAPs on mother
	if("mother"%in%ut & "mother"%in%param$type){
		m.param<-param[param$type=="mother",setdiff(names(param),"type")]
		output<-data.frame(FUN(expData=pairExpData$mother,param=m.param,p=p,...),type="mother")	
	}

	#FRAPs on daughter
	if("daughter"%in%ut & "daughter"%in%param$type){
		d.param<-param[param$type=="daughter",setdiff(names(param),"type")]
		output<-rbind(output,data.frame(FUN(expData=pairExpData$daughter,param=d.param,p=p,...),type="daughter"))	
	}

	return(output)
}

#' calculate initial values of parameters and upper/lower bounds
#' 
#' This funcion has build-in estimations that can be overwriten either by p$VALUE.PARAMETERS or by EXPDATA$expParam. 
#' if 'tot0_mean' and 'tot0_sd' are in p$MODEL.PARAMETERS, and have NA values in p$VALUE.PARAMETERS, they are calculated from the fluorescence variables of the cell and geomVfrac
#' 
#' @param EXPDATA list as returned by \code{makeEXPDATA}, which contains the data for a cell (required if db=NULL). 
#' @param  p list of options for the FRAP experiment
#' @return a list with elements x0, LB, UB and names
#' @export
calcInitBound<-function(EXPDATA,p=param.FRAP.pairs){

	#checking that model parameters in p are consistent
	if(!setequal(p$MODEL.PARAMETERS,names(p$FIXED.PARAMETERS))) stop("p$FIXED.PARAMETERS do not coincide with p$MODEL.PARAMETERS")
	if(!setequal(p$MODEL.PARAMETERS,names(p$VALUE.PARAMETERS))) stop("p$FIXED.PARAMETERS do not coincide with p$VALUE.PARAMETERS")
	if(!setequal(p$MODEL.PARAMETERS,names(p$LOWER.PARAMETERS))) stop("p$FIXED.PARAMETERS do not coincide with p$LOWER.PARAMETERS")
	if(!setequal(p$MODEL.PARAMETERS,names(p$UPPER.PARAMETERS))) stop("p$FIXED.PARAMETERS do not coincide with p$UPPER.PARAMETERS")

	#estimating initial guess
	nFrames <- length(EXPDATA$time)
	mif <- apply(EXPDATA$y[EXPDATA$PBindex==0,],2,mean) #mean initial fluorescence 
	mff <- apply(EXPDATA$y[(nFrames-5):(length(EXPDATA$time)),],2,mean) #mean final fluorescence 
	mif_sd <- apply(EXPDATA$std[EXPDATA$PBindex==0,],2,mean) #mean initial fluorescence standard deviation
	PB.tFrame<-which(EXPDATA$PBindex[2:(nFrames)]-EXPDATA$PBindex[1:(nFrames-1)]==1)
	n.bff <- EXPDATA$y[PB.tFrame,1] #before frap fluorescence 
	n.bff_sd <- EXPDATA$std[PB.tFrame,1] #before frap fluorescence standard deviation
	n.aff <- EXPDATA$y[PB.tFrame+1,1] #after frap fluorescence 
	n.aff_sd <- EXPDATA$std[PB.tFrame+1,1] #after frap fluorescence standard deviation
	c.bff <- EXPDATA$y[PB.tFrame,2] #before frap fluorescence 
	c.bff_sd <- EXPDATA$std[PB.tFrame,2] #before frap fluorescence standard deviation
	c.aff <- EXPDATA$y[PB.tFrame+1,2] #after frap fluorescence 
	c.aff_sd <- EXPDATA$std[PB.tFrame+1,2] #after frap fluorescence standard deviation

	PBn <- mean((n.bff - n.aff) / n.bff)
	PBfrac <- mean((c.bff - c.aff) / c.bff)
	auto1n <- mif[1] * 0.1
	autoPBn <- 1
	auto0n <- 0	
	an1 <- 0
	an2 <- 0
	auto1c <- 0
	autoPBc <- 0
	auto0c <- mif[2] * 0.1
	ac0 <- 1
	ac1 <- 0
	ac2 <- 0
	deltaPBn1 <- 0
	log2_geomVfrac <- 3
	if(!is.na(p$VALUE.PARAMETERS["log2_geomVfrac"])) log2_geomVfrac <- p$VALUE.PARAMETERS["log2_geomVfrac"]
	log2_geomVfrac_sd <- log2_geomVfrac * p$GEOMVFRAC.REL.SD
	if(!is.null(EXPDATA$expParam)){
		#using experimental geomVfrac if available
		if("log2_geomVfrac" %in% names(EXPDATA$expParam)) log2_geomVfrac <- EXPDATA$expParam["log2_geomVfrac"]
		#using experimental geomVfrac_sd if available
		if("log2_geomVfrac_sd" %in% names(EXPDATA$expParam)) log2_geomVfrac_sd <- EXPDATA$expParam["log2_geomVfrac_sd"]
	}

	#calculating tot0 based on geomVfrac
	tot0 <- mif[1] + 2^log2_geomVfrac * mif[2]
	tot0_sd <- mif_sd[1] + 2^log2_geomVfrac * mif_sd[2] + 2^log2_geomVfrac * log2_geomVfrac_sd * mif[2]
	log10_kEVfrac <- 0
	log10_kI <- 0

	#creating vectos with best guess and upper and lower bounds
	#kEVfrac, kI, tot0, PBn, PBfrac, geomVfrac, auto1n, autoPBn, auto0n, an1, an2, auto1c, autoPBc, auto0c, ac0, ac1, ac2
	x0 <- c(0,0 ,tot0 , PBn, PBfrac, log2_geomVfrac, auto1n, autoPBn, auto0n, an1, an2, auto1c, autoPBc, auto0c, ac0, ac1, ac2, deltaPBn1, tot0, tot0_sd)
	UB <- c( 1, 1,tot0 + p$TOT0.SIGMA.FACTOR * tot0_sd,1,1,log2_geomVfrac + p$GEOMVFRAC.SIGMA.FACTOR * log2_geomVfrac_sd,mif[1],1,mif[1], 0.1,  0.001, mif[2], 1, mif[2], 1.5, 0.1, 0.001,1, NA, NA)
	LB <- c(-3,-3,tot0 - p$TOT0.SIGMA.FACTOR * tot0_sd,0,0,log2_geomVfrac - p$GEOMVFRAC.SIGMA.FACTOR * log2_geomVfrac_sd,     0,0,     0,-0.1, -0.001,      0, 0,      0, 0.5,-0.1,-0.001,0, NA, NA)
	varNames <-c("log10_kEVfrac", "log10_kI", "tot0", "PBn", "PBfrac", "log2_geomVfrac", "auto1n", "autoPBn", "auto0n", "an1", "an2", "auto1c", "autoPBc", "auto0c", "ac0", "ac1", "ac2", "deltaPBn1", "tot0_mean", "tot0_sd")
	names(x0) <- varNames
	names(UB) <- varNames
	names(LB) <- varNames

	#reordering x0, UB and LB based on p$MODEL.PARAMETERS and intial guesses
	x0 <- sapply(p$MODEL.PARAMETERS,function(x){y <- try(eval(as.name(x),envir=as.list(x0)),silent=TRUE); ifelse(class(y)=="try-error",NA,y)})
	UB <- sapply(p$MODEL.PARAMETERS,function(x){y <- try(eval(as.name(x),envir=as.list(UB)),silent=TRUE); ifelse(class(y)=="try-error",NA,y)})
	LB <- sapply(p$MODEL.PARAMETERS,function(x){y <- try(eval(as.name(x),envir=as.list(LB)),silent=TRUE); ifelse(class(y)=="try-error",NA,y)})

	#replacing x0, UB, LB by predifined values when available
	x0 <- ifelse(is.na(p$VALUE.PARAMETERS),x0,p$VALUE.PARAMETERS)
	UB <- ifelse(is.na(p$UPPER.PARAMETERS),UB,p$UPPER.PARAMETERS)
	LB <- ifelse(is.na(p$LOWER.PARAMETERS),LB,p$LOWER.PARAMETERS)	

	#replacing x0, LB, UB by experimental parameters 
	if(!is.null(EXPDATA$expParam)){
		UBindex <- pmatch("UB_",names(EXPDATA$expParam),duplicates.ok=TRUE)
		LBindex <- pmatch("LB_",names(EXPDATA$expParam),duplicates.ok=TRUE)	
		X0index <- setdiff(seq_along(EXPDATA$expParam),c(UBindex,LBindex))
		x0[which(p$MODEL.PARAMETERS %in% names(EXPDATA$expParam))] <- EXPDATA$expParam[X0index]
		UB[which(paste0("UB_",p$MODEL.PARAMETERS) %in% names(EXPDATA$expParam))] <- EXPDATA$expParam[UBindex]
		LB[which(paste0("LB_",p$MODEL.PARAMETERS) %in% names(EXPDATA$expParam))] <- EXPDATA$expParam[LBindex]
	}

	return(list(x0=t(x0),LB=t(LB),UB=t(UB),names=p$MODEL.PARAMETERS))
}

#' Makes a list with a pair of EXPDATA 
#'
#' @param db a data.frame with the data of a single pair. 'time', 'FRAP', 'type', nuclear and cytoplasmic fluorescence variables for mother and daughter required
#' @param  p list of options for the FRAP experiment
#' @return a list with a 'mother' and a 'daughter' element, each of which contain the correspondent EXPDATA
#' @export
makePairEXPDATA <- function(db,p=param.FRAP.pairs){
	ut<-unique(db$type)
	if(p$SIGMA.METHOD=="SINGLE.CELL.PAIR"){
		if(!all(c("mother","daughter")%in%ut)) stop("Complete pair required for SIGMA.METHOD=SINGLE.CELL.PAIR")
		p$SIGMA.FUN<-I
		p$SIGMA.METHOD<-"SIGMA.FUN"

		m.expData<-makeEXPDATA(db[db$type=="mother",c("pos","time","f.nuc.m","f.cyt.m","FRAP","type")],p=p)
		d.expData<-makeEXPDATA(db[db$type=="daughter",c("pos","time","f.nuc.d","f.cyt.d","FRAP","type")],p=p)

		#crossData has the variables of the mother in the daughter FRAPs and vice-versa
		m.crossData<-makeEXPDATA(db[db$type=="daughter",c("pos","time","f.nuc.m","f.cyt.m","FRAP","type")],p=p)
		d.crossData<-makeEXPDATA(db[db$type=="mother",c("pos","time","f.nuc.d","f.cyt.d","FRAP","type")],p=p)

		#calculating sd for each variable from the "unused" data 
		#standard deviation from the "loess" model calculated by makeEXPDATA with p$SIGMA.FUN<-I
		m.sd<-matrix(apply(m.crossData$y-m.crossData$std,2,sd),nrow=dim(m.expData$std)[1],ncol=2,byrow=TRUE)
		d.sd<-matrix(apply(d.crossData$y-d.crossData$std,2,sd),nrow=dim(d.expData$std)[1],ncol=2,byrow=TRUE)

		#calculating mean value to use in rescaling
		m.mean<-matrix(apply(m.crossData$y,2,mean),nrow=dim(m.expData$std)[1],ncol=2,byrow=TRUE)
		d.mean<-matrix(apply(d.crossData$y,2,mean),nrow=dim(d.expData$std)[1],ncol=2,byrow=TRUE)

		#rescaling by fluorescence level (m.expData$std is just the smoothed verion of the data p$SIGMA.FUN<-I)
		I1<-function(x) ifelse(x>=1, x, 1)
		meanRep<-function(x) rep(mean(x),times=length(x))

		m.expData$std<-I1(meanRep(m.expData$std/m.mean)) * m.sd	
		d.expData$std<-I1(meanRep(d.expData$std/d.mean)) * d.sd	

		# m.expData$std<-sqrt(m.expData$std/m.mean) * m.sd	
		# d.expData$std<-sqrt(d.expData$std/d.mean) * d.sd	

		# m.expData$std<-(m.expData$std/m.mean) * m.sd	
		# d.expData$std<-(d.expData$std/d.mean) * d.sd	

	} else if(p$SIGMA.METHOD=="SIGMA.FUN")	{
		#FRAPs on mother
		m.expData <- NULL
		if("mother"%in%ut)
			m.expData<-makeEXPDATA(db[db$type=="mother",c("pos","time","f.nuc.m","f.cyt.m","FRAP","type")],p=p)

		#FRAPs on daughter
		d.expData <- NULL
		if("daughter"%in%ut)
			d.expData<-makeEXPDATA(db[db$type=="daughter",c("pos","time","f.nuc.d","f.cyt.d","FRAP","type")],p=p)
	}

	return(list(mother=m.expData,daughter=d.expData))
}

#' Makes a EXPDATA list, as required by matlab's functions
#'
#' @param db a data.frame with the data of a single cell. 'time', 'FRAP', nuclear and cytoplasmic fluorescence variables required
#' @param nuc string that defines the pattern that identifies the nuclear fluorescence variable
#' @param cyt string that defines the pattern that identifies the cytoplasmic fluorescence variable
#' @param  p list of options for the FRAP experiment
#' @return a EXPDATA list with elements time, y, std, FRAP, PBindex, timeFocusDrift and sF
#' @export
makeEXPDATA<-function(db,nuc="f.nuc*",cyt="f.cyt*",sigma.prefix="sd.",p=param.FRAP.pairs){
	if(nrow(db)==0) stop("no data to create EXPDATA")

	#finding nuc and cyt variables
	cyt.index<-grep(glob2rx(cyt),names(db))
	nuc.index<-grep(glob2rx(nuc),names(db))
	if(length(cyt.index)!=1) stop(length(cyt.index)," matches of cyt pattern to db names")
	if(length(nuc.index)!=1) stop(length(nuc.index)," matches of cyt pattern to db names")
	cyt<-names(db)[cyt.index]
	nuc<-names(db)[nuc.index]

	if(!"FRAP"%in%names(db)) stop("FRAP variable required in db")
	if(p$SIGMA.METHOD=="SIGMA.FUN")	{
		db <- ddply(db,.(FRAP),function(df){
			paf<-1:p$PRE.ACTIVATION.FRAMES #pre activation frames	
			df[paf,"sigma.nuc"] <- rep(p$SIGMA.FUN(mean(df[1:p$PRE.ACTIVATION.FRAMES,nuc])),p$PRE.ACTIVATION.FRAMES)  
			df[paf,"sigma.cyt"] <- rep(p$SIGMA.FUN(mean(df[1:p$PRE.ACTIVATION.FRAMES,cyt])),p$PRE.ACTIVATION.FRAMES)  
			paf<-(p$PRE.ACTIVATION.FRAMES+1):(dim(df)[1]) #post activation frames
			df[paf,"sigma.nuc"] <- p$SIGMA.FUN(predict(loess(as.formula(paste0(nuc,"~time")),data=df[paf,])))
			df[paf,"sigma.cyt"] <- p$SIGMA.FUN(predict(loess(as.formula(paste0(cyt,"~time")),data=df[paf,])))
			return(df)
		})
	} else if(p$SIGMA.METHOD=="FROM.DATA"){
		sigma.cyt.index<-grep(glob2rx(paste0(sigma.prefix,cyt)),names(db))
		sigma.nuc.index<-grep(glob2rx(paste0(sigma.prefix,nuc)),names(db))
		if(length(sigma.cyt.index)!=1) stop(length(sigma.cyt.index)," matches of sigma.prefix-cyt pattern to db names")
		if(length(sigma.nuc.index)!=1) stop(length(sigma.nuc.index)," matches of sigma.prefix-cyt pattern to db names")
		db$sigma.cyt <- db[,sigma.cyt.index]
		db$sigma.nuc <- db[,sigma.nuc.index]
	}

	db <- transformBy(db,.(FRAP),t.index=seq_along(FRAP), timeFocusDrift=min(time))
	db <- transform(db,PBindex=(t.index>p$PRE.ACTIVATION.FRAMES) + (FRAP - min(FRAP))
					  ,timeFocusDrift=timeFocusDrift-min(timeFocusDrift))

	#Experimentally determined parameters
	expParam <- NULL
	if(!is.null(p$EXPERIMENTAL.PARAMETERS.DB)){
		if("type" %in% names(db)){ #pairs
			colIndex <- which(names(p$EXPERIMENTAL.PARAMETERS.DB) %in% c(
								paste0(p$MODEL.PARAMETERS,"_",unique(db$type))
								,paste0("LB_",p$MODEL.PARAMETERS,"_",unique(db$type))
								,paste0("UB_",p$MODEL.PARAMETERS,"_",unique(db$type))))
			expParam <- unlist(subset(p$EXPERIMENTAL.PARAMETERS.DB,pos == unique(db$pos),select=colIndex))
			names(expParam) <- substr(names(expParam),1,nchar(names(expParam)) - nchar(as.character(unique(db$type))) - 1)
		} else { #single cell
			expParam <- unlist(subset(p$EXPERIMENTAL.PARAMETERS.DB,pos == unique(db$pos),select=-pos))
		}
	}

	EXPDATA <- list(time=		db$time
				    ,y=			cbind(db[[nuc]],db[[cyt]])
					,std=		cbind(db[["sigma.nuc"]],db[["sigma.cyt"]])
					,FRAP=		db$FRAP
					,PBindex= 	db$PBindex
					,timeFocusDrift= 	db$timeFocusDrift)
					#,sF= p$TOT0.SIGMA.FACTOR
					#,expParam = expParam)
	if(!is.null(expParam)) EXPDATA[["expParam"]] <- expParam

	return(EXPDATA)
}
