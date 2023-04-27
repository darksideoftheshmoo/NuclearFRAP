library('stats')

#----------------------------------------------------
modelFit <- function( EXPDATA, x0, LB, UB, fixed, kEVfrac_kI_guess )
{
	# Multi start from regular grid
	startPts <- customStartPointSet( x0, fixed, guess=kEVfrac_kI_guess )

	optparam <- runOptimProblem( EXPDATA, x0, LB, UB, fixed, startPts )
	optparam <- cbind( c(1:dim(optparam)[1]), optparam )

	return( optparam )
}

#----------------------------------------------------
costFun <- function(EXPDATA,x)
{
	#running model 
	simData <- nuclearFrapModel(EXPDATA,x[1:15])

	# Calculating residuals
	resid <- (EXPDATA$y-simData) / EXPDATA$std

	# Adding penalization for geomVfrac
	geomVfrac <- x[8]
	geomVfrac_mean <- x[16]
	geomVfrac_sd <- x[17]

	# ridge regression cost
	geomVfrac_cost <- ((geomVfrac - geomVfrac_mean)^2)/(geomVfrac_sd^2)

	# Returning the residuals (for nucleus and cytoplasm)
	F <- cbind(t(resid[,1]),t(resid[,2]),geomVfrac_cost)

	return( F )
}

#----------------------------------------------------
nuclearFrapModel <- function(EXPDATA,x)
{
	#Seperating mechanistic parameters from "microscopy" parameters
	kEVfrac <- 10 ^ x[1]
	kI <- 10 ^ x[2]
	Stot0 <- x[3]
	Fn0 <- x[4]
	Fc0 <- x[5]
	PBn <- x[6]
	PBfrac <- x[7]
	geomVfrac <- 2 ^ x[8]
	auto0n <- x[9]
	auto0frac <- x[10]
	an1 <- x[11]
	an2 <- x[12]
	ac0 <- x[13]
	ac1 <- x[14]
	ac2 <- x[15]

	# Simulating the model
	Sn0 <- Stot0 * kI / (kEVfrac + kI) #assuming first points are in ss
	Sc0 <- Stot0 * kEVfrac / (kEVfrac + kI) #assuming first points are in ss
	uPBi <- unique(EXPDATA$PBindex[EXPDATA$PBindex>0])
	Sdata <- kronecker(matrix(1,sum(EXPDATA$PBindex==0),1),cbind(Sn0,Sc0))
	Fdata <- kronecker(matrix(1,sum(EXPDATA$PBindex==0),1),cbind(Fn0,Fc0))
	for( i in 1:length(uPBi)) {
		sF <- simpleFrap(EXPDATA$time[EXPDATA$PBindex==uPBi[i]],kEVfrac,kI,Sdata[nrow(Sdata),1]*(1-PBn),Sdata[nrow(Sdata),2]*(1-(PBn * PBfrac)))
		Sdata <- rbind(Sdata, sF)
		Fdata <- rbind(Fdata, kronecker( matrix( 1,sum(EXPDATA$PBindex==uPBi[i]),1), cbind(Fn0 * (1-PBn)^i,Fc0 * (1 - (PBn * PBfrac))^i)))
	}
	simData <- Sdata + Fdata

	#transforming to fluorescence levels
	simData[,1] <- (simData[,1] + auto0n) 
	simData[,2] <- (simData[,2] / geomVfrac + auto0n * auto0frac)

	#modeling multiplicative focus drift
	simData[,1] <- simData[,1] * (1 + an1 * EXPDATA$timeFocusDrift + (an2 * (EXPDATA$timeFocusDrift ^ 2)))
	simData[,2] <- simData[,2] * (ac0 + ac1 * EXPDATA$timeFocusDrift + (ac2 * (EXPDATA$timeFocusDrift ^ 2)))

	return(simData)
}

#----------------------------------------------------
simpleFrap <- function( t, kEVfrac, kI, n0, c0 )
{
	t <- t - t[1]
	Nuc <- (kI * (c0 + n0) + exp(t * (-kI-kEVfrac)) * (-(c0*kI) + n0*kEVfrac))/(kI + kEVfrac)
	Cyt <- c0 + n0 - Nuc
	NucCyt <- cbind(Nuc, Cyt)
	return(NucCyt)
}

#----------------------------------------------------
makeParamVector <- function( x, x0, fixed )
{
	x0[fixed==0] <- x
	paramVector <- x0
	return(paramVector)
}

#----------------------------------------------------
# Run porfile likelihood analysis to find a 2D contour
#
# EXPDATA
#	time
#	y
#	std
#	FRAP
#	PBindex
#	timeFocusDrift
#	sF

profileLikelihoodContour <- function( EXPDATA, x0, LB, UB, fixed, scan, cost, theta, maxLambda, precisionLambda, nStarts )
{	
	# Esto debería ser la "Floating Point Relative Accuracy". (Se puede usar eps(x=1.0) del package Pragma.)
	eps <- 2.2204e-016 

	TolX <- precisionLambda 
	getLast <- function(x) x[length(x)]
	fixedScan <- fixed | scan

	for( i in 1:dim(theta)[2] ) 
	{
		plCostTheta <- function(lambda) (getLast(profileLikelihoodCost(EXPDATA,makeParamVector(x0[scan==1] + abs(lambda) * cbind(cos(theta[i]),sin(theta[i])),x0,1-scan),LB,UB,fixedScan,nStarts)) - cost)
		extremePoint <- x0[scan==1] + abs(maxLambda) * cbind(cos(theta[i]),sin(theta[i]))
		startLambda <- maxLambda
		if( any(cbind(extremePoint > UB[scan==1], extremePoint < LB[scan==1])) ) {
			UBx0 <- (UB[scan==1] - x0[scan==1])
			UBx0[abs(UBx0)<sqrt(eps)] <- sqrt(eps)
			LBx0 <- (LB[scan==1] - x0[scan==1])
			LBx0[abs(LBx0)<sqrt(eps)] <- -sqrt(eps)
			startLambda <- 0.99 *(maxLambda / max(cbind((extremePoint - x0[scan==1]) / UBx0, (extremePoint - x0[scan==1]) / LBx0)))
		}

		lambda1 <- tryCatch(
			{
				sol <- uniroot(plCostTheta,lower=min(0,startLambda),upper=max(0,startLambda),tol=TolX)
				#alt#sol <- uniroot(plCostTheta,cbind(min(0,startLambda),max(0,startLambda)),tol=TolX)
				abs( sol$root )
			},
			error = function(cond) {
				return( -startLambda )
			}
		)
		
		temp <- cbind( x0[scan==1] + abs(lambda1) * cbind(cos(theta[i]),sin(theta[i])),theta[i],lambda1,plCostTheta(lambda1))
		if( exists('boundary') ) {
			boundary <- rbind( boundary, temp )
		} else {
			boundary <- temp
		}
	}
	
	boundary <- as.matrix(boundary)
	colnames( boundary ) <- NULL
	return( boundary )
}

#----------------------------------------------------
profileLikelihoodCost <- function( EXPDATA, x0, LB, UB, fixed, nStarts )
{
	starts <- x0[fixed==0]
	for( i in 1:(nStarts-1) ) {
		starts <- rbind( starts, runif( length(x0[fixed==0]), min=LB[fixed==0], max=UB[fixed==0] ) )
	}

	solutions <- runOptimProblem( EXPDATA, x0, LB, UB, fixed, startPts=starts )
	
	objFun <- function(x) sum( costFun(EXPDATA,makeParamVector(x,x0,fixed))^2 )
	solution <- nlminb( x0[fixed==0], objFun, gradient=NULL, hessian=NULL, lower=LB[fixed==0], upper=UB[fixed==0] )
	x <- unlist( solution$par )
	resnorm <- solution$objective

	paramCost <- FALSE
	if( (dim(solutions)[1] > 0) && (solutions[1,(dim(solutions)[2])] < resnorm) )
	{
		paramCost <- c( solutions[1,] )
	}
	else 
	{
		#print('*** Uso x0')
		paramCost <- c( makeParamVector(x,x0,fixed), resnorm )
		#if( !(solution$convergence == 0) ) {
		#	print(sprintf('*** Uso x0, que no converge, con costo: %f',resnorm))
		#}
	}
	
	paramCost <- t( as.matrix( paramCost ) )
	
	return( paramCost )
}

#----------------------------------------------------
profileLikelihoodGrid <- function( EXPDATA, x0, LB, UB, fixed, scan, gridParam, nStarts )
{
	TolFun <- 1e-2
	TolX <- 1e3

	#%Profile likelihood on scan gridParam
	n <- length(gridParam)
	startPts <- matrix( rep( x0, n ), nrow=n, byrow=TRUE)
	startPts[,scan==1] <- gridParam

	fixedScan <- fixed | scan
	#rm('profL')

	nRows <- dim(startPts)[1]
	nCols <- length(x0)+1
	profL <- matrix( data=NA_real_, nrow=nRows, ncol=nCols)
	for( i in 1:dim(startPts)[1]) 
	{
		starts <- x0[fixedScan==0]
		for( j in 1:(nStarts-1) ) {
			starts <- rbind( starts, runif( length(x0[fixedScan==0]), min=LB[fixedScan==0], max=UB[fixedScan==0] ) )
		}
		solutions <- runOptimProblem( EXPDATA, startPts[i,], LB, UB, fixedScan, startPts=starts )

		if( dim(solutions)[1] > 0 ) {
			x <- solutions[1,1:(dim(solutions)[2]-1)]
			fx <- solutions[1,(dim(solutions)[2])]
			profL[i,] <- c( makeParamVector( x[fixedScan==0], startPts[i,], fixedScan ), fx )
		}
		else 
		{
			z = rep( 0, length(fixedScan)-sum(fixedScan) )
			profL[i,] <- c( makeParamVector( z, startPts[i,], fixedScan ), -1 )
		}
	}

	colnames( profL ) <- NULL
	return( profL )
}

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
runOptimProblem <- function( EXPDATA, x0, LB, UB, fixed, startPts )
{
	upper <- UB[fixed==0]
	lower <- LB[fixed==0]

	objFun <- function(x) sum( costFun(EXPDATA,makeParamVector(x,x0,fixed))^2 )

	d <- length( x0 )
	solutions <- matrix( nrow=nrow(startPts), ncol=d+1 )
	for( i in 1:nrow(solutions) ) 
	{
		xi <- startPts[i,]
		solution <- nlminb( xi, objFun, gradient=NULL, hessian=NULL, lower=lower, upper=upper )
		xSol <- unlist( solution$par )
		vSol <- unlist( solution$objective )
		solutions[i,] <- c( makeParamVector( xSol, x0, fixed ), vSol )
	}

	# ordena
	optparam <- solutions[ sort.list( solutions[,ncol(solutions)], decreasing=FALSE ), ]

	# filtra
	optparam <- filterSolutions( optparam, 1e-2, 1e3 )

	return( optparam )
}

#----------------------------------------------------
filterSolutions <- function( solutions, TolFun, TolX ){

	norm <- function(x) sqrt(sum(x^2))
	iCost <- ncol(solutions)
	minF <- sqrt( solutions[1,iCost] )
	minX <- solutions[1,2:(iCost-1)]
	tolF <- TolFun * max( 1, abs( minF ) )
	tolX <- TolX * max( 1, norm( minX ) )

	#print( minF )
	#print( norm(minX) )
	#print( t( apply( solutions, 1, function(x) c(abs(sqrt(x[iCost])-minF) ,norm(x[2:(iCost-1)]-minX)) )))

	filterF <- function(x) (abs(sqrt(x[iCost])-minF) <= tolF)
	filterX <- function(x) (norm(x[2:(iCost-1)]-minX) <= tolX)
	filter <- function(x) (filterF(x) && filterX(x))
	solutions <- matrix( solutions[ apply( solutions, 1, filter ), ], ncol=iCost )

	return( solutions )
}

#----------------------------------------------------
customStartPointSet <- function( x0, fixed, guess )
{
	n <- dim( guess )[1]
	startPts <- matrix( rep( x0[fixed==0], n ), nrow=n, byrow=TRUE )
	startPts[,c(1,2)] <- guess
	
	return( startPts )
}
