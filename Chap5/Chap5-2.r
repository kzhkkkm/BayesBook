rm( list = ls(  ) )
set.seed( 123456 )
source( "../MCMC.r" )

## set timer
ttt <- proc.time(  )

## load data
sData <- read.csv( "../Data/mi.csv" )
vZ.0 <- sData$Y
vZ.1 <- sData$K
vZ.2 <- sData$L
vY    <- log( vZ.0 ) - log( vZ.2 )
vX.1  <- log( vZ.1 ) - log( vZ.2 )
vX.2  <- log( vZ.2 )
mX    <- cbind( 1.0, vX.1, vX.2 )

sName.P <- sData$X2009年

## set dimensions
cN <- nrow( mX )
cK <- ncol( mX )

## set hyper-parameters
vBeta.0  <- rep( 0.0, cK )
mSigma.0 <- diag( 100.0, cK )
dNu.0    <- 2.0
dLamb.0  <- 0.01
dEta.0   <- 2.0
dPsi.0   <- 0.01

## set initial values
vBeta   <- rep( 0.0, cK )
dSigma2 <- 1.0
dTau2   <- 1.0
vU      <- rep( 0.0, cN )

# set biru-in & draws
iBurn <- 10000
iDraw <- 40000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

# set scalars, vectors & matrices
mInvSigma.0 <- solve( mSigma.0 )
vISB.0      <- mInvSigma.0 %*% vBeta.0

t.mX <- t( mX )
mXX  <- t.mX %*% mX

dNu.1  <- cN + dNu.0
dEta.1 <- cN + dEta.0

## set posterior matrices
mPostParam <- matrix( 0, iDraw, cK + 2 )
mPostU     <- matrix( 0, iDraw, cN )
mPostL     <- matrix( 0, iDraw, 2 )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", vBeta[1] = ", sprintf( "%7.3f", vBeta[1] ), sep = "" ) )
	}

	## sampling vBeta
	mSigma.1 <- solve( mXX / dSigma2 + mInvSigma.0 )
	vBeta.1  <- mSigma.1 %*% ( t.mX %*% ( vY + vU )  / dSigma2 + vISB.0 )

	vBeta <- vBeta.1 + t( chol( mSigma.1 ) ) %*% rnorm( cK )

	## sampling dSigma2
	vE <- vY - mX %*% vBeta + vU
	dLamb.1 <- sum( vE ^ 2 ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## sampling dTau2
	dPsi.1 <- sum( vU ^ 2 ) + dPsi.0
	
	dTau2 <- 1.0 / rgamma( 1, 0.5 * dEta.1, 0.5 * dPsi.1 )

	## sampling vU
	dTau2.1 <- 1.0 / ( 1.0 / dSigma2 + 1.0 / dTau2 )
	vMu.1   <- dTau2.1 * ( mX %*% vBeta - vY ) / dSigma2
	for ( i in 1 : cN )
	{
		vU[i] <- rtnorm( vMu.1[i], dTau2.1, 0.0, Inf )
	}
	
	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, ] <- c( vBeta, dSigma2, dTau2 )
		mPostU[it, ]     <- exp( - vU )
		mPostL[it, 1]    <- prod( dnorm( ( vY - mX %*% vBeta + vU ) / sqrt( dSigma2 ) ) / sqrt( dSigma2 ) )
		mPostL[it, 2]    <- prod( sqrt( 2.0 ) * exp( - 0.5 * ( vY - mX %*% vBeta ) ^ 2 / ( dSigma2 + dTau2 ) ) / sqrt( pi * ( dSigma2 + dTau2 ) ) * pnorm( - sqrt( dTau2 ) * ( vY - mX %*% vBeta ) / ( sqrt( dSigma2 ) * sqrt( dSigma2 + dTau2 ) ) ) )
	}
}

## set parameter names
sName <- NULL
for ( i in 1 : ( cK + 2 ) )
{
	if ( i <= cK )
	{
		sName <- c( sName, bquote( beta[.( i )] ) )
	}
	else if ( i == cK + 1 )
	{
		sName <- c( sName, bquote( sigma^2 ) )
	}
	else
	{
		sName <- c( sName, bquote( tau^2 ) )
	}
}
sName.U <- NULL
for( i in 1 : cN )
{
	sName.U <- c( sName.U, bquote( u[.( i )] ) )
}	

## posterior analysis
mSummary   <- fSummary( mPostParam, sName, "SF-Param", vLik = mPostL[, 2] )
mSummary.U <- fSummary( mPostU, sName.U, "SF-U" )
fDraw( mPostParam, sName, "SF-Param-" )
#fDraw( mPostU, sName.U, "SF-U-" )

print( c( min( mSummary.U[,1] ), max( mSummary.U[,1] ) ) )
print( paste( "log marginal likelihood = ", log( 1.0 / apply( 1.0 / mPostL, 2, mean ) ), sep = "" ) )

pdf( file = "Inefficiency.pdf", width = 7, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2.9, 2.1, 0.5, 0.5 ), family = "Japan1Ryumin" )
boxplot( mPostU[, order( mSummary.U[,1], decreasing = T )], names = sName.P[order( mSummary.U[,1], decreasing = T )], cex = 0.5, las = 2, cex.axis = 0.6, frame = F )
dev.off(  )

print( paste( "Execution time: ", ( proc.time(  ) - ttt )[3], sep = "" ) )
