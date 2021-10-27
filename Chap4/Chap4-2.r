rm( list = ls(  ) )
set.seed( 123456 )
source( "../MCMC.r" )

fLogN <- function( dX, dMu, dSigma2 )
{
	return( - 0.5 * ( dX - dMu ) ^ 2 / dSigma2 )
}

fLogIG <- function( dX, dNu, dLamb )
{
	return( - ( dNu + 1.0 ) * log( dX ) - dLamb / dX )
}

fLogL <- function( vTheta, vX, vN )
{
	dMu     <- vTheta[1]
	dSigma2 <- vTheta[2]
	dSigma <- sqrt( dSigma2 )
	
	cK <- length( vX )
	vZ <- ( log( vX ) - dMu ) / dSigma

	dLogL <- 0.0
	for ( i in 1 : ( cK + 1 ) )
	{
		if ( i == 1 )
		{
			dLogL <- log( dnorm( vZ[i] ) ) - log( dSigma ) + ( vN[1] - 1 ) * log( pnorm( vZ[i] ) )
		}
		else if ( i < cK + 1 )
		{
			dLogL <- dLogL + log( dnorm( vZ[i] ) ) - log( dSigma ) + ( vN[i] - vN[( i - 1 )] - 1 ) * log( pnorm( vZ[i] ) - pnorm( vZ[( i - 1 )] ) )
		}
		else
		{
			dLogL <- dLogL + ( vN[i] - vN[( i - 1 )] ) * log( 1.0 - pnorm( vZ[( i - 1 )] ) )
		}
	}

	return( dLogL )
}

## set timer
ttt <- proc.time(  )

## load data
sData <- read.csv( "../Data/Income.csv" )
vX <- sData$I2005 / 1000
vN <- seq( 1000, 10000, 1000 )

## set prior
dMu.0   <- 0.0
dTau2.0 <- 100.0
dNu.0   <- 2.0
dLamb.0 <- 0.01

## set initial values
vTheta <- c( 8.4, 0.4 )

## set biru-in & draws
iBurn <- 20000
iDraw <- 20000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
res <- optim( vTheta, fLogL, method = "L-BFGS-B", lower = c( -Inf, 0.0 ), hessian = TRUE, control = list( fnscale=-1 ), vX = vX, vN = vN )
mSigma <- t( chol( solve( -res$hessian ) ) )
vTheta <- res$par

dC      <- 0.5
dCount  <- 0
dAccept <- 0

## set posterior matrices
mPostParam.05 <- matrix( 0, iDraw, 2 )
vPostGini.05  <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dC = ", sprintf( "%7.3f", dC ), sep = "" ) )
	}
	dCount <- dCount + 1

	## draw vTheta
	vTheta.n <- vTheta + dC * mSigma %*% rnorm( 2 )
	if ( vTheta.n[2] > 0.0 )
	{
		if ( runif( 1 ) < exp( fLogN( vTheta.n[1], dMu.0, dTau2.0 ) + fLogIG( vTheta.n[2], 0.5 * dNu.0, 0.5 * dLamb.0 ) + fLogL( vTheta.n, vX, vN ) - fLogN( vTheta[1], dMu.0, dTau2.0 ) - fLogIG( vTheta[2], 0.5 * dNu.0, 0.5 * dLamb.0 ) - fLogL( vTheta, vX, vN ) ) )
		{
			vTheta  <- vTheta.n
			dAccept <- dAccept + 1
		}
	}
	if ( iter <= iBurn )
	{
		if ( dAccept / dCount < 0.4 )
		{
			dC <- dC / 1.1
		}
		if ( dAccept / dCount > 0.6 )
		{
			dC <- 1.1 * dC
		}
	}
	if ( iter == iBurn )
	{
		dCount  <- 0
		dAccept <- 0
	}

	## save parameters
	if ( iter > iBurn	)
	{
		it <- iter - iBurn

		mPostParam.05[it, ] <- vTheta
		vPostGini.05[it]    <- 2.0 * pnorm( sqrt( vTheta[2] / 2.0 ) ) - 1
	}
}

print( paste( "acceptance rate = ", dAccept / dCount, sep = "" ) )

## set parameter names
sName <- c( expression( mu ), expression( sigma^2 ) )

## Posterior analysis
mSummary      <- fSummary( mPostParam.05, sName, "LN05" )
vSummary.Gini <- fSummary( vPostGini.05, "Gini", "Gini05" )
fDraw( mPostParam.05, sName, "LN05" )
fDraw( vPostGini.05, "Gini", "Gini05" )

print( paste( "Execution time: ", ( proc.time(  ) - ttt )[3], sep = "" ) )

## load data
vX <- sData$I2010 / 1000

## set initial values
vTheta <- c( 8.4, 0.4 )

## set some matrices, vectors and scalars
res <- optim( vTheta, fLogL, method = "L-BFGS-B", lower = c( -Inf, 0.0 ), hessian = TRUE, control = list( fnscale=-1 ), vX = vX, vN = vN )
mSigma <- t( chol( solve( -res$hessian ) ) )
vTheta <- res$par

dC      <- 0.5
dCount  <- 0
dAccept <- 0

## set posterior matrices
mPostParam.10 <- matrix( 0, iDraw, 2 )
vPostGini.10  <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dC = ", sprintf( "%7.3f", dC ), sep = "" ) )
	}
	dCount <- dCount + 1

	## draw vTheta
	vTheta.n <- vTheta + dC * mSigma %*% rnorm( 2 )
	if ( vTheta.n[2] > 0.0 )
	{
		if ( runif( 1 ) < exp( fLogN( vTheta.n[1], dMu.0, dTau2.0 ) + fLogIG( vTheta.n[2], 0.5 * dNu.0, 0.5 * dLamb.0 ) + fLogL( vTheta.n, vX, vN ) - fLogN( vTheta[1], dMu.0, dTau2.0 ) - fLogIG( vTheta[2], 0.5 * dNu.0, 0.5 * dLamb.0 ) - fLogL( vTheta, vX, vN ) ) )
		{
			vTheta  <- vTheta.n
			dAccept <- dAccept + 1
		}
	}
	if ( iter <= iBurn )
	{
		if ( dAccept / dCount < 0.4 )
		{
			dC <- dC / 1.1
		}
		if ( dAccept / dCount > 0.6 )
		{
			dC <- 1.1 * dC
		}
	}
	if ( iter == iBurn )
	{
		dCount  <- 0
		dAccept <- 0
	}

	## save parameters
	if ( iter > iBurn	)
	{
		it <- iter - iBurn

		mPostParam.10[it, ] <- vTheta
		vPostGini.10[it]    <- 2.0 * pnorm( sqrt( vTheta[2] / 2.0 ) ) - 1
	}
}

print( paste( "acceptance rate = ", dAccept / dCount, sep = "" ) )

## set parameter names
sName <- c( expression( mu ), expression( sigma^2 ) )

## Posterior analysis
mSummary      <- fSummary( mPostParam.10, sName, "LN10" )
vSummary.Gini <- fSummary( vPostGini.10, "Gini", "Gini10" )
fDraw( mPostParam.10, sName, "LN10" )
fDraw( vPostGini.10, "Gini", "Gini10" )

print( paste( "Execution time: ", ( proc.time(  ) - ttt )[3], sep = "" ) )

## draw figure
mG.05 <- density( vPostGini.05 )
mG.10 <- density( vPostGini.10 )
vG    <- density( vPostGini.10 - vPostGini.05 )
y.max <- max( mG.05$y, mG.10$y )
x.min <- min( mG.05$x, mG.10$x )
x.max <- max( mG.05$x, mG.10$x )
pdf( file = "gini-comparison.pdf", width = 8, height = 5 )
par( mfrow = c( 1, 2 ), mar = c( 2, 2, 0.5, 0.5 ), family = "Japan1Ryumin" )
plot( mG.05, lty = 1, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "", xlab = "", ylab = "", bty = "n" )
par( new = T )
plot( mG.10, lty = 2, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n" )
legend( "topright", legend = c( "2005年", "2010年" ), lty = c( 1, 2 ), bty = "n" )
plot( vG, lty = 1, main = "", xlab = "", ylab = "", bty = "n" )
dev.off(  )

print( sum( vPostGini.10 - vPostGini.05 > 0.0 ) / iDraw )
