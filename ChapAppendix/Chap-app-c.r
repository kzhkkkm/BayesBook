rm( list = ls(  ) )
source( "../MCMC.r" )
set.seed( 123456 )

SimData <- function( cT )
{
	## set true parameters
	dPhi <- 0.8

	## generate vY
	vY <- rep( 0.0, cT )
	for ( tt in 1 : cT )
	{
		if ( tt == 1 )
		{
			vY[tt] <- sqrt( 1.0 / ( 1.0 - dPhi ^ 2 ) ) * rnorm( 1 )
		}
		else
		{
			vY[tt] <- dPhi * vY[( tt - 1 )] + rnorm( 1 )
		}
	}

	return( list( Y = vY ) )
}

fLogB <- function( dX, dAlpha, dBeta )
{
	dY <- 0.5 * ( dX + 1.0 )

	return( ( dAlpha - 1.0 ) * log( dY ) + ( dBeta - 1.0 ) * log( 1.0 - dY ) )
}

fLogL <- function( vY, dPhi )
{
	cT <- length( vY )

	dLogL <- 0.5 * log( 1.0 - dPhi ^ 2 ) - 0.5 * ( ( 1.0 - dPhi ^ 2 ) * vY[1] ^ 2 + sum( ( vY[2 : cT] - dPhi * vY[1 : ( cT - 1 )] ) ^ 2 ) )

	return( dLogL )
}

fLogPi.d <- function( dPhi, vY, dAlpha.0, dBeta.0 )
{
	cT <- length( vY )

	dLogL.d <- ( dAlpha.0 - 1.0 ) / ( dPhi + 1.0 ) - ( dBeta.0 - 1.0 ) / ( 1.0 - dPhi ) - dPhi / ( 1.0 - dPhi ^ 2 ) + dPhi * vY[1] ^ 2 + sum( vY[1 : ( cT - 1 )] * ( vY[2 : cT] - dPhi * vY[1 : ( cT - 1 )] ) )

	return( dLogL.d )
}

## simulate data
sData <- SimData( 300 )
vY    <- sData$Y

## set dimensions
cT <- length( vY )

## set hyper-parameters
dAlpha.0 <- 1.0
dBeta.0  <- 1.0

## set biru-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dS2    <- 1.0 / sum( vY[2 : ( cT - 1 )] ^ 2 )
dPhi.h <- dS2 * sum( vY[2 : cT] * vY[1 : ( cT - 1 )] )

## set posterior matrices
mPostParam <- matrix( 0.0, iDraw, 4 )
vTime      <- rep( 0.0, 4 )
vAccept    <- rep( 0.0, 5 )

## random-walk chain MH
## set initial parameters
dPhi <- 0.0

## set some matrices, vectors and scalars
dLogPi.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )

dC      <- 1.0
dCount  <- 0
dAccept <- 0

## set timer
ttt <- proc.time(  )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ), ", dC = ", sprintf( "%5.3f", dC ), sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dPhi
	dPhi.n <- dPhi + dC * rnorm( 1 )
	if ( abs( dPhi.n ) < 1.0 )
	{
		dLogPi.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
		if ( runif( 1 ) <= exp( dLogPi.n - dLogPi.o ) )
		{
			dPhi     <- dPhi.n
			dAccept  <- dAccept + 1
			dLogPi.o <- dLogPi.n
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
	if ( iter - iBurn == 0 )
	{
	  dCount  <- 0
	  dAccept <- 0
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, 1] <- dPhi
	}
}
vTime[1]   <- ( proc.time(  ) - ttt )[3]
vAccept[1] <- dAccept / dCount

## independence-chain MH
## set initial parameters
dPhi <- 0.0

## set some matrices, vectors and scalars
dLogPi.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )
dQ.o     <- -0.5 * ( dPhi - dPhi.h ) ^ 2 / dS2

dCount  <- 0
dAccept <- 0

## set timer
ttt <- proc.time(  )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ), sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dPhi
	dPhi.n <- rtnorm( dPhi.h, dS2, -1.0, 1.0 )
	
	dLogPi.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
	dQ.n     <- -0.5 * ( dPhi.n - dPhi.h ) ^ 2 / dS2
	if ( runif( 1 ) <= exp( dLogPi.n + dQ.o - dLogPi.o - dQ.n ) )
	{
		dPhi     <- dPhi.n
		dAccept  <- dAccept + 1
		dLogPi.o <- dLogPi.n
		dQ.o     <- dQ.n

	}
	if ( iter - iBurn == 0 )
	{
	  dCount  <- 0
	  dAccept <- 0
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, 2] <- dPhi
	}
}
vTime[2]   <- ( proc.time(  ) - ttt )[3]
vAccept[2] <- dAccept / dCount

## ARMH
## set initial parameters
dPhi <- 0.0

## set some matrices, vectors and scalars
dLogKappa <- - 0.5 * ( dS2 * sum( vY ^ 2 ) - dPhi.h ^ 2 ) / dS2
dLogPi.o  <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )	
dQ.o      <- dLogKappa - 0.5 * sqrt( dS2 ) - 0.5 * ( dPhi - dPhi.h ) ^ 2 / dS2

dCount.AR <- 0
dCount    <- 0
dAccept   <- 0

## set timer
ttt <- proc.time(  )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ), sep = "" ) )
	}

	## sampling dPhi
	dAlpha.AR <- -100.0
	while ( runif( 1 ) > dAlpha.AR )
	{
		dPhi.n    <- rtnorm( dPhi.h, dS2, -1.0, 1.0 )
		dLogPi.n  <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
		dQ.n      <- dLogKappa - 0.5 * ( dPhi.n - dPhi.h ) ^ 2 / dS2
		dAlpha.AR <- exp( dLogPi.n - dQ.n )
		dCount.AR <- dCount.AR + 1
	}
	dCount   <- dCount + 1
	if ( runif( 1 ) <= exp( dLogPi.n + min( dLogPi.o, dQ.o ) - dLogPi.o - min( dLogPi.n, dQ.n ) ) )
	{
		dPhi     <- dPhi.n
		dAccept  <- dAccept + 1
		dLogPi.o <- dLogPi.n	
		dQ.o     <-	dQ.n
	}
	if ( iter - iBurn == 0 )
	{
		dCount.AR <- 0
	  dCount    <- 0
	  dAccept   <- 0
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, 3] <- dPhi
	}
}
vTime[3]   <- ( proc.time(  ) - ttt )[3]
vAccept[3] <- dCount / dCount.AR
vAccept[4] <- dAccept / dCount

## ランジェヴァン MH
## set initial parameters
dPhi <- 0.0

## set some matrices, vectors and scalars
dLogPi.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )
dQ.o     <- -0.5 * ( dPhi - dPhi.n - 0.5 * dC ^ 2 * fLogPi.d( dPhi.n, vY, dAlpha.0, dBeta.0 ) ) ^ 2 / dC ^ 2

dC      <- 0.1
dCount  <- 0
dAccept <- 0

## set timer
ttt <- proc.time(  )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ), sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dPhi
	dPhi.n <- dPhi + 0.5 * dC ^ 2 * fLogPi.d( dPhi, vY, dAlpha.0, dBeta.0 ) + dC * rnorm( 1 )

	if ( abs( dPhi.n ) < 1.0  )
	{
		dLogPi.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
		dQ.n     <- -0.5 * ( dPhi.n - dPhi - 0.5 * dC ^ 2 * fLogPi.d( dPhi, vY, dAlpha.0, dBeta.0 ) ) ^ 2 / dC ^ 2
		if ( runif( 1 ) <= exp( dLogPi.n + dQ.o - dLogPi.o - dQ.n ) )
		{
			dPhi     <- dPhi.n
			dAccept  <- dAccept + 1
			dLogPi.o <- dLogPi.n
			dQ.o     <- dQ.n
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
	if ( iter - iBurn == 0 )
	{
	  dCount  <- 0
	  dAccept <- 0
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, 4] <- dPhi
	}
}
vTime[4]   <- ( proc.time(  ) - ttt )[3]
vAccept[5] <- dAccept / dCount


## set parameter names
sName <- c( "酔歩連鎖", "独立連鎖", "採択・棄却連鎖", "ランジェヴァン連鎖" )

## posterior analysis
vSummary <- fSummary( mPostParam, sName, "mh-comparison" )
print( vSummary )
print( vAccept )
print( vTime )

y.min <- 0.01 * floor( 100 * min( mPostParam ) )
y.max <- 0.01 * ceiling( 100 * max( mPostParam ) )
pdf( file = "mh-comparison.pdf", width = 8, height = 5 )
par( mfrow = c( 2, 2 ), mar = c( 2, 2, 3, 0.5 ), family = "Japan1Ryumin" )
for ( i in 1 : 4 )
{
  acf( mPostParam[, i], main = sName[i], xlab = "", ylab = "", bty = "n", lag = 30 )
}
dev.off(  )

x.min <- y.min
x.max <- y.max
y.max <- 0.0
for ( i in 1 : 4 )
{
	y.max <- max( y.max, density( mPostParam[, i] )$y )
}
pdf( file = "mh-density.pdf" )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
for ( i in 1 : 4 )
{
	if ( i > 1 )
	{
		par( new = T )
		plot( density( mPostParam[, i] ), main = "", xlab = "", ylab = "", lty = i, xlim = c( x.min, x.max ), ylim = c( 0.0, y.max ), bty = "n", xaxt = "n", yaxt = "n" )
	}
	else
	{
		plot( density( mPostParam[, i] ), main = "", xlab = "", ylab = "", lty = i, xlim = c( x.min, x.max ), ylim = c( 0.0, y.max ), bty = "n" )
	}
}
legend( "topright", sName, lty = 1 : 4, bty = "n" )
dev.off(  )

