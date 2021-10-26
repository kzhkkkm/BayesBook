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

## simulate data
sData <- SimData( 300 )
vY    <- sData$Y
pdf( file = "ar1-tsplot.pdf", width = 7, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
ts.plot( vY, main = "", xlab = "", ylab = "",  gpars = list( bty = "n" ) )
dev.off(  )

## set timer
ttt <- proc.time(  )

## set dimensions
cT <- length( vY )

## set hyper-parameters
dAlpha.0 <- 1.0
dBeta.0  <- 1.0

## set initial parameters
dPhi    <- 0.0

## set biru-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dC      <- 1.0
dCount  <- 0
dAccept <- 0
vC      <- rep( 0.0, iBurn )

## set posterior matrices
vPostParam <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ), "dC = ", sprintf( "%5.3f", dC ), sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dPhi
	dPhi.n <- dPhi + dC * rnorm( 1 )
	if ( abs( dPhi.n ) < 1.0 )
	{
		dLogL.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )
		dLogL.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
		if ( runif( 1 ) <= exp( dLogL.n - dLogL.o ) )
		{
			dPhi    <- dPhi.n
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
	if ( iter - iBurn == 0 )
	{
	  dCount  <- 0
	  dAccept <- 0
	}

	## save parameters
	if ( iter <= iBurn )
	{
		vC[iter] <- dC
	}
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		vPostParam[it] <- dPhi
	}
}

## set parameter names
sName <- c( bquote( phi ) )

## posterior analysis
vSummary <- fSummary( vPostParam, sName, "AR1" )
print( paste( "Acceptance rate =", sprintf( "%6.3f", dAccept / dCount ), sep = "" ) )

pdf(file = "ar1-tuned.pdf", width = 6, height = 2 )
par( mfrow = c( 1, 3 ), family = "Times" )
par( mar = c( 2, 2, 1.5, 0.5 ) )
ts.plot( vPostParam, main = expression( phi ), xlab = "", ylab = "", gpars = list( bty = "n" ) )
plot( density( vPostParam ), main = expression( phi ), xlab = "", ylab = "", bty = "n" )
par( mar = c( 2, 2, 2.8, 0.5 ) )
acf( vPostParam, main = expression( phi ), xlab = "", ylab = "", bty = "n", lag = 300 )
dev.off(  )

pdf( file = "ar1-C.pdf", width = 7, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Times" )
ts.plot( vC[1 : 500], main = "",  xlab = "", ylab = "", gpars = list( bty = "n" ) )
dev.off(  )

print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## Case I (c = 0.005)
## set initial parameters
dPhi <- 0.8

## set biru-in & draws
iBurn <- 0
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dC      <- 0.005
dCount  <- 0
dAccept <- 0

## set posterior matrices
vPostParam <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ),  sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dPhi
	dPhi.n <- dPhi + dC * rnorm( 1 )
	if ( abs( dPhi.n ) < 1.0 )
	{
		dLogL.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )
		dLogL.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
		if ( runif( 1 ) <= exp( dLogL.n - dLogL.o ) )
		{
			dPhi    <- dPhi.n
			dAccept <- dAccept + 1
		}
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		vPostParam[it] <- dPhi
	}
}

mPostParam <- vPostParam
vAccept    <- round( dAccept / dCount, digits = 3 )
print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## Case II (c = 0.05)
## set initial parameters
dPhi <- 0.8

## set biru-in & draws
iBurn <- 0
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dC      <- 0.05
dCount  <- 0
dAccept <- 0

## set posterior matrices
vPostParam <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
  ## monitoring
  if ( iter %% nstep == 0 )
  {
    print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ),  sep = "" ) )
  }
  dCount <- dCount + 1
  
  ## sampling dPhi
  dPhi.n <- dPhi + dC * rnorm( 1 )
  if ( abs( dPhi.n ) < 1.0 )
  {
    dLogL.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )
    dLogL.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
    if ( runif( 1 ) <= exp( dLogL.n - dLogL.o ) )
    {
      dPhi    <- dPhi.n
      dAccept <- dAccept + 1
    }
  }
  
  ## save parameters
  if ( iter > iBurn )
  {
    it <- iter - iBurn
    
    vPostParam[it] <- dPhi
  }
}

mPostParam <- cbind( mPostParam, vPostParam )
vAccept    <- cbind( vAccept, round( dAccept / dCount, digits = 3 ) )
print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## Case III ( dC = 1.0 )
## set initial parameters
dPhi    <- 0.8

## set biru-in & draws
iBurn <- 0
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dC      <- 1.0
dCount  <- 0
dAccept <- 0

## set posterior matrices
vPostParam <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dPhi = ", sprintf( "%7.3f", dPhi ),  sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dPhi
	dPhi.n <- dPhi + dC * rnorm( 1 )
	if ( abs( dPhi.n ) < 1.0 )
	{
		dLogL.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi )
		dLogL.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vY, dPhi.n )
		if ( runif( 1 ) <= exp( dLogL.n - dLogL.o ) )
		{
			dPhi    <- dPhi.n
			dAccept <- dAccept + 1
		}
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		vPostParam[it] <- dPhi
	}
}

mPostParam <- cbind( mPostParam, vPostParam )
vAccept    <- cbind( vAccept, round( dAccept / dCount, digits = 3 ) )
print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## comparison
sName <- c( expression( italic( c ) == 0.005 ), expression( italic( c ) == 0.05 ), expression( italic( c ) == 1.0 ) )
y.max <- 0.01 * ceiling( 100 * max( mPostParam ) )
y.min <- 0.01 * floor( 100 * min( mPostParam ) )
pdf( file = "ar1-comparison.pdf", width = 7, height = 5 )
par( mfrow = c( 3, 3 ), family = "Times" )
par( mar = c( 2, 2, 1.5, 0.5 ) )
for ( j in 1 : 3 )
{
  ts.plot( mPostParam[1 : 2000, j], main = sName[j], xlab = "", ylab = "", ylim = c( y.min, y.max ), gpars = list( bty = "n" ) )
}
for ( j in 1 : 3 )
{
  plot( density( mPostParam[, j] ), main = sName[j], xlab = "", ylab = "", xlim = c( y.min, y.max ), bty = "n" )
}
par( mar = c( 2, 2, 2.8, 0.5 ) )
for ( j in 1 : 3 )
{
  acf( mPostParam[, j], main = sName[j], xlab = "", ylab = "", bty = "n", lag = 300 )
}
dev.off(  )
mSummary <- fSummary( mPostParam, sName, "comparison")
print( vAccept )
