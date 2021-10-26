rm( list = ls(  ) )
source( "../MCMC.r" )
set.seed( 123456 )

SimData <- function( cN )
{
	## set true parameters
	dMu     <- 2.0
	dSigma2 <- 1.5

	## generate vY
	vY <- dMu + sqrt( dSigma2 ) * rnorm( cN )

	return( list( Y = vY ) )
}

fLogN <- function( dX, dMu, dSigma2 )
{
	return( - 0.5 * ( dX - dMu ) ^ 2 / dSigma2 )
}

fLogIG <- function( dX, dNu, dLamb )
{
	return( - ( dNu + 1.0 ) * log( dX ) - dLamb / dX )
}

fLogL <- function( vY, dMu, dSigma2 )
{
	cN <- length( vY )

	dLogL <- - 0.5 * cN * log( dSigma2 ) - 0.5 * sum( ( vY - dMu ) ^ 2 ) / dSigma2

	return( dLogL )
}

## simulate data
sData <- SimData( 100 )
vY    <- sData$Y

## Gibbs sampler
## set timer
ttt <- proc.time(  )

## set dimensions
cN <- length( vY )

## set hyper-parameters
dMu.0   <- 0.0
dTau2.0 <- 100.0
dNu.0   <- 2.0
dLamb.0 <- 0.01

## set initial parameters
dMu     <- mean( vY )
dSigma2 <- var( vY )

## set biru-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dSumY <- sum( vY )
dNu.1 <- cN + dNu.0

## set posterior matrices
mPostParam <- matrix( 0, iDraw, 2 )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dMu = ", sprintf( "%7.3f", dMu ),  sep = "" ) )
	}

	## sampling dMu
	dTau2.1 <- 1.0 / ( cN / dSigma2 + 1.0 / dTau2.0 )
	dMu.1   <- dTau2.1 * ( dSumY / dSigma2 + dMu.0 / dTau2.0 )

	dMu <- dMu.1 + sqrt( dTau2.1 ) * rnorm( 1 )

	## sampling dSigma2
	dLamb.1 <- sum( ( vY - dMu ) ^ 2 ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, ] <- c( dMu, dSigma2 )
	}
}
mPostPara <- mPostParam

## set parameter names
sName <- c( bquote( mu ), bquote( sigma^2 ) )

## posterior analysis
mSummary <- fSummary( mPostParam, sName, "Normal-Gibbs" )

print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## set timer
ttt <- proc.time(  )

## set dimensions
cN <- length( vY )

## set hyper-parameters
dMu.0   <- 0.0
dTau2.0 <- 100.0
dNu.0   <- 2.0
dLamb.0 <- 0.01

## set initial parameters
dMu     <- mean( vY )
dSigma2 <- var( vY )

## set biru-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dCount  <- 0
vAccept <- rep( 0, 2 )
vC      <- rep( 0.1, 2 )

## set posterior matrices
mPostParam <- matrix( 0, iDraw, 2 )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dMu = ", sprintf( "%7.3f", dMu ),  sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dMu
	dMu.n   <- dMu + vC[1] * rnorm( 1 )
	dLogL.o <- fLogN( dMu, dMu.0, dTau2.0 ) + fLogL( vY, dMu, dSigma2 )
	dLogL.n <- fLogN( dMu.n, dMu.0, dTau2.0 ) + fLogL( vY, dMu.n, dSigma2 )
	if ( runif( 1 ) <= exp( dLogL.n - dLogL.o ) )
	{
		dMu        <- dMu.n
		vAccept[1] <- vAccept[1] + 1
	}
	if ( iter <= iBurn )
	{
		if( vAccept[1] / dCount < 0.4 )
		{
			vC[1] <- vC[1] / 1.1
		}
		if ( vAccept[1] / dCount > 0.6 )
		{
			vC[1] <- 1.1 * vC[1]
		}
	}

	## sampling dSigma2
	dSigma2.n <- dSigma2 + vC[2] * rnorm( 1 )
	if ( dSigma2.n > 0.0 )
	{
		dLogL.o <- fLogIG( dSigma2, 0.5 * dNu.0, 0.5 * dLamb.0 ) + fLogL( vY, dMu, dSigma2 )
		dLogL.n <- fLogIG( dSigma2.n, 0.5 * dNu.0, 0.5 * dLamb.0 ) + fLogL( vY, dMu, dSigma2.n )
		if ( runif( 1 ) <= exp( dLogL.n - dLogL.o ) )
		{
			dSigma2    <- dSigma2.n
			vAccept[2] <- vAccept[2] + 1
		}
	}
	if ( iter <= iBurn )
	{
		if( vAccept[2] / dCount < 0.4 )
		{
			vC[2] <- vC[2] / 1.1
		}
		if ( vAccept[2] / dCount > 0.6 )
		{
			vC[2] <- 1.1 * vC[2]
		}
	}
	if ( iter - iBurn == 0 )
	{
	  dCount  <- 0
	  vAccept <- rep( 0, 2 )
	}
	
	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, ] <- c( dMu, dSigma2 )
	}
}
mPostPara <- cbind( mPostPara, mPostParam )

## set parameter names
sName <- c( bquote( mu ), bquote( sigma^2 ) )

## posterior analysis
mSummary <- fSummary( mPostParam, sName, "Normal-MH" )
print( paste( "acceptance rate for dMu     =", sprintf( "%6.3f", vAccept[1] / dCount ), sep = "" ) )
print( paste( "acceptance rate for dSigma2 =", sprintf( "%6.3f", vAccept[2] / dCount ), sep = "" ) )

print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

x.min <- 0.1 * floor( 10 * min( mPostPara[, c( 1, 3 )] ) )
x.max <- 0.1 * ceiling( 10 * max( mPostPara[, c( 1, 3 )] ) )
y.max <- 0.1 * ceiling( 10 * ( max( density( mPostPara[, 1] )$y, density( mPostPara[, 3] )$y ) ) )
sName <- c( "ギブズ・サンプラー", "ギブズ・サンプラー", "MHアルゴリズム", "MHアルゴリズム" )
pdf( file = "normal-comparison-mu.pdf", width = 7, height = 5 )
par( mfrow = c( 2, 3 ), family = "Japan1Ryumin" )
for ( i in c( 1, 3 ) )
{
  par( mar = c( 2, 2, 1.5, 0.5 ) )
  ts.plot( mPostPara[, i], main = sName[i], xlab = "", ylab = "", ylim = c( x.min, x.max ), gpars = list( bty = "n" ) )
  plot( density( mPostPara[, i] ), main = sName[i], xlab = "", ylab = "", xlim = c( x.min, x.max ), ylim = c( 0, y.max ), bty = "n" )
  par( mar = c( 2, 2, 3, 0.5 ) )
  acf( mPostPara[, i], main = sName[i], xlab = "", ylab = "", bty = "n", lag = 50 )
}
dev.off(  )
