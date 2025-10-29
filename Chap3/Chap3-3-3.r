rm( list = ls(  ) )
source( "../MCMC.r" )
set.seed( 123456 )

SimData <- function( cN, dNu )
{
	## set true parameters
	dNu     <- 5.0

	## generate vY
	vY <- rt( cN, dNu )

	return( list( Y = vY ) )
}

fLogIG <- function( dX, dNu, dLamb )
{
	return( - ( dNu + 1.0 ) * log( dX ) - dLamb / dX )
}

fLogL <- function( vY, dNu, dSigma2 )
{
	cN <- length( vY )

	dLogL <- - 0.5 * cN * log( dSigma2 ) - 0.5 * ( dNu + 1 ) * sum( log( 1.0 +  vY ^ 2 / ( dNu * dSigma2 ) ) )

	return( dLogL )
}

## simulate data
sData <- SimData( 300, 5.0 )
vY    <- sData$Y

## Data Augmentation
## set timer
ttt <- proc.time(  )

## set dimensions
cN <- length( vY )

## set hyper-parameters
dNu.0   <- 2.0
dLamb.0 <- 0.001
dNu     <- 5.0

## set initial parameters
dSigma2 <- 0.1
vW      <- rep( 1.0, cN )

## set burn-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dSumY <- sum( vY )
dNu.1 <- cN + dNu.0
dNu.2 <- dNu + 1.0

## set posterior matrices
vPostParam <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dSigma2 = ", sprintf( "%7.3f", dSigma2 ),  sep = "" ) )
	}

	## sampling dSigma2
	dLamb.1 <- sum( vY ^ 2 / vW ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## sampling vW
	vLamb.2 <- vY ^ 2 / dSigma2 + dNu
	vW <- 1.0 / rgamma( cN, 0.5 * dNu.2, 0.5 * vLamb.2 )

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		vPostParam[it] <- dSigma2
	}
}
mPostParam <- vPostParam

print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## MH-algorithm
## set timer
ttt <- proc.time(  )

## set dimensions
cN <- length( vY )

## set hyper-parameters
dNu.0   <- 2.0
dLamb.0 <- 0.001
dNu     <- 5.0

## set initial parameters
dSigma2 <- 0.1

## set burn-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dLogL.o <- fLogIG( dSigma2, 0.5 * dNu.0, 0.5 * dLamb.0 ) + fLogL( vY, dNu, dSigma2 )
dCount  <- 0
dAccept <- 0
dC      <- 0.1

## set posterior matrices
vPostParam <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dSigma2 = ", sprintf( "%7.3f", dSigma2 ),  sep = "" ) )
	}
	dCount <- dCount + 1

	## sampling dSigma2
	dSigma2.n <- dSigma2 + dC * rnorm( 1 )
	if ( dSigma2.n > 0.0 )
	{
		dLogL.n <- fLogIG( dSigma2.n, 0.5 * dNu.0, 0.5 * dLamb.0 ) + fLogL( vY, dNu, dSigma2.n )
		if ( runif( 1 ) <= exp( dLogL.n - dLogL.o ) )
		{
			dSigma2 <- dSigma2.n
			dAccept <- dAccept + 1
			dLogL.o <- dLogL.n
		}
	}
	if ( iter <= iBurn )
	{
		if( dAccept / dCount < 0.4 )
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

		vPostParam[it] <- dSigma2
	}
}
mPostParam <- cbind( mPostParam, vPostParam )

print( paste( "acceptance rate for dSigma2 =", sprintf( "%6.3f", dAccept / dCount ), sep = "" ) )
print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## Gibbs Normal-model
## set timer
ttt <- proc.time(  )

## set dimensions
cN <- length( vY )

## set hyper-parameters
dNu.0   <- 2.0
dLamb.0 <- 0.001

## set initial parameters
dSigma2 <- 0.1

## set burn-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dNu.1   <- cN + dNu.0
dLamb.1 <- sum( vY ^ 2 ) + dLamb.0

## set posterior matrices
vPostParam <- rep( 0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dSigma2 = ", sprintf( "%7.3f", dSigma2 ),  sep = "" ) )
	}

	## sampling dSigma2
	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		vPostParam[it] <- dSigma2
	}
}
mPostParam <- cbind( mPostParam, vPostParam )

print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

sName <- c( "データ拡大法", "MHアルゴリズム", "ギブズ・サンプラー (正規分布)" )

mSummary  <- fSummary( mPostParam, sName, "t-DA" )

x.min <- 0.1 * floor( 10 * min( mPostParam[, 1 : 2] ) )
x.max <- 0.1 * ceiling( 10 * max( mPostParam[, 1 : 2] ) )
y.max <- 0.1 * ceiling( 10 * ( max( density( mPostParam[, 1] )$y, density( mPostParam[, 2] )$y ) ) )
pdf( file = "t-comparison.pdf", width = 7, height = 5 )
par( mfrow = c( 2, 3 ), family = "Japan1Ryumin" )

for ( i in 1 : 2 )
{
	par( mar = c( 2, 2, 1.5, 0.5 ) )
	ts.plot( mPostParam[, i], main = sName[i], xlab = "", ylab = "", ylim = c( x.min, x.max ), gpars = list( bty = "n" ) )
	plot( density( mPostParam[, i] ), main = sName[i], xlab = "", ylab = "", xlim = c( x.min, x.max ), ylim = c( 0, y.max ), bty = "n" )
	par( mar = c( 2, 2, 3, 0.5 ) )
	acf( mPostParam[, i], main = sName[i], xlab = "", ylab = "", bty = "n", lag = 50 )
}
dev.off(  )
