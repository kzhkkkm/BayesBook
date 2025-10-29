rm( list = ls(  ) )
source( "../MCMC.r" )
set.seed( 123456 )

fLogB <- function( dX, dAlpha, dBeta )
{
	dY <- 0.5 * ( dX + 1.0 )

	return( ( dAlpha - 1.0 ) * log( dY ) + ( dBeta - 1.0 ) * log( 1.0 - dY ) )
}

fLogL <- function( vH, dMu, dPhi, dSigma2 )
{
	cT <- length( vH )

	vHH <- vH[2 : cT]
	vHL <- vH[1 : ( cT - 1 )]
	
	dLogL <- - 0.5 * ( 1.0 - dPhi ^ 2 ) * ( vH[1] - dMu ) ^ 2 / dSigma2
	dLogL <- dLogL - 0.5 * sum( ( vHH - dMu - dPhi * ( vHL - dMu ) ) ^ 2 ) / dSigma2

	return( dLogL )
}

fLogH <- function( tt, vY, vH, dMu, dPhi, dSigma2 )
{
	cT <- length( vY )

	dLogH <- - 0.5 * vH[tt] - 0.5 * vY[tt] ^ 2 / exp( vH[tt] )
	if ( tt == 1 )
	{
		dLogH <- dLogH - 0.5 * ( 1.0 - dPhi ^ 2 ) * ( vH[tt] - dMu ) ^ 2 / dSigma2
	}
	if ( tt < cT )
	{
		dLogH <- dLogH - 0.5 * ( vH[(tt + 1)] - dMu - dPhi * ( vH[tt] - dMu ) ) ^ 2 / dSigma2
	}
	if ( tt > 1 )
	{
		dLogH <- dLogH - 0.5 * ( vH[tt] - dMu - dPhi * ( vH[( tt - 1 )] - dMu ) ) ^ 2 / dSigma2
	}

	return( dLogH )
}

## set timer
ttt <- proc.time(  )

## load data
sData <- read.csv( "../Data/ExRate.csv", header = T )
vY    <- 100 * diff( log( sData$ドル ) )

## set dimensions
cT <- length( vY )

## set hyper-parameters
dMu.0    <- 0.0
dTau2.0  <- 100.0
dNu.0    <- 2.0
dLamb.0  <- 0.01
dAlpha.0 <- 25.0
dBeta.0  <- 1.5

## set initial parameters
dMu     <- 0.0
dPhi    <- 0.0
dSigma2 <- 0.1
vH      <- rep( 0.0, cT )

## set burn-in & draws
iBurn <- 10000
iDraw <- 1000000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )
k     <- seq( 1, iDraw, by = 100 )

## set some matrices, vectors and scalars
dNu.1   <- cT + dNu.0
dCount  <- 0
dAccept <- 0
dC      <- 0.1
vAccept <- rep( 0, cT )
vC      <- rep( 0.1, cT )
it      <- 1

## set posterior matrices
mPostParam.sv <- matrix( 0, length( k ), 3 )
mPostH.sv     <- matrix( 0, length( k ), cT )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dMu = ", sprintf( "%7.3f", dMu ), ", dC = ", sprintf( "%7.3f", dC ), sep = "" ) )
	}
	dCount <- dCount + 1

	## draw dMu
	dTau2.1 <- 1.0 / ( ( 1.0 - dPhi ^ 2 + ( cT - 1 ) * ( 1.0 - dPhi ) ^ 2 ) / dSigma2 + 1.0 / dTau2.0 )
	dMu.1   <- dTau2.1 * ( ( ( 1.0 - dPhi ^ 2 ) * vH[1] + ( 1.0 - dPhi ) * sum( vH[2 : cT] - dPhi * vH[1 : ( cT - 1 )] ) ) / dSigma2 + dMu.0 / dTau2.0 )

	dMu <- dMu.1 + sqrt( dTau2.1 ) * rnorm( 1 )

	## draw dPhi
	dPhi.n <- dPhi + dC * rnorm( 1 )
	if ( abs( dPhi.n ) < 1.0 )
	{
		dLogL.o <- fLogB( dPhi, dAlpha.0, dBeta.0 ) + fLogL( vH, dMu, dPhi, dSigma2 )
		dLogL.n <- fLogB( dPhi.n, dAlpha.0, dBeta.0 ) + fLogL( vH, dMu, dPhi.n, dSigma2 )
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

	## draw dSigma2
	dLamb.1 <- ( 1.0 - dPhi ^ 2 ) * ( vH[1] - dMu ) ^ 2 + sum( ( vH[2 : cT] - dMu - dPhi * ( vH[1 : ( cT - 1 )] - dMu ) ) ^ 2 ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## draw vH
	for ( tt in 1 : cT )
	{
		vH.n     <- vH
		vH.n[tt] <- vH[tt] + vC[tt] * rnorm( 1 )
		dLogH.o  <- fLogH( tt, vY, vH, dMu, dPhi, dSigma2 )
		dLogH.n  <- fLogH( tt, vY, vH.n, dMu, dPhi, dSigma2 )
		if ( runif( 1 ) <= exp( dLogH.n - dLogH.o ) )
		{
			vH          <- vH.n
			vAccept[tt] <- vAccept[tt] + 1
		}
		if ( iter <= iBurn )
		{
			if ( vAccept[tt] / dCount < 0.2 )
			{
				vC[tt] <- vC[tt] / 1.1
			}
			if ( vAccept[tt] / dCount > 0.4 )
			{
				vC[tt] <- 1.1 * vC[tt]
			}
		}
	}
	if ( iter == iBurn )
	{
		dCount  <- 0
		dAccept <- 0
		vAccept <- rep( 0, cT )
	}

	## save parameters
	if ( iter == k[it] + iBurn )
	{
		mPostParam.sv[it, ] <- c( dMu, dPhi, dSigma2 )
		mPostH.sv[it, ]     <- vH
		if ( it == length( k ) )
		{
			it <- 1
		}
		else
		{
			it <- it + 1
		}
	}
}

## set parameter names
sName.sv <- c( bquote( mu ), bquote( phi ), bquote( sigma^2 ) )
sName.H  <- NULL
for ( tt in 1 : cT )
{
	sName.H <- c( sName.H, bquote( h[.( tt )] ) )
}

print( paste( "acceptance rate =", sprintf( "%6.3f", dAccept / dCount ), sep = "" ) )
print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## set timer
ttt <- proc.time(  )

## set dimensions
cT <- length( vY )

## set hyper-parameters
dMu.0   <- 0.0
dTau2.0 <- 100.0
dNu.0   <- 2.0
dLamb.0 <- 0.01

## set initial parameters
dMu     <- mean( vY )
dSigma2 <- var( vY )

## set burn-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set some matrices, vectors and scalars
dSumY <- sum( vY )
dNu.1 <- cT + dNu.0

## set posterior matrices
mPostParam.hv <- matrix( 0, iDraw, 3 )

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
	dTau2.1 <- 1.0 / ( cT / dSigma2 + 1.0 / dTau2.0 )
	dMu.1   <- dTau2.1 * ( dSumY / dSigma2 + dMu.0 / dTau2.0 )

	dMu <- dMu.1 + sqrt( dTau2.1 ) * rnorm( 1 )

	## sampling dSigma2
	dLamb.1 <- sum( ( vY - dMu ) ^ 2 ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam.hv[it, ] <- c( dMu, dSigma2, sqrt( dSigma2 ) )
	}
}

## set parameter names
sName.hv   <- c( bquote( mu ), bquote( sigma^2 ), bquote( sigma ) )

print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )

## posterior analysis
mSummary.sv   <- fSummary( mPostParam.sv, sName.sv, "SV" )
mSummary.H.sv <- fSummary( mPostH.sv, sName.H, "SV-H" )
mSummary.hv   <- fSummary( mPostParam.hv, sName.hv, "HV" )
fDraw( mPostParam.sv, sName.sv, "SV" )
# fDraw( mPostH.sv, sName.H, "SV-H-" )
fDraw( mPostParam.hv, sName.hv, "HV" )

## draw figures
mSigma.sv <- exp( 0.5 * mSummary.H.sv[, c( 1, 3, 4 )] )
y.min.sv  <- 0.1 * floor( 10 * min( mSigma.sv[, 2] ) )
y.max.sv  <- 0.1 * ceiling( 10 * max( mSigma.sv[, 3] ) )
vSigma.hv <- mSummary.hv[3, c( 1, 3, 4 )]
inter  <- seq( 1, cT, by = 200 )
pdf( file = "sv-dollar-tsplot.pdf", width = 7, height = 5 )
par( mfrow = c( 2, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
ts.plot( 100 * diff( log( sData$ドル ) ), main = "対数変化率",  xlab = "", ylab = "", gpars = list( xaxt = "n", bty = "n" ) )
axis( 1, at = inter, labels = sData$日付[( 1 + inter )] )

ts.plot( mSigma.sv[, 1], main = "ボラティリティ", xlab = "", ylab = "", ylim = c( y.min.sv, y.max.sv ), gpars = list( xaxt = "n", bty = "n" ) )
polygon( c( 1, 1, cT, cT ), c( vSigma.hv[2], vSigma.hv[3], vSigma.hv[3], vSigma.hv[2] ), col = grey( 0.8 ), lty = 0 )
lines( c( 1, cT ), c( vSigma.hv[1], vSigma.hv[1] ), col = "red" )
par( new = T )
ts.plot( mSigma.sv[, 1], main = "", xlab = "", ylab = "", ylim = c( y.min.sv, y.max.sv ), gpars = list( xaxt = "n", bty = "n" ) )
par( new = T )
ts.plot( mSigma.sv[, 2], lty = 3, main = "", xlab = "", ylab = "", ylim = c( y.min.sv, y.max.sv ), gpars = list( xaxt = "n", yaxt = "n", bty = "n" ) )
par( new = T )
ts.plot( mSigma.sv[, 3], lty = 3, main = "", xlab = "", ylab = "", ylim = c( y.min.sv, y.max.sv ), gpars = list( xaxt = "n", yaxt = "n", bty = "n" ) )
axis( 1, at = inter, labels = sData$日付[( 1 + inter )] )
legend( "topright", legend = c( "SVの事後平均", "SVの95%信用区間", "HVの事後平均" ), lty = c( 1, 3, 1 ), col = c( 1, 1, "red" ), bty = "n" )
dev.off(  )
