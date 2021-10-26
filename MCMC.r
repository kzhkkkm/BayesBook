##
## Last Change: 10 06 14:38 PM 2021.
##
#install.packages( "coda", dep = T )
library( coda )

rtnorm <- function( dMu, dSigma2, dLeft, dRight )
{
	dSigma <- sqrt( dSigma2 )
	dP.l   <- ifelse( dLeft == -Inf, 0.0, pnorm( ( dLeft  - dMu ) / dSigma ) )
	dP.r   <- ifelse( dRight ==  Inf, 1.0, pnorm( ( dRight - dMu ) / dSigma ) )
	dU     <- runif( 1, dP.l, dP.r )
	dX     <- qnorm( dU, dMu, dSigma )

	return( dX )
}

fHarmonicMean <- function( vX )
{
	cN <- length( vX )
	vZ <- 1.0 / vX
	
	dMu <- 1.0 / mean( vZ )
	dS  <- dMu ^ 2 * sd( vZ ) / sqrt( cN )

	mY <- cbind( rbind( dMu, dS ), rbind( log( dMu ), dS / dMu ) )
	rownames( mY ) <- c( "mean", "standard error" )
	colnames( mY ) <- c( "standard", "log" )

	return( mY )
}

fSummary <- function( mX, sName = NULL, sFileName = "Output", vProb = c( 0.025, 0.975 ), vLik = NULL )
{
	cK <- NCOL( mX )	

	mSummary <- NULL
	for ( i in 1 : cK )
	{
		if ( cK == 1 )
		{
			vX <- mX
		}
		else
		{
			vX <- mX[, i]
		}
		aa <- round( mean( vX ), digits = 5 )
		bb <- round( sd( vX ), digits = 5 )
		cc <- round( quantile( vX, probs = vProb ), digits = 5 )
		dd <- round( 2.0 * ( 1.0 - pnorm( abs( geweke.diag( vX, frac1 = 0.1, frac2 = 0.5 )$z ) ) ), digits = 5 )
		ee <- round( max( 1.0, spectrum0.ar( vX )$spec / var( vX ) ), digits = 5 )

		mSummary <- rbind( mSummary, c( aa, bb, cc, dd, ee ) )
	}

	sLabel <- c( "Post Mean", "Std Dev", paste( 100 * vProb, "%", sep = ""  ), "Geweke CD", "IF" )
	if ( is.null( sName ) )
	{
		for ( i in 1 : cK )
		{
			sName <- c( sName, paste( "param", i, sep = "" ) )
		}
	}
	colnames( mSummary ) <- sLabel
	rownames( mSummary ) <- sName

	sink( paste( sFileName, ".txt", sep = "" ) )
	print( mSummary )
	print( paste( "# of no-convergence:", sprintf( "%3d", sum( mSummary[!is.na( mSummary[, 5] ), 5] < 0.01 ) ), sep = "" ) )
	if ( !is.null( vLik ) )
	{
		hm <- fHarmonicMean( vLik )
		print( "marginal likelihood:" )
		print( hm )
	}
	sink(  )

	return( mSummary )
}

fDraw <- function( mX, sName = NULL, sFileName = "Fig" )
{
	cK <- NCOL( mX )
	j  <- 1
	aa <- cK %/% 3

	if ( is.null( sName ) )
	{
		for ( i in 1 : cK )
		{
			sName <- c( sName, paste( "param", i, sep = "" ) )
		}
	}
	
	for ( i in 1 : cK )
	{
		if ( i %% 3 == 1 )
		{
			if ( cK <= 3 )
			{
				pdf( file = paste( sFileName, ".pdf", sep = "" ) )
			}
			else
			{
				pdf( file = paste( sFileName, j, ".pdf", sep = "" ) )
			}
			bb <- ifelse( i %/% 3 < aa, 3, cK %% 3 )
			par( mfrow = c( bb, 3 ), family = "Times" )
		}
		if ( cK == 1 )
		{
			vX <- mX
		}
		else
		{
			vX <- mX[, i]
		}
		par( mar = c( 2, 2, 1.5, 0.5 ) )
		ts.plot( vX, main = sName[i], xlab = "", ylab = "", gpars = list( bty = "n" ) )
		plot( density( vX ), main = sName[i], xlab = "", ylab = "", bty = "n" )
		par( mar = c( 2, 2, 2.8, 0.5 ) )
		acf( vX, main = sName[i], xlab = "", ylab = "", bty = "n", lag = 300 )
		if ( i %% 3 == 0 | i == cK )
		{
			dev.off(  )
			j <- j + 1
		}
	}
}
