rm( list = ls(  ) )
source( "../MCMC.r" )
set.seed( 123456 )

library( boot )

## set timer
ttt <- proc.time(  )

## load data
data( coal )
vY <- tabulate( floor( coal[[1]] ) )
vY <- vY[1851 : length( vY )]

## set dimensions
cT <- length( vY )

## set hyper-parameters
dA.0 <- 0.01
dB.0 <- 0.01
dC.0 <- 0.01
dD.0 <- 0.01

## set initial parameters
dLamb <- 1.0
dMu   <- 1.0
dK    <- as.integer( 0.5 * cT )

## set burn-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set posterior matrices
mPostParam <- matrix( 0, iDraw, 3 )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", dK = ", dK,  sep = "" ) )
	}

	## sampling dLamb
	dA.1 <- sum( vY[1 : dK] ) + dA.0
	dB.1 <- dK + dB.0

	dLamb <- rgamma( 1, dA.1, dB.1 )

	## sampling dMu
	dC.1 <- dC.0
	if ( dK < cT )
	{
		dC.1 <- dC.1 + sum( vY[( dK + 1 ) : cT] )
	}
	dD.1 <- cT - dK + dD.0

	dMu <- rgamma( 1, dC.1, dD.1 )

	## sampling dK
	vProb <- rep( 0, cT )
	for ( tt in 1 : cT )
	{
		dLP <- log( dLamb ) * sum( vY[1 : tt] ) - dLamb * tt
		if ( tt < cT )
		{
			 dLP <- dLP + log( dMu ) * sum( vY[( tt + 1 ) : cT] ) - dMu * ( cT - tt )
		}
		vProb[tt] <- exp( dLP )
	}
	vProb <- vProb / sum( vProb )
	
	dK <- sample( 1 : cT, 1, 0, vProb )

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, ] <- c( dLamb, dMu, 1850 + dK )
	}
}

## set parameter names
sName <- c( bquote( lambda ), bquote( mu ), bquote( k ) )

## posterior analysis
fDraw( mPostParam, sName, "Poisson" )
mSummary <- fSummary( mPostParam, sName, "Poisson" )

## draw figures
pdf( file = "poi-tsplot.pdf", width = 7, height = 4 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Palatino" )
ts.plot( vY, type = "b", main = "", xlab = "", ylab = "", gpars = list( cex = 0.7, xaxt = "n", bty = "n" ) )
axis( 1, at = seq( 10, cT, by = 20 ), labels = seq( 1860, 1850 + cT, by = 20 ) )
dev.off(  )

sX.1  <- density( mPostParam[, 1] )
x.1.L <- quantile( mPostParam[, 1], probs = 0.025 )
x.1.U <- quantile( mPostParam[, 1], probs = 0.975 )
sX.2  <- density( mPostParam[, 2] )
x.2.L <- quantile( mPostParam[, 2], probs = 0.025 )
x.2.U <- quantile( mPostParam[, 2], probs = 0.975 )
x.min <- 0.01 * floor( 100 * min( sX.1$x, sX.2$x ) )
x.max <- 0.01 * ceiling( 100 * max( sX.1$x, sX.2$x ) )
y.max <- 0.01 * ceiling( 100 * max( sX.1$y, sX.2$y ) )

pdf( file = "poi-lm.pdf", width = 7, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Palatino" )
plot( sX.1, main = "", xlab = "", ylab = "", lwd = 1, lty = 1, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), bty = "n" )
polygon( c( x.1.L, sX.1$x[sX.1$x > x.1.L & sX.1$x < x.1.U], x.1.U ), c( 0, sX.1$y[sX.1$x > x.1.L & sX.1$x < x.1.U], 0 ), col = grey( 0.8 ), lty = 0 )
polygon( c( x.2.L, sX.2$x[sX.2$x > x.2.L & sX.2$x < x.2.U], x.2.U ), c( 0, sX.2$y[sX.2$x > x.2.L & sX.2$x < x.2.U], 0 ), col = grey( 0.8 ), lty = 0 )
par( new = T )
plot( sX.2, main = "", xlab = "", ylab = "", lwd = 1, lty = 2, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), xaxt = "n", yaxt = "n", bty = "n" )
legend( "topright",  legend = c( expression( ),  sName[1 : 2] ), lty = c( 1, 2 ), lwd = 1, bty = "n" )
dev.off(  )

pdf( file = "poi-k.pdf", width = 7, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Palatino" )
barplot( table( factor( mPostParam[, 3], levels = min( mPostParam[, 3] ) : max( mPostParam[, 3] ) ) ), bty = "n" )
dev.off(  )

print( paste( "Execution time: ", sprintf( "%7.3f", ( proc.time(  ) - ttt )[3] ), sep = "" ) )
