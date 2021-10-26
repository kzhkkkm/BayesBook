rm( list = ls(  ) )
source( "../MCMC.r" )
set.seed( 123456 )

fKernel <- function( vX )
{
	cN <- length( vX )

	dPhi    <- 0.9
	dSigma2 <- 1.0 - 0.9 ^ 2

	return( dPhi * vX + sqrt( dSigma2 ) * rnorm( cN ) )
}

## number of draws
iDraw <- 50000

## multiple chain
cK <- 50
vX <- runif( iDraw, -3, 3 )
mParam.M <- matrix( 0, iDraw, 4 )
mParam.M[,1] <- vX
for ( iter in 1 : cK )
{
	vX <- fKernel( vX )
	if ( iter == 1 )
	{
		mParam.M[, 2] <- vX
	}
	if ( iter == 10 )
	{
		mParam.M[, 3] <- vX
	}
	if ( iter == cK )
	{
		mParam.M[, 4] <- vX
	}
}

## single chain
dX <- runif( 1, -3, 3 )
vParam.S <- rep( 0, iDraw )
for ( iter in 1 : iDraw)
{
	dX <- fKernel( dX )
	vParam.S[iter] <- dX
}

## thinning
k <- seq( 1, iDraw, by = 10 )
vParam.T <- vParam.S[k]

## draw figure
sDensity.M1 <- density( mParam.M[, 2] )
sDensity.M2 <- density( mParam.M[, 3] )
sDensity.M3 <- density( mParam.M[, 4] )
sDensity.S <- density( vParam.S )
x.max <- ceiling( max( mParam.M, vParam.S ) )
x.min <- floor( min( mParam.M, vParam.S ) )
y.max <- 0.01 * ceiling( 100 * max( sDensity.M1$y, sDensity.M2$y, sDensity.M3$y, sDensity.S$y ) )
vX <- seq( x.min, x.max, by = 0.01 )
vY <- dnorm( vX )
pdf( file = "mc-sc-dist.pdf", height = 5, width = 10 )
par( mfrow = c( 1, 2 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
hist( mParam.M[, 1], freq = F, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "", xlab = "", ylab = "", bty = "n" )
par( new = T )
plot( sDensity.M1, lty = 3, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n" )
par( new = T )
plot( sDensity.M2, lty = 2, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n" )
par( new = T )
plot( sDensity.M3, lty = 1, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "多重連鎖", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n" )
par( new = T )
plot( vX, vY, type = "l", lty = 4, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n" )
legend( "topright", legend = c( "1回", "10回", "50回", "不変分布" ), lty = c( 3, 2, 1, 4 ), bty = "n" )
plot( sDensity.S, lty = 1, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "単一連鎖", xlab = "", ylab = "", bty = "n" )
par( new = T )
plot( vX, vY, type = "l", lty = 4, xlim = c( x.min, x.max ), ylim = c( 0, y.max ), main = "", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n" )
legend( "topright", legend = c( "単一連鎖", "不変分布" ), lty = c( 1, 4 ), bty = "n" )
dev.off(  )

pdf( file = "mc-sc-acf.pdf" )
par( mfrow = c( 2, 2 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
plot( mParam.M[1 : ( iDraw - 1 ), 4], mParam.M[2 : iDraw, 4], pch = 20, cex = 0.2, main = "多重連鎖", xlab = "", ylab = "", bty = "n" )
plot( vParam.S[1 : ( iDraw - 1 )], vParam.S[2 : iDraw], pch = 20, cex = 0.2, main = "単一連鎖", xlab = "", ylab = "", bty = "n" )
acf( mParam.M[4,], main = "", xlab = "", ylab = "", xlim = c(0, 50), bty = "n", lag = 50 )
acf( vParam.S, main = "", xlab = "", ylab = "", bty = "n", lag = 50 )
dev.off(  )

y.max <- ceiling( max( vParam.S ) )
y.min <- floor( min( vParam.S ) )
pdf( file = "thinning.pdf" )
par( mfrow = c( 2, 2 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
ts.plot( vParam.S, ylim = c( y.min, y.max ), main = "全てのMCMC標本", xlab = "", ylab = "", gpars = list( bty = "n" ) )
ts.plot( vParam.T, ylim = c( y.min, y.max ), main = "シニングしたMCMC標本", xlab = "", ylab = "", gpars = list( bty = "n" ) )
acf( vParam.S, ylim = c( -0.02, 1.0 ), main = "", xlab = "", ylab = "", bty = "n", lag = 50 )
acf( vParam.T, ylim = c( -0.02, 1.0 ), main = "", xlab = "", ylab = "", bty = "n", lag = 50 )
dev.off(  )

fSummary( cbind( mParam.M[, 2 : 4], vParam.S ), c( "多重連鎖 ($t=1$)", "多重連鎖 ($t=10$)", "多重連鎖 ($t=50$)", "単一連鎖" ), "full" )
fSummary( vParam.T, "シニング", "thinn" )
