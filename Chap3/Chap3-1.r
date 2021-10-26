rm( list = ls(  ) )
source( "../MCMC.r" )
set.seed( 123456 )

dtnorm <- function( x, dMu, dSigma2, dA, dB )
{
	dSigma <- sqrt( dSigma2 )
	if ( dA == -Inf )
	{
		p.a <- 0.0
	}
	else
	{
		p.a <- pnorm( ( dA - dMu ) / dSigma )
	}
	if ( dB == Inf )
	{
		p.b <- 1.0
	}
	else
	{
		p.b <- pnorm( ( dB - dMu ) / dSigma )
	}
	z <- dnorm( ( x - dMu ) / dSigma ) / ( dSigma * ( p.b - p.a ) )

	return( ifelse( x >= dA & x <= dB, z, 0.0 ) )
}

## number of draws
iDraw <- 50000

## Monte Carlo samples from truncated normal distribution
mu     <- 1.0
sigma2 <- 1.0
aa     <- -1.0
bb     <- 2.0
vParam <- rep( 0, iDraw )
for ( iter in 1 : iDraw )
{
	vParam[iter] <- rtnorm( mu, sigma2, aa, bb )
}
x <- seq( -1.0, 2.0, by = 0.01 )
y <- dtnorm( x, 1.0, 1.0, -1.0, 2.0 )
pdf( file = "mc-tnorm.pdf" )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 0.5, 0.5 ), family = "Japan1Ryumin" )
hist( vParam, freq = F, ylim = c( 0.0, 0.5 ), main = "", xlab = "", ylab = "", bty = "n" )
par( new = T )
plot( x, y, type = "l", ylim = c( 0.0, 0.5 ), main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n" )
dev.off(  )
tn.mean <- mu + sqrt( sigma2 ) * ( dnorm( ( aa - mu ) / sqrt( sigma2 ) ) - dnorm( ( bb - mu ) / sqrt( sigma2 ) ) ) / ( pnorm( ( bb - mu ) / sqrt( sigma2 ) ) - pnorm( ( aa - mu ) / sqrt( sigma2 ) ) )
tn.var  <- ( 1.0 + ( ( ( aa - mu ) / sqrt( sigma2 ) ) * dnorm( ( aa - mu ) / sqrt( sigma2 ) ) - ( ( bb - mu ) / sqrt( sigma2 ) ) * dnorm( ( bb - mu ) / sqrt( sigma2 ) ) ) / ( pnorm( ( bb - mu ) / sqrt( sigma2 ) ) - pnorm( ( aa - mu ) / sqrt( sigma2 ) ) ) - ( ( dnorm( ( aa - mu ) / sqrt( sigma2 ) ) - dnorm( ( bb - mu ) / sqrt( sigma2 ) ) ) / ( pnorm( ( bb - mu ) / sqrt( sigma2 ) ) - pnorm( ( aa - mu ) / sqrt( sigma2 ) ) ) ) ^ 2 )
sink( "tnorm.txt" )
print( round( c( tn.mean, mean( vParam ), tn.var, var( vParam ) ), digits = 3 ) )
sink(  )
