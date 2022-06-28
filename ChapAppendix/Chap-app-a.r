rm( list = ls(  ) )
set.seed( 123456 )

fLogL <- function( vParam, vY, mX )
{
	nn <- length( vParam )
	vBeta   <- vParam[1 : ( nn - 1 )]
	dSigma2 <- vParam[nn]

	cN    <- length( vY )
	vE    <- vY - mX %*% vBeta
	dLogL <- -0.5 * cN * log( dSigma2 ) - 0.5 * sum( vE ^ 2 ) / dSigma2

	return( dLogL )
}

## load data
sData <- read.csv( "../Data/mi.csv" )
vZ.0 <- sData$Y
vZ.1 <- sData$K
vZ.2 <- sData$L
vY    <- log( vZ.0 ) - log( vZ.2 )
vX.1  <- log( vZ.1 ) - log( vZ.2 )
vX.2  <- log( vZ.2 )
mX    <- cbind( 1.0, vX.1, vX.2 )

cK <- ncol( mX )

vParam <- c( 2, 0.6, 0.1, 0.02 )
vL     <- c( rep( -Inf, cK ), 0.0 )
res    <- optim( vParam, fLogL, method = "L-BFGS-B", lower = vL, hessian = TRUE, control = list( fnscale=-1 ), vY = vY, mX = mX )
print( res$par )
