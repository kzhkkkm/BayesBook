rm( list = ls(  ) )
set.seed( 123456 )
source( "../MCMC.r" )

## set timer
ttt <- proc.time(  )

## load data
sData <- read.csv( "../Data/mi.csv" )
vZ.0 <- sData$Y
vZ.1 <- sData$K
vZ.2 <- sData$L
vY    <- log( vZ.0 ) - log( vZ.2 )
vX.1  <- log( vZ.1 ) - log( vZ.2 )
vX.2  <- log( vZ.2 )
mX    <- cbind( 1.0, vX.1, vX.2 )

## set dimensions
cN <- nrow( mX )
cK <- ncol( mX )

## set hyper-parameters
vBeta.0  <- rep( 0.0, cK )
mSigma.0 <- diag( 100.0, cK )
dNu.0    <- 2.0
dLamb.0  <- 0.01

## set initial values
vBeta   <- rep( 0.0, cK )
dSigma2 <- 1.0

## set biru-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set scalars, vectors & matrices
mInvSigma.0 <- solve( mSigma.0 )
vISB.0      <- mInvSigma.0 %*% vBeta.0

t.mX <- t( mX )
mXX  <- t.mX %*% mX
vXY  <- t.mX %*% vY

dNu.1 <- cN + dNu.0

## set posterior matrices
mPostParam.1 <- matrix( 0, iDraw, cK + 1 )
vPostL <- rep( 0.0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", vBeta[1] = ", sprintf( "%7.3f", vBeta[1] ), sep = "" ) )
	}

	## sampling vBeta
	mSigma.1 <- solve( mXX / dSigma2 + mInvSigma.0 )
	vBeta.1  <- mSigma.1 %*% ( vXY / dSigma2 + vISB.0 )

	vBeta <- vBeta.1 + t( chol( mSigma.1 ) ) %*% rnorm( cK )

	## sampling dSigma2
	vE <- vY - mX %*% vBeta
	dLamb.1 <- sum( vE ^ 2 ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam.1[it, ] <- c( vBeta, dSigma2 )
		vPostL[it]         <- prod( dnorm( ( vY - mX %*% vBeta ) / sqrt( dSigma2 ) ) / sqrt( dSigma2 ) )
	}
}

## set parameter names
sName <- NULL
for ( i in 1 : ( cK + 1 ) )
{
	if ( i <= cK )
	{
		sName <- c( sName, bquote( beta[.( i )] ) )
	}
	else
	{
		sName <- c( sName, bquote( sigma^2 ) )
	}
}

## posterior analysis
fDraw( mPostParam.1, sName, "LRM-m1-" )
mSummary <- fSummary( mPostParam.1, sName, "LRM-m1", vLik = vPostL )
dLogL.1 <- log( 1.0 / mean( 1.0 / vPostL ) )
print( paste( "Execution time: ", ( proc.time(  ) - ttt )[3], sep = "" ) )

## set timer
ttt <- proc.time(  )

## load data
vY    <- log( vZ.0 ) - log( vZ.2 )
vX.1  <- log( vZ.1 ) - log( vZ.2 )
mX    <- cbind( 1.0, vX.1 )

## set dimensions
cN <- nrow( mX )
cK <- ncol( mX )

## set hyper-parameters
vBeta.0  <- rep( 0.0, cK )
mSigma.0 <- diag( 100.0, cK )
dNu.0    <- 2.0
dLamb.0  <- 0.01

## set initial values
vBeta   <- rep( 0.0, cK )
dSigma2 <- 1.0

## set biru-in & draws
iBurn <- 10000
iDraw <- 10000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set scalars, vectors & matrices
mInvSigma.0 <- solve( mSigma.0 )
vISB.0      <- mInvSigma.0 %*% vBeta.0

t.mX <- t( mX )
mXX  <- t.mX %*% mX
vXY  <- t.mX %*% vY

dNu.1 <- cN + dNu.0

## set posterior matrices
mPostParam.2 <- matrix( 0, iDraw, cK + 1 )
vPostL <- rep( 0.0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", vBeta[1] = ", sprintf( "%7.3f", vBeta[1] ), sep = "" ) )
	}

	## sampling vBeta
	mSigma.1 <- solve( mXX / dSigma2 + mInvSigma.0 )
	vBeta.1  <- mSigma.1 %*% ( vXY / dSigma2 + vISB.0 )

	vBeta <- vBeta.1 + t( chol( mSigma.1 ) ) %*% rnorm( cK )

	## sampling dSigma2
	vE <- vY - mX %*% vBeta
	dLamb.1 <- sum( vE ^ 2 ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam.2[it, ] <- c( vBeta, dSigma2 )
		vPostL[it]       <- prod( dnorm( ( vY - mX %*% vBeta ) / sqrt( dSigma2 ) ) / sqrt( dSigma2 ) )
	}
}

## set parameter names
sName <- NULL
for ( i in 1 : ( cK + 1 ) )
{
	if ( i <= cK )
	{
		sName <- c( sName, bquote( beta[.( i )] ) )
	}
	else
	{
		sName <- c( sName, bquote( sigma^2 ) )
	}
}

## posterior analysis
fDraw( mPostParam.2, sName, "LRM-m2" )
mSummary <- fSummary( mPostParam.2, sName, "LRM-m2", vLik = vPostL )
dLogL.2 <- log( 1.0 / mean( 1.0 / vPostL ) )
print( paste( "Execution time: ", ( proc.time(  ) - ttt )[3], sep = "" ) )
print( paste( "Bayes Factor = ", exp( dLogL.1 - dLogL.2 ), sep = "" ) )

## draw figures
dens.1.1 <- density( mPostParam.1[, 1] )
dens.1.2 <- density( mPostParam.1[, 2] )
dens.2.1 <- density( mPostParam.2[, 1] )
dens.2.2 <- density( mPostParam.2[, 2] )
min.1.x <- min( dens.1.1$x, dens.2.1$x )
max.1.x <- max( dens.1.1$x, dens.2.1$x )
min.2.x <- min( dens.1.2$x, dens.2.2$x )
max.2.x <- max( dens.1.2$x, dens.2.2$x )
min.1.y <- min( dens.1.1$y, dens.2.1$y )
max.1.y <- max( dens.1.1$y, dens.2.1$y )
min.2.y <- min( dens.1.2$y, dens.2.2$y )
max.2.y <- max( dens.1.2$y, dens.2.2$y )
pdf( file = "posterior.pdf", width = 10, height = 5 )
par( mfrow = c( 1, 2 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
plot( dens.1.1, main = expression( paste( beta[1], "の事後分布", sep = "" ) ), xlab = "", ylab = "", bty = "n", xlim = c( min.1.x, max.1.x ), ylim = c( min.1.y, max.1.y ), lty = 1 )
par( new = T )
plot( dens.2.1, main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", xlim = c( min.1.x, max.1.x ), ylim = c( min.1.y, max.1.y ), lty = 2 )
legend( "topright", legend = c( "(5.5)式", "(5.7)式" ), lty = c( 1, 2 ), bty = "n" )
plot( dens.1.2, main = expression( paste( beta[2], "の事後分布", sep = "" ) ), xlab = "", ylab = "", bty = "n", xlim = c( min.2.x, max.2.x ), ylim = c( min.2.y, max.2.y ), lty = 1 )
par( new = T )
plot( dens.2.2, main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", xlim = c( min.2.x, max.2.x ), ylim = c( min.2.y, max.2.y ), lty = 2 )
legend( "topright", legend = c( "(5.5)式", "(5.7)式" ), lty = c( 1, 2 ), bty = "n" )
dev.off(  )
