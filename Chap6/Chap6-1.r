rm( list = ls(  ) )
set.seed( 123456 )
source( "../MCMC.r" )

## set timer
ttt <- proc.time(  )

sData <- read.csv( "../Data/bc.csv" )
vY <- sData$BC
vX <- sData$CI
mX <- cbind( 1.0, sData$C1, sData$C2, sData$C3, sData$C4, sData$C5, sData$C6, sData$C7, sData$C8, sData$C9, sData$C10 )

## set dimensions
cT <- nrow( mX )
cK <- ncol( mX )

## set hyper-parameters
vBeta.0  <- rep( 0.0, cK )
mSigma.0 <- diag( 100.0, cK )

## set some matrices, vectors and scalars
mInvSigma.0 <- solve( mSigma.0 )
vISB.0      <- mInvSigma.0 %*% vBeta.0

t.mX <- t( mX )
mXX  <- t.mX %*% mX

## set initial parameters
vZ    <- vY
vBeta <- solve( mXX ) %*% t( mX ) %*% vZ

## set burn-in & draws
iBurn <- 50000
iDraw <- 50000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set posterior matrices
mPostParam <- matrix( 0, iDraw, cK )
mPostZ     <- matrix( 0, iDraw, cT )
mPostY     <- matrix( 0, iDraw, cT )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	# monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", vBeta[1] = ", sprintf( "%7.3f", vBeta[1] ),  sep = "" ) )
	}

	## draw vBeta
	mSigma.1 <- solve( mXX + mInvSigma.0 )
	vBeta.1  <- mSigma.1 %*% ( t.mX %*% vZ + vISB.0 )

	vBeta <- vBeta.1 + t( chol( mSigma.1 ) ) %*% rnorm( cK )

	## draw vZ
	vMu <- mX %*% vBeta
	for ( tt in 1 : cT )
	{
		if ( vY[tt] == 1 )
		{
			vZ[tt] <- rtnorm( vMu[tt], 1, 0, Inf )
		}
		else
		{
			vZ[tt] <- rtnorm( vMu[tt], 1, -Inf, 0 )
		}
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		mPostParam[it, ] <- vBeta
		mPostZ[it, ]     <- vZ
		mPostY[it, ]     <- ifelse( mX %*% vBeta + rnorm( cT ) >= 0.0, 1, 0 )
	}
}

## set parameter names
sName <- NULL
for ( i in 1 : cK )
{
	sName <- c( sName, bquote( beta[.( i )] ) )
}
sName.Z <- NULL
for ( i in 1 : cT )
{
	sName.Z <- c( sName.Z, bquote( z[.( i )] ) )
}

## posterior analysis
mSummary   <- fSummary( mPostParam, sName, "Probit" )
mSummary.Z <- fSummary( mPostZ, sName.Z, "Z" )
mSummary.Y <- fSummary( mPostY, sName.Z, "Y" )
fDraw( mPostParam, sName, "Probit" )
#fDraw( mPostZ, sName.Z, "ProbitZ" )

print( paste( "Execution time: ", ( proc.time(  ) - ttt )[3], sep = "" ) )

## draw figures
vP <- NULL
for ( tt in 2 : cT )
{
  if ( vY[tt] - vY[( tt - 1 )] == 1 )
  {
    vP <- c( vP, tt )
  }
  if ( vY[tt] - vY[( tt - 1 )] == -1 )
  {
    vP <- c( vP, tt - 1 )
  }
}
vP <- c( vP, cT )
cK <- length( vP ) / 2
y.min <- 10 * floor( 0.1 * min( vX ) )
y.max <- 10 * ceiling( 0.1 * max( vX ) )
inter <- seq( 1, cT, by = 12)
pdf( file = "bc.pdf", width = 10, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
plot.new(  )
plot.window( xlim = c( 1, cT ), ylim = c( y.min, y.max ) )
for ( k in 1 : cK )
{
  dS <- vP[( 2 * k - 1 )]
  dE <- vP[( 2 * k )]
  polygon( c( dS, dS, dE, dE), c( y.min, y.max, y.max, y.min ), col="gray", border = NA )
}
par( new = T )
plot( vX, type = "l", ylim = c( y.min, y.max ), main = "", xlab = "", ylab = "", bty = "n", xaxt = "n" )
axis( 1, at = inter, labels = sData$年月[inter] )
dev.off(  )

pdf( file = "prediction.pdf", width = 10, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
plot.new(  )
plot.window( xlim = c( 1, cT ), ylim = c( 0, 1 ) )
for ( k in 1 : cK )
{
  dS <- vP[( 2 * k - 1 )]
  dE <- vP[( 2 * k )]
  polygon( c( dS, dS, dE, dE), c( 0, 1, 1, 0 ), col="gray", border = NA )
}
par( new = T )
plot( mSummary.Y[, 1], type = "l", ylim = c( 0, 1 ), main = "", xlab = "", ylab = "", bty = "n", xaxt = "n" )
axis( 1, at = inter, labels = sData$年月[inter] )
dev.off(  )
