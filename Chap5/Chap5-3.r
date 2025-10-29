rm( list = ls(  ) )
set.seed( 123456 )
source( "../MCMC.r" )

fLogL <- function( vY, mX, mW, vBeta, dRho, dSigma2 )
{
	cN <- length( vY )

	vE    <- vY - dRho * mW %*% vY - mX %*% vBeta
	dLogL <- log( det( diag( cN ) - dRho * mW ) ) - 0.5 * sum( vE ^ 2 ) / dSigma2

	return( dLogL )
}

fEffect <- function( vGamma, vDelta, dRho, mW )
{
	cN <- nrow( mW )
	cK <- length( vGamma )

	mI <- diag( cN )

	vEffect <- NULL
	for ( i in 1 : cK )
	{
		mS        <- solve( mI - dRho * mW ) %*% ( vGamma[i] * mI + vDelta[i] * mW )
		dTotal    <- sum( mS ) / cN
		dDirect   <- sum( diag( mS ) ) / cN
		dIndirect <-  dTotal - dDirect
		vEffect   <- c( vEffect, dTotal, dDirect, dIndirect )
	}

	return( vEffect )
}

## set timer
ttt <- proc.time(  )

## load data
sData <- read.csv( "../Data/mi.csv" )
vZ.0  <- sData$Y
vZ.1  <- sData$K
vZ.2  <- sData$L
mC    <- cbind( sData$北海道, sData$青森, sData$岩手, sData$宮城, sData$秋田, sData$山形, sData$福島, sData$茨城, sData$栃木, sData$群馬, sData$埼玉, sData$千葉, sData$東京, sData$神奈川, sData$新潟, sData$富山, sData$石川, sData$福井, sData$山梨, sData$長野, sData$岐阜, sData$静岡, sData$愛知, sData$三重, sData$滋賀, sData$京都, sData$大阪, sData$兵庫, sData$奈良, sData$和歌山, sData$鳥取, sData$島根, sData$岡山, sData$広島, sData$山口, sData$徳島, sData$香川, sData$愛媛, sData$高知, sData$福岡, sData$佐賀, sData$長崎, sData$熊本, sData$大分, sData$宮崎, sData$鹿児島, sData$沖縄 )
vS    <- apply( mC, 1, sum )
vSum  <- ifelse( vS == 0, 1, vS )
mW    <- mC / vSum
vY    <- log( vZ.0 ) - log( vZ.2 )
vX.1  <- log( vZ.1 ) - log( vZ.2 )
vX.2  <- log( vZ.2 )
mZ    <- cbind( vX.1, vX.2 )
mX    <- cbind( 1, mZ, mW %*% mZ )

## set dimensions
cN  <- nrow( mX )
cK  <- ncol( mX )
cKK <- ncol( mZ )

## set initial values
vBeta   <- rep( 0.0, cK )
dRho    <- 0.6
dSigma2 <- 1.0

## set hyper-parameters
vBeta.0  <- rep( 0.0, cK )
mSigma.0 <- diag( 100.0, cK )
dNu.0    <- 2.0
dLamb.0  <- 0.01

## set burn-in & draws
iBurn <- 1000000
iDraw <- 4000000
iIter <- iBurn + iDraw
nstep <- as.integer( 0.1 * iIter )

## set scalars, vectors & matrices
mInvSigma.0 <- solve( mSigma.0 )
vISB.0      <- mInvSigma.0 %*% vBeta.0

mI   <- diag( cN )
t.mX <- t( mX )
mXX  <- t.mX %*% mX
vWY  <- mW %*% vY

dNu.1 <- cN + dNu.0

dC      <- 1.0
dCount  <- 0
dAccept <- 0
dMin    <- 1.0 / min( eigen( mW )$values )

## set posterior matrices
mPostParam  <- matrix( 0, iDraw, cK + 2 )
mPostEffect <- matrix( 0, iDraw, 3 * cKK )
vPostL      <- rep( 0.0, iDraw )

## MCMC main loop
print( "main loop:" )
for ( iter in 1 : iIter )
{
	## monitoring
	if ( iter %% nstep == 0 )
	{
		print( paste( "iter = ", sprintf( "%7d", iter - iBurn ), ", vBeta[1] = ", sprintf( "%7.3f", vBeta[1] ), ", dC = ", sprintf( "%7.3f", dC ), sep = "" ) )
	}

	## sampling vBeta
	mSigma.1 <- solve( mXX / dSigma2 + mInvSigma.0 )
	vBeta.1  <- mSigma.1 %*% ( t.mX %*% ( vY - dRho * vWY ) / dSigma2 + vISB.0 )

	vBeta <- vBeta.1 + t( chol( mSigma.1 ) ) %*% rnorm( cK )

	## sampling dSigma2
	dLamb.1 <- sum( ( vY - dRho * vWY - mX %*% vBeta ) ^ 2 ) + dLamb.0

	dSigma2 <- 1.0 / rgamma( 1, 0.5 * dNu.1, 0.5 * dLamb.1 )

	## sampling dRho
	dCount <- dCount + 1
	dRho.n <- dRho + dC * rnorm( 1 )
	if ( dRho.n > dMin & dRho.n < 1.0 )
	{
		dLogL.o <- fLogL( vY, mX, mW, vBeta, dRho, dSigma2 )
		dLogL.n <- fLogL( vY, mX, mW, vBeta, dRho.n, dSigma2 )
		if ( runif( 1 ) < exp( dLogL.n - dLogL.o ) )
		{
			dRho    <- dRho.n
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
	if ( iter - iBurn == 0 )
	{
	  dCount  <- 0
	  dAccept <- 0
	}

	## save parameters
	if ( iter > iBurn )
	{
		it <- iter - iBurn

		vGamma <- vBeta[2 : ( cKK + 1 )]
		vDelta <- vBeta[( cKK + 2 ) : cK]

		mPostParam[it, ]  <- c( vBeta, dRho, dSigma2 )
		mPostEffect[it, ] <- fEffect( vGamma, vDelta, dRho, mW )
		vPostL[it]        <- det( diag( cN ) - dRho * mW ) * prod( dnorm( ( vY - dRho * vWY - mX %*% vBeta ) / sqrt( dSigma2 ) ) / sqrt( dSigma2 ) )
	}
}
k <- seq( 1, iDraw, by = 100 )
mPostParam  <- mPostParam[k, ]
mPostEffect <- mPostEffect[k, ]
vPostL      <- vPostL[k]

## set parameter names
sName   <- NULL
for ( i in 1 : ( cK + 2 ) )
{
	if ( i == 1 )
	{
		sName <- c( sName, bquote( alpha ) )
	}
	else if ( i <= cKK + 1 )
	{
		sName <- c( sName, bquote( beta[.( i - 1 )] ) )
	}
	else if ( i <= cK )
	{
		sName <- c( sName, bquote( gamma[.( i - cKK - 1 )] ) )
	}
	else if ( i == cK + 1 )
	{
		sName <- c( sName, bquote( rho ) )
	}
	else
	{
		sName <- c( sName, bquote( sigma ^ 2 ) )
	}
}
sName.E <- NULL
j <- 1
for ( i in 1 : ( 3 * cKK ) )
{
	if ( i %% 3 == 1 )
	{
		sName.E <- c( sName.E, bquote( Total[.( j )] ) )
	}
	else if ( i %% 3 == 2 )
	{
		sName.E <- c( sName.E, bquote( Direct[.( j )] ) )
	}
	else
	{
		sName.E <- c( sName.E, bquote( Indirect[.( j )] ) )
		j <- j + 1
	}
}

## posterior analysis
mSummary   <- fSummary( mPostParam, sName, "SDM" )
mSummary.E <- fSummary( mPostEffect, sName.E, "SDM-E" )
print( paste( "acceptance rate: ", sprintf( "%6.3f", dAccept / dCount ), sep = "" ) )
dLogL <- log( 1.0 / mean( 1.0 / vPostL ) )
print( paste( "log marginal likelihood = ", dLogL, sep = "" ) )
print( fHarmonicMean( vPostL ) )
fDraw( mPostParam, sName, "SDM-" )
fDraw( mPostEffect, sName.E, "SDM-E-" )

print( paste( "Execution time: ", ( proc.time(  ) - ttt )[3], sep = "" ) )
