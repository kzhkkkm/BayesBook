rm( list = ls(  ) )
set.seed( 123456 )

## read csv-data
sData <- read.csv( "../Data/ExRate.csv", header = T )
vY    <- 100 * diff( log( sData$ドル ) )

## set dimension
cT <- length( vY )

print( paste( "hv = ", sd( vY ), sep = "" ) )

## Posterior Summary
dNu.1   <- 0.5 * cT
dLamb.1 <- 0.5 * sum( ( vY - mean( vY ) ) ^ 2 )
dX.l <- sqrt( 1.0 / qgamma( 0.975, dNu.1, dLamb.1 ) )
dX.u <- sqrt( 1.0 / qgamma( 0.025, dNu.1, dLamb.1 ) )
print( paste( "cT = ", cT, sep = "" ) )
print( paste( "dLamb.1 = ", sum( ( vY - mean( vY ) ) ^ 2 ), sep = "" ) )
print( c( dX.l, dX.u ) )
print( exp( 0.5 * log( dLamb.1 ) + lgamma( dNu.1 - 0.5 ) - lgamma( dNu.1 ) ) )

## draw figures
inter  <- seq( 1, cT, by = 200 )
pdf( file = "hv-tsplot.pdf", width = 8, height = 5 )
par( mfrow = c( 2, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
y.min <- 0.1 * floor( 10 * min( sData$ドル ) )
y.max <- 0.1 * ceiling( 10 * max( sData$ドル ) )
ts.plot( sData$ドル, main = "終値", col = grey( 0.2 ), xlab = "", ylab = "", ylim = c( y.min, y.max ), gpars = list( xaxt = "n", bty = "n" ) )
axis( 1, at = ( 1 + inter ), labels = sData$日付[( 1 + inter )] )
y.min <- 0.1 * floor( 10 * min( 100 * diff( log( sData$ドル ) ) ) )
y.max <- 0.1 * ceiling( 10 * max( 100 * diff( log( sData$ドル ) ) ) )
ts.plot( 100 * diff( log( sData$ドル ) ), col = grey( 0.2 ), main = "対数変化率",  xlab = "", ylab = "", ylim = c( y.min, y.max ), gpars = list( xaxt = "n", bty = "n" ) )
axis( 1, at = inter, labels = sData$日付[( 1 + inter )] )
dev.off(  )

pdf( file = "hv.pdf", width = 8, height = 5 )
par( mfrow = c( 1, 1 ), mar = c( 2, 2, 1.5, 0.5 ), family = "Japan1Ryumin" )
x <- seq( 0.5, 0.6, by = 0.00001 )
y <- 2 * x ^ ( - 3 ) * dgamma( ( 1 / ( x ^ 2 ) ), dNu.1, dLamb.1 )
plot.new(  )
plot.window( xlim = c( 0.5, 0.6 ), ylim = c( 0, 36 ) )
polygon( c( dX.l, x[x >= dX.l & x <= dX.u], dX.u ), c( 0, y[x >= dX.l & x <= dX.u], 0 ), col = grey( 0.8 ), lty = 0 )
par( new = T )
plot( x, y, type = "l", xlim = c( 0.5, 0.6 ), ylim = c( 0, 36 ), main = "", xlab = "", ylab = "", bty = "n" )
lines( c( 0.5, 0.6 ), c( 0, 0 ) )
lines( c( dX.l, dX.l ), c( 0, 2 * dX.l ^ ( - 3 ) * dgamma( 1 / ( dX.l ^ 2 ), dNu.1, dLamb.1 ) ) )
lines( c( dX.u, dX.u ), c( 0, 2 * dX.u ^ ( - 3 ) * dgamma( 1 / ( dX.u ^ 2 ), dNu.1, dLamb.1 ) ) )
arrows( 0.53, 14, dX.l, 0, length = 0.05 )
arrows( 0.59, 14, dX.u, 0, length = 0.05 )
text( 0.53, 15, round( dX.l, digit = 3 ) )
text( 0.59, 15, round( dX.u, digit = 3 ) )
dev.off(  )
