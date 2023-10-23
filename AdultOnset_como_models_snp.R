gen_cor_norm <- function( n, mu, sigma, rho, pop=F ) {
	
	t <- length( mu )
	X <- matrix( rnorm( t*n ), n, t )
	X.sig <- var( scale( X ) )
	X.sig.C.inv <- solve( chol( X.sig ) )
	X.clean <- scale( X %*% X.sig.C.inv )
	rho.C <- chol( rho )
	Y <- matrix( rep( mu, n ), n, t, byrow=T ) + 
			( X.clean %*% rho.C ) * matrix( rep( sigma, n ), n, t, byrow=T )

	return( Y )
	
}

plot.profiles.age <- function( l.value, p.value, plotName, traits=c("Trait 1", "Trait 2", "Trait 3") ) {
	
	l.1.con <- mean( l.value[ p.value[,1] == 0, 1 ] )
	l.2.con <- mean( l.value[ p.value[,1] == 0, 2 ] )
	l.3.con <- mean( l.value[ p.value[,1] == 0, 3 ] )
			
	l.1.ch <- mean( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 1 ] )
	l.2.ch <- mean( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 2 ] )
	l.3.ch <- mean( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 3 ] )
	
	l.1.adult <- mean( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 1 ] ) 
	l.2.adult <- mean( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 2 ] )
	l.3.adult <- mean( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 3 ] )

	sd.1.con <- sd( l.value[ p.value[,1] == 0, 1 ] ) / 
						sqrt( sum( p.value[,1] == 0 ) )
	sd.2.con <- sd( l.value[ p.value[,1] == 0, 2 ] ) / 
						sqrt( sum( p.value[,1] == 0 ) )
	sd.3.con <- sd( l.value[ p.value[,1] == 0, 3 ] ) / 
						sqrt( sum( p.value[,1] == 0 ) )
			
	sd.1.ch <- sd( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 1 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,4] < 18 ) )
	sd.2.ch <- sd( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 2 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,4] < 18 ) )
	sd.3.ch <- sd( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 3 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,4] < 18 ) )
	
	sd.1.adult <- sd( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 1 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,4] > 18 ) ) 
	sd.2.adult <- sd( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 2 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,4] > 18 ) )
	sd.3.adult <- sd( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 3 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,4] > 18 ) )


	plot.dat <- rbind( 	c( l.1.con, l.1.adult, l.1.ch ),
						c( l.2.con, l.2.adult, l.2.ch ),
						c( l.3.con, l.3.adult, l.3.ch ) )
	plot.sd <- rbind( 	c( sd.1.con, sd.1.adult, sd.1.ch ),
						c( sd.2.con, sd.2.adult, sd.2.ch ),
						c( sd.3.con, sd.3.adult, sd.3.ch ) )
	
	col.cats <- c( "#FF62BC", "#FF62BC", "#00BFC4" )
	dense.cats <- c( 0, NA, 40 )
	angle.cats <- c( NA, NA, 0 )
	
	pdf( plotName )	
	
	b.plot <- print( barplot( t( plot.dat ),
					names.arg=traits, 
					beside=T,
					angle=angle.cats, density=dense.cats, 
					ylab="Expected Genetic Value",
					col=rep( c( "#FF62BC", "#FF62BC", "#00BFC4" ), each=3 ),
#					ylim=c( min( plot.dat ) - 0.1*( max(plot.dat)-min(plot.dat) ),
#							max( plot.dat ) + 0.1*( max(plot.dat)-min(plot.dat) ) )
					ylim=c(-0.75,2)
				) )
	print( segments( b.plot, t( plot.dat )-1.96*t( plot.sd ), 
			b.plot, t( plot.dat )+1.96*t( plot.sd ) ) )
	print( segments( b.plot-0.25, t( plot.dat )+1.96*t( plot.sd ), 
			b.plot+0.25, t( plot.dat )+1.96*t( plot.sd ) ) )
	print( segments( b.plot-0.25, t( plot.dat )-1.96*t( plot.sd ), 
			b.plot+0.25, t( plot.dat )-1.96*t( plot.sd ) ) )

	print( legend( 'topright', legend=c("ADHD-", "1st dx Adult", "1 dx Child"),
					angle=angle.cats, 
					density=dense.cats, 
					col=rep( c( "#FF62BC", "#FF62BC", "#00BFC4" ), each=3 ) ) )

	dev.off()
	
	out.dat <- list( plot.dat, plot.sd, b.plot )	
	
	return( out.dat )
	
}


plot.profiles.age.bnw <- function( l.value, p.value, plotName, traits=c("Trait 1", "Trait 2", "Trait 3") ) {
	
	l.1.con <- mean( l.value[ p.value[,1] == 0, 1 ] )
	l.2.con <- mean( l.value[ p.value[,1] == 0, 2 ] )
	l.3.con <- mean( l.value[ p.value[,1] == 0, 3 ] )
			
	l.1.ch <- mean( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 1 ] )
	l.2.ch <- mean( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 2 ] )
	l.3.ch <- mean( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 3 ] )
	
	l.1.adult <- mean( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 1 ] ) 
	l.2.adult <- mean( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 2 ] )
	l.3.adult <- mean( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 3 ] )

	q05.1.con <- quantile( l.value[ p.value[,1] == 0, 1 ], 0.05 )
	q05.2.con <- quantile( l.value[ p.value[,1] == 0, 2 ], 0.05 )
	q05.3.con <- quantile( l.value[ p.value[,1] == 0, 3 ], 0.05 )
			
	q05.1.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 1 ], 0.05 )
	q05.2.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 2 ], 0.05 )
	q05.3.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 3 ], 0.05 )
	
	q05.1.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 1 ], 0.05 ) 
	q05.2.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 2 ], 0.05 )
	q05.3.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 3 ], 0.05 )

	q25.1.con <- quantile( l.value[ p.value[,1] == 0, 1 ], 0.25 )
	q25.2.con <- quantile( l.value[ p.value[,1] == 0, 2 ], 0.25 )
	q25.3.con <- quantile( l.value[ p.value[,1] == 0, 3 ], 0.25 )
			
	q25.1.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 1 ], 0.25 )
	q25.2.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 2 ], 0.25 )
	q25.3.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 3 ], 0.25 )
	
	q25.1.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 1 ], 0.25 ) 
	q25.2.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 2 ], 0.25 )
	q25.3.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 3 ], 0.25 )

	q75.1.con <- quantile( l.value[ p.value[,1] == 0, 1 ], 0.75 )
	q75.2.con <- quantile( l.value[ p.value[,1] == 0, 2 ], 0.75 )
	q75.3.con <- quantile( l.value[ p.value[,1] == 0, 3 ], 0.75 )
			
	q75.1.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 1 ], 0.75 )
	q75.2.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 2 ], 0.75 )
	q75.3.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 3 ], 0.75 )
	
	q75.1.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 1 ], 0.75 ) 
	q75.2.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 2 ], 0.75 )
	q75.3.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 3 ], 0.75 )

	q95.1.con <- quantile( l.value[ p.value[,1] == 0, 1 ], 0.95 )
	q95.2.con <- quantile( l.value[ p.value[,1] == 0, 2 ], 0.95 )
	q95.3.con <- quantile( l.value[ p.value[,1] == 0, 3 ], 0.95 )
			
	q95.1.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 1 ], 0.95 )
	q95.2.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 2 ], 0.95 )
	q95.3.ch <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] < 18, 3 ], 0.95 )
	
	q95.1.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 1 ], 0.95 ) 
	q95.2.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 2 ], 0.95 )
	q95.3.adult <- quantile( l.value[ p.value[,1] == 1 & p.value[,4] > 18, 3 ], 0.95 )
	
	plot.dat <- rbind( 	c( l.1.con, l.1.adult, l.1.ch ),
						c( l.2.con, l.2.adult, l.2.ch ),
						c( l.3.con, l.3.adult, l.3.ch ) )
	plot.q05 <- rbind( 	c( q05.1.con, q05.1.adult, q05.1.ch ),
						c( q05.2.con, q05.2.adult, q05.2.ch ),
						c( q05.3.con, q05.3.adult, q05.3.ch ) )
	plot.q25 <- rbind( 	c( q25.1.con, q25.1.adult, q25.1.ch ),
						c( q25.2.con, q25.2.adult, q25.2.ch ),
						c( q25.3.con, q25.3.adult, q25.3.ch ) )
	plot.q75 <- rbind( 	c( q75.1.con, q75.1.adult, q75.1.ch ),
						c( q75.2.con, q75.2.adult, q75.2.ch ),
						c( q75.3.con, q75.3.adult, q75.3.ch ) )	
	plot.q95 <- rbind( 	c( q95.1.con, q95.1.adult, q95.1.ch ),
						c( q95.2.con, q95.2.adult, q95.2.ch ),
						c( q95.3.con, q95.3.adult, q95.3.ch ) )
																		
	col.cats <- c( "#FF62BC", "#FF62BC", "#00BFC4" )
	dense.cats <- c( 0, NA, 40 )
	angle.cats <- c( NA, NA, 0 )
	
	pdf( plotName )	
	
	b.plot <- print( barplot( t( plot.dat ),
					names.arg=traits, 
					beside=T,
					angle=angle.cats, density=dense.cats, 
					ylab="Expected Genetic Value",
					col='white', border='white',
#					ylim=c( min( plot.dat ) - 0.1*( max(plot.dat)-min(plot.dat) ),
#							max( plot.dat ) + 0.1*( max(plot.dat)-min(plot.dat) ) )
					ylim=c(-2,3), type='n'		
				) )

	abline( h=0 )

	for ( i in 1:3 ) {
	for ( g in 1:3 ) {
		
		t.col=col.cats[ i ]
		t.border=col.cats[ i ]
		t.angle=angle.cats[ g ]
		t.density=dense.cats[ g ]	

		tempMat <- t( plot.dat )
		tempMat.q05 <- t( plot.q05 )
		tempMat.q25 <- t( plot.q25 )
		tempMat.q75 <- t( plot.q75 )
		tempMat.q95 <- t( plot.q95 )
				
		polygon(	c( b.plot[ g,i ]-0.4, b.plot[ g,i ]+0.4, b.plot[ g,i ]+0.4, b.plot[ g,i ]-0.4 ), 
					c( tempMat.q25[ g,i ], tempMat.q25[ g,i ], tempMat.q75[ g,i ], tempMat.q75[ g,i ] ), 
					col=t.col,
					border=t.border,
					angle=t.angle,
					density=t.density
					)
		segments(	b.plot[ g,i ], tempMat.q75[ g,i ], 
					b.plot[ g,i ], tempMat.q95[ g,i ], 
					col=t.col
					)
		segments(	b.plot[ g,i ], tempMat.q25[ g,i ], 
					b.plot[ g,i ], tempMat.q05[ g,i ], 
					col=t.col
					)							
		segments(	b.plot[ g,i ]+0.3, tempMat.q05[ g,i ], 
					b.plot[ g,i ]-0.3, tempMat.q05[ g,i ], 
					col=t.col
					)	
		segments(	b.plot[ g,i ]+0.3, tempMat.q95[ g,i ], 
					b.plot[ g,i ]-0.3, tempMat.q95[ g,i ], 
					col=t.col
					)
		segments(	b.plot[ g,i ]+0.4, tempMat[ g,i ], 
					b.plot[ g,i ]-0.4, tempMat[ g,i ], 
					col=1
					)					
	}
	}


	print( legend( 'topright', legend=c("ADHD-", "1st dx Adult", "1 dx Child"),
					angle=angle.cats, 
					density=dense.cats, 
					col=rep( c( "#FF62BC", "#FF62BC", "#00BFC4" ), each=3 ) ) )

	dev.off()
	
	out.dat <- list( plot.dat, plot.sd, b.plot )	
	
	return( out.dat )
	
}

## Basic parameters

traits <- c( "ADHD", "MDD", "EA" )
n=100000
m=250
k=c( 0.05, 0.15, NA )
h2=c( 0.2, 0.2, 0.2 )
rg=matrix( c( 1, 0.4, -0.5, 0.4, 1, -0.25, -0.5, -0.25, 1), 3, 3 )
re=diag( length(h2) )
onset <- c( 5, 15, NA )		#ADHD uniform: 5 to 20, MDD: uniform 15 to 65
offset <- c( 20, 65, NA)
t <- length( h2 )

## Model 0: Add random, uniform age of onset of ADHD

# Genotypes
g <- matrix( rnorm( n*m ), n, m )
beta <- gen_cor_norm( m, rep(0,t), sqrt(h2/m), rg, pop=F )
# Genetic values
g.value.0 <- g %*% beta
# Environmental values				
e.value.0 <- gen_cor_norm( n, rep(0,t), sqrt(1-h2), re, pop=F )
# Liability values
l.value.0 <- g.value.0 + e.value.0
# Convert liability to 'disorder' phenotypes
p.value.0 <- l.value.0
p.value.0[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.0[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )

# Add onset
p.value.0 <- cbind( p.value.0, 0 )
p.value.0 <- cbind( p.value.0, 0 )
p.value.0[ , 4 ] <- sample( seq( onset[1], offset[1], by=5 ),
											n, replace=T )
p.value.0[ , 5 ] <- sample( seq( onset[2], offset[2], by=5 ),
											n, replace=T )

#Plot baseline
plot.0 <- plot.profiles.age.bnw( g.value.0, p.value.0,"AdultOnset_Model0_snp.bnw.pdf", traits )

# Model 1a: Add age of onset of ADHD associated with ADHD genetics ( younger = more )
p.value.1a <- p.value.0
p.value.1a[ ,4 ] <- 0
p.value.1a[ g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ),1 ], 0.25 ), 4 ] <- 
			sample( c( 20,20,20,15,10,5), 
#			sample( c( 20,20 ), 			
					replace=T,
					size=sum( g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.25 ) ) 
			)
p.value.1a[ g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.25 ) & 
			g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ), 4 ] <- 
			sample( c( 20,15,15,15,10,5), 
#			sample( c( 15,15 ), 
					replace=T,
					size=sum( g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.25 ) & 
							g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ) 
					) 
			)
p.value.1a[ g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ) & 
			g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ), 4 ] <-
			sample( c( 20,15,10,10,10,5), 
#			sample( c( 10,10 ), 
					replace=T,
					sum( g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ) & 
							g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ) 
					) 
			)
p.value.1a[ g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ), 4 ] <-  
			sample( c( 20,15,10,5,5,5), 
#			sample( c( 5,5 ),
					replace=T,
					size=sum( g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ) ) 
			)

plot.1a <- plot.profiles.age.bnw( g.value.0, p.value.1a,"AdultOnset_Model1a_snp.bnw.pdf", traits )

# Model 1b: Add age of onset of ADHD associated with ADHD genetics ( younger = less )
p.value.1b <- p.value.0
p.value.1b[ ,4 ] <- 0
p.value.1b[ g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.25 ), 4 ] <- 
			sample( c( 20,15,10,5,5,5), 
#			sample( c( 5,5 ),
					replace=T,
					size=sum( g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.25 ) ) 
			)
p.value.1b[ g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.25 ) & 
			g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ), 4 ] <- 
			sample( c( 20,15,10,10,10,5), 
#			sample( c( 10,10 ), 
					replace=T,
					size=sum( g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.25 ) & 
							g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ) 
					) 
			)
p.value.1b[ g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ) & 
			g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ), 4 ] <-
			sample( c( 20,15,15,15,10,5), 
#			sample( c( 15,15 ), 
					replace=T,
					sum( g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.5 ) & 
							g.value.0[,1] < quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ) 
					) 
			)
p.value.1b[ g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ), 4 ] <-  
			sample( c( 20,20,20,15,10,5), 
#			sample( c( 20,20 ), 			
					replace=T,
					size=sum( g.value.0[,1] >= quantile( g.value.0[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ], 0.75 ) ) 
			)

plot.1b <- plot.profiles.age.bnw( g.value.0, p.value.1b,"AdultOnset_Model1b_snp.bnw.pdf", traits )

# Model 2: Add childhood onset, and adult including misdiagnosed

#rg.2=matrix( c( 1, 0.4, -0.5, 0.4, 1, -0.25, -0.5, -0.25, 1), 3, 3 )
#rg.2=matrix( c( 1, 0.4, -0.33, 0.4, 1, -0.33, -0.33, -0.33, 1), 3, 3 )
#rg.2=matrix( c( 1, 0.4, -0.25, 0.4, 1, -0.5, -0.25, -0.5, 1), 3, 3 )
#g.value.2 <- gen_cor_norm( n, rep(0,t), sqrt(h2), rg.2, pop=F )
#e.value.2 <- gen_cor_norm( n, rep(0,t), sqrt(1-h2), re, pop=F )
#l.value.2 <- g.value.2 + e.value.2
#p.value.2 <- l.value.2
#p.value.2[ ,1 ] <- 1*( l.value.2[,1] >= qnorm( 1-k[1] ) )
#p.value.2[ ,2 ] <- 1*( l.value.2[,2] >= qnorm( 1-k[2] ) )

p.value.2 <- l.value.0
p.value.2[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.2[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )

p.value.2 <- cbind( p.value.2, 0 )
p.value.2 <- cbind( p.value.2, 0 )
p.value.2[ ,4 ] <- sample( seq( onset[1], offset[1], by=5 ),
											n, replace=T )
p.value.2[ ,5 ] <- sample( seq( onset[2], offset[2], by=5 ),
											n, replace=T )
											
misDx <- sample( which( p.value.2[,2] == 1 & p.value.2[,4] == 20 ), 
					sum( p.value.2[,2] == 1 & p.value.2[,4] == 20 )/10 )

p.value.2[ misDx, 1 ] <- 1
p.value.2[ misDx, 4 ] <- 20

#plot.2 <- plot.profiles.age( g.value.2, p.value.2,"AdultOnset_Model2.pdf" )
plot.2 <- plot.profiles.age.bnw( g.value.0, p.value.2,"AdultOnset_Model2_snp.bnw.pdf", traits )

# Model 3: Adult have more environmental risk, that is related to EA (active GxE correlation?)
# assign age of onset "pre-diagnosis"
l.value.3 <- l.value.0
l.value.3 <- cbind( l.value.3,0 )
l.value.3 <- cbind( l.value.3,0 )
l.value.3[ ,4 ] <- sample( seq( onset[1], offset[1], by=5 ),
											n, replace=T )
l.value.3[ ,5 ] <- sample( seq( onset[2], offset[2], by=5 ),
											n, replace=T )
											
# Add opposite of EA liability to ADHD, as a function of age 
l.value.3[ l.value.3[,4] == 5, 1 ] <- l.value.3[ l.value.3[,4] == 5, 1 ] - 
										0*l.value.3[ l.value.3[,4] == 5, 3 ]
l.value.3[ l.value.3[,4] == 10, 1 ] <- l.value.3[ l.value.3[,4] == 10, 1 ] - 
										0*l.value.3[ l.value.3[,4] == 10, 3 ]
l.value.3[ l.value.3[,4] == 15, 1 ] <- l.value.3[ l.value.3[,4] == 15, 1 ] - 
										0*l.value.3[ l.value.3[,4] == 15, 3 ]
l.value.3[ l.value.3[,4] == 20, 1 ] <- l.value.3[ l.value.3[,4] == 20, 1 ] - 
										1*l.value.3[ l.value.3[,4] == 20, 3 ]

# Convert to dx							
p.value.3 <- l.value.3
p.value.3[ ,1 ] <- 1*( l.value.3[,1] >= qnorm( 1-k[1] ) )
p.value.3[ ,2 ] <- 1*( l.value.3[,2] >= qnorm( 1-k[2] ) )

plot.3 <- plot.profiles.age.bnw( g.value.0, p.value.3,"AdultOnset_Model3_snp.bnw.pdf", traits )

# Model 4: Adult have more random environmental effects (i.e., onset function of e.value)
p.value.4 <- p.value.0
p.value.4[ ,4 ] <- 0
p.value.4[ e.value.0[,1] < quantile( e.value.0[,1], 0.25 ), 4 ] <- 5
p.value.4[ e.value.0[,1] >= quantile( e.value.0[,1], 0.25 ) & 
			e.value.0[,1] < quantile( e.value.0[,1], 0.5 ), 4 ] <- 10
p.value.4[ e.value.0[,1] >= quantile( e.value.0[,1], 0.5 ) & 
			e.value.0[,1] < quantile( e.value.0[,1], 0.75 ), 4 ] <- 15
p.value.4[ e.value.0[,1] >= quantile( e.value.0[,1], 0.75 ) & 
			e.value.0[,1] < quantile( e.value.0[,1], 1 ), 4 ] <- 20

plot.4 <- plot.profiles.age.bnw( g.value.0, p.value.4,"AdultOnset_Model4_snp.bnw.pdf", traits )




