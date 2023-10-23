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

plot.profiles.como <- function( l.value, p.value, plotName, 
								traits.txt=c( "Trait 1","Trait 2","Trait 3" ),
								legend.txt=c( "Group 0-0", "Group 1-1", "Group 1-0", "Group 0-1") ) {
	
	l.1.00 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 1 ] )
	l.2.00 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 2 ] )
	l.3.00 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 3 ] )
			
	l.1.10 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 1 ] )
	l.2.10 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 2 ] )
	l.3.10 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 3 ] )
	
	l.1.01 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 1 ] )
	l.2.01 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 2 ] )
	l.3.01 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 3 ] )

	l.1.11 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 1 ] )
	l.2.11 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 2 ] )
	l.3.11 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 3 ] )

	sd.1.00 <- sd( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 1 ] ) / 
						sqrt( sum( p.value[,1] == 0 & p.value[,2] == 0 ) )
	sd.2.00 <- sd( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 2 ] ) / 
						sqrt( sum( p.value[,1] == 0 & p.value[,2] == 0 ) )
	sd.3.00 <- sd( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 3 ] ) / 
						sqrt( sum( p.value[,1] == 0 & p.value[,2] == 0 ) )
			
	sd.1.10 <- sd( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 1 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,2] == 0 ) )
	sd.2.10 <- sd( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 2 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,2] == 0 ) )
	sd.3.10 <- sd( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 3 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,2] == 0 ) )
						
	sd.1.01 <- sd( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 1 ] ) / 
						sqrt( sum( p.value[,1] == 0 & p.value[,2] == 1 ) )
	sd.2.01 <- sd( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 2 ] ) / 
						sqrt( sum( p.value[,1] == 0 & p.value[,2] == 1 ) )
	sd.3.01 <- sd( l.value[ p.value[,1] == 0 & p.value[,2] == 1 , 3 ]) / 
						sqrt( sum( p.value[,1] == 0 & p.value[,2] == 1 ) )

	sd.1.11 <- sd( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 1 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,2] == 1 ) )
	sd.2.11 <- sd( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 2 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,2] == 1 ) )
	sd.3.11 <- sd( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 3 ] ) / 
						sqrt( sum( p.value[,1] == 1 & p.value[,2] == 1 ) )

	plot.dat <- rbind( 	c( l.1.00, l.1.11, l.1.10, l.1.01 ),
						c( l.2.00, l.2.11, l.2.10, l.2.01 ),
						c( l.3.00, l.3.11, l.3.10, l.3.01 ) )
	plot.sd <- rbind( 	c( sd.1.00, sd.1.11, sd.1.10, sd.1.01 ),
						c( sd.2.00, sd.2.11, sd.2.10, sd.2.01 ),
						c( sd.3.00, sd.3.11, sd.3.10, sd.3.01 ) )
	
	col.cats <- c( "#FF62BC", "#FF62BC", "#00BFC4" )
	dense.cats <- c( 0, NA, 40, 40 )
	angle.cats <- c( NA, NA, 0, 135 )
	
	pdf( plotName )	
	
	b.plot <- print( barplot( t( plot.dat ),
					names.arg=traits.txt, 
					beside=T,
					angle=angle.cats, density=dense.cats, 
					ylab="Expected Genetic Value",
					col=rep( col.cats, each=4 ),
					ylim=c( -1,3 )
				) )
	print( segments( b.plot, t( plot.dat )-1.96*t( plot.sd ), 
			b.plot, t( plot.dat )+1.96*t( plot.sd ) ) )
	print( segments( b.plot-0.25, t( plot.dat )+1.96*t( plot.sd ), 
			b.plot+0.25, t( plot.dat )+1.96*t( plot.sd ) ) )
	print( segments( b.plot-0.25, t( plot.dat )-1.96*t( plot.sd ), 
			b.plot+0.25, t( plot.dat )-1.96*t( plot.sd ) ) )

	print( legend( 'topright', 
					legend=legend.txt,
					angle=angle.cats, 
					density=dense.cats, 
					col=1 ) )

	dev.off()
	
	out.dat <- list( plot.dat, plot.sd, b.plot )	
	
	return( out.dat )
	
}


plot.profiles.como.bnw <- function( l.value, p.value, plotName, 
								traits.txt=c( "Trait 1","Trait 2","Trait 3" ),
								legend.txt=c( "Group 0-0", "Group 1-1", "Group 1-0", "Group 0-1") ) {

	#l.value <- g.value.0
	#p.value <- p.value.0 
	#traits.txt <- traits
	#legend.txt <- legend.txt
	
	l.1.00 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 1 ] )
	l.2.00 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 2 ] )
	l.3.00 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 3 ] )
			
	l.1.10 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 1 ] )
	l.2.10 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 2 ] )
	l.3.10 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 3 ] )
	
	l.1.01 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 1 ] )
	l.2.01 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 2 ] )
	l.3.01 <- mean( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 3 ] )

	l.1.11 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 1 ] )
	l.2.11 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 2 ] )
	l.3.11 <- mean( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 3 ] )

	q05.1.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 1 ], 0.05 )
	q05.2.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 2 ], 0.05 )
	q05.3.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 3 ], 0.05 )
				
	q05.1.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 1 ], 0.05 )
	q05.2.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 2 ], 0.05 )
	q05.3.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 3 ], 0.05 )
						
	q05.1.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 1 ], 0.05 )
	q05.2.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 2 ], 0.05 )
	q05.3.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 3 ], 0.05 )

	q05.1.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 1 ], 0.05 )
	q05.2.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 2 ], 0.05 )
	q05.3.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 3 ], 0.05 )

	q25.1.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 1 ], 0.25 )
	q25.2.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 2 ], 0.25 )
	q25.3.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 3 ], 0.25 )
				
	q25.1.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 1 ], 0.25 )
	q25.2.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 2 ], 0.25 )
	q25.3.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 3 ], 0.25 )
						
	q25.1.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 1 ], 0.25 )
	q25.2.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 2 ], 0.25 )
	q25.3.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 3 ], 0.25 )

	q25.1.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 1 ], 0.25 )
	q25.2.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 2 ], 0.25 )
	q25.3.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 3 ], 0.25 )

	q75.1.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 1 ], 0.75 )
	q75.2.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 2 ], 0.75 )
	q75.3.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 3 ], 0.75 )
				
	q75.1.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 1 ], 0.75 )
	q75.2.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 2 ], 0.75 )
	q75.3.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 3 ], 0.75 )
						
	q75.1.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 1 ], 0.75 )
	q75.2.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 2 ], 0.75 )
	q75.3.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 3 ], 0.75 )

	q75.1.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 1 ], 0.75 )
	q75.2.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 2 ], 0.75 )
	q75.3.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 3 ], 0.75 )

	q95.1.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 1 ], 0.95 )
	q95.2.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 2 ], 0.95 )
	q95.3.00 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 0, 3 ], 0.95 )
				
	q95.1.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 1 ], 0.95 )
	q95.2.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 2 ], 0.95 )
	q95.3.10 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 0, 3 ], 0.95 )
						
	q95.1.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 1 ], 0.95 )
	q95.2.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 2 ], 0.95 )
	q95.3.01 <- quantile( l.value[ p.value[,1] == 0 & p.value[,2] == 1, 3 ], 0.95 )

	q95.1.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 1 ], 0.95 )
	q95.2.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 2 ], 0.95 )
	q95.3.11 <- quantile( l.value[ p.value[,1] == 1 & p.value[,2] == 1, 3 ], 0.95 )

	plot.dat <- rbind( 	c( l.1.00, l.1.11, l.1.10, l.1.01 ),
						c( l.2.00, l.2.11, l.2.10, l.2.01 ),
						c( l.3.00, l.3.11, l.3.10, l.3.01 ) )
	plot.q05 <- rbind( 	c( q05.1.00, q05.1.11, q05.1.10, q05.1.01 ),
						c( q05.2.00, q05.2.11, q05.2.10, q05.2.01 ),
						c( q05.3.00, q05.3.11, q05.3.10, q05.3.01 ) )
	plot.q25 <- rbind( 	c( q25.1.00, q25.1.11, q25.1.10, q25.1.01 ),
						c( q25.2.00, q25.2.11, q25.2.10, q25.2.01 ),
						c( q25.3.00, q25.3.11, q25.3.10, q25.3.01 ) )
	plot.q75 <- rbind( 	c( q75.1.00, q75.1.11, q75.1.10, q75.1.01 ),
						c( q75.2.00, q75.2.11, q75.2.10, q75.2.01 ),
						c( q75.3.00, q75.3.11, q75.3.10, q75.3.01 ) )
	plot.q95 <- rbind( 	c( q95.1.00, q95.1.11, q95.1.10, q95.1.01 ),
						c( q95.2.00, q95.2.11, q95.2.10, q95.2.01 ),
						c( q95.3.00, q95.3.11, q95.3.10, q95.3.01 ) )					
			
	col.cats <- c( "#FF62BC", "#FF62BC", "#00BFC4" )
	dense.cats <- c( 0, NA, 40, 40 )
	angle.cats <- c( NA, NA, 0, 135 )
	
	pdf( plotName )	
	
	b.plot <- print( barplot( t( plot.dat ),
					names.arg=traits.txt, 
					beside=T,
					angle=angle.cats, density=dense.cats,
					col='white', border='white',
					ylab="Expected Genetic Value",
					ylim=c( -2,4 )
				) )
	abline( h=0 )

	for ( i in 1:3 ) {
	for ( g in 1:4 ) {
		
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

	print( legend( 'topright', 
					legend=legend.txt,
					angle=angle.cats, 
					density=dense.cats, 
					col=1 ) )

	dev.off()
	
	out.dat <- list( plot.dat, plot.sd, b.plot )	
	
	return( out.dat )
	
}



## Basic parameters 

traits <- c( "ADHD", "ASD", "EA" )
legend.txt <- c( "ADHD-, ASD-", "ADHD+, ASD+", "ADHD+, ASD-", "ADHD-, ASD+" )
n=100000
m=250
k=c( 0.05, 0.01, NA )
h2=c( 0.2, 0.2, 0.2 )
rg=matrix( c( 1, 0.33, -0.5, 0.33, 1, 0.2, -0.5, 0.2, 1), 3, 3 )
re=diag( length(h2) )
onset <- c( 5, 0, NA )		#ADHD uniform: 5 to 20, ASD: uniform 0 to 15
offset <- c( 20, 15, NA)
t <- length( h2 )

## Model 0: Comorbidity

# Genotypes
g <- matrix( rnorm( n*m ), n, m )
beta <- gen_cor_norm( m, rep(0,t), sqrt(h2/m), rg, pop=F )
 
# Genetic values
#g.value.0 <- gen_cor_norm( n, rep(0,t), sqrt(h2), rg, pop=F )
g.value.0 <- g %*% beta
# Environmental values				
e.value.0 <- gen_cor_norm( n, rep(0,t), sqrt(1-h2), re, pop=F )
# Liability values
l.value.0 <- g.value.0 + e.value.0
# Convert liability to 'disorder' phenotypes
p.value.0 <- l.value.0
p.value.0[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.0[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )

# Add random onset
p.value.0 <- cbind( p.value.0, 0 )
p.value.0 <- cbind( p.value.0, 0 )
p.value.0[ , 4 ] <- sample( seq( onset[1], offset[1], by=5 ),
											n, replace=T )
p.value.0[ , 5 ] <- sample( seq( onset[2], offset[2], by=5 ),
											n, replace=T )

#Plot baseline COMO model
plot.0 <- plot.profiles.como.bnw( g.value.0, p.value.0,"ASD_Model0_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1a: Misdiagnosis one-way (ADHD -> ASD)

misDx.r <- 0.10

p.value.1a <- l.value.0
p.value.1a[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1a[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1a[,1] == 1 ), 
					sum( p.value.1a[,1] == 1 )*misDx.r )

p.value.1a[ misDx, 1 ] <- 0
p.value.1a[ misDx, 2 ] <- 1

plot.1a <- plot.profiles.como.bnw( g.value.0, p.value.1a, "ASD_Model1a_flip_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1b: Misdiagnosis one-way (ASD -> ADHD)

misDx.r <- 0.80

p.value.1b <- l.value.0
p.value.1b[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1b[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1b[,2] == 1 ), 
					sum( p.value.1b[,2] == 1 )*misDx.r )

p.value.1b[ misDx, 1 ] <- 1
p.value.1b[ misDx, 2 ] <- 0

plot.1b <- plot.profiles.como.bnw( g.value.0, p.value.1b, "ASD_Model1b_flip_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1c: Misdiagnosis two-way (ASD -> ADHD & ADHD -> ASD) 

misDx.r <- 0.08

p.value.1c <- l.value.0
p.value.1c[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1c[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx1 <- sample( which( p.value.1c[,1] == 1 ), 
					sum( p.value.1c[,1] == 1 )*misDx.r )
misDx2 <- sample( which( p.value.1c[,2] == 1 ), 
					sum( p.value.1c[,2] == 1 )*misDx.r )

p.value.1c[ misDx1, 1 ] <- 0
p.value.1c[ misDx1, 2 ] <- 1
p.value.1c[ misDx2, 1 ] <- 1
p.value.1c[ misDx2, 2 ] <- 0

plot.1c <- plot.profiles.como.bnw( g.value.0, p.value.1c, "ASD_Model1c_flip_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1a: Misdiagnosis one-way (ADHD -> ASD)

misDx.r <- 0.10

p.value.1a <- l.value.0
p.value.1a[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1a[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1a[,1] == 1 ), 
					sum( p.value.1a[,1] == 1 )*misDx.r )

p.value.1a[ misDx, 2 ] <- 1

plot.1a <- plot.profiles.como.bnw( g.value.0, p.value.1a, "ASD_Model1a_add_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1b: Misdiagnosis one-way (ASD -> ADHD)

misDx.r <- 0.80

p.value.1b <- l.value.0
p.value.1b[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1b[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1b[,2] == 1 ), 
					sum( p.value.1b[,2] == 1 )*misDx.r )

p.value.1b[ misDx, 1 ] <- 1

plot.1b <- plot.profiles.como.bnw( g.value.0, p.value.1b, "ASD_Model1b_add_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1c: Misdiagnosis two-way (ASD -> ADHD & ADHD -> ASD) 

misDx.r <- 0.08

p.value.1c <- l.value.0
p.value.1c[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1c[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx1 <- sample( which( p.value.1c[,1] == 1 ), 
					sum( p.value.1c[,1] == 1 )*misDx.r )
misDx2 <- sample( which( p.value.1c[,2] == 1 ), 
					sum( p.value.1c[,2] == 1 )*misDx.r )

p.value.1c[ misDx1, 2 ] <- 1
p.value.1c[ misDx2, 1 ] <- 1

plot.1c <- plot.profiles.como.bnw( g.value.0, p.value.1c, "ASD_Model1c_add_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 2a: One-way causation (being above the threshold for ADHD, lowers the threshold for ASD)

effect <- 0.5

l.value.2a <- l.value.0
l.value.2a[ l.value.0[,1] >= qnorm( 1-k[1] ), 2 ] <- 
				l.value.2a[ l.value.0[,1] >= qnorm( 1-k[1] ), 2 ] + 
				effect*l.value.2a[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ]

p.value.2a <- l.value.2a
p.value.2a[ ,1 ] <- 1*( l.value.2a[,1] >= qnorm( 1-k[1] ) )
p.value.2a[ ,2 ] <- 1*( l.value.2a[,2] >= qnorm( 1-k[2] ) )
											
plot.2a <- plot.profiles.como.bnw( g.value.0, p.value.2a, "ASD_Model2a_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 2b: One-way causation (being above the threshold for ASD, lowers the threshold for ADHD)

effect <- 0.5

l.value.2b <- l.value.0
l.value.2b[ l.value.0[,1] >= qnorm( 1-k[1] ), 2 ] <- 
				l.value.2b[ l.value.0[,1] >= qnorm( 1-k[1] ), 2 ] + 
				effect*l.value.2b[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ]
l.value.2b[ l.value.0[,2] >= qnorm( 1-k[2] ), 1 ] <- 
				l.value.2b[ l.value.0[,2] >= qnorm( 1-k[2] ), 1 ] + 
				effect*l.value.2b[ l.value.0[,2] >= qnorm( 1-k[2] ), 2 ]

p.value.2b <- l.value.2b
p.value.2b[ ,1 ] <- 1*( l.value.2b[,1] >= qnorm( 1-k[1] ) )
p.value.2b[ ,2 ] <- 1*( l.value.2b[,2] >= qnorm( 1-k[2] ) )
											
plot.2b <- plot.profiles.como.bnw( g.value.0, p.value.2b, "ASD_Model2b_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 2c: Mutual causation (being above the threshold for one, lowers the threshold for the other)

effect <- 0.5

l.value.2c <- l.value.0
l.value.2c[ l.value.0[,1] >= qnorm( 1-k[1] ), 2 ] <- 
				l.value.2c[ l.value.0[,1] >= qnorm( 1-k[1] ), 2 ] + 
				effect*l.value.2c[ l.value.0[,1] >= qnorm( 1-k[1] ), 1 ]
l.value.2c[ l.value.0[,2] >= qnorm( 1-k[2] ), 1 ] <- 
				l.value.2c[ l.value.0[,2] >= qnorm( 1-k[2] ), 1 ] + 
				effect*l.value.2c[ l.value.0[,2] >= qnorm( 1-k[2] ), 2 ]

p.value.2c <- l.value.2c
p.value.2c[ ,1 ] <- 1*( l.value.2c[,1] >= qnorm( 1-k[1] ) )
p.value.2c[ ,2 ] <- 1*( l.value.2c[,2] >= qnorm( 1-k[2] ) )
											
plot.2c <- plot.profiles.como.bnw( g.value.0, p.value.2c, "ASD_Model2c_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 2d: Mutual causation (being above the threshold for one, lowers the threshold for the other)

effect <- 5

l.value.2d <- l.value.0

p.value.2d <- l.value.2d
p.value.2d[ ,1 ] <- 1*( l.value.2d[,1] >= qnorm( 1-k[1] ) )
p.value.2d[ ,2 ] <- 1*( l.value.2d[,2] >= qnorm( 1-k[2] ) )

incl.adhd <- which( l.value.2d[ ,1 ] >= qnorm( 1-effect*k[1] ) & p.value.2d[ ,2 ] == 1 )
incl.asd <- which( l.value.2d[ ,2 ] >= qnorm( 1-effect*k[2] ) & p.value.2d[ ,1 ] == 1 )

p.value.2d[ incl.adhd, 1 ] <- 1
p.value.2d[ incl.asd, 2 ] <- 1
											
plot.2d <- plot.profiles.como.bnw( g.value.0, p.value.2d, "ASD_Model2d_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )



## Model 3a: One way severity ('severe ADHD' -> ASD)

effect <- 0.05

l.value.3a <- l.value.0

p.value.3a <- l.value.3a
p.value.3a[ ,1 ] <- 1*( l.value.3a[,1] >= qnorm( 1-k[1] ) )
p.value.3a[ ,2 ] <- 1*( l.value.3a[,2] >= qnorm( 1-k[2] ) )

#Add extreme multiformity
p.value.3a[ p.value.3a[ ,2 ] == 0, 2 ] <- 1*( l.value.3a[ p.value.3a[ ,2 ] == 0, 1 ] >= qnorm( 1-effect*k[1] ) )

plot.3a <- plot.profiles.como.bnw( g.value.0, p.value.3a, "ASD_Model3a_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 3b: One way severity ('severe ASD' -> ADHD)

effect <- 0.1

l.value.3b <- l.value.0

p.value.3b <- l.value.3b
p.value.3b[ ,1 ] <- 1*( l.value.3b[,1] >= qnorm( 1-k[1] ) )
p.value.3b[ ,2 ] <- 1*( l.value.3b[,2] >= qnorm( 1-k[2] ) )

#Add extreme multiformity
p.value.3b[ p.value.3b[ ,1 ] == 0, 1 ] <- 1*( l.value.3b[ p.value.3b[ ,1 ] == 0, 2 ] >= qnorm( 1-effect*k[2] ) )

plot.3b <- plot.profiles.como.bnw( g.value.0, p.value.3b, "ASD_Model3b_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 3c: Two way severity ('severe ADHD' -> ASD & 'severe ASD' -> ADHD)
effect1 <- 0.05
effect2 <- 0.1

l.value.3c <- l.value.0

p.value.3c <- l.value.3c
p.value.3c[ ,1 ] <- 1*( l.value.3c[,1] >= qnorm( 1-k[1] ) )
p.value.3c[ ,2 ] <- 1*( l.value.3c[,2] >= qnorm( 1-k[2] ) )

#Add extreme multiformity
p.value.3c[ p.value.3c[ ,2 ] == 0, 2 ] <- 1*( l.value.3c[ p.value.3c[ ,2 ] == 0, 1 ] >= qnorm( 1-effect1*k[1] ) )
p.value.3c[ p.value.3c[ ,1 ] == 0, 1 ] <- 1*( l.value.3c[ p.value.3c[ ,1 ] == 0, 2 ] >= qnorm( 1-effect2*k[2] ) )

plot.3c <- plot.profiles.como.bnw( g.value.0, p.value.3c, "ASD_Model3c_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 4a: Biological heirarchy ( ADHD > ASD ), with random comorbidity
p.asd <- 0.032

l.value.4a <- l.value.0

p.value.4a <- l.value.4a
p.value.4a[ ,1 ] <- 1*( l.value.4a[,1] >= qnorm( 1-k[1] ) )
p.value.4a[ ,2 ] <- 1*( l.value.4a[,2] >= qnorm( 1-k[2] ) )

# Add biological hierarchy
p.value.4a[ p.value.4a[ ,1 ] == 1 & p.value.4a[ ,2 ] == 1, 2 ] <- 0
p.value.4a[ sample( which( p.value.4a[ ,1 ] == 1 ), p.asd*sum( p.value.4a[ ,1 ] == 1 ) ), 2 ] <- 1

plot.4a <- plot.profiles.como.bnw( g.value.0, p.value.4a, "ASD_Model4a_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 4b: Biological heirarchy ( ASD > ADHD ), with random comorbidity
p.adhd <- 0.166

l.value.4b <- l.value.0

p.value.4b <- l.value.4b
p.value.4b[ ,1 ] <- 1*( l.value.4b[,1] >= qnorm( 1-k[1] ) )
p.value.4b[ ,2 ] <- 1*( l.value.4b[,2] >= qnorm( 1-k[2] ) )

# Add biological hierarchy
p.value.4b[ p.value.4b[ ,1 ] == 1 & p.value.4b[ ,2 ] == 1, 1 ] <- 0
p.value.4b[ sample( which( p.value.4b[ ,2 ] == 1 ), p.adhd*sum( p.value.4b[ ,2 ] == 1 ) ), 1 ] <- 1

plot.4b <- plot.profiles.como.bnw( g.value.0, p.value.4b, "ASD_Model4b_snp.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )


