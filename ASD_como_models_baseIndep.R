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

he_reg <- function( g, p, k, m1, m2, misType='flip', misdir=0 ) {
	
	# Using the same people in both studies breaks estimators.  
	#	Sampling needs to be done carefully to ensure E( y1y2 | s1,s2,g12 ) can
	#	be factored under some assumptions of independence (my guess), or 
	#	it doesn't like within and across study similarities to be mixed. 
	#	Weissbrod, Flint, Rosset. 2018. AJHG. Supp. Note, eq. 5. 
	# Derived for Gij ~ 0, so when subjects are compared with themselves it breaks 
	#	this assumption.  With m=1000, anything > ~0.8 is a self-self correlations.
	#	With m=250, unrel x < 0.375, within person x > 0.675 seems ~reasonable. 
	#	With m=100, unrel x and within person x overlap - too much variance in (1/m)( X'X )
	#	I remove those pairs by checking x > 0.5		

	if ( misdir %in% c( 1,3 ) ) {
		m1.p <- sample( which( p[,1] == 1 ), floor( m1*sum( p[,1] ) ) ) 
		p[ m1.p, 2 ] <- 1 
	}
	if ( misdir %in% c( 2,3 ) ) { 
		m2.p <- sample( which( p[,2] == 1 ), floor( m2*sum( p[,2] ) ) )
		p[ m2.p, 1 ] <- 1 
	}

	if ( misType == 'flip' ) {
		if ( misdir %in% c( 1,3 ) ) { p[ m1.p, 1 ] <- 0 }	
		if ( misdir %in% c( 2,3 ) ) { p[ m2.p, 2 ] <- 0 }
	}
		
	# Initial sampling
	g.keep.1 <- sort( unique( c( sample( which( p[,1] == 0 ), 10000 ), 
							which( p[,1] == 1 ) ) ) )
	g.keep.2 <- sort( unique( c( sample( which( p[,2] == 0 ), 10000 ), 
							which( p[,2] == 1 ) ) ) )
	g.keep.3 <- sort( unique( c( sample( 
							which( !is.na( p[,3] ) ), 10000 ) ) ) )
								
	n.1 <- length( g.keep.1 )
	n.2 <- length( g.keep.2 )
	n.3 <- length( g.keep.3 )
		
	g.temp.1 <- g[ g.keep.1, ]
	g.temp.2 <- g[ g.keep.2, ]
	g.temp.3 <- g[ g.keep.3, ]	
	p.temp.1 <- p[ g.keep.1, ]	
	p.temp.2 <- p[ g.keep.2, ]
	p.temp.3 <- p[ g.keep.3, ]
		
	pt.1 <- sum( p.temp.1[ ,1 ] )/n.1
	p.temp.1[ ,1 ] <- ( p.temp.1[ ,1 ] - pt.1 ) / sqrt( pt.1*(1-pt.1) )	
	pt.2 <- sum( p.temp.2[ ,2 ] )/n.2
	p.temp.2[ ,2 ] <- ( p.temp.2[ ,2 ] - pt.2 ) / sqrt( pt.2*(1-pt.2) )

	x.temp <- ( 1/dim(g)[2] ) * g.temp.1 %*% t( g.temp.1 )
	x.temp[ upper.tri( x.temp, diag=T ) ] <- NA
	x.temp.l <- c( x.temp[ !is.na( x.temp ) ] )
	
	y.temp <- p.temp.1[,1] %*% t( p.temp.1[,1] )
	y.temp[ upper.tri( y.temp, diag=T ) ] <- NA
	y.temp.l <- c( y.temp[ !is.na( y.temp ) ] )
	
	c.1.1.p <- sqrt( pt.1*(1-pt.1)*pt.1*(1-pt.1) )
	c.1.1.n <- dnorm( qnorm(1-k[1]) )*dnorm( qnorm(1-k[1]) )
	c.1.1.k <- k[1]*( 1-k[1] )*k[1]*( 1-k[1] )
	c.1.1 <- ( c.1.1.p*c.1.1.n ) / ( c.1.1.k )
	y.1.1.h2 <- ( ( c.1.1*x.temp.l ) %*% y.temp.l ) / 
				( ( c.1.1*x.temp.l ) %*% ( c.1.1*x.temp.l ) )						

	x.temp <- ( 1/dim(g)[2] ) * g.temp.2 %*% t( g.temp.2 )
	x.temp[ upper.tri( x.temp, diag=T ) ] <- NA
	x.temp.l <- c( x.temp[ !is.na( x.temp ) ] )

	y.temp <- p.temp.2[,2] %*% t( p.temp.2[,2] )
	y.temp[ upper.tri( y.temp, diag=T ) ] <- NA
	y.temp.l <- c( y.temp[ !is.na( y.temp ) ] )
	
	c.2.2.p <- sqrt( pt.2*(1-pt.2)*pt.2*(1-pt.2) )
	c.2.2.n <- dnorm( qnorm(1-k[2]) )*dnorm( qnorm(1-k[2]) )
	c.2.2.k <- k[2]*( 1-k[2] )*k[2]*( 1-k[2] )
	c.2.2 <- ( c.2.2.p*c.2.2.n ) / ( c.2.2.k )
	y.2.2.h2 <- ( ( c.2.2*x.temp.l ) %*% y.temp.l ) / 
				( ( c.2.2*x.temp.l ) %*% ( c.2.2*x.temp.l ) )

	x.temp <- ( 1/dim(g)[2] ) * g.temp.3 %*% t( g.temp.3 )
	x.temp[ upper.tri( x.temp, diag=T ) ] <- NA
	x.temp.l <- c( x.temp[ !is.na( x.temp ) ] )

	y.temp <- p.temp.3[,3] %*% t( p.temp.3[,3] )
	y.temp[ upper.tri( y.temp, diag=T ) ] <- NA
	y.temp.l <- c( y.temp[ !is.na( y.temp ) ] )

	# Biased if studies with ascertainment on secondary trait
	y.3.3.h2 <- ( ( x.temp.l ) %*% y.temp.l ) / 
				( ( x.temp.l ) %*% ( x.temp.l ) )

	x.temp <- ( 1/dim(g)[2] ) * g.temp.1 %*% t( g.temp.2 )
	x.temp[ upper.tri( x.temp, diag=T ) ] <- NA
	x.temp.l <- c( x.temp[ !is.na( x.temp ) ] )
	rem.x.temp.l <- which( x.temp.l > 0.5 )

	y.temp <- p.temp.1[,1] %*% t( p.temp.2[,2] )
	y.temp[ upper.tri( y.temp, diag=T ) ] <- NA
	y.temp.l <- c( y.temp[ !is.na( y.temp ) ] )

	c.1.2.p <- sqrt( pt.1*(1-pt.1)*pt.2*(1-pt.2) )
	c.1.2.n <- dnorm( qnorm(1-k[1]) )*dnorm( qnorm(1-k[2]) )
	c.1.2.k <- k[1]*( 1-k[1] )*k[2]*( 1-k[2] )
	c.1.2 <- ( c.1.2.p*c.1.2.n ) / ( c.1.2.k )
	y.1.2.cov <- ( ( c.1.2*x.temp.l[ -rem.x.temp.l ] ) %*% y.temp.l[ -rem.x.temp.l ] ) / 
				( ( c.1.2*x.temp.l[ -rem.x.temp.l ] ) %*% ( c.1.2*x.temp.l[ -rem.x.temp.l ] ) )
	y.1.2.rg <- y.1.2.cov / sqrt( y.1.1.h2*y.2.2.h2 )

	x.temp <- ( 1/dim(g)[2] ) * g.temp.1 %*% t( g.temp.3 )
	x.temp[ upper.tri( x.temp, diag=F ) ] <- NA
	x.temp.l <- c( x.temp[ !is.na( x.temp ) ] )
	rem.x.temp.l <- which( x.temp.l > 0.5 )

	y.temp <- p.temp.1[,1] %*% t( p.temp.3[,3] )
	y.temp[ upper.tri( y.temp, diag=F ) ] <- NA
	y.temp.l <- c( y.temp[ !is.na( y.temp ) ] )

	c.1.3.p <- sqrt( pt.1*(1-pt.1) )
	c.1.3.n <- dnorm( qnorm(1-k[1]) )
	c.1.3.k <- k[1]*( 1-k[1] )
	c.1.3 <- ( c.1.3.p*c.1.3.n ) / ( c.1.3.k )
	y.1.3.cov <- ( ( c.1.3*x.temp.l[ -rem.x.temp.l ] ) %*% y.temp.l[ -rem.x.temp.l ] ) / 
				( ( c.1.3*x.temp.l[ -rem.x.temp.l ] ) %*% ( c.1.3*x.temp.l[ -rem.x.temp.l ] ) )
	y.1.3.rg <- y.1.3.cov / sqrt( y.1.1.h2*y.3.3.h2 )

	x.temp <- ( 1/dim(g)[2] ) * g.temp.2 %*% t( g.temp.3 )
	x.temp[ upper.tri( x.temp, diag=F ) ] <- NA
	x.temp.l <- c( x.temp[ !is.na( x.temp ) ] )
	rem.x.temp.l <- which( x.temp.l > 0.5 )
		
	y.temp <- p.temp.2[,2] %*% t( p.temp.3[,3] )
	y.temp[ upper.tri( y.temp, diag=F ) ] <- NA
	y.temp.l <- c( y.temp[ !is.na( y.temp ) ] )

	c.2.3.p <- sqrt( pt.2*(1-pt.2) )
	c.2.3.n <- dnorm( qnorm(1-k[2]) )
	c.2.3.k <- k[2]*( 1-k[2] )
	c.2.3 <- ( c.2.3.p*c.2.3.n ) / ( c.2.3.k )
	y.2.3.cov <- ( ( c.2.3*x.temp.l[ -rem.x.temp.l ] ) %*% y.temp.l[ -rem.x.temp.l ] ) / 
				( ( c.2.3*x.temp.l[ -rem.x.temp.l ] ) %*% ( c.2.3*x.temp.l[ -rem.x.temp.l ] ) )
	y.2.3.rg <- y.2.3.cov / sqrt( y.2.2.h2*y.3.3.h2 )
	
	output <- list( y.1.1.h2=y.1.1.h2, y.2.2.h2=y.2.2.h2, y.3.3.h2=y.3.3.h2, 
						y.1.2.rg=y.1.2.rg, y.1.3.rg=y.1.3.rg, y.2.3.rg=y.2.3.rg,
						m1=m1, m2=m2 )
	
	return( output )
					
}

## Basic parameters 

traits <- c( "ADHD", "ASD", "EA" )
legend.txt <- c( "ADHD-, ASD-", "ADHD+, ASD+", "ADHD+, ASD-", "ADHD-, ASD+" )
n=100000
m=250
k=c( 0.05, 0.01, NA )
h2=c( 0.7, 0.8, 0.4 )
rg=matrix( c( 1, 0, -0.5, 0, 1, 0.2, -0.5, 0.2, 1), 3, 3 )
re=diag( length(h2) )
onset <- c( 5, 0, NA )		#ADHD uniform: 5 to 20, ASD: uniform 0 to 15
offset <- c( 20, 15, NA)
t <- length( h2 )

## Model 0: Uncorrelated liability and Comorbidity

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
plot.0 <- plot.profiles.como.bnw( g.value.0, p.value.0,"ASD_v2_Model0.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Learn reasonable misdiagnosis rates (to create an rg of 0.3)

# # output <- NULL
# for ( i in 0:15 ) {
	# print( i )
	# m1 <- i/100
	# m2 <- i/100
	# output <- rbind( output, he_reg( g, p.value.0, k, m1, m2, misType='flip', misdir=3 ) )
	# print( output )
# }
# save( output, file="misDx_famh2_rg0_flipMisDx_01to10.RData")

# output <- NULL
# for ( i in 0:15 ) {
	# print( i )
	# m1 <- i/100
	# m2 <- i/100
	# output <- rbind( output, he_reg( g, p.value.0, k, m1, m2, misType='add', misdir=3 ) )
	# print( output )
# }
# save( output, file="misDx_famh2_rg0_addMisDx_01to11.RData")

# output <- NULL
# for ( i in 0:10 ) {
	# print( i )
	# m1 <- i/50
	# m2 <- i/50
	# output <- rbind( output, he_reg( g, p.value.0, k, m1, m2, misType='add', misdir=1 ) )
	# print( output )
# }
# save( output, file="misDx_famh2_rg0_addMisDx_01to11_adhdOnly.RData")

# output <- NULL
# for ( i in 0:10 ) {
	# print( i )
	# m1 <- i/50
	# m2 <- 0
	# output <- rbind( output, he_reg( g, p.value.0, k, m1, m2, misType='flip', misdir=1 ) )
	# print( output )
# }
# save( output, file="misDx_famh2_rg0_flipMisDx_01to10_adhdOnly.RData")

# output <- NULL
# for ( i in 0:10 ) {
	# print( i )
	# m1 <- 0
	# m2 <- i/10
	# output <- rbind( output, he_reg( g, p.value.0, k, m1, m2, misType='flip', misdir=2 ) )
	# print( output )
# }
# save( output, file="misDx_famh2_rg0_flipMisDx_01to10_asdOnly.RData")

 # output <- NULL
 # for ( i in 19:0 ) {
	# print( i )
	# m1 <- 0
	# m2 <- i/20
	# output <- rbind( output, he_reg( g, p.value.0, k, m1, m2, misType='add', misdir=2 ) )
	# print( output )
 # }
 # save( output, file="misDx_famh2_rg0_addMisDx_01to11_asdOnly.RData")

# Model 1a (flip): Misdiagnosis one-way, flipped (ADHD -> ASD) up to rG of 0.3.  10% rate does it.

misDx.r <- 0.10

p.value.1a <- l.value.0
p.value.1a[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1a[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1a[,1] == 1 ), 
					floor( sum( p.value.1a[,1] == 1 )*misDx.r ) )

p.value.1a[ misDx, 1 ] <- 0
p.value.1a[ misDx, 2 ] <- 1

plot.1a <- plot.profiles.como.bnw( g.value.0, p.value.1a, "ASD_v2_Model1a_flip.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1b (flip): Model 1a: Misdiagnosis one-way, flipped (ASD -> ADHD) up to rG of 0.3.  80% rate probably still not enough.

misDx.r <- 0.8

p.value.1b <- l.value.0
p.value.1b[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1b[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1b[,2] == 1 ), 
					floor( sum( p.value.1b[,2] == 1 )*misDx.r ) )

p.value.1b[ misDx, 1 ] <- 1
p.value.1b[ misDx, 2 ] <- 0

plot.1b <- plot.profiles.como.bnw( g.value.0, p.value.1b, "ASD_v2_Model1b_flip.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1c (flip): Misdiagnosis two-way, flipped (ASD -> ADHD & ADHD -> ASD) up to rG of 0.3.  8% rate enough.

misDx.r <- 0.08

p.value.1c <- l.value.0
p.value.1c[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1c[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx1 <- sample( which( p.value.1c[,1] == 1 ), 
					floor( sum( p.value.1c[,1] == 1 )*misDx.r ) )
misDx2 <- sample( which( p.value.1c[,2] == 1 ), 
					floor( sum( p.value.1c[,2] == 1 )*misDx.r ) )

p.value.1c[ misDx1, 1 ] <- 0
p.value.1c[ misDx1, 2 ] <- 1
p.value.1c[ misDx2, 1 ] <- 1
p.value.1c[ misDx2, 2 ] <- 0

plot.1c <- plot.profiles.como.bnw( g.value.0, p.value.1c, "ASD_v2_Model1c_flip.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1a (add): Misdiagnosis one-way, add COMO (ADHD -> ASD) up to rG of 0.3.  10% rate does it.

misDx.r <- 0.10

p.value.1a <- l.value.0
p.value.1a[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1a[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1a[,1] == 1 ), 
					floor( sum( p.value.1a[,1] == 1 )*misDx.r ) )

p.value.1a[ misDx, 2 ] <- 1

plot.1a <- plot.profiles.como.bnw( g.value.0, p.value.1a, "ASD_v2_Model1a_add.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1b (add): Model 1a: Misdiagnosis one-way, add COMO (ASD -> ADHD) up to rG of 0.3.  80% rate probably still not enough.

misDx.r <- 0.8

p.value.1b <- l.value.0
p.value.1b[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1b[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx <- sample( which( p.value.1b[,2] == 1 ), 
					floor( sum( p.value.1b[,2] == 1 )*misDx.r ) )

p.value.1b[ misDx, 1 ] <- 1

plot.1b <- plot.profiles.como.bnw( g.value.0, p.value.1b, "ASD_v2_Model1b_add.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

# Model 1c (add): Misdiagnosis two-way, add COMO (ASD -> ADHD & ADHD -> ASD) up to rG of 0.3.  8% rate enough.

misDx.r <- 0.08

p.value.1c <- l.value.0
p.value.1c[ ,1 ] <- 1*( l.value.0[,1] >= qnorm( 1-k[1] ) )
p.value.1c[ ,2 ] <- 1*( l.value.0[,2] >= qnorm( 1-k[2] ) )
											
misDx1 <- sample( which( p.value.1c[,1] == 1 ), 
					floor( sum( p.value.1c[,1] == 1 )*misDx.r ) )
misDx2 <- sample( which( p.value.1c[,2] == 1 ), 
					floor( sum( p.value.1c[,2] == 1 )*misDx.r ) )

p.value.1c[ misDx1, 2 ] <- 1
p.value.1c[ misDx2, 1 ] <- 1

plot.1c <- plot.profiles.como.bnw( g.value.0, p.value.1c, "ASD_v2_Model1c_add.bnw.pdf", 
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
											
plot.2a <- plot.profiles.como.bnw( g.value.0, p.value.2a, "ASD_v2_Model2a.bnw.pdf", 
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
											
plot.2b <- plot.profiles.como.bnw( g.value.0, p.value.2b, "ASD_v2_Model2b.bnw.pdf", 
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
											
plot.2c <- plot.profiles.como.bnw( g.value.0, p.value.2c, "ASD_v2_Model2c.bnw.pdf", 
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
											
plot.2d <- plot.profiles.como.bnw( g.value.0, p.value.2d, "ASD_v2_Model2d.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 3a: One way severity ('severe ADHD' -> ASD)

effect <- 0.05

l.value.3a <- l.value.0

p.value.3a <- l.value.3a
p.value.3a[ ,1 ] <- 1*( l.value.3a[,1] >= qnorm( 1-k[1] ) )
p.value.3a[ ,2 ] <- 1*( l.value.3a[,2] >= qnorm( 1-k[2] ) )

#Add extreme multiformity
p.value.3a[ p.value.3a[ ,2 ] == 0, 2 ] <- 1*( l.value.3a[ p.value.3a[ ,2 ] == 0, 1 ] >= qnorm( 1-effect*k[1] ) )

plot.3a <- plot.profiles.como.bnw( g.value.0, p.value.3a, "ASD_v2_Model3a.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 3b: One way severity ('severe ASD' -> ADHD)

effect <- 0.1

l.value.3b <- l.value.0

p.value.3b <- l.value.3b
p.value.3b[ ,1 ] <- 1*( l.value.3b[,1] >= qnorm( 1-k[1] ) )
p.value.3b[ ,2 ] <- 1*( l.value.3b[,2] >= qnorm( 1-k[2] ) )

#Add extreme multiformity
p.value.3b[ p.value.3b[ ,1 ] == 0, 1 ] <- 1*( l.value.3b[ p.value.3b[ ,1 ] == 0, 2 ] >= qnorm( 1-effect*k[2] ) )

plot.3b <- plot.profiles.como.bnw( g.value.0, p.value.3b, "ASD_v2_Model3b.bnw.pdf", 
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

plot.3c <- plot.profiles.como.bnw( g.value.0, p.value.3c, "ASD_v2_Model3c.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 4a: Biological heirarchy ( ADHD > ASD ), with random comorbidity
p.asd <- 0.01

l.value.4a <- l.value.0

p.value.4a <- l.value.4a
p.value.4a[ ,1 ] <- 1*( l.value.4a[,1] >= qnorm( 1-k[1] ) )
p.value.4a[ ,2 ] <- 1*( l.value.4a[,2] >= qnorm( 1-k[2] ) )

# Add biological hierarchy
p.value.4a[ p.value.4a[ ,1 ] == 1 & p.value.4a[ ,2 ] == 1, 2 ] <- 0
p.value.4a[ sample( which( p.value.4a[ ,1 ] == 1 ), p.asd*sum( p.value.4a[ ,1 ] == 1 ) ), 2 ] <- 1

plot.4a <- plot.profiles.como.bnw( g.value.0, p.value.4a, "ASD_v2_Model4a.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )

## Model 4b: Biological heirarchy ( ASD > ADHD ), with random comorbidity
p.adhd <- 0.05

l.value.4b <- l.value.0

p.value.4b <- l.value.4b
p.value.4b[ ,1 ] <- 1*( l.value.4b[,1] >= qnorm( 1-k[1] ) )
p.value.4b[ ,2 ] <- 1*( l.value.4b[,2] >= qnorm( 1-k[2] ) )

# Add biological hierarchy
p.value.4b[ p.value.4b[ ,1 ] == 1 & p.value.4b[ ,2 ] == 1, 1 ] <- 0
p.value.4b[ sample( which( p.value.4b[ ,2 ] == 1 ), p.adhd*sum( p.value.4b[ ,2 ] == 1 ) ), 1 ] <- 1

plot.4b <- plot.profiles.como.bnw( g.value.0, p.value.4b, "ASD_v2_Model4b.bnw.pdf", 
									traits.txt=traits, legend.txt=legend.txt )


