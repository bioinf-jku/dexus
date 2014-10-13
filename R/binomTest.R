# binomtest by Simon Anders as implemented in the Bioconductor package
# DESEQ.
nbinomTestForMatrices <- function( countsA, countsB, sizeFactorsA, sizeFactorsB, 
		dispsA, dispsB )
{
	
	kAs <- rowSums( cbind(countsA) )
	kBs <- rowSums( cbind(countsB) )
	
	mus <- rowMeans( cbind(      
					t( t( countsA ) / sizeFactorsA ),
					t( t( countsB ) / sizeFactorsB ) ) )      
	
	fullVarsA <- pmax( mus * sum( sizeFactorsA ) + dispsA * mus^2 * sum(sizeFactorsA^2), 
			mus * sum( sizeFactorsA ) * (1+1e-8) )
	fullVarsB <- pmax( mus * sum( sizeFactorsB ) + dispsB * mus^2 * sum(sizeFactorsB^2), 
			mus * sum( sizeFactorsB ) * (1+1e-8) )
	
	sumDispsA <- ( fullVarsA - mus * sum( sizeFactorsA ) ) / ( mus * sum( sizeFactorsA ) )^2
	sumDispsB <- ( fullVarsB - mus * sum( sizeFactorsB ) ) / ( mus * sum( sizeFactorsB ) )^2
	
	sapply( 1:length(kAs), function(i) {
				
				if( kAs[i] == 0 & kBs[i] == 0 )
					return( NA )
				
				# probability of all possible counts sums with the same total count:
				ks <- 0 : ( kAs[i] + kBs[i] )
				ps <- dnbinom(                   ks, mu = mus[i] * sum( sizeFactorsA ), size = 1/sumDispsA[i] ) * 
						dnbinom( kAs[i] + kBs[i] - ks, mu = mus[i] * sum( sizeFactorsB ), size = 1/sumDispsB[i] )
				
				# probability of observed count sums:
				pobs <- dnbinom( kAs[i], mu = mus[i] * sum( sizeFactorsA ), size = 1/sumDispsA[i] ) * 
						dnbinom( kBs[i], mu = mus[i] * sum( sizeFactorsB ), size = 1/sumDispsB[i] )
				
				stopifnot( pobs == ps[ kAs[i]+1 ] )
				if( kAs[i] * sum( sizeFactorsB ) < kBs[i] * sum( sizeFactorsA ) )
					numer <- ps[ 1 : (kAs[i]+1) ]
				else 
					numer <- ps[ (kAs[i]+1) : length(ps) ]
				min( 1, 2 * sum(numer) / sum(ps) )
			} )
}

# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>
# calls the C re-implementation of nbinomTestForMatrices
nbinomTestForMatricesC <- function(countsA, countsB, sizeFactorsA, sizeFactorsB, dispsA, dispsB) {
	kAs <- rowSums(cbind(countsA))
	kBs <- rowSums(cbind(countsB))
	
	mus <- rowMeans(cbind(
					sweep(countsA, 2, sizeFactorsA, "/"),
					sweep(countsB, 2, sizeFactorsB, "/")))
	
	muA = mus*sum(sizeFactorsA)
	muB = mus*sum(sizeFactorsB)
	
	fullVarsA <- pmax(muA + dispsA * mus^2 * sum(sizeFactorsA^2), muA * (1+1e-8) )
	fullVarsB <- pmax(muB + dispsB * mus^2 * sum(sizeFactorsB^2), muB * (1+1e-8) )
	sumDispsA <- (fullVarsA - muA) / (muA)^2
	sumDispsB <- (fullVarsB - muB) / (muB)^2
	
	pval = .Call("dexus_pval_calculation", 
			as.integer(kAs), as.integer(kBs),
			muA, muB,
			1 / sumDispsA, 1 / sumDispsB, 
			sum(sizeFactorsA), sum(sizeFactorsB))
	return (pval)
}