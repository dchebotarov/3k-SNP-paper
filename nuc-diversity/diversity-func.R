library(plyr)
options(digits=4)


#'  Nucleotide diversity (heterozygosity) for a biallelic SNP
#'  calculated by the formula
#'  pi = n/(n-1) * 2pq
#'  
#'  @param frq Vector of SNP frequencies
#'  @param N Vector of numbers of nonmissing calls for SNPs. Same length as frq.
#'  @return Vector of pi values
per_base_nucleotide_diversity = function(frq, N){
    Pi = 2 * frq * (1-frq) * N / (N-1)   
    Pi[ N == 1 ] = NA
    Pi
}


#' Compute theta, a nucleotide diversity estimate from number of polymorphic sites
#' Note: Not used in the paper.
#' @param poly Vector of number of polymorphic SNPs in the windows
#' @param nsites Number (or vector of numbers per window) of positions in a window
#' @param nsamples A single number, number of samples
#' @return A vector of values of theta
theta = function(poly, nsites, nsamples){
    s = poly / nsites;
    a = sum( 1/(1:(nsamples-1)))
    s/a
}


#'  Tajima's D calculation
#' @param Pi mean number of sites that differ in a window
#' @param S  mean number of segregating sites in a window
#' @param n  number of samples
tD = function(Pi, S, n){
    if(length(S) != length(Pi)) stop("pi and s must be same length\n")
    a1 = sum( 1/(1:(n-1)) )
    a2 = sum( 1/(1:(n-1))^2 )
    b1 = (n+1)/(3*(n-1))
    b2 = 2/9*(n^2 + n + 3)/(n*(n-1))
    c1 = b1 - 1/a1
    c2 = b2 - (n+2)/(a1*n)  + a2 / (a1^2)
    e1 = c1 / a1
    e2 = c2 / ( a1^2 + a2 )
    # cat("a1=", a1, "\n a2=", a2, "\nb1=",b1, "\nb2=", b2, "\n c1=", c1, "\n c2=",c2, "\n e1=", e1, 
    #    "\n e2=", e2, "\n")
    d = (Pi - S/a1 ) / ( S * e1 + e2 * S *(S-1))^0.5
    d
}



#' Compute two measures of nucleotide diversity, pi and theta, per window chromosomewide
#' @param window_bp window size in basepairs
#' @param pos list of positions in basepairs
#' @param ndiv Nucleotide diversity (pi) per site
#' @param nsamples Number of samples in the dataset
#' @return A data frame with columns: avg_pi, theta, tajD, pos
one_chr_window_nuc_div = function(window_bp, pos, ndiv, nsamples){
    if(length(window_bp)!=1) stop("window_bp should be a scalar value, length of window in basepairs")
    notNA = is.finite(ndiv)
    pos = pos[notNA]
    ndiv = ndiv[notNA]
    cat("The dataset has", length(pos), " positions after removing those with N-1 missing\n")
    df = cbind.data.frame(pos=pos, ndiv=ndiv)
    win_start = seq(0, max(pos), by=window_bp)
    df$win_ix = findInterval(pos, win_start, all.inside=T)
    # Window data
    cat("Computing pi and theta per window...\n")
    wdata = ddply(df, .(win_ix), summarize, 
                  total_pi = sum(ndiv, na.rm=TRUE) , 
                  num_poly_sites=length(ndiv) )
    wdata$avg_pi = wdata$total_pi / window_bp
    wdata$theta  = theta( wdata$num_poly_sites, window_bp, nsamples )
    wdata$tajD = tD(Pi=wdata$total_pi, S=wdata$num_poly_sites,n=nsamples)
    wdata$pos = as.integer(  win_start[ wdata$win_ix ]  + window_bp / 2 ) # middle positioin per window
    wdata
}


#' Wrapper around previous functions to compute nuc.diversoty from frequency and position data
nuc_diversity_from_freqx_pos = function(window_bp, chr, pos, frq, miss, nsamples ){
    nucdiv = per_base_nucleotide_diversity(frq, 2*(nsamples - miss) ) # N is the num of chromosomes, so mult by 2
    ans = NULL
    for(chr_num in unique(chr)){
        ss = (chr_num == chr)
        #ss[ is.na(nucdiv) ] = FALSE  # ignore
        df = one_chr_window_nuc_div(window_bp, pos[ss], nucdiv[ss], nsamples )
        df = cbind.data.frame(chr=chr_num, df)
        if(is.null(ans)){ 
            ans = df
        } else {
            ans = rbind(ans, df)
        }
    }
    ans
}


plot_win_data = function(wdata, what=c("pi", "theta"), ...){
    
    if( all( c("theta", "pi") %in% what)){
        with(wdata, 
             {
                 m = max( c(theta,pi) )
                 plot( pos/1000, theta, type="l", xlab="Position, kb", 
                       ylab=expression( paste("Nucleotide diversity (", pi, " and ", theta, ") per window" )) 
                       ,...)
                 lines( pos/1000, avg_pi	, col="red", ...)			
                 legend("topright", legend=c(expression("theta"), "pi"), col=c("black", "red"), lty=1)
             }
        )
    } else if( what == "pi"){
        with( wdata, plot( pos/1000, avg_pi, type="l", xlab="Position, kb", 
                           ylab=expression( paste("Mean ", pi, " per window")) 
                           , ...) )        
    } else if(what == "theta"){
        with( wdata, plot( pos/1000, theta, type="l", xlab="Position, kb", 
                           ylab=expression( paste( theta, " per window")) 
                           , ...) )
    } else  stop("No options given. Use what=c('pi', 'theta') argument " )
}



