
library(plyr)
library(ggplot2)
library(tools)
source("diversity-func.R")

usage = "Usage: $1: PLINK frqx file, $2: position file (snpid, pos) \n"

args = commandArgs(TRUE)

# Modify window size as necessary
win_size = 1e5
win_size_str = "100kb"

MAKE_PLOTS = TRUE

## Functions
read.frqx = function(filename){
    fx = read.table(filename, header=FALSE, skip=1, stringsAsFactors = FALSE)
    names(fx)=c("CHR", "SNP", "A1", "A2", "HOM_A1", "HET", "HOM_A2", "HAP_A1", "HAP_A2", "MISS")
    fx
}

frqx_add_basic_stats <- function(fx){
    if(nrow(fx)==0) stop("No data\n")
    # Number of samples
    N = with(fx[1,,drop=F], HOM_A1 + HET + HOM_A2 + MISS)
    attr(fx, "N") = N
    # Frequency of the first allele
    fA1 = with(fx, (2*HOM_A1 + HET) / (2*(N-MISS)))
    fx$frq1 = fA1
    fmiss = with(fx, MISS/N)
    fx$FMISS = fmiss
    
    fhet = with(fx, HET/ ( N - MISS) )
    fx$fhet = fhet
    fx
}


##### Main #####

if(length(args)<2) stop(usage)
frqx_file = args[[1]]
pos_file = args[[2]]

#CHR_NUM = args[[3]]

cat("Loading frequencies...\n")
fx = read.frqx( frqx_file )
fv = frqx_add_basic_stats(fx)
cat("done\n")
N_samples = sum( fx[1,5:10] )


cat("Loading positions...")
pos_df = read.table( pos_file , h=F)
cat("done\n")
names(pos_df) = c("id", "pos")  

for(chr_num in 1:12){
    cat("Chr", chr_num , "\n")
    wh = which(fv$CHR == chr_num )
    if(length(wh)==0) next
    fv_chr  = merge( fv[ fv$CHR == chr_num ,], pos_df, by.x = "SNP", by.y = "id", sort=F)
    
    cat("Calculating heterozygosity per position..\n")
    fv_chr$nucdiv = per_base_nucleotide_diversity(fv_chr$frq1, 2*(N_samples - fv_chr$MISS))
    
    
    print(head(fv_chr))
    nd_ans = one_chr_window_nuc_div( window_bp = win_size, 
                                     pos = fv_chr$pos, 
                                     ndiv = fv_chr$nucdiv,
                                     nsamples = N_samples)  
    cat("\tSaving results...\n")
    table_outfile = paste0( basename(file_path_sans_ext(frqx_file)), "-chr", chr_num, ".nucdiv.txt")
    write.table(nd_ans, table_outfile, quote=F, row.names=F)
    
    if(MAKE_PLOTS){
        plot_outfile = paste0( basename(file_path_sans_ext(frqx_file)), "-chr", chr_num, ".nucdiv.pdf")
        pdf(file=plot_outfile, width=12, height = 4)
        main_title = paste0("Nucleotide diversity (theta and pi) per ", win_size_str, " window, chr ", chr_num, ", using all SNPs")
        plot.sub = paste0("Chr ", chr_num)
        ymax = max( c(nd_ans$avg_pi, nd_ans$theta))
        plot_win_data(nd_ans, what=c("theta", "pi"), 
                      main=main_title, 
                      ylim=c(0,ymax), sub=plot.sub)
        dev.off()
    }

    
}






