suppressPackageStartupMessages({
	library(argparse)
	library(Rmpfr)
	library(data.table)
})
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--chisq_per_nampc_file",type="character")
parser$add_argument("--ks_file",type="character")
parser$add_argument("--outfile",type="character")
args <- parser$parse_args()

# Import
all_res = fread(args$chisq_per_nampc_file, header=FALSE) # T values
all_res = all_res**2 # T-squared
ks = read.table(args$ks_file, header=FALSE)[,1]

# Evaluate multi-NAM-PC associations across SNPs
for(k in ks){
    all_res[[paste0("k",k,"_P")]] = apply(as.matrix(rowSums(all_res[,1:k, drop=FALSE]), ncol=1), 1, pchisq, df = k, lower = F)
}
col_names <- paste0("k",ks,"_P")
all_res[['P']] = apply(as.matrix(all_res[,..col_names, drop=FALSE], ncol=length(ks)), 1, min)
small_vals = which(all_res[,'P']<1e-15) # If p<2.2e-16, need small-values handling
all_res[['P']] = 1-(1-all_res[['P']])**length(ks) # Multiple testing correction across values for k
all_res[['k']]= apply(as.matrix(all_res[,..col_names, drop=FALSE], ncol=length(ks)), 1, which.min)
all_res$k <- as.numeric(all_res$k)
all_res[['k']] = ks[all_res$k]

# Small-values handling with precision
for(k in ks){ # Define sum of Tsq for each value of k considered
    all_res[[paste0('k',k,'_sumTsq')]] = apply(as.matrix(rowSums(all_res[,1:k, drop=FALSE]), ncol=1), 1, sum)
}

p_mht_corrected <-function(sumTsq_vec, ks, prec_bits = 10000){
    # Use Rmpfr to evaluate associations and perform multiple-testing correction
    uncorr_ps=c(exp(mpfr(pchisq(sumTsq_vec[1], df = ks[1], lower = F, log.p=TRUE),precBits = prec_bits)))
    if(length(ks)>1){ 
	for(i in c(2:length(ks))){
	      new_p = exp(mpfr(pchisq(sumTsq_vec[i], df = ks[i], lower = F, log.p=TRUE),precBits = prec_bits))
	      uncorr_ps = c(uncorr_ps, new_p)
	}
    }
    sel_k = as.character(ks[which.min(uncorr_ps)])
    # Multiple testing correction across values for k
    corr_p = formatDec(mpfr(1, prec_bits)-(mpfr(1, prec_bits)-mpfr(min(uncorr_ps), prec_bits))**length(ks), precBits=64)
    return(c(corr_p, sel_k))
}

col_names <- paste0("k",ks,"_sumTsq")
mtx <- t(apply(as.matrix(all_res[small_vals,..col_names, drop=FALSE], ncol=length(ks)), 1, p_mht_corrected, ks))
all_res[small_vals, P := mtx[, 1]]
all_res[small_vals, k := mtx[, 2]]

fwrite(all_res[,c("P", "k")], args$outfile, quote=FALSE, row.names=FALSE, sep = "\t")
