#' MLE estimate for SNP marker
#' 
#' Estimates allele frequency of first allele in specified marker by maximising the likelihood.
#' 
#' @param x	A ped object.
#' @param mark Index or name of marker
#' 
#' @return Estimate of frequency of first allele.
#' 
#' @examples 
#' x = halfSibPedSim(seed = 1729, nCh = 100, nMark = 2)
#' mleSNP(x, mark = 1)
#' mleSNP(x, mark = NULL)
#' 
#' @export

mleSNP = function(x, mark = 1, alleles = 1:2){

# Help function, log likelihood on format suitable for optimisation  
  loglik = function(p, x, mark, alleles = 1:2){
    afreq = c(p, 1 - p)
    names(afreq) = alleles
    ped = setAfreq(x, marker = mark, afreq = afreq)
    -log(likelihood(ped, mark))
  }
  
  if(is.null(mark)){ # calculate for all markers
    nM = nMarkers(x)
    p.hat = rep(NA, nM)
    names(p.hat) = name(x,markers = 1:nM)
    for (i in 1:nM)
      p.hat[i] =  optimize(loglik, c(0,1), x, i, alleles)$minimum
    } 
  else{
    p.hat = optimize(loglik, c(0,1), x, mark, alleles)$minimum
    names(p.hat) = name(x,markers = mark)
  }
  p.hat
}


