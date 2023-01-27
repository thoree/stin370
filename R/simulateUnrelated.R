#' Simulates marker data for unrelated individuals
#' 
#' A wrapper for `markerSim`.
#' 
#' @param nMark Integer. No of markers
#' @param nInd Integer. No of individuals
#' @param alleles designations
#' @param afreq frequencies
#' @param seed Integer. 
#' 
#' @return Simulated marker data
#' 
#' @examples 
#' simulateUnrelated(seed = 1729)
#' 
#' @export


simulateUnrelated = function(nMark = 5, nInd = 2, alleles = NULL,
                             afreq = NULL, seed = NULL){
  ind = lapply(1:nInd, function(x) singleton(x))
  db = markerSim(ind, N = nMark, alleles = alleles, afreq = afreq, 
                 seed = seed, verbose = F)
  db
}
