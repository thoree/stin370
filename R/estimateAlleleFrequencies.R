#' Estimates allele frequencies
#' 
#' A wrapper for genetics:::summary.genotype.
#' 
#' @param x a ped object or a list of such
#' @param markers If `NULL`, results are given for all markers

#' @return A list, eachelement containing allele frequencies
#' 
#' @examples 
#' db = simulateUnrelated(nInd = 100, nMark = 2, seed = 1729)
#' estimateAlleleFrequencies(db)
#' estimateAlleleFrequencies(db, markers = 1)
#' @export


estimateAlleleFrequencies = function(x, markers = NULL){
  g = getGenotypes(x)
  nMark = dim(g)[2]
  if(nMark == 0)
    stop("No marker data found")
  if(is.null(markers))
    loopOver = 1:nMark
  else
    loopOver = markers
  if(!all(loopOver %in% 1:nMark))
    stop("Result requested for non-existing markers")
  res = list()
  j = 0
  for(i in loopOver){
    g1 = genotype(g[,i])
    j = j+1
    res[[j]]= genetics:::summary.genotype(g1)$allele.freq
  }
  names(res) = paste("Marker", loopOver)
  res
}
