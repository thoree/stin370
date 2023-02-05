#' Simulates data for male fathering children with multiple females
#' 
#' Simulates SNP data for father and children
#' 
#' @param nMark Integer. No of markers.
#' @param nCh Integer. No of children.
#' @param plotte Logical. Plots if TRUE.
#' @param write Logical. Writes to file if TRUE.
#' @param seed Integer. Required.
#' 
#' @return Simulated pedigre
#' 
#' @examples 
#' halfSibPedSim(seed = 1729)
#' 
#' @export

halfSibPedSim = function(nMark = 5, nCh = 2, plotte = F, 
                      write = F, seed = NULL){
  if(is.null(seed))
    stop("seed must be specified")
  set.seed(seed)
  
  # Define pedigree with nCh paternal half sibs
  x = nuclearPed(father = "F1", mother = "M1", child ="C1")
  mothers = paste0("M", 1:nCh)
  children = paste0("C", 1:nCh)
  if(nCh > 1){
    for(i in 1:(nCh-1))
      x = addChildren(x, father = "F1", sex = sample(1:2,1), 
                      ids = paste0("C", i+1), verbose = FALSE)
    x = relabel(x, mothers,c("M1", paste0("NN_", 1:(nCh-1))))
  }
  # Define SNP markers
  snps = data.frame(
    CHROM  = 1,
    MARKER = paste0("M", 1:nMark),
    MB     = 1:nMark,
    A1     = rep(1, nMark),
    A2     = rep(2, nMark),
    FREQ1  = rep(0.5, nMark))
  z = setSNPs(x, snpData = snps)
  
  # Simulate
  s = suppressWarnings(
    profileSim(z, simplify1 = T, verbose = F, ids = c("F1", children)))
  # Plot with typed hatched
  if(plotte) plot(s, arrow = T, hatched = typedMembers)
  if(write)
    writePed(s, paste0("simDataSeed", seed))
  s
}
