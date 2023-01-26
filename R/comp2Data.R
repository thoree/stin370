#' Simulates data for Compulsory 2
#' 
#' Simulates data for half sibs
#' 
#' @param nMark Integer. No of markers
#' @param nCh Integer. No of children, at least 2
#' @param plotte Logical. Plots if TRUE
#' @param write Logical. Writes to file if TRUE
#' @param seed Integer. Required
#' 
#' @return Simulated half sib ped.
#' 
#' @examples 
#' comp2Data(seed = 1729)
#' 
#' @export



comp2Data = function(nMark = 5, nCh = 2, plotte = F, 
                     write = F, seed = NULL){
  if (nCh < 2)
    stop("Require at least 2 children (variable nCh)")
  if(is.null(seed))
    stop("seed must be specified")
  set.seed(seed)
  
  # Define pedigree with nCh paternal half sibs
  x = nuclearPed(father = "F1", mother = "M1", child ="C1")
  mothers = paste0("M", 1:nCh)
  children = paste0("C", 1:nCh)
  for(i in 1:(nCh-1))
    x = addChildren(x, father = "F1", sex = sample(1:2,1), 
                    ids = paste0("C", i+1), verbose = FALSE)
  x = relabel(x, mothers,c("M1", paste0("NN_", 1:(nCh-1))))
  
  # Define SNP markers
  snps = data.frame(
    CHROM  = 1,
    MARKER = paste0("M", 1:nMark),
    MB     = 1:nMark,
    A1     = rep(1, nMark),
    A2     = rep(2, nMark),
    FREQ1  = rep(0.5, nMark))
  z = setSNPs(x, snpData = snps)
  
  # Simulate. Introduce one false paternity
  skip = sample(2:nCh,1)
  s = suppressWarnings(
    profileSim(z, simplify1 = T, verbose = F, ids = c("F1", children[-skip])))
  id = paste0("C",skip)
  unr = singleton(paste0("C",skip), sex = getSex(s, id))
  unr = setSNPs(x, snpData = snps)
  unr = suppressWarnings(
    profileSim(unr, simplify1 = T, verbose = F, ids = id))
  s = transferMarkers(unr, s, erase = F, ids = id)
  # Plot with typed hatched
  if(plotte) plot(s, arrow = T, title = "Presumed pedigree", hatched = typedMembers)
  if(write)
    writePed(s, paste0("simDataSeed", seed))
  s
}
