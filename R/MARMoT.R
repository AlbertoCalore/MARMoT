#' MARMoT balancing method
#'
#' @description
#' Matching on poset-based average rank for multiple treatments (MARMoT)
#'
#' @param data
#' a dataframe or equivalent
#' @param confounders
#' a vector containing the column names of the confounders to balance by
#' @param treatment
#' a string indicating the column name of the treatment variable
#' @param reference
#' "median" default, the statistic used to determine the reference frequencies in the balancing process
#' @param n.cores
#' 1 default, number of cores to be used (Linux and Mac systems only!); if a number grater than 1 is specified the function will use a parallelized version of the deloof approximation
#'
#' @return
#' a list of objects, also containing the balanced dataset with the same structure of the input dataset
#' @export
#'
#' @examples
#' out = MARMoT(data = MARMoT_data, confounders = c("race", "age"),
#'       treatment = "hospital", n.cores = 1)
#'
MARMoT = function(data, confounders, treatment, reference = "median", n.cores = 1){


  # Convert factor to ordinal -----------------------------------------------


  comparable_confounders = comparables(data = data, confounders = confounders)


  # ASB pre balancing -------------------------------------------------------


  print("Absolute standardized bias - before balancing")
  ASB_pre = ASB(data = data, confounders = confounders, treatment = treatment)


  # Average rank with deloof approx -----------------------------------------


  print("... Computing average rank ...")
  if(n.cores == 1){
    AR = deloof(comparable_confounders)}
  if(n.cores > 1){
    AR = mcdeloof(comparable_confounders, n.cores)}


  # Normalize AR ------------------------------------------------------------


  AR = normalize.AR(AR)


  # Balancing ---------------------------------------------------------------


  print("... Balancing ...")
  balanced_data = balancing(data = data, treatment = treatment, AR = AR, reference)


  # ASB post ----------------------------------------------------------------


  print("Absolute standardized bias - after balancing")
  ASB_post = ASB(data = balanced_data, confounders = confounders, treatment = treatment)


  # Output ------------------------------------------------------------------


  output = output.maker(balanced_data, ASB_pre, ASB_post)


  return(output)
}
