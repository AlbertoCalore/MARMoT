#' MARMoT balancing method
#'
#' @description
#' Matching on poset-based average rank for multiple treatments (MARMoT).
#'
#' @details
#' There are many scenarios where classic propensity score techniques are not
#' applicable (e.g. there are many treatments). In a multiple-treatment
#' framework, MARMoT is a method to balance the distribution of covariates among
#' several treatment groups. MARMoT introduces a method for achieving balance
#' among treatment groups by utilizing partially ordered set (poset) theory.
#' This approach focuses on equalizing individual characteristics without
#' relying on propensity score techniques or a dependent variable. Unlike
#' propensity score methods, poset theory doesn't require assumptions about
#' model specifications for treatment allocation. Each subject is represented by
#' a profile of their characteristics, and an average rank approximation is
#' associated with each profile. This value represents the significance of
#' individual characteristics for treatment allocation and can be normalized for
#' better interpretability.
#'
#' @param data
#' A dataframe or equivalent.
#' @param confounders
#' A vector containing the column names of the confounders to balance by.
#' @param treatment
#' A string indicating the column name of the treatment variable.
#' @param reference
#' The statistic used to determine the reference frequencies in the balancing
#' process. Default is median.
#' @param n.cores
#' Number of cores to be used (Linux and Mac systems only!); if a
#' number grater than 1 is specified the function will use a parallelized
#' version of the deloof approximation. Default set to 1.
#'
#' @return
#' A list of objects, also containing the balanced dataset with the same
#' structure of the input dataset.
#' @export
#'
#' @references
#' Silan, Margherita, Arpino, Bruno and Boccuzzo, Giovanna (2021), "Evaluating
#' inverse propensity score weighting in the presence of many treatments. An
#' application to the estimation of the neighbourhood efect.", Journal of
#' Statistical Computation and Simulation, 91(4), 836–859.
#' https://doi.org/10.1080/00949655.2020.1832092
#'
#' Silan, Margherita, Arpino, Bruno and Boccuzzo, Giovanna (2021), "Matching
#' on poset-based average rank for multiple treatments to compare many
#' unbalanced groups.", Statistics in Medicine, 40(28), 6443–6458.
#' https://doi.org/10.1002/sim.9192
#'
#' @examples
#' MARMoT(data = MARMoT_data, confounders = c("race", "age"), treatment = "hospital", n.cores = 1)
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
