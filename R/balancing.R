# MARMoT balancing --------------------------------------------------------


# -------------------------------------------------------------------------

bal.table = function(data, treatment){
  tabella = table(data[, "AR"], data[, treatment])
  return(tabella)
}

# -------------------------------------------------------------------------

bal.caliper = function(data){
  caliper = stats::sd(data[, "AR"])/4
  return(caliper)
}

# -------------------------------------------------------------------------

bal.median = function(tabella){
  median = ifelse(ceiling(apply(tabella, 1, median)) == 0, 1,
                  ceiling(apply(tabella, 1, median)))
  return(median)
}

# -------------------------------------------------------------------------

#alternative statistics

# -------------------------------------------------------------------------

bal.ref.selector = function(tabella, reference){
  if(reference == "median"){ref = bal.median(tabella)}
  return(ref)
}

# -------------------------------------------------------------------------

bal.sort.AR = function(AR){
  AR = sort(unique(AR))
  return(AR)
}

# -------------------------------------------------------------------------

bal.treatment.data = function(data, treatment, trt){
  data_trt = data[data[, treatment] == trt, ]
  return(data_trt)
}

# -------------------------------------------------------------------------

bal.zero.match = function(tabella, trt, AR){
  zeroes_match = AR[tabella[, trt] == 0]
  return(zeroes_match)
}

# -------------------------------------------------------------------------

bal.single.ref= function(ref, single_x_match){
  single_ref = as.numeric(ref[names(ref) == single_x_match])
  return(single_ref)
}

# -------------------------------------------------------------------------

bal.AR.ceiling = function(single_zero_match, caliper){
  top_AR = single_zero_match + caliper
  return(top_AR)
}

# -------------------------------------------------------------------------

bal.AR.floor = function(single_zero_match, caliper){
  bottom_AR = single_zero_match - caliper
  return(bottom_AR)
}

# -------------------------------------------------------------------------

bal.closest.elegible = function(elegible, single_zero_match){
  AR_elegible = unique(elegible[, "AR"])
  closest_elegible_AR = AR_elegible[which.min(abs(AR_elegible -
                                                 single_zero_match))]
  closest_elegible = elegible[elegible[, "AR"] == closest_elegible_AR, ]
  return(closest_elegible)
}

# -------------------------------------------------------------------------

bal.elegible.substitutes = function(data_trt, top_AR, bottom_AR){
  elegible = data_trt[data_trt[, "AR"] <= top_AR &
                        data_trt[, "AR"] >= bottom_AR,]
  return(elegible)
}

# -------------------------------------------------------------------------

bal.elegible.data = function(ref, single_zero_match, caliper, data_trt){
  top_AR = bal.AR.ceiling(single_zero_match, caliper)
  bottom_AR = bal.AR.floor(single_zero_match, caliper)
  elegible = bal.elegible.substitutes(data_trt, top_AR, bottom_AR)
  return(elegible)
}

# -------------------------------------------------------------------------

bal.drop.in.table = function(tabella, single_zero_match){
  tabella = tabella[rownames(tabella) != single_zero_match, ]
  return(tabella)
}

# -------------------------------------------------------------------------

bal.drop.in.ref = function(ref, single_zero_match){
  ref = ref[names(ref) != single_zero_match]
  return(ref)
}

# -------------------------------------------------------------------------

bal.drop.in.AR = function(AR, single_zero_match){
  AR = AR[AR != single_zero_match]
  return(AR)
}

# -------------------------------------------------------------------------

bal.exact.match = function(AR, tabella, trt, ref){
  exact_match = AR[tabella[, trt] == ref]
  return(exact_match)
}

# -------------------------------------------------------------------------

bal.exact.data = function(data_trt, exact_match){
  exact = data_trt[data_trt[, "AR"] %in% exact_match,]
  return(exact)
}

# -------------------------------------------------------------------------

bal.exact = function(data_trt, AR, tabella, trt, ref, no_match){
  exact_match = bal.exact.match(AR, tabella, trt, ref)
  exact_match = bal.rm.no.match(exact_match, no_match)
  exact = bal.exact.data(data_trt, exact_match)
  return(exact)
}

# -------------------------------------------------------------------------

bal.inexact.match = function(AR, tabella, trt, ref){
  inexact_match = AR[tabella[, trt] != ref & tabella[, trt] != 0]
  return(inexact_match)
}

# -------------------------------------------------------------------------

bal.single.inexact.data = function(data_trt, single_inexact_match){
  single_inexact_data = data_trt[data_trt[, "AR"] == single_inexact_match,]
  return(single_inexact_data)
}

# -------------------------------------------------------------------------

bal.single.inexact.sampling = function(single_inexact_data, single_ref){
  single_inexact = single_inexact_data[sample(nrow(single_inexact_data),
                                              single_ref, replace = T), ]
  return(single_inexact)
}

# -------------------------------------------------------------------------

bal.inexact = function(data, data_trt, AR, tabella, trt, ref, no_match){
  inexact_match = bal.inexact.match(AR, tabella, trt, ref)
  inexact = data[0,]

  for (single_inexact_match in inexact_match) {

    if(single_inexact_match %in% no_match){next}
    single_ref = bal.single.ref(ref, single_inexact_match)
    single_inexact_data = bal.single.inexact.data(data_trt,
                                                  single_inexact_match)
    single_inexact = bal.single.inexact.sampling(single_inexact_data,
                                                 single_ref)
    inexact = rbind(inexact, single_inexact)
  }
  return(inexact)
}

# -------------------------------------------------------------------------

bal.no.elegible = function(single_zero_match, elegible, no_match){
  if(nrow(elegible) == 0){
    no_match = c(no_match, single_zero_match)
  }
  return(no_match)
}

# -------------------------------------------------------------------------

bal.rm.no.match = function(exact_match, no_match){
  exact_match = exact_match[!exact_match %in% no_match]
  return(exact_match)
}

# -------------------------------------------------------------------------

bal.zeroes = function(data, data_trt, AR, tabella, trt, ref, caliper, no_match){
  zeroes_match = bal.zero.match(tabella, trt, AR)
  zeroes = data[0, ]

  for(single_zero_match in zeroes_match){
    if(single_zero_match %in% no_match){next}
    single_ref = bal.single.ref(ref, single_zero_match)
    elegible = bal.elegible.data(ref, single_zero_match, caliper, data_trt)
    elegible = bal.closest.elegible(elegible, single_zero_match)
    single_zeroes = elegible[sample(nrow(elegible),
                                    single_ref, replace = T), ]
    zeroes = rbind(zeroes, single_zeroes)
  }
  return(zeroes)
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


balancing = function(data, treatment, AR, reference){

  data$AR = AR

  tabella = bal.table(data, treatment)

  caliper = bal.caliper(data)

  ref = bal.ref.selector(tabella, reference)

  # Non match removing -----------------------------------------------------

  AR = bal.sort.AR(AR)
  no_match = c()

  for (trt in unique(data[, treatment])) {

    data_trt = bal.treatment.data(data, treatment, trt)

    zeroes_match = bal.zero.match(tabella, trt, AR)

    for(single_zero_match in zeroes_match){

      elegible = bal.elegible.data(ref, single_zero_match, caliper, data_trt)
      no_match = bal.no.elegible(single_zero_match, elegible, no_match)

    }
  }

  no_match = unique(no_match)

  # Balancing -----------------------------------------------------------

  balanced_data = data[0,]

  for (trt in unique(data[, treatment])) {

    data_trt = bal.treatment.data(data, treatment, trt)

    # Zeroes ----------------------------------------------------------------

    zeroes = bal.zeroes(data, data_trt, AR, tabella, trt, ref, caliper, no_match)

    # Exact match -----------------------------------------------------------

    exact = bal.exact(data_trt, AR, tabella, trt, ref, no_match)

    # Inexact match (but non zero) ------------------------------------------

    inexact = bal.inexact(data, data_trt, AR, tabella, trt, ref, no_match)

    # Merge -----------------------------------------------------------------

    balanced_data = rbind(balanced_data, exact, inexact, zeroes)
  }

  return(balanced_data)
}


