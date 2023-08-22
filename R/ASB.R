#' Absolute standardized bias
#'
#' @description
#' Compute the absolute standardized bias of given confounders and return some
#' useful statistics.
#'
#' @param data
#' A dataframe or equivalent.
#' @param confounders
#' A vector with the column names of the confounders to balance by.
#' @param treatment
#' A string with the column name of the treatment variable.
#'
#' @return
#' A list of objects, containing the ASB matrix and some summary statistics.
#' @export
#'
#' @examples
#' ASB(data = MARMoT_data, confounders = c("race", "age"), treatment = "hospital")
#'
ASB = function(data, confounders, treatment){

  ASB_data = ASB.table.init(data, confounders, treatment)

  ASB_data_single = ASB.confounders.tables(confounders, data, treatment, ASB_data)

  ASB_data = ASB.table(ASB_data_single, ASB_data)

  nms_list = ASB.tab.names(data, confounders)

  ASB_data = ASB.tab.give.rownames(ASB_data, nms_list)

  n_trt = ASB.n.treatment(data, treatment)

  all_pop = ASB.pop.confounders(ASB_data)

  ASB_data = ASB.prop.confounders(ASB_data, n_trt)

  all_pop = ASB.all.pop.prop(all_pop, data)

  all_pop_matrix = ASB.pop.matrix(ASB_data, all_pop)

  numerator = ASB.numerator(ASB_data, all_pop_matrix)

  denominator = ASB.denominator(ASB_data, all_pop_matrix)

  ASB_result = ASB.result.matrix(numerator, denominator)

  ASB_statistics = ASB.statistics(ASB_result)
  print(round(ASB_statistics, digits = 3))

  output = ASB.out(ASB_result, ASB_statistics)

  return(output)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

ASB.table.init = function(data, confounders, treatment){
  ASB_data = table(data[, confounders[1]], data[, treatment])[0,]
  return(ASB_data)
}

# -------------------------------------------------------------------------

ASB.single.table.data = function(c, data, treatment){
  ASB_single_table_data = table(data[, c], data[, treatment])
  return(ASB_single_table_data)
}

# -------------------------------------------------------------------------

ASB.single.table = function(ASB_single_table_data, ASB_data){
  if(nrow(ASB_single_table_data) == 2){
    ASB_data_single = rbind(ASB_data, ASB_single_table_data[1, ])
  } else(ASB_data_single = rbind(ASB_data, ASB_single_table_data))
  return(ASB_data_single)
}

# -------------------------------------------------------------------------

ASB.binding.tables = function(c, data, treatment, ASB_data){
  ASB_single_table_data = ASB.single.table.data(c, data, treatment)
  ASB_data_single = ASB.single.table(ASB_single_table_data, ASB_data)
  return(ASB_data_single)
}

# -------------------------------------------------------------------------

ASB.confounders.tables = function(confounders, data, treatment, ASB_data){
  ASB_data_single = lapply(confounders, ASB.binding.tables, data, treatment, ASB_data)
  return(ASB_data_single)
}


# -------------------------------------------------------------------------

ASB.table = function(ASB_data_single, ASB_data){
  for (i in 1:length(ASB_data_single)) {
    ASB_data = rbind(ASB_data, ASB_data_single[[i]])
  }
  return(ASB_data)
}

# -------------------------------------------------------------------------

ASB.tab.single.name = function(c, data){
  if(is.factor(data[, c]))
    {lvs = levels(data[, c])}else{lvs = unique(data[, c])}
  nms = paste0(c, sep = "_", lvs)
  if(length(nms) == 2){nms = nms[1]}
  return(nms)
}

# -------------------------------------------------------------------------

ASB.tab.names = function(data, confounders){
  nms_list = lapply(confounders, ASB.tab.single.name, data)
  nms_list = unlist(nms_list)
  return(nms_list)
}

# -------------------------------------------------------------------------

ASB.tab.give.rownames = function(ASB_data, nms_list){
  rownames(ASB_data) = nms_list
  return(ASB_data)
}


# -------------------------------------------------------------------------

ASB.n.treatment = function(data, treatment){
  n_trt = as.numeric(table(data[, treatment]))
  return(n_trt)
}

# -------------------------------------------------------------------------

ASB.pop.confounders = function(ASB_data){
  all_pop = rowSums(ASB_data)
  return(all_pop)
}


# -------------------------------------------------------------------------

ASB.prop.confounders = function(ASB_data, n_trt){
  ASB_data = t(t(ASB_data) / n_trt)
  return(ASB_data)
}

# -------------------------------------------------------------------------

ASB.all.pop.prop = function(all_pop, data){
  all_pop = all_pop / nrow(data)
  return(all_pop)
}


# -------------------------------------------------------------------------

ASB.pop.matrix.single.col = function(n, all_pop){
  n = all_pop
  return(all_pop)
}

# -------------------------------------------------------------------------

ASB.pop.matrix = function(ASB_data, all_pop){
  all_pop_matrix = apply(ASB_data, 2, ASB.pop.matrix.single.col, all_pop)
  return(all_pop_matrix)
}

# -------------------------------------------------------------------------

ASB.numerator = function(ASB_data, all_pop_matrix){
  numerator = abs(ASB_data - all_pop_matrix)
  return(numerator)
}

# -------------------------------------------------------------------------

ASB.var.trt = function(ASB_data){
  var_trt = ASB_data*(1-ASB_data)
  return(var_trt)
}

# -------------------------------------------------------------------------

ASB.var.pop = function(all_pop_matrix){
  var_pop = all_pop_matrix*(1 - all_pop_matrix)
  return(var_pop)
}

# -------------------------------------------------------------------------

ASB.denominator = function(ASB_data, all_pop_matrix){
  var_trt = ASB.var.trt(ASB_data)
  var_pop = ASB.var.pop(all_pop_matrix)
  denominator = sqrt(var_pop/2 + var_trt/2)
  return(denominator)
}

# -------------------------------------------------------------------------

ASB.result.matrix = function(numerator, denominator){
  ASB_result = numerator / denominator *100
  return(ASB_result)
}

# -------------------------------------------------------------------------

ASB.quantiles = function(ASB_result){
  quant = stats::quantile(ASB_result, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(quant) = c("Min", "1st quartile", "Median", "3rd quartile", "Max")
  return(quant)
}

# -------------------------------------------------------------------------

ASB.mean = function(ASB_result){
  mn = mean(ASB_result)
  names(mn) = "Mean"
  return(mn)
}

# -------------------------------------------------------------------------

ASB.over5 = function(ASB_result){
  over5 = length(which(ASB_result>5))
  names(over5) = "Over 5%"
  return(over5)
}

# -------------------------------------------------------------------------

ASB.over10 = function(ASB_result){
  over10 = length(which(ASB_result>10))
  names(over10) = "Over 10%"
  return(over10)
}

# -------------------------------------------------------------------------

ASB.statistics = function(ASB_result){
  quant = ASB.quantiles(ASB_result)
  mn = ASB.mean(ASB_result)
  over5 = ASB.over5(ASB_result)
  over10 = ASB.over10(ASB_result)
  ASB_statistics = c(quant, mn, over5, over10)
  return(ASB_statistics)
}

# -------------------------------------------------------------------------

ASB.out = function(ASB_result, ASB_statistics){
  output = list("Matrix" =  ASB_result, "Stat" = ASB_statistics)
  return(output)
}


