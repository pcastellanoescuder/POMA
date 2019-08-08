
PomaImpute <- function(data,
                       ZerosAsNA = TRUE,
                       RemoveNA = TRUE,
                       cutoff = 20,
                       method = c("none", "half_min", "median", "mean", "min", "knn")){

  data <- as.data.frame(data)
  samples_groups <- data[, 1:2]
  to_imp_data <- data[, c(3:ncol(data))]
  colnames(samples_groups)[2] <- "Group"

  if (ZerosAsNA == TRUE){
    to_imp_data[to_imp_data == 0] <- NA
    to_imp_data <- cbind(samples_groups, to_imp_data)
  }

  else{
    to_imp_data <- cbind(samples_groups, to_imp_data)
  }

  if (RemoveNA == TRUE){
    count_NA <- aggregate(. ~ Group, data = to_imp_data[, 2:ncol(to_imp_data)],
                          function(x) {100*(sum(is.na(x))/(sum(is.na(x))+sum(!is.na(x))))},
                          na.action = NULL)
    count_NA$Group <- NULL
    supress <- as.data.frame(lapply(count_NA, function(x) all(x > cutoff)))
    supress <- unlist(supress)
    depurdata <- to_imp_data[, 3:ncol(to_imp_data)][!supress]

    if (method == "none"){
      depurdata[is.na(depurdata)] <- 0
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "half_min"){
      depurdata <- apply(depurdata, 2, function(x) "[<-"(x, !x | is.na(x),
                                                         min(x[x >= 0], na.rm = TRUE) / 2))
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "median"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), median(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "mean"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), mean(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "min"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), min(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "knn"){
      depurdata <- t(depurdata)
      datai <- impute::impute.knn(as.matrix(depurdata))
      depurdata <- t(datai$data)
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }
    return (depurdata)

  }

  else{

    depurdata <- to_imp_data[, 3:ncol(to_imp_data)]

    if (method == "none"){
      depurdata[is.na(depurdata)] <- 0
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "half_min"){
      depurdata <- apply(depurdata, 2, function(x) "[<-"(x, !x | is.na(x),
                                                         min(x[x >= 0], na.rm = TRUE) / 2))
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "median"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), median(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "mean"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), mean(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "min"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), min(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "knn"){
      depurdata <- t(depurdata)
      datai <- impute::impute.knn(as.matrix(depurdata))
      depurdata <- t(datai$data)
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }
    return (depurdata)
  }

}

