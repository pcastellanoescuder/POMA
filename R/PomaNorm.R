
PomaNorm <- function(data,
                     method = c("none", "auto_scaling", "level_scaling", "log_scaling",
                                "log_transformation", "vast_scaling","log_pareto"),
                     round = 3){

  data <- as.data.frame(data)
  samples_groups <- data[, 1:2]
  to_norm_data <- data[, c(3:ncol(data))]

  # remove columns that only have zeros
  to_norm_data <- to_norm_data[, apply(to_norm_data, 2, function(x) !all(x==0))]

  if (method == "none"){
    not_norm_data <- round(to_norm_data, round)
    normalized <- cbind(samples_groups, not_norm_data)
    return (normalized)
  }

  else if (method == "auto_scaling"){
    auto_scaling_data <- round(apply(to_norm_data, 2, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)), round)
    normalized <- cbind(samples_groups, auto_scaling_data)
    return (normalized)
  }

  else if (method == "level_scaling"){
    level_scaling_data <- round(apply(to_norm_data, 2, function(x) (x-mean(x,na.rm=T))/mean(x,na.rm=T)), round)
    normalized <- cbind(samples_groups, level_scaling_data)
    return (normalized)
  }

  else if (method == "log_scaling"){
    log_scaling_data <- round(apply(to_norm_data, 2, function(x) (log10(x+1)-mean(log10(x+1),na.rm=T))/sd(log10(x+1),na.rm=T)), round)
    normalized <- cbind(samples_groups, log_scaling_data)
    return (normalized)
  }

  else if (method == "log_transformation"){
    log_transformation_data <- round(apply(to_norm_data, 2, function(x) (log10(x+1))), round)
    normalized <- cbind(samples_groups, log_transformation_data)
    return (normalized)
  }

  else if (method == "vast_scaling"){
    vast_scaling_data <- round(apply(to_norm_data, 2, function(x) ((x-mean(x,na.rm=T))/sd(x,na.rm=T))*(mean(x,na.rm=T)/sd(x,na.rm=T))), round)
    normalized <- cbind(samples_groups, vast_scaling_data)
    return (normalized)
  }

  else if (method == "log_pareto"){
    log_pareto_data <- round(apply(to_norm_data, 2, function(x) (log10(x+1)-mean(log10(x+1),na.rm=T))/sqrt(sd(log10(x+1),na.rm=T))), round)
    normalized <- cbind(samples_groups, log_pareto_data)
    return (normalized)
  }
}

