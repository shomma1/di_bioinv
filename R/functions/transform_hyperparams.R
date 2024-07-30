
# MLCAR: Multivariate Leroux CAR
transform.hyper.mlerouxcar <- function (res) {
  
  margs_internal <- res$marginals.hyperpar
  
  par.names <- c("tau1", "tau2", "rho","lambda")
  n_tau <- 1:2
  n_rho <- 3
  n_lambda <- 4
  
  ## transform internal scale marginals
  margs1 <- lapply(margs_internal[n_tau],
                   function(xx) {
                     inla.tmarginal(function(x) {exp(x)}, xx)})
  
  margs2 <- lapply(margs_internal[n_rho],
                   function(xx) {
                     inla.tmarginal(function(x) {
                       ((2 * exp(x))/(1 + exp(x)) - 1)}, xx)})
  
  margs3 <- lapply(margs_internal[n_lambda],
                   function(xx) {
                     inla.tmarginal(function(x) {
                       (exp(x) / (1 + exp(x)))}, xx)})
  
  
  margs <- c(margs1, margs2, margs3)
  names(margs) <- par.names
  
  # summary statics
  zmargs <- lapply(margs, inla.zmarginal, silent = T)
  zmargs <- lapply(zmargs, unlist)
  zmargs <- do.call(rbind, zmargs)
  rownames(zmargs) <- par.names
  zmargs <- as.data.frame(zmargs)
  zmargs <- zmargs[c("mean", "sd", "quant0.025", "quant0.5","quant0.975")]
  
  colnames(zmargs) <- 
    c("mean", "sd", "0.025quant", "0.5quant","0.975quant")
  zmargs <- cbind(zmargs, mode = NA)
  
  return(list("zmargs" = zmargs, "margs" = margs))
}

# PCAR: Leroux
transform.hyper.lerouxcar <- function (res) {
  
  margs_internal <- res$marginals.hyperpar
  
  par.names <- c("tau","lambda")
  n_tau <- 1
  n_lambda <- 2
  
  ## transform internal scale marginals
  margs1 <- lapply(margs_internal[n_tau],
                   function(xx) {
                     inla.tmarginal(function(x) {exp(x)}, xx)})
  
  margs2 <- lapply(margs_internal[n_lambda],
                   function(xx) {
                     inla.tmarginal(function(x) {
                       (exp(x) / (1 + exp(x)))}, xx)})
  
  margs <- c(margs1, margs2)
  names(margs) <- par.names
  
  # summary statics
  zmargs <- lapply(margs, inla.zmarginal, silent = T)
  zmargs <- lapply(zmargs, unlist)
  zmargs <- do.call(rbind, zmargs)
  rownames(zmargs) <- par.names
  zmargs <- as.data.frame(zmargs)
  zmargs <- zmargs[c("mean", "sd", "quant0.025", "quant0.5","quant0.975")]
  
  colnames(zmargs) <- 
    c("mean", "sd", "0.025quant", "0.5quant","0.975quant")
  zmargs <- cbind(zmargs, mode = NA)
  
  return(list("zmargs" = zmargs, "margs" = margs))
}


