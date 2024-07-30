
# Leroux Proper CAR model
PCAR_Leroux <- function (
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", 
            "log.prior", "quit"), theta = NULL) {
  interpret.theta <- function() {

    return(
      list(prec = exp(theta[1L]), # tau = exp(theta) <=> theta = log(tau)
           lambda = exp(theta[2L]) / ( 1 + exp(theta[2L]) ) 
           # lambda = logistic(theta) <=> theta = logit(lambda)
           )
    )
  
  }
  graph <- function() {

    #Graph of precision function; i.e., a 0/1 representation of precision matrix
    G <- diag(nrow(W), x = 1) + W
    
    # need when W is not binary
    G <- replace(G, G > 0, 1) 
    return(G)
  }
  Q <- function() {
    param <- interpret.theta()

    # tau * (1-lambda) * I + lambda * (D - W) 
    R <- (1 - param$lambda) * (Diagonal(nrow(W), 1)) +
      (param$lambda * (Diagonal(nrow(W), x = rowSums(W)) - W))
    Q <- param$prec * R
    
    return(Q)
  }
  
  mu <- function() {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    val <- numeric(0)
    return(val)
  }
  
  log.prior <- function() {
    param <- interpret.theta()
    
    # tau
    val <- 
      dgamma(param$prec, shape = 1, rate = 0.01, log = TRUE) +
      log( exp(theta[1L]) ) # Jacobian
    
    # lambda ~ uniform(a = 0, b = 1)
    val <- val + log(1) + 
      log(param$lambda) + log(1 - param$lambda) # Jacobian
    
    return(val)
  }
  
  # initial in internal scale
  # tau_0 = exp(0) = 1
  initial <- function() {
    return(c(0, 0)) 
  }
  quit <- function() {
    return(invisible())
  }
  if (as.integer(R.version$major) > 3) {
    if (!length(theta)) 
      theta <- initial()
  }
  else {
    if (is.null(theta)) {
      theta <- initial()
    }
  }
  val <- do.call(match.arg(cmd), args = list())
  return(val)
}
