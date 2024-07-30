### 
# from INLAMSM package

###

PMCAR_Leroux <- function (cmd = c("graph", "Q", "mu", "initial", "log.norm.const", 
                  "log.prior", "quit"), theta = NULL) 
{
  interpret.theta <- function() {
    # theta(internal scale) to parameter.
    
    # mprec: precision vector.
    mprec <- sapply(theta[1:2], 
                    function(x) { exp(x) })
    
    # corre: correlation. theta = logit( (rho + 1) / 2)
    corre <- sapply(theta[3],
                    function(x) { (2 * exp(x))/(1 + exp(x)) - 1 })
    
    # lambda
    lambda <- exp(theta[4]) / (1 + exp(theta[4]))
    
    param <- c(mprec, corre, lambda)
    
    # M: variance-covariance matrix
    n <- (k - 1) * k/2
    M <- diag(1, k)
    M[lower.tri(M)] <- param[k+1] # corre
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
    st.dev <- 1/sqrt(param[1:k]) # sd = 1/√tau
    st.dev.mat <- 
      matrix(st.dev, ncol = 1) %*% matrix(st.dev, nrow = 1)
    M <- M * st.dev.mat
    
    # PREC: precision matrix of k by k
    PREC <- solve(M)
    return(list(param = param, PREC = PREC))
  }
  graph <- function() {
    # a 0/1 representation of precision matrix
    PREC <- matrix(1, ncol = k, nrow = k)
    G <- kronecker(PREC, (Matrix::Diagonal(nrow(W), 1) + W))
    G <- replace(G, G > 0, 1)
    return(G)
  }
  Q <- function() {
    param <- interpret.theta()
    require(Matrix)
    # kronecker(PREC, (1-λ)I + λ(D-W) )
    Q_Leroux <- (1 - param$param[4]) * (Diagonal(nrow(W), 1)) +
      (param$param[4] * (Diagonal(nrow(W), x = rowSums(W)) - W))
    Q <- kronecker(param$PREC, Q_Leroux)
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
    
    # tau, rho
    val <- 
      log(MCMCpack::dwish(W = param$PREC, v = 4, S = diag(rep(1, 2)))) +
      sum(theta[1:2]) + # Jacobian
      log(2) + log( exp(theta[3]) ) - 2 * log( 1 + exp(theta[3]) ) # Jacobian

    # lambda ~ uniform(a = 0, b = 1)
    # lambda
    val <- val + log(1) + 
      log(param$param[4]) + log(1 - param$param[4]) # Jacobian
      
    return(val)
  }
  initial <- function() {
    return( c(0, 0, 0, 0)) # internal scale
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
