library(doParallel)

learn.model.laplace <- function(X, psi, Phi, alpha, m, lambda, a, b, l, mprior, eps, niter)
{
  #Function to calculate log probability of normal gamma
  loglapgam <- function(m, lambda, mprior, l, a, b) {
    return(
      log(b) * a +
        log(sqrt(l)) -
        lgamma(a) -
        log(2*lambda) * (a - 0.5) -
        b / lambda -
        sqrt(l) * abs(m - mprior)  / (lambda)
    )
  }
  
  # Function to calculate the expected -log(probabity(Normal gamma))
  NLentropy <- function(lambda, bhat, ahat) {
    return(
      log((2 * exp(1) * lambda))  +
        0.5 * (bhat / ahat) +
        lgamma(ahat) -
        log(bhat) +
        ahat * (1 - 1 / bhat) +
        (1 - ahat) * digamma(ahat)
    )
  }
  
  # Function to calculate Expectation of log phi
  E_logphi <- function(alpha, phi_i, exp_gamma) {
    return(lgamma(sum(alpha))  - sum(lgamma(alpha)) + #expectation of log multinomial dirichlet prior
             sum((alpha - 1) * exp_gamma) +
             sum(phi_i * exp_gamma)) # Expectation of y given expected log phi
  }
  
  # Function to calculate Phi entropy
  phi_ientropy <- function(phi_i, var_gamma, exp_gamma) {
    return(sum(lgamma(var_gamma)) - lgamma(sum(var_gamma)) - # Entropy of Dirichlet prior for expected phi
             sum((var_gamma - 1) * exp_gamma) -
             sum(phi_i * log(phi_i))                # Entropy of log multinomial for expected y
    )
  }
  
  
  #Function to calculate f(x_i | m_k, lambda_k, phi_k_)
  Qfx_yi <- function(X_i,
                     psi_i,
                     theta_i,
                     m,
                     lambda,
                     exp_gamma)
  {
    return(sum(vapply(
      1:2,
      FUN.VAL = vector(mode = 'double', length = 1),
      FUN = function(j)
        -abs(X_i[j]-theta_i)/sqrt(psi[j]) - log(2*sqrt(psi_i[j]))
    )) + 
      -abs(theta_i-m)/(lambda) - log(2*lambda) +
      exp_gamma)
    # Vectorised Calculation
  }
  
  
  # Function to calculate E(theta_bar |  X, psi, Phi, m, lambda)
  thetabarupdate <-
    function(phi_i,
             X_i,
             psi_i,
             m,
             lambda) {
      w <- phi_i / lambda
      return((sum(w * m) + sum(X_i / psi_i)) / (sum(w) + sum(1 / psi_i)))
      # Vectorised Calculation
    }
  
  # Function to calculate document level ELBO
  doc.lik <-
    function(X_i,
             psi_i,
             theta_i,
             m,
             lambda,
             phi_i,
             a,
             b,
             l,
             mprior,
             exp_gamma,
             var_gamma) {
      ahat <-
        a  + 0.5 * phi_i
      
      bhat <-
        b + 0.5 * (phi_i * (abs(theta_i - m)    + (sqrt(l) * phi_i * abs(m - mprior)) / (phi_i + sqrt(l))))
      
      return(
        sum(phi_i * (
          sum(vapply(
            1:2,
            FUN.VAL = vector(mode = 'double', length = 1.0),
            FUN = function(j)
             -abs(X_i[j]-theta_i)/sqrt(psi[j]) - log(2*sqrt(psi_i[j]))
          )) - abs(theta_i - m)/(lambda) - log(2*lambda) +
            loglapgam(m, lambda, mprior, l, a, b)
        )) +
          E_logphi(alpha, phi_i, exp_gamma)  +
          phi_ientropy(phi_i, var_gamma, exp_gamma) +
          sum((NLentropy(lambda, bhat, ahat)))
      )
    }
  
  # SNP level expectation or parameters
  SNP.expectation.laplace <-
    function(X_i, psi_i, m, lambda, a, b, l, mprior) {
      iter <- 1
      converged <- 1
      old.ll <- 1
      doc.ll <- 0
      phi_i <-  rep(1 / 3, 3)
      var_gamma <-
        alpha +  rep(1 / 3, 3)
      exp_gamma <-
        digamma(var_gamma) - digamma(sum((var_gamma)))
      log_phi <-
        log(phi_i)
      theta_i <-
        thetabarupdate(phi_i, X_i, psi_i, m, lambda)
      while (converged > eps & iter < niter)  {
        old_phi <-
          phi_i
        
        fx_yi <-
          Qfx_yi(X_i, psi_i, theta_i, m, lambda, exp_gamma)   #derivative of expected probability of y
        
        log_phi <-
          (fx_yi)  -
          (max(fx_yi) +
             log(sum(exp(
               fx_yi - max(fx_yi)
             )))) #normalise in log space
        
        log_phi[log_phi < (-200)] <-
          (-200)  # avoid NAN with zero probability
        phi_i <-
          exp(log_phi)
        var_gamma <-
          var_gamma + (phi_i - old_phi)
        exp_gamma <-
          digamma(var_gamma) - digamma(sum((var_gamma)))
        
        theta_i <- thetabarupdate(phi_i, X_i, psi_i, m, lambda)
        
        doc.ll <- doc.lik(X_i,
                          psi_i,
                          theta_i,
                          m,
                          lambda,
                          phi_i,
                          a,
                          b,
                          l,
                          mprior,
                          exp_gamma,
                          var_gamma)
        
        converged <- abs((doc.ll - old.ll) / old.ll)
        old.ll <- doc.ll
        iter <- iter + 1
        cat(phi_i, theta_i, doc.ll, '\n')
      }
      return(c(phi_i, theta_i , doc.ll))
    }
  ncores <- parallel::detectCores(logical = T)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  
  converged <- 1
  oldll <- 1
  ll <- 0.5
  ll.hist <- c()
  EM.iter <- 0
  
  while (abs(converged) > eps & EM.iter < niter)
  {
    ####----Expectation SNP level
    SNP.expectation.results <-
      foreach::foreach(
        position = 1:nrow(X),
        .combine = rbind,
        .inorder = T,
        .multicombine = T
      ) %dopar% {
        SNP.expectation.laplace(X_i = X[position,],
                        psi_i = psi[position,],
                        m,
                        lambda,
                        a,
                        b,
                        l,
                        mprior)
      }
    ####---- Overall Maximation
    gamma <- colSums(SNP.expectation.results[, 1:3]) #/nrow(X)
    Phi <-
      (gamma) / sum((gamma)) # exp(digamma(gamma) - digamma(sum(gamma)))
    Phi_var <- Phi * (1 - Phi) / (sum(gamma + alpha) + 1)
    m <-
      ((
        t(SNP.expectation.results[, 4]) %*% (SNP.expectation.results[, 1:3])
      ) + (sqrt(l) * mprior)) /
      (colSums(SNP.expectation.results[, 1:3]) + l)
    m[2] <- 0
    bhat <-
      b + 0.5 * (colSums((abs(
        SNP.expectation.results[, 4] - matrix(rep(m, nrow(X)), byrow = T , nrow = nrow(X))
      ) 
      ) * SNP.expectation.results[, 1:3])  +
        (sqrt(l) * colSums(SNP.expectation.results[, 1:3]) * (abs(m - mprior))) /
        (colSums(SNP.expectation.results[, 1:3]) + sqrt(l)))
    ahat <- (colSums(SNP.expectation.results[, 1:3]) / 2) + a
    lambda <-  bhat / ahat
    cat(gamma, Phi, m, lambda, '\n')
    ll <- sum(SNP.expectation.results[5, ])
    cat(ll, '\n')
    converged <- (ll - oldll) / oldll
    oldll <- ll
    ll.hist <- c(ll.hist, ll)
    EM.iter <- EM.iter + 1
  }
  results <- list(Phi = Phi, 
                  Phi_var = Phi_var,
                  theta = SNP.expectation.results[, 4],
                  m = m,
                  lambda = lambda,
                  phi_i = SNP.expectation.results[, 1:3],
                  ll.hist = ll.hist)
  return(results)
}

setwd(
  "C:/Users/mczcjc/OneDrive - The University of Nottingham/Courses/MSc Statistics/Summer project"
)
load("final_data.RData")
X <- as.matrix(cbind(final_data$ic_beta, final_data$gw_beta))
psi <- as.matrix(cbind(final_data$ic_se ^ 2, final_data$gw_se ^ 2))

# log(1.05)
theta <- as.matrix(rowMeans(X))

#Initialisation
niter <- 100
Phi <- rep(1 / 3, 3)

m <- c(-0.02, 0, 0.02)
lambda = c(0.01, 0.01, 0.01)

# hyperparameters
alpha = rep(1 / 2, 3)  
a <- c(4, 4, 4)
b <- c(0.15, 0.35, 0.29)
l <- c(10, 10, 10)
mprior <- c(-0.08, 0, 0.14)
eps <- 1e-6

base.model.laplace <- learn.model.laplace(
  l = l,
  X = X,
  psi = psi,
  Phi = Phi,
  alpha = alpha,
  m = m,
  lambda = lambda,
  a = a,
  b = b,
  mprior = mprior,
  eps = eps, 
  niter = niter
)


dlaplace <- function(x,m,s) {
  return(  exp(-abs(x - m)/(s))/(2*s))
}

# Summary measures
png("Thetaposterior_laplace.png")
hist(base.model.laplace$theta,
     breaks = 1000,
     freq = F,
     main = "",
     xlim = c(-0.3,0.3),
     ylim = c(0,60),
     xlab = expression('SNP posterior mean log odds ratio'(theta[i])),
     ylab = 'Density')
lines(
  x = seq(-0.3, 0.3, 0.001),
  y =dlaplace(seq(-0.3, 0.3, 0.001), m = base.model.laplace$m[1], s =base.model.laplace$lambda[1]),
  type = 'l',
  col = 'green'
)
lines(
  x = seq(-0.3, 0.3, 0.001),
  y = dlaplace(seq(-0.3, 0.3, 0.001), m = base.model.laplace$m[2], s = base.model.laplace$lambda[2]),
  col = 'blue'
)
lines(
  x = seq(-0.3, 0.3, 0.001),
  y = dlaplace(seq(-0.3, 0.3, 0.001), m = base.model.laplace$m[3], s = base.model.laplace$lambda[3]),
  col = 'red'
)
legend(
  "topright",
  c(
    expression("Negative posterior m"[1]),
    expression("Null posterior m"[2]),
    expression("Positive posterior m"[3])
  ),
  col = c("green", "blue", "red"),
  lty = c(1, 1, 1)
)
dev.off()



base.model.laplace$Phi
base.model.laplace$Phi + sqrt(base.model.laplace$Phi_var) * 1.96
base.model.laplace$Phi - sqrt(base.model.laplace$Phi_var) * 1.96

base.model.laplace$Phi * nrow(X)
(base.model.laplace$Phi + sqrt(base.model.laplace$Phi_var) * 1.96) * nrow(X)
(base.model.laplace$Phi - sqrt(base.model.laplace$Phi_var) * 1.96) * nrow(X)


  sum(base.model.laplace$phi_i[,1] > 0.90)
  sum(base.model.laplace$phi_i[,2] > 0.90)
  sum(base.model.laplace$phi_i[,3] > 0.90)
  
  sum(base.model.laplace$phi_i[,1] > 0.975)
  sum(base.model.laplace$phi_i[,2] > 0.975)
  sum(base.model.laplace$phi_i[,3] > 0.975)

png('VBll_laplace.png')
plot(base.model.laplace$ll.hist,
     type = 'l',
     xlab = 'iterations',
     ylab = 'ELBO')
dev.off()
