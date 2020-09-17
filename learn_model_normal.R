
library(doParallel)

learn.model <- function(X, psi, Phi, alpha, m, lambda, a, b, l, mprior, eps, niter)
{
  #Function to calculate log probability of normal gamma
  lognormgam <- function(m, lambda, mprior, l, a, b) {
    return(
      log(b) * a +
        log(sqrt(l)) -
        lgamma(a) -
        log(sqrt(2 * pi)) -
        log(lambda) * (a - 0.5) -
        b / lambda -
        l * (m - mprior) ^ 2 / (2 * lambda)
    )
  }
  
  # Function to calculate the expected -log(probabity(Normal gamma))
  NGentropy <- function(lambda, bhat, ahat) {
    return(
      log(sqrt(2 * pi * lambda))  +
        0.5 * (bhat / ahat) +
        lgamma(ahat) -
        log(bhat) +
        ahat * (1 - 1 / bhat) +
        (1 - ahat) * digamma(ahat)
    )
  }
  
  # Function to calculate the expectation of log Phi from Dirichlet multinomial
  E_logphi <- function(alpha, phi_i, exp_gamma) {
    return(lgamma(sum(alpha))  - sum(lgamma(alpha)) + #expectation of log multinomial dirichlet prior
             sum((alpha - 1) * exp_gamma) +
             sum(phi_i * exp_gamma)) # Expectation of y given expected log phi
  }
  
  # Entropy of Dirichlet multinomial
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
        dnorm(
          x = X_i[j],
          mean = theta_i,
          sd = sqrt(psi[j]),
          log = T
        )
    )) +
      dnorm(
        x = theta_i,
        mean = m,
        sd = sqrt(lambda),
        log = T
      )  +
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
  
  # Function to calculate document component of ELBO
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
        b + 0.5 * (phi_i * ((theta_i - m) ^ 2   + (l * phi_i * (m - mprior) ^ 2) / (phi_i + l)))
      
      return(
        sum(phi_i * (
          sum(vapply(
            1:2,
            FUN.VAL = vector(mode = 'double', length = 1.0),
            FUN = function(j)
              dnorm(
                x = X_i[j],
                mean = theta_i,
                sd = sqrt(psi[j]),
                log = T
              )
          )) +
            dnorm(
              x = theta_i,
              mean = m,
              sd = sqrt(lambda),
              log = T
            ) +
            lognormgam(m, lambda, mprior, l, a, b)
        )) +
          E_logphi(alpha, phi_i, exp_gamma)  +
          phi_ientropy(phi_i, var_gamma, exp_gamma) +
          sum((NGentropy(lambda, bhat, ahat)))
      )
    }
  
  # Expectation of parameters at document level
  SNP.expectation <-
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
        SNP.expectation(X_i = X[position,],
                        psi_i = psi[position,],
                        m,
                        lambda,
                        a,
                        b,
                        l,
                        mprior)
      }
    ####---- Overall Maximation of normal gamma prior
    gamma <- colSums(SNP.expectation.results[, 1:3]) 
    Phi <-
      (gamma) / sum((gamma))
    Phi_var <- Phi * (1 - Phi) / (sum(gamma + alpha) + 1)
    m <-
      ((
        t(SNP.expectation.results[, 4]) %*% (SNP.expectation.results[, 1:3])
      ) + (l * mprior)) /
      (colSums(SNP.expectation.results[, 1:3]) + l)
    m[2] <- 0
    bhat <-
      b + 0.5 * (colSums(((
        SNP.expectation.results[, 4] - matrix(rep(m, nrow(X)), byrow = T , nrow = nrow(X))
      ) ^ 2
      ) * SNP.expectation.results[, 1:3])  +
        (l * colSums(SNP.expectation.results[, 1:3]) * ((m - mprior) ^ 2)) /
        (colSums(SNP.expectation.results[, 1:3]) + l))
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

#Load data
setwd(
  "C:/Users/mczcjc/OneDrive - The University of Nottingham/Courses/MSc Statistics/Summer project"
)
load("final_data.RData")
X <- as.matrix(cbind(final_data$ic_beta, final_data$gw_beta))
psi <- as.matrix(cbind(final_data$ic_se ^ 2, final_data$gw_se ^ 2))


theta <- as.matrix(rowMeans(X))

# Initialise parameters
niter <- 100
Phi <- rep(1 / 3, 3)
m <- c(-0.02, 0, 0.02)
lambda = c(0.01, 0.01, 0.01)

# Hyperparameters
alpha = rep(1 / 2, 3)  

a <- c(4, 4, 4)
b <- c(0.008, 0.052, 0.036)
l <- c(10, 10, 10)
mprior <- c(-0.08, 0, 0.14)
eps <- 1e-6

base.model <- learn.model(
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


# Summary statistics
base.model$Phi
base.model$Phi + sqrt(base.model$Phi_var) * 1.96
base.model$Phi - sqrt(base.model$Phi_var) * 1.96

base.model$Phi * nrow(X)
(base.model$Phi + sqrt(base.model$Phi_var) * 1.96) * nrow(X)
(base.model$Phi - sqrt(base.model$Phi_var) * 1.96) * nrow(X)


par(mar = c(0, 0, 0, 0))
png("NGposterior.png")
plot(
  x = seq(-1, 1, 0.001),
  y = dnorm(seq(-1, 1, 0.001), mean = base.model$m[2], sd = sqrt(base.model$lambda[2])),
  type = 'l',
  xlab = 'Value of posterior m',
  ylab = 'Probability',
  col = 'blue'
)
lines(
  x = seq(-1, 1, 0.001),
  y = dnorm(seq(-1, 1, 0.001), mean = base.model$m[1], sd = sqrt(base.model$lambda[1])),
  col = 'green'
)
lines(
  x = seq(-1, 1, 0.001),
  y = dnorm(seq(-1, 1, 0.001), mean = base.model$m[3], sd = sqrt(base.model$lambda[3])),
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

png("Thetaposterior.png")
hist(base.model$theta,
     breaks = 1000,
     freq = F,
     main = "",
     xlim = c(-0.1,0.1),  
     xlab = expression('SNP posterior mean log odds ratio'(theta[i])),
     ylab = 'Density')
lines(
  x = seq(-0.1, 0.1, 0.001),
  y = dnorm(seq(-0.1, 0.1, 0.001), mean = base.model$m[1], sd = sqrt(base.model$lambda[1])),
  type = 'l',
  col = 'green'
)
lines(
  x = seq(-0.1, 0.1, 0.001),
  y = dnorm(seq(-0.1, 0.1, 0.001), mean = base.model$m[2], sd = sqrt(base.model$lambda[2])),
  col = 'blue'
)
lines(
  x = seq(-0.1, 0.1, 0.001),
  y = dnorm(seq(-0.1, 0.1, 0.001), mean = base.model$m[3], sd = sqrt(base.model$lambda[3])),
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

base.model$m
base.model$lambda

png('VBll.png')
plot(-base.model$ll.hist,
     type = 'l',
     xlab = 'iterations',
     ylab = 'ELBO')
dev.off()

png('phi_i.png')
par(mfrow = c(2,2))
hist(base.model$phi_i[,1],
     xlab = 'Probability of negative SNP association',
     ylab = 'Frequency',
     breaks = 100,
     main = 'Negative association',
     freq = T)
hist(base.model$phi_i[,3],
     xlab = 'Probability of positive SNP association',
     ylab = 'Frequency',
     breaks = 100,
     main = 'Positive association',
     freq = T)
hist(base.model$phi_i[base.model$phi_i[,1]>0.9,1],
     xlab = 'Probability of negative SNP association',
     ylab = 'Frequency',
     breaks = 100,
     main = 'SNPs probability >0.9',
     freq = T)
hist(base.model$phi_i[base.model$phi_i[,3]>0.9,3],
     xlab = 'Probability of positive SNP association',
     ylab = 'Frequency',
     breaks = 100,
     main = 'SNPs probability >0.9',
     freq = T)
dev.off()
par(mfrow = c(1,1))


# Summary number of SNPs allocated at different thresholds
sum(base.model$phi_i[,1] > 0.90)
sum(base.model$phi_i[,2] > 0.90)
sum(base.model$phi_i[,3] > 0.90)

sum(base.model$phi_i[,1] > 0.975)
sum(base.model$phi_i[,2] > 0.975)
sum(base.model$phi_i[,3] > 0.975)


# sensitivity analyses:

l <- c(0.001,0.01,0.1,1,10,100,1000,10000)

l.sens <- vapply(1:8, FUN.VALUE = vector(mode = 'list',length = 7),FUN = function(x) learn.model(l = l[x], X=X,psi=psi,Phi=Phi, alpha=alpha, m=m, lambda=lambda,a=a,b=b,mprior=mprior,
                                                                                                 eps = eps, 
                                                                                                 niter = niter) )

l <- c(10,10,10)
b <- c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10)


b.sens <- vapply(1:8, FUN.VALUE = vector(mode = 'list',length = 7),FUN = function(x) learn.model(l = l, X=X,psi=psi,Phi=Phi, alpha=alpha, m=m, lambda=lambda,a=a,b=b[x],mprior=mprior,
                                                                                                 eps = eps, 
                                                                                                 niter = niter) )

b <- c(0.002,0.013,0.009)
mprior <- matrix (c(-0.5,-0.4,-0.3,-0.2,-0.1,-0.01,0,0,0,0,0,0,0.5,0.4,0.3,0.2,0.1,0.01),byrow = F,nrow = 6)

mprior.sens <- vapply(1:6, FUN.VALUE = vector(mode = 'list',length = 7),FUN = function(x) learn.model(l = l, X=X,psi=psi,Phi=Phi, alpha=alpha, m=m, lambda=lambda,a=a,b=b,mprior=mprior[x,],
                                                                                                      eps = eps, 
                                                                                                      niter = niter) )

mprior <- c(-0.08, 0, 0.14)
a <- c(0.001,0.01,0.1,1,10,100,1000)

a.sens <- vapply(1:6, FUN.VALUE = vector(mode = 'list',length = 7),FUN = function(x) learn.model(l = l, X=X,psi=psi,Phi=Phi, alpha=alpha, m=m, lambda=lambda,a=a[x],b=b,mprior=mprior,
                                                                                                 eps = eps, 
                                                                                                 niter = niter) )

a <- c(4,4,4)
alpha <- c(5,2,1,1/2,1/3,1/5)

alpha.sens <- vapply(1:6, FUN.VALUE = vector(mode = 'list',length = 7),FUN = function(x) learn.model(l = l, X=X,psi=psi,Phi=Phi, alpha=alpha[x], m=m, lambda=lambda,a=a,b=b,mprior=mprior,
                                                                                                     eps = eps, 
                                                                                                     niter = niter) )

alpha <- rep(1/2,3)

# Summary of sensitivity analyses
library(tidyr)
library(dplyr)
library(ggplot2)
lmat <- (matrix(unlist(l.sens[1,1:8]), nrow = 8, byrow = T))
rownames(lmat) <- c('0.001','0.01','0.1','1','10','100','1000','10000')
colnames(lmat) <- c('Negative','Null','Positive')
lplot <- data.frame(lmat) %>% tibble::rownames_to_column("l.prior") %>%
  pivot_longer(-l.prior,names_to = "Phi") %>%
  ggplot() + geom_bar(aes(x=Phi,y=value,fill=l.prior), position = 'dodge', stat = 'identity') +
  ylab("Posterior probability") 
print(lplot)

ggplot2::ggsave('lsensplot.png',plot =lplot, dpi = 300, height = 6.5, units = 'cm')

bmat <- (matrix(unlist(b.sens[1,1:8]), nrow = 8, byrow = T))
rownames(bmat) <- c('0.000001','0.00001','0.0001','0.001','0.01','0.1','1','10')
colnames(bmat) <- c('Negative','Null','Positive')
bplot <- data.frame(bmat) %>% tibble::rownames_to_column("b.prior") %>%
  pivot_longer(-b.prior,names_to = "Phi") %>%
  ggplot() + geom_bar(aes(x=Phi,y=value,fill=b.prior), position = 'dodge', stat = 'identity')+
  ylab("Posterior probability")
print(bplot)

ggplot2::ggsave('bsensplot.png',plot =bplot, dpi = 300, height = 6.5, units = 'cm')


alphamat <- (matrix(unlist(alpha.sens[1,1:6]), nrow = 6, byrow = T))
rownames(alphamat) <- c('5','2','1','1/2','1/3','1/5')
colnames(alphamat) <- c('Negative','Null','Positive')
alphaplot <- data.frame(alphamat) %>% tibble::rownames_to_column("alpha.prior") %>%
  pivot_longer(-alpha.prior,names_to = "Phi") %>% mutate(alpha.prior = factor(alpha.prior,levels = c('5','2','1','1/2','1/3','1/5'),labels = c('5','2','1','1/2','1/3','1/5'))) %>%
  ggplot() + geom_bar(aes(x=Phi,y=value,fill=alpha.prior), position = 'dodge', stat = 'identity') +
  ylab("Posterior probability")
print(alphaplot)
ggplot2::ggsave('alphasensplot.png',plot =alphaplot, dpi = 300, height = 8, units = 'cm')

alphamat*9716

vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(alpha.sens[[6,x]][,1] > 0.90))
vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(alpha.sens[[6,x]][,3] > 0.90))


vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(alpha.sens[[6,x]][,1] > 0.975))
vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(alpha.sens[[6,x]][,3] > 0.975))

mmat <- (matrix(unlist(mprior.sens[1,1:6]), nrow = 6, byrow = T))
rownames(mmat) <- c('0.5','0.4','0.3','0.2','0.1','0.01')
colnames(mmat) <- c('Negative','Null','Positive')
mplot <- data.frame(mmat) %>% tibble::rownames_to_column("m.prior") %>%
  pivot_longer(-m.prior,names_to = "Phi") %>%
  ggplot() + geom_bar(aes(x=Phi,y=value,fill=m.prior), position = 'dodge', stat = 'identity')+
  ylab("Posterior probability")
print(mplot)
ggplot2::ggsave('msensplot.png',plot =mplot, dpi = 300, height = 5.5, units = 'cm')



amat <- (matrix(unlist(a.sens[1,1:6]), nrow = 6, byrow = T))
rownames(amat) <- c('0.001','0.01','0.1','1','10','100')
colnames(amat) <- c('Negative','Null','Positive')
aplot <- data.frame(amat) %>% tibble::rownames_to_column("a.prior") %>%
  pivot_longer(-a.prior,names_to = "Phi") %>%
  ggplot() + geom_bar(aes(x=Phi,y=value,fill=a.prior), position = 'dodge', stat = 'identity')+
  ylab("Posterior probability")
print(aplot)
ggplot2::ggsave('asensplot.png',plot =aplot, dpi = 300,  height = 5.5, units = 'cm')

vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(mprior.sens[[6,x]][,1] > 0.90))
vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(mprior.sens[[6,x]][,3] > 0.90))

vapply(1:7, FUN.VAL = 1.0, FUN = function(x) sum(l.sens[[6,x]][,1] > 0.90))
vapply(1:7, FUN.VAL = 1.0, FUN = function(x) sum(l.sens[[6,x]][,3] > 0.90))

vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(a.sens[[6,x]][,1] > 0.90))
vapply(1:6, FUN.VAL = 1.0, FUN = function(x) sum(a.sens[[6,x]][,3] > 0.90))

vapply(1:7, FUN.VAL = 1.0, FUN = function(x) sum(b.sens[[6,x]][,1] > 0.90))
vapply(1:7, FUN.VAL = 1.0, FUN = function(x) sum(b.sens[[6,x]][,3] > 0.90))



save(l.sens,a.sens, b.sens, alpha.sens, a.sens,mprior.sens, file = "Sensitivity.Rdata")
