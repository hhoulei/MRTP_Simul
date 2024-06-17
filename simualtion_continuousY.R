# Note that the "MRTP" has been changed to "TEMR".

library(MASS)
library(doMC)
library(doParallel)
library(foreach)

Wald_Ratio <- function(betaXG,betaYG,sebetaXG,sebetaYG){
  WR <- betaYG/betaXG
  varWR_2 <- (betaYG^2)*(sebetaXG^2)/(betaXG^4)+
    (sebetaYG^2)/(betaXG^2)
  
  return(data.frame(WR=WR,
                    varWR_2=varWR_2))
}

DataGenerator_wald <- function(N1,N2,N3,N4,g,p1,p2,beta1,beta2,sigma_x1,sigma_x2,
                               sigma2_xi1,sigma2_xi2,sigma2_epi1,sigma2_epi2){
  
  G1 <- NULL
  for(i in 1:g){
    Gg <- rbinom(N1,2,p1)
    G1 <- cbind(G1,Gg)
  }
  
  G2 <- NULL
  for(i in 1:g){
    Gg <- rbinom(N2,2,p2)
    G2 <- cbind(G2,Gg)
  }
  
  sigma_alpha <- matrix(c(sigma_x1,0,0,sigma_x2),nrow=2)
  alpha <- mvrnorm(n=g, rep(0, 2),  sigma_alpha)
  
  U1 <- rnorm(N1,0,1)
  U2 <- rnorm(N2,0,1)
  
  X1 <- G1 %*% alpha[,1] + rnorm(N1,0,1) + U1
  X2 <- G2 %*% alpha[,2] + rnorm(N2,0,1) + U2
  
  G3 <- NULL
  for(i in 1:g){
    Gg <- rbinom(N3,2,p1)
    G3 <- cbind(G3,Gg)
  }
  
  G4 <- NULL
  for(i in 1:g){
    Gg <- rbinom(N4,2,p2)
    G4 <- cbind(G4,Gg)
  }
  
  U3 <- rnorm(N3,0,1)
  U4 <- rnorm(N4,0,1)
  
  X3 <- G3 %*% alpha[,1] + rnorm(N3,0,1) + U3
  X4 <- G4 %*% alpha[,2] + rnorm(N4,0,1) + U4
  
  Y1 <-   beta1*X3 + rnorm(N3,0,1) + U3 #+ G1 %*% gamma1
  Y2 <-   beta2*X4 + rnorm(N4,0,1) + U4 #+ G2 %*% gamma2 
  
  datasimul <- list(G1=G1,
                    G2=G2,
                    X1=X1,
                    X2=X2,
                    G3=G3,
                    G4=G4,
                    Y1=Y1,
                    Y2=Y2)
  return(datasimul)
}

Comp_wald <- function(N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                      sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b){
  
  fdata1 <- DataGenerator_wald(N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                               sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b)
  
  betaX1G <- NULL
  betaX2G <- NULL
  betaY1G <- NULL
  betaY2G <- NULL
  seX1G <- NULL
  seX2G <- NULL
  seY1G <- NULL
  seY2G <- NULL
  for(j in 1:gb){
    FF1 <- lm(fdata1$X1~fdata1$G1[,j])
    betaX1G <- c(betaX1G,FF1$coef[2])
    seX1G <- c(seX1G,summary(FF1)$coef[2,2])
    FF2 <- lm(fdata1$X2~fdata1$G2[,j])
    betaX2G <- c(betaX2G,FF2$coef[2])
    seX2G <- c(seX2G,summary(FF2)$coef[2,2])
    FF3 <- lm(fdata1$Y1~fdata1$G3[,j])
    betaY1G <- c(betaY1G,FF3$coef[2])
    seY1G <- c(seY1G,summary(FF3)$coef[2,2])
    FF4 <- lm(fdata1$Y2~fdata1$G4[,j])
    betaY2G <- c(betaY2G,FF4$coef[2])
    seY2G <- c(seY2G,summary(FF4)$coef[2,2])
  }
  
  res_wald1 <- Wald_Ratio(betaX1G,
                          betaY1G,
                          seX1G,
                          seY1G)
  res_wald2 <- Wald_Ratio(betaX2G,
                          betaY2G,
                          seX2G,
                          seY2G)
  res_wald1 <- res_wald1[order(res_wald1$varWR_2),]
  res_wald2 <- res_wald2[order(res_wald2$varWR_2),]
  sigma_b1 <- sqrt(res_wald1$varWR_2)
  sigma_b2 <- sqrt(res_wald2$varWR_2)
  
  return(list(sigma_b1=sigma_b1,
              sigma_b2=sigma_b2))
}

DataGeneration <- function(g,b1,b2,rho_b,sigma_b1,sigma_b2,siga1,siga2){
  
  gamma1 <- runif(g,siga1,siga2)
  gamma2 <- runif(g,siga1,siga2)
  
  beta1 <- NULL
  beta2 <- NULL
  for(j in 1:g){
    beta <- mvrnorm(1,mu = c(b1,b2),
                    Sigma = matrix(c(sigma_b1[j]^2,rho_b*sigma_b1[j]*sigma_b2[j],
                                     rho_b*sigma_b1[j]*sigma_b2[j],sigma_b2[j]^2),
                                   nrow=2))
    beta1 <- c(beta1,beta[1])
    beta2 <-c(beta2,beta[2]) 
  }
  beta1 <- beta1 + gamma1
  beta2 <- beta2 + gamma2
  
  data_sum <- data.frame(gamma1=gamma1,
                         gamma2=gamma2,
                         beta1=beta1,
                         beta2=beta2,
                         sebeta1=sigma_b1,
                         sebeta2=sigma_b2,
                         rho_b=rho_b)
  return(data_sum)
}

#########IVW############
Trad_IVW <- function(beta,sigma_b){
  
  beta_IVW <- sum(beta*(1/(sigma_b^2)))/sum(1/(sigma_b^2))
  var_IVW <- 1/sum(1/(sigma_b^2))
  pvalue_IVW <- 2 * stats::pnorm(abs(beta_IVW/((var_IVW)^0.5)), lower.tail = FALSE)
  
  resi <- c(beta_IVW, var_IVW, pvalue_IVW)
  return(resi)
}

############### MRTP #############
MLikelihood <- function(beta1,beta2,b2,sigma_b1,sigma_b2,rho_b,gamma1,gamma2){
  
  bb2=b2
  ggamma1=gamma1
  ggamma2=gamma2
  
  nll <- function(bb1){
    det_sigma = (1-rho_b^2)*(sigma_b1^2)
    pp = (beta1-bb1-ggamma1-rho_b*sigma_b1*(beta2-bb2-ggamma2)/sigma_b2)^2
    sum(log(2*pi)+0.5*log(det_sigma)+0.5*pp/det_sigma)
  }
  
  fit <- stats4::mle(minuslog=nll, 
                     start=list(bb1=0))
  res <- fit@coef
  
  fit1 <- stats4::mle(minuslog=nll, 
                      start=list(bb1=0),
                      fixed = list(bb1=0) )
  
  stat_beta1=  2 * (fit@min - fit1@min)
  pvalue_beta1 = pchisq(-stat_beta1,1,lower.tail=F)
  
  return(c(res,pvalue_beta1))
}

########mr-median##########
weighted_median <- function (b_iv, weights) 
{
  betaIV.order <- b_iv[order(b_iv)]
  weights.order <- weights[order(b_iv)]
  weights.sum <- cumsum(weights.order) - 0.5 * weights.order
  weights.sum <- weights.sum/sum(weights.order)
  below <- max(which(weights.sum < 0.5))
  b = betaIV.order[below] + (betaIV.order[below + 1] - betaIV.order[below]) * 
    (0.5 - weights.sum[below])/(weights.sum[below + 1] - 
                                  weights.sum[below])
  return(b)
}

weighted_median_bootstrap <- function (beta1, sigma_b1, weights, nboot) 
{
  med <- rep(0, nboot)
  for (i in 1:nboot) {
    betaIV.boot = stats::rnorm(length(beta1), mean = beta1, 
                               sd = sigma_b1)
    med[i] = weighted_median(betaIV.boot, weights)
  }
  return(stats::sd(med))
}

MR_MEDIAN <- function(beta1,sigma_b1){
  
  penk = 20
  nboot=1000
  beta_median <- NULL
  se_median <- NULL
  P_median <- NULL
  
  #simple
  weights=rep(1/length(beta1), length(beta1))
  b <- weighted_median(beta1,weights )
  se <- weighted_median_bootstrap(beta1, sigma_b1,  weights,nboot)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail = FALSE)
  beta_median <- c(beta_median,b)
  se_median <- c(se_median,se)
  P_median <- c(P_median,pval)
  
  #weighted
  weights=1/(sigma_b1^2)
  b <- weighted_median(beta1, weights)
  se <- weighted_median_bootstrap(beta1, sigma_b1, weights, nboot)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail = FALSE)
  beta_median <- c(beta_median,b)
  se_median <- c(se_median,se)
  P_median <- c(P_median,pval)
  
  #penalised
  weights=1/(sigma_b1^2)
  bwm <- weighted_median(beta1, weights)
  penalty <- stats::pchisq(weights * (beta1 - bwm)^2, df = 1, 
                           lower.tail = FALSE)
  pen.weights <- weights * pmin(1, penalty * penk)
  b <- weighted_median(beta1, pen.weights)
  se <- weighted_median_bootstrap(beta1, sigma_b1,pen.weights,nboot)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail = FALSE)
  beta_median <- c(beta_median,b)
  se_median <- c(se_median,se)
  P_median <- c(P_median,pval)
  
  Method <- c("Simple median", "Weighted median", "Penalised median")
  Results <- data.frame(method = Method,b = beta_median, 
                        se = se_median, pval = P_median, stringsAsFactors = FALSE)
  
  return(Results)
}

#############Mode-based method#############
MR_MODE <- function(beta1,sigma_b){
  
  phi = 1
  penk = 20
  nboot=1000
  alpha=0.05
  # b_exp <- dat$beta.exposure
  # b_out <- dat$beta.outcome
  # se_exp <- dat$se.exposure
  # se_out <- dat$se.outcome
  beta <- function(BetaIV.in, seBetaIV.in, phi) {
    s <- 0.9 * (min(stats::sd(BetaIV.in), stats::mad(BetaIV.in)))/length(BetaIV.in)^(1/5)
    weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)
    beta <- NULL
    for (cur_phi in phi) {
      h <- max(1e-08, s * cur_phi)
      densityIV <- stats::density(BetaIV.in, weights = weights, 
                                  bw = h)
      beta[length(beta) + 1] <- densityIV$x[densityIV$y == 
                                              max(densityIV$y)]
    }
    return(beta)
  }
  boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in, nboot) {
    beta.boot <- matrix(nrow = nboot, ncol = length(beta_Mode.in))
    for (i in 1:nboot) {
      BetaIV.boot <- stats::rnorm(length(BetaIV.in), mean = BetaIV.in, 
                                  sd = seBetaIV.in)
      # BetaIV.boot_NOME <- stats::rnorm(length(BetaIV.in), 
      #                                  mean = BetaIV.in, sd = seBetaIV.in[, 2])
      beta.boot[i, 1:length(phi)] <- beta(BetaIV.in = BetaIV.boot, 
                                          seBetaIV.in = rep(1, length(BetaIV)), phi = phi)
      beta.boot[i, (length(phi) + 1):(2 * length(phi))] <- beta(BetaIV.in = BetaIV.boot, 
                                                                seBetaIV.in = seBetaIV.in, phi = phi)
      weights <- 1/seBetaIV.in^2
      penalty <- stats::pchisq(weights * (BetaIV.boot - 
                                            beta.boot[i, (length(phi) + 1):(2 * length(phi))])^2, 
                               df = 1, lower.tail = FALSE)
      pen.weights <- weights * pmin(1, penalty * penk)
      beta.boot[i, (2 * length(phi) + 1):(3 * length(phi))] <- beta(BetaIV.in = BetaIV.boot, 
                                                                    seBetaIV.in = sqrt(1/pen.weights), phi = phi)
      # beta.boot[i, (3 * length(phi) + 1):(4 * length(phi))] <- beta(BetaIV.in = BetaIV.boot_NOME, 
      #                                                               seBetaIV.in = rep(1, length(BetaIV)), phi = phi)
      # beta.boot[i, (4 * length(phi) + 1):(5 * length(phi))] <- beta(BetaIV.in = BetaIV.boot_NOME, 
      #                                                               seBetaIV.in = seBetaIV.in[, 2], phi = phi)
    }
    return(beta.boot)
  }
  # phi <- parameters$phi
  # nboot <- parameters$nboot
  # alpha <- parameters$alpha
  # BetaIV <- b_out/b_exp
  # seBetaIV <- cbind(sqrt((se_out^2)/(b_exp^2) + ((b_out^2) * 
  #                                                  (se_exp^2))/(b_exp^4)), se_out/abs(b_exp))
  BetaIV <- beta1
  seBetaIV <- sigma_b
  beta_SimpleMode <- beta(BetaIV.in = BetaIV, seBetaIV.in = rep(1,length(BetaIV)), phi = phi)
  beta_WeightedMode <- beta(BetaIV.in = BetaIV, seBetaIV.in = seBetaIV, phi = phi)
  weights <- 1/seBetaIV^2
  penalty <- stats::pchisq(weights * (BetaIV - beta_WeightedMode)^2, 
                           df = 1, lower.tail = FALSE)
  pen.weights <- weights * pmin(1, penalty * penk)
  beta_PenalisedMode <- beta(BetaIV.in = BetaIV, seBetaIV.in = sqrt(1/pen.weights), 
                             phi = phi)
  #beta_WeightedMode_NOME <- beta(BetaIV.in = BetaIV, seBetaIV.in = seBetaIV[,2], phi = phi)
  beta_Mode <- rep(c(beta_SimpleMode, beta_WeightedMode, beta_PenalisedMode))
  beta_Mode.boot <- boot(BetaIV.in = BetaIV, seBetaIV.in = seBetaIV, 
                         beta_Mode.in = beta_Mode, nboot = nboot)
  se_Mode <- apply(beta_Mode.boot, 2, stats::mad)
  CIlow_Mode <- beta_Mode - stats::qnorm(1 - alpha/2) * se_Mode
  CIupp_Mode <- beta_Mode + stats::qnorm(1 - alpha/2) * se_Mode
  P_Mode <- stats::pt(abs(beta_Mode/se_Mode), df = length(beta1) - 
                        1, lower.tail = FALSE) * 2
  Method <- rep(c("Simple mode", "Weighted mode", "Penalised mode"), each = length(phi))
  Results <- data.frame(method = Method,b = beta_Mode, 
                        se = se_Mode,  pval = P_Mode, stringsAsFactors = FALSE)
  
  return(Results)
  
}

######COMP###########
Comp <- function(x,N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                 sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b,
                 b1,b2,rho_b,siga1,siga2){
  cat('x=',x,'\n')
  sigma_b <- Comp_wald(N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                       sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b)
  sigma_b1 <- sigma_b[[1]]
  sigma_b2 <- sigma_b[[2]]
  rm(sigma_b)
  
  fdata <- DataGeneration(gb,b1,b2,rho_b,sigma_b1,sigma_b2,siga1,siga2)
  
  beta1 = fdata$beta1
  beta2 = fdata$beta2
  gamma1 = fdata$gamma1
  gamma2 = fdata$gamma2
  obeta1 <- beta1-gamma1
  
  res_ivw <- Trad_IVW(obeta1,sigma_b1)
  
  res_mle <- MLikelihood(beta1,beta2,b2,sigma_b1,sigma_b2,rho_b,gamma1,gamma2)
  
  res_median <- MR_MEDIAN(obeta1,sigma_b1)
  res_MBE <- MR_MODE(obeta1,sigma_b1)
  
  res_all <- c(res_ivw,res_mle,
               unlist(c(res_median[,-1])),unlist(c(res_MBE[,-1])),
               N1b,N2b,N3b,N4b,gb,b1,b2,rho_b,siga1,siga2)
  res_all <- matrix(res_all ,nrow=1)
  
  return(res_all)
}

#############MAIN############
main <- function(NN,x,N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                 sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b,
                 b1,b2,rho_b,siga1,siga2,mc){
  
  registerDoMC(mc)
  
  tt1 <- foreach(x=1:NN,
                 .combine=rbind,
                 .packages = c("MASS")) %dopar% {
                   Comp(x,N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                        sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b,
                        b1,b2,rho_b,siga1,siga2)
                 }
  return(tt1)
}


##########SIMULATION##########
ccname <- c('IVW_beta1','IVW_var1','IVW_pval1',
            
            'MRTP_beta1','MRTP_pvalue_beta1',
            
            'simple_median_beta1','weighted_median_beta1','penalised_median_beta1',
            'simple_median_se1','weighted_median_se1','penalised_median_se1',
            'simple_median_pval1','weighted_median_pval1','penalised_median_pval1',
            
            'simple_mode_beta1','weighted_mode_beta1','penalised_mode_beta1',
            'simple_mode_se1','weighted_mode_se1','penalised_mode_se1',
            'simple_mode_pval1','weighted_mode_pval1','penalised_mode_pval1',
            
            'N1b','N2b','N3b','N4b','g','b1','b2','rho_b',"siga1","siga2")

NN=1000
gb=100
mc=25
N1b=N3b=3000
N2bb=c(5000)
N4bb=c(500000)
p1b=0.2
p2b=0.3
beta1b=beta2b=0.5
sigma_x1b=0.15
sigma_x2b=0.09
sigma2_xi1b=0.15
sigma2_xi2b=0.09
sigma2_epi1b=0.1
sigma2_epi2b=0.09

bb <- c(0,0.1)
rho_bb <- c(0.2,0.8)
siga1b <- c(0,0)
siga2b <- c(0,0.2)


#m=n=t=p=1
result_mean <- NULL
for (q in 1:length(N2bb)){
  for(n in 1:length(N4bb)){
    for(m in 1:length(bb)){
      for(t in 1:length(rho_bb)){
        for(p in 1:length(siga1b)){
          N2b= N2bb[q]
          N4b= N4bb[n]
          b1=bb[m]
          b2=bb[m]
          rho_b=rho_bb[t]
          siga1=siga1b[p]
          siga2=siga2b[p]
          result <- main(NN,x,N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                         sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b,
                         b1,b2,rho_b,siga1,siga2,mc)
          #result <- Comp(x,g,b1,b2,rho_b,sigma_b1,sigma_b2)
          result <- as.data.frame(result)
          colnames(result) <- ccname
          write.csv(result,paste0("MC","_N2_",N2b,"_N4_",N4b,"_g_",gb,
                                  "_beta1_",b1,"_beta2_",b2,"_rho_b_",rho_b,"_siga1_",siga1,"_siga2_",siga2,
                                  "_0120.csv"))
        }
      }
    }
  }
}

NN=1000
gb=100
mc=25
N1b=N3b=3000
N2b=300000
N4b=300000
p1b=0.2
p2b=0.3
beta1b=beta2b=0.5
sigma_x1b=0.15
sigma_x2b=0.09
sigma2_xi1b=0.15
sigma2_xi2b=0.09
sigma2_epi1b=0.1
sigma2_epi2b=0.09

bb1 <- c(0,0.05,0.1,0.15,0.2)
bb2 <- c(0,0.05,0.1,0.15,0.2)
rho_bb <- c(0.1,0.2,0.4,0.6,0.8,0.9)
siga1b <- c(0,0,-0.2)
siga2b <- c(0,0.2,0.2)


#m=n=t=p=z=1
result_mean <- NULL
for(n in 1:length(bb1)){
  for(m in 1:length(bb2)){
    for(t in 1:length(rho_bb)){
      for(p in 1:length(siga1b)){
        b1=bb1[n]
        b2=bb2[m]
        rho_b=rho_bb[t]
        siga1=siga1b[p]
        siga2=siga2b[p]
        result <- main(NN,x,N1b,N2b,N3b,N4b,gb,p1b,p2b,beta1b,beta2b,sigma_x1b,sigma_x2b,
                       sigma2_xi1b,sigma2_xi2b,sigma2_epi1b,sigma2_epi2b,
                       b1,b2,rho_b,siga1,siga2,mc)
        #result <- Comp(x,g,b1,b2,rho_b,sigma_b1,sigma_b2)
        result <- as.data.frame(result)
        colnames(result) <- ccname
        write.csv(result,paste0("MC","_N2_",N2b,"_N4_",N4b,"_g_",gb,
                                "_beta1_",b1,"_beta2_",b2,"_rho_b_",rho_b,"_siga1_",siga1,"_siga2_",siga2,
                                "_0120.csv"))
      }
    }
  }
}
