
library(MASS)
library(glmnet)
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

DataGenerator_wald <- function(N, g, p, beta, sigma_x) {
  
  G <- list()
  X <- list()
  Y <- list()
  U <- list()
  
  sigma_alpha <- diag(sigma_x)
  
  alpha <-mvrnorm(n = g, mu = rep(0, length(sigma_x)), Sigma = sigma_alpha)
  
  for(i in 1:length(N)) {
    
    G[[i]] <- matrix(rbinom(N[i]*g, 2, p[i]), ncol = g)
    
    U[[i]] <- rnorm(N[i],0,1)
    
    X[[i]] <- G[[i]] %*% alpha[,i] + rnorm(N[i], 0, 1) + U[[i]]
    
    Y[[i]] <- beta[i] * X[[i]] + rnorm(N[i], 0, 1) + U[[i]]
  }
  
  datasimul <- c(G, X, Y)
  names(datasimul) <- c(paste0("G", 1:length(N)), paste0("X", 1:length(N)), 
                        paste0("Y", 1:length(N)))
  
  return(datasimul)
}

Comp_wald <- function(Nb, gb, pb, betab, sigma_xb) {

  fdata <- DataGenerator_wald(Nb, gb, pb, betab, sigma_xb)
  
  results <- list()
  
  for(race in 1:length(Nb)) {
    betaXG <- list()
    betaYG <- list()
    seXG <- list()
    seYG <- list()
    
    for(j in 1:gb) {

      FF_X <- lm(fdata[[paste0("X", race)]] ~ fdata[[paste0("G", race)]][,j])
      betaXG[[j]] <- FF_X$coef[2]
      seXG[[j]] <- summary(FF_X)$coef[2,2]

      FF_Y <- lm(fdata[[paste0("Y", race)]] ~ fdata[[paste0("G", race)]][,j])
      betaYG[[j]] <- FF_Y$coef[2]
      seYG[[j]] <- summary(FF_Y)$coef[2,2]
    }

    res_wald <- Wald_Ratio(unlist(betaXG),
                           unlist(betaYG),
                           unlist(seXG),
                           unlist(seYG))
    res_wald <- res_wald[order(res_wald$varWR_2),]
    sigma_b <- sqrt(res_wald$varWR_2)

    results[[paste0("sigma_b_race", race)]] <- sigma_b
  }
  
  return(results)
}

DataGeneration <- function(g, b, rhob, sigma_b, siga1, siga2, race) {
  gamma <- lapply(1:race, function(x) runif(g, siga1, siga2))
  
  beta <- matrix(0, nrow = g, ncol = race)
  sigma_matrices <- list() 
  
  rho_index <- 0
  rho_matrix <- matrix(0, nrow = race, ncol = race)
  for(r1 in 1:(race-1)){
    for(r2 in (r1+1):race){
      rho_index <- rho_index + 1
      rho_matrix[r1, r2] <- rhob[rho_index]
      rho_matrix[r2, r1] <- rhob[rho_index] 
    }
  }
  
  for(j in 1:g){
    Sigma <- matrix(0, nrow = race, ncol = race)
    for(r1 in 1:race){
      for(r2 in 1:race){
        if(r1 == r2){
          Sigma[r1, r2] <- sigma_b[[r1]][j]^2
        } else {
          Sigma[r1, r2] <- rho_matrix[r1, r2] * sigma_b[[r1]][j] * sigma_b[[r2]][j]
        }
      }
    }
    
    beta_j <- mvrnorm(1, mu = b, Sigma = Sigma)
    beta[j, ] <- beta_j
    sigma_matrices[[j]] <- Sigma 
  }
  
  for(i in 1:race){
    beta[, i] <- beta[, i] + gamma[[i]]
  }
  
  data_list <- list()
  for(i in 1:ncol(beta)) {
    data_list[[paste0("beta", i)]] <- beta[, i]
    data_list[[paste0("sebeta", i)]] <- sigma_b[[i]]
    data_list[[paste0("gamma", i)]] <- gamma[[i]]
  }
  
  data_sum <- as.data.frame(data_list)
  
  return(list(data_sum = data_sum, sigma_matrices = sigma_matrices))
}

#########IVW############
Trad_IVW <- function(beta,sigma_b){
  
  beta_IVW <- sum(beta*(1/(sigma_b^2)))/sum(1/(sigma_b^2))
  var_IVW <- 1/sum(1/(sigma_b^2))
  se_IVW <- sqrt(var_IVW)
  pvalue_IVW <- 2 * stats::pnorm(abs(beta_IVW/((var_IVW)^0.5)), lower.tail = FALSE)
  
  resi <- c(beta_IVW, se_IVW, pvalue_IVW)
  return(resi)
}

############### MRTP #############
MLikelihood <- function(beta, beta_know,known_b, gamma, gamma_know, 
                        sigma_b_all, rhob, race, gb){
  
  nll <- function(bbe1) {
    
    bbe <- c(bbe1, known_b[1], known_b[2])
    sigma_bb <- list()
    pp <- NULL
    det_sigma <- NULL
    
    for (i in 1:gb) {
      sigma_tt <- sigma_b_all[[i]][1:(race-1), 1:(race-1)]
      sigma_t1 <- sigma_b_all[[i]][1:(race-1), race]
      sigma_e2 <- 1/as.numeric(sigma_b_all[[i]][race, race])
      sigma_bb <- sigma_tt - sigma_e2 * sigma_t1 %*% t(sigma_t1)
      
      inv_sigma <- ginv(sigma_bb)
      diff <- beta[, i] - bbe - gamma[, i] - sigma_e2 * sigma_t1 * (beta_know[i] - known_b[3] - gamma_know[i])
      pp[i] <- as.numeric(t(diff) %*% inv_sigma %*% diff)
      det_sigma[i] <- det(sigma_bb)
    }
    
    sum(gb * log(2 * pi) + 0.5 * log(abs(det_sigma)) + 0.5 * pp)
  }
  
  tryFirstMethod <- function() {
    
    fit <- stats4::mle(minuslog=nll, 
                       start=list(bbe1=0),
                       method = "CG")
    res <- fit@coef
    
    fit1 <- stats4::mle(minuslog=nll, 
                        start=list(bbe1=0),
                        fixed = list(bbe1=0),
                        method = "CG" )
    
    stat_beta1=  2 * (fit@min - fit1@min)
    pvalue_beta1 = pchisq(-stat_beta1,1,lower.tail=F)
    
    return(c(res, pvalue_beta1))
  }

  result <- tryCatch(
    tryFirstMethod(),
    error = function(e) NULL  
  )
  
  if (is.null(result)) {
    cat("Both methods failed for the current iteration.\n")
    return(NULL) 
  }
  
  return(result)
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

    }
    return(beta.boot)
  }

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
Comp <- function(x, Nb, gb, pb, betab, sigma_xb, bb, rhob, siga1, siga2, race, maxRetries) {
  cat('x=', x, '\n')
  
  success <- FALSE
  attempt <- 0
  
  while (!success && attempt < maxRetries) {
    attempt <- attempt + 1
    nn <- race - 1
    
    sigma_b <- Comp_wald(Nb, gb, pb, betab, sigma_xb)
    fdata <- DataGeneration(gb, bb, rhob, sigma_b, siga1, siga2, race)
    
    beta <- t(fdata$data_sum[grep("^beta\\d+$", colnames(fdata$data_sum), value = TRUE)[-race]])
    gamma <- t(fdata$data_sum[grep("^gamma\\d+$", colnames(fdata$data_sum), value = TRUE)[-race]])
    beta_know <- fdata$data_sum[[paste0("beta", race)]]
    gamma_know <- fdata$data_sum[[paste0("gamma", race)]]
    sigma_b_all <- fdata$sigma_matrices
    known_b <- bb[2:race]
    
    res_mle <- MLikelihood(beta, beta_know, known_b, gamma, gamma_know, sigma_b_all, rhob, race, gb)
    
    if (!is.null(res_mle)) {
      success <- TRUE
      
      obetas <- fdata$data_sum[paste0("beta", 1:nn)] - fdata$data_sum[paste0("gamma", 1:nn)]
      res_ivw <- lapply(1:nn, function(i) Trad_IVW(obetas[[i]], sigma_b[[i]]))
      res_median <- lapply(1:nn, function(i) MR_MEDIAN(obetas[[i]], sigma_b[[i]]))
      res_MBE <- lapply(1:nn, function(i) MR_MODE(obetas[[i]], sigma_b[[i]]))
      
      extractNumericValues <- function(listResult) {
        allValues <- unlist(listResult)
        numericValues <- as.numeric(allValues)
        validValues <- numericValues[!is.na(numericValues)]
        return(validValues)
      }
      
      numericResMedian <- extractNumericValues(res_median)
      numericResMBE <- extractNumericValues(res_MBE)

      res_all <- c(
        matrix(unlist(res_ivw), nrow = 1),
        res_mle,
        numericResMedian,
        numericResMBE,
        bb, rhob, race, siga1, siga2
      )
      
      res_all <- matrix(res_all, nrow = 1)
    } else {
      cat("Both methods failed for the current iteration, retrying...\n")
    }
  }
  
  if (!success) {
    return(NULL) 
  }
  
  return(res_all)
}


#############MAIN############
main <- function(NN,x,Nb, gb, pb, betab, sigma_xb, bb, rhob, siga1, siga2, race, maxRetries, mc){
  
  registerDoMC(mc)
  
  tt1 <- foreach(x=1:NN,
                 .combine=rbind,
                 .errorhandling = "remove",
                 .packages = c("MASS","glmnet")) %dopar% {
                   Comp(x, Nb, gb, pb, betab, sigma_xb, bb, rhob, siga1, siga2, race, maxRetries)
                 }
  return(tt1)
}


##########SIMULATION##########

NN=1000
maxRetries=20
mc=50
race = 4
Nb <- c(rep(3000,(race-1)),300000) # N1-N4
gb <- 100
pb <- c(rep(0.2,(race-1)),0.3) # p1-p4
betab <- rep(0.5,race) # beta1-beta4
sigma_xb <- c(runif((race-1),0.12,0.15),0.09) # sigma_x1-sigma_x4

siga1b <- c(0,0,-0.2)
siga2b <- c(0,0.2,0.2)
bb_all <- list(rep(0.2,race),
               rep(0.15,race),
               rep(0.1,race),
               rep(0.05,race),
               rep(0,race)) 
rho_bb <- list(rep(0.1,sum(1:(race-1))),
               rep(0.2,sum(1:(race-1))),
               rep(0.4,sum(1:(race-1))),
               rep(0.6,sum(1:(race-1))),
               rep(0.8,sum(1:(race-1))),
               rep(0.9,sum(1:(race-1))),
               seq(0.15,0.4,0.05),
               seq(0.65,0.9,0.05)
)


n_beta <- race - 1
n_rho <- sum(1:(race - 1))

b_strings <- paste0("b", 1:race)
rho_strings <- paste0("rho", 1:n_rho)

base_names <- c('IVW',
                'MRTP',
                'simple_median', 'weighted_median', 'penalised_median',
                'simple_mode', 'weighted_mode', 'penalised_mode')

ccname <- unlist(sapply(base_names, function(base) {
  if (base == "MRTP") {
    c(paste0(base, "_beta", 1), paste0(base, "_pvalue_beta", 1))
  } else {
    c(
      paste0(base, "_beta", 1:n_beta),
      paste0(base, "_se", 1:n_beta),
      paste0(base, "_pval", 1:n_beta)
    )
  }
}), use.names = FALSE)

ivw_indices <- which(grepl("IVW", ccname))
ordered_ccname <- ccname
for (i in 1:n_beta) {
  indices <- c(i, i + n_beta, i + 2*n_beta)
  ordered_ccname[indices] <- c(
    ccname[((i-1)*n_beta + 1):((i-1)*n_beta + 3)]
  )
}

ccname <- ordered_ccname

ccname <- c(ccname, b_strings, rho_strings, 'race', "siga1", "siga2")


result_mean <- NULL
for(m in 1:length(bb_all)){
  for(t in 1:length(rho_bb)){
    for(p in 1:length(siga1b)){
      bb=bb_all[[m]]
      rhob=rho_bb[[t]]
      siga1=siga1b[p]
      siga2=siga2b[p]
      result <- main(NN, x, Nb, gb, pb, betab, 
                     sigma_xb, bb, rhob, siga1, siga2, race, maxRetries, mc)
      result <- as.data.frame(result)
      colnames(result) <- ccname
      write.csv(result,paste0("MC","_g_",gb,"_bb_",bb[1],"_rho_b_",rhob[1],
                              "_siga1_",siga1,"_siga2_",siga2,"_race_",race,
                              "_0222.csv"))
      
    }
  }
}  

