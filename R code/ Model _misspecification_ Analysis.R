

# True Data Generating Model: Burr XII
# Fitted (Misspecified) Model: Weibull
rm(list=ls())
# Fixed Parameters & Design Setup
 
beta0 <- -0.2
alpha <- 2.5
beta1 <- 0.5

K <- 3 # Number of stress levels
I <- 3 # Number of inspection intervals
x <- matrix(c(
  0.20, 0.55, 1.00,
  0.25, 0.60, 1.00,
  0.30, 0.65, 1.00
), nrow = K, byrow = TRUE)

tau <- matrix(c(
  0.5, 1.0, 1.5,
  1.0, 1.5, 2.0,
  2, 2.5, 3
), nrow = K, byrow = TRUE)

R.reps <- 1000   

q_vec <- c(0.10, 0.50, 0.95)
Nk_vec <- c(50, 100)

# TRUE MODEL: Burr XII Log-likelihood
loglik_burr <- function(par, n, N, x, tau){  
  b0 <- par[1]; a <- par[2]; b1 <- par[3]  
  ll <- 0  
  for(k in 1:nrow(x)){  
    Hcum <- 0; Rprev <- 1; nk <- 0  
    for(i in 1:ncol(x)){  
      lambda <- exp(b0 + b1 * x[k, i])  
      if(i == 1){ Hcum <- lambda * log(1 + tau[k, i]^a) } 
      else { Hcum <- Hcum + lambda * (log(1 + tau[k, i]^a) - log(1 + tau[k, i-1]^a)) }  
      
      Rcurr <- exp(-Hcum)  
      D <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)  
      ll <- ll + n[k, i] * log(D)  
      nk <- nk + n[k, i]
      Rprev <- Rcurr  
    }  
    Rsurv <- pmin(pmax(Rprev, 1e-10), 1 - 1e-10)  
    ll <- ll + (N[k] - nk) * log(Rsurv)  
  }  
  return(-ll)  
}

# MISSPECIFIED MODEL: Weibull Log-likelihood
loglik_weibull <- function(par, n, N, x, tau){  
  b0 <- par[1]; a <- par[2]; b1 <- par[3]  
  ll <- 0  
  for(k in 1:nrow(x)){  
    Hcum <- 0; Rprev <- 1; nk <- 0  
    for(i in 1:ncol(x)){  
      lambda <- exp(b0 + b1 * x[k, i])  
      if(i == 1){ Hcum <- lambda * tau[k, i]^a } 
      else { Hcum <- Hcum + lambda * (tau[k, i]^a - tau[k, i-1]^a) }  
      
      Rcurr <- exp(-Hcum)  
      D <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)  
      ll <- ll + n[k, i] * log(D)  
      nk <- nk + n[k, i]
      Rprev <- Rcurr  
    }  
    Rsurv <- pmin(pmax(Rprev, 1e-10), 1 - 1e-10)  
    ll <- ll + (N[k] - nk) * log(Rsurv)  
  }  
  return(-ll)  
}

Results_List <- list()

# Initial values for optimization 
start_B <- c(-1.0, 2.5, 0.4) 
start_W <- c(-1.0, 2.0, 0.5) 

cat("Starting Simulations...\n")

for(Nk_val in Nk_vec) {
  for(q_val in q_vec) {
    
    cat(sprintf("Processing Nk = %d, q = %.2f...\n", Nk_val, q_val))
    
    # Set current design parameters
    q <- q_val
    Nk <- rep(Nk_val, K)
    
    # True Theoretical Quantile (Burr XII Baseline)
    Q.true <- ( (1 - q)^(-exp(-beta0)) - 1 )^(1/alpha)
    
    Q.hat.true  <- rep(NA, R.reps) 
    Q.hat.wrong <- rep(NA, R.reps) 
    valid <- rep(FALSE, R.reps)
    
    for(r in 1:R.reps){  
      
      #  Generate Data based on the True Model (Burr XII)
      failure <- matrix(0, K, I)  
      for(k in 1:K){  
        Hcum <- 0; Rprev <- 1; D <- numeric(I)  
        for(i in 1:I){  
          lambda <- exp(beta0 + beta1 * x[k, i])  
          
          if(i == 1){ Hcum <- lambda * log(1 + tau[k, i]^alpha) } 
          else { Hcum <- Hcum + lambda * (log(1 + tau[k, i]^alpha) - log(1 + tau[k, i-1]^alpha)) }  
          
          Rcurr <- exp(-Hcum)  
          D[i] <- Rprev - Rcurr  
          Rprev <- Rcurr  
        }  
        probs <- c(D, Rprev)  
        failure[k, ] <- rmultinom(1, Nk[k], probs)[1:I]  
      }  
      
      # Fit the True Model (Burr XII)
      fitB <- try(optim(start_B, loglik_burr, n = failure, N = Nk, x = x, tau = tau, 
                        method = "L-BFGS-B", lower = c(-10, 0.01, -10), upper = c(10, 10, 10)), silent = TRUE)  
      
      #  Fit the Misspecified Model (Weibull)
      fitW <- try(optim(start_W, loglik_weibull, n = failure, N = Nk, x = x, tau = tau, 
                        method = "L-BFGS-B", lower = c(-10, 0.01, -10), upper = c(10, 10, 10)), silent = TRUE)  
      
      # Convergence Check
      if(inherits(fitB, "try-error") || inherits(fitW, "try-error")) next
      if(fitB$convergence != 0 || fitW$convergence != 0) next
      
      valid[r] <- TRUE
      
      #Calculate Estimated Baseline Quantiles for this iteration
      Q.hat.true[r] <- ( (1 - q)^(-exp(-fitB$par[1])) - 1 )^(1/fitB$par[2])
      Q.hat.wrong[r] <- exp(-fitW$par[1]/fitW$par[2]) * (-log(1-q))^(1/fitW$par[2])
    }
    
    valid_idx <- which(valid)
    Q_true_est  <- Q.hat.true[valid_idx]
    Q_wrong_est <- Q.hat.wrong[valid_idx]
    
    # Averages and MSE
    mean_Q_wrong <- mean(Q_wrong_est)
    mean_Q_true  <- mean(Q_true_est)
    MSE_wrong    <- mean((Q_wrong_est - Q.true)^2)
    
    # RB and RV
    RB <- (mean_Q_wrong - mean_Q_true) / mean_Q_true
    num_RV <- sum((Q_wrong_est - mean_Q_true)^2)
    den_RV <- sum((Q_true_est - mean_Q_true)^2)
    RV <- num_RV / den_RV
    
    # results
    Results_List[[length(Results_List) + 1]] <- data.frame(
      Nk = Nk_val,
      q = q_val,
      Valid_Reps = length(valid_idx),
      True_Quantile = Q.true,
      Mean_Est_Wrong = mean_Q_wrong,
      MSE_Wrong = MSE_wrong,
      RB = RB,
      RV = RV
    )
  }
} 
cat("\n############################################################\n")
cat("  MISSPECIFICATION ANALYSIS SUMMARY\n")
cat(" True Model: Burr XII | Fitted Model: Weibull\n")
cat("############################################################\n\n")

Final_Results <- do.call(rbind, Results_List)
Final_Results[, 4:8] <- round(Final_Results[, 4:8], 5)

print(Final_Results, row.names = FALSE)


#########################################

rm(list=ls())
# True Data Generating Model: Weibull
# Fitted (Misspecified) Model: Burr XII

beta0 <- -0.2
alpha <- 2.5
beta1 <- 0.5

K <- 3 # Number of stress levels
I <- 3 # Number of inspecti
x <- matrix(c(
  0.20, 0.55, 1.00,
  0.25, 0.60, 1.00,
  0.30, 0.65, 1.00
), nrow = K, byrow = TRUE)

tau <- matrix(c(
  0.5, 1.0, 1.5,
  1.0, 1.5, 2.0,
  2, 2.5, 3
), nrow = K, byrow = TRUE)

R.reps <- 1000   

q_vec <- c(0.10, 0.50, 0.95)
Nk_vec <- c(50, 100)

# MISSPECIFIED MODEL: Burr XII Log-likelihood
loglik_burr <- function(par, n, N, x, tau){  
  b0 <- par[1]; a <- par[2]; b1 <- par[3]  
  ll <- 0  
  for(k in 1:nrow(x)){  
    Hcum <- 0; Rprev <- 1; nk <- 0  
    for(i in 1:ncol(x)){  
      lambda <- exp(b0 + b1 * x[k, i])  
      if(i == 1){ Hcum <- lambda * log(1 + tau[k, i]^a) } 
      else { Hcum <- Hcum + lambda * (log(1 + tau[k, i]^a) - log(1 + tau[k, i-1]^a)) }  
      
      Rcurr <- exp(-Hcum)  
      D <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)  
      ll <- ll + n[k, i] * log(D)  
      nk <- nk + n[k, i]
      Rprev <- Rcurr  
    }  
    Rsurv <- pmin(pmax(Rprev, 1e-10), 1 - 1e-10)  
    ll <- ll + (N[k] - nk) * log(Rsurv)  
  }  
  return(-ll)  
}

# TRUE MODEL: Weibull Log-likelihood
loglik_weibull <- function(par, n, N, x, tau){  
  b0 <- par[1]; a <- par[2]; b1 <- par[3]  
  ll <- 0  
  for(k in 1:nrow(x)){  
    Hcum <- 0; Rprev <- 1; nk <- 0  
    for(i in 1:ncol(x)){  
      lambda <- exp(b0 + b1 * x[k, i])  
      if(i == 1){ Hcum <- lambda * tau[k, i]^a } 
      else { Hcum <- Hcum + lambda * (tau[k, i]^a - tau[k, i-1]^a) }  
      
      Rcurr <- exp(-Hcum)  
      D <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)  
      ll <- ll + n[k, i] * log(D)  
      nk <- nk + n[k, i]
      Rprev <- Rcurr  
    }  
    Rsurv <- pmin(pmax(Rprev, 1e-10), 1 - 1e-10)  
    ll <- ll + (N[k] - nk) * log(Rsurv)  
  }  
  return(-ll)  
}

Results_List <- list()

start_B <- c(-1.0, 2.5, 0.4) 
start_W <- c(-1.0, 2.0, 0.5) 

cat("Starting Simulations...\n")

for(Nk_val in Nk_vec) {
  for(q_val in q_vec) {
    
    cat(sprintf("Processing Nk = %d, q = %.2f...\n", Nk_val, q_val))
    
    q <- q_val
    Nk <- rep(Nk_val, K)
    
    # True Theoretical Quantile (Weibull Baseline)
    Q.true <- exp(-beta0/alpha) * (-log(1-q))^(1/alpha)
    
    Q.hat.true  <- rep(NA, R.reps) 
    Q.hat.wrong <- rep(NA, R.reps) 
    valid <- rep(FALSE, R.reps)
    
    for(r in 1:R.reps){  
      
      # Generate Data based on the True Model (Weibull)
      failure <- matrix(0, K, I)  
      for(k in 1:K){  
        Hcum <- 0; Rprev <- 1; D <- numeric(I)  
        for(i in 1:I){  
          lambda <- exp(beta0 + beta1 * x[k, i])  
          
          if(i == 1){ Hcum <- lambda * tau[k, i]^alpha } 
          else { Hcum <- Hcum + lambda * (tau[k, i]^alpha - tau[k, i-1]^alpha) }  
          
          Rcurr <- exp(-Hcum)  
          D[i] <- Rprev - Rcurr  
          Rprev <- Rcurr  
        }  
        probs <- c(D, Rprev)  
        failure[k, ] <- rmultinom(1, Nk[k], probs)[1:I]  
      }  
      
      #  Fit the Misspecified Model (Burr XII)
      fitB <- try(optim(start_B, loglik_burr, n = failure, N = Nk, x = x, tau = tau, 
                        method = "L-BFGS-B", lower = c(-10, 0.01, -10), upper = c(10, 5, 10)), silent = TRUE)  
      
      #  Fit the True Model (Weibull)
      fitW <- try(optim(start_W, loglik_weibull, n = failure, N = Nk, x = x, tau = tau, 
                        method = "L-BFGS-B", lower = c(-10, 0.01, -10), upper = c(10, 5, 10)), silent = TRUE)  
      
      if(inherits(fitB, "try-error") || inherits(fitW, "try-error")) next
      if(fitB$convergence != 0 || fitW$convergence != 0) next
      
      valid[r] <- TRUE
      
    
      Q.hat.true[r] <- exp(-fitW$par[1]/fitW$par[2]) * (-log(1-q))^(1/fitW$par[2])
      
      Q.hat.wrong[r] <- ( (1 - q)^(-exp(-fitB$par[1])) - 1 )^(1/fitB$par[2])
    }
   
    valid_idx <- which(valid)
    Q_true_est  <- Q.hat.true[valid_idx]
    Q_wrong_est <- Q.hat.wrong[valid_idx]
    
    # Averages and MSE
    mean_Q_wrong <- mean(Q_wrong_est)
    mean_Q_true  <- mean(Q_true_est)
    MSE_wrong    <- mean((Q_wrong_est - Q.true)^2)
    
    # RB and RV
    RB <- (mean_Q_wrong - mean_Q_true) / mean_Q_true
    num_RV <- sum((Q_wrong_est - mean_Q_true)^2)
    den_RV <- sum((Q_true_est - mean_Q_true)^2)
    RV <- num_RV / den_RV
    
    Results_List[[length(Results_List) + 1]] <- data.frame(
      Nk = Nk_val,
      q = q_val,
      Valid_Reps = length(valid_idx),
      True_Quantile = Q.true,
      Mean_Est_Wrong = mean_Q_wrong,
      MSE_Wrong = MSE_wrong,
      RB = RB,
      RV = RV
    )
  }
} 
cat("\n############################################################\n")
cat(" MISSPECIFICATION ANALYSIS SUMMARY\n")
cat(" True Model: Weibull | Fitted Model: Burr XII\n")
cat("############################################################\n\n")

Final_Results <- do.call(rbind, Results_List)

Final_Results[, 4:8] <- round(Final_Results[, 4:8], 5)

print(Final_Results, row.names = FALSE)

#########################################################################

# 1. Weibull (True) -> Chen (Misspecified)
# 2. Chen (True) -> Weibull (Misspecified)
#########################################################################

rm(list=ls())  
K <- 3         # Number of stress levels
I <- 3         # Number of inspection intervals
R.reps <- 1000 

x <- matrix(c(
  0.20, 0.55, 1.00,
  0.25, 0.60, 1.00,
  0.30, 0.65, 1.00
), nrow = K, byrow = TRUE)

tau <- matrix(c(
  0.5, 1.0, 1.5,
  1.0, 1.5, 2.0,
  2, 2.5, 3
), nrow = K, byrow = TRUE)

q_vec <- c(0.10, 0.50, 0.95)
Nk_vec <- c(50, 100)

# Weibull Log-likelihood
loglik_weibull <- function(par, n, N, x, tau) {  
  b0 <- par[1]; a <- par[2]; b1 <- par[3]  
  ll <- 0  
  for(k in 1:nrow(x)){  
    Hcum <- 0; Rprev <- 1; nk <- 0  
    for(i in 1:ncol(x)){  
      lambda <- exp(b0 + b1 * x[k, i])  
      deltaH <- if(i == 1) tau[k, i]^a else (tau[k, i]^a - tau[k, i-1]^a)
      Hcum <- Hcum + lambda * deltaH
      Rcurr <- exp(-Hcum)  
      D <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)  
      ll <- ll + n[k, i] * log(D)  
      nk <- nk + n[k, i]; Rprev <- Rcurr  
    }  
    Rsurv <- pmin(pmax(Rprev, 1e-10), 1 - 1e-10)  
    ll <- ll + (N[k] - nk) * log(Rsurv)  
  }  
  return(-ll)  
}

# Chen Log-likelihood
loglik_chen <- function(par, n, N, x, tau) {
  b0 <- par[1]; a <- par[2]; b1 <- par[3]
  ll <- 0
  for(k in 1:nrow(x)){
    Hcum <- 0; Rprev <- 1; nk <- 0
    for(i in 1:ncol(x)){
      lambda <- exp(b0 + b1 * x[k, i])
      deltaH <- if(i == 1) (exp(tau[k, i]^a) - 1) else (exp(tau[k, i]^a) - exp(tau[k, i-1]^a))
      Hcum <- Hcum + lambda * deltaH
      Rcurr <- exp(-Hcum)
      D <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)
      ll <- ll + n[k, i] * log(D)
      nk <- nk + n[k, i]; Rprev <- Rcurr
    }
    Rsurv <- pmin(pmax(Rprev, 1e-10), 1 - 1e-10)
    ll <- ll + (N[k] - nk) * log(Rsurv)
  }
  return(-ll)
}

# Data Generator
generate_data <- function(model, par, K, I, x, tau, Nk) {
  b0 <- par[1]; a <- par[2]; b1 <- par[3]
  failure <- matrix(0, K, I)
  for(k in 1:K){
    Hcum <- 0; Rprev <- 1; D <- numeric(I)
    for(i in 1:I){
      lambda <- exp(b0 + b1 * x[k, i])
      if(model == "Weibull") {
        deltaH <- if(i==1) tau[k, i]^a else (tau[k, i]^a - tau[k, i-1]^a)
      } else if (model == "Chen") {
        deltaH <- if(i==1) (exp(tau[k, i]^a) - 1) else (exp(tau[k, i]^a) - exp(tau[k, i-1]^a))
      }
      Hcum <- Hcum + lambda * deltaH
      Rcurr <- exp(-Hcum)
      D[i] <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)
      Rprev <- Rcurr
    }
    probs <- c(D, Rprev)
    failure[k, ] <- rmultinom(1, Nk[k], probs)[1:I]
  }
  return(failure)
}

# Baseline Quantile 
calc_Q <- function(model, b0, a, q) {
  if(model == "Weibull") {
    return(exp(-b0/a) * (-log(1-q))^(1/a))
  } else if (model == "Chen") {
    return( ( log(1 - log(1-q)*exp(-b0)) )^(1/a) )
  }
}

get_loglik <- function(model) {
  if(model == "Weibull") return(loglik_weibull)
  if(model == "Chen") return(loglik_chen)
}

# Optimization Starting Guesses
get_start_par <- function(model) {
  if(model == "Weibull") return(c(-1.0, 2.0, 0.5))
  if(model == "Chen") return(c(-2.0, 1.6, 0.5))
}
scenarios <- list(
  list(ID = 1, TrueModel = "Weibull",  WrongModel = "Chen",     TruePar = c(-0.2, 2.0, 0.5)),
  list(ID = 2, TrueModel = "Chen",     WrongModel = "Weibull",  TruePar = c(-0.2, 2.0, 0.5))
)

Results_Table <- list()

cat("Running Model Misspecification Analysis (Weibull & Chen)...\n\n")

for(sc in scenarios) {
  
  cat(sprintf("Processing Scenario %d: True = %s | Fitted (Wrong) = %s\n", 
              sc$ID, sc$TrueModel, sc$WrongModel))
  
  # Set up functions and parameters for current scenario
  ll_true  <- get_loglik(sc$TrueModel)
  ll_wrong <- get_loglik(sc$WrongModel)
  
  start_true  <- get_start_par(sc$TrueModel)
  start_wrong <- get_start_par(sc$WrongModel)
  
  for(Nk_val in Nk_vec) {
    Nk <- rep(Nk_val, K)
    
    for(q_val in q_vec) {
      q <- q_val
      
      cat(sprintf("  -> Nk = %d, q = %.2f...\n", Nk_val, q_val))
      
      Q.true <- calc_Q(sc$TrueModel, sc$TruePar[1], sc$TruePar[2], q)
      
      Q.hat.true  <- rep(NA, R.reps) 
      Q.hat.wrong <- rep(NA, R.reps) 
      valid <- rep(FALSE, R.reps)
      
      for(r in 1:R.reps){  
        
        # Generate Data from True Model
        failure <- generate_data(sc$TrueModel, sc$TruePar, K, I, x, tau, Nk)
        
        # Fit True Model
        fitT <- try(optim(start_true, ll_true, n = failure, N = Nk, x = x, tau = tau, 
                          method = "L-BFGS-B", lower = c(-10, 1e-3, -10), upper = c(10, 10, 10)), silent = TRUE)  
        
        # Fit Misspecified Model
        fitW <- try(optim(start_wrong, ll_wrong, n = failure, N = Nk, x = x, tau = tau, 
                          method = "L-BFGS-B", lower = c(-10, 1e-3, -10), upper = c(10, 10, 10)), silent = TRUE)  
        
       
        if(inherits(fitT, "try-error") || inherits(fitW, "try-error")) next
        if(fitT$convergence != 0 || fitW$convergence != 0) next
        
        valid[r] <- TRUE
        
        # Compute baseline quantiles
        Q.hat.true[r]  <- calc_Q(sc$TrueModel,  fitT$par[1], fitT$par[2], q)
        Q.hat.wrong[r] <- calc_Q(sc$WrongModel, fitW$par[1], fitW$par[2], q)
      }
      
      valid_idx <- which(valid)
      
      if(length(valid_idx) == 0) {
        cat("     [!] ERROR: No models converged for this iteration.\n")
        next
      }
      
      Q_true_est  <- Q.hat.true[valid_idx]
      Q_wrong_est <- Q.hat.wrong[valid_idx]
      
      mean_Q_wrong <- mean(Q_wrong_est)
      mean_Q_true  <- mean(Q_true_est)
      MSE_wrong    <- mean((Q_wrong_est - Q.true)^2)
      
      # RB & RV
      RB <- (mean_Q_wrong - mean_Q_true) / mean_Q_true
      RV <- sum((Q_wrong_est - mean_Q_true)^2) / sum((Q_true_est - mean_Q_true)^2)
      
      Results_Table[[length(Results_Table) + 1]] <- data.frame(
        Scenario     = paste0(sc$TrueModel, " -> ", sc$WrongModel),
        Nk           = Nk_val,
        q            = q_val,
        Valid_Reps   = length(valid_idx),
        True_Quantile= Q.true,
        Est_True     = mean_Q_true,
        Est_Wrong    = mean_Q_wrong,
        MSE          = MSE_wrong,
        RB           = RB,
        RV           = RV
      )
    }
  }
}
cat("\n###################################################################################\n")
cat(" MISSPECIFICATION ANALYSIS SUMMARY (WEIBULL AND CHEN)\n")
cat("###################################################################################\n")

Final_Summary <- do.call(rbind, Results_Table)

# Format numerical columns for cleaner terminal output
Final_Summary[, 5:10] <- round(Final_Summary[, 5:10], 5)

print(Final_Summary, row.names = FALSE)
