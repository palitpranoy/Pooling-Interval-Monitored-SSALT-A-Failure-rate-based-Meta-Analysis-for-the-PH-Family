

rm(list = ls())

# True parameters
beta0 <- -2
alpha <- 1
beta1 <- 3


# Design
K <- 3
I <- 3

x <- matrix(c(
  0.2, 0.55, 1,
  0.25, 0.60, 1,
  0.3, 0.65, 1
),
nrow = 3, byrow = TRUE)

tau <- matrix(c(
  0.5, 1.0, 1.5,
  1, 1.5, 2,
  2, 2.5, 3
),
nrow = 3, byrow = TRUE)
Nk <- rep(50, K)

# Baseline cumulative hazard
H0 <- function(t, alpha) {
  t^alpha
}

# Negative log-likelihood

loglik <- function(par, n, N, x, tau) {

  beta0 <- par[1]
  beta1 <- par[2]

  ll <- 0

  for(k in 1:nrow(x)) {

    Hcum <- 0
    Rprev <- 1
    nk <- 0

    for(i in 1:ncol(x)) {

      lambda <- exp(beta0 + beta1 * x[k, i])

      if(i == 1) {

        Hcum <- lambda * tau[k, i]

      } else {

        Hcum <- Hcum +
          lambda * (tau[k, i] - tau[k, i - 1])

      }

      Rcurr <- exp(-Hcum)

      D <- Rprev - Rcurr

      D <- pmin(
        pmax(D, 1e-10),
        1 - 1e-10
      )

      ll <- ll + n[k, i] * log(D)

      nk <- nk + n[k, i]

      Rprev <- Rcurr
    }

    Rsurv <- pmin(
      pmax(Rprev, 1e-10),
      1 - 1e-10
    )

    ll <- ll +
      (N[k] - nk) * log(Rsurv)
  }

  return(-ll)
}


start <- c(-3, 4)

# Monte Carlo settings

R.reps <- 1000
t0 <- 3.5

Estimate <- matrix(
  NA,
  R.reps,
  2
)

colnames(Estimate) <- c(
  "beta0",
  "beta1"
)

Rhat <- rep(NA, R.reps)

valid <- rep(FALSE, R.reps)

Lower.theta <- matrix(
  NA,
  R.reps,
  2
)

Upper.theta <- matrix(
  NA,
  R.reps,
  2
)

Lower.R0 <- rep(NA, R.reps)
Upper.R0 <- rep(NA, R.reps)

CP.store <- matrix(
  NA,
  R.reps,
  5
)

AL.store <- matrix(
  NA,
  R.reps,
  5
)

z <- qnorm(0.975)


R0.true <- exp(
  -exp(beta0) * t0
)


q <- 0.5

MTTF.true <- exp(-beta0)

Q.true <- exp(-beta0) *
  (-log(1 - q))


MTTF.hat <- rep(
  NA,
  R.reps
)

Q.hat <- rep(
  NA,
  R.reps
)

Lower.MTTF <- rep(
  NA,
  R.reps
)

Upper.MTTF <- rep(
  NA,
  R.reps
)

Lower.Q <- rep(
  NA,
  R.reps
)

Upper.Q <- rep(
  NA,
  R.reps
)

SE.MTTF.store <- rep(
  NA,
  R.reps
)

SE.Q.store <- rep(
  NA,
  R.reps
)

discard <- 0


# Monte Carlo loop
SE.R0.store <- rep(
  NA,
  R.reps
)


for(r in 1:R.reps) {

  # Generate data
  
  failure <- matrix(
    0,
    K,
    I
  )

  for(k in 1:K) {

    Hcum <- 0
    Rprev <- 1

    D <- numeric(I)

    for(i in 1:I) {

      lambda <- exp(
        beta0 + beta1 * x[k, i]
      )

      if(i == 1) {

        Hcum <- lambda * tau[k, i]

      } else {

        Hcum <- Hcum +
          lambda *
          (tau[k, i] - tau[k, i - 1])
      }

      Rcurr <- exp(-Hcum)

      D[i] <- Rprev - Rcurr

      Rprev <- Rcurr
    }

    probs <- c(
      D,
      Rprev
    )

    draw <- rmultinom(
      1,
      Nk[k],
      probs
    )

    failure[k, ] <- draw[1:I]
  }



  # Maximum likelihood estimation

  fit <- try(

    optim(
      start,
      loglik,
      n = failure,
      N = Nk,
      x = x,
      tau = tau,
      method = "L-BFGS-B",

      lower = c(
        -10,
        -10
      ),

      upper = c(
        10,
        10
      )
    ),

    silent = TRUE
  )


  if(inherits(fit, "try-error"))
    next

  if(fit$convergence != 0)
    next


  # Analytical Fisher Information Matrix
  # For theta = (beta0, beta1)
  

  Info <- matrix(
    0,
    2,
    2
  )


  for(k in 1:K) {

    Hcum <- 0
    Rprev <- 1

    numB <- 0

    for(i in 1:I) {

      lambda <- exp(
        fit$par[1] +
        fit$par[2] * x[k, i]
      )

      if(i == 1) {

        deltaH <- tau[k, i]

      } else {

        deltaH <-
          tau[k, i] -
          tau[k, i - 1]
      }

      Hcum <- Hcum +
        lambda * deltaH

      numB <- numB +
        x[k, i] *
        lambda *
        deltaH

      B <- numB / Hcum

      Rcurr <- exp(-Hcum)

      D <- pmax(
        Rprev - Rcurr,
        1e-10
      )

      U_b0 <-
        Rcurr * Hcum -
        Rprev * ifelse(i == 1, 0, Hprev)

      U_b1 <-
        Rcurr * Hcum * B -
        Rprev *
        ifelse(i == 1, 0, Hprev * Bprev)

      Info[1, 1] <-
        Info[1, 1] +
        Nk[k] * U_b0^2 / D

      Info[1, 2] <-
        Info[1, 2] +
        Nk[k] * U_b0 * U_b1 / D

      Info[2, 2] <-
        Info[2, 2] +
        Nk[k] * U_b1^2 / D


      Hprev <- Hcum
      Rprev <- Rcurr
      Bprev <- B
    }

    Rsurv <- pmax(
      Rprev,
      1e-10
    )

    V_b0 <-
      -Rprev * Hprev

    V_b1 <-
      -Rprev * Hprev * Bprev


    Info[1, 1] <-
      Info[1, 1] +
      Nk[k] * V_b0^2 / Rsurv

    Info[1, 2] <-
      Info[1, 2] +
      Nk[k] * V_b0 * V_b1 / Rsurv

    Info[2, 2] <-
      Info[2, 2] +
      Nk[k] * V_b1^2 / Rsurv
  }


  Info[2, 1] <- Info[1, 2]

  # Variance-Covariance Matrix

  V <- try(
    solve(Info),
    silent = TRUE
  )


  if(inherits(V, "try-error")) {

    discard <- discard + 1

    next
  }


  if(any(!is.finite(V))) {

    discard <- discard + 1

    next
  }


  valid[r] <- TRUE

  Estimate[r, ] <- fit$par


  SE <- sqrt(
    pmax(diag(V), 0)
  )

  # CI for beta0 and beta1

  Lower.theta[r, ] <-
    fit$par - z * SE

  Upper.theta[r, ] <-
    fit$par + z * SE

  # Coverage Probability

  CP.store[r, 1] <-
    (beta0 >= Lower.theta[r, 1]) &
    (beta0 <= Upper.theta[r, 1])

  CP.store[r, 2] <-
    (beta1 >= Lower.theta[r, 2]) &
    (beta1 <= Upper.theta[r, 2])

  # Average Length

  AL.store[r, 1] <-
    Upper.theta[r, 1] -
    Lower.theta[r, 1]

  AL.store[r, 2] <-
    Upper.theta[r, 2] -
    Lower.theta[r, 2]


  Rhat[r] <-
    exp(
      -exp(fit$par[1]) * t0
    )

  MTTF.hat[r] <-
    exp(
      -fit$par[1]
    )

  Q.hat[r] <-
    exp(
      -fit$par[1]
    ) *
    (-log(1 - q))


  # Delta method for R0(t)

  grad <- matrix(
    c(
      -Rhat[r] *
        exp(fit$par[1]) *
        t0,

      0
    ),
    nrow = 1
  )


  Var.R0 <-
    grad %*%
    V %*%
    t(grad)

  SE.R0 <-
    sqrt(
      as.numeric(Var.R0)
    )

  SE.R0.store[r] <-
    SE.R0

  # Delta Method : MTTF

  grad.MTTF <- matrix(
    c(
      -MTTF.hat[r],
      0
    ),
    nrow = 1
  )


  Var.MTTF <-
    grad.MTTF %*%
    V %*%
    t(grad.MTTF)

  SE.MTTF <-
    sqrt(
      as.numeric(Var.MTTF)
    )

  SE.MTTF.store[r] <-
    SE.MTTF

  # Delta Method : Quantile
  grad.Q <- matrix(
    c(
      -Q.hat[r],
      0
    ),
    nrow = 1
  )

  Var.Q <-
    grad.Q %*%
    V %*%
    t(grad.Q)

  SE.Q <-
    sqrt(
      as.numeric(Var.Q)
    )

  SE.Q.store[r] <-
    SE.Q

  # Logit CI for R0(t)

  S <-
    exp(
      z * SE.R0 /
        (
          Rhat[r] *
          (1 - Rhat[r])
        )
    )


  Lower.R0[r] <-
    Rhat[r] /
    (
      Rhat[r] +
        (1 - Rhat[r]) * S
    )


  Upper.R0[r] <-
    Rhat[r] /
    (
      Rhat[r] +
        (1 - Rhat[r]) / S
    )

  # CP and AL for R0(t)
 

  CP.store[r, 3] <-
    (R0.true >= Lower.R0[r]) &
    (R0.true <= Upper.R0[r])


  AL.store[r, 3] <-
    Upper.R0[r] -
    Lower.R0[r]

  # Log CI : MTTF
  w <-
    z *
    SE.MTTF /
    MTTF.hat[r]


  Lower.MTTF[r] <-
    MTTF.hat[r] *
    exp(-w)

  Upper.MTTF[r] <-
    MTTF.hat[r] *
    exp(w)

  # Log CI : Quantile
  w <-
    z *
    SE.Q /
    Q.hat[r]


  Lower.Q[r] <-
    Q.hat[r] *
    exp(-w)

  Upper.Q[r] <-
    Q.hat[r] *
    exp(w)
  
  # CP and AL : MTTF
  CP.store[r, 4] <-
    MTTF.true >= Lower.MTTF[r] &
    MTTF.true <= Upper.MTTF[r]

  AL.store[r, 4] <-
    Upper.MTTF[r] -
    Lower.MTTF[r]

  # CP and AL : Quantile
  CP.store[r, 5] <-
    Q.true >= Lower.Q[r] &
    Q.true <= Upper.Q[r]

  AL.store[r, 5] <-
    Upper.Q[r] -
    Lower.Q[r]

}


valid_idx <- which(valid)

# Average asymptotic SD of R0(t)

ASD.R0 <-
  mean(
    SE.R0.store[valid_idx],
    na.rm = TRUE
  )
###########################

AE <-
  colMeans(
    Estimate[valid_idx, ],
    na.rm = TRUE
  )

AE.R0 <-
  mean(
    Rhat[valid_idx],
    na.rm = TRUE
  )

AE.MTTF <-
  mean(
    MTTF.hat[valid_idx],
    na.rm = TRUE
  )

AE.Q <-
  mean(
    Q.hat[valid_idx],
    na.rm = TRUE
  )

True <- c(
  beta0,
  beta1
)


MSE <-
  colMeans(
    (
      Estimate[valid_idx, ] -
        matrix(
          True,
          nrow = length(valid_idx),
          ncol = 2,
          byrow = TRUE
        )
    )^2,
    na.rm = TRUE
  )


MSE.R0 <-
  mean(
    (
      Rhat[valid_idx] -
        R0.true
    )^2,
    na.rm = TRUE
  )


MSE.MTTF <-
  mean(
    (
      MTTF.hat[valid_idx] -
        MTTF.true
    )^2,
    na.rm = TRUE
  )


MSE.Q <-
  mean(
    (
      Q.hat[valid_idx] -
        Q.true
    )^2,
    na.rm = TRUE
  )

RMSE <-
  sqrt(MSE)

RMSE.R0 <-
  sqrt(MSE.R0)

RMSE.MTTF <-
  sqrt(MSE.MTTF)

RMSE.Q <-
  sqrt(MSE.Q)

CP <-
  colMeans(
    CP.store[valid_idx, ],
    na.rm = TRUE
  )

AL <-
  colMeans(
    AL.store[valid_idx, ],
    na.rm = TRUE
  )

# Final Table

Result <- data.frame(

  Parameter = c(
    "beta0",
    "beta1",
    paste0("R0(", t0, ")"),
    "MTTF",
    "q(0.5)"
  ),

  True = c(
    beta0,
    beta1,
    R0.true,
    MTTF.true,
    Q.true
  ),

  AE = c(
    AE,
    AE.R0,
    AE.MTTF,
    AE.Q
  ),

  MSE = c(
    MSE,
    MSE.R0,
    MSE.MTTF,
    MSE.Q
  ),

  RMSE = c(
    RMSE,
    RMSE.R0,
    RMSE.MTTF,
    RMSE.Q
  ),

  CP = CP,

  AL = AL
)

Result[, -1] <-
  round(
    Result[, -1],
    5
  )

print(Result)

cat(
  "\nDiscarded (singular Fisher matrix):",
  discard,
  "\n"
)

cat(
  "Valid replications:",
  length(valid_idx),
  "\n"
)

cat(
  "\nAverage asymptotic SD of R0(t):",
  round(ASD.R0, 5),
  "\n"
)

# KERNEL DENSITY VS STANDARD NORMAL
# FOUR PANELs

valid_idx <- which(valid)

beta0_hat <-
  Estimate[
    valid_idx,
    "beta0"
  ]

beta1_hat <-
  Estimate[
    valid_idx,
    "beta1"
  ]


R0_hat <-
  Rhat[valid_idx]

MTTF_hat <-
  MTTF.hat[valid_idx]

Q_hat <-
  Q.hat[valid_idx]

SE_beta0 <- (

  Upper.theta[
    valid_idx,
    1
  ] -

    Lower.theta[
      valid_idx,
      1
    ]

) / (2 * z)


SE_beta1 <- (

  Upper.theta[
    valid_idx,
    2
  ] -

    Lower.theta[
      valid_idx,
      2
    ]

) / (2 * z)

SE_R0 <-
  SE.R0.store[
    valid_idx
  ]

SE_MTTF <-
  SE.MTTF.store[
    valid_idx
  ]

SE_Q <-
  SE.Q.store[
    valid_idx
  ]

ok_beta0 <-
  is.finite(beta0_hat) &
  is.finite(SE_beta0) &
  SE_beta0 > 0


ok_beta1 <-
  is.finite(beta1_hat) &
  is.finite(SE_beta1) &
  SE_beta1 > 0


ok_R0 <-
  is.finite(R0_hat) &
  is.finite(SE_R0) &
  SE_R0 > 0


ok_MTTF <-
  is.finite(MTTF_hat) &
  is.finite(SE_MTTF) &
  SE_MTTF > 0


ok_Q <-
  is.finite(Q_hat) &
  is.finite(SE_Q) &
  SE_Q > 0


#  STANDARDIZED ESTIMATORS
Z_beta0 <- (

  beta0_hat[ok_beta0] -
    beta0

) /

  SE_beta0[ok_beta0]


Z_beta1 <- (

  beta1_hat[ok_beta1] -
    beta1

) /

  SE_beta1[ok_beta1]


Z_R0 <- (

  R0_hat[ok_R0] -
    R0.true

) /

  SE_R0[ok_R0]


Z_MTTF <- (

  MTTF_hat[ok_MTTF] -
    MTTF.true

) /

  SE_MTTF[ok_MTTF]


Z_Q <- (

  Q_hat[ok_Q] -
    Q.true

) /

  SE_Q[ok_Q]


Z_beta0 <-
  Z_beta0[
    is.finite(Z_beta0)
  ]

Z_beta1 <-
  Z_beta1[
    is.finite(Z_beta1)
  ]

Z_R0 <-
  Z_R0[
    is.finite(Z_R0)
  ]

Z_MTTF <-
  Z_MTTF[
    is.finite(Z_MTTF)
  ]

Z_Q <-
  Z_Q[
    is.finite(Z_Q)
  ]

#KDE plotting function
plot_kde_normal <-
  function(
    z,
    main_title
  ) {

    kd <- density(
      z,
      kernel = "gaussian",
      bw = "SJ",
      n = 1024,
      from = -4,
      to = 4
    )
    
    plot(
      kd,
      type = "l",
      lwd = 2,
      col = "blue",
      xlim = c(-4, 4),
      ylim = c(0, 0.6),
      yaxt = "n",
      main = main_title,
      xlab = "Value",
      ylab = "Density",
      font.lab = 2,
      font.axis = 2
    )

    # Y-axis
    axis(
      2,
      at = c(
        0.2,
        0.4,
        0.6
      ),
      labels = c(
        "0.2",
        "0.4",
        "0.6"
      ),
      font = 2
    )
    # Standard normal curve
    xx <- seq(
      -4,
      4,
      length.out = 1024
    )
    lines(
      xx,
      dnorm(xx),
      col = "red",
      lwd = 2,
      lty = 2
    )

    legend(
      "topright",
      legend = c(
        "Kernel Density",
        "Standard Normal"
      ),
      col = c(
        "blue",
        "red"
      ),
      lty = c(
        1,
        2
      ),
      lwd = c(
        2,
        1.9
      ),
      bty = "o",
      cex = 1.3,
      text.font = 2
    )
  }

layout(
  matrix(
    c(1, 1,2,2,3,3,
      0, 4, 4,5,5,0),
    nrow = 2,
    byrow = TRUE
  ),
  widths = rep(1,6),
  heights = c(1, 1)
)


par(
  mar = c(4.5, 4.8, 2.2, 1),
  cex.axis = 1.8,
  cex.lab = 1.8,
  cex.main = 1.8
)

plot_kde_normal(
  Z_beta0,
  expression(hat(beta)[0])
)

plot_kde_normal(
  Z_beta1,
  expression(hat(beta)[1])
)

plot_kde_normal(
  Z_R0,
  expression(hat(R)[NOC](t[M]))
)

plot_kde_normal(
  Z_MTTF,
  "MTTF"
)

plot_kde_normal(
  Z_Q,
  "Median"
)

layout(1)
par(mfrow = c(1, 1))