# LED data: pooled (K = 2) vs each experiment alone (K = 1)
#  Chen baseline with interval monitoring

rm(list = ls())
set.seed(2026)

B  <- 1000        # bootstrap replications
t0 <- 5           # mission time (hours)
z  <- qnorm(0.975)

# ---- data ----
Nk.all  <- c(30, 30)
x.all   <- rbind(c(0,     1, 1),      # Exp I : 7 V -> 8 V
                 c(0.517, 1, 1))      # Exp II: 7.5 V -> 8 V
tau.all <- rbind(c(100/60, 600/60, 1340/60),
                 c(100/60, 600/60, 1137/60))
n.all   <- rbind(c(11, 10, 7),
                 c(16,  8, 4))

# Chen baseline cumulative hazard
H0 <- function(t, alpha) exp(t^alpha) - 1

# negative log-likelihood, par = (beta0, log(alpha), beta1) 
negloglik <- function(par, n, N, x, tau) {
  beta0 <- par[1]
  alpha <- exp(par[2])
  beta1 <- par[3]
  ll <- 0
  for (k in 1:nrow(x)) {
    Hcum <- 0; Rprev <- 1; nk <- 0
    for (i in 1:ncol(x)) {
      lambda <- exp(beta0 + beta1 * x[k, i])
      if (i == 1) {
        Hcum <- lambda * (exp(tau[k, i]^alpha) - 1)
      } else {
        Hcum <- Hcum + lambda * (exp(tau[k, i]^alpha) - exp(tau[k, i - 1]^alpha))
      }
      Rcurr <- exp(-Hcum)
      if (!is.finite(Rcurr)) return(1e10)
      D <- pmin(pmax(Rprev - Rcurr, 1e-10), 1 - 1e-10)
      ll <- ll + n[k, i] * log(D)
      nk <- nk + n[k, i]
      Rprev <- Rcurr
    }
    Rsurv <- pmin(pmax(Rprev, 1e-10), 1 - 1e-10)
    ll <- ll + (N[k] - nk) * log(Rsurv)
  }
  if (!is.finite(ll)) return(1e10)
  -ll
}

# MLE
fit_mle <- function(n, N, x, tau, start) {
  fit <- try(optim(c(start[1], log(start[2]), start[3]), negloglik,
                   n = n, N = N, x = x, tau = tau, method = "BFGS",
                   hessian = TRUE, control = list(maxit = 2000)),
             silent = TRUE)
  if (inherits(fit, "try-error") || fit$convergence != 0) return(NULL)
  
  alpha <- exp(fit$par[2])
  theta <- c(fit$par[1], alpha, fit$par[3])
  
  # Hessian from (beta0, log(alpha), beta1) to (beta0, alpha, beta1)
  J <- diag(c(1, 1 / alpha, 1))
  H_theta <- J %*% fit$hessian %*% J
  
  gr <- get_grad(function(pp) negloglik(pp, n, N, x, tau), fit$par)
  list(theta = theta, hessian = H_theta,
       grad_norm = sqrt(sum(gr^2)),
       kappa = tryCatch(kappa(H_theta), error = function(e) NA))
}

#  numerical gradient for trapezoidal rule
get_grad <- function(func, x, h = 1e-5) {
  n <- length(x)
  g <- numeric(n)
  for (i in 1:n) {
    x1 <- x; x2 <- x
    x1[i] <- x[i] - h; x2[i] <- x[i] + h
    g[i] <- (func(x2) - func(x1)) / (2 * h)
  }
  g
}

# trapezoidal rule on [0, Tmax]
trap <- function(f, Tmax, n = 20000) {
  u  <- seq(0, 1, length.out = n + 1)
  tt <- Tmax * u^3
  y  <- f(tt)
  sum(diff(tt) * (head(y, -1) + tail(y, -1)) / 2)
}

# estimates and CIs 
qnames <- c("beta0", "alpha", "beta1", "R_NOC(t_M)",
            "MTTF", "t_0.10", "t_0.50")

all_quantities <- function(theta, V) {
  b0 <- theta[1]; a <- theta[2]; b1 <- theta[3]
  se <- sqrt(diag(V))
  
  ci.b0 <- b0 + c(-1, 1) * z * se[1]    
  ci.a  <- a + c(-1, 1) * z * se[2]   
  ci.b1 <- b1 + c(-1, 1) * z * se[3]    
  
  # reliability at t0, logit-scale CI
  calc_R <- function(p) exp(-exp(p[1]) * (exp(t0^p[2]) - 1))
  R   <- calc_R(theta)
  gR  <- get_grad(calc_R, theta)
  seR <- sqrt(as.numeric(t(gR) %*% V %*% gR))
  S   <- exp(z * seR / (R * (1 - R)))
  ci.R <- c(R / (R + (1 - R) * S), R / (R + (1 - R) / S))
  
  # MTTF = integral of R(t) over (0, Inf)
  Rf   <- function(t) exp(-exp(b0) * (exp(t^a) - 1))
  Tmax <- (log(1 + (-log(1e-12)) * exp(-b0)))^(1 / a)
  if (!is.finite(Tmax) || Tmax <= 0) {
    M <- NA; ci.M <- c(NA, NA)
  } else {
    M   <- trap(Rf, Tmax)
    # gradient by differentiating under the integral
    dMb <- -trap(function(t) Rf(t) * exp(b0) * (exp(t^a) - 1), Tmax)
    dMa <- -trap(function(t) ifelse(t > 0,
                                    Rf(t) * exp(b0) * exp(t^a) * t^a * log(t), 0), Tmax)
    gM   <- c(dMb, dMa, 0)
    seM  <- sqrt(as.numeric(t(gM) %*% V %*% gM))
    ci.M <- M * exp(c(-1, 1) * z * seM / M)    # log scale
  }
  
  # quantiles, log-scale CI
  qfun <- function(q) {
    calc_Q <- function(p) (log(1 + (-log(1 - q)) * exp(-p[1])))^(1 / p[2])
    Q   <- calc_Q(theta)
    gQ  <- get_grad(calc_Q, theta)
    seQ <- sqrt(as.numeric(t(gQ) %*% V %*% gQ))
    c(Q, Q * exp(c(-1, 1) * z * seQ / Q))
  }
  Q1 <- qfun(0.10); Q5 <- qfun(0.50)
  
  list(est = c(b0, a, b1, R, M, Q1[1], Q5[1]),
       low = c(ci.b0[1], ci.a[1], ci.b1[1], ci.R[1], ci.M[1], Q1[2], Q5[2]),
       upp = c(ci.b0[2], ci.a[2], ci.b1[2], ci.R[2], ci.M[2], Q1[3], Q5[3]))
}

# pooled fit on the real data (reference values)
fit_ref   <- fit_mle(n.all, Nk.all, x.all, tau.all, start = c(-1, 0.5, 0.3))
theta_ref <- fit_ref$theta
V_ref     <- solve(fit_ref$hessian)
real      <- all_quantities(theta_ref, V_ref)
ref_value <- real$est

real_tab <- data.frame(Quantity = qnames, Estimate = real$est,
                       Lower95 = real$low, Upper95 = real$upp)
write.csv(real_tab, "LED_real_data_estimates.csv", row.names = FALSE)

# LR test of beta1 = 0
nll_full <- negloglik(c(theta_ref[1], log(theta_ref[2]), theta_ref[3]),
                      n.all, Nk.all, x.all, tau.all)
nll_null <- optim(c(theta_ref[1], log(theta_ref[2])),
                  function(p) negloglik(c(p, 0), n.all, Nk.all, x.all, tau.all),
                  method = "BFGS")$value
LR <- 2 * (nll_null - nll_full)
cat("Pooled LED fit:\n"); print(real_tab, digits = 4)
cat("\nLR test of beta1 = 0: LR =", round(LR, 3),
    ", p-value =", round(pchisq(LR, 1, lower.tail = FALSE), 3), "\n")

# AIC / BIC
logLik_full  <- -nll_full
n_parameters <- 3
n_units      <- sum(Nk.all)

AIC_pooled <- -2 * logLik_full + 2 * n_parameters
BIC_pooled <- -2 * logLik_full + log(n_units) * n_parameters

cat("\nPooled model information criteria\n")
cat("Log-likelihood =", round(logLik_full, 4), "\n")
cat("AIC =", round(AIC_pooled, 4), "\n")
cat("BIC =", round(BIC_pooled, 4), "\n")

#  Pearson chi-square statistic 
calc_chi2 <- function(theta, n, N, x, tau) {
  b0 <- theta[1]; a <- theta[2]; b1 <- theta[3]
  chi2 <- 0
  for (k in seq_len(nrow(x))) {
    Hcum <- 0; Rprev <- 1; nk <- 0
    for (i in seq_len(ncol(x))) {
      lambda <- exp(b0 + b1 * x[k, i])
      deltaH <- if (i == 1) (exp(tau[k, i]^a) - 1) else
        (exp(tau[k, i]^a) - exp(tau[k, i - 1]^a))
      Hcum  <- Hcum + lambda * deltaH
      Rcurr <- exp(-Hcum)
      
      D        <- pmax(Rprev - Rcurr, 1e-10)
      expected <- N[k] * D
      observed <- n[k, i]
      chi2 <- chi2 + (observed - expected)^2 / expected
      
      nk    <- nk + observed
      Rprev <- Rcurr
    }
    Rsurv         <- pmax(Rprev, 1e-10)
    expected_surv <- N[k] * Rsurv
    observed_surv <- N[k] - nk
    chi2 <- chi2 + (observed_surv - expected_surv)^2 / expected_surv
  }
  chi2
}

obs_chi2 <- calc_chi2(theta_ref, n.all, Nk.all, x.all, tau.all)

#  parametric bootstrap GOF test 
boot_chi2 <- rep(NA_real_, B)

for (b in seq_len(B)) {
  boot_fail <- matrix(0L, nrow(x.all), ncol(x.all))
  
  for (k in seq_len(nrow(x.all))) {
    Hcum <- 0; Rprev <- 1; D <- numeric(ncol(x.all))
    for (i in seq_len(ncol(x.all))) {
      lambda <- exp(theta_ref[1] + theta_ref[3] * x.all[k, i])
      deltaH <- if (i == 1) H0(tau.all[k, i], theta_ref[2]) else
        (H0(tau.all[k, i], theta_ref[2]) - H0(tau.all[k, i - 1], theta_ref[2]))
      Hcum  <- Hcum + lambda * deltaH
      Rcurr <- exp(-Hcum)
      D[i]  <- pmax(Rprev - Rcurr, 1e-10)
      Rprev <- Rcurr
    }
    probs <- c(D, max(Rprev, 1e-10))
    draw  <- rmultinom(1, Nk.all[k], probs)
    boot_fail[k, ] <- draw[seq_len(ncol(x.all))]
  }
  
  boot_fit <- fit_mle(n = boot_fail, N = Nk.all, x = x.all,
                      tau = tau.all, start = theta_ref)
  
  boot_chi2[b] <- if (!is.null(boot_fit))
    calc_chi2(boot_fit$theta, boot_fail, Nk.all, x.all, tau.all) else NA
}

valid_boot  <- boot_chi2[is.finite(boot_chi2)]
p_value_gof <- mean(boot_chi2 > obs_chi2, na.rm = TRUE)

cat("\nParametric bootstrap chi-square GOF\n")
cat("Observed chi-square =", round(obs_chi2, 4), "\n")
cat("Valid bootstrap replicates =", length(valid_boot), "of", B, "\n")
cat("Bootstrap GOF p-value =", round(p_value_gof, 4), "\n")

gof_ic_tab <- data.frame(
  Model = "Pooled Chen FRM",
  LogLik = logLik_full,
  Parameters = n_parameters,
  AIC = AIC_pooled,
  BIC = BIC_pooled,
  ChiSquare = obs_chi2,
  Bootstrap_GOF_P_value = p_value_gof,
  Bootstrap_replicates_valid = length(valid_boot)
)

write.csv(gof_ic_tab, "LED_GOF_AIC_BIC_Chen.csv", row.names = FALSE)
cat("\nSaved: LED_GOF_AIC_BIC_Chen.csv\n")

# parametric bootstrap for one design 
run_design <- function(idx) {
  Nk  <- Nk.all[idx]
  x   <- x.all[idx, , drop = FALSE]
  tau <- tau.all[idx, , drop = FALSE]
  K <- nrow(x); I <- ncol(x)
  b0 <- theta_ref[1]; a <- theta_ref[2]; b1 <- theta_ref[3]
  
  EST <- LOW <- UPP <- matrix(NA, B, 7, dimnames = list(NULL, qnames))
  status <- character(B)
  DIAG <- matrix(NA, B, 2, dimnames = list(NULL, c("grad_norm", "kappa")))
  
  for (r in 1:B) {
    # simulate interval counts
    failure <- matrix(0, K, I)
    for (k in 1:K) {
      Hcum <- 0; Rprev <- 1; D <- numeric(I)
      for (i in 1:I) {
        lambda <- exp(b0 + b1 * x[k, i])
        if (i == 1) Hcum <- lambda * H0(tau[k, i], a)
        else Hcum <- Hcum + lambda * (H0(tau[k, i], a) - H0(tau[k, i - 1], a))
        Rcurr <- exp(-Hcum); D[i] <- Rprev - Rcurr; Rprev <- Rcurr
      }
      failure[k, ] <- rmultinom(1, Nk[k], c(D, Rprev))[1:I]
    }
    
    fit_b <- fit_mle(failure, Nk, x, tau, start = theta_ref)
    if (is.null(fit_b)) { status[r] <- "no convergence"; next }
    th <- fit_b$theta
    
    V <- try(solve(fit_b$hessian), silent = TRUE)
    if (inherits(V, "try-error") || any(!is.finite(V)) || any(diag(V) <= 0) ||
        any(eigen(fit_b$hessian, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
      status[r] <- "singular information"; next
    }
    
    qs <- suppressWarnings(all_quantities(th, V))
    if (any(!is.finite(qs$est))) { status[r] <- "non-finite estimate"; next }
    
    lo <- qs$low; up <- qs$upp
    lo[!is.finite(lo)] <- NA; up[!is.finite(up)] <- NA
    status[r] <- "ok"
    EST[r, ] <- qs$est; LOW[r, ] <- lo; UPP[r, ] <- up
    DIAG[r, ] <- c(fit_b$grad_norm, fit_b$kappa)
  }
  
  ok <- status == "ok"
  list(est = EST[ok, , drop = FALSE],
       diag = DIAG[ok, , drop = FALSE],
       low = LOW[ok, , drop = FALSE],
       upp = UPP[ok, , drop = FALSE],
       status = table(factor(status, levels = c("ok", "no convergence",
                                                "singular information",
                                                "non-finite estimate"))))
}

summarise_design <- function(res, name) {
  E <- res$est; L <- res$low; U <- res$upp
  covered <- sweep(L, 2, ref_value, "<=") & sweep(U, 2, ref_value, ">=")
  data.frame(
    Design    = name,
    Quantity  = qnames,
    Reference = ref_value,
    AE        = colMeans(E),
    Bias      = colMeans(E) - ref_value,
    RMSE      = sqrt(colMeans(sweep(E, 2, ref_value)^2)),
    CP        = colMeans(covered, na.rm = TRUE),
    AL        = colMeans(U - L, na.rm = TRUE),
    row.names = NULL
  )
}

designs <- list("Exp I" = 1, "Exp II" = 2, "Pooled" = 1:2)
mc <- lapply(designs, run_design)

summary_tab <- do.call(rbind, Map(summarise_design, mc, names(mc)))

#  % reduction in RMSE and AL, pooled vs single experiment 
# reduction = 100 * (1 - pooled / single).
get_val <- function(design, quantity, what)
  summary_tab[summary_tab$Design == design & summary_tab$Quantity == quantity, what]

reduction <- function(pooled, single) {
  if (!is.finite(single) || !is.finite(pooled) || single <= 0) return(NA)
  if (single > 1e6 * pooled) return(NA)
  round(100 * (1 - pooled / single), 1)
}

reduction_tab <- do.call(rbind, lapply(qnames, function(qn) {
  rm1 <- get_val("Exp I", qn, "RMSE");  rm2 <- get_val("Exp II", qn, "RMSE")
  rmp <- get_val("Pooled", qn, "RMSE")
  al1 <- get_val("Exp I", qn, "AL");    al2 <- get_val("Exp II", qn, "AL")
  alp <- get_val("Pooled", qn, "AL")
  data.frame(
    Quantity                = qn,
    RMSE_ExpI               = signif(rm1, 4),
    RMSE_ExpII              = signif(rm2, 4),
    RMSE_Pooled             = signif(rmp, 4),
    RMSE_reduction_vs_ExpI  = reduction(rmp, rm1),
    RMSE_reduction_vs_ExpII = reduction(rmp, rm2),
    AL_ExpI                 = signif(al1, 4),
    AL_ExpII                = signif(al2, 4),
    AL_Pooled               = signif(alp, 4),
    AL_reduction_vs_ExpI    = reduction(alp, al1),
    AL_reduction_vs_ExpII   = reduction(alp, al2),
    row.names = NULL)
}))

cat("\nPercentage reduction (pooled vs single experiment)\n")
cat("NA = single-experiment value not finite or dominated by extreme replications\n")
print(reduction_tab)
write.csv(reduction_tab, "LED_reduction_percent.csv", row.names = FALSE)

#  diagnostics 
diagnose <- function(res, name) {
  E <- res$est; st <- res$status
  data.frame(
    Design                 = name,
    Excluded_pct           = round(100 * (1 - st["ok"] / sum(st)), 1),
    NoConvergence_pct      = round(100 * st["no convergence"] / sum(st), 1),
    Singular_pct           = round(100 * st["singular information"] / sum(st), 1),
    NonFinite_pct          = round(100 * st["non-finite estimate"] / sum(st), 1),
    Alpha_below_0.05_pct   = round(100 * mean(E[, "alpha"] < 0.05), 1),
    Abs_beta0_above_10_pct = round(100 * mean(abs(E[, "beta0"]) > 10), 1),
    Abs_beta1_above_10_pct = round(100 * mean(abs(E[, "beta1"]) > 10), 1),
    Median_kappa           = signif(median(res$diag[, "kappa"], na.rm = TRUE), 3),
    Max_grad_norm          = signif(max(res$diag[, "grad_norm"], na.rm = TRUE), 3),
    row.names = NULL
  )
}
diag_tab <- do.call(rbind, Map(diagnose, mc, names(mc)))
ci_missing <- do.call(rbind, lapply(names(mc), function(d)
  data.frame(Design = d, t(round(100 * colMeans(is.na(mc[[d]]$low)), 1)),
             check.names = FALSE)))
cat("\nDiagnostics per design\n");       print(diag_tab)
cat("\nMissing CIs (%) per quantity\n"); print(ci_missing)
write.csv(diag_tab,   "LED_diagnostics.csv",        row.names = FALSE)
write.csv(ci_missing, "LED_CI_missing_percent.csv", row.names = FALSE)
excluded_tab <- do.call(rbind, lapply(names(mc), function(d)
  data.frame(Design = d, t(as.matrix(mc[[d]]$status)), check.names = FALSE)))

write.csv(summary_tab,  "LED_simulation_summary.csv",    row.names = FALSE)
write.csv(excluded_tab, "LED_excluded_replications.csv", row.names = FALSE)

cat("\nReplications by status:\n"); print(excluded_tab)
cat("\nSummary (reference = pooled estimate from LED data):\n")
print(summary_tab, digits = 4)

# plot settings 
col_line <- c("Exp I" = "grey55", "Exp II" = "#1F3A68", "Pooled" = "black")
lty_line <- c("Exp I" = 1,        "Exp II" = 2,         "Pooled" = 1)
lwd_line <- c("Exp I" = 2,        "Exp II" = 2,         "Pooled" = 2.5)
col_fill <- c("Exp I" = "yellow", "Exp II" = "#20B2E8", "Pooled" = "#C8A2C8")
col_med  <- c("Exp I" = "black",  "Exp II" = "black",   "Pooled" = "black")
col_ref  <- "red"

panels <- list(
  list(j = 1, main = expression(beta[0]),
       lab = "Estimate", lim = c(-3, 1), pos = FALSE),
  list(j = 2, main = expression (alpha),
       lab = "Estimate", lim = c(0, 2.5), pos = TRUE),
  list(j = 3, main = expression(beta[1]),
       lab = "Estimate", lim = c(-5, 4), pos = FALSE),
  list(j = 4, main = expression(R[NOC](5)),
       lab = "Estimate", lim = c(0, 1), pos = TRUE),
  list(j = 5, main = "MTTF",
       lab = "Estimate", lim = c(0, 60), pos = TRUE),
  list(j = 6, main = expression(bold(t[0.10])),
       lab = "Estimate", lim = c(0, 0.5), pos = TRUE),
  list(j = 7, main = "Median lifetime",
       lab = "Estimate", lim = c(0, 15), pos = TRUE)
)

share_outside <- function(v, lim) mean(v < lim[1] | v > lim[2])


# boxplots 
outside_tab <- NULL
pdf("LED_boxplots.pdf", width = 12, height = 7)
par(mfrow = c(2, 4), mar = c(5, 5, 3, 3.2), font.axis = 2, font.lab = 2,
    font.main = 2, cex.axis = 1.4, cex.lab = 1.8, cex.main = 1.8)
for (p in panels) {
  vals <- lapply(mc, function(m) m$est[, p$j])
  
 
  st   <- lapply(vals, function(v) boxplot.stats(v)$stats)
  ylo  <- min(sapply(st, `[`, 1)); yhi <- max(sapply(st, `[`, 5))
  pad  <- 0.05 * (yhi - ylo)
  ylim <- c(ylo - pad, yhi + pad)
  if (p$pos) ylim[1] <- max(0, ylim[1])
  
  boxplot(vals, ylim = ylim, main = p$main, ylab = p$lab,
          names = names(mc), col = col_fill, border = "black",
          medcol = col_med, outpch = 16, outcex = 0.5, outcol = "grey40",
          whisklty = 1)
  abline(h = ref_value[p$j], col = col_ref, lty = 3, lwd = 3)
  out <- sapply(vals, share_outside, lim = ylim)
  outside_tab <- rbind(outside_tab,
                       data.frame(Quantity = qnames[p$j], Lower = ylim[1], Upper = ylim[2],
                                  t(round(100 * out, 1)), check.names = FALSE))
}
plot.new()

legend("center",
       legend = c("Exp I", "Exp II", "Pooled Exp", "Pooled estimate"),
       fill = c(col_fill, NA),
       border = c("black", "black", "black", NA),
       lty = c(NA, NA, NA, 3),
       lwd = c(1.5, 1.5, 1.5, 3),
       col = c(NA, NA, NA, col_ref),
       box.col = "black",
       box.lwd = 1.5,
       cex = 1.5,
       text.font = 2,
       pt.cex = 1.6,
       x.intersp = 1.2,
       y.intersp = 1.2,
       xpd = NA)
dev.off()

# share of estimates outside the plotted range
write.csv(outside_tab,  "LED_boxplot_outside_range.csv", row.names = FALSE)
cat("\nPercentage of estimates outside the plotted boxplot range:\n")
print(outside_tab)
