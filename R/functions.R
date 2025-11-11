# ------------------------------------------------------------------------------
# Core least-squares loss (sum of squared residuals)
# ------------------------------------------------------------------------------


ls_loss <- function(beta, X, y) {
  # X: n x p (already includes intercept if desired)
  # y: n-vector
  # beta: p-vector
  res <- y - as.vector(X %*% beta)
  sum(res^2)
}

# ------------------------------------------------------------------------------
# Fitter: my_lm()
#   - Minimizes ls_loss with optim()
#   - Pairs bootstrap of rows to estimate vcov(N2L)
# ------------------------------------------------------------------------------

#' @title Linear regression via optimization
#'
#' @description Uses a least-square loss to estimate coefficients and provide inference tools
#' @param X A \code{matrix} (not a data frame) of dimension nxp representing numerical predictors
#' @param y A \code{vector} that represents the numeric response.
#' @param B A \code{numeric} (integer) used to denote the number of bootstrap replications.
#' @param seed A \code{numeric} used to control the seed of the random number
#' generator used by this function.
#' @return A \code{list} of class ls_optim:
#' \describe{
#'      \item{I}{Estimated value of the integral}
#'      \item{var}{Estimated variance of the estimator}
#' }
#' @author Roberto Molinari
#' @importFrom stats runif
#' @export
#' @examples
#' set.seed(1)
#' n <- 120; p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("x", 1:p)
#' beta_true <- c(2, -1.5, 0, 0, 1, rep(0, p - 5))
#' y <- as.vector(X %*% beta_true + rnorm(n, sd = 1))
#' fit <- my_lm(X, y, B = 200, add_intercept = TRUE)
my_lm <- function(X, y, B = 100, add_intercept = TRUE, optim_method = "BFGS", optim_control = list(maxit = 1000), ...) {

  stopifnot(is.matrix(X), is.numeric(y), nrow(X) == length(y))

  if (add_intercept) {
    X <- cbind(`(Intercept)` = 1, X)
  }

  n <- nrow(X)
  p <- ncol(X)

  beta0 <- rep(0, p)

  fit_opt <- optim(
    par     = beta0,
    fn      = ls_loss,
    X       = X,
    y       = y,
    method  = optim_method,
    control = optim_control,
    ...
  )

  beta_hat <- as.numeric(fit_opt$par)
  conv <- fit_opt$convergence
  if (conv != 0) warning("optim() did not report successful convergence (code ", conv, ").")

  fitted_vals <- as.vector(X %*% beta_hat)
  resid_vals  <- as.vector(y - fitted_vals)
  rss <- sum(resid_vals^2)
  df_resid <- n - p
  sigma2_hat <- rss / df_resid

  betas_boot <- matrix(NA_real_, nrow = B, ncol = p)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    Xb  <- X[idx, , drop = FALSE]
    yb  <- y[idx]

    fit_b <- try(
      optim(
        par     = beta0,
        fn      = ls_loss,
        X       = Xb,
        y       = yb,
        method  = optim_method,
        control = optim_control
      ),
      silent = TRUE
    )
    if (inherits(fit_b, "try-error") || fit_b$convergence != 0) {
      next
    }
    betas_boot[b, ] <- fit_b$par
  }

  keep <- stats::complete.cases(betas_boot)
  betas_boot <- betas_boot[keep, , drop = FALSE]

  if (nrow(betas_boot) < 2)
    warning("Too few successful bootstrap replicates to estimate covariance.")

  vcov_boot <- if (nrow(betas_boot) >= 2) stats::cov(betas_boot) else matrix(NA_real_, p, p)

  out <- list(
    coefficients = setNames(beta_hat, colnames(X)),
    vcov         = vcov_boot,
    betas_boot   = betas_boot,
    fitted       = fitted_vals,
    residuals    = resid_vals,
    sigma2_hat   = sigma2_hat,
    df_residual  = df_resid,
    call         = match.call(),
    X            = X,
    y            = y,
    converged    = (conv == 0),
    B            = B,
    add_intercept = add_intercept
  )
  class(out) <- "ls_optim"
  out
}

# ------------------------------------------------------------------------------
# Methods: print(), coef(), vcov(), fitted(), residuals(), predict()
# ------------------------------------------------------------------------------
print.ls_optim <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nConverged:", x$converged, " |  Bootstrap reps (kept):", nrow(x$betas_boot), "\n")
  invisible(x)
}

coef.ls_optim <- function(object, ...) object$coefficients
vcov.ls_optim <- function(object, ...) object$vcov
fitted.ls_optim <- function(object, ...) object$fitted
residuals.ls_optim <- function(object, ...) object$residuals

predict.ls_optim <- function(object, newdata = NULL, ...) {
  X <- object$X
  beta <- object$coefficients

  if (is.null(newdata)) {
    return(as.vector(X %*% beta))
  }

  Xnew <- as.matrix(newdata)
  if (object$add_intercept) {
    Xnew <- cbind(`(Intercept)` = 1, Xnew)
  }
  as.vector(Xnew %*% beta)
}

# ------------------------------------------------------------------------------
# Summary method: coef table with bootstrap SEs, z/t stats, p-values, CIs
# ------------------------------------------------------------------------------
summary.ls_optim <- function(object, conf.level = 0.95, use.t = TRUE, ...) {
  beta <- object$coefficients
  V    <- object$vcov
  se   <- sqrt(diag(V))

  if (any(!is.finite(se))) {
    warning("Bootstrap vcov not available; using homoskedastic OLS-style SEs from sigma^2 (X'X)^{-1}.")
    XtX_inv <- try(solve(crossprod(object$X)), silent = TRUE)
    if (inherits(XtX_inv, "try-error")) {
      se <- rep(NA_real_, length(beta))
    } else {
      se <- sqrt(diag(XtX_inv) * object$sigma2_hat)
    }
  }

  df <- object$df_residual
  alpha <- 1 - conf.level
  crit <- if (use.t) stats::qt(1 - alpha/2, df = max(df, 1)) else stats::qnorm(1 - alpha/2)

  stat <- beta / se
  pval <- if (use.t) 2 * stats::pt(-abs(stat), df = max(df, 1)) else 2 * stats::pnorm(-abs(stat))

  ci_l <- beta - crit * se
  ci_u <- beta + crit * se

  tab <- cbind(Estimate = beta, `Std. Error` = se, `Statistic` = stat, `Pr(>|t|)` = pval,
               `CI Lower` = ci_l, `CI Upper` = ci_u)
  class(tab) <- c("coef_table_ls_optim", class(tab))

  ans <- list(
    call        = object$call,
    coefficients= tab,
    df_residual = df,
    conf.level  = conf.level,
    use.t       = use.t,
    B           = object$B
  )
  class(ans) <- "summary.ls_optim"
  ans
}

print.summary.ls_optim <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n"); print(x$call)
  cat("\nBootstrap-based coefficient table:\n")
  printCoefmat(x$coefficients[, c("Estimate","Std. Error","Statistic","Pr(>|t|)")],
               digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  cli <- paste0(round(100 * x$conf.level), "%")
  cat(sprintf("\nConfidence intervals (%s):\n", cli))
  print(x$coefficients[, c("CI Lower","CI Upper")], digits = digits)
  cat("\nResidual df:", x$df_residual, " | Bootstrap reps used:", attr(x$coefficients, "nboot") %||% NA, "\n")
  invisible(x)
}

# Helper to avoid importing rlang
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------------
# Plot marginals:
#  - Scatter of y vs x_j
#  - Marginal fitted line: alpha + beta_j * x_j + sum_{k!=j} beta_k * mean(X_k)
# ------------------------------------------------------------------------------
plot_marginals <- function(object, variables = NULL, point_alpha = 0.6, point_cex = 0.8,
                           line_lwd = 2, npoints_line = 200) {
  stopifnot(inherits(object, "ls_optim"))
  X <- object$X
  y <- object$y
  beta <- object$coefficients

  cn <- colnames(X)
  has_int <- any(cn == "(Intercept)")
  idx_vars <- if (has_int) seq_len(ncol(X))[-1] else seq_len(ncol(X))
  var_names <- cn[idx_vars]

  if (!is.null(variables)) {
    if (is.numeric(variables)) {
      idx_vars <- intersect(idx_vars, variables)
      var_names <- cn[idx_vars]
    } else {
      idx_vars <- which(cn %in% variables)
      var_names <- cn[idx_vars]
    }
  }

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar), add = TRUE)
  k <- length(idx_vars)
  nrow_plot <- ceiling(k^0.5)
  ncol_plot <- ceiling(k / nrow_plot)
  par(mfrow = c(nrow_plot, ncol_plot), mar = c(4,4,2,1))

  Xbar <- colMeans(X)
  for (j in idx_vars) {
    xj <- X[, j]
    plot(xj, y,
         xlab = cn[j],
         ylab = "y",
         main = paste("Marginal fit for", cn[j]),
         cex = point_cex, col = rgb(0,0,0,point_alpha), pch = 16)

    xg <- seq(min(xj), max(xj), length.out = npoints_line)

    baseline <- 0
    if (has_int) baseline <- baseline + beta[cn == "(Intercept)"] * 1
    others <- setdiff(seq_along(beta), c(which(cn == "(Intercept)"), j))
    if (length(others)) {
      baseline <- baseline + sum(beta[others] * Xbar[others])
    }
    yg <- baseline + beta[j] * xg

    lines(xg, yg, lwd = line_lwd)
  }

  invisible(NULL)
}



################################################################################
# Example usage
################################################################################
set.seed(1)
n <- 120; p <- 8
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
beta_true <- c(2, -1.5, 0, 0, 1, rep(0, p - 5))
y <- as.vector(X %*% beta_true + rnorm(n, sd = 1))

fit <- my_lm(X, y, B = 200, add_intercept = TRUE)
print(fit)
s <- summary(fit)
print(s)
plot_marginals(fit)   # base R panel of marginal fits

# Predictions:
preds <- predict(fit, newdata = X[1:5, ])
preds
