#' Least-squares loss for linear regression
#'
#' Computes the sum of squared residuals \eqn{\sum_i (y_i - x_i^\top \beta)^2}
#' for a given coefficient vector \code{beta}, design matrix \code{X}, and
#' response vector \code{y}. Intended for use as an objective in \code{\link[stats]{optim}}.
#'
#' @param beta Numeric vector of length \eqn{p}; regression coefficients.
#' @param X Numeric matrix \eqn{n \times p}. Should already include an intercept
#'   column if desired.
#' @param y Numeric vector of length \eqn{n}; responses.
#'
#' @return A single numeric value: the sum of squared residuals.
ls_loss <- function(beta, X, y) {
  # X: n x p (already includes intercept if desired)
  # y: n-vector
  # beta: p-vector
  res <- y - as.vector(X %*% beta)
  sum(res^2)
}

#' Least-squares via \code{optim()} with bootstrap covariance
#'
#' Fits a linear model by minimizing the least-squares loss using
#' \code{\link[stats]{optim}}, then estimates \eqn{\mathrm{Cov}(\hat\beta)} via a
#' pairs (row) bootstrap of size \code{B}.
#'
#' @param X Numeric matrix \eqn{n \times p} of predictors. If \code{add_intercept = TRUE},
#'   an intercept column is prepended internally.
#' @param y Numeric response vector of length \eqn{n}.
#' @param B Integer; number of bootstrap resamples (default \code{100}).
#' @param add_intercept Logical; if \code{TRUE} (default) add an intercept column to \code{X}.
#' @param optim_method Character; passed to \code{\link[stats]{optim}}'s \code{method} (default \code{"BFGS"}).
#' @param optim_control List; control list passed to \code{\link[stats]{optim}} (default \code{list(maxit = 1000)}).
#' @param ... Additional arguments forwarded to \code{\link[stats]{optim}}.
#'
#' @details
#' The function performs:
#' \enumerate{
#' \item Deterministic fit of \eqn{\hat\beta} by minimizing \code{\link{ls_loss}} with \code{\link[stats]{optim}}.
#' \item Pairs bootstrap: resample rows of \code{(X, y)} \code{B} times, refit, and compute
#'   the sample covariance matrix of the bootstrap coefficients.
#' }
#' If too few bootstrap fits succeed, the covariance matrix is filled with \code{NA}.
#'
#' @return An object of class \code{"ls_optim"} with elements:
#' \itemize{
#' \item \code{coefficients} — numeric named vector of \eqn{\hat\beta}.
#' \item \code{vcov} — bootstrap covariance matrix of \eqn{\hat\beta}.
#' \item \code{betas_boot} — matrix of bootstrap coefficient draws (rows).
#' \item \code{fitted}, \code{residuals}, \code{sigma2_hat}, \code{df_residual}.
#' \item \code{X}, \code{y}, \code{call}, \code{converged}, \code{B}, \code{add_intercept}.
#' }
#'
#' @seealso \code{\link[stats]{optim}}, \code{\link{summary.ls_optim}},
#'   \code{\link{predict.ls_optim}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 100; p <- 3
#' X <- matrix(rnorm(n*p), n, p); colnames(X) <- paste0("x",1:p)
#' beta <- c(2, -1, 0.5)
#' y <- as.vector(cbind(1, X) %*% c(1, beta) + rnorm(n))
#' fit <- my_lm(X, y, B = 50)
#' print(fit)
#' }
#'
#' @importFrom stats optim cov complete.cases
#' @export
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

#' Print method for \code{ls_optim}
#'
#' @param x An object of class \code{"ls_optim"}.
#' @param ... Unused; included for S3 consistency.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print ls_optim
#' @export
print.ls_optim <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nConverged:", x$converged, " |  Bootstrap reps (kept):", nrow(x$betas_boot), "\n")
  invisible(x)
}

#' Coefficients for \code{ls_optim}
#'
#' @param object An \code{"ls_optim"} fit.
#' @param ... Unused.
#' @return Named numeric vector of coefficients.
#'
#' @method coef ls_optim
#' @export
coef.ls_optim <- function(object, ...) object$coefficients

#' Variance-covariance matrix for \code{ls_optim}
#'
#' @param object An \code{"ls_optim"} fit.
#' @param ... Unused.
#' @return Numeric covariance matrix of coefficients (bootstrap estimate).
#'
#' @method vcov ls_optim
#' @export
vcov.ls_optim <- function(object, ...) object$vcov

#' Fitted values for \code{ls_optim}
#'
#' @param object An \code{"ls_optim"} fit.
#' @param ... Unused.
#' @return Numeric vector of fitted values.
#'
#' @method fitted ls_optim
#' @export
fitted.ls_optim <- function(object, ...) object$fitted

#' Residuals for \code{ls_optim}
#'
#' @param object An \code{"ls_optim"} fit.
#' @param ... Unused.
#' @return Numeric vector of residuals.
#'
#' @method residuals ls_optim
#' @export
residuals.ls_optim <- function(object, ...) object$residuals

#' Predict method for \code{ls_optim}
#'
#' @param object An \code{"ls_optim"} fit.
#' @param newdata Optional matrix or data frame of predictors (without intercept).
#'   If omitted, in-sample fitted values are returned.
#' @param ... Unused.
#'
#' @return Numeric vector of predictions.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(50*2), 50, 2); colnames(X) <- c("x1","x2")
#' y <- as.vector(cbind(1, X) %*% c(1, 2, -1) + rnorm(50))
#' fit <- my_lm(X, y, B = 20)
#' predict(fit, newdata = X[1:3, ])
#' }
#'
#' @method predict ls_optim
#' @export
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

#' Summary for \code{ls_optim} fits
#'
#' Produces a summary object including estimates, bootstrap standard errors,
#' test statistics, p-values, and confidence intervals.
#'
#' If the bootstrap covariance is unavailable (e.g., too few successful
#' resamples), a fallback OLS-style variance based on \eqn{\hat\sigma^2 (X^\top X)^{-1}}
#' is used when invertible.
#'
#' @param object An \code{"ls_optim"} fit.
#' @param conf.level Confidence level for two-sided intervals (default \code{0.95}).
#' @param use.t Logical; if \code{TRUE} use \eqn{t}-reference with residual df; otherwise normal.
#' @param ... Unused.
#'
#' @return An object of class \code{"summary.ls_optim"} containing:
#' \itemize{
#' \item \code{call}, \code{coefficients} (matrix with estimate, SE, stat, p, CI),
#' \item \code{df_residual}, \code{conf.level}, \code{use.t}, \code{B}.
#' }
#'
#' @importFrom stats qt pt pnorm qnorm printCoefmat
#' @export
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
  attr(tab, "nboot") <- nrow(object$betas_boot)
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

#' Print method for \code{summary.ls_optim}
#'
#' @param x An object of class \code{"summary.ls_optim"}.
#' @param digits Integer; number of significant digits to print. Defaults to
#'   two fewer than \code{getOption("digits")} but at least 3.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print summary.ls_optim
#' @export
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

#' Marginal plots for \code{ls_optim} fits (base graphics)
#'
#' Produces scatter plots of \code{y} versus each predictor and overlays a
#' “marginal” fitted line for variable \eqn{x_j} defined as
#' \eqn{\alpha + \beta_j x_j + \sum_{k \neq j} \beta_k \bar X_k}, i.e., holding
#' other covariates at their sample means.
#'
#' @param object An \code{"ls_optim"} fit returned by \code{\link{my_lm}}.
#' @param variables Optional subset of variables to plot, as names or integer indices.
#'   Defaults to all non-intercept columns.
#' @param point_alpha Point opacity (0–1) passed to \code{\link[grDevices]{rgb}} for plotting.
#' @param point_cex Point size (passed to \code{\link[graphics]{plot}}).
#' @param line_lwd Line width for the marginal fit lines.
#' @param npoints_line Number of points for the line grid along each predictor.
#'
#' @return \code{NULL} (invisibly). Called for its side effect of plotting.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 80; p <- 2
#' X <- matrix(rnorm(n*p), n, p); colnames(X) <- c("x1","x2")
#' y <- as.vector(cbind(1, X) %*% c(1, 2, -1) + rnorm(n))
#' fit <- my_lm(X, y, B = 30)
#' plot_marginals(fit)
#' }
#'
#' @importFrom graphics par plot lines
#' @importFrom grDevices rgb
#' @export
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
