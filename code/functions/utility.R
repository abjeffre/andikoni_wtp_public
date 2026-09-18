###################################################################
############ SHARED HELPERS FOR THE WTP REPLICATION ###############
# Sourced by every MAKE script after its library() calls.
# Assumes the working directory is the repository root.

# Test mode: WTP_TEST=1 runs every model with a handful of iterations so the
# whole pipeline can be exercised in minutes. Unset it for the paper's settings.
WTP_TEST <- nzchar(Sys.getenv("WTP_TEST"))
ITER_FULL <- 250
ITER_TEST <- 20
SAMPLER_ITER <- if (WTP_TEST) ITER_TEST else ITER_FULL
SAMPLER_CHAINS <- 4
SAMPLER_REFRESH <- if (WTP_TEST) 0 else 100

# The SI sweep is a lambda x N grid; test mode keeps only the corners.
SWEEP_LAMBDA_FULL <- seq(0, 0.9, by = 0.1)
SWEEP_N_FULL <- c(250, 500, 1000)
SWEEP_LAMBDA_TEST <- c(0, 0.9)
SWEEP_N_TEST <- c(250)

# Output directories are not all in git; every step writes into them.
dir.create("figures", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

###################################################################
############ POSTERIOR HELPERS ####################################

create_empty_like <- function(obj) {
  if (is.vector(obj)) {
    return(vector(mode = typeof(obj), length = length(obj)))
  } else if (is.matrix(obj)) {
    return(matrix(NA, nrow = nrow(obj), ncol = ncol(obj)))
  } else if (is.data.frame(obj)) {
    return(obj[0, ])
  } else if (is.array(obj)) {
    return(array(NA, dim = dim(obj)))
  } else if (is.list(obj)) {
    empty_list <- vector("list", length(obj))
    names(empty_list) <- names(obj)
    return(empty_list)
  } else {
    stop("Unsupported object type")
  }
}

dims <- function(x) {
  z <- dim(x)
  if (is.null(z)) {
    return(1)
  }
  return(length(z))
}

get_init_list <- function(pf) {
  stanfit <- posterior::as_draws_rvars(pf)
  init_list <- list()
  cnt <- 1
  for (i in names(stanfit)) {
    means <- posterior::summarise_draws(posterior::as_draws_array(stanfit[[i]]))
    temp <- create_empty_like(stanfit[[i]])
    for (j in 1:length(means$mean)) {
      temp[j] <- c(means$mean[j])
    }
    init_list[[cnt]] <- as.array(temp)
    cnt <- cnt + 1
  }
  names(init_list) <- names(stanfit)
  for (i in names(init_list)) {
    if (dims(init_list[[i]]) == 1 & length(init_list[[i]]) == 1) {
      init_list[[i]] <- init_list[[i]][1]
    }
  }
  return(init_list)
}

# Returns a plain named list of posterior draws, including log_lik. The MAKE
# steps save this rather than the cmdstanr fit object: it is far smaller and
# carries everything the figures and the model comparison need.
extract.samples2 <- function(x) {
  output <- list()
  stanfit <- posterior::as_draws_rvars(x)
  for (i in names(stanfit)) {
    output[[i]] <- posterior::draws_of(stanfit[[i]])
  }
  return(output)
}

# One place for the sampler settings, so test mode reaches every fit.
fit_model <- function(stan_file, stan_data, seed = 1) {
  model <- cmdstanr::cmdstan_model(stan_file)
  model$sample(
    data = stan_data,
    chains = SAMPLER_CHAINS,
    parallel_chains = SAMPLER_CHAINS,
    iter_sampling = SAMPLER_ITER,
    iter_warmup = SAMPLER_ITER,
    refresh = SAMPLER_REFRESH,
    seed = seed
  )
}

###################################################################
############ RESIDUAL AND COVARIANCE STRUCTURE ####################

# Posterior draws of the residual covariance implied by the unrestricted MVN
# model: Sigma = (diag(tau) L)(diag(tau) L)'.
get_cov_mvn <- function(post) {
  S <- nrow(post$tau)
  cov_list <- vector("list", S)
  for (s in 1:S) {
    L <- post$L_Omega[s, , ]
    tau <- post$tau[s, ]
    L_Sigma <- diag(tau) %*% L
    cov_list[[s]] <- L_Sigma %*% t(L_Sigma)
  }
  cov_list
}

# The same object for the one-factor model, where the loadings (-lambda, +lambda)
# induce Sigma = lambda^2 * [[1,-1],[-1,1]] + diag(sigma^2).
get_cov_latent <- function(post) {
  S <- length(post$lambda)
  cov_list <- vector("list", S)
  for (s in 1:S) {
    lambda <- post$lambda[s]
    Sigma_lat <- lambda^2 * matrix(c(1, -1, -1, 1), 2, 2)
    Sigma_noise <- diag(c(post$sigma[s, 1]^2, post$sigma[s, 2]^2))
    cov_list[[s]] <- Sigma_lat + Sigma_noise
  }
  cov_list
}

# Share of residual covariance carried by the dominant eigenvalue.
rank1_post <- function(cov_list) {
  sapply(cov_list, function(Sigma) {
    ev <- eigen(Sigma, symmetric = TRUE)$values
    ev[1] / sum(ev)
  })
}

rank1_share <- function(r1, r2) {
  sv <- svd(cov(cbind(r1, r2)))$d
  sv[1] / sum(sv)
}

# Posterior-mean residuals for each outcome.
#
# The "mvn" branch rebuilds base_model.stan's linear predictor in full:
#   mu[o][i] = a[o] + bS[o, sex[i]] + bW[o]*wealth[i]
#              + bA[o]*age[i] + bA2[o]*age[i]^2 + edu_effect[o][edu[i]]
# edu_effect is a transformed parameter of the Stan program (a monotonic
# effect built from an ordered simplex), so it is read from the posterior
# rather than recomputed here.
extract_mu_residuals <- function(post, df, model_type = c("latent", "mvn")) {
  model_type <- match.arg(model_type)

  if (model_type == "latent") {
    mu_p <- apply(post$mu_priv_hat, 2, mean)
    mu_c <- apply(post$mu_pub_hat, 2, mean)
    return(list(
      res_p = df$priv - mu_p,
      res_c = df$pub - mu_c,
      T_hat = apply(post$T, 2, mean)
    ))
  }

  a <- apply(post$a, 2, mean)
  bW <- apply(post$bW, 2, mean)
  bA <- apply(post$bA, 2, mean)
  bA2 <- apply(post$bA2, 2, mean)
  bS <- apply(post$bS, c(2, 3), mean)                 # outcome x sex
  edu_effect <- apply(post$edu_effect, c(2, 3), mean) # outcome x education level

  mu_o <- function(o) {
    a[o] + bS[o, df$sex] + bW[o] * df$wealth +
      bA[o] * df$age + bA2[o] * df$age^2 + edu_effect[o, df$edu]
  }

  list(res_p = df$priv - mu_o(1), res_c = df$pub - mu_o(2), T_hat = NULL)
}

leading_eigvec <- function(cov_list) {
  lapply(cov_list, function(Sigma) {
    v <- eigen(Sigma, symmetric = TRUE)$vectors[, 1]
    if (v[2] < 0) v <- -v
    v
  })
}

###################################################################
############ PLOTTING #############################################

# Covariance ellipse with its principal axis, drawn over the residuals.
plot_cov_with_data <- function(Sigma, x, y, main, col_ellipse = "black") {
  e <- eigen(Sigma)
  theta <- seq(0, 2 * pi, length.out = 300)
  ellipse <- e$vectors %*% diag(sqrt(e$values)) %*% rbind(cos(theta), sin(theta))

  plot(y, x,
    pch = 16, col = rgb(0, 0, 0, 0.25),
    ylab = "Private WTP", xlab = "Community WTP", main = main, asp = 1
  )
  lines(ellipse[1, ], ellipse[2, ], lwd = 2, col = col_ellipse)
  arrows(0, 0,
    e$vectors[1, 1] * sqrt(e$values[1]), e$vectors[2, 1] * sqrt(e$values[1]),
    col = "red", lwd = 2, length = 0.1
  )
  abline(h = 0, v = 0, lty = 3, col = "grey80")
}

# dens2(), my_dotchart() and bookcase() below are the figure primitives used by
# 2_MAKE_figures.R. dens2 is a variant of rethinking::dens that shades an HPDI
# and marks zero; my_dotchart is a coefficient plot that preserves input order.

dens2 <- function(
  x, adj = .5, main = "", show.HPDI = FALSE, show.zero = FALSE,
  rm.na = TRUE, add = FALSE, col_fill = col.alpha("#005555", .2),
  norm.comp = FALSE, ...
) {
  if (inherits(x, "data.frame")) {
    n <- ncol(x)
    cnames <- colnames(x)
    set_nice_margins()
    par(mfrow = make.grid(n))
    for (i in 1:n) {
      dens(x[, i],
        adj = adj, norm.comp = norm.comp, show.HPDI = show.HPDI,
        show.zero = TRUE, xlab = cnames[i], ...
      )
    }
  } else {
    if (rm.na == TRUE) {
      x <- x[!is.na(x)]
    }
    thed <- density(x, adjust = adj)
    if (add == FALSE) {
      set_nice_margins()
      plot(thed, main = main, ...)
    } else {
      lines(thed$x, thed$y, ...)
    }
    if (show.HPDI != FALSE) {
      hpd <- HPDI(x, prob = show.HPDI)
      rethinking::shade(thed, hpd, col = col_fill)
    }
    if (norm.comp == TRUE) {
      mu <- mean(x)
      sigma <- sd(x)
      curve(dnorm(x, mu, sigma),
        col = "white", lwd = 2,
        add = TRUE
      )
      curve(dnorm(x, mu, sigma), add = TRUE)
    }
    if (show.zero == TRUE) {
      lines(c(0, 0), c(0, max(thed$y) * 2), lty = 2)
    }
  }
}
my_dotchart <- function(
  x,
  error = NULL,
  labels = NULL,
  xlab = "Value",
  main = "",
  ylab = "",
  pch = 19,
  col = "black",
  buffer = 0.5,
  err.col = "gray",
  err.lwd = 1,
  err.length = 0.05,
  label.cex = 1,
  xlab.cex = 1,
  ...
) {
  # If labels are not provided, try to use names from x; otherwise, use sequence numbers.
  if (is.null(labels)) {
    labels <- names(x)
    if (is.null(labels)) {
      labels <- as.character(seq_along(x))
    }
  }
  # Preserve the original order of x.
  ord <- seq_along(x)
  x_ordered <- x[ord]
  labels_ordered <- labels[ord]
  # Create y-positions corresponding to each observation in the original order.
  n <- length(x_ordered)
  y_positions <- seq(from = 1, to = n)
  # Set up the plot with additional y-axis buffer at the top and bottom.
  # cex.lab controls the size of the x- and y-axis labels (here we target xlab via xlab.cex).
  plot(x_ordered,
    y_positions,
    type = "n",
    xlab = xlab,
    ylab = ylab,
    main = main,
    yaxt = "n",
    cex.lab = xlab.cex,
    ylim = c(min(y_positions) - buffer, max(y_positions) + buffer),
    ...
  )
  # Optionally add horizontal gridlines for clarity.
  abline(h = y_positions, col = "gray90", lty = "dotted")
  # If error data is provided, plot error bars.
  if (!is.null(error)) {
    error <- as.matrix(error)
    if (ncol(error) != length(x)) {
      stop("The number of columns in 'error' must match the length of 'x'.")
    }
    # Determine the lower and upper error bounds.
    if (!is.null(rownames(error))) {
      lower_index <- which(rownames(error) %in% c("5%", "2.5%", "lower"))
      upper_index <- which(rownames(error) %in% c("94%", "97.5%", "upper"))
      if (length(lower_index) == 1 && length(upper_index) == 1) {
        lower_vals <- error[lower_index, ]
        upper_vals <- error[upper_index, ]
      } else {
        lower_vals <- error[1, ]
        upper_vals <- error[2, ]
      }
    } else {
      lower_vals <- error[1, ]
      upper_vals <- error[2, ]
    }
    # Preserve the original order for the error bounds.
    lower_vals_ordered <- lower_vals[ord]
    upper_vals_ordered <- upper_vals[ord]
    # Draw horizontal error bars.
    arrows(
      x0 = lower_vals_ordered,
      y0 = y_positions,
      x1 = upper_vals_ordered,
      y1 = y_positions,
      angle = 90, code = 3,
      length = err.length,
      col = err.col, lwd = err.lwd
    )
  }
  # Plot the points.
  points(x_ordered, y_positions, pch = pch, col = col)
  # Add a custom y-axis with the original order labels; cex.axis controls the tick label size.
  axis(side = 2, at = y_positions, labels = labels_ordered, las = 1, cex.axis = label.cex)
  # Return (invisibly) the ordered x-values, labels, and y positions.
  invisible(list(x_ordered = x_ordered, labels_ordered = labels_ordered, y_positions = y_positions))
}
bookcase <- function(x) {
  x <- gsub("_", " ", x)
  x <- tolower(trimws(x))
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
