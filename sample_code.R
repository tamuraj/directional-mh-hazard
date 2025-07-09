  # ============================================================
  # Directional Measure φ and Its Asymptotic Confidence Interval
  # ============================================================
  # This R script defines a directional measure φ for square r×r 
  # contingency tables, its analytical gradient, and computes 
  # its asymptotic confidence interval via the delta method.
  # ============================================================
  
  #' Calculate the directional measure φ from a square contingency table
  #' @param table Square contingency table (r × r)
  #' @return Scalar φ value
  calculate_phi <- function(table) {
    r  <- nrow(table)
    ta <- table / sum(table)  # Normalize to probability matrix
    
    p_row <- rowSums(ta)
    p_col <- colSums(ta)
    
    omega_x <- omega_y <- numeric(r)
    omega_x[1] <- p_row[1]
    omega_y[1] <- p_col[1]
    for (i in 2:r) {
      omega_x[i] <- p_row[i] / (1 - sum(p_row[1:(i - 1)]))
      omega_y[i] <- p_col[i] / (1 - sum(p_col[1:(i - 1)]))
    }
    omega_x <- omega_x[-r]
    omega_y <- omega_y[-r]
    
    W1 <- omega_x * (1 - omega_y)
    W2 <- omega_y * (1 - omega_x)
    
    Delta <- sum(W1 + W2)
    W1s   <- W1 / Delta
    W2s   <- W2 / Delta
    theta <- acos(W1 / sqrt(W1^2 + W2^2))
    
    4 / pi * sum((W1s + W2s) * (theta - pi / 4))
  }
  
  #' Compute the gradient vector of φ with respect to π_ab
  #' @param p A square probability matrix (r × r)
  #' @return Gradient vector (length r^2, column-major)
  phi_grad_theory <- function(p) {
    r  <- nrow(p)
    rs <- rowSums(p)
    cs <- colSums(p)
    Fx <- c(0, cumsum(rs))[1:r]
    Fy <- c(0, cumsum(cs))[1:r]
    
    omegaX <- rs[1:(r - 1)] / (1 - Fx[1:(r - 1)])
    omegaY <- cs[1:(r - 1)] / (1 - Fy[1:(r - 1)])
    W1 <- omegaX * (1 - omegaY)
    W2 <- omegaY * (1 - omegaX)
    Delta <- sum(W1 + W2)
    Wst1 <- W1 / Delta
    Wst2 <- W2 / Delta
    theta <- acos(W1 / sqrt(W1^2 + W2^2))
    
    grad <- matrix(0, r, r)
    for (a in 1:r) for (b in 1:r) {
      dW1 <- dW2 <- numeric(r - 1)
      for (i in 1:(r - 1)) {
        denomX <- 1 - Fx[i]
        denomY <- 1 - Fy[i]
        dOmX <- (a == i) / denomX + omegaX[i] * (a < i) / denomX
        dOmY <- (b == i) / denomY + omegaY[i] * (b < i) / denomY
        dW1[i] <- (1 - omegaY[i]) * dOmX - omegaX[i] * dOmY
        dW2[i] <- (1 - omegaX[i]) * dOmY - omegaY[i] * dOmX
      }
      dDelta <- sum(dW1 + dW2)
      
      d_phi <- 0
      for (i in 1:(r - 1)) {
        dWst1 <- (dW1[i] - Wst1[i] * dDelta) / Delta
        dWst2 <- (dW2[i] - Wst2[i] * dDelta) / Delta
        dtheta <- (W1[i] * dW2[i] - W2[i] * dW1[i]) / (W1[i]^2 + W2[i]^2)
        d_phi <- d_phi +
          (theta[i] - pi / 4) * (dWst1 + dWst2) +
          (Wst1[i] + Wst2[i]) * dtheta
      }
      grad[b, a] <- 4 / pi * d_phi  # column-major gradient
    }
    as.vector(grad)
  }
  
  # ============================================================
  # Example: Estimating φ and its asymptotic 95% confidence interval
  # ============================================================
  
  # Load required package
  library(Matrix)
  
  # Set parameters
  set.seed(123)
  n <- 500                          # Sample size
  conf.z <- qnorm(0.975)            # z-value for 95% CI
  
  # True probability matrix (4×4 example)
  p <- matrix(
    c(0.30, 0.05, 0.03, 0.02,
      0.04, 0.25, 0.06, 0.05,
      0.01, 0.05, 0.20, 0.06,
      0.02, 0.03, 0.07, 0.15),
    nrow = 4, byrow = TRUE
  )
  
  # Generate a sample from multinomial distribution
  a <- rmultinom(1, n, prob = c(t(p)))
  pa <- matrix(a, 4, 4, byrow = TRUE) / n  # Empirical frequencies
  
  # Estimate φ from observed sample
  phi_hat <- calculate_phi(pa)
  
  # Compute gradient vector (row-major format)
  grad_vec <- c(t(matrix(phi_grad_theory(pa), 4, 4)))
  
  # Multinomial covariance matrix for vectorized table (row-major)
  Sigma <- diag(c(pa)) - tcrossprod(c(pa))  # 16 × 16
  
  # Delta method: standard error
  se <- sqrt(as.numeric(t(grad_vec) %*% Sigma %*% grad_vec) / n)
  
  # 95% confidence interval
  ci <- phi_hat + c(-1, 1) * conf.z * se
  
  # Output results
  cat(sprintf("Estimated φ: %.4f\n", phi_hat))
  cat(sprintf("Standard Error: %.4f\n", se))
  cat(sprintf("95%% Confidence Interval: (%.4f, %.4f)\n", ci[1], ci[2]))
