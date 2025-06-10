estimadores_SpadeR <- function(DATOS, k = 10, conf = 0.95) {
  z <- -qnorm((1 - conf) / 2)
  
  resultados <- lapply(seq_len(nrow(DATOS)), function(i) {
    x <- as.numeric(DATOS[i, ])
    x <- x[x > 0]
    n <- sum(x)
    D <- length(x)
    f <- function(i) sum(x == i)
    f1 <- f(1)
    f2 <- f(2)
    
    # Datos básicos para ACE
    x_rare <- x[x <= k]
    n_rare <- sum(x_rare)
    D_rare <- length(x_rare)
    D_abun <- sum(x > k)
    C_rare <- ifelse(n_rare > 0, 1 - f1 / n_rare, 1)
    
    j <- 1:k
    a1 <- sum(j * (j - 1) * sapply(j, f))
    a2 <- sum(j * sapply(j, f))
    
    gamma_hat <- if (a2 > 1 && C_rare != 0) {
      max(D_rare / C_rare * a1 / a2 / (a2 - 1) - 1, 0)
    } else 0
    
    gamma_bc <- if (C_rare != 0 && a2 > 1) {
      max(gamma_hat * (1 + (1 - C_rare) / C_rare * a1 / (a2 - 1)), 0)
    } else 0
    
    # ---- Chao1 ----
    if (f1 > 0 && f2 > 0) {
      S_Chao1 <- D + (n - 1) / n * f1^2 / (2 * f2)
      var_Chao1 <- f2 * ((n - 1) / n * (f1 / f2)^2 / 2 +
                         ((n - 1) / n)^2 * (f1 / f2)^3 +
                         ((n - 1) / n)^2 * (f1 / f2)^4 / 4)
    } else if (f1 > 1 && f2 == 0) {
      S_Chao1 <- D + (n - 1) / n * f1 * (f1 - 1) / 2
      var_Chao1 <- (n - 1) / n * f1 * (f1 - 1) / 2 +
                   ((n - 1) / n)^2 * f1 * (2 * f1 - 1)^2 / 4 -
                   ((n - 1) / n)^2 * f1^4 / (4 * S_Chao1)
    } else {
      S_Chao1 <- D
      var_Chao1 <- 0
    }
    t1 <- S_Chao1 - D
    K1 <- exp(z * sqrt(log(1 + var_Chao1 / t1^2)))
    CI_Chao1 <- c(D + t1 / K1, D + t1 * K1)

    # ---- Chao1-bc ----
    S_Chao1_bc <- D + (n - 1) / n * f1 * (f1 - 1) / (2 * (f2 + 1))
    if (f2 > 0) {
      var_Chao1_bc <- (n - 1)/n * f1 * (f1 - 1)/2/(f2 + 1) +
        ((n - 1)/n)^2 * f1 * (2 * f1 - 1)^2 / (4 * (f2 + 1)^2) +
        ((n - 1)/n)^2 * f1^2 * f2 * (f1 - 1)^2 / (4 * (f2 + 1)^4)
    } else {
      var_Chao1_bc <- (n - 1)/n * f1 * (f1 - 1)/2 +
        ((n - 1)/n)^2 * f1 * (2 * f1 - 1)^2 / 4 -
        ((n - 1)/n)^2 * f1^4 / (4 * S_Chao1_bc)
    }
    t2 <- S_Chao1_bc - D
    K2 <- exp(z * sqrt(log(1 + var_Chao1_bc / t2^2)))
    CI_Chao1_bc <- c(D + t2 / K2, D + t2 * K2)

    # ---- ACE ----
    S_ACE <- D_abun + D_rare / C_rare + f1 / C_rare * gamma_hat
    s_ace <- S_ACE
    diff_ace <- function(q) {
      if (q == 1) {
        return((1 - f1 / n_rare + D_rare * (n_rare - f1) / n_rare^2) / C_rare^2)
      }
      return(1)
    }
    COV.f <- function(i, j) {
      if (i == j) f(i) * (1 - f(i) / s_ace) else -f(i) * f(j) / s_ace
    }
    i <- rep(unique(x), each = length(unique(x)))
    j <- rep(unique(x), length(unique(x)))
    var_ace <- sum(mapply(function(i, j) diff_ace(i) * diff_ace(j) * COV.f(i, j), i, j))
    t_ace <- s_ace - D
    K_ace <- exp(z * sqrt(log(1 + var_ace / t_ace^2)))
    CI_ACE <- c(D + t_ace / K_ace, D + t_ace * K_ace)

    # ---- ACE-bc ----
    S_ACE_bc <- D_abun + D_rare / C_rare + f1 / C_rare * gamma_bc
    t_ace_bc <- S_ACE_bc - D
    K_ace_bc <- exp(z * sqrt(log(1 + var_ace / t_ace_bc^2)))
    CI_ACE_bc <- c(D + t_ace_bc / K_ace_bc, D + t_ace_bc * K_ace_bc)

    # ---- Jackknife 1 ----
    jack1 <- D + (n - 1) / n * f1
    var_jack1 <- f1 * (1 - f1 / jack1) * (1 + (n - 1) / n)^2
    K_j1 <- exp(z * sqrt(log(1 + var_jack1 / (jack1 - D)^2)))
    CI_j1 <- c(D + (jack1 - D) / K_j1, D + (jack1 - D) * K_j1)

    # ---- Jackknife 2 ----
    jack2 <- D + (2 * n - 3)/n * f1 - ((n - 2)^2)/(n * (n - 1)) * f2
    var_jack2 <- f1 + 4 * f2
    K_j2 <- exp(z * sqrt(log(1 + var_jack2 / (jack2 - D)^2)))
    CI_j2 <- c(D + (jack2 - D) / K_j2, D + (jack2 - D) * K_j2)

    # Tabla final con redondeo
    out <- rbind(
      "Chao1"     = c(S_Chao1, sqrt(var_Chao1), CI_Chao1),
      "Chao1-bc"  = c(S_Chao1_bc, sqrt(var_Chao1_bc), CI_Chao1_bc),
      "ACE"       = c(S_ACE, sqrt(var_ace), CI_ACE),
      "ACE-bc"    = c(S_ACE_bc, sqrt(var_ace), CI_ACE_bc),
      "Jack1"     = c(jack1, sqrt(var_jack1), CI_j1),
      "Jack2"     = c(jack2, sqrt(var_jack2), CI_j2)
    )
    
    colnames(out) <- c("Estimado", "s.error", "Inferior", "Superior")
    round(out, 3)
  })

  names(resultados) <- rownames(DATOS)
  return(resultados)
}