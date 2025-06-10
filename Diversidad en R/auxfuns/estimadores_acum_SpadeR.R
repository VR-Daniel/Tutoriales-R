# --- Carga de paquetes necesarios ---
paquetes <- c("parallel", "SpadeR")

faltantes <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]

if (length(faltantes) > 0) {
  stop(paste("Faltan paquetes necesarios:", paste(faltantes, collapse = ", ")))
}

# --- Carga de función estimadores_SpadeR desde GitHub (modifica la URL) ---
tryCatch({
  source("D:/Daniel/Desktop/estimadores_SpadeR.R")  # <-- CAMBIA ESTA URL
}, error = function(e) {
  stop("No se pudo cargar la función estimadores_SpadeR desde GitHub")
})


# --- Función principal: estimadores_acum_SpadeR ---
estimadores_acum_SpadeR <- function(DATOS, k = 10, conf = 0.95, perm = 100) {
  z <- qnorm(1 - (1 - conf) / 2)
  nsites <- nrow(DATOS)
  nspp   <- ncol(DATOS)

  # Obtener nombres desde una fila
  nombres <- rownames(estimadores_SpadeR(DATOS[1, , drop = FALSE], k = k, conf = conf)[[1]])
  n_estim <- length(nombres)

  # Inicializar arrays/matrices
  est_array   <- array(NA_real_, dim = c(n_estim, nsites, perm))
  sobs_mat    <- matrix(NA_integer_, nrow = nsites, ncol = perm)
  abun_mat    <- matrix(NA_integer_, nrow = nsites, ncol = perm)
  single_mat  <- matrix(NA_integer_, nrow = nsites, ncol = perm)
  double_mat  <- matrix(NA_integer_, nrow = nsites, ncol = perm)

  # Define worker
  worker <- function(p) {
    datos_perm <- DATOS[sample(nsites), , drop = FALSE]
    datos_acum <- matrix(0L, nrow = nsites, ncol = nspp)
    datos_acum[1, ] <- as.numeric(datos_perm[1, ])
    for (i in 2:nsites) {
      datos_acum[i, ] <- datos_acum[i - 1, ] + as.numeric(datos_perm[i, ])
    }

    res <- estimadores_SpadeR(datos_acum, k = k, conf = conf)
    est_mat <- do.call(cbind, lapply(res, function(x) x[, "Estimado"]))

    list(
      est   = est_mat,
      sobs  = rowSums(datos_acum > 0L),
      abun  = rowSums(datos_acum),
      sing  = rowSums(datos_acum == 1L),
      doub  = rowSums(datos_acum == 2L)
    )
  }

  # Detecta número de núcleos
  ncores <- max(1, parallel::detectCores() - 1)

  # Ejecutar en paralelo según sistema operativo
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl, varlist = c("DATOS", "k", "conf", "nombres", "nsites", "nspp", "estimadores_SpadeR"), envir = environment())
    parallel::clusterEvalQ(cl, library(SpadeR))
    res <- parallel::parLapply(cl, seq_len(perm), worker)
    parallel::stopCluster(cl)
  } else {
    res <- parallel::mclapply(seq_len(perm), worker, mc.cores = ncores)
  }

  # Volcar resultados
  for (p in seq_len(perm)) {
    est_array[, , p] <- res[[p]]$est
    sobs_mat[, p]    <- res[[p]]$sobs
    abun_mat[, p]    <- res[[p]]$abun
    single_mat[, p]  <- res[[p]]$sing
    double_mat[, p]  <- res[[p]]$doub
  }

  # Función auxiliar para resumen
  calcular_estad <- function(m) {
    mu <- rowMeans(m)
    se <- apply(m, 1, sd)
    list(media = mu,
         s.error = se,
         inferior = mu - z * se / sqrt(perm),
         superior = mu + z * se / sqrt(perm))
  }

  # Estadísticas para matrices simples
  ic_sobs   <- calcular_estad(sobs_mat)
  ic_abun   <- calcular_estad(abun_mat)
  ic_single <- calcular_estad(single_mat)
  ic_double <- calcular_estad(double_mat)

  # Estadísticas para estimadores
  media_est   <- apply(est_array, c(1, 2), mean)
  desv_est    <- apply(est_array, c(1, 2), sd)
  ic_low_est  <- media_est - z * desv_est / sqrt(perm)
  ic_high_est <- media_est + z * desv_est / sqrt(perm)

  rownames(media_est)   <- nombres
  rownames(desv_est)    <- nombres
  rownames(ic_low_est)  <- nombres
  rownames(ic_high_est) <- nombres

  # Ensamblar tablas
  promedio <- rbind("Sobs"       = ic_sobs$media,
                    media_est,
                    "Singletons" = ic_single$media,
                    "Doubletons" = ic_double$media,
                    "Abundancia" = ic_abun$media)
  serror <- rbind("Sobs"       = ic_sobs$s.error,
                  desv_est,
                  "Singletons" = ic_single$s.error,
                  "Doubletons" = ic_double$s.error,
                  "Abundancia" = ic_abun$s.error)
  inferior <- rbind("Sobs"       = ic_sobs$inferior,
                    ic_low_est,
                    "Singletons" = ic_single$inferior,
                    "Doubletons" = ic_double$inferior,
                    "Abundancia" = ic_abun$inferior)
  superior <- rbind("Sobs"       = ic_sobs$superior,
                    ic_high_est,
                    "Singletons" = ic_single$superior,
                    "Doubletons" = ic_double$superior,
                    "Abundancia" = ic_abun$superior)

  # Asignar nombres de columnas
  nombres_cols <- as.character(1:nsites)
  colnames(promedio) <- colnames(serror) <- colnames(inferior) <- colnames(superior) <- nombres_cols

  # Reordenar filas
  orden_filas <- c("Sobs", nombres, "Singletons", "Doubletons", "Abundancia")
  promedio  <- promedio[orden_filas, , drop = FALSE]
  serror    <- serror[orden_filas, , drop = FALSE]
  inferior  <- inferior[orden_filas, , drop = FALSE]
  superior  <- superior[orden_filas, , drop = FALSE]

  # Redondear y devolver
  list(
    promedio = as.data.frame(round(promedio, 3)),
    s.error  = as.data.frame(round(serror, 3)),
    inferior = as.data.frame(round(inferior, 3)),
    superior = as.data.frame(round(superior, 3))
  )
}
