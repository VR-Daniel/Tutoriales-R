# Definimos una semilla (para la reproducibilidad del código)
set.seed(42)

# Vector de abundancias para la Comunidad 1
C1.Jan <- c(rep(1, 70), rep(2, 17), rep(3, 4), rep(4, 5), rep(5, 5), rep(6, 5), rep(7, 5), rep(8, 3), rep(9, 1), rep(10, 2), rep(11, 3), rep(12, 2), rep(14, 2), rep(17, 1), rep(19, 2), rep(20, 3), rep(21, 1), rep(24, 1), rep(26, 1), rep(40, 1), rep(57, 2), rep(60, 1), rep(64, 1), rep(71, 1), rep(77, 1))

# Aleatorizamos el orden de los datos
C1.Jan <- sample(C1.Jan, 140)

# Vector de abundancias para la Comunidad 2
C2.Jan <- c(rep(1, 84), rep(2, 10), rep(3, 4), rep(4, 3), rep(5, 5), rep(6, 1), rep(7, 2), rep(8, 1), rep(14, 1), rep(42, 1), rep(0, 28))

# Aleatorizamos el orden de los datos
C2.Jan <- sample(C2.Jan, 140)

# Añadimos ambos vectores a una tabla
Janzen73 <- as.data.frame(rbind(C1 = C1.Jan, C2 = C2.Jan))

# Eliminamos las variables intermedias
rm(C1.Jan, C2.Jan)