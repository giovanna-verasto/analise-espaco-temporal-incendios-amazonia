# =========================================
# 1. PACOTES
# =========================================
# Base metodológica: González & Moraga (2023), The R Journal (RJ-2023-025)
library(sf)
library(readr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(stpp)
library(KernSmooth)
library(plot3D)

# =========================================
# 2. PREPARAÇÃO DOS DADOS ESPAÇO-TEMPORAIS
# =========================================
ti <- st_read("dados/Poligonos de Terra indigenas", quiet = TRUE)
javari <- filter(ti, terrai_nom == "Vale do Javari")
javari <- st_transform(javari, 5880)
javari <- st_simplify(javari, dTolerance = 500)
janela_A <- as.owin(st_geometry(javari))

dados <- read_csv("dados/dados_javari_final.csv", show_col_types = FALSE)
dados <- dados[is.finite(dados$x_utm) & is.finite(dados$y_utm) &
                 is.finite(dados$tempo_continuo), , drop = FALSE]

ppp_st <- ppp(
  x      = dados$x_utm,
  y      = dados$y_utm,
  window = janela_A,
  marks  = dados$tempo_continuo,
  check  = TRUE
)

tempos   <- as.numeric(marks(ppp_st))
X        <- unmark(ppp_st)
n        <- npoints(X)
W        <- Window(X)
janela_t <- c(0, 132) # janeiro/2015 a janeiro/2026

cat("Número de eventos:", n, "\n")
cat("Intervalo temporal:", round(min(tempos), 2), "a", round(max(tempos), 2), "meses\n")

# =========================================
# 3. INTENSIDADE ESPAÇO-TEMPORAL NÃO SEPARÁVEL (Figura 8.5)
# =========================================
h_s_km <- 29.87551   # critério de Turner (2009)
h_s_m  <- h_s_km * 1000
h_t <- bw.nrd0(tempos)   # regra de Silverman (1986), como em González e Moraga (2023)

corte_temporal <- 3
resolucao      <- 256
paleta_tg      <- hcl.colors(120, "YlOrRd", rev = TRUE)

pasta_saida <- "resultados_espaco_temporal"
dir.create(pasta_saida, showWarnings = FALSE, recursive = TRUE)

anos <- 2015:2025
tempos_referencia <- 8 + 12 * (anos - 2015) # 1º de setembro de cada ano

estimar_snapshot_ht <- function(t_ref, ht_val, res_grid = resolucao, corte = corte_temporal) {
  et_val <- pnorm((janela_t[2] - tempos) / ht_val) -
    pnorm((janela_t[1] - tempos) / ht_val)
  idx <- which(abs(tempos - t_ref) <= corte * ht_val)

  pesos_t <- dnorm(tempos[idx], mean = t_ref, sd = ht_val) / et_val[idx]
  mapa_m2 <- density.ppp(
    X[idx], sigma = h_s_m, weights = pesos_t, kernel = "gaussian",
    edge = TRUE, diggle = TRUE, at = "pixels", dimyx = res_grid,
    positive = TRUE
  )
  mapa_m2 * 1e6 # eventos/m²/mês -> eventos/km²/mês
}

mapas_setembro <- lapply(tempos_referencia, estimar_snapshot_ht, ht_val = h_t)
names(mapas_setembro) <- as.character(anos)

z_max <- max(vapply(mapas_setembro, function(im) max(im$v, na.rm = TRUE), numeric(1)))

png(file.path(pasta_saida, "fig_8_5_intensidade_01_setembro.png"),
    width = 2600, height = 1900, res = 220)
op <- par(mfrow = c(3, 4), mar = c(1.5, 1.5, 2.8, 4.2))
for (ano in anos) {
  plot(mapas_setembro[[as.character(ano)]],
       main = paste("01/09/", ano),
       col = paleta_tg, zlim = c(0, z_max))
  plot(W, add = TRUE, border = "black", lwd = 0.8)
}
plot.new()
par(op)
dev.off()

# -----------------------------------------
# 3.1 Análise de sensibilidade da largura de banda temporal h_t (Figura 8.6)
# -----------------------------------------
ht_candidatos <- c(2.0, h_t, 14.0)
names(ht_candidatos) <- c("ht_2.0m", "ht_7.18m", "ht_14.0m")

anos_comp   <- c(2017, 2018, 2023, 2024)
t_refs_comp <- 8 + 12 * (anos_comp - 2015)

mapas_comp <- list()
zmax_comp  <- 0
for (nome_ht in names(ht_candidatos)) {
  ht_val <- ht_candidatos[[nome_ht]]
  mapas_comp[[nome_ht]] <- lapply(t_refs_comp, function(t_ref) {
    im <- estimar_snapshot_ht(t_ref, ht_val, res_grid = 200, corte = 4)
    zmax_comp <<- max(zmax_comp, max(im$v, na.rm = TRUE))
    im
  })
}

png(file.path(pasta_saida, "fig_8_6_sensibilidade_ht_intensidade.png"),
    width = 2600, height = 1800, res = 220)
op <- par(mfrow = c(length(ht_candidatos), length(anos_comp)),
          mar = c(1.5, 1.5, 2.8, 3.8))
for (nome_ht in names(ht_candidatos)) {
  ht_val <- ht_candidatos[[nome_ht]]
  for (j in seq_along(anos_comp)) {
    plot(mapas_comp[[nome_ht]][[j]],
         main = paste0("01/09/", anos_comp[j], " | h_t = ", round(ht_val, 2), "m"),
         col = paleta_tg, zlim = c(0, zmax_comp))
    plot(W, add = TRUE, border = "black", lwd = 0.8)
  }
}
par(op)
dev.off()

# =========================================
# 4. ESTIMAÇÃO DE LAMBDA(s_i, t_i) NOS EVENTOS OBSERVADOS
# =========================================
calcular_lambda_nos_eventos <- function(X, tempos, hs_m, ht_meses,
                                        janela_t, corte = 3) {
  n_eventos <- npoints(X)
  et <- pnorm((janela_t[2] - tempos) / ht_meses) -
    pnorm((janela_t[1] - tempos) / ht_meses)

  lambda <- rep(NA_real_, n_eventos)
  tempos_unicos <- sort(unique(tempos))

  for (k in seq_along(tempos_unicos)) {
    t_ref <- tempos_unicos[k]
    idx_contrib <- if (is.infinite(corte)) seq_len(n_eventos) else which(abs(tempos - t_ref) <= corte * ht_meses)
    idx_alvo <- which(tempos == t_ref)

    pesos_t <- dnorm(tempos[idx_contrib], mean = t_ref, sd = ht_meses) / et[idx_contrib]

    dens_pontos <- density.ppp(
      X[idx_contrib], sigma = hs_m, weights = pesos_t,
      kernel = "gaussian", edge = TRUE, diggle = TRUE,
      at = "points", positive = TRUE
    )
    lambda[idx_alvo] <- as.numeric(dens_pontos)[match(idx_alvo, idx_contrib)] * 1e6
  }
  lambda
}

lambda_st_eventos <- calcular_lambda_nos_eventos(
  X, tempos, h_s_m, h_t, janela_t, corte_temporal
)

# =========================================
# 5. PARÂMETROS E ESTIMAÇÃO DE SEGUNDA ORDEM (stpp)
# =========================================
FMD      <- as.3dpoints(X$x / 1000, X$y / 1000, tempos)
s.region <- cbind(W$bdry[[1]]$x / 1000, W$bdry[[1]]$y / 1000)
t.region <- janela_t

# Seleção automática das larguras de banda de segunda ordem
sigmaD <- bw.scott(X)
MP     <- density.ppp(X, diggle = TRUE, positive = TRUE, sigma = sigmaD, at = "points")
bwG    <- bw.pcf(X, cv.method = "compLik", divisor = "d", lambda = MP)
hs_pcf_km <- as.numeric(bwG)[1] / 1000

dt     <- as.numeric(dist(sort(unique(tempos))))
ht_pcf <- 2 * dpik(dt, kernel = "epanech", gridsize = 50, scalest = "iqr")

intensidade_media <- n / area.owin(W)
rmax_m <- rmax.rule("K", W, intensidade_media)
r0_km  <- 0.7 * rmax_m / 1000
tau0   <- 0.15 * max(dt)
ngrid  <- 71

dist_grid <- seq(0, r0_km, length.out = ngrid)
tau_grid  <- seq(0, tau0,  length.out = ngrid)

# -----------------------------------------
# 5.1 Função K espaço-temporal não homogênea centrada D(r, tau)
# -----------------------------------------
K_st <- STIKhat(
  xyt      = FMD,
  s.region = s.region,
  t.region = t.region,
  lambda   = lambda_st_eventos,
  dist     = dist_grid,
  times    = tau_grid,
  infectious = TRUE # unilateral: 0 < t_j - t_i <= tau
)

# D(r, tau) = K_inhom(r, tau) - pi * r^2 * tau
ref_cstr   <- outer(K_st$dist^2, K_st$times, function(r2, tau) pi * r2 * tau)
K_centrado <- K_st$Khat - ref_cstr

dist_K_plot <- K_st$dist
tau_K_plot  <- K_st$times

# -----------------------------------------
# 5.2 Função de correlação de pares espaço-temporal g(r, tau)
# -----------------------------------------
g_st <- PCFhat(
  xyt      = FMD,
  s.region = s.region,
  t.region = t.region,
  lambda   = lambda_st_eventos,
  dist     = dist_grid,
  times    = tau_grid,
  ks       = "epanech",
  hs       = hs_pcf_km,
  kt       = "epanech",
  ht       = ht_pcf
)

dist_g <- g_st$dist
tau_g  <- g_st$times

# =========================================
# 6. VISUALIZAÇÃO DE D(r, tau) (Figura 8.7)
# =========================================
zmax_D <- as.numeric(quantile(K_centrado, probs = 0.99, na.rm = TRUE))
zmin_D <- -0.15 * zmax_D
zlim_K <- c(zmin_D, zmax_D)
K_cor  <- pmax(pmin(K_centrado, zlim_K[2]), zlim_K[1])

png(file.path(pasta_saida, "fig_8_7_K_espaco_temporal_painel.png"),
    width = 3200, height = 1400, res = 220)
op <- par(mfrow = c(1, 2), mar = c(4.5, 5.5, 3.5, 5.2))

persp3D(
  x = dist_K_plot, y = tau_K_plot, z = K_centrado,
  colvar = K_cor, clim = zlim_K,
  facets = NA, curtain = FALSE,
  col = paleta_tg,
  colkey = list(width = 0.8, clab = "D(r,tau)"),
  bty = "g", theta = -30, phi = 30, border = NA,
  ticktype = "detailed", cex.axis = 0.75,
  xlab = "r (km)", ylab = "tau (meses)", zlab = "",
  main = "(a) Superfície 3D de D(r, tau)"
)

image2D(
  x = dist_K_plot, y = tau_K_plot, z = K_cor,
  col = paleta_tg, zlim = zlim_K,
  colkey = list(width = 0.8, clab = "D(r,tau)"),
  xlab = "Distância espacial r (km)",
  ylab = "Defasagem temporal tau (meses)",
  main = "(b) Mapa de calor 2D de D(r, tau) e isolinha 0"
)
contour(
  x = dist_K_plot, y = tau_K_plot, z = K_centrado,
  levels = 0, add = TRUE, col = "black", lwd = 2, lty = 2, labcex = 0.9
)
par(op)
dev.off()

# =========================================
# 7. VISUALIZAÇÃO DE g(r, tau) COMPLETA E TRUNCADA EM 5 (Figuras 8.8 e 8.9)
# =========================================
# -----------------------------------------
# 7.1 Escala completa (Figura 8.8)
# -----------------------------------------
zlim_g_full <- c(0, max(g_st$pcf, na.rm = TRUE))

png(file.path(pasta_saida, "fig_8_8_g_espaco_temporal_completo.png"),
    width = 3200, height = 1400, res = 220)
op <- par(mfrow = c(1, 2), mar = c(4.5, 5.5, 3.5, 5.2))

persp3D(
  x = dist_g, y = tau_g, z = g_st$pcf,
  facets = NA, curtain = FALSE,
  col = paleta_tg, zlim = zlim_g_full, clim = zlim_g_full,
  colkey = list(width = 0.8, clab = "g(r,tau)"),
  bty = "g", theta = -30, phi = 30, border = NA,
  ticktype = "detailed", cex.axis = 0.8,
  xlab = "r (km)", ylab = "tau (meses)", zlab = "",
  main = "(a) Superfície 3D completa de g(r, tau)"
)

image2D(
  x = dist_g, y = tau_g, z = g_st$pcf,
  col = paleta_tg, zlim = zlim_g_full,
  colkey = list(width = 0.8, clab = "g(r,tau)"),
  xlab = "Distância espacial r (km)",
  ylab = "Defasagem temporal tau (meses)",
  main = "(b) Mapa de calor 2D completo e isolinha 1"
)
contour(
  x = dist_g, y = tau_g, z = g_st$pcf,
  levels = 1, add = TRUE, col = "black", lwd = 2, lty = 2, labcex = 0.9
)
par(op)
dev.off()

# -----------------------------------------
# 7.2 Zoom truncado em g(r, tau) = 5 (Figura 8.9)
# -----------------------------------------
teto_zoom_g <- 5
zlim_g_zoom <- c(0, teto_zoom_g)
g_zoom      <- pmin(g_st$pcf, teto_zoom_g)

png(file.path(pasta_saida, "fig_8_9_g_espaco_temporal_zoom.png"),
    width = 3200, height = 1400, res = 220)
op <- par(mfrow = c(1, 2), mar = c(4.5, 5.5, 3.5, 5.2))

persp3D(
  x = dist_g, y = tau_g, z = g_zoom,
  facets = NA, curtain = FALSE,
  col = paleta_tg, zlim = zlim_g_zoom, clim = zlim_g_zoom,
  colkey = list(width = 0.8, clab = "g(r,tau)"),
  bty = "g", theta = -30, phi = 30, border = NA,
  ticktype = "detailed", cex.axis = 0.8,
  xlab = "r (km)", ylab = "tau (meses)", zlab = "",
  main = paste0("(a) Zoom 3D de g(r, tau) (truncado em ", teto_zoom_g, ")")
)

image2D(
  x = dist_g, y = tau_g, z = g_zoom,
  col = paleta_tg, zlim = zlim_g_zoom,
  colkey = list(width = 0.8, clab = "g(r,tau)"),
  xlab = "Distância espacial r (km)",
  ylab = "Defasagem temporal tau (meses)",
  main = paste0("(b) Zoom 2D (teto = ", teto_zoom_g, ") e isolinha 1")
)
contour(
  x = dist_g, y = tau_g, z = g_st$pcf,
  levels = 1, add = TRUE, col = "black", lwd = 2, lty = 2, labcex = 0.9
)
par(op)
dev.off()

# =========================================
# 8. SENSIBILIDADE DE g(r, tau) À LARGURA DE BANDA TEMPORAL h_t (Figura 8.10)
# =========================================
g_comp_list <- lapply(ht_candidatos, function(ht_val) {
  lam_ht <- calcular_lambda_nos_eventos(
    X, tempos, h_s_m, ht_val, janela_t, corte = 4
  )
  PCFhat(
    xyt      = FMD,
    s.region = s.region,
    t.region = t.region,
    lambda   = lam_ht,
    dist     = dist_grid,
    times    = tau_grid,
    ks       = "epanech",
    hs       = hs_pcf_km,
    kt       = "epanech",
    ht       = ht_pcf
  )$pcf
})

png(file.path(pasta_saida, "fig_8_10_sensibilidade_ht_g.png"),
    width = 3000, height = 1000, res = 220)
op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 5.0))
for (i in seq_along(ht_candidatos)) {
  g_comp_zoom <- pmin(g_comp_list[[i]], teto_zoom_g)
  image2D(
    x = dist_g, y = tau_g, z = g_comp_zoom,
    col = paleta_tg, zlim = c(0, teto_zoom_g),
    colkey = list(width = 0.8, clab = "g(r,tau)"),
    xlab = "r (km)", ylab = "tau (meses)",
    main = paste0("g(r, tau) com h_t = ", round(ht_candidatos[i], 2), " meses")
  )
  contour(x = dist_g, y = tau_g, z = g_comp_list[[i]],
          levels = 1, add = TRUE, col = "black", lwd = 1.8, lty = 2)
}
par(op)
dev.off()
