# =========================================
# 1. PACOTES
# =========================================
library(sf)
library(spatstat)
library(spatstat.geom)
library(readr)
library(dplyr)

# =========================================
# 2. JANELA DE OBSERVAÇÃO (POLÍGONO REAL)
# =========================================
ti <- st_read("dados/Poligonos de Terra indigenas", quiet = TRUE)

javari <- ti %>%
  filter(terrai_nom == "Vale do Javari")
javari_utm <- st_transform(javari, 5880)

# Simplificar o polígono para estabilidade computacional
javari_utm  <- st_simplify(javari_utm, dTolerance = 500)
javari_owin <- as.owin(st_geometry(javari_utm))

dados_csv <- read_csv("dados/dados_javari_final.csv", show_col_types = FALSE)
dados_sf  <- st_as_sf(
  dados_csv,
  coords = c("x_utm", "y_utm"),
  crs    = 5880
)

# =========================================
# 3. OBJETO PPP ESPACIAL
# =========================================
coords <- st_coordinates(dados_sf)

ppp_focos <- ppp(
  x      = coords[, 1],
  y      = coords[, 2],
  window = javari_owin
)

# =========================================
# 4. INTENSIDADE ESPACIAL POR KERNEL (Figura 5.1)
# =========================================
# Critério de Turner (2009): sigma = 0.10 * sqrt(|A|) -> h_s = 29,9 km
escala_W <- sqrt(area.owin(javari_owin))
h_turner <- 0.10 * escala_W

kernel_focos <- density(
  ppp_focos,
  sigma  = h_turner,
  edge   = TRUE,
  kernel = "gaussian"
) * 1e6 # conversão de focos/m² para focos/km²

cores <- hcl.colors(100, "YlOrRd", rev = TRUE)

plot(kernel_focos,
     col  = cores,
     main = paste0("turner_010 (", round(h_turner / 1000, 1), " km)"))

plot(st_geometry(javari_utm),
     add    = TRUE,
     border = "black",
     lwd    = 1.5)

# =========================================
# 5. FUNÇÃO L CENTRALIZADA E ENVELOPES DE MONTE CARLO (Figura 5.2)
# =========================================
set.seed(123)
env_focos <- envelope(
  ppp_focos, Lest,
  nsim       = 999,
  nrank      = 25,
  correction = "best",
  rmax       = 50000
)

plot(env_focos, . - r ~ r,
     main   = "Função L - Focos de calor",
     xlab   = "r (metros)",
     ylab   = "L(r) - r",
     legend = FALSE)

# =========================================
# 6. FUNÇÃO DE CORRELAÇÃO DE PARES g(r) (Figura 5.3)
# =========================================
pcf_focos <- pcf(
  ppp_focos,
  correction = "best",
  rmax       = 50000
)

r_cruzamento <- pcf_focos$r[which(pcf_focos$iso < 1)[1]]
cat("g(r) cruza 1 em aproximadamente:", round(r_cruzamento), "metros\n")

# (a) Visão geral até 50 km
plot(pcf_focos,
     main   = "",
     xlab   = "r (metros)",
     ylab   = "g(r)",
     legend = FALSE)

# (b) Detalhe até 25 km
plot(pcf_focos,
     main   = "",
     xlab   = "r (metros)",
     ylab   = "g(r)",
     legend = FALSE,
     xlim   = c(0, 25000),
     ylim   = c(0, 10))

# =========================================
# 7. ANÁLISE DE ISOTROPIA - FUNÇÃO K DIRECIONAL (Figura 5.4)
# =========================================
K_0   <- Ksector(ppp_focos, begin = 0,   end = 45,  rmax = 50000)
K_45  <- Ksector(ppp_focos, begin = 45,  end = 90,  rmax = 50000)
K_90  <- Ksector(ppp_focos, begin = 90,  end = 135, rmax = 50000)
K_135 <- Ksector(ppp_focos, begin = 135, end = 180, rmax = 50000)

plot(K_0$r, K_0$iso,
     type = "l", lwd = 2,
     col  = "blue",
     main = "K direcional - Focos de calor",
     xlab = "r (metros)",
     ylab = "K(r)",
     ylim = range(c(K_0$iso, K_45$iso, K_90$iso, K_135$iso)))

lines(K_45$r,  K_45$iso,  col = "red",       lwd = 2, lty = 2)
lines(K_90$r,  K_90$iso,  col = "darkgreen", lwd = 2, lty = 3)
lines(K_135$r, K_135$iso, col = "orange",    lwd = 2, lty = 4)

legend("topleft",
       legend = c("0-45°", "45-90°", "90-135°", "135-180°"),
       col    = c("blue", "red", "darkgreen", "orange"),
       lwd    = 2,
       lty    = 1:4,
       bty    = "n")
