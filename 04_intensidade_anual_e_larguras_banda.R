# =========================================
# 1. PACOTES
# =========================================
library(sf)
library(spatstat)
library(spatstat.geom)
library(readr)
library(dplyr)
library(ggplot2)

# =========================================
# 2. CONFIGURAÇÃO E CARREGAMENTO DOS DADOS
# =========================================
nx <- 4
ny <- 3
min_pontos <- 5
q_zmax <- 1
KM2 <- 1e6   # conversão de focos/m² para focos/km²

ti <- st_read("dados/Poligonos de Terra indigenas", quiet = TRUE)
javari <- ti %>% filter(terrai_nom == "Vale do Javari")
javari_utm  <- st_transform(javari, 5880)
javari_utm  <- st_simplify(javari_utm, dTolerance = 500)
javari_owin <- as.owin(st_geometry(javari_utm))

dados_csv <- read_csv("dados/dados_javari_final.csv", show_col_types = FALSE)
dados_sf  <- st_as_sf(dados_csv, coords = c("x_utm", "y_utm"), crs = 5880)
coords    <- st_coordinates(dados_sf)

ppp_st <- ppp(
  x      = coords[, 1],
  y      = coords[, 2],
  window = javari_owin,
  marks  = data.frame(
    tempo = dados_csv$tempo_continuo,
    ano   = factor(dados_csv$ano)
  )
)

anos    <- levels(marks(ppp_st)$ano)
n_anos  <- length(anos)
bb_glob <- Frame(javari_owin)
X_geral <- unmark(ppp_st)

cores_kernel <- hcl.colors(100, "YlOrRd", rev = TRUE)

# =========================================
# 3. CANDIDATOS DE LARGURA DE BANDA ESPACIAL (sigma)
# =========================================
sigma_da_particao <- function(janela, nx, ny, modo = c("geom", "media", "min")) {
  modo <- match.arg(modo)
  bb <- Frame(janela)
  largura_x <- diff(bb$xrange) / nx
  largura_y <- diff(bb$yrange) / ny
  switch(modo,
    geom  = sqrt(largura_x * largura_y),
    media = mean(c(largura_x, largura_y)),
    min   = min(largura_x, largura_y)
  )
}

tentar <- function(expr) {
  tryCatch(as.numeric(expr), error = function(e) {
    message("Seletor falhou: ", conditionMessage(e)); NA_real_
  })
}

escala_W      <- sqrt(area.owin(javari_owin))
w_celula_geom <- sigma_da_particao(javari_owin, nx, ny, "geom")
w_celula_min  <- sigma_da_particao(javari_owin, nx, ny, "min")

candidatos <- c(
  diggle             = tentar(bw.diggle(X_geral)),
  ppl                = tentar(bw.ppl(X_geral)),
  cvl                = tentar(bw.CvL(X_geral)),
  scott              = tentar(bw.scott.iso(X_geral)),
  turner_0025        = 0.025 * escala_W,
  turner_005         = 0.05  * escala_W,
  turner_010         = 0.10  * escala_W,
  particao_equiv     = w_celula_geom / sqrt(12),
  particao_equiv_min = w_celula_min  / sqrt(12)
)

tab_sigma <- data.frame(
  metodo       = names(candidatos),
  sigma_km     = round(candidatos / 1000, 1),
  fator_turner = round(candidatos / escala_W, 3)
)
tab_sigma <- tab_sigma[order(tab_sigma$sigma_km), ]
rownames(tab_sigma) <- NULL
print(tab_sigma)

# =========================================
# 4. FUNÇÕES AUXILIARES DE ESTIMAÇÃO E PLOTAGEM
# =========================================
calcular_kernels_ano <- function(ppp_obj, sigma, min_pontos = 5) {
  Tj <- split(ppp_obj %mark% marks(ppp_obj)$ano)
  Tj <- Tj[sapply(Tj, npoints) >= min_pontos]
  lapply(Tj, function(a) {
    density(unmark(a), sigma = sigma, edge = TRUE, kernel = "gaussian") * KM2
  })
}

zmax_comum <- function(lista_kernels, q_zmax = 1) {
  valores <- unlist(lapply(lista_kernels, function(k) as.vector(k$v)))
  as.numeric(quantile(valores, probs = q_zmax, na.rm = TRUE))
}

camadas_base <- function(poligono_sf, nx, ny) {
  xbreaks <- seq(bb_glob$xrange[1], bb_glob$xrange[2], length.out = nx + 1)
  ybreaks <- seq(bb_glob$yrange[1], bb_glob$yrange[2], length.out = ny + 1)
  list(
    geom_sf(data = poligono_sf, fill = NA, color = "black",
            linewidth = 0.7, inherit.aes = FALSE),
    geom_vline(xintercept = xbreaks, color = "black",
               linetype = "dashed", linewidth = 0.25),
    geom_hline(yintercept = ybreaks, color = "black",
               linetype = "dashed", linewidth = 0.25),
    coord_sf(xlim = bb_glob$xrange, ylim = bb_glob$yrange, expand = FALSE),
    theme_minimal(),
    theme(axis.text = element_blank(), panel.grid = element_blank())
  )
}

escala_fill <- function(zmax, nome = "focos/km²") {
  scale_fill_gradientn(colors = cores_kernel, limits = c(0, zmax),
                       oob = scales::squish, name = nome)
}

painel_kernels_ano <- function(kernels, poligono_sf, nx, ny, zmax, titulo) {
  df_long <- bind_rows(lapply(names(kernels), function(a) {
    d <- as.data.frame(kernels[[a]]); d$ano <- a; d
  }))
  ggplot() +
    geom_raster(data = df_long, aes(x = x, y = y, fill = value)) +
    escala_fill(zmax) +
    camadas_base(poligono_sf, nx, ny) +
    facet_wrap(~ano) +
    labs(title = titulo, x = NULL, y = NULL)
}

# =========================================
# 5. COMPARAÇÃO ENTRE LARGURAS DE BANDA ESPACIAIS
# =========================================
metodos_testar <- c("ppl", "turner_0025", "turner_005", "turner_010", "particao_equiv")
metodos_testar <- metodos_testar[!is.na(candidatos[metodos_testar])]

df_rel <- bind_rows(lapply(metodos_testar, function(m) {
  k <- density(X_geral, sigma = candidatos[[m]], edge = TRUE, kernel = "gaussian") * KM2
  d <- as.data.frame(k / max(k$v, na.rm = TRUE))
  d$metodo <- paste0(m, " (", round(candidatos[[m]] / 1000, 1), " km)")
  d
}))
df_rel$metodo <- factor(df_rel$metodo, levels = unique(df_rel$metodo))

p_comp_forma <- ggplot() +
  geom_raster(data = df_rel, aes(x = x, y = y, fill = value)) +
  scale_fill_gradientn(colors = cores_kernel, limits = c(0, 1),
                       name = "intensidade\nrelativa") +
  camadas_base(javari_utm, nx, ny) +
  facet_wrap(~metodo) +
  labs(title = "Kernel geral, intensidade relativa por método", x = NULL, y = NULL)

print(p_comp_forma)

# =========================================
# 6. INTENSIDADE POR ANO (ESCALA COMUM E NORMALIZADA) - turner_010 (29,9 km)
# =========================================
sigma_escolhido <- candidatos[["turner_010"]]
sufixo <- paste0("turner_010 | sigma = ", round(sigma_escolhido / 1000, 1), " km")

# -----------------------------------------
# 6.1 Kernel por ano em escala comum (Figura 8.1)
# -----------------------------------------
kernels_ano <- calcular_kernels_ano(ppp_st, sigma_escolhido, min_pontos)
zmax_abs    <- zmax_comum(kernels_ano, q_zmax)

p_anual <- painel_kernels_ano(
  kernels_ano, javari_utm, nx, ny, zmax_abs,
  titulo = paste0("Kernel por ano (escala comum) — ", sufixo)
)
print(p_anual)

# -----------------------------------------
# 6.2 Kernel por ano normalizado (cada ano soma 1) (Figura 8.2)
# -----------------------------------------
kernels_norm <- lapply(kernels_ano, function(k) k / integral(k) * KM2)
zmax_norm    <- zmax_comum(kernels_norm, q_zmax)

df_norm_long <- bind_rows(lapply(names(kernels_norm), function(a) {
  d <- as.data.frame(kernels_norm[[a]]); d$ano <- a; d
}))

p_normalizado <- ggplot() +
  geom_raster(data = df_norm_long, aes(x = x, y = y, fill = value)) +
  scale_fill_gradientn(colors = cores_kernel, limits = c(0, zmax_norm),
                       oob = scales::squish,
                       name = "fração dos\nfocos do ano\npor km²") +
  camadas_base(javari_utm, nx, ny) +
  facet_wrap(~ano) +
  labs(title = paste0("Kernel por ano normalizado (cada ano soma 1) — ", sufixo),
       x = NULL, y = NULL)

print(p_normalizado)
