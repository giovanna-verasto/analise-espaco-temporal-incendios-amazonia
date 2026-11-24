# =========================================
# 1. PACOTES
# =========================================
library(sf)
library(spatstat)
library(spatstat.geom)
library(readr)
library(dplyr)
library(kernstadapt)
library(ggplot2)

# =========================================
# 2. CARREGAR JANELA E CRIAR OBJETO PPP ESPAÇO-TEMPORAL
# =========================================
ti <- st_read("dados/Poligonos de Terra indigenas", quiet = TRUE)

javari <- ti %>%
  filter(terrai_nom == "Vale do Javari")
javari_utm  <- st_transform(javari, 5880)
javari_utm  <- st_simplify(javari_utm, dTolerance = 500)
javari_owin <- as.owin(st_geometry(javari_utm))

dados_csv <- read_csv("dados/dados_javari_final.csv", show_col_types = FALSE)
dados_sf  <- st_as_sf(
  dados_csv,
  coords = c("x_utm", "y_utm"),
  crs    = 5880
)
coords <- st_coordinates(dados_sf)

ppp_st <- ppp(
  x      = coords[, 1],
  y      = coords[, 2],
  window = javari_owin,
  marks  = data.frame(
    tempo = dados_csv$tempo_continuo,
    ano   = factor(dados_csv$ano)
  )
)

# =========================================
# 3. PARÂMETROS DA GRADE ESPACIAL E VISUALIZAÇÃO (Figura 8.3)
# =========================================
nx <- 4
ny <- 3
nt <- 11 # 11 anos civis (2015 a 2025)

bb      <- Frame(javari_owin)
xbreaks <- seq(bb$xrange[1], bb$xrange[2], length.out = nx + 1)
ybreaks <- seq(bb$yrange[1], bb$yrange[2], length.out = ny + 1)

qc <- quadratcount(unmark(ppp_st), nx = nx, ny = ny)

parse_mid <- function(s) {
  mean(as.numeric(regmatches(s, gregexpr("[-0-9.eE+]+", s))[[1]]))
}
y_mid <- vapply(rownames(qc), parse_mid, numeric(1))
x_mid <- vapply(colnames(qc), parse_mid, numeric(1))

centros <- expand.grid(y = y_mid, x = x_mid)
centros$letra <- LETTERS[seq_len(nrow(centros))]

p_grade <- ggplot() +
  geom_sf(data = javari_utm, fill = NA, color = "black", linewidth = 0.9) +
  geom_vline(xintercept = xbreaks, color = "firebrick", linetype = "dashed") +
  geom_hline(yintercept = ybreaks, color = "firebrick", linetype = "dashed") +
  geom_point(data = dados_csv, aes(x = x_utm, y = y_utm), size = 0.3, alpha = 0.3) +
  geom_text(data = centros, aes(x = x, y = y, label = letra),
            color = "firebrick", fontface = "bold", size = 5) +
  coord_sf() +
  labs(title = paste0("Grade espacial nx=", nx, " × ny=", ny,
                      " sobre o retângulo envolvente"),
       x = NULL, y = NULL) +
  theme_minimal()

print(p_grade)

# =========================================
# 4. TESTE DE SEPARABILIDADE (9.999 SIMULAÇÕES MONTE CARLO)
# =========================================
# Referência: Ghorbani et al. (2021) / González & Moraga (2023)
set.seed(123)
sep_test <- separability.test(
  X     = unmark(ppp_st),
  t     = marks(ppp_st)$tempo,
  nx    = nx,
  ny    = ny,
  nt    = nt,
  nperm = 9999
)
print(sep_test)

# =========================================
# 5. MATRIZ DE RESÍDUOS PADRONIZADOS DE PEARSON (Figura 8.4)
# =========================================
nomes_celulas <- LETTERS[1:(nx * ny)]

celula_x <- cut(ppp_st$x, breaks = xbreaks, include.lowest = TRUE, labels = FALSE)
celula_y <- cut(ppp_st$y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
idx_celula <- (celula_x - 1) * ny + (ny - celula_y + 1)
celula_fator <- factor(nomes_celulas[idx_celula], levels = nomes_celulas)

tab_obs <- table(Celula = celula_fator, Ano = marks(ppp_st)$ano)
teste_chisq <- chisq.test(tab_obs, simulate.p.value = TRUE, B = 9999)

df_diag <- as.data.frame(tab_obs)
names(df_diag) <- c("Celula", "Ano", "Obs")
df_diag$Esp     <- round(as.data.frame(as.table(teste_chisq$expected))$Freq, 2)
df_diag$Residuo <- round(as.data.frame(as.table(teste_chisq$residuals))$Freq, 2)

# Top 12 combinações com maior |Resíduo|
top_residuos <- df_diag %>%
  arrange(desc(abs(Residuo))) %>%
  head(12)
print(top_residuos)

# Heatmap dos resíduos padronizados de Pearson (exibe rótulos para |R_ij| >= 2.5)
df_diag$Celula <- factor(df_diag$Celula, levels = rev(nomes_celulas))

p_residuos <- ggplot(df_diag, aes(x = factor(Ano), y = Celula, fill = Residuo)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(abs(Residuo) >= 2.5, sprintf("%.1f", Residuo), "")),
            size = 3.0, color = "black") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, name = "Resíduo\nde Pearson") +
  labs(title = "Resíduos padronizados de Pearson — dependência espaço x ano",
       subtitle = "Valores numéricos exibidos para |Resíduo| >= 2,5",
       x = "Ano", y = "Célula espacial") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10))

print(p_residuos)
