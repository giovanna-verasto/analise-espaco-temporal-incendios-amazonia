# =========================================
# 1. PACOTES
# =========================================
library(sf)
library(dplyr)
library(readr)
library(lubridate)
library(corrplot)
library(spatstat.geom)

# =========================================
# 2. CARREGAR DADOS JÁ TRATADOS
# =========================================
dados_javari_final <- read_csv("dados/dados_javari_final.csv") %>%
  mutate(
    data_pas = with_tz(data_pas, tzone = "America/Manaus")
  )

# =========================================
# 3. SÉRIE TEMPORAL MENSAL (Figura 4.4)
# =========================================
serie_mensal <- dados_javari_final %>%
  count(ano_mes)

plot(serie_mensal$ano_mes, serie_mensal$n,
     type = "l",
     lwd = 2,
     col = "darkred",
     main = "Focos de calor por mês",
     xlab = "Ano",
     ylab = "Número de focos",
     xaxt = "n")

anos <- unique(year(serie_mensal$ano_mes))
datas_ano <- unique(floor_date(serie_mensal$ano_mes, "year"))

abline(v = datas_ano, col = "gray70", lty = 2)

axis(1,
     at = datas_ano,
     labels = anos,
     las = 2,
     cex.axis = 0.8)

points(serie_mensal$ano_mes, serie_mensal$n,
       pch = 16,
       cex = 0.5,
       col = "darkred")

# =========================================
# 4. GRÁFICOS TEMPORAIS E COVARIÁVEIS (Figuras 4.6, 4.9, 4.10)
# =========================================
x <- dados_javari_final$precipitacao

hist(x[x <= 20],
     breaks = 30,
     col = "lightblue",
     border = "white",
     main = "Distribuição da precipitação (até 20)",
     xlab = "Precipitação")

dados_javari_final %>%
  count(mes) %>%
  with(plot(mes, n, type = "b", pch = 16, col = "steelblue",
            main = "Focos por mês",
            xlab = "Mês", ylab = "Quantidade"))

hist(dados_javari_final$numero_dias_sem_chuva,
     breaks = 30,
     xlim = c(0, 30),
     col = "gray70",
     border = "white",
     main = "Dias sem chuva",
     xlab = "Dias")

dados_javari_final %>%
  count(hora) %>%
  with(plot(hora, n, type = "h", lwd = 3, col = "orange",
            main = "Focos por hora",
            xlab = "Hora", ylab = "Quantidade"))

boxplot(precipitacao ~ mes,
        data = dados_javari_final,
        col = "lightblue",
        border = "darkblue",
        main = "Precipitação por mês",
        xlab = "Mês", ylab = "Precipitação")

dados_javari_final %>%
  mutate(ano = year(data_pas),
         mes = factor(month(data_pas, label = TRUE))) %>%
  count(ano, mes) %>%
  boxplot(n ~ mes, data = .,
          col = "lightcoral",
          border = "darkred",
          main = "Distribuição de focos por mês",
          xlab = "Mês",
          ylab = "Número de focos")

# =========================================
# 5. FRP E RISCO DE FOGO (Figuras 4.12 a 4.15)
# =========================================
frp_log <- log1p(dados_javari_final$frp)

hist(frp_log,
     breaks = 40,
     col = "salmon",
     border = "white",
     main = "Distribuição do log(FRP + 1)",
     xlab = "log(FRP + 1)")

boxplot(frp_log,
        col = "salmon",
        main = "log(FRP + 1)",
        ylab = "log(FRP + 1)")

hist(dados_javari_final$risco_fogo,
     breaks = 30,
     col = "tomato",
     border = "white",
     main = "Distribuição do risco de fogo",
     xlab = "Risco de fogo")

boxplot(risco_fogo ~ mes,
        data = dados_javari_final,
        col = "tomato",
        main = "Risco de fogo por mês",
        xlab = "Mês", ylab = "Risco de fogo")

# =========================================
# 6. CORRELAÇÕES ENTRE COVARIÁVEIS
# =========================================
vars <- dados_javari_final %>%
  select(frp, precipitacao, numero_dias_sem_chuva, risco_fogo)

cor_mat_spearman <- cor(vars, use = "complete.obs", method = "spearman")
cor_mat_pearson  <- cor(vars, use = "complete.obs")

corrplot(cor_mat_spearman,
         method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         number.cex = 0.7,
         title = "Spearman", mar = c(0, 0, 1, 0))

corrplot(cor_mat_pearson,
         method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         number.cex = 0.7,
         title = "Pearson", mar = c(0, 0, 1, 0))

# =========================================
# 7. DADOS ESPACIAIS (POLÍGONO IBGE - SIRGAS 2000 / EPSG:5880)
# =========================================
ti <- st_read("dados/Poligonos de Terra indigenas", quiet = TRUE)

javari <- ti %>%
  filter(terrai_nom == "Vale do Javari")

javari_utm <- st_transform(javari, 5880)

dados_sf <- st_as_sf(
  dados_javari_final,
  coords = c("x_utm", "y_utm"),
  crs = 5880
)

# Mapa geral dos focos (Figura 4.3)
plot(st_geometry(javari_utm),
     col = "gray95", border = "black",
     main = "Focos de incêndio - Vale do Javari")
plot(st_geometry(dados_sf), add = TRUE,
     col = adjustcolor("red", alpha.f = 0.5),
     pch = 16, cex = 0.6)

# =========================================
# 8. VARIÁVEIS AGRUPADAS
# =========================================
dados_sf$periodo <- cut(dados_sf$ano,
                        breaks = c(2015, 2018, 2021, 2025),
                        labels = c("2015-18", "2019-21", "2022-25"),
                        include.lowest = TRUE)

dados_sf$periodo_dia <- cut(dados_sf$hora,
                            breaks = c(-1, 6, 12, 18, 24),
                            labels = c("Madrugada", "Manhã", "Tarde", "Noite"))

dados_sf$estacao <- cut(dados_sf$mes,
                        breaks = c(0, 3, 6, 9, 12),
                        labels = c("Jan-Mar", "Abr-Jun", "Jul-Set", "Out-Dez"))

cores_periodo <- c("#2c7bb6", "#fdae61", "#d7191c")
cores_hora    <- c("#4575b4", "#91cf60", "#fdae61", "#d73027")
cores_estacao <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")

# =========================================
# 9. ANO + MAPA + LEGENDA (Figura 4.7)
# =========================================
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(2, 2, 1))

dados_javari_final %>%
  count(ano) %>%
  with(plot(ano, n, type = "b", pch = 16, col = "darkred",
            main = "Focos por ano",
            xlab = "Ano", ylab = "Quantidade"))

plot(st_geometry(javari_utm),
     col = "gray95", border = "black",
     main = "Distribuição espacial por período")

plot(st_geometry(dados_sf), add = TRUE,
     col = cores_periodo[as.factor(dados_sf$periodo)],
     pch = 16, cex = 0.6)

par(mar = c(5, 1, 4, 1))
plot.new()
legend("center",
       legend = levels(dados_sf$periodo),
       col = cores_periodo,
       pch = 16, title = "Período", bty = "n")
par(mar = c(5, 4, 4, 2))

# =========================================
# 10. HORA + MAPA + LEGENDA (Figura 4.8)
# =========================================
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(2, 2, 1))

dados_javari_final %>%
  count(hora) %>%
  with(plot(hora, n, type = "h", lwd = 3, col = "orange",
            main = "Focos por hora",
            xlab = "Hora", ylab = "Quantidade"))

plot(st_geometry(javari_utm),
     col = "gray95", border = "black",
     main = "Distribuição por período do dia")

plot(st_geometry(dados_sf), add = TRUE,
     col = cores_hora[as.factor(dados_sf$periodo_dia)],
     pch = 16, cex = 0.6)

par(mar = c(5, 1, 4, 1))
plot.new()
legend("center",
       legend = levels(dados_sf$periodo_dia),
       col = cores_hora,
       pch = 16, title = "Período do dia", bty = "n")
par(mar = c(5, 4, 4, 2))

# =========================================
# 11. MÊS + MAPA + LEGENDA (Figura 4.5)
# =========================================
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(2, 2, 1))

dados_javari_final %>%
  count(mes) %>%
  with(plot(mes, n, type = "b", pch = 16, col = "steelblue",
            main = "Focos por mês",
            xlab = "Mês", ylab = "Quantidade"))

plot(st_geometry(javari_utm),
     col = "gray95", border = "black",
     main = "Distribuição por períodos do ano")

plot(st_geometry(dados_sf), add = TRUE,
     col = cores_estacao[as.factor(dados_sf$estacao)],
     pch = 16, cex = 0.6)

par(mar = c(5, 1, 4, 1))
plot.new()
legend("center",
       legend = levels(dados_sf$estacao),
       col = cores_estacao,
       pch = 16, title = "Período do ano", bty = "n")
par(mar = c(5, 4, 4, 2))

layout(1)

# =========================================
# 12. MAPAS DE COVARIÁVEIS E OBJETOS PPP MARCADOS (Figura 4.11)
# =========================================
plotar_mapa_covariavel <- function(val, titulo, paleta_cores) {
  breaks  <- pretty(val[!is.na(val)], n = 5)
  classes <- cut(val, breaks = breaks, include.lowest = TRUE)
  paleta  <- colorRampPalette(paleta_cores)(length(levels(classes)))
  cores   <- paleta[as.numeric(classes)]

  plot(st_geometry(javari_utm),
       col = "gray95", border = "black", main = titulo)
  plot(st_geometry(dados_sf), add = TRUE,
       col = cores, pch = 16, cex = 0.6)
  legend("bottomleft",
         legend = levels(classes),
         fill = paleta, title = titulo, bty = "n")
}

plotar_mapa_covariavel(dados_sf$frp,
                       "Intensidade do fogo (FRP)",
                       c("#ffcccc", "#990000"))

plotar_mapa_covariavel(dados_sf$numero_dias_sem_chuva,
                       "Dias sem chuva",
                       c("#ffe6cc", "#cc5500"))

plotar_mapa_covariavel(dados_sf$risco_fogo,
                       "Risco de fogo",
                       c("#ffe5e5", "#990000"))

# Objetos ppp marcados
coords <- st_coordinates(dados_sf)
janela <- as.owin(st_bbox(javari_utm))
dados_sf$risco_fogo[dados_sf$risco_fogo < 0] <- NA

criar_ppp <- function(valores) {
  ok <- !is.na(valores)
  ppp(x = coords[ok, 1], y = coords[ok, 2],
      window = janela, marks = valores[ok])
}

ppp_frp   <- criar_ppp(dados_sf$frp)
ppp_dias  <- criar_ppp(dados_sf$numero_dias_sem_chuva)
ppp_prec  <- criar_ppp(dados_sf$precipitacao)
ppp_risco <- criar_ppp(dados_sf$risco_fogo)

# Plot final de risco de fogo com marcas proporcionais (Figura 4.11)
plot(ppp_risco,
     use.marks = TRUE,
     cols = hcl.colors(100, palette = "Reds", rev = TRUE),
     main = "Risco de fogo",
     symbolmap = symbolmap("circles", size = 0.015))
plot(st_geometry(javari_utm), add = TRUE, border = "black", lwd = 1.5)
