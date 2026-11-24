# =========================================
# 1. PACOTES
# =========================================
library(dplyr)
library(readr)
library(lubridate)

# =========================================
# 2. CARREGAR BASE DE DADOS BRUTA
# =========================================
# Ajuste o caminho caso os arquivos estejam em outra pasta local
dados <- read_csv("dados/dados_javari_utm.csv")

# -----------------------------------------
# 2.1 Remover NA e duplicatas exatas
# -----------------------------------------
dados <- dados %>%
  filter(!is.na(data_pas))

dados_unique <- dados %>%
  distinct(x_utm, y_utm, data_pas, .keep_all = TRUE)

# -----------------------------------------
# 2.2 Converter fuso: UTC -> America/Manaus
#     (Vale do Javari está em UTC-4)
# -----------------------------------------
dados_unique <- dados_unique %>%
  mutate(
    data_pas = with_tz(data_pas, tzone = "America/Manaus")
  )

# -----------------------------------------
# 2.3 Filtro temporal (diferença > 24h por local)
#     Mantém apenas a primeira ocorrência de um mesmo evento
# -----------------------------------------
dados_time <- dados_unique %>%
  arrange(x_utm, y_utm, data_pas) %>%
  group_by(x_utm, y_utm) %>%
  mutate(
    diff_tempo = as.numeric(difftime(data_pas, lag(data_pas), units = "hours"))
  )

dados_final <- dados_time %>%
  mutate(
    manter = is.na(diff_tempo) | diff_tempo > 24
  ) %>%
  filter(manter) %>%
  ungroup()

# -----------------------------------------
# 2.4 Variável contínua de tempo
#     t0 = 1 de janeiro do primeiro ano dos dados (2015)
# -----------------------------------------
ano_inicial <- year(min(dados_final$data_pas))
t0 <- ymd(paste0(ano_inicial, "-01-01"), tz = "America/Manaus")

dados_final <- dados_final %>%
  mutate(
    tempo_continuo = time_length(interval(t0, data_pas), "month"),
    ano     = year(data_pas),
    mes     = month(data_pas),
    hora    = hour(data_pas),
    ano_mes = floor_date(data_pas, "month")
  )

# -----------------------------------------
# 2.5 Tratar valores inconsistentes em risco de fogo
# -----------------------------------------
dados_final <- dados_final %>%
  mutate(
    risco_fogo = ifelse(risco_fogo == -999, NA, risco_fogo)
  )

# =========================================
# 3. SALVAR BASE TRATADA (n = 1.975)
# =========================================
write_csv(dados_final, "dados/dados_javari_final.csv")
