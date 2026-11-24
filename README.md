# Análise Estatística de Focos de Incêndio na Amazônia Brasileira como Padrões Pontuais Espaciais e Espaço-Temporais

Repositório contendo os scripts em **R** desenvolvidos para o Trabalho de Conclusão de Curso (TCC) em Estatística pela **Universidade Federal de São Carlos (UFSCar)**.

---

## Sobre o Projeto

Este trabalho analisa a dinâmica espacial e espaço-temporal dos focos de calor registrados na **Terra Indígena Vale do Javari** (Amazonas) entre **fevereiro de 2015 e dezembro de 2025** ($n = 1.975$ eventos válidos), sob a perspectiva da teoria de **Processos Pontuais Espaciais e Espaço-Temporais**.

As análises estão estruturadas nas seguintes etapas metodológicas:

* **Tratamento e padronização geoespacial:** conversão de fuso horário (UTC para horário local `America/Manaus`), projeção métrica SIRGAS 2000 / Brasil Policônica (`EPSG:5880`), filtragem de eventos redundantes em menos de 24 horas e construção da escala temporal contínua em meses.
* **Análise exploratória:** séries temporais mensais e anuais, distribuição horária, mapas temáticos e análise das covariáveis ambientais (precipitação acumulada, dias consecutivos sem chuva, índice de risco de fogo e potência radiativa do fogo — FRP).
* **Análise puramente espacial:** estimação de densidade por kernel gaussiano com correção de borda ($h_s = 29,9\text{ km}$, critério de Turner, 2009), função $K$ de Ripley e função $L(r) - r$ com envelopes de Monte Carlo (999 simulações sob CSR), função de correlação de pares $g(r)$ e análise de isotropia via função $K$ direcional.
* **Análise espaço-temporal:** estimação de intensidade espacial por ano civil e comparação de larguras de banda; teste formal de separabilidade de primeira ordem via simulações de Monte Carlo ($B = 9.999$) e diagnóstico por resíduos padronizados de Pearson; estimação da superfície de intensidade conjunta não separável $\hat{\lambda}(s, t)$ com análise de sensibilidade para $h_t$; e investigação da estrutura de dependência de segunda ordem via função $K$ espaço-temporal não homogênea centrada $\hat{D}(r, \tau)$ e função de correlação de pares $\hat{g}(r, \tau)$.

---

## Estrutura do Repositório

Os scripts estão numerados na ordem sequencial de execução e correspondem diretamente aos capítulos da monografia:

* `01_tratamento_dados.R` — Limpeza da base de focos de calor, remoção de duplicatas, conversão de fuso horário (`America/Manaus`), filtro temporal de 24 horas, tratamento de valores ausentes (`-999` em risco de fogo) e criação da variável `tempo_continuo` (Seção 4.3).
* `02_analise_descritiva.R` — Gráficos temporais (mensal, anual e horário), histogramas e boxplots das covariáveis, matrizes de correlação e mapas exploratórios sobre o polígono da Terra Indígena Vale do Javari (Seção 4.4, Figuras 4.3 a 4.15).
* `03_analise_espacial.R` — Estimação da intensidade espacial global por kernel ($h_s = 29,9\text{ km}$), função $L(r) - r$ com 999 simulações de Monte Carlo, função de correlação de pares $g(r)$ e função $K$ direcional em quatro setores de $45^\circ$ (Capítulo 5, Figuras 5.1 a 5.4).
* `04_intensidade_anual_e_larguras_banda.R` — Avaliação comparativa de seletores de largura de banda espacial ($h_s$) e estimação dos mapas de intensidade por ano civil (2015–2025), tanto em escala absoluta comum quanto normalizados para integrar a 1 (Seção 8.1, Figuras 8.1 e 8.2).
* `05_teste_separabilidade_e_residuos.R` — Construção da grade espacial regular $4 \times 3$ (células A a L), aplicação do teste de separabilidade de primeira ordem com 9.999 simulações de Monte Carlo (`kernstadapt::separability.test`) e geração da matriz/heatmap de resíduos padronizados de Pearson com destaque para $|R_{ij}| \ge 2,5$ (Seção 8.2, Figuras 8.3 e 8.4).
* `06_analise_espaco_temporal.R` — Estimação da superfície de intensidade não separável $\hat{\lambda}(s, t)$ avaliada em $1^\circ$ de setembro de cada ano ($h_s = 29,88\text{ km}$, $h_t = 7,18\text{ meses}$), análise de sensibilidade para $h_t \in \{2,0;\; 7,18;\; 14,0\}\text{ meses}$, cálculo da função $\hat{D}(r, \tau)$ e estimação da função $\hat{g}(r, \tau)$ em escala completa e truncada em $\hat{g}(r, \tau) = 5$ (Seções 8.3 a 8.5, Figuras 8.5 a 8.10).

---

## Fonte dos Dados e Reprodutibilidade

* **Focos de calor:** Instituto Nacional de Pesquisas Espaciais (INPE) — Programa Queimadas (período de 2015 a 2025).
* **Malha territorial:** Instituto Brasileiro de Geografia e Estatística (IBGE) — Bases Cartográficas Contínuas (polígono oficial da Terra Indígena Vale do Javari).

Para reproduzir os resultados localmente, crie uma pasta `dados/` na raiz do repositório contendo o arquivo tabular tratado (`dados_javari_final.csv`) e o shapefile dos polígonos de Terras Indígenas (`Poligonos de Terra indigenas/`), e execute os scripts na ordem numérica (`01` a `06`).

---

## Pacotes Utilizados

```r
install.packages(c(
  "dplyr", "readr", "lubridate", "ggplot2", "corrplot",
  "sf", "spatstat", "spatstat.geom", "spatstat.explore",
  "kernstadapt", "stpp", "KernSmooth", "plot3D"
))
```

---

## Referências Metodológicas Principais

* Baddeley, A., Rubak, E., & Turner, R. (2016). *Spatial Point Patterns: Methodology and Applications with R*. CRC Press.
* Diggle, P. J. (2014). *Statistical Analysis of Spatial and Spatio-Temporal Point Patterns* (3rd ed.). CRC Press.
* Ghorbani, M., Vafaei, N., Dvořák, J., & Myllymäki, M. (2021). Testing the first-order separability hypothesis for spatio-temporal point patterns. *Computational Statistics & Data Analysis*, 161, 107245.
* González, J. A., & Moraga, P. (2023). Non-parametric analysis of spatial and spatio-temporal point patterns. *The R Journal*, 15(1), 65–82.
* Turner, R. (2009). Point patterns of forest fire locations. *Environmental and Ecological Statistics*, 16, 197–223.

---

## Autoria

* **Giovanna Martins Verasto** — Bacharelanda em Estatística (UFSCar)
* **Orientação:** Prof. Dr. Danilo Lourenço Lopes
