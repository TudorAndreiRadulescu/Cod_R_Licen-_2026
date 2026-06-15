# Documentație Detaliată — Pipeline FRM Stocuri Bancare
## Structura completă a codului și fluxul de date între fișiere

---

## ARHITECTURA GENERALĂ

```
FRM_Stable_Banci.R  ← FIȘIERUL PRINCIPAL (2.781 linii, conține toate blocurile inline)
        │
        ├── source("FRM_Statistics_Algorithm.R")  ← algoritmul matematic LASSO
        ├── source("FRM_SC_config.R")              ← configurație (la fiecare bloc)
        └── source("FRM_SC_utils.R")              ← helpers vizuali (la fiecare bloc)
```

**Regula de bază:** `FRM_Stable_Banci.R` nu soursiază celelalte fișiere de pipeline
(`FRM_SC_load_data.R`, `FRM_SC_estimation_varying.R` etc.) — le conține **inline** ca blocuri
numerotate. Fișierele separate există pentru rulare independentă.

**Ordinea obligatorie de source:** întotdeauna `config` ÎNAINTE de `utils`, deoarece
`config` conține `rm(list = ls(all = TRUE))` care șterge tot mediul, inclusiv
`theme_transparent_bottom` dacă `utils` a fost deja soursat.

---

## FIȘIER 1: `FRM_Statistics_Algorithm.R`
### Rol: Algoritmul matematic de bază (nu se modifică niciodată)

**Autori originali:** Youjuan Li & Ji Zhu, University of Michigan (2006),
adaptat de Andrija Mihoci pentru Financial Risk Meter.

**Funcții exportate:**

### `FRM_Quantile_Regression(Xd, j, a, max.steps, eps1, eps2)`

**Parametri de intrare:**
| Parametru | Tip | Descriere |
|---|---|---|
| `Xd` | matrix n×m | Matricea de date: n=obs (90 zile), m=variabile (7 bănci + macro) |
| `j` | integer | Indexul coloanei care devine variabila dependentă (Y) |
| `a` | numeric | Cuantila: 0.05 în pipeline-ul tău (coada de 5%) |
| `max.steps` | integer | Numărul maxim de iterații LASSO: I=25 |
| `eps1` | numeric | Toleranță numerică internă (implicit 1e-10) |
| `eps2` | numeric | Prag minim pentru lambda (implicit 1e-6) |

**Ce face intern:**
1. Separă Y = coloana j, X = restul coloanelor
2. Calculează SVD pentru condiționarea numerică (`FRM_Condition`)
3. Perturbă ușor Y dacă există valori duplicate (evită degenerarea)
4. Rulează algoritmul path LASSO: pornind de la λ=∞ (zero variabile active),
   adaugă câte o variabilă pe rând reducând penalizarea
5. La fiecare pas k calculează GACV și SIC pentru selecția modelului

**Output returnat (listă):**
| Câmp | Descriere |
|---|---|
| `beta` | Matrice (max.steps+1) × (m-1): coeficienții la fiecare pas k |
| `beta0` | Vector: interceptul la fiecare pas k |
| `lambda` | Vector: valoarea λ la fiecare pas k |
| `Cgacv` | Vector: criteriul GACV la fiecare pas k → folosit pentru selecție |
| `Csic` | Vector: criteriul SIC la fiecare pas k (verificare secundară) |
| `FRM_Condition` | Scalar: numărul de condiție SVD (max/min valori singulare) |
| `V` | Listă: setul de variabile active la fiecare pas |
| `Elbow` | Listă: setul de observații la "cot" la fiecare pas |

**Cum se folosește în pipeline:**
```r
est    <- FRM_Quantile_Regression(as.matrix(X), k, tau=0.05, I=25)
gacv   <- est$Cgacv
k_best <- which.min(gacv[is.finite(gacv)])   # pasul optim
lambda_optim <- abs(est$lambda[k_best])       # intensitatea riscului de coadă
beta_optim   <- est$beta[k_best, ]            # coeficienții de contagiune
```

### `qrL1Ini(x, y, a)` — inițializarea (intern, pas 0)
Calculează cuantila empirică a y și setul inițial de indici pentru iterație.
Nu se cheamă direct din afara algoritmului.

### `pf(beta0, y0, tau)` — check function (intern)
Funcția de pierdere asimetrică pentru regresia cantilă:
```
pf = τ·Σ(yᵢ - ŷᵢ)₊  +  (τ-1)·Σ(yᵢ - ŷᵢ)₋
```
Penalizează asimetric erorile: erorile pozitive (subestimare) cu τ=0.05,
erorile negative (supraestimare) cu 1-τ=0.95.

**Când este `source`-at în fișierul principal:**
```r
# Block 1 (Config) — ultima linie:
source("FRM_Statistics_Algorithm.R")
# Disponibil în tot restul fișierului principal după acest punct.
```

---

## FIȘIER 2: `FRM_SC_config.R`
### Rol: Configurație globală — primul care rulează

**Când este `source`-at:** La începutul fiecărui bloc major din `FRM_Stable_Banci.R`
și la începutul fiecărui fișier separat.

**Ce face la rulare:**
1. `rm(list = ls(all = TRUE))` — ȘTERGE TOT mediul R curent
2. `graphics.off()` — închide toate dispozitivele grafice
3. `setwd(wdir)` — schimbă directorul de lucru

**Variabile globale create:**

| Variabilă | Valoare | Folosită în |
|---|---|---|
| `channel` | `"Banci"` | toate numele de fișiere și foldere output |
| `wdir` | calea locală | `setwd()` |
| `date_start` / `date_end` | `20010103` / `20251231` | Block 4 (varying estimation) |
| `date_start_fixed` / `date_end_fixed` | idem | Block Fixed estimation |
| `date_start_source` / `date_end_source` | idem | `input_path` |
| `s` | `90` | fereastra rulantă în zile |
| `tau` | `0.05` | cuantila pentru regresia cantilă |
| `I` | `25` | iterații maxime LASSO |
| `J` | `7` | numărul de bănci |
| `stock_main` | `"BT"` | nodul principal în GIF și heatmap-uri |
| `output_path` | `Output/Banci/` | toate fișierele de ieșire |
| `website_path` | `Website/Banci/` | PNG-urile finale |
| `input_path` | `Input/Banci/...` | (neutilizat în pipeline actual) |

**Logica `output_path`:**
```r
output_path <- if (tau == 0.05 & s == 90) file.path("Output", channel) else
  file.path("Output", channel, paste0("Sensitivity/tau=", 100*tau, "/s=", s))
```
Dacă schimbi τ sau s pentru analiză de sensibilitate, output-urile merg automat
într-un subfolder separat fără să suprascrie rezultatele de bază.

**Foldere create automat:** 13 subfoldere sub `Output/Banci/` și `Website/Banci/`

**Funcții helper definite:**
```r
idx_start_at_or_after(key, target)  # primul index din 'key' >= target
idx_end_at_or_before(key, target)   # ultimul index din 'key' <= target
```

**Pachete încărcate:** readr, dplyr, tidyr, lubridate, zoo, ggplot2, data.table,
igraph, magick, scales, stringr, plotly, reshape2, quadprog, MASS

**ATENȚIE — Efectul `rm(list=ls())`:**
Când `source("FRM_SC_config.R")` rulează, șterge orice obiect din mediu,
inclusiv `theme_transparent_bottom` din utils. De aceea ordinea corectă este
**întotdeauna config înainte de utils**, nu invers.

---

## FIȘIER 3: `FRM_SC_utils.R`
### Rol: Helpers vizuali și utilitare — al doilea care rulează

**Când este `source`-at:** Imediat după config în fiecare bloc.

**Obiecte create:**

### `theme_transparent_bottom`
Temă ggplot2 cu fundal alb, fără grilă, legendă jos.
```r
# Folosit în fiecare grafic ggplot2:
ggplot(...) + theme_transparent_bottom
```

### `theme_transparent_min`
Variantă fără legendă, pentru grafice simple.

### `risk_colors`
```r
c("1. Low risk"="green", "2. General risk"="blue",
  "3. Elevated risk"="yellow", "4. High risk"="orange", "5. Severe risk"="red")
```
Folosit în `FRMColor_Banci.png` și `FRM_SC_history_outputs.R`.

### `BANK_NAME_MAP` / `pretty_bank(x)`
Mapează numele intern → nume de afișare:
```r
pretty_bank("BT")        # → "Banca Transilvania"
pretty_bank("BNPPariBas") # → "BNP Paribas"
```

### `pick_file_flexible(path, strong_regex, weak_regex, role)`
Găsește fișierul cu data cea mai recentă în nume. Încearcă regex strict,
apoi regex slab. Folosit la încărcarea lambdas_fixed CSV în GIF.

### `parse_date_robust(x)`
Parsează orice format de dată: Date, POSIXt, numeric (zile de la 1970),
"YYYY-MM-DD", "MM/DD/YYYY". Folosit în Block Centrality și HHI.

### `idx_start_at_or_after(key, target)` / `idx_end_at_or_before(key, target)`
Identice cu cele din config — redefinite și aici pentru când utils e soursat
independent fără config.

### `clamp01(v)`
```r
clamp01(v) = pmin(pmax(v, 0), 1)  # limitează la [0,1]
```
Folosit în HHI pentru a asigura că valorile sunt valide.

---

## FIȘIER 4: `FRM_SC_load_data.R`
### Rol: Citirea datelor brute și construirea Stage 20

**Când rulează în fișierul principal:** Block 3

**Dependențe necesare înainte de rulare:**
```r
source("FRM_SC_config.R")  # pentru output_path, channel
source("FRM_SC_utils.R")   # (opțional, nu folosește direct din utils)
```

**Ce citește:** `Baza_FRM_Banci.csv` din working directory

**Procesare:**
1. `ticker` ← coloana `Date` ca integer YYYYMMDD
2. `dates` ← `ticker` convertit la Date R
3. `banci_cols` ← hardcodat: `c("AlphaBank","BNPPariBas","BRD","BT","Erste","ING","Unicredit")`
4. `macro_cols` ← tot ce nu e `Date` sau bancă (auto-detectat prin `setdiff`)
5. `stock_return` ← matrice 5456×7 din coloanele banci_cols
6. `macro_return` ← matrice 5456×n_macro din coloanele macro_cols
7. `mktcap_index` ← matrice sintetică 5456×8:
   - Col 1 = ticker
   - Col 2-8 = ranking static 1,2,3,4,5,6,7 (același pentru fiecare zi)
   - Rațiune: în versiunea crypto era ranking real de market cap zilnic;
     pentru bănci nu avem date zilnice de market cap, deci toate băncile
     sunt mereu incluse cu rang egal

**Obiecte salvate în `stage_20_loaded.RData`:**
```r
save(ticker, stock_return, macro_return, mktcap_index, M_stock, M_macro, M, ...)
```

**Cum ajută fișierul principal:**
Când Block 3 din `FRM_Stable_Banci.R` rulează, creează exact aceste obiecte
în memorie și le salvează. Blocurile 4, 5, HHI, GMV, GIF toate fac
`load(stage_20_loaded.RData)` sau presupun că obiectele sunt deja în memorie.

---

## FIȘIER 5: `FRM_SC_estimation_varying.R`
### Rol: Bucla principală LASSO — Stage 30

**Când rulează în fișierul principal:** Block 4

**Dependențe:**
```r
load(file.path(output_path, "stage_20_loaded.RData"))
# Necesită: ticker, stock_return, macro_return, mktcap_index,
#           M_stock, M_macro, s, tau, I, J
```

**Logica completă a buclei:**

```
pentru fiecare zi t de la N0 la N1:
  dacă t < N0 + 90: skip (nu avem încă 90 de zile de istoric)
  
  1. biggest_index ← top-J bănci din mktcap_index[t, ] (mereu 1:7)
  2. rows_ret ← (t-90):(t-1) — fereastra de 90 de zile
  3. X ← cbind(stock_return[rows_ret, biggest_index],
               macro_return[rows_ret, ])
  4. X[!is.finite(X)] ← 0
  5. Elimină coloanele cu toate valorile zero (colSums(X!=0) > 0)
  6. pentru fiecare coloană k din X:
       est ← FRM_Quantile_Regression(X, k, tau=0.05, I=25)
       k_best ← which.min(est$Cgacv)   ← pasul cu GACV minim
       lambda[k] ← abs(est$lambda[k_best])
       A[k, -k] ← est$beta[k_best, ]   ← coeficienții de contagiune
  7. Salvează A ca adj_matrix_YYYYMMDD.csv
  8. Adaugă lambda-urile băncilor (nu ale macro) în FRM_individ[[t]]
```

**Obiecte produse:**
- `FRM_individ` — listă cu ~5.366 elemente, fiecare `[[date]] = matrix(1×J_t)`
  cu lambda-urile băncilor prezente în ziua respectivă
- `J_dynamic` — vector cu numărul de bănci valide per zi
- ~5.366 fișiere CSV în `Output/Banci/Adj_Matrices/`

**Salvate în `stage_30_varying.RData`:** `FRM_individ`, `J_dynamic`

**Timp de rulare:** 20-40 minute pe un singur core pentru 25 de ani.
Cu paralelizare `doParallel`: ~5-10 minute.

**Cum ajută fișierul principal:**
Block 4 din `FRM_Stable_Banci.R` conține exact această logică inline.
Fișierul separat există pentru rulare independentă sau înlocuire cu
versiunea paralelizată.

---

## FIȘIER 6: `FRM_SC_history_outputs.R`
### Rol: Construiește indicele FRM și scatter-ul de risc — Stage 50

**Când rulează în fișierul principal:** Block 5

**Dependențe:**
```r
load(file.path(output_path, "stage_20_loaded.RData"))
load(file.path(output_path, "stage_30_varying.RData"))
# Necesită: FRM_individ
```

**Ce face pas cu pas:**

**1. Gestionarea istoricului (RDS):**
```r
rds_path <- "Output/Banci/Lambda/FRM_Banci.rds"
FRM_history <- c(readRDS(rds_path), FRM_individ)  # combină cu istoricul
FRM_history <- FRM_history[!duplicated(names, fromLast=TRUE)]  # deduplicare
FRM_history <- FRM_history[order(as.Date(names))]              # sortare cronologică
saveRDS(FRM_history, rds_path)
```
Permite rulări incrementale: poți adăuga zile noi fără să reestimezi tot.

**2. `lambdas_wide.csv`:**
Matrice (date × bănci) cu toate lambda-urile individuale.
Rânduri = zile, coloane = bănci. Folosit în blocul Granger.

**3. Indicele FRM zilnic:**
```r
FRM_index$frm[t] = mean(lambda_AlphaBank[t], lambda_BNP[t], ..., lambda_Unicredit[t])
```
Media simplă a lambda-urilor băncilor disponibile în ziua t.

**4. Colorarea ECDF:**
```r
risk_ecdf <- ecdf(FRM_index$frm)  # funcția de distribuție empirică cumulativă
risk_pct   <- 100 * risk_ecdf(frm)  # percentila fiecărei valori
# 0-20%ile → Low risk (verde)
# 20-40%ile → General risk (albastru)
# 40-60%ile → Elevated risk (galben)
# 60-80%ile → High risk (portocaliu)
# 80-100%ile → Severe risk (roșu)
```
Colorarea e relativă la distribuția istorică, nu la valori absolute.

**5. Normalizarea [0,1]** (adăugată ulterior pentru FRMColor):
```r
frm_norm = (frm - min) / (max - min)
```
Folosită în `FRMColor_Banci.png` pentru scala Y lizibilă.

**Fișiere produse:**
- `Output/Banci/Lambda/FRM_Banci_index.csv`
- `Output/Banci/Lambda/lambdas_wide.csv`
- `Output/Banci/Lambda/FRM_Banci.rds`
- `Website/Banci/20251231/FRMColor_Banci.png`
- `Output/Banci/stage_50_index.RData`

---

## FIȘIER 7: `FRM_SC_frm_plot_stable.R`
### Rol: Graficul liniei albastre FRM

**Când rulează în fișierul principal:** Block 6

**Dependențe:**
```r
load(file.path(output_path, "stage_50_index.RData"))
# Necesită: FRM_index (data.frame cu coloanele 'date' și 'frm')
```

**Ce produce:**
PNG 1200×700px cu seria temporală FRM ca linie albastră continuă.
- `scale_x_date(date_breaks="1 year", date_labels="%Y")` — doar anii pe axa X
- `theme_transparent_bottom` + `angle=90` pe etichetele X

**Output:** `Website/Banci/20251231/FRM_Banci_Index.png`

**Cum ajută fișierul principal:**
Block 6 din fișierul principal conține exact aceeași logică.
Diferența: în fișierul principal e `FRM_Banci_Index.png` hardcodat;
în fișierul separat e `paste0("FRM_", channel, "_Index.png")`.

---

## FIȘIER 8: `FRM_SC_adj_heatmap.R`
### Rol: Funcțiile pentru heatmap-urile de adiacență

**Când rulează în fișierul principal:** Block 7

**Dependențe:**
```r
source("FRM_SC_config.R")
source("FRM_SC_utils.R")
library(ggplot2); library(reshape2)
```

**Funcții exportate și cum funcționează:**

### `plot_adj_matrix(adj_file, out_dir, digits, zero_eps, text_size)`
Funcția de bază. Citește un CSV de adiacență și desenează heatmap-ul.

Pași interni:
1. Citește matricea din CSV cu `read.csv(..., row.names=1)`
2. Calculează `L = max(abs(mat))` — scala pentru gradientul de culoare
3. Convertește la format lung (`reshape2::melt`)
4. Celulele cu `|valoare| < zero_eps` (implicit 0.001) → NA → afișate albe
5. Etichetele coloanelor sunt `toupper()` — majuscule
6. Gradientul: albastru (#3B82F6) = negativ, alb = zero, roșu (#EF4444) = pozitiv
7. `coord_fixed()` — celule pătrate
8. `ggsave()` la 300 DPI, 8×6 inch

### `plot_adj_for_date(date_input, fixed=FALSE, nearest_ok=TRUE, out_dir)`
Wrapper public. Convertește data la YYYYMMDD, găsește fișierul CSV,
cheamă `plot_adj_matrix`.

Cu `nearest_ok=TRUE` (implicit): dacă data exactă nu există,
folosește cea mai apropiată dată disponibilă și afișează un mesaj.

### `find_crash_day_frm()`
Citește `FRM_Banci_index.csv`, calculează `diff(frm)` și returnează
ziua cu cea mai mare scădere zilnică a FRM-ului.

### `find_crash_day_bank(bank="BT")`
Citește `stage_20_loaded.RData`, găsește ziua cu cel mai negativ
randament pentru banca specificată.

### `plot_adj_for_crash_day()` / `plot_adj_for_bank_crash()`
Combină funcțiile de găsire și de plotare.

**Cum ajută fișierul principal:**
Block 7 din fișierul principal conține toate aceste funcții inline,
plus apelurile directe pentru zilele de criză:
```r
plot_adj_for_date("2020-03-16", fixed=TRUE)
plot_adj_for_date("2009-03-04", fixed=TRUE)
# etc.
```

---

## BLOCURI SUPLIMENTARE DIN `FRM_Stable_Banci.R`
### (nu au fișiere separate corespunzătoare)

### Block Fixed Estimation (rulează de două ori în fișier!)
Construiește matricele de adiacență cu univers fix (aceleași 7 bănci mereu).
**Problemă identificată:** Există două versiuni ale acestui bloc în fișier —
una ca funcție `build_fixed_all()` și una inline. Rulează de fapt cea inline.
Diferența față de varying: `biggest_index_fixed` e calculat O SINGURĂ DATĂ
la `N0_fixed` și nu se mai schimbă.

### Block Centrality
Iterează prin toate adj CSVs din `Fixed/`, construiește grafuri `igraph` și
calculează per fiecare zi:
- `outdegree`, `indegree` (normalizate)
- `closeness` (cu ponderi = 1/|weight|)
- `betweenness` (cu ponderi = 1/|weight|)
- `eigenvector` (cu ponderi = |weight|)
- `out_strength`, `in_strength` (cu ponderi = |weight|)

Salvează `Centrality_ByNode_Fixed.csv`, `Centrality_Averages_Fixed.csv`,
`Centrality_FRM_Corr_Fixed.csv`, `Centrality_FRM_CorrP_Fixed.csv`
și 7 PNG-uri dual-axis (FRM vs fiecare centralitate).

### Block HHI
Calculează Herfindahl-Hirschman Index:
- **Lambda HHI:** concentrarea riscului între bănci
  `HHI(t) = Σ(λₖ/Σλ)²` → 1/7 = dispersat, 1 = concentrat în o singură bancă
- **MarketCap HHI:** (valori sintetice, de limitat utilitatea)

### Block GMV Portfolio
Optimizare Global Minimum Variance cu fereastră rulantă de 90 zile.
Produce grafice de wealth, volatilitate, Sharpe, boxplot-uri și
stacked area chart al ponderilor.

### Block Granger (Block 5.1)
Cauzalitate Granger FRM ↔ lambda individual per bancă.
La lag-uri: 1, 5, 10, 22 zile de tranzacționare.

### Block Network GIF
Animație a rețelei de contagiune pentru o perioadă aleasă (implicit 2020).
Limitat la o perioadă scurtă (max 500 zile) din motive de performanță.

---

## FLUXUL COMPLET DE DATE

```
Baza_FRM_Banci.csv
       │
       ▼ Block 3 (load_data)
stage_20_loaded.RData
  ├── ticker (int YYYYMMDD)
  ├── dates (Date)
  ├── stock_return (5456×7)
  ├── macro_return (5456×n)
  └── mktcap_index (5456×8)
       │
       ▼ Block 4 (estimation_varying) ──→ Adj_Matrices/adj_matrix_YYYYMMDD.csv (×5366)
stage_30_varying.RData
  ├── FRM_individ [[date]] → matrix(1×J_t)
  └── J_dynamic
       │
       ▼ Block 5 (history_outputs)
stage_50_index.RData          lambdas_wide.csv
  └── FRM_index (date, frm)   FRM_Banci_index.csv
       │                      FRMColor_Banci.png
       ▼ Block 6
FRM_Banci_Index.png
       │
       ▼ Block Fixed ──────→ Adj_Matrices/Fixed/adj_matrix_YYYYMMDD.csv (×5366)
lambdas_fixed_*.csv
       │
       ▼ Block Centrality
Centrality_ByNode_Fixed.csv
Centrality_Averages_Fixed.csv
Centrality_FRM_Corr_Fixed.csv    FRM_vs_*_Fixed_DualAxis.png (×7)
       │
       ▼ Block HHI
HHI_lambda.csv
HHI_Lambda.png
       │
       ▼ Block Granger
Granger_FRM_vs_Banks.csv
Granger_Heatmap_lag5.png
Granger_Heatmap_lag22.png
       │
       ▼ Block GMV
GMV_Weights_Daily_wide.csv
CAPM_Portfolio_vs_Individual_Banks_Wealth.png
CAPM_Boxplot_RollVol_90d_ByBank.png
GMV_Weights_AllBanks_StackedArea.png
       │
       ▼ Block 7 (adj_heatmap)
AdjMatrix_20020706.png
AdjMatrix_20090304.png
AdjMatrix_20200316.png  (etc.)
       │
       ▼ Block GIF
Network_20010103_20251231_Banci.gif
```

---

## PROBLEME CUNOSCUTE ÎN FIȘIERUL PRINCIPAL

| Linie | Problemă | Status |
|---|---|---|
| ~200 | `build_fixed_all()` definită și chemată, dar există și o versiune inline mai jos | Rulează cea inline |
| ~50 | `date_banci <- read.csv(...)` — load redundant în Block 1 | Inofensiv, de curățat |
| ~800 | `!all(colSums(X != 0) > 0)` în loc de `keep <- colSums...` | Corectat în versiunile recente |
| ~2200 | `GMV_Weights_AllCoins_StackedArea.png` — "Coins" în loc de "Banks" | Corectat |
| Multiplu | `source("FRM_SC_utils.R"); source("FRM_SC_config.R")` — ordinea greșită | De verificat per bloc |
| GIF | Bucla pe toți cei 5456 ani durează ore | Limitat la 2020 cu `gif_start`/`gif_end` |

