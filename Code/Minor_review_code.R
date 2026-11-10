###############################################################################
# SWI anomaly: probit M0-M4, robustness, multicollinearity, marginal effects
# Script consolide (ordre d'execution: de haut en bas)
###############################################################################
library(readxl); library(lubridate); library(dplyr)
library(fixest); library(pROC); library(ggplot2); library(patchwork)
library(marginaleffects)

###############################################################################
# 1. Import et fusion
###############################################################################
df <- read_excel(choose.files())      # base habituelle
ev <- read_excel(choose.files())      # CatNat_drought_events_SWI_19902025anom.xlsx

pad5 <- function(x) { x <- as.character(x); paste0(strrep("0", pmax(0, 5 - nchar(x))), x) }

df$code_insee  <- pad5(df$code_insee)
df$d_start     <- as.Date(df$event_start)
df$d_end       <- as.Date(df$event_end)
df$year        <- year(df$d_start)
df$departement <- substr(df$code_insee, 1, 2)

ev2 <- ev %>%
  mutate(code_insee = pad5(code_insee),
         d_start    = as.Date(`Début d'Evénement`),
         d_end      = as.Date(`Fin d'Evénement`)) %>%
  select(code_insee, d_start, d_end, swi_min, swi_anom_min) %>%
  distinct(code_insee, d_start, d_end, .keep_all = TRUE)

d_all <- left_join(df, ev2, by = c("code_insee", "d_start", "d_end"))
cat("Sans anomalie :", sum(is.na(d_all$swi_anom_min)), "\n")

d <- subset(d_all, year >= 2002)
cat("Échantillon 2002+ :", nrow(d), "(attendu 64389)\n")

d$dur     <- as.numeric(difftime(d$d_end, d$d_start, units = "days")) / 30.44
d$start_m <- factor(as.integer(format(d$d_start, "%m")))

# Echantillon commun : departements sans variation de la decision exclus (2B)
bad_dept <- d %>% group_by(departement) %>%
  summarise(n_outcomes = n_distinct(decision), .groups = "drop") %>%
  filter(n_outcomes < 2) %>% pull(departement)
d_common <- d %>% filter(!departement %in% bad_dept)
cat("Départements exclus :", bad_dept, "| N commun :", nrow(d_common), "(attendu 64386)\n")

# Niveaux de CSS pour les effets conditionnels
css_pts <- c(Q1 = unname(quantile(d$CSS, 0.25)),
             Mean = mean(d$CSS),
             Median = median(d$CSS))
print(round(css_pts, 2))   # 70.13 / 81.36 / 95.35

ord <- c("^swi_anom_min$", "^CSS$", ":", "SPI", "SPEI", "CDD", "dur")

###############################################################################
# 2. Probit M0 a M4 (anomalie SWI) : + SPI, puis + SPEI, puis + CDD
###############################################################################
h  <- "+ old_housing_share + share_individual_houses + log(n_dwellings)"
fe <- " | year + departement"

a0 <- feglm(as.formula(paste("decision ~ swi_anom_min + CSS + swi_anom_min:CSS", h)),
            family = binomial("probit"), cluster = ~code_insee, data = d)
mk <- function(rhs) feglm(as.formula(paste("decision ~", rhs, h, fe)),
                          family = binomial("probit"), cluster = ~code_insee, data = d)
a1 <- mk("swi_anom_min + CSS + swi_anom_min:CSS")
a2 <- mk("swi_anom_min + CSS + swi_anom_min:CSS + SPI")
a3 <- mk("swi_anom_min + CSS + swi_anom_min:CSS + SPI + SPEI")
a4 <- mk("swi_anom_min + CSS + swi_anom_min:CSS + SPI + SPEI + CDD")

etable(a0, a1, a2, a3, a4,
       headers = c("M0 (no FE)", "M1", "M2 (+SPI)", "M3 (+SPEI)", "M4 (+CDD)"),
       fitstat = ~ n + pr2 + bic, keep = "swi|CSS|SPI|SPEI|CDD", order = ord)

# Controles de bati dans les cinq modeles (Panel C du tableau 5.1)
etable(a0, a1, a2, a3, a4,
       headers = c("M0", "M1", "M2", "M3", "M4"),
       keep = "old_housing|share_individual|n_dwellings", fitstat = ~ n)

auc_m <- function(m) as.numeric(auc(roc(m$y, as.numeric(predict(m, type = "response")), quiet = TRUE)))
round(c(M1 = auc_m(a1), M2 = auc_m(a2), M3 = auc_m(a3), M4 = auc_m(a4)), 3)

###############################################################################
# 3. Avec duree et mois de debut
###############################################################################
hb <- "+ log(dur + 1) + start_m + old_housing_share + share_individual_houses + log(n_dwellings) | year + departement"
mkb <- function(rhs) feglm(as.formula(paste("decision ~", rhs, hb)),
                           family = binomial("probit"), cluster = ~code_insee, data = d)
b1 <- mkb("swi_anom_min + CSS + swi_anom_min:CSS")
b2 <- mkb("swi_anom_min + CSS + swi_anom_min:CSS + SPI")
b3 <- mkb("swi_anom_min + CSS + swi_anom_min:CSS + SPI + SPEI")
b4 <- mkb("swi_anom_min + CSS + swi_anom_min:CSS + SPI + SPEI + CDD")

etable(b1, b2, b3, b4, headers = c("M1", "M2 (+SPI)", "M3 (+SPEI)", "M4 (+CDD)"),
       fitstat = ~ n + pr2 + bic, keep = "swi|CSS|SPI|SPEI|CDD|dur", order = ord)

# Poids de la duree et du mois de debut (M4)
wald(b4, keep = "start_m")

h0 <- "+ old_housing_share + share_individual_houses + log(n_dwellings) | year + departement"
mkh <- function(rhs) feglm(as.formula(paste("decision ~", rhs, h0)),
                           family = binomial("probit"), cluster = ~code_insee, data = d)
base  <- "swi_anom_min + CSS + swi_anom_min:CSS + SPI + SPEI + CDD"
c_dur  <- mkh(paste(base, "+ log(dur + 1)"))
c_mois <- mkh(paste(base, "+ start_m"))
c_les2 <- mkh(paste(base, "+ log(dur + 1) + start_m"))

etable(a4, c_dur, c_mois, c_les2,
       headers = c("M4", "+ durée", "+ mois", "+ les deux"),
       fitstat = ~ n + pr2 + bic, keep = "swi|CSS|SPI|SPEI|CDD|dur", order = ord)

###############################################################################
# 4. Effets conditionnels sur l'indice latent (tab:me)
###############################################################################
me_c <- function(m, css = mean(d$CSS)) {
  b <- coef(m); V <- vcov(m)
  est <- b["swi_anom_min"] + b["swi_anom_min:CSS"] * css
  se  <- sqrt(V["swi_anom_min", "swi_anom_min"] +
                css^2 * V["swi_anom_min:CSS", "swi_anom_min:CSS"] +
                2 * css * V["swi_anom_min", "swi_anom_min:CSS"])
  c(estimate = est, se = se, z = est / se)
}
me_tab <- function(models) {
  do.call(rbind, lapply(names(css_pts), function(nm) {
    r <- t(sapply(models, function(m) me_c(m, css = css_pts[[nm]])))
    colnames(r) <- c("estimate", "se", "z")
    data.frame(CSS = nm, Model = rownames(r), round(r, 3), row.names = NULL)
  }))
}
me_tab(list(M1 = a1, M2 = a2, M3 = a3, M4 = a4))   # Panel A : sans duree
me_tab(list(M1 = b1, M2 = b2, M3 = b3, M4 = b4))   # Panel B : avec duree et mois

###############################################################################
# 5. Effets marginaux moyens en probabilite (tableau 5.1, Panel D)
#    Modeles avec indicatrices (marginaleffects ne gere pas les FE absorbes)
###############################################################################
mk_dum <- function(rhs) feglm(
  as.formula(paste("decision ~", rhs,
                   "+ old_housing_share + share_individual_houses + log(n_dwellings)",
                   "+ factor(year) + factor(departement)")),
  family = binomial("probit"), cluster = ~code_insee, data = d_common)

f0 <- feglm(decision ~ swi_anom_min * CSS + old_housing_share +
              share_individual_houses + log(n_dwellings),
            family = binomial("probit"), cluster = ~code_insee, data = d_common)
f1 <- mk_dum("swi_anom_min * CSS")
f2 <- mk_dum("swi_anom_min * CSS + SPI")
f3 <- mk_dum("swi_anom_min * CSS + SPI + SPEI")
f4 <- mk_dum("swi_anom_min * CSS + SPI + SPEI + CDD")

# Controle : coefficients SWI (doivent redonner -1.270 pour M1 et -1.236 pour M4)
c(coef(f1)["swi_anom_min"], coef(f4)["swi_anom_min"])

# AME globaux (attendu : M1 -0.400 (0.0042) ; M4 -0.374 (0.0049))
for (nm in c("f0", "f1", "f2", "f3", "f4")) {
  m <- get(nm)
  cat("\n==", nm, "==\n")
  print(avg_slopes(m, variables = "swi_anom_min", vcov = vcov(m)))
}

# AME aux trois niveaux de CSS (M1 et M4)
ame_at <- function(m, css) {
  nd <- transform(d_common, CSS = css)
  avg_slopes(m, variables = "swi_anom_min", newdata = nd, vcov = vcov(m))
}
lapply(css_pts, function(x) ame_at(f1, x))
lapply(css_pts, function(x) ame_at(f4, x))

###############################################################################
# 6. Figure : ROC + calibration
###############################################################################
roc1 <- roc(a1$y, as.numeric(predict(a1, type = "response")), quiet = TRUE)
roc2 <- roc(a2$y, as.numeric(predict(a2, type = "response")), quiet = TRUE)
roc3 <- roc(a3$y, as.numeric(predict(a3, type = "response")), quiet = TRUE)
roc4 <- roc(a4$y, as.numeric(predict(a4, type = "response")), quiet = TRUE)

mk_roc <- function(r, nm) data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Model = nm)
roc_df <- bind_rows(mk_roc(roc1, "Baseline"), mk_roc(roc2, "+ SPI"),
                    mk_roc(roc3, "+ SPEI"),   mk_roc(roc4, "+ CDD"))
roc_df$Model <- factor(roc_df$Model, levels = c("Baseline", "+ SPI", "+ SPEI", "+ CDD"))

legend_labels <- c(paste0("Baseline (AUC = ", round(as.numeric(roc1$auc), 3), ")"),
                   paste0("+ SPI (AUC = ",    round(as.numeric(roc2$auc), 3), ")"),
                   paste0("+ SPEI (AUC = ",   round(as.numeric(roc3$auc), 3), ")"),
                   paste0("+ CDD (AUC = ",    round(as.numeric(roc4$auc), 3), ")"))

p1 <- ggplot(roc_df, aes(FPR, TPR, colour = Model)) +
  geom_line(linewidth = 1.15) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey70") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(labels = legend_labels,
                      values = c("#4E79A7", "#F28E2B", "#59A14F", "#B22222")) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "False Positive Rate", y = "True Positive Rate", colour = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = c(.67, .22),
        legend.background = element_rect(fill = "white", colour = "grey75"),
        legend.key = element_blank(), legend.text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = .6))

calib <- data.frame(Predicted = as.numeric(predict(a4, type = "response")), Observed = a4$y) |>
  mutate(Decile = ntile(Predicted, 10)) |>
  group_by(Decile) |>
  summarise(Predicted = mean(Predicted), Observed = mean(Observed), .groups = "drop")

p2 <- ggplot(calib, aes(Predicted, Observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey70") +
  geom_line(linewidth = 1, colour = "#B22222") +
  geom_point(size = 3, colour = "#1B4F72") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Predicted Probability", y = "Observed Probability") +
  theme_bw(base_size = 13) +
  theme(panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = .6))

figure_validation_anom <- p1 + p2 + plot_layout(ncol = 2)
figure_validation_anom
ggsave("Figure_Model_Validation_anomaly.pdf", figure_validation_anom,
       width = 11, height = 5, dpi = 600)

###############################################################################
# 7. Contribution statistique propre des indicateurs (tests de Wald dans M4)
###############################################################################
wald(a4, keep = "^SPI$")
wald(a4, keep = "^SPEI$")
wald(a4, keep = "^CDD$")
wald(a4, keep = "^(SPI|SPEI|CDD)$")      # conjoint

wald(b4, keep = "^SPI$")                 # avec duree et mois de debut
wald(b4, keep = "^SPEI$")
wald(b4, keep = "^CDD$")
wald(b4, keep = "^(SPI|SPEI|CDD)$")

###############################################################################
# 8. Correlations et multicolinearite
###############################################################################
vars   <- c("swi_anom_min", "SPI", "SPEI", "CDD", "CSS")
cor_mat <- cor(d[, vars], use = "complete.obs", method = "pearson")
cat("\n=== CORRELATIONS DE PEARSON ===\n"); print(round(cor_mat, 3))

# VIF : variables centrees, effets fixes annee et departement neutralises
v <- d %>%
  select(swi_anom_min, CSS, SPI, SPEI, CDD, old_housing_share,
         share_individual_houses, n_dwellings, year, departement) %>%
  filter(complete.cases(.), n_dwellings > 0) %>%
  mutate(swi_c = swi_anom_min - mean(swi_anom_min),
         css_c = CSS - mean(CSS),
         interaction = swi_c * css_c,
         log_dwellings = log(n_dwellings))

X <- v %>% select(swi_c, css_c, interaction, SPI, SPEI, CDD,
                  old_housing_share, share_individual_houses, log_dwellings)

X_res <- as.data.frame(sapply(names(X), function(nm)
  resid(feols(as.formula(paste(nm, "~ 1 | year + departement")), data = v))))

vif_fun <- function(variable, data) {
  r2 <- summary(lm(reformulate(setdiff(names(data), variable), response = variable),
                   data = data))$r.squared
  1 / (1 - r2)
}
vif_results <- data.frame(Variable = names(X_res),
                          VIF = sapply(names(X_res), vif_fun, data = X_res))
cat("\n=== VIF APRES NEUTRALISATION DES FE ===\n")
print(vif_results %>% mutate(VIF = round(VIF, 3)), row.names = FALSE)

# Correlations brute et partielle SPEI / decision
cat("\n=== CORRELATION BRUTE SPEI / DECISION ===\n")
print(cor(d$SPEI, d$decision, use = "complete.obs"))

controls <- paste("swi_anom_min * CSS + SPI + CDD +",
                  "old_housing_share + share_individual_houses +",
                  "log(n_dwellings) | year + departement")
r_spei <- resid(feols(as.formula(paste("SPEI ~", controls)), data = d))
r_dec  <- resid(feols(as.formula(paste("decision ~", controls)), data = d))
cat("\n=== CORRELATION PARTIELLE SPEI / DECISION ===\n")
print(cor(r_spei, r_dec))

# Figure : matrice de correlation (sans diagonale)
labels <- c("SWI", "SPI", "SPEI", "CDD", "CSS")
cm <- cor_mat; rownames(cm) <- colnames(cm) <- labels

cor_df <- as.data.frame(as.table(cm)) |>
  setNames(c("Var1", "Var2", "r")) |>
  mutate(i = match(Var1, labels), j = match(Var2, labels)) |>
  filter(i > j) |>
  mutate(Var1 = factor(Var1, levels = rev(labels[-1])),
         Var2 = factor(Var2, levels = labels[-length(labels)]),
         txt  = sprintf("%.2f", r))

fig_cor <- ggplot(cor_df, aes(Var2, Var1, fill = r)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = txt), size = 4, colour = "grey15") +
  scale_fill_gradient2(low = "#4E79A7", mid = "#F4F4F4", high = "#B22222",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson's r") +
  coord_equal() + labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = .6),
        legend.position = "right", legend.key.height = unit(1.4, "cm"))
fig_cor
ggsave("Figure_Correlation_matrix.pdf", fig_cor, width = 6.5, height = 5, dpi = 600)

# VIF sans effets fixes (colonne 1)
X_df <- as.data.frame(X)
vif_no_fe <- data.frame(Variable = names(X_df),
                        VIF = sapply(names(X_df), vif_fun, data = X_df))
print(vif_no_fe %>% mutate(VIF = round(VIF, 2)), row.names = FALSE)

###############################################################################
# 9. Robustesse : heterogeneite territoriale (effets fixes annee / departement)
###############################################################################
rhs <- paste("swi_anom_min + CSS + swi_anom_min:CSS +",
             "old_housing_share + share_individual_houses + log(n_dwellings)")

fit_fe <- function(fe = NULL) {
  f <- if (is.null(fe)) as.formula(paste("decision ~", rhs))
  else as.formula(paste("decision ~", rhs, "|", fe))
  feglm(f, family = binomial("probit"), cluster = ~code_insee, data = d_common)
}
r0     <- fit_fe()
r_year <- fit_fe("year")
r_dept <- fit_fe("departement")
r_both <- fit_fe("year + departement")

c(No_FE = nobs(r0), Year_FE = nobs(r_year), Department_FE = nobs(r_dept), Both_FE = nobs(r_both))

etable(r0, r_year, r_dept, r_both,
       headers = c("No FE", "Year FE", "Department FE", "Year + Department FE"),
       fitstat = ~ n + pr2 + bic)

# Wald conjoint des effets departementaux (annee en effet fixe)
r_dept_test <- feglm(as.formula(paste("decision ~", rhs, "+ i(departement) | year")),
                     family = binomial("probit"), cluster = ~code_insee, data = d_common)
wald_dept <- wald(r_dept_test, keep = "departement::")
print(wald_dept)

# Meme test avec le modele complet M4
r_full <- feglm(as.formula(paste("decision ~", rhs, "+ SPI + SPEI + CDD + i(departement) | year")),
                family = binomial("probit"), cluster = ~code_insee, data = d_common)
wald_full <- wald(r_full, keep = "departement::")
print(wald_full)

# Test du rapport de vraisemblance : departements
lr  <- 2 * (as.numeric(logLik(r_both)) - as.numeric(logLik(r_year)))
ddl <- r_both$fixef_sizes["departement"] - 1
cat("LR =", round(lr, 1), " ddl =", ddl, " p =", pchisq(lr, ddl, lower.tail = FALSE), "\n")

###############################################################################
# 10. Duree de la periode declaree et decision (5.4)
###############################################################################
# Taux de reconnaissance par classe de duree (en mois)
d %>%
  mutate(dur_class = cut(dur, breaks = c(0, 3, 6, 9, 12, Inf), right = TRUE)) %>%
  group_by(dur_class) %>%
  summarise(n = n(), rate = round(mean(decision), 3), .groups = "drop")

# Duree moyenne selon la decision (1 = acceptee)
d %>% group_by(decision) %>%
  summarise(mean_dur = round(mean(dur), 2), median_dur = round(median(dur), 2), .groups = "drop")

# Part des periodes qui se terminent le dernier jour d'un trimestre
d %>%
  mutate(q_end = month(d_end) %in% c(3, 6, 9, 12) &
           d_end == ceiling_date(d_end, "month") - days(1)) %>%
  group_by(decision) %>%
  summarise(share_quarter_end = round(mean(q_end), 3), .groups = "drop")

########################### Accesories 
# Figure : matrice de correlation (sans diagonale)
labels <- c("SWI", "SPI", "SPEI", "CDD", "CSS")
cm <- cor_mat; rownames(cm) <- colnames(cm) <- labels

cor_df <- as.data.frame(as.table(cm)) |>
  setNames(c("Var1", "Var2", "r")) |>
  mutate(i = match(Var1, labels), j = match(Var2, labels)) |>
  filter(i > j) |>
  mutate(Var1 = factor(Var1, levels = rev(labels[-1])),
         Var2 = factor(Var2, levels = labels[-length(labels)]),
         txt  = sprintf("%.2f", r))

fig_cor <- ggplot(cor_df, aes(Var2, Var1, fill = r)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = txt), size = 5.5, colour = "grey15") +
  scale_fill_gradient2(low = "#4E79A7", mid = "#F4F4F4", high = "#B22222",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson's r") +
  coord_equal() + labs(x = NULL, y = NULL) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = .6),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "right",
        legend.key.height = unit(1.3, "cm"))
fig_cor
ggsave("Figure_Correlation_matrix.pdf", fig_cor, width = 6.4, height = 4.4)
