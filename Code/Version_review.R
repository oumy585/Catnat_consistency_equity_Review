
###############################################################################
# 0. Packages
###############################################################################

library(readxl)
library(lubridate)
library(fixest)
library(marginaleffects)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(pROC)

###############################################################################
# 1. Import des données
###############################################################################

# C:\Users\OumyDIONE\Downloads\Drought_risk_review\DATA\Data_merging"
df <- read_excel(choose.files())

names(df)
summary(df$decision)


###############################################################################
# 2. Construction des effets fixes
###############################################################################

# Effet fixe année
df$year <- year(as.Date(df$`event_start`))
table(df$year)

# Effet fixe département
df$departement <- substr(df$code_insee, 1, 2)
table(df$departement)[1:10]



# Sample
###############################################################

df2002 <- subset(df, year >= 2002)

###############################################################
# M0 : Hydro-geotechnical baseline (No fixed effects)
###############################################################

m0 <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings),
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df2002
)

###############################################################
# M1 : Hydro-geotechnical baseline
# + Year & Department fixed effects
###############################################################

m1 <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    year + departement,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df2002
)

###############################################################
# M2 : + SPEI
###############################################################

m2 <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    SPEI +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    year + departement,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df2002
)

###############################################################
# M3 : + SPI
###############################################################

m3 <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    SPEI +
    SPI +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    year + departement,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df2002
)

###############################################################
# M4 : Full model (+ CDD)
###############################################################

m4 <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    SPEI +
    SPI +
    CDD +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    year + departement,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df2002
)

###############################################################
# Regression table
###############################################################

etable(
  m0,
  m1,
  m2,
  m3,
  m4,
  headers = c(
    "Baseline\n(No FE)",
    "+ Year & Department FE",
    "+ SPEI",
    "+ SPI",
    "+ CDD"
  ),
  fitstat = ~ n + pr2 + bic
)


############################################################
## Figure : ROC curves + Calibration
############################################################
############################################################
## Predictive performance
############################################################

library(pROC)
library(dplyr)

model_perf <- function(model){
  
  p <- predict(model, type = "response")
  y <- model$y
  
  auc <- as.numeric(
    roc(
      y,
      p,
      quiet = TRUE
    )$auc
  )
  
  logloss <- -mean(
    y * log(pmax(p,1e-15)) +
      (1-y) * log(pmax(1-p,1e-15))
  )
  
  brier <- mean(
    (p-y)^2
  )
  
  tibble(
    AUC = auc,
    LogLoss = logloss,
    Brier = brier
  )
  
}

############################################################
# Performance table
############################################################

perf <- bind_rows(
  
  model_perf(m1),
  model_perf(m2),
  model_perf(m3),
  model_perf(m4)
  
)

perf$Model <- c(
  
  "Baseline",
  "+ SPEI",
  "+ SPI",
  "+ CDD"
  
)

perf

############################################################
## Figure : ROC curves + Calibration
############################################################

library(ggplot2)
library(patchwork)

############################################################
# ROC curves
############################################################

roc1 <- roc(m1$y, predict(m1, type = "response"), quiet = TRUE)
roc2 <- roc(m2$y, predict(m2, type = "response"), quiet = TRUE)
roc3 <- roc(m3$y, predict(m3, type = "response"), quiet = TRUE)
roc4 <- roc(m4$y, predict(m4, type = "response"), quiet = TRUE)

roc_df <- bind_rows(
  
  data.frame(
    FPR = 1 - roc1$specificities,
    TPR = roc1$sensitivities,
    Model = "Baseline"
  ),
  
  data.frame(
    FPR = 1 - roc2$specificities,
    TPR = roc2$sensitivities,
    Model = "+ SPEI"
  ),
  
  data.frame(
    FPR = 1 - roc3$specificities,
    TPR = roc3$sensitivities,
    Model = "+ SPI"
  ),
  
  data.frame(
    FPR = 1 - roc4$specificities,
    TPR = roc4$sensitivities,
    Model = "+ CDD"
  )
  
)

############################################################
# Legend labels
############################################################

legend_labels <- c(
  
  paste0("Baseline (AUC = ", round(as.numeric(roc1$auc),3),")"),
  paste0("+ SPEI (AUC = ", round(as.numeric(roc2$auc),3),")"),
  paste0("+ SPI (AUC = ", round(as.numeric(roc3$auc),3),")"),
  paste0("+ CDD (AUC = ", round(as.numeric(roc4$auc),3),")")
  
)

############################################################
# ROC panel
############################################################

p1 <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = TPR,
    colour = Model
  )
) +
  
  geom_line(linewidth = 1.15) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    colour = "grey70"
  ) +
  
  coord_equal(
    xlim = c(0,1),
    ylim = c(0,1)
  ) +
  
  scale_colour_manual(
    
    labels = legend_labels,
    
    values = c(
      "#4E79A7",
      "#F28E2B",
      "#59A14F",
      "#B22222"
    )
    
  ) +
  
  scale_x_continuous(
    breaks = seq(0,1,0.25)
  ) +
  
  scale_y_continuous(
    breaks = seq(0,1,0.25)
  ) +
  
  labs(
    
    x = "False Positive Rate",
    
    y = "True Positive Rate",
    
    colour = NULL
    
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    
    legend.position = c(.67,.22),
    
    legend.background = element_rect(
      fill = "white",
      colour = "grey75"
    ),
    
    legend.key = element_blank(),
    
    legend.text = element_text(size = 10),
    
    panel.grid.major = element_line(
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      fill = NA,
      linewidth = .6
    )
    
  )

############################################################
# Calibration
############################################################

pred <- predict(
  m4,
  type = "response"
)

calib <- data.frame(
  
  Predicted = pred,
  
  Observed = m4$y
  
) |>
  
  mutate(
    
    Decile = ntile(Predicted,10)
    
  ) |>
  
  group_by(Decile) |>
  
  summarise(
    
    Predicted = mean(Predicted),
    
    Observed = mean(Observed),
    
    .groups = "drop"
    
  )

############################################################
# Calibration panel
############################################################

p2 <- ggplot(
  
  calib,
  
  aes(
    Predicted,
    Observed
  )
  
) +
  
  geom_abline(
    
    intercept = 0,
    
    slope = 1,
    
    linetype = 2,
    
    colour = "grey70"
    
  ) +
  
  geom_line(
    
    linewidth = 1,
    
    colour = "#B22222"
    
  ) +
  
  geom_point(
    
    size = 3,
    
    colour = "#1B4F72"
    
  ) +
  
  coord_equal(
    
    xlim = c(0,1),
    
    ylim = c(0,1)
    
  ) +
  
  scale_x_continuous(
    breaks = seq(0,1,0.25)
  ) +
  
  scale_y_continuous(
    breaks = seq(0,1,0.25)
  ) +
  
  labs(
    
    x = "Predicted Probability",
    
    y = "Observed Probability"
    
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    
    panel.grid.major = element_line(
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      fill = NA,
      linewidth = .6
    )
    
  )

############################################################
# Final figure
############################################################

figure_validation <-
  
  p1 + p2 +
  
  plot_layout(
    ncol = 2
  )

figure_validation

ggsave(
  
  "Figure_Model_Validation.pdf",
  
  figure_validation,
  
  width = 11,
  
  height = 5,
  
  dpi = 600
  
)

########################################################################################################################
#####################################################################################################################

###############################################################################
# ROBUSTNESS: TERRITORIAL HETEROGENEITY
# Contribution of year and department fixed effects
###############################################################################

library(dplyr)
library(fixest)

###############################################################################
# 1. IDENTIFY DEPARTMENTS WITH NO VARIATION IN THE OUTCOME
###############################################################################

dept_check <- df2002 %>%
  group_by(departement) %>%
  summarise(
    N = n(),
    Mean_decision = mean(decision, na.rm = TRUE),
    N_outcomes = n_distinct(decision),
    .groups = "drop"
  )

bad_dept <- dept_check %>%
  filter(N_outcomes < 2)

bad_dept


###############################################################################
# 2. CONSTRUCT A COMMON SAMPLE
###############################################################################

bad_dept_id <- bad_dept$departement

df_common <- df2002 %>%
  filter(!departement %in% bad_dept_id)

# Checks
nrow(df2002)
nrow(df_common)

table(df_common$departement == "2B")


###############################################################################
# 3. R0: NO FIXED EFFECTS
###############################################################################

r0_common <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings),
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df_common
)

###############################################################################
# 4. R1: YEAR FIXED EFFECTS ONLY
###############################################################################

r_year_common <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    year,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df_common
)

###############################################################################
# 5. R2: DEPARTMENT FIXED EFFECTS ONLY
###############################################################################

r_dept_common <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    departement,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df_common
)

###############################################################################
# 6. R3: YEAR + DEPARTMENT FIXED EFFECTS
###############################################################################

r_both_common <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    year + departement,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df_common
)

###############################################################################
# 7. VERIFY THAT ALL MODELS USE EXACTLY THE SAME SAMPLE
###############################################################################

c(
  No_FE = nobs(r0_common),
  Year_FE = nobs(r_year_common),
  Department_FE = nobs(r_dept_common),
  Year_Department_FE = nobs(r_both_common)
)



###############################################################################
# 8. MODEL COMPARISON
###############################################################################

etable(
  r0_common,
  r_year_common,
  r_dept_common,
  r_both_common,
  headers = c(
    "No FE",
    "Year FE",
    "Department FE",
    "Year + Department FE"
  ),
  fitstat = ~ n + pr2 + bic
)


###############################################################################
# 9. JOINT TEST OF DEPARTMENT EFFECTS
###############################################################################

glm_year <- glm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) +
    factor(year),
  family = binomial(link = "probit"),
  data = df_common
)


glm_year_dept <- glm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) +
    factor(year) +
    factor(departement),
  family = binomial(link = "probit"),
  data = df_common
)


###############################################################################
# 10. LIKELIHOOD-RATIO TEST:
# ARE DEPARTMENT EFFECTS JOINTLY SIGNIFICANT?
###############################################################################

dept_LR_test <- anova(
  glm_year,
  glm_year_dept,
  test = "Chisq"
)

dept_LR_test

print(dept_LR_test)


###############################################################################
# JOINT SIGNIFICANCE OF DEPARTMENT FIXED EFFECTS
###############################################################################

# Log-likelihoods
LL_year <- as.numeric(logLik(glm_year))
LL_year_dept <- as.numeric(logLik(glm_year_dept))

# Number of estimated parameters
K_year <- length(coef(glm_year))
K_year_dept <- length(coef(glm_year_dept))

# Likelihood-ratio statistic
LR <- 2 * (LL_year_dept - LL_year)

# Degrees of freedom = number of additional department parameters
DF <- K_year_dept - K_year

# p-value
P_value <- pchisq(
  LR,
  df = DF,
  lower.tail = FALSE
)

###############################################################################
# RESULTS
###############################################################################

data.frame(
  Test = "Joint significance of department fixed effects",
  LogLik_Year_FE = LL_year,
  LogLik_Year_Department_FE = LL_year_dept,
  LR_statistic = LR,
  df = DF,
  p_value = P_value
)

length(coef(glm_year))
length(coef(glm_year_dept))

logLik(glm_year)
logLik(glm_year_dept)


####################################
fitstat(r0_common, "pr2")
fitstat(r_year_common, "pr2")
fitstat(r_dept_common, "pr2")
fitstat(r_both_common, "pr2")

BIC(r0_common)
BIC(r_year_common)
BIC(r_dept_common)
BIC(r_both_common)


############################################################################# avec M4
glm_year_full <- glm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    SPEI +
    SPI +
    CDD +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) +
    factor(year),
  family = binomial(link = "probit"),
  data = df_common
)

glm_year_dept_full <- glm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    SPEI +
    SPI +
    CDD +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) +
    factor(year) +
    factor(departement),
  family = binomial(link = "probit"),
  data = df_common
)

LL_year_full <- as.numeric(logLik(glm_year_full))
LL_year_dept_full <- as.numeric(logLik(glm_year_dept_full))

K_year_full <- length(coef(glm_year_full))
K_year_dept_full <- length(coef(glm_year_dept_full))

LR_full <- 2 * (LL_year_dept_full - LL_year_full)
DF_full <- K_year_dept_full - K_year_full

P_value_full <- pchisq(
  LR_full,
  df = DF_full,
  lower.tail = FALSE
)

data.frame(
  Test = "Joint significance of department FE - Full model",
  LogLik_Year_FE = LL_year_full,
  LogLik_Year_Department_FE = LL_year_dept_full,
  LR_statistic = LR_full,
  df = DF_full,
  p_value = P_value_full
)


#######################################################################################
############################################################
# ROBUSTNESS CHECK 4
# Joint significance of Department Fixed Effects
############################################################

library(fixest)
library(dplyr)

# ==========================================================
# 1. COMMON SAMPLE
# ==========================================================

# Department 2B contains only recognized observations
# and is dropped by the FE estimation.
# We remove it so that all models use exactly the same sample.

df_common <- df2002 %>%
  filter(departement != "2B")

nrow(df_common)
# Expected: 64386


# ==========================================================
# 2. MODEL WITH YEAR FE ONLY
# ==========================================================

r_year_test <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) |
    year,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df_common
)


# ==========================================================
# 3. MODEL WITH YEAR FE + EXPLICIT DEPARTMENT DUMMIES
# ==========================================================

# i(departement) creates explicit department coefficients,
# which allows us to test them jointly.

r_dept_test <- feglm(
  decision ~
    SWI +
    CSS +
    SWI:CSS +
    old_housing_share +
    share_individual_houses +
    log(n_dwellings) +
    i(departement) |
    year,
  family = binomial("probit"),
  cluster = ~code_insee,
  data = df_common
)


# ==========================================================
# 4. CHECK IDENTICAL SAMPLE SIZE
# ==========================================================

c(
  Year_FE = nobs(r_year_test),
  Year_Department_FE = nobs(r_dept_test)
)


# ==========================================================
# 5. JOINT WALD TEST OF DEPARTMENT EFFECTS
# ==========================================================

wald_dept <- wald(
  r_dept_test,
  keep = "departement::"
)

wald_dept

