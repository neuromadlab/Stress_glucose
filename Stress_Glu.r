# =============================================================================
# Higher glucose levels are associated with lower everyday stress load
# =============================================================================
# Description:  This script performs linear mixed-effects modeling and creates
#               publication-ready figures to examine the association between
#               continuous glucose monitoring (CGM) and ecological momentary
#               assessments (EMA) of stress and mood.
#
# Authors:     Madeleine Kördel, Anne Kühnel, Kristin Kaduk, Alica Guzman, Melanie Henes, Melina Grahlow, Birgit Derntl, Nils B. Kroemer
# Last updated: 2026-05-21
#
# Structure:
#   1. Setup & Package Loading
#   2. Data Import
#   3. Data Preprocessing
#   4. Sample Description & Demographic Tests
#   5. Correlational Analyses
#   6. Linear Mixed-Effects Models
#      6.1  Stress-Mood Association (Figure 2)
#      6.2  Glucose-Stress Association (Figure 3)
#      6.3  Temporal Glucose-Stress Dynamics (Figure 4)
#      6.4  Sensitivity Analyses
#      6.5  Age-Related Analyses
#   7. Distribution Comparisons (Kolmogorov-Smirnov Tests)
#   8. Figures (Main)
#      8.1  Figure 2 - Individual Stress Traces & Heatmap
#      8.2  Figure 3 - Glucose & Stress Density Plot
#      8.3  Figure 4 - Temporal Dynamics
#      8.4  Figure 5 - BMI & Stress Dynamics
#   9. Supplementary Figures
#      9.1  Figure S1 - Stress and Glucose across Meal Times
#      9.2  Figure S2 - Sensitivity: Pre-EMA Window Comparison
# =============================================================================


# =============================================================================
# 1. SETUP & PACKAGE LOADING
# =============================================================================

library(ggplot2)
library(foreign)
library(MASS)
library(cowplot)
library(viridis)
library(readxl)
library(tidyverse)
library(ggridges)
library(ggforce)
library(ggdist)
library(Hmisc)
library(ggside)
library(dplyr)
library(hexbin)
library(ggeffects)
library(lmerTest)
library(emmeans)
library(forcats)

# Set global plotting theme
theme_set(theme_cowplot(font_size = 12))


# =============================================================================
# 2. DATA IMPORT
# =============================================================================

# Main dataframe: EMA and CGM data merged at the observation level
d <- read_excel("TUE009_Glucose_VAS_Blood.xlsx")

# CGM data time-locked to EMA assessments (12h pre-EMA window, 1.5h bins)
dT <- read.csv("TUE009_GCM_VAS_Stress_12h.csv")

# EMA dataframe used for plotting only (complete EMA runs with glucose linkage)
dEMA <- read_csv("TUE009_Influenca_data_runs_clean_Glu_LR.csv")

# Sensitivity analysis: CGM data with alternative (variable) pre-EMA binning
dT_sens <- read.csv("TUE009_CGM_VAS_Stress.csv")

# Stress-eating questionnaire data (Stress and Snacking Effects Scale, SSES)
dQ <- read_excel("TUE009_SSES_questionnaire.xlsx")


# =============================================================================
# 3. DATA PREPROCESSING
# =============================================================================

# --- 3.1 Main dataframe (d) --------------------------------------------------

# Exclude participant with fewer than 20 EMA runs (insufficient data)
d <- d %>% filter(ID != 87)

# Grand-mean centering of continuous covariates
d$cHappy    <- d$Happy    - mean(d$Happy,    na.rm = TRUE)
d$cSad      <- d$Sad      - mean(d$Sad,      na.rm = TRUE)
d$cBMI      <- d$BMI      - mean(d$BMI,      na.rm = TRUE)
d$cAge      <- d$Age      - mean(d$Age,      na.rm = TRUE)
d$cMetState <- d$MetState - mean(d$MetState, na.rm = TRUE)  # hunger-satiety; used in sensitivity analyses

# Effect-code binary variables (subtract 0.5 to center around zero)
d$cSex <- d$Sex_male - 0.5

# Recent eating: L_Meal coding is 1=<0.5h, 2=~1h, 3=~1.5h, 4=~2h, 5=~2.5h, 6=>3h
# Classify as recent eating if last meal was within ~2h (still in postprandial window)
d$RecentEating  <- ifelse(d$L_Meal <= 4, 1, 0)
d$cRecentEating <- d$RecentEating - 0.5

# Within-person (group-mean) centering of glucose: deviation from participant's mean
d$M_Glu <- d$Glu - d$gcGlu

# BMI categories
d$fBMI <- cut(d$BMI,
              breaks = c(-100, 25, 30, Inf),
              labels = c("Normal weight", "Overweight", "Obese"))

# Temporal derivatives: change in stress, mood, and time between consecutive EMA runs
d <- d %>%
  arrange(ID, run_ind) %>%
  group_by(ID) %>%
  mutate(
    delStress = c(NA, diff(Stress)),      # change in stress (current - previous)
    delMood   = c(NA, diff(MoodState)),   # change in mood
    delTime   = c(NA, diff(timestamp))    # time elapsed (raw units -> converted below)
  ) %>%
  ungroup()

d$delTime <- d$delTime / 60000  # convert milliseconds to minutes

# Classify stress changes into categories for visualisation
d$fdelStress <- cut(d$delStress,
                    breaks = c(-100, -5, 5, Inf),
                    labels = c("Stress relief", "No change", "Stress onset"))

# Participant-level summaries: mean mood and stress, mood rank
dRank <- d %>%
  group_by(ID) %>%
  summarise(
    M_MoodState = mean(MoodState),
    M_Stress    = mean(Stress)
  ) %>%
  mutate(R_MoodState = rank(M_MoodState))  # rank participants by average mood

d <- merge(d, dRank, by = "ID")

# Scale mood variables to a 0-100 range and grand-mean center
d$MoodState   <- d$MoodState   * 100
d$cMoodState  <- d$MoodState   - mean(d$MoodState,  na.rm = TRUE)
d$M_MoodState <- d$M_MoodState * 100

# Z-score glucose and stress for supplementary Figure S1
d <- d %>%
  mutate(
    zGlu    = as.numeric(scale(Glu)),
    zStress = as.numeric(scale(Stress))
  )


# --- 3.2 CGM-EMA time-locked dataframe (dT) ----------------------------------

# Participants excluded for fewer than 20 EMA runs
exclude_ids <- c(9, 11, 13, 25, 27, 73, 87)
dT <- dT %>% filter(!id %in% exclude_ids)

# Rename variables for clarity
dT <- dT %>%
  rename(gl       = glucose,
         estIVglu = estIVglucose)

# Classify stress derivatives
dT$fdelStress <- cut(dT$StressDeriv,
                     breaks = c(-100, -5, 5, Inf),
                     labels = c("Stress relief", "No change", "Stress onset"))

# Grand-mean center the pre-EMA CGM index
dT$cIndEMA <- dT$IndexBefore - mean(dT$IndexBefore, na.rm = TRUE)

# Bin the pre-EMA CGM index into 8 quantile-based time bins (~1.5h each over 12h)
dT <- dT %>%
  mutate(bTime = ntile(IndexBefore, n = 8))

dT$fbTime <- factor(dT$bTime,
                    labels = c("b-8", "b-7", "b-6", "b-5", "b-4", "b-3", "b-2", "b-1"))

# Convert CGM index to hours before EMA
dT$iTime <- dT$IndexBefore / 12

# Extract hour of day from EMA timestamp and grand-mean center
dT$EMA_tT <- dmy_hms(dT$time)
dT$EMA_h  <- hour(dT$EMA_tT)
dT$cEMA_h <- dT$EMA_h - mean(dT$EMA_h, na.rm = TRUE)

# Collapse the two earliest time bins (b-8 and b-7) to ensure sufficient data
# density and model stability for the time bin furthest from EMA
dT$fbTime_col <- fct_collapse(dT$fbTime, "b-87" = c("b-8", "b-7"))


# --- 3.3 Sensitivity CGM-EMA dataframe (dT_sens) -----------------------------

dT_sens <- dT_sens %>% filter(!id %in% exclude_ids)

dT_sens$fdelStress <- cut(dT_sens$StressDeriv,
                          breaks = c(-100, -5, 5, Inf),
                          labels = c("Stress relief", "No change", "Stress onset"))

dT_sens$cIndEMA <- dT_sens$IndexBefore - mean(dT_sens$IndexBefore, na.rm = TRUE)

# Alternative binning: 6 quantile-based bins (~2h each over 12h)
dT_sens <- dT_sens %>%
  mutate(bTime = ntile(IndexBefore, n = 6))

dT_sens$fbTime <- factor(dT_sens$bTime,
                         labels = c("b-12", "b-10", "b-8", "b-6", "b-4", "b-2"))

dT_sens$iTime  <- dT_sens$IndexBefore / 12

dT_sens$EMA_tT <- dmy_hms(dT_sens$time)
dT_sens$EMA_h  <- hour(dT_sens$EMA_tT)
dT_sens$cEMA_h <- dT_sens$EMA_h - mean(dT_sens$EMA_h, na.rm = TRUE)


# --- 3.4 SSES (stress eating) questionnaire dataframe (dQ) -----------------------------------

# Retain only participants present in the main dataframe
dQ <- dQ %>% filter(ID %in% unique(d$ID))

# Rename for clarity
dQ <- dQ %>% rename(SSES = SSES_mean)

# Descriptive: categorise participants by stress-eating behaviour
dQ %>%
  mutate(eating_behavior = case_when(
    SSES < 3  ~ "eats less",
    SSES == 3 ~ "no change",
    SSES > 3  ~ "eats more"
  )) %>%
  count(eating_behavior) %>%
  mutate(percent = n / sum(n) * 100)

# Merge SSES scores into main dataframe (inner join -> only IDs present in dQ)
# Note: run this after all d preprocessing steps are complete
d_sses <- d %>% inner_join(dQ %>% select(ID, SSES), by = "ID")

# Grand-mean center SSES
d_sses$cSSES <- d_sses$SSES - mean(d_sses$SSES, na.rm = TRUE)


# --- 3.5 EMA plotting dataframe (dEMA) ---------------------------------------

# Retain only participants present in the main dataframe
dEMA <- dEMA %>% filter(glucose_id %in% unique(d$ID))
dEMA <- dEMA %>% rename(ID = glucose_id)

# Composite mood state: positive minus negative affect
dEMA <- dEMA %>% mutate(MoodState = Happy - Sad)

# Temporal derivatives
dEMA <- dEMA %>%
  arrange(ID, run_ind) %>%
  group_by(ID) %>%
  mutate(
    delStress = c(NA, diff(Stress)),
    delMood   = c(NA, diff(MoodState)),
    delTime   = c(NA, diff(timestamp))
  ) %>%
  ungroup()

dEMA$delTime <- dEMA$delTime / 60000

# Participant-level mood ranks
dRank_EMA <- dEMA %>%
  group_by(ID) %>%
  summarise(
    M_MoodState = mean(MoodState),
    M_Stress    = mean(Stress)
  ) %>%
  mutate(R_MoodState = rank(M_MoodState))

dEMA <- merge(dEMA, dRank_EMA, by = "ID")

# Scale mood to 0-100 and grand-mean center
dEMA$MoodState   <- dEMA$MoodState   * 100
dEMA$cMoodState  <- dEMA$MoodState   - mean(dEMA$MoodState,  na.rm = TRUE)
dEMA$M_MoodState <- dEMA$M_MoodState * 100

# =============================================================================
# 4. SAMPLE DESCRIPTION & DEMOGRAPHIC TESTS
# =============================================================================

# Overall sample characteristics (one row per participant)
d %>%
  group_by(ID) %>%
  summarise(
    Age      = first(Age),
    BMI      = first(BMI),
    Sex_male = first(Sex_male),
    N_runs   = max(run_ind, na.rm = TRUE)
  ) %>%
  summarise(
    Mean_Age  = round(mean(Age,    na.rm = TRUE), 2),
    SD_Age    = round(sd(Age,      na.rm = TRUE), 2),
    Min_Age   = round(min(Age,     na.rm = TRUE), 2),
    Max_Age   = round(max(Age,     na.rm = TRUE), 2),
    Mean_BMI  = round(mean(BMI,    na.rm = TRUE), 2),
    SD_BMI    = round(sd(BMI,      na.rm = TRUE), 2),
    Min_BMI   = round(min(BMI,     na.rm = TRUE), 2),
    Max_BMI   = round(max(BMI,     na.rm = TRUE), 2),
    Mean_runs = round(mean(N_runs, na.rm = TRUE), 2),
    SD_runs   = round(sd(N_runs,   na.rm = TRUE), 2),
    N_female  = sum(Sex_male == 0),
    N_male    = sum(Sex_male == 1)
  )

# Demographic differences between sexes
d_subj_demo <- d %>%
  group_by(ID) %>%
  summarise(
    Age      = first(Age),
    BMI      = first(BMI),
    Sex_male = first(Sex_male)
  )

t.test(Age ~ Sex_male, data = d_subj_demo)
t.test(BMI ~ Sex_male, data = d_subj_demo)


# =============================================================================
# 5. CORRELATIONAL ANALYSES
# =============================================================================
# Related to Figure 2: everyday stress is associated with mood, but
# acute stress spikes occur independently of momentary mood

cor.test(d$Stress,    d$MoodState)    # momentary stress ~ momentary mood
cor.test(d$delStress, d$M_MoodState)  # stress change ~ person-mean mood
cor.test(d$delStress, d$delMood)      # stress change ~ mood change
cor.test(d$M_Stress,  d$M_MoodState)  # person-mean stress ~ person-mean mood


# =============================================================================
# 6. LINEAR MIXED-EFFECTS MODELS
# =============================================================================

# =============================================================================
# 6.1 Stress-Mood Association (Figure 2)
# =============================================================================

# Model 1a: Composite mood state (Happy - Sad) predicts momentary stress
fm_1a <- lmer(Stress ~ cMoodState + cBMI + cAge + cSex +
                (1 + cMoodState | ID),
              data = d, REML = TRUE)
summary(fm_1a)

# Model 1b: Happy and Sad rated separately to decompose composite mood
fm_1b <- lmer(Stress ~ cHappy + cSad + cBMI + cAge + cSex +
                (1 + cHappy + cSad | ID),
              data = d, REML = TRUE)
summary(fm_1b)

# Model 1c: Stress derivative predicted by person-mean mood and mood change
fm_1c <- lmer(delStress ~ scale(M_MoodState) + delMood + cBMI + cAge + cSex +
                (1 + scale(M_MoodState) + delMood | ID),
              data = d, REML = TRUE)
summary(fm_1c)

# Model 1d: As Model 1c but using centered momentary mood instead of person-mean
fm_1d <- lmer(delStress ~ cMoodState + delMood + cBMI + cAge + cSex +
                (1 + cMoodState + delMood | ID),
              data = d, REML = TRUE)
summary(fm_1d)


# =============================================================================
# 6.2 Glucose-Stress Association (Figure 3)
# =============================================================================

# Model 2res: Residualise stress for mood (residuals used in Figure 3a)
fm_2res <- lmer(Stress ~ (cHappy + cSad) + cBMI + cAge +
                  (1 + cHappy + cSad | ID),
                data = d)
summary(fm_2res)
d$lmeResStress <- residuals(fm_2res)

# Model 2a: Glucose (within-person and between-person) x mood interaction predicts stress
fm_2a <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex +
                (1 + cHappy + cSad + gcGlu | ID),
              data = d)
summary(fm_2a)

# Effect size and confidence interval for the within-person glucose effect
confint(fm_2a, parm = "gcGlu", method = "Wald")
t  <- -2.141
df <- 720
r  <- sign(t) * sqrt(t^2 / (t^2 + df))
r

# Model-predicted glucose-stress association for plotting
dPfm2a <- ggpredict(fm_2a, terms = "gcGlu")

# Extract fixed-effect estimates and CIs
coef_tab <- coef(summary(fm_2a))
coef_tab["gcGlu", ]
confint(fm_2a, parm = "gcGlu", method = "Wald")

# Empirical Bayes (EB) estimates of participant-level glucose slopes (Figure 3b)
d_fm2a <- coef(fm_2a)$ID
sum(d_fm2a$gcGlu < 0)           # number of participants with negative slopes
nrow(d_fm2a)                    # total participants
mean(d_fm2a$gcGlu < 0) * 100   # percentage with negative slopes

# Individual OLS slopes (unshrunken) for comparison with EB estimates
ols_fits  <- lme4::lmList(Stress ~ (cHappy + cSad) * gcGlu | ID, data = d)
ols_coefs <- coef(ols_fits)
mean(ols_coefs[, "gcGlu"] < 0, na.rm = TRUE)  # proportion with negative OLS slopes

# Model 2b: Add sex x glucose interaction
fm_2b <- lmer(Stress ~ (cHappy + cSad + cSex) * gcGlu + M_Glu + cBMI + cAge + cSex +
                (1 + cHappy + cSad + gcGlu | ID),
              data = d)
summary(fm_2b)


# =============================================================================
# 6.3 Temporal Glucose-Stress Dynamics (Figure 4)
# =============================================================================

# Set the collapsed earliest bin (b-87) as the reference category
dT <- within(dT, fbTime_col <- relevel(fbTime_col, ref = 7))

# Model 3a: Glucose x time bin interaction predicting stress derivative (12h window, 1.5h bins)
fm_3a <- lmer(estIVglu ~ fbTime_col * scale(StressDeriv) +
                (1 + fbTime_col + scale(StressDeriv) | id),
              data = dT)
summary(fm_3a)

# Model 3b: As Model 3a, additionally controlling for time of day
fm_3b <- lmer(estIVglu ~ fbTime_col * scale(StressDeriv) + cEMA_h +
                (1 + fbTime_col + scale(StressDeriv) + cEMA_h | id),
              data = dT,
              control = lmerControl(
                optimizer    = "nloptwrap",
                calc.derivs  = FALSE,
                check.conv.grad     = .makeCC("warning", tol = 2e-3),
                check.conv.singular = .makeCC("ignore",  tol = 1e-4),
                optCtrl      = list(maxfun = 1e6)
              ))
summary(fm_3b)

# Model 3c: Test whether time of day is associated with stress derivative (Figure 4d)
# Restricted to the measurement immediately before EMA (IndexBefore == -1)
fm_3c <- lmer(EMA_h ~ scale(StressDeriv) +
                (1 + scale(StressDeriv) | id),
              data = subset(dT, IndexBefore == -1))
summary(fm_3c)

# Control analysis: test whether effects are independent of the actual
# time distance between CGM measurements and EMA assessments
fm_4a <- lmer(estIVglu ~ cIndEMA * scale(StressDeriv) +
                (1 + cIndEMA | id) + (1 | Run),
              data = dT)
summary(fm_4a)


# =============================================================================
# 6.4 Sensitivity Analyses
# =============================================================================

# --- Metabolic state (hunger-satiety) as alternative predictor ---------------

# Model 2c: Add metabolic state as additional covariate alongside glucose
fm_2c <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + cMetState +
                (1 + cHappy + cSad + gcGlu + cMetState | ID),
              data = d)
summary(fm_2c)

# Model 2d: Replace glucose with metabolic state as the within-person predictor
fm_2d <- lmer(Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cBMI + cAge + cSex +
                (1 + cHappy + cSad + gcMetState | ID),
              data = d)
summary(fm_2d)

# Model comparison: glucose vs. metabolic state
data.frame(
  model     = c("fm_2a", "fm_2d"),
  df        = AIC(fm_2a, fm_2d)$df,
  AIC       = AIC(fm_2a, fm_2d)$AIC,
  delta_AIC = AIC(fm_2a, fm_2d)$AIC - min(AIC(fm_2a, fm_2d)$AIC),
  BIC       = BIC(fm_2a, fm_2d)$BIC,
  delta_BIC = BIC(fm_2a, fm_2d)$BIC - min(BIC(fm_2a, fm_2d)$BIC)
)

# Check whether glucose and metabolic state are correlated
cor.test(d$Glu, d$MetState)

# --- Recent eating -----------------------------------------------------------

# Model 2a_eating: Recent eating as additional covariate
fm_2a_eating <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + cRecentEating +
                       (1 + cHappy + cSad + gcGlu | ID),
                     data = d)
summary(fm_2a_eating)

# Model 2b_eating: Recent eating x glucose interaction
fm_2b_eating <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + cRecentEating * gcGlu +
                       (1 + cHappy + cSad + gcGlu | ID),
                     data = d)
summary(fm_2b_eating)

# Model 2c_eating: Recent eating added to metabolic state model
fm_2c_eating <- lmer(Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cBMI + cAge + cSex + cRecentEating +
                       (1 + cHappy + cSad + gcMetState | ID),
                     data = d)
summary(fm_2c_eating)

# Model 2d_eating: Recent eating x metabolic state interaction
fm_2d_eating <- lmer(Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cBMI + cAge + cSex + cRecentEating * gcMetState +
                       (1 + cHappy + cSad + gcMetState | ID),
                     data = d)
summary(fm_2d_eating)

# Likelihood-ratio tests for recent eating
anova(fm_2a, fm_2a_eating)
anova(fm_2c, fm_2c_eating)

# --- Stress-eating behavior (SSES) -------------------------------------------

# Model 2a_sses: Replicate Model 2a on the SSES subsample
fm_2a_sses <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex +
                     (1 + cHappy + cSad + gcGlu | ID),
                   data = d_sses)
summary(fm_2a_sses)

# Model 2b_sses: Add centered SSES as covariate
fm_2b_sses <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + cSSES +
                     (1 + cHappy + cSad + gcGlu | ID),
                   data = d_sses)
summary(fm_2b_sses)

# Model 2c_sses: SSES x glucose interaction
fm_2c_sses <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + cSSES * gcGlu +
                     (1 + cHappy + cSad + gcGlu | ID),
                   data = d_sses)
summary(fm_2c_sses)

# Model 2d_sses: Replicate Model 2d (metabolic state) on SSES subsample
fm_2d_sses <- lmer(Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cBMI + cAge + cSex +
                     (1 + cHappy + cSad + gcMetState | ID),
                   data = d_sses)
summary(fm_2d_sses)

# Model 2e_sses: Add centered SSES to metabolic state model
fm_2e_sses <- lmer(Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cBMI + cAge + cSex + cSSES +
                     (1 + cHappy + cSad + gcMetState | ID),
                   data = d_sses)
summary(fm_2e_sses)

# Model 2f_sses: SSES x metabolic state interaction
fm_2f_sess <- lmer(Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cBMI + cAge + cSex + cSSES * gcMetState +
                     (1 + cHappy + cSad + gcMetState | ID),
                   data = d_sses)
summary(fm_2f_sess)

# Likelihood-ratio tests for SSES
anova(fm_2a_sses, fm_2b_sses)
anova(fm_2d_sses, fm_2e_sses)

# Cross-model comparisons
AIC(fm_2b_sses, fm_2e_sses)
BIC(fm_2b_sses, fm_2e_sses)

# Person-level residual stress ~ SSES correlation
d_sses$lmeResStress_sses <- residuals(fm_2a_sses)

res_person <- d_sses %>%
  group_by(ID) %>%
  summarise(M_ResStress = mean(lmeResStress_sses, na.rm = TRUE)) %>%
  left_join(dQ %>% select(ID, SSES), by = "ID")

cor.test(res_person$M_ResStress, res_person$SSES)

# --- Metabolic parameters as moderators --------------------------------------

# HOMA-IR as moderator of glucose x stress coupling
fm_HOMA <- lmer(Stress ~ (cHappy + cSad) * gcGlu * log_HOMA + M_Glu + cBMI + cAge + cSex +
                  (1 + cHappy + cSad + gcGlu | ID),
                data = d)
summary(fm_HOMA)

# Triglyceride-glucose index (TyG) as moderator
fm_TyG <- lmer(Stress ~ (cHappy + cSad) * gcGlu * cTyG + M_Glu + cBMI + cAge + cSex +
                 (1 + cHappy + cSad + gcGlu | ID),
               data = d)
summary(fm_TyG)

# Between-person associations: do metabolic parameters explain glucose-stress slopes?
d_subj <- d %>%
  group_by(ID) %>%
  summarise(
    cBMI     = first(cBMI),
    log_HOMA = first(log_HOMA),
    cTyG     = first(cTyG)
  )

d_subj$gcGlu_slope <- d_fm2a$gcGlu  # EB glucose slopes from Model 2a

cor.test(d_subj$gcGlu_slope, d_subj$cBMI)      # BMI
cor.test(d_subj$gcGlu_slope, d_subj$log_HOMA)  # HOMA-IR
cor.test(d_subj$gcGlu_slope, d_subj$cTyG)      # TyG index

# --- Alternative pre-EMA window (sensitivity for Figure S2) ------------------

# Set the closest bin (0-2h before EMA) as reference
dT_sens <- within(dT_sens, fbTime <- relevel(fbTime, ref = 6))

# Model 3a_sens: Glucose x time bin interaction (variable 2h bins)
fm_3a_sens <- lmer(estIVglu ~ fbTime * scale(StressDeriv) +
                     (1 + fbTime + scale(StressDeriv) | id),
                   data = dT_sens)
summary(fm_3a_sens)

# Model 3b_sens: As above, controlling for time of day
fm_3b_sens <- lmer(estIVglu ~ fbTime * scale(StressDeriv) + cEMA_h +
                     (1 + fbTime + scale(StressDeriv) + cEMA_h | id),
                   data = dT_sens)
summary(fm_3b_sens)

# Model 3c_sens: Time of day ~ stress derivative
fm_3c_sens <- lmer(EMA_h ~ scale(StressDeriv) +
                     (1 + scale(StressDeriv) | id),
                   data = subset(dT_sens, IndexBefore == -1))
summary(fm_3c_sens)


# =============================================================================
# 6.5 Age-Related Analyses
# =============================================================================

# --- Mood-stress association moderated by age --------------------------------

fm_1a_age <- lmer(Stress ~ cMoodState * cAge + cBMI + cSex +
                    (1 + cMoodState | ID),
                  data = d, REML = TRUE)
summary(fm_1a_age)

fm_1b_age <- lmer(Stress ~ cHappy * cAge + cSad * cAge + cBMI + cSex +
                    (1 + cHappy + cSad | ID),
                  data = d, REML = TRUE)
summary(fm_1b_age)

fm_1c_age <- lmer(delStress ~ scale(M_MoodState) * cAge + delMood * cAge + cBMI + cSex +
                    (1 + scale(M_MoodState) + delMood | ID),
                  data = d, REML = TRUE)
summary(fm_1c_age)

fm_1d_age <- lmer(delStress ~ (cMoodState + delMood) * cAge + cBMI + cSex +
                    (1 + cMoodState + delMood | ID),
                  data = d)
summary(fm_1d_age)

# --- Glucose-stress coupling moderated by age --------------------------------

fm_2a_age <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cAge * gcGlu + cBMI + cSex +
                    (1 + cHappy + cSad + gcGlu | ID),
                  data = d)
summary(fm_2a_age)

fm_2d_age <- lmer(Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cAge * gcMetState + cBMI + cAge + cSex +
                    (1 + cHappy + cSad + gcMetState | ID),
                  data = d)
summary(fm_2d_age)

# --- Between-person correlations with age ------------------------------------

d_subj <- d %>%
  group_by(ID) %>%
  summarise(
    Age       = first(Age),
    cAge      = first(cAge),
    BMI       = first(BMI),
    M_Glu     = mean(Glu,    na.rm = TRUE),
    SD_Glu    = sd(Glu,      na.rm = TRUE),
    M_Stress  = mean(Stress, na.rm = TRUE),
    SD_Stress = sd(Stress,   na.rm = TRUE),
    log_HOMA  = first(log_HOMA),
    cTyG      = first(cTyG)
  )

cor.test(d_subj$Age, d_subj$M_Glu)      # age ~ mean glucose
cor.test(d_subj$Age, d_subj$M_Stress)   # age ~ mean stress
cor.test(d_subj$Age, d_subj$SD_Stress)  # age ~ stress variability
cor.test(d_subj$Age, d_subj$SD_Glu)     # age ~ glucose variability
cor.test(d_subj$Age, d_subj$log_HOMA)   # age ~ HOMA-IR
cor.test(d_subj$Age, d_subj$cTyG)       # age ~ TyG index

# Between-person models: age and BMI predicting stress/glucose variability
fm_age_stressSD <- lm(SD_Stress ~ Age + BMI, data = d_subj)
summary(fm_age_stressSD)

fm_age_gluSD <- lm(SD_Glu ~ Age + BMI, data = d_subj)
summary(fm_age_gluSD)

# Does stress variability moderate the glucose-stress coupling?
d <- d %>%
  left_join(d_subj %>% select(ID, SD_Stress), by = "ID") %>%
  mutate(cSD_Stress = SD_Stress - mean(SD_Stress, na.rm = TRUE))

fm_stressvar_glu <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu +
                           cSD_Stress * gcGlu +
                           cBMI + cAge + cSex +
                           (1 + cHappy + cSad + gcGlu | ID),
                         data = d)
summary(fm_stressvar_glu)

# Does age moderate the effect of insulin resistance on glucose-stress coupling?
fm_HOMA_age <- lmer(Stress ~ (cHappy + cSad) * gcGlu * log_HOMA + cAge * log_HOMA +
                      M_Glu + cBMI + cSex +
                      (1 + cHappy + cSad + gcGlu | ID),
                    data = d)
summary(fm_HOMA_age)

fm_TyG_age <- lmer(Stress ~ (cHappy + cSad) * gcGlu * cTyG + cAge * cTyG +
                     M_Glu + cBMI + cSex +
                     (1 + cHappy + cSad + gcGlu | ID),
                   data = d)
summary(fm_TyG_age)


# =============================================================================
# 7. DISTRIBUTION COMPARISONS (Kolmogorov-Smirnov Tests)
# =============================================================================
# Related to Figure 5: stress dynamics across BMI categories

# Subset by BMI category
d_NW <- filter(d, fBMI == "Normal weight")
d_OW <- filter(d, fBMI == "Overweight")
d_OB <- filter(d, fBMI == "Obese")

# Number of participants per category
d %>%
  group_by(fBMI) %>%
  summarise(n = n_distinct(ID))

# Pairwise KS tests on stress-change distributions
ks.test(d_NW$delStress, d_OW$delStress,
        exact = TRUE, simulate.p.value = TRUE, B = 2000)  # Normal weight vs. Overweight
ks.test(d_NW$delStress, d_OB$delStress,
        exact = TRUE, simulate.p.value = TRUE, B = 2000)  # Normal weight vs. Obese
ks.test(d_OW$delStress, d_OB$delStress,
        exact = TRUE, simulate.p.value = TRUE, B = 2000)  # Overweight vs. Obese

# =============================================================================
# 8. FIGURES (MAIN)
# =============================================================================
# Note: Ensure all models in Section 6 have been fitted before running figures.
# Figures are saved to the /figures sub-directory.

# =============================================================================
# 8.1 Figure 2 - Individual Stress Traces & Heatmap
# =============================================================================

# Panel a: Ridgeline plot of individual stress traces, coloured by mean mood
plot2a <-
  ggplot(aes(x = run_ind, y = R_MoodState), data = d) +
  geom_ridgeline(aes(height = Stress / 15, color = M_MoodState, group = ID),
                 linewidth = 0.7, alpha = 0.01) +
  scale_color_gradient2("Mood",
                        low  = "steelblue",
                        mid  = "grey80",
                        high = "tomato",
                        midpoint = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(20, 80, 20)) +
  theme(
    text              = element_text(face = "bold",  size = 20.0),
    axis.text         = element_text(face = "plain", size = 18.0),
    axis.text.x       = element_text(size = 18.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 20.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 18.0, hjust = 0.5),
    legend.key.height = unit(2, "cm")
  ) +
  ggtitle("Individual traces of stress ratings", "EMA, rescaled for display") +
  xlab("Runs") +
  ylab("Mood rank [VAS happy - sad]")

# Panel b: Heatmap of stress ratings per run, sorted by mood rank
plot2b <-
  ggplot(dEMA, aes(x = run_ind, y = R_MoodState, fill = Stress)) +
  geom_tile() +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    limits    = c(0, 110),
    values    = scales::rescale(c(10, 30, 65, 100)),
    name      = "Stress",
    na.value  = "grey80"
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 60, 10)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(20, 80, 20)) +
  theme(
    text              = element_text(face = "bold",  size = 20.0),
    axis.text         = element_text(face = "plain", size = 18.0),
    axis.text.x       = element_text(size = 18.0),
    axis.line.x       = element_line(color = "black", linewidth = 0.4),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 20.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 18.0, hjust = 0.5),
    legend.key.height = unit(2, "cm"),
    panel.grid        = element_blank()
  ) +
  ggtitle("EMA stress ratings per run") +
  xlab("Runs") +
  ylab("")

plot2 <- plot_grid(plot2a, plot2b,
                   labels = c("a", "b"), label_size = 18,
                   ncol = 2, rel_widths = c(1.5, 1.5), align = "h")

ggsave("figures/Figure_2_Stress_Traces_Heatmap.png",
       plot = plot2, height = 7, width = 12, units = "in", dpi = 600, bg = "white")


# =============================================================================
# 8.2 Figure 3 - Glucose & Stress Density Plot
# =============================================================================

# Panel a: Density plot with model predictions
plot3a <-
  ggplot(dPfm2a, aes(x, predicted)) +
  stat_density_2d(aes(x = gcGlu, y = lmeResStress + 35.72,
                      fill = after_stat(ndensity)),
                  geom = "raster", contour = FALSE, n = 500, data = d) +
  geom_line(linewidth = 2, color = "slateblue") +
  scale_fill_viridis_c(guide = guide_legend(title = "Density"), option = "F") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "grey80", alpha = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(30, 45), xlim = c(-30, 40)) +
  theme(
    text         = element_text(face = "bold",  size = 14.0),
    axis.text    = element_text(face = "plain", size = 14.0),
    axis.text.x  = element_text(size = 12.0),
    strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))
  ) +
  ylab("Stress [residualised for mood]") +
  xlab("Glucose levels [group centered]")

# Panel b: Distribution of EB glucose-stress slopes across participants
plot3b <-
  ggplot(d_fm2a, aes(x = 1, y = gcGlu)) +
  geom_bar(linewidth = 1, color = "grey80", fill = "slateblue", stat = "summary") +
  geom_sina(size = 2.5, adjust = 0.5) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 1) +
  theme(
    text          = element_text(face = "bold",  size = 14.0),
    axis.text     = element_text(face = "plain", size = 14.0),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    strip.text.x  = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))
  ) +
  ylab("EB estimate") +
  xlab("Glu slope")

plot3 <- plot_grid(plot3a, plot3b,
                   labels = c("a", "b"), label_size = 12,
                   ncol = 2, rel_widths = c(1.40, 0.60), align = "h")

ggsave("figures/Figure_3_Glucose_Stress.png",
       plot = plot3, height = 4.5, width = 6.5, units = "in", dpi = 600, bg = "white")


# =============================================================================
# 8.3 Figure 4 - Temporal Dynamics of Glucose before Stress Changes
# =============================================================================

# Panel a: GAM-smoothed glucose trajectories stratified by stress-change category
plot4a <-
  ggplot(aes(y = estIVglu, x = iTime),
         data = subset(dT, !is.na(fdelStress))) +
  # Highlight significant interaction bin (-6 to -7.5h)
  annotate("rect", xmin = -7.5, xmax = -6, ymin = 103, ymax = 109,
           alpha = 0.06, fill = "grey70", color = "grey90") +
  # Highlight reference bin (-1.5 to 0h; main effect)
  annotate("rect", xmin = -1.5, xmax = 0, ymin = 103, ymax = 109,
           alpha = 0.8, fill = "#e8f4f8") +
  stat_smooth(aes(fill = fdelStress), linewidth = 1.2,
              method = "gam", alpha = 0.8, color = "grey30") +
  stat_summary_bin(aes(group = fdelStress, color = fdelStress),
                   fun.data = "mean_cl_boot",
                   geom = "linerange", linewidth = 1.5, alpha = 0.8, bins = 12) +
  scale_fill_manual(guide  = guide_legend(title = "Stress bin"),
                    breaks = c("Stress onset", "No change", "Stress relief"),
                    values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_color_manual(guide  = guide_legend(title = "Stress bin"),
                     breaks = c("Stress onset", "No change", "Stress relief"),
                     values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_x_continuous(breaks = seq(-12, 0, by = 1.5), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-12, 0), ylim = c(103, 109)) +
  annotate("text",    x = -6.75, y = 108.5, label = "*",           size = 6) +
  annotate("segment", x = -1.6,  y = 108.6, xend = -3.6, yend = 108.6,
           arrow = arrow(type = "closed", length = unit(0.03, "npc"))) +
  annotate("text",    x = -3.6,  y = 109,   label = "Interaction", size = 4) +
  annotate("text",    x = -0.75, y = 104.5, label = "*",           size = 6) +
  theme(
    text              = element_text(face = "bold",  size = 12.0),
    axis.text         = element_text(face = "plain", size = 12.0),
    axis.text.x       = element_text(size = 12.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 13.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 12.0, hjust = 0.5),
    legend.key.height = unit(1, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.text       = element_text(face = "bold", size = 10.0)
  ) +
  ggtitle("Differences in glucose", "before changes in stress (controlled for time)") +
  xlab("Time before EMA [h]") +
  ylab("Glucose level [mg/dl]")

# Panel b: Density distribution of stress-change categories
plot4b <-
  ggplot(aes(y = delStress, fill = fdelStress), data = d) +
  stat_slabinterval(aes(fill = fdelStress),
                    alpha = 1, color = "grey20", adjust = 2) +
  scale_fill_manual(guide  = guide_legend(title = "Stress cat"),
                    values = c("#96f5c3", "#faed7c", "#f5c396")) +
  coord_cartesian(ylim = c(-100, 100), xlim = c(0, 0.75)) +
  theme(
    legend.position   = "none",
    text              = element_text(face = "bold",  size = 12.0),
    axis.text         = element_text(face = "plain", size = 12.0),
    axis.text.x       = element_text(size = 10.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 13.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 12.0, hjust = 0.5),
    legend.key.height = unit(1, "cm")
  ) +
  ggtitle("Density") +
  xlab("Density") +
  ylab("Change in stress ratings")

# Panel c: Early glucose differences predict later stress changes
plot4c <-
  ggplot(aes(y = StressDeriv, x = estIVglu),
         data = subset(dT, bTime == 1)) +
  stat_density_2d(geom = "raster", aes(fill = after_stat(ndensity)),
                  contour = FALSE, n = 500) +
  stat_smooth(linewidth = 0.5, aes(group = id),
              method = "rlm", alpha = 0.01, color = "grey70") +
  stat_smooth(linewidth = 2, method = "rlm", alpha = 0.05) +
  scale_fill_viridis_c(guide = guide_legend(title = "Density"), option = "F") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-25, 25), xlim = c(75, 135)) +
  theme(
    legend.position   = "none",
    text              = element_text(face = "bold",  size = 12.0),
    axis.text         = element_text(face = "plain", size = 12.0),
    axis.text.x       = element_text(size = 10.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 13.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 12.0, hjust = 0.5),
    legend.key.height = unit(1, "cm")
  ) +
  ggtitle("Early differences in glucose", "precipitate changes in stress") +
  ylab("Changes in stress [deriv]") +
  xlab("Glucose levels [12 to 10h before EMA]")

# Panel d: Time of day associated with stress derivative
plot4d <-
  ggplot(aes(x = StressDeriv, y = EMA_h),
         data = subset(dT, IndexBefore == -1)) +
  stat_smooth(color = "cornflowerblue", fill = "#90D5FF", linewidth = 2,
              method = "gam", alpha = 0.5) +
  stat_summary_bin(color = "cornflowerblue", fun.data = "mean_cl_boot",
                   geom = "linerange", linewidth = 1.5, alpha = 0.8, bins = 10) +
  scale_x_continuous(expand = c(0.01, 5)) +
  coord_cartesian(ylim = c(10, 20), xlim = c(-100, 100)) +
  theme(
    legend.position   = "none",
    text              = element_text(face = "bold",  size = 12.0),
    axis.text         = element_text(face = "plain", size = 12.0),
    axis.text.x       = element_text(size = 10.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 13.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 12.0, hjust = 0.5),
    legend.key.height = unit(1, "cm")
  ) +
  ggtitle("Time of day") +
  ylab("Time [h]") +
  xlab("Changes in stress")

plot4 <- plot_grid(plot4a, plot4b, plot4c, plot4d,
                   labels = c("a", "b", "c", "d"), label_size = 12,
                   ncol = 2, rel_widths = c(1.39, 0.68), align = "h")

ggsave("figures/Figure_4.png",
       plot = plot4, height = 7.3, width = 7.3, units = "in", dpi = 600, bg = "white")

# =============================================================================
# 8.4 Figure 5 - BMI Categories & Stress Dynamics
# =============================================================================

plot5 <-
  ggplot(aes(y = delStress, x = fBMI), data = d) +
  stat_slabinterval(aes(group = fBMI, fill = fBMI), alpha = 0.9) +
  scale_fill_manual(guide  = guide_legend(title = "BMI cat"),
                    values = c("#CFE0C3", "#9EC1A3", "#70A9A1")) +
  theme(
    text              = element_text(face = "bold",  size = 14.0),
    axis.text         = element_text(face = "plain", size = 14.0),
    axis.text.x       = element_text(size = 12.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 16.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 14.0, hjust = 0.5),
    legend.key.height = unit(1.5, "cm")
  ) +
  ggtitle("Altered distribution of changes in stress", "Temporal derivative") +
  xlab("BMI categories") +
  ylab("Delta stress [VAS]")

ggsave("figures/Figure_5_BMI_Stress.png",
       plot = plot5, height = 4, width = 7, units = "in", dpi = 600, bg = "white")


# =============================================================================
# 9. SUPPLEMENTARY FIGURES
# =============================================================================

# =============================================================================
# 9.1 Figure S1 - Stress and Glucose across Meal Times (z-scores)
# =============================================================================

meal_labels <- c("1" = "<0.5h", "2" = "~1h", "3" = "~1.5h",
                 "4" = "~2h",   "5" = "~2.5h", "6" = ">3h")

# Person-level means per meal-timing category
df_person <- d %>%
  filter(!is.na(L_Meal), !is.na(zStress), !is.na(zGlu)) %>%
  group_by(ID, L_Meal) %>%
  summarise(
    person_zStress = mean(zStress, na.rm = TRUE),
    person_zGlu    = mean(zGlu,    na.rm = TRUE),
    .groups = "drop"
  )

# Grand means and SEM across participants
df_sum <- df_person %>%
  group_by(L_Meal) %>%
  summarise(
    mean_zStress = mean(person_zStress, na.rm = TRUE),
    se_zStress   = sd(person_zStress,  na.rm = TRUE) / sqrt(n()),
    mean_zGlu    = mean(person_zGlu,   na.rm = TRUE),
    se_zGlu      = sd(person_zGlu,     na.rm = TRUE) / sqrt(n()),
    n_persons    = n(),
    .groups = "drop"
  ) %>%
  mutate(L_Meal = factor(L_Meal, levels = 1:6, labels = meal_labels))

plot_S1 <-
  ggplot(df_sum, aes(x = as.numeric(L_Meal))) +
  # Shade postprandial window (<=2h since last meal)
  annotate("rect", xmin = 1, xmax = 4, ymin = -Inf, ymax = Inf,
           fill = "#b5e0b5", alpha = 0.25) +
  annotate("text", x = 2.5, y = Inf, vjust = 1.5,
           label = "Recent eating (<=2h)", size = 3.2, colour = "#4a7a4a") +
  # Glucose line
  geom_line(aes(y = mean_zGlu, group = 1, colour = "Glucose"),
            linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = mean_zGlu, colour = "Glucose"), size = 2.5) +
  geom_errorbar(aes(ymin = mean_zGlu - se_zGlu,
                    ymax = mean_zGlu + se_zGlu, colour = "Glucose"),
                width = 0.15, linewidth = 0.7) +
  # Stress line
  geom_line(aes(y = mean_zStress, group = 1, colour = "Stress"),
            linewidth = 0.9) +
  geom_point(aes(y = mean_zStress, colour = "Stress"), size = 2.5) +
  geom_errorbar(aes(ymin = mean_zStress - se_zStress,
                    ymax = mean_zStress + se_zStress, colour = "Stress"),
                width = 0.15, linewidth = 0.7) +
  scale_colour_manual(values = c("Glucose" = "#2980B9", "Stress" = "#C0392B"),
                      name = NULL) +
  scale_y_continuous(
    name     = "Glucose level\nMean +/- SEM across participants (z-score)",
    sec.axis = dup_axis(name = "Stress rating\nMean +/- SEM across participants (z-score)"),
    expand   = c(0.0, 0.01)
  ) +
  scale_x_continuous(breaks = 1:6, labels = levels(df_sum$L_Meal),
                     expand = c(0.0, 0.00)) +
  coord_cartesian(xlim = c(0.5, 6.5)) +
  labs(x = "Time since last meal",
       title = "Stress and glucose across meal times") +
  theme(
    text               = element_text(face = "bold",  size = 12.0),
    axis.text          = element_text(face = "plain", size = 12.0),
    axis.text.x        = element_text(size = 12.0),
    axis.title.y.left  = element_text(colour = "#2980B9", face = "bold", size = 12.0,
                                      margin = margin(r = 8)),
    axis.title.y.right = element_text(colour = "#C0392B", face = "bold", size = 12.0,
                                      margin = margin(l = 8)),
    strip.text.x       = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title         = element_text(size = 13.0, hjust = 0.5),
    legend.key.width   = unit(1.2, "cm"),
    legend.text        = element_text(face = "bold", size = 10.0),
    legend.position    = "top",
    legend.direction   = "horizontal"
  )

ggsave("figures/Figure_S1_meal_timing.png",
       plot = plot_S1, height = 4.5, width = 6.5, units = "in", dpi = 600, bg = "white")


# =============================================================================
# 9.2 Figure S2 - Sensitivity: Pre-EMA Window Comparison
# =============================================================================
# Panels a and b replicate Figure 4a using two different CGM binning approaches
# to demonstrate robustness of the temporal glucose-stress interaction.

# Compute actual bin extents for the sensitivity dataset (dT_sens)
bin_info_S2a <- dT_sens %>%
  filter(!is.na(fdelStress), !is.na(fbTime), !is.na(iTime)) %>%
  group_by(fbTime) %>%
  summarise(
    xmin    = min(iTime, na.rm = TRUE),
    xmax    = max(iTime, na.rm = TRUE),
    xmid    = mean(c(xmin, xmax)),
    width_h = xmax - xmin,
    .groups = "drop"
  )

# Identify significant interaction bin and reference bin for annotation
sig_bins_S2a <- bin_info_S2a %>%
  filter(fbTime == "b-12") %>%
  mutate(star = "***", alpha_val = 0.12)

ref_bin_S2a <- bin_info_S2a %>%
  filter(fbTime == "b-2") %>%
  mutate(star = "*")

# Add bin midpoints to data for accurate error-bar placement
dT_plot_S2a <- dT_sens %>%
  left_join(bin_info_S2a, by = "fbTime")

# Panel a: Variable pre-EMA window (2h bins, dT_sens)
plot_S2a <-
  ggplot(data = subset(dT_plot_S2a, !is.na(fdelStress)),
         aes(y = estIVglu, x = iTime)) +
  geom_rect(data = sig_bins_S2a,
            aes(xmin = xmin, xmax = xmax, ymin = 103, ymax = 109, alpha = alpha_val),
            inherit.aes = FALSE, fill = "grey70", color = "grey90") +
  geom_rect(data = ref_bin_S2a,
            aes(xmin = xmin, xmax = xmax, ymin = 103, ymax = 109),
            inherit.aes = FALSE, alpha = 0.8, fill = "#e8f4f8") +
  scale_alpha_identity() +
  stat_smooth(aes(fill = fdelStress), linewidth = 1.2,
              method = "gam", alpha = 0.8, color = "grey30") +
  stat_summary(
    data = subset(dT_plot_S2a, !is.na(fdelStress) & !is.na(xmid)),
    aes(x = xmid, y = estIVglu, group = fdelStress, color = fdelStress),
    fun.data = "mean_cl_boot", geom = "linerange",
    linewidth = 1.5, alpha = 0.8, inherit.aes = FALSE
  ) +
  scale_fill_manual(guide  = guide_legend(title = "Stress bin"),
                    breaks = c("Stress onset", "No change", "Stress relief"),
                    values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_color_manual(guide  = guide_legend(title = "Stress bin"),
                     breaks = c("Stress onset", "No change", "Stress relief"),
                     values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_x_continuous(breaks = seq(-12, 0, by = 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-12, 0), ylim = c(103, 109)) +
  geom_text(data = sig_bins_S2a,
            aes(x = xmid, y = 108.5, label = star),
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  geom_text(data = ref_bin_S2a,
            aes(x = xmid, y = 104.5, label = "*"),
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  annotate("segment", x = -2, y = 108.6, xend = -3.5, yend = 108.6,
           arrow = arrow(type = "closed", length = unit(0.03, "npc"))) +
  annotate("text", x = -3.8, y = 109, label = "Interaction", size = 4) +
  theme(
    text              = element_text(face = "bold",  size = 12.0),
    axis.text         = element_text(face = "plain", size = 12.0),
    axis.text.x       = element_text(size = 12.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 13.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 12.0, hjust = 0.5),
    legend.key.height = unit(1, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.text       = element_text(face = "bold", size = 10.0),
    legend.position   = "none"
  ) +
  ggtitle("Glucose differences before stress changes",
          "with variable pre-EMA window") +
  xlab("Time before EMA [h]") +
  ylab("Glucose level [mg/dl]")

# Panel b: Fixed 12h pre-EMA window (1.5h bins, dT) - mirrors Figure 4a
plot_S2b <-
  ggplot(aes(y = estIVglu, x = iTime),
         data = subset(dT, !is.na(fdelStress))) +
  annotate("rect", xmin = -7.5, xmax = -6, ymin = 103, ymax = 109,
           alpha = 0.06, fill = "grey70", color = "grey90") +
  annotate("rect", xmin = -1.5, xmax = 0,  ymin = 103, ymax = 109,
           alpha = 0.8,  fill = "#e8f4f8") +
  stat_smooth(aes(fill = fdelStress), linewidth = 1.2,
              method = "gam", alpha = 0.8, color = "grey30") +
  stat_summary_bin(aes(group = fdelStress, color = fdelStress),
                   fun.data = "mean_cl_boot",
                   geom = "linerange", linewidth = 1.5, alpha = 0.8, bins = 12) +
  scale_fill_manual(guide  = guide_legend(title = "Stress bin"),
                    breaks = c("Stress onset", "No change", "Stress relief"),
                    values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_color_manual(guide  = guide_legend(title = "Stress bin"),
                     breaks = c("Stress onset", "No change", "Stress relief"),
                     values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_x_continuous(breaks = seq(-12, 0, by = 1.5), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-12, 0), ylim = c(103, 109)) +
  annotate("text",    x = -6.75, y = 108.5, label = "*",           size = 6) +
  annotate("segment", x = -1.6,  y = 108.6, xend = -3.6, yend = 108.6,
           arrow = arrow(type = "closed", length = unit(0.03, "npc"))) +
  annotate("text",    x = -3.6,  y = 109,   label = "Interaction", size = 4) +
  annotate("text",    x = -0.75, y = 104.5, label = "*",           size = 6) +
  theme(
    text              = element_text(face = "bold",  size = 12.0),
    axis.text         = element_text(face = "plain", size = 12.0),
    axis.text.x       = element_text(size = 12.0),
    strip.text.x      = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
    plot.title        = element_text(size = 13.0, hjust = 0.5),
    plot.subtitle     = element_text(size = 12.0, hjust = 0.5),
    legend.key.height = unit(1, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.text       = element_text(face = "bold", size = 10.0),
    axis.title.y      = element_blank()
  ) +
  ggtitle("Glucose differences before stress changes",
          "with fixed 12h pre-EMA window") +
  xlab("Time before EMA [h]") +
  ylab("Glucose level [mg/dl]")

plot_S2 <- plot_grid(plot_S2a, plot_S2b,
                     labels = c("a", "b"), label_size = 12,
                     ncol = 2, rel_widths = c(1.2, 1.5), align = "h")

ggsave("figures/Figure_S2.png",
       plot = plot_S2, height = 5.5, width = 9.5, units = "in", dpi = 600, bg = "white")


# =============================================================================
# END OF SCRIPT
# =============================================================================
