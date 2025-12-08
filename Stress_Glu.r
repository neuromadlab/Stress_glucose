# =============================================================================
# Higher glucose levels buffer against everyday stress load
# Analysis Script 
# Creating figures
# =============================================================================
# Description: This script performs linear mixed-effects modeling and creates 
#              figures to examine the relationship between continuous glucose 
#              monitoring and ecological momentary assessments of stress & mood
# =============================================================================

# Load required packages -------------------------------------------------------
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

# Set plotting theme
theme_set(theme_cowplot(font_size = 12))

# =============================================================================
# 1. DATA IMPORT
# =============================================================================

# Main dataframe with EMA and glucose data
d <- read_excel("TUE009_Glucose_VAS_Blood.xlsx")

# Dataframe with glucose data time-locked to EMA
dT <- read.csv("TUE009_CGM_VAS_Stress.csv")

# =============================================================================
# 2. DATA PREPROCESSING
# =============================================================================

# Grand mean centering for continuous variables -------------------------------
d$cHappy <- d$Happy - mean(d$Happy, na.rm = TRUE)
d$cSad <- d$Sad - mean(d$Sad, na.rm = TRUE)
d$cBMI <- d$BMI - mean(d$BMI, na.rm = TRUE)
d$cAge <- d$Age - mean(d$Age, na.rm = TRUE)

# Center binary variable around 0
d$cSex <- d$Sex_male - 0.5

# Group-center glucose (within-person centering)
d$M_Glu <- d$Glu - d$gcGlu

# Create BMI categories --------------------------------------------------------
d$fBMI <- cut(d$BMI, 
              breaks = c(-100, 25, 30, 'Inf'), 
              labels = c("Normal weight", "Overweight", "Obese"))

# Calculate temporal derivatives -----------------------------------------------
d <- d %>%
  arrange(ID, run_ind) %>%
  group_by(ID) %>%
  mutate(
    delStress = c(NA, diff(Stress)),  # Change in stress
    delMood = c(NA, diff(MoodState)), # Change in mood
    delTime = c(NA, diff(timestamp))  # Time difference
  ) %>%
  ungroup()

# Convert time difference to minutes
d$delTime <- d$delTime / 60000

# Classify stress changes
d$fdelStress <- cut(d$delStress, 
                    breaks = c(-100, -5, 5, 'Inf'), 
                    labels = c("Stress relief", "No change", "Stress onset"))

# Aggregate mood rank data -----------------------------------------------------
dRank <- d %>%  
  group_by(ID) %>%
  summarise(
    M_MoodState = mean(MoodState),
    M_Stress = mean(Stress)
  ) %>%
  mutate(R_MoodState = rank(M_MoodState)) # adds new variable that ranks participants based on their average mood state

# Merge rank data
d <- merge(d, dRank, by.x = "ID", by.y = "ID")

# Scale mood variables to 0-100
d$MoodState <- d$MoodState * 100
d$cMoodState <- d$MoodState - mean(d$MoodState, na.rm = TRUE)
d$M_MoodState <- d$M_MoodState * 100


# Process dT dataframe ---------------------------------------------------------
# Classify stress derivatives
dT$fdelStress <- cut(dT$StressDeriv, 
                     breaks = c(-100,-5,5,'Inf'), 
                     labels = c("Stress relief", "No change", "Stress onset"))

# Center index before EMA
dT$cIndEMA <- dT$IndexBefore - mean(dT$IndexBefore, na.rm = TRUE)

# Bin time data
dT <- dT %>% 
  mutate(bTime = ntile(IndexBefore, n = 6))

# Create factor labels for time bins
dT$fbTime <- factor(dT$bTime, 
                    labels = c("b-12", "b-10", "b-8", "b-6", "b-4", "b-2"))
# Convert index before to hours
dT$iTime <- dT$IndexBefore / 12

bin <- hexbin(d$Stress, d$Happy, xbins=40) #  create hexbin object for plotting
plot(bin, main="" , legend=F ) #  plot hexbins

# Extract time of day
dT$EMA_tT <- dmy_hms(dT$time)
dT$EMA_h <- hour(dT$EMA_tT)
dT$cEMA_h <- dT$EMA_h - mean(dT$EMA_h, na.rm = TRUE)

# =============================================================================
# 3. CORRELATIONAL ANALYSES
# =============================================================================

# Everyday stress is associated with mood, but stress spikes occur independently (Figure 2) -------------------
#### correlations ####
cor.test(d$Stress,d$MoodState) # correlation stress and mood state
cor.test(d$delStress,d$M_MoodState) # correlation change in stress and mean mood state
cor.test(d$delStress,d$delMood) # correlation change in stress and change in mood state
cor.test(d$M_Stress,d$M_MoodState) # correlation mean stress and mean mood state

# =============================================================================
# 4. LINEAR MIXED-EFFECTS MODELS
# =============================================================================

# Everyday stress is associated with mood, but stress spikes occur independently (Figure 2) -------------------
# Model 1a
fm_1a <- lmer(Stress ~ cMoodState + cBMI + cAge + cSex + (1 + cMoodState|ID), d, REML = T)
summary(fm_1a)

# Model 1b
fm_1b <- lmer(Stress ~ cHappy + cSad + cBMI + cAge + cSex + (1 + cHappy + cSad|ID), d, REML = T)
summary(fm_1b)

# Model 1c
fm_1c <- lmer(delStress ~ scale(M_MoodState) + delMood + cBMI + cAge + cSex + (1 + scale(M_MoodState) + delMood|ID), d, REML = T)
summary(fm_1c)

# Model 1d
fm_1d <- lmer(delStress ~ cMoodState + delMood + cBMI + cAge + cSex + (1 + cMoodState + delMood|ID), d, REML = T)
summary(fm_1d)



# Higher glucose levels are associated with lower stress ratings (Figure 3) -------------------------
# Model 2res: Collect stress residuals from mood predictors for Figure 3
fm_2res <- lmer(Stress ~ (cHappy + cSad) + cBMI + cAge + (1 + cHappy + cSad|ID), data = d)
summary(fm_2res)
d$lmeResStress <- residuals(fm_2res)

# Tests whether glucose levels (group-centered and mean) interact with mood (Happy, Sad) to predict Stress (Figure 3)
# Model 2a
fm_2a <- lmer(Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + (1 + cHappy + cSad + gcGlu|ID), d)
summary(fm_2a)

# Test whether glucose levels (group-centered and mean) interact with mood and sex to predict Stress
fm_2b <- lmer(Stress ~ (cHappy + cSad + cSex) * gcGlu + M_Glu + cBMI + cAge + cSex + (1 + cHappy + cSad + gcGlu + cSex|ID), d)
summary(fm_2b)


#Test whether metabolic parameters moderate the glucose stress effects
# Model for HOMA-IR
fm_HOMA <- lmer(Stress ~ (cHappy + cSad) * gcGlu * log_HOMA + M_Glu + cBMI + cAge + cSex + (1 + cHappy + cSad + gcGlu|ID), data = d)
summary(fm_HOMA)

# Model for triglyceride-glucose index (TyG)
fm_TyG <- lmer(Stress ~ (cHappy + cSad) * gcGlu * cTyG + M_Glu + cBMI + cAge + cSex + (1 + cHappy + cSad + gcGlu|ID), data = d)
summary(fm_TyG)

# Extract empirical Bayes estimates for Figure 3b
d_fm2a <- coef(fm_2a)$ID

# Are metabolic parameters associated with glucose slopes? 
# Create subject-level dataset
d_subj <- d %>%
  group_by(ID) %>%
  summarise(
    cBMI = first(cBMI),
    log_HOMA = first(log_HOMA),
    cTyG = first(cTyG)
  )

# Add glucose slopes (EB estimates)
d_subj$gcGlu_slope <- d_fm2a$gcGlu

# Correlations with metabolic parameters
cor.test(d_subj$gcGlu_slope, d_subj$cBMI) # BMI

cor.test(d_subj$gcGlu_slope, d_subj$log_HOMA) # HOMA-IR

cor.test(d_subj$gcGlu_slope, d_subj$cTyG) # TyG

# Differences in glucose levels precipitate changes in stress reported during EMA (Figure 4) -------------------------
# Set last time bin as reference 
dT <- within(dT, fbTime <- relevel(fbTime, ref = 6)) # 0-2h as reference category

# Model 3a: Test whether temporal derivatives of stress (change in stress) predict glucose levels at different binned time points before EMA assessments
fm_3a <- lmer(estIVglu ~ fbTime * scale(StressDeriv) + (1 + fbTime + scale(StressDeriv)|id), data = dT)
summary(fm_3a)

# Model 3b: Including time of day as covariate
fm_3b <- lmer(estIVglu ~ fbTime * scale(StressDeriv) + cEMA_h + (1 + fbTime + scale(StressDeriv) + cEMA_h|id), data = dT)
summary(fm_3b)

# Model 3c: Time of day effects (Figure 4d)
fm_3c <- lmer(EMA_h ~ scale(StressDeriv) + (1 + scale(StressDeriv)|id), data = subset(dT, IndexBefore == -1))
summary(fm_3c)

# =============================================================================
# 5. DISTRIBUTION COMPARISONS (Kolmogorov-Smirnov Tests)
# =============================================================================
# NW = Normal weight, OW = Overweight, OB = Obese

# Subset data by BMI category
d_NW <- filter(d, fBMI == "Normal weight")
d_OW <- filter(d, fBMI == "Overweight")
d_OB <- filter(d, fBMI == "Obese")

# BMI category comparisons (Figure 5) -----------------------------------------

# Normal weight vs. Overweight
ks.test(d_NW$delStress, d_OW$delStress, 
              exact = TRUE, simulate.p.value = TRUE, B = 2000)
# Normal weight vs. Obese
ks.test(d_NW$delStress, d_OB$delStress, 
              exact = TRUE, simulate.p.value = TRUE, B = 2000)
# Overweight vs. Obese
ks.test(d_OW$delStress, d_OB$delStress, 
              exact = TRUE, simulate.p.value = TRUE, B = 2000)


# =============================================================================
# 6. Create FIGURES
# =============================================================================
# =============================================================================
# Figure Generation Script
# =============================================================================
# Note: Run main analysis script (4. LINEAR MIXED-EFFECTS MODELS) first
# =============================================================================

# FIGURE 2: Individual traces of stress ratings (Ridgeline Plot) -----------------------
plot2 <- 
  ggplot(aes(x = run_ind, y = R_MoodState), data = d) +
  geom_ridgeline(aes(height = Stress/15, color = M_MoodState, group = ID), 
                 linewidth = 0.7, alpha = 0.01) +
  scale_color_gradient2("Mood", 
                        low = "steelblue", 
                        mid = "grey80", 
                        high = "tomato", 
                        midpoint = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(20, 80, 20)) +
  theme(text = element_text(face = 'bold', size = 20.0),
        axis.text = element_text(face = 'plain', size = 18.0),
        axis.text.x = element_text(size = 18.0), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
        plot.title = element_text(size = 20.0, hjust = 0.5),
        plot.subtitle = element_text(size = 18.0, hjust = 0.5),
        legend.key.height = unit(2, 'cm')) +
  ggtitle('Individual traces of stress ratings',
          'EMA, rescaled for display') +
  xlab(label = 'Runs') +
  ylab(label = 'Mood rank [VAS happy - sad]')

# save plot2
ggsave("figures/Figure_2_Stress_Traces.png", 
       plot = plot2, 
       height = 9, width = 8, units = "in", 
       dpi = 300, bg = "white")

# FIGURE 3: Glucose & Stress plot (density Plot) -----------------------

# Generate predicted values from mixed model
dPfm2a <- ggpredict(fm_2a, terms = "gcGlu")

# FIGURE 3a: Density plot with model predictions
plot3a <-
  ggplot(dPfm2a, aes(x, predicted)) +
  stat_density_2d(aes(x = gcGlu, y = lmeResStress + 35.72, 
                      fill = after_stat(ndensity)), 
                  geom = "raster",
                  contour = FALSE, n = 500, data = d) +
  geom_line(linewidth = 2, color = 'slateblue') +
  scale_fill_viridis_c(guide = guide_legend(title = "Density"), 
                       option = "F") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              fill = 'grey80', alpha = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(30, 45), xlim = c(-30, 40)) +
  theme(text = element_text(face = 'bold', size = 14.0),
        axis.text = element_text(face = 'plain', size = 14.0),
        axis.text.x = element_text(size = 12.0), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  ylab(label = 'Stress [residualized for mood]') +
  xlab(label = 'Glucose levels [group centered]')

# FIGURE 3b: Distribution of EB estimates
plot3b <-
  ggplot(d_fm2a, aes(x = 1, y = gcGlu)) +
  geom_bar(linewidth = 1, color = 'grey80', fill = 'slateblue', 
           stat = "summary") +
  geom_sina(size = 2.5, adjust = 0.5) +
  geom_hline(yintercept = 0, color = 'grey20', size = 1) +
  theme(text = element_text(face = 'bold', size = 14.0),
        axis.text = element_text(face = 'plain', size = 14.0),
        axis.text.x = element_blank(), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
        axis.ticks.x = element_blank()) +
  ylab(label = 'EB estimate') +
  xlab(label = 'Glu slope')


# Combine panels
plot3 <- plot_grid(plot3a, plot3b, 
                  labels = c("a", "b"), 
                  label_size = 12, 
                  ncol = 2, 
                  rel_widths = c(1.40, 0.60), 
                  align = "h")
# save plot3
ggsave("figures/Figure_3_Glucose_Stress.png", 
       plot = plot3, 
       height = 4.5, width = 6.5, units = "in", 
       dpi = 300, bg = "white")


# FIGURE 4: Temporal Dynamics --------------------------------

# FIGURE 4a: Differences in glucose before changes in stress
plot4a <- 
  ggplot(aes(y = estIVglu, x = iTime), 
         data = subset(dT, !is.na(fdelStress))) +
  annotate("rect", xmin = -12, xmax = -10, ymin = 103, ymax = 109,
           alpha = 0.1, color = "grey90") +
  annotate("rect", xmin = -8, xmax = -6, ymin = 103, ymax = 109,
           alpha = 0.1, color = "grey90") +
  annotate("rect", xmin = -2, xmax = 0, ymin = 103, ymax = 109,
           alpha = 0.8, fill = "#e8f4f8") +
  stat_smooth(aes(fill = fdelStress), size = 1.2, 
              method = 'gam', alpha = 0.8, color = 'grey30') +
  stat_summary_bin(aes(group = fdelStress, color = fdelStress), 
                   fun.data = "mean_cl_boot", 
                   geom = "linerange", size = 1.5, alpha = 0.8, bins = 12) + 
  scale_fill_manual(guide = guide_legend(title = "Stress bin"),
                    breaks = c('Stress onset', 'No change', 'Stress relief'),
                    values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_color_manual(guide = guide_legend(title = "Stress bin"),
                     breaks = c('Stress onset', 'No change', 'Stress relief'),
                     values = c("#f5c396", "#faed7c", "#96f5c3")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-12, 0), ylim = c(103, 109)) +
  annotate("text", x = -11, y = 108.5, label = "***", size = 6) +
  annotate("text", x = -7, y = 108.5, label = "*", size = 6) +
  annotate("text", x = -1, y = 104.5, label = "*", size = 6) +
  annotate("segment", x = -2, y = 108.6, xend = -4, yend = 108.6,
           arrow = arrow(type = "closed", length = unit(0.03, "npc"))) +
  annotate("text", x = -3.8, y = 109, label = "Interaction", size = 4) +
  theme(text = element_text(face = 'bold', size = 12.0),
        axis.text = element_text(face = 'plain', size = 12.0),
        axis.text.x = element_text(size = 12.0), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
        plot.title = element_text(size = 13.0, hjust = 0.5),
        plot.subtitle = element_text(size = 12.0, hjust = 0.5),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text = element_text(face = 'bold', size = 10.0)) +
  ggtitle('Differences in glucose',
          'before changes in stress (deriv)') +
  xlab(label = 'Time before EMA [h]') +
  ylab(label = 'Glucose level [mg/dl]')

# FIGURE 4b: Density distribution of stress changes
plot4b <- 
  ggplot(aes(y = delStress, fill = fdelStress), data = d) +
  stat_slabinterval(aes(fill = fdelStress), 
                    alpha = 1, color = 'grey20', adjust = 2) +
  scale_fill_manual(guide = guide_legend(title = "Stress cat"),
                    values = c("#96f5c3", "#faed7c", "#f5c396")) +
  coord_cartesian(ylim = c(-100, 100), xlim = c(0, 0.75)) +
  theme(legend.position = 'none', 
        text = element_text(face = 'bold', size = 12.0),
        axis.text = element_text(face = 'plain', size = 12.0),
        axis.text.x = element_text(size = 10.0), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
        plot.title = element_text(size = 13.0, hjust = 0.5),
        plot.subtitle = element_text(size = 12.0, hjust = 0.5),
        legend.key.height = unit(1, 'cm')) +
  ggtitle('Density') +
  xlab(label = 'Density') +
  ylab(label = 'Change in stress ratings')

# FIGURE 4c: Early differences in glucose percipitate changes in stress
plot4c <- 
  ggplot(aes(y = StressDeriv, x = estIVglu), 
         data = subset(dT, bTime == 1)) +
  stat_density_2d(geom = "raster", aes(fill = after_stat(ndensity)),
                  contour = FALSE, n = 500) +
  stat_smooth(size = 0.5, aes(group = id), 
              method = 'rlm', alpha = 0.01, color = 'grey70') +
  stat_smooth(size = 2, method = 'rlm', alpha = 0.05) +
  scale_fill_viridis_c(guide = guide_legend(title = "Density"), 
                       option = "F") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-25, 25), xlim = c(75, 135)) +
  theme(legend.position = 'none', 
        text = element_text(face = 'bold', size = 12.0),
        axis.text = element_text(face = 'plain', size = 12.0),
        axis.text.x = element_text(size = 10.0), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
        plot.title = element_text(size = 13.0, hjust = 0.5),
        plot.subtitle = element_text(size = 12.0, hjust = 0.5),
        legend.key.height = unit(1, 'cm')) +
  ggtitle('Early differences in glucose', 
          'precipitate changes in stress') +
  ylab(label = 'Changes in stress [deriv]') +
  xlab(label = 'Glucose levels [12 to 10h before EMA]')

# FIGURE 4d: Time of day effects
plot4d <- 
  ggplot(aes(x = StressDeriv, y = EMA_h), 
         data = subset(dT, IndexBefore == -1)) +
  stat_smooth(color = 'cornflowerblue', fill = '#90D5FF', size = 2, 
              method = 'gam', alpha = 0.5) +
  stat_summary_bin(color = 'cornflowerblue', fun.data = "mean_cl_boot", 
                   geom = "linerange", size = 1.5, alpha = 0.8, bins = 10) +
  scale_x_continuous(expand = c(0.01, 5)) +
  coord_cartesian(ylim = c(10, 20), xlim = c(-100, 100)) +
  theme(legend.position = 'none', 
        text = element_text(face = 'bold', size = 12.0),
        axis.text = element_text(face = 'plain', size = 12.0),
        axis.text.x = element_text(size = 10.0), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
        plot.title = element_text(size = 13.0, hjust = 0.5),
        plot.subtitle = element_text(size = 12.0, hjust = 0.5),
        legend.key.height = unit(1, 'cm')) +
  ggtitle('Time of day') +
  ylab(label = 'Time [h]') +
  xlab(label = 'Changes in stress')

# Combine all panels
plot4 <- plot_grid(plot4a, plot4b, plot4c, plot4d, 
                    labels = c("a", "b", "c", "d"), 
                    label_size = 12, 
                    ncol = 2, 
                    rel_widths = c(1.32, 0.68), 
                    align = "h")
# save plot4
ggsave("figures/Figure_4_Temporal_Dynamics.png", 
       plot = plot4, 
       height = 7.6, width = 6.8, units = "in", 
       dpi = 300, bg = "white")

# FIGURE 5: BMI and Stress Dynamics --------------------------------
plot5 <- 
  ggplot(aes(y = delStress, x = fBMI), data = d) +
  stat_slabinterval(aes(group = fBMI, fill = fBMI), alpha = 0.9) +
  scale_fill_manual(guide = guide_legend(title = "BMI cat"),
                    values = c("#CFE0C3", "#9EC1A3", "#70A9A1")) +
  theme(text = element_text(face = 'bold', size = 14.0),
        axis.text = element_text(face = 'plain', size = 14.0),
        axis.text.x = element_text(size = 12.0), 
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm")),
        plot.title = element_text(size = 16.0, hjust = 0.5),
        plot.subtitle = element_text(size = 14.0, hjust = 0.5),
        legend.key.height = unit(1.5, 'cm')) +
  ggtitle('Altered distribution of changes in stress',
          'Temporal derivative') +
  xlab(label = 'BMI categories') +
  ylab(label = 'Delta stress [VAS]')

# save plot5
ggsave("figures/Figure_5_BMI_Stress.png", 
       plot = plot5, 
       height = 6, width = 7, units = "in", 
       dpi = 600, bg = "white")

# =============================================================================
# END OF SCRIPT
# =============================================================================