# Higher glucose levels buffer against everyday stress load

## Authors
- Madeleine Kördel
- Anne Kühnel
- Kristin Kaduk
- Alica Guzman
- Melanie Henes
- Melina Grahlow
- Birgit Derntl
- Nils B. Kroemer

## Overview

This repository contains the R workflow used to compute the linear mixed-effects models for the Glucose ↔︎ Stress project, as well as the corresponding figures of the manuscript. The analysis examines the relationship between continuous glucose monitoring data and ecological momentary assessments (EMA) of stress and glucose.

---

## Repository Structure

### Data Files

- `TUE009_Glucose_VAS_Blood.xlsx` - Main dataframe with EMA and glucose data
- `TUE009_CGM_VAS_Stress.csv` - Dataframe with glucose data time-locked to EMA assessments

### Scripts

- `Stress_Glu.R` - Main analysis script containing:
  - Data preprocessing and variable creation
  - Linear mixed-effects models
  - Statistical tests (correlations, K-S tests)
  - Publication-ready figure generation

### Output
- `figures/` - Directory containing generated publication figures:
  - `Figure_2_Stress_Traces.png` - Ridgeline plot of stress trajectories
  - `Figure_3_Glucose_Stress.png` - Glucose-stress relationship panels
  - `Figure_4_Temporal_Dynamics.png` - Temporal dynamics analysis
  - `Figure_5_BMI_Stress.png` - BMI category comparisons

---

## Data Preprocessing

### Variable Centering

The script creates standardized variables to properly partition variance:

**Grand-mean centering** (continuous predictors):

- `cHappy`, `cSad` - Centered mood items
- `cBMI`, `cAge` - Centered demographic variables
- `cMoodState` - Centered overall mood state (scaled 0-100)
- `cTyG` - Centered triglyceride-glucose index

**Binary variable centering**:

- `cSex` - Sex coded as -0.5 and 0.5 (centered around 0)

**Group-centering** (within-person):

- `gcGlu` - Glucose deviations from each participant's mean
- `M_Glu` - Between-person glucose (participant means)

### Derived Variables

**Temporal derivatives** (run-to-run changes):

- `delStress` - Change in stress between consecutive assessments
- `delMood` - Change in mood state between assessments  
- `delTime` - Time difference between assessments (converted to minutes)

**Categorical variables**:

- `fBMI` - BMI categories: "Normal weight" (<25), "Overweight" (25-30), "Obese" (>30)
- `fdelStress` - Stress change categories: "Stress relief", "No change", "Stress onset"
- `fbTime` - Time bins before EMA labeled as factors (b-12 to b-2)

**Aggregate measures**:

- `M_MoodState`, `M_Stress` - Participant-level averages
- `R_MoodState` - Mood state rank across participants

**Time variables**:

- `bTime` - Binned time periods (6 bins) before EMA
- `iTime` - Index before EMA in hours
- `cIndEMA` - Centered index before EMA
- `EMA_h` - Hour of day for EMA assessment
- `cEMA_h` - Centered hour of day

**Residualized stress**:

- `lmeResStress` - Stress residuals after controlling for mood (Happy, Sad) for visualization

---

## Correlational Analyses

The script computes Pearson correlations to examine:

1. **Stress-mood associations**:
   - Stress with overall mood state
   - Mean stress with mean mood state

2. **Independence of stress changes**:
   - Change in stress (delStress) with mean mood state
   - Change in stress with change in mood (delMood)

3. **Metabolic associations** (subject-level):
   - Glucose-stress slopes (EB estimates) with BMI
   - Glucose-stress slopes with HOMA-IR (insulin resistance)
   - Glucose-stress slopes with TyG index (tryglyceride-glucose)

---

## Linear Mixed-Effects Models

### Model 1: Everyday Stress and Mood Associations

**fm_1a** - Overall mood state predicting stress:
```r
Stress ~ cMoodState + cBMI + cAge + cSex + (1 + cMoodState|ID)
```

**fm_1b** - Happy/sad mood items predicting stress:
```r
Stress ~ cHappy + cSad + cBMI + cAge + cSex + (1 + cHappy + cSad|ID)
```

**fm_1c** - Average mood predicting stress changes:
```r
delStress ~ scale(M_MoodState) + delMood + cBMI + cAge + cSex + 
            (1 + scale(M_MoodState) + delMood|ID)
```

**fm_1d** - Concurrent mood predicting stress changes:
```r
delStress ~ cMoodState + delMood + cBMI + cAge + cSex + 
            (1 + cMoodState + delMood|ID)
```

### Model 2: Glucose-Stress Relationships

**fm_2res** - Residualizing stress for visualization:
```r
Stress ~ (cHappy + cSad) + cBMI + cAge + (1 + cHappy + cSad|ID)
```
*Residuals stored as `lmeResStress` for Figure 3a*

**fm_2a** - Main glucose-stress model (Figure 3):
```r
Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + 
         (1 + cHappy + cSad + gcGlu|ID)
```
*Tests mood × glucose interaction, controls for between-person glucose*

**fm_2b** - Glucose-stress with sex interaction:
```r
Stress ~ (cHappy + cSad + cSex) * gcGlu + M_Glu + cBMI + cAge + cSex + 
         (1 + cHappy + cSad + gcGlu + cSex|ID)
```

### Model 3: Metabolic Moderation

**fm_HOMA** - Insulin resistance moderation:
```r
Stress ~ (cHappy + cSad) * gcGlu * log_HOMA + M_Glu + cBMI + cAge + cSex + 
         (1 + cHappy + cSad + gcGlu|ID)
```
*Three-way interaction: mood × glucose × HOMA-IR*

**fm_TyG** - Triglyceride-glucose index moderation:
```r
Stress ~ (cHappy + cSad) * gcGlu * cTyG + M_Glu + cBMI + cAge + cSex + 
         (1 + cHappy + cSad + gcGlu|ID)
```
*Three-way interaction: mood × glucose × TyG*

### Model 4: Temporal Dynamics (Figure 4)

**fm_3a** - Time-binned glucose predicting stress changes:
```r
estIVglu ~ fbTime * scale(StressDeriv) + (1 + fbTime + scale(StressDeriv)|id)
```
*Tests whether glucose at different time points before EMA predicts subsequent stress changes*

**fm_3b** - Adding time-of-day control:
```r
estIVglu ~ fbTime * scale(StressDeriv) + cEMA_h + 
           (1 + fbTime + scale(StressDeriv) + cEMA_h|id)
```
*Last time bin (0-2h before EMA) set as reference using `relevel()`*

**fm_3c** - Time-of-day effects on stress changes:
```r
EMA_h ~ scale(StressDeriv) + (1 + scale(StressDeriv)|id)
```
*Uses data from IndexBefore == -1 only*

---

## Distribution Comparisons

### Kolmogorov-Smirnov Tests
Tests whether stress change distributions differ across BMI categories:

1. Normal weight vs. Overweight
2. Normal weight vs. Obese  
3. Overweight vs. Obese

Parameters: `exact = TRUE, simulate.p.value = TRUE, B = 2000`

---

## Key Variables

### Glucose measures
- `Glu` - Raw glucose levels (mg/dL)
- `gcGlu` - Group-centered glucose (within-person deviations)
- `M_Glu` - Mean glucose per participant (between-person)
- `estIVglu` - Estimated glucose in time bins before EMA

### Mood/Stress measures
- `Stress` - Stress VAS (0-100 scale)
- `Happy`, `Sad` - Individual mood items
- `MoodState` - Composite mood state (happy - sad, scaled 0-100)
- `delStress` - Temporal derivative of stress
- `delMood` - Temporal derivative of mood
- `lmeResStress` - Residualized stress (controlling for Happy, Sad)

### Metabolic parameters
- `HOMA_IR` - Homeostatic Model Assessment for Insulin Resistance
- `log_HOMA` - Log-transformed HOMA-IR
- `TyG` - Triglyceride-glucose index
- `cTyG` - Centered TyG index
- `BMI` - Body Mass Index
- `fBMI` - BMI categories

### Temporal variables
- `run_ind` - Run/assessment index
- `timestamp` - Time of assessment
- `delTime` - Time between assessments (minutes)
- `IndexBefore` - Index of time points before EMA
- `iTime` - Time before EMA in hours
- `bTime` - Binned time periods (1-6)
- `fbTime` - Labeled time bins (b-12 through b-2)
- `EMA_h` - Hour of day (0-23)
- `cEMA_h` - Centered hour of day

### Demographics
- `ID` - Participant identifier
- `Age`, `cAge` - Age variables
- `Sex_male` - Binary sex variable
- `cSex` - Centered sex (-0.5, 0.5)

### Derived categorical variables
- `fdelStress` - Stress change categories
- `R_MoodState` - Mood state rank

---

## Visualizations

The script generates publication-ready figures using `ggplot2` and `cowplot`:

### Figure 2: Stress Trajectories (`plot2`)

**Type**: Ridgeline plot

- Individual stress traces across runs
- Y-axis: Participants ranked by mean mood state
- Color: Mean mood state
- **Output**: `Figure_2_Stress_Traces.png`

### Figure 3: Glucose-Stress Relationship (`plot3`)

**Panel a** (`plot3a`): Density plot with model predictions

- X-axis: Group-centered glucose
- Y-axis: Residualized stress
- Overlays: 2D density, model predictions with confidence intervals
- Uses `ggpredict()` to extract model predictions from `fm_2a`

**Panel b** (`plot3b`): Distribution of individual slopes

- Y-axis: Empirical Bayes estimates for glucose slopes
- Sina plot shows individual participant estimates
- Bar shows mean effect
- Uses `coef(fm_2a)$ID` to extract EB estimates

**Output**: `Figure_3_Glucose_Stress.png` 

### Figure 4: Temporal Dynamics (`plot4`)

**Panel a** (`plot4a`): Glucose trajectories before stress changes

- Stratified by stress change categories (onset, relief, no change)
- GAM smoothing with confidence bands
- Highlights three time windows with significance annotations
- Uses data from `dT` dataframe

**Panel b** (`plot4b`): Stress change distribution

- Density distribution
- Colored by stress change category

**Panel c** (`plot4c`): Early glucose predicting stress

- 2D density plot with individual regression lines
- Uses glucose 12-10h before EMA (`bTime == 1`)
- Robust linear regression (`method = 'rlm'`)

**Panel d** (`plot4d`): Time-of-day effects

- Relationship between stress changes and EMA timing
- GAM smoothing with confidence intervals
- Uses `IndexBefore == -1` data

**Output**: `Figure_4_Temporal_Dynamics.png` 

### Figure 5: BMI Categories (`plot5`)

**Type**: Slab-interval plot

- Distribution of stress changes by BMI category
- Shows full distribution with intervals
- Color-coded by BMI group
- **Output**: `Figure_5_BMI_Stress.png`

---

## Key Findings

1. **Stress-mood relationship**: Stress and mood state are robustly correlated, but run-to-run stress changes are independent of average mood state

2. **Glucose buffering effect**: Higher glucose levels were associated with lower self-reported stress

3. **Temporal dynamics**: Dynamic changes in glucose levels precipitated subsequent increases in stress levels, with lower glucose levels associated with subsequent stress increases and elevated glucose with stress relief

4. **BMI effects**: Stress change distributions differ by BMI category, with overweight and obese groups showing altered dynamics compared to normal weight

---

## Requirements

### R Packages

```r
# Data manipulation and visualization
library(tidyverse)
library(dplyr)
library(readxl)
library(foreign)

# Statistical modeling
library(lmerTest)      # Linear mixed-effects models
library(MASS)          # Robust regression
library(emmeans)       # Post-hoc comparisons
library(ggeffects)     # Model predictions
library(Hmisc)         # Statistical functions

# Visualization
library(ggplot2)
library(cowplot)       # Plot composition
library(viridis)       # Color palettes
library(ggridges)      # Ridgeline plots
library(ggforce)       # Extended ggplot2
library(ggdist)        # Distribution visualizations
library(ggside)        # Marginal plots
library(hexbin)        # Hexagonal binning
```


### Theme Settings
The script sets a global plotting theme:
```r
theme_set(theme_cowplot(font_size = 12))
```

---



## Citation

If you use this code or data, please cite the original publication.

---

## Contact

For questions or issues, please contact the corresponding authors or open an issue in this repository.

---

## Notes

- All glucose measurements are in mg/dL
- VAS scores are scaled 0-100
- Time indices are in hours relative to EMA assessments  
- Group centering (gcGlu) removes between-subject variance in glucose levels
- Missing data handled via `na.rm = TRUE` in calculations
