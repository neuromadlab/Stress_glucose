# Higher glucose levels are associated with lower everyday stress load

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

This repository contains the R workflow used to compute the linear mixed-effects models for the Glucose ↔ Stress project, as well as the corresponding figures of the manuscript. The analysis examines the association between continuous glucose monitoring (CGM) data and ecological momentary assessments (EMA) of stress and mood.

---

## Repository Structure

### Data Files

- `TUE009_Glucose_VAS_Blood.xlsx` — Main dataframe: EMA and CGM data merged at the observation level
- `TUE009_GCM_VAS_Stress_12h.csv` — CGM data time-locked to EMA assessments (12h pre-EMA window, 1.5h bins)
- `TUE009_Influenca_data_runs_clean_Glu_LR.csv` — EMA dataframe used for plotting only
- `TUE009_CGM_VAS_Stress.csv` — Sensitivity analysis: CGM data with alternative (variable) pre-EMA binning
- `TUE009_SSES_questionnaire.xlsx` — Stress-eating questionnaire data (Stress and Snacking Effects Scale, SSES)

### Scripts

- `TUE009_analysis_clean.R` — Main analysis script containing:
  - Data preprocessing and variable creation
  - Sample description and demographic tests
  - Correlational analyses
  - Linear mixed-effects models (main and sensitivity)
  - Age-related analyses
  - Distribution comparisons (Kolmogorov–Smirnov tests)
  - Publication-ready figure generation

### Output

- `figures/` — Directory containing generated publication figures:
  - `Figure_2_Stress_Traces_Heatmap.png` — Ridgeline plot of individual stress traces and heatmap
  - `Figure_3_Glucose_Stress.png` — Glucose–stress relationship panels
  - `Figure_4.png` — Temporal dynamics of glucose before stress changes
  - `Figure_5_BMI_Stress.png` — BMI category comparisons
  - `Figure_S1_meal_timing.png` — Stress and glucose across meal times (z-scores)
  - `Figure_S2.png` — Sensitivity: pre-EMA window comparison (two binning approaches)

---

## Data Preprocessing

### Exclusion Criteria

Participant ID 87 was excluded from the main dataframe `d` for having fewer than 20 EMA runs. IDs 9, 11, 13, 25, 27, 73, and 87 were excluded from the CGM–EMA time-locked dataframes (`dT`, `dT_sens`) for the same reason.

### Variable Centering

The script creates standardized variables to properly partition variance:

**Grand-mean centering** (continuous predictors):
- `cHappy`, `cSad` — Centered mood items
- `cBMI`, `cAge` — Centered demographic variables
- `cMoodState` — Centered composite mood state (scaled 0–100)
- `cMetState` — Centered metabolic state (hunger–satiety); used in sensitivity analyses
- `cTyG` — Centered triglyceride–glucose index
- `cSSES` — Centered SSES stress-eating score

**Effect-coded binary variables** (centered at 0 by subtracting 0.5):
- `cSex` — Sex coded as −0.5 and +0.5
- `cRecentEating` — Recent eating (last meal within ~2h) coded as −0.5 and +0.5

**Within-person (group-mean) centering**:
- `gcGlu` — Glucose deviations from each participant's mean
- `M_Glu` — Between-person glucose (participant means); computed as `Glu - gcGlu`

### Derived Variables

**Temporal derivatives** (run-to-run changes):
- `delStress` — Change in stress between consecutive assessments
- `delMood` — Change in mood state between assessments
- `delTime` — Time elapsed between assessments (converted from milliseconds to minutes)

**Categorical variables**:
- `fBMI` — BMI categories: "Normal weight" (<25), "Overweight" (25–30), "Obese" (>30)
- `fdelStress` — Stress change categories: "Stress relief" (<−5), "No change" (−5 to 5), "Stress onset" (>5)
- `fbTime` — Quantile-based time bins before EMA (8 bins of ~1.5h in `dT`; 6 bins of ~2h in `dT_sens`)
- `fbTime_col` — Collapsed version of `fbTime` in `dT`: bins b-8 and b-7 merged into "b-87" to ensure sufficient data density

**Aggregate measures**:
- `M_MoodState`, `M_Stress` — Participant-level averages
- `R_MoodState` — Mood state rank across participants

**Time variables** (in `dT` and `dT_sens`):
- `iTime` — Pre-EMA CGM index converted to hours
- `cIndEMA` — Grand-mean centered pre-EMA index
- `EMA_h` — Hour of day for EMA assessment
- `cEMA_h` — Grand-mean centered hour of day

**Residualized stress**:
- `lmeResStress` — Stress residuals after controlling for Happy and Sad; used in Figure 3a

**Z-scored variables** (for Figure S1):
- `zGlu` — Z-scored raw glucose levels
- `zStress` — Z-scored stress ratings

---

## Sample Description & Demographic Tests

The script computes descriptive statistics (mean, SD, min, max) for Age, BMI, and number of EMA runs at the participant level, and reports counts by sex. Independent-samples t-tests are used to test demographic differences between male and female participants.

---

## Correlational Analyses

The script computes Pearson correlations to examine:

1. **Stress–mood associations**:
   - Momentary stress with momentary mood state
   - Person-mean stress with person-mean mood state

2. **Independence of stress changes**:
   - Change in stress (`delStress`) with person-mean mood state
   - Change in stress with concurrent change in mood (`delMood`)

3. **Metabolic associations** (subject-level):
   - Glucose–stress slopes (EB estimates from `fm_2a`) with BMI, HOMA-IR, and TyG index
   - Person-mean residual stress (`lmeResStress`) with SSES scores

4. **Age-related associations** (subject-level):
   - Age with mean glucose, mean stress, stress variability, glucose variability, HOMA-IR, and TyG index

5. **Metabolic state**:
   - Raw glucose with metabolic state (hunger–satiety) to assess collinearity

---

## Linear Mixed-Effects Models

All models were fit using `lmerTest::lmer()`. REML estimation is used for Models 1a–1d; ML estimation (default) is used for all other models to allow likelihood-ratio testing.

### 6.1 Stress–Mood Association (Figure 2)

**fm_1a** — Composite mood state predicting momentary stress:
```r
Stress ~ cMoodState + cBMI + cAge + cSex + (1 + cMoodState | ID)
```

**fm_1b** — Happy and Sad items separately to decompose composite mood:
```r
Stress ~ cHappy + cSad + cBMI + cAge + cSex + (1 + cHappy + cSad | ID)
```

**fm_1c** — Person-mean mood and mood change predicting stress derivative:
```r
delStress ~ scale(M_MoodState) + delMood + cBMI + cAge + cSex +
            (1 + scale(M_MoodState) + delMood | ID)
```

**fm_1d** — As fm_1c but using centered momentary mood instead of person-mean:
```r
delStress ~ cMoodState + delMood + cBMI + cAge + cSex +
            (1 + cMoodState + delMood | ID)
```

### 6.2 Glucose–Stress Association (Figure 3)

**fm_2res** — Residualise stress for mood (residuals used in Figure 3a):
```r
Stress ~ (cHappy + cSad) + cBMI + cAge + (1 + cHappy + cSad | ID)
```

**fm_2a** — Main glucose–stress model:
```r
Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex +
         (1 + cHappy + cSad + gcGlu | ID)
```
*Tests mood × within-person glucose interaction; controls for between-person glucose (M_Glu)*

**fm_2b** — Add sex × glucose interaction:
```r
Stress ~ (cHappy + cSad + cSex) * gcGlu + M_Glu + cBMI + cAge + cSex +
         (1 + cHappy + cSad + gcGlu | ID)
```

### 6.3 Temporal Glucose–Stress Dynamics (Figure 4)

**fm_3a** — Glucose × time bin interaction predicting stress derivative (12h window):
```r
estIVglu ~ fbTime_col * scale(StressDeriv) +
           (1 + fbTime_col + scale(StressDeriv) | id)
```
*Reference category: collapsed earliest bin "b-87" (set via `relevel()`)*

**fm_3b** — As fm_3a, additionally controlling for time of day:
```r
estIVglu ~ fbTime_col * scale(StressDeriv) + cEMA_h +
           (1 + fbTime_col + scale(StressDeriv) + cEMA_h | id)
```
*Uses `nloptwrap` optimizer with relaxed convergence tolerances*

**fm_3c** — Time of day associated with stress derivative (Figure 4d):
```r
EMA_h ~ scale(StressDeriv) + (1 + scale(StressDeriv) | id)
```
*Restricted to `IndexBefore == -1` (measurement immediately before EMA)*

**fm_4a** — Control analysis: effects independent of actual CGM–EMA time distance:
```r
estIVglu ~ cIndEMA * scale(StressDeriv) + (1 + cIndEMA | id) + (1 | Run)
```

### 6.4 Sensitivity Analyses

#### Metabolic state as alternative predictor

**fm_2c** — Add metabolic state alongside glucose:
```r
Stress ~ (cHappy + cSad) * gcGlu + M_Glu + cBMI + cAge + cSex + cMetState +
         (1 + cHappy + cSad + gcGlu + cMetState | ID)
```

**fm_2d** — Replace glucose with metabolic state as within-person predictor:
```r
Stress ~ (cHappy + cSad) * gcMetState + M_MetState + cBMI + cAge + cSex +
         (1 + cHappy + cSad + gcMetState | ID)
```
*AIC/BIC model comparison performed between fm_2a and fm_2d*

#### Recent eating

**fm_2a_eating** — Recent eating as additional covariate  
**fm_2b_eating** — Recent eating × glucose interaction  
**fm_2c_eating** — Recent eating added to metabolic state model  
**fm_2d_eating** — Recent eating × metabolic state interaction  
*Likelihood-ratio tests compare models with and without recent eating*

#### Stress-eating behavior (SSES)

Models 2a–2f_sses test whether SSES scores (as covariate or interaction term) explain the glucose–stress and metabolic state–stress associations. Run on the `d_sses` subsample (inner join with `dQ`).

#### Metabolic parameters as moderators

**fm_HOMA** — HOMA-IR as three-way moderator:
```r
Stress ~ (cHappy + cSad) * gcGlu * log_HOMA + M_Glu + cBMI + cAge + cSex +
         (1 + cHappy + cSad + gcGlu | ID)
```

**fm_TyG** — TyG index as three-way moderator:
```r
Stress ~ (cHappy + cSad) * gcGlu * cTyG + M_Glu + cBMI + cAge + cSex +
         (1 + cHappy + cSad + gcGlu | ID)
```

#### Alternative pre-EMA window (Figure S2)

Models fm_3a_sens, fm_3b_sens, fm_3c_sens replicate the temporal dynamics models using `dT_sens` (6 bins of ~2h over 12h) to verify robustness of results.

### 6.5 Age-Related Analyses

Models fm_1a_age through fm_1d_age test whether age moderates the stress–mood association. Models fm_2a_age and fm_2d_age test whether age moderates the glucose–stress and metabolic state–stress associations. Between-person models (`lm`) test whether age predicts stress variability and glucose variability.

---

## Distribution Comparisons

### Kolmogorov–Smirnov Tests

Tests whether stress-change distributions differ across BMI categories (Figure 5):

1. Normal weight vs. Overweight
2. Normal weight vs. Obese
3. Overweight vs. Obese

Parameters: `exact = TRUE, simulate.p.value = TRUE, B = 2000`

---

## Key Variables

### Glucose measures
| Variable | Description |
|---|---|
| `Glu` | Raw glucose levels (mg/dL) |
| `gcGlu` | Within-person (group-mean centered) glucose |
| `M_Glu` | Between-person glucose (participant mean) |
| `estIVglu` | Estimated glucose in time bins before EMA |
| `zGlu` | Z-scored glucose (for Figure S1) |

### Mood & Stress measures
| Variable | Description |
|---|---|
| `Stress` | Stress VAS (0–100) |
| `Happy`, `Sad` | Individual mood items |
| `MoodState` | Composite mood (Happy – Sad, scaled 0–100) |
| `delStress` | Temporal derivative of stress |
| `delMood` | Temporal derivative of mood |
| `lmeResStress` | Stress residualized for Happy and Sad |
| `zStress` | Z-scored stress (for Figure S1) |

### Metabolic parameters
| Variable | Description |
|---|---|
| `BMI` | Body Mass Index |
| `fBMI` | BMI categories (Normal weight / Overweight / Obese) |
| `log_HOMA` | Log-transformed HOMA-IR |
| `cTyG` | Centered triglyceride–glucose index |
| `MetState` | Metabolic state (hunger–satiety) |
| `gcMetState` | Within-person centered metabolic state |
| `L_Meal` | Time since last meal (1=<0.5h … 6=>3h) |
| `RecentEating` | Binary: last meal within ~2h (postprandial window) |

### Temporal variables
| Variable | Description |
|---|---|
| `run_ind` | EMA run/assessment index |
| `timestamp` | Time of assessment |
| `delTime` | Time between assessments (minutes) |
| `IndexBefore` | Pre-EMA CGM index |
| `iTime` | Time before EMA in hours |
| `bTime` | Quantile-based time bin (integer) |
| `fbTime` | Labelled time bins (factor) |
| `fbTime_col` | Collapsed time bins (b-8 and b-7 merged into b-87) |
| `EMA_h` | Hour of day (0–23) |
| `cEMA_h` | Grand-mean centered hour of day |

### Demographics
| Variable | Description |
|---|---|
| `ID` | Participant identifier |
| `Age`, `cAge` | Age (raw and grand-mean centered) |
| `Sex_male` | Binary sex variable (0 = female, 1 = male) |
| `cSex` | Effect-coded sex (−0.5 / +0.5) |

### SSES
| Variable | Description |
|---|---|
| `SSES` | Mean stress-eating score (renamed from `SSES_mean`) |
| `cSSES` | Grand-mean centered SSES |

---

## Visualizations

Figures are generated with `ggplot2` and assembled with `cowplot`. All figures are saved to `figures/` at 600 dpi.

### Figure 2: Stress Traces & Heatmap

**Panel a** (`plot2a`): Ridgeline plot of individual stress traces across runs, ranked by mean mood state and coloured by mean mood.  
**Panel b** (`plot2b`): Heatmap of EMA stress ratings per run, sorted by participant mood rank.  
**Output**: `Figure_2_Stress_Traces_Heatmap.png`

### Figure 3: Glucose–Stress Relationship

**Panel a** (`plot3a`): 2D density plot with GAM model prediction overlay. X-axis: group-centered glucose; Y-axis: residualized stress. Model predictions and confidence intervals extracted via `ggpredict(fm_2a)`.  
**Panel b** (`plot3b`): Sina plot of empirical Bayes (EB) glucose–stress slopes per participant, extracted via `coef(fm_2a)$ID`.  
**Output**: `Figure_3_Glucose_Stress.png`

### Figure 4: Temporal Dynamics

**Panel a** (`plot4a`): GAM-smoothed glucose trajectories 0–12h before EMA, stratified by stress-change category (onset / relief / no change). Significant interaction and main-effect bins are annotated.  
**Panel b** (`plot4b`): Slab-interval density of stress-change categories.  
**Panel c** (`plot4c`): 2D density with individual and group robust regression lines for the earliest CGM window (12–10h before EMA, `bTime == 1`).  
**Panel d** (`plot4d`): GAM of time of day as a function of stress derivative, restricted to `IndexBefore == -1`.  
**Output**: `Figure_4.png`

### Figure 5: BMI & Stress Dynamics

Slab-interval plot of stress-change distributions by BMI category.  
**Output**: `Figure_5_BMI_Stress.png`

### Figure S1: Stress and Glucose across Meal Times

Line plot of z-scored glucose and stress as a function of time since last meal (6 categories). Person-level means are computed first; grand means ± SEM are plotted. The postprandial window (≤2h) is shaded.  
**Output**: `Figure_S1_meal_timing.png`

### Figure S2: Sensitivity – Pre-EMA Window Comparison

**Panel a**: Temporal dynamics using the alternative binning approach (`dT_sens`, 6 bins of ~2h).  
**Panel b**: Temporal dynamics using the main binning approach (`dT`, 8 bins of ~1.5h); mirrors Figure 4a.  
**Output**: `Figure_S2.png`

---

## Key Findings

1. **Stress–mood relationship**: Stress and mood state are robustly correlated at the momentary level, but run-to-run stress changes are independent of a person's average mood state.

2. **Glucose buffering effect**: Higher within-person glucose levels were associated with lower self-reported stress, an effect that was consistent across the majority of participants (negative EB slopes).

3. **Temporal dynamics**: Lower glucose levels up to 12h before EMA were associated with subsequent increases in stress, while elevated glucose was associated with stress relief. This interaction was robust to control for time of day and replicated across two binning approaches (Figure S2).

4. **BMI effects**: Stress-change distributions differed by BMI category, with overweight and obese groups showing altered dynamics compared to the normal-weight group.

---

## Requirements

### R Packages

```r
# Data manipulation
library(tidyverse)
library(dplyr)
library(readxl)
library(foreign)
library(forcats)

# Statistical modeling
library(lmerTest)      # Linear mixed-effects models with p-values
library(MASS)          # Robust regression (rlm)
library(emmeans)       # Post-hoc comparisons
library(ggeffects)     # Model predictions for plotting
library(Hmisc)         # Bootstrap CIs (mean_cl_boot)

# Visualization
library(ggplot2)
library(cowplot)       # Plot composition and theming
library(viridis)       # Perceptually uniform color palettes
library(ggridges)      # Ridgeline plots
library(ggforce)       # Sina plots
library(ggdist)        # Distribution visualizations (stat_slabinterval)
library(ggside)        # Marginal plots
library(hexbin)        # Hexagonal binning
```

### Theme Settings

The script sets a global plotting theme at the top:
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
- VAS stress scores are scaled 0–100
- Mood items (Happy, Sad) are combined into a composite `MoodState = Happy – Sad`, then scaled to 0–100
- Time indices are in hours relative to EMA assessments
- Within-person centering (`gcGlu`) removes between-subject variance in glucose levels
- Missing data are handled via `na.rm = TRUE` in all summary calculations
- The `d_sses` dataframe is a subset of `d` restricted to participants with available SSES data; all sensitivity models using SSES are run on this subset
