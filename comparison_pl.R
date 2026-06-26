
# Load Library
library(tidyverse)
library(afex)
options(scipen = "999")

# Read the data
data <- read.csv('pl.csv')

# Convert relevant columns to factors
data <- data %>%
  mutate(
    PID = as.factor(PID),
    Group = as.factor(Group),
    Session = as.factor(Session),
    Stimuli = as.factor(Stimuli)
  )

# Reshape the data: pivot electrode columns into long format
df <- data %>%
  pivot_longer(
    cols = starts_with("PeakLat_E"),  # all electrode columns
    names_to = "Electrode",
    values_to = "Latency"
  ) %>%
  mutate(
    Electrode = gsub("PeakLat_", "", Electrode)  # clean electrode names
  )

# --------------------------------
# Model 1 : Individual Electrode
# --------------------------------
mod1_electrode <- c("E6", "E11", "E55", "E62",
                    "E72", "E75", "E129")

mod1_df <- df %>% 
  filter(Electrode %in% mod1_electrode) 

mod1 <- aov_car(
  Latency ~ Group * Session * Stimuli * Electrode +
    Error(PID/(Session * Stimuli * Electrode)), 
  data = mod1_df, factorize = FALSE
)

mod1_summary <- mod1$anova_table

# -----------------------
# Model 2: Averaged ROIs
#------------------------
roi_lookup <- tibble(
  Electrode = c("E6", "E7", "E13", "E106", "E112", "E129",
                "E31", "E37", "E54", "E55", "E79", "E80", "E87", "E129",
                "E37", "E42", "E53", "E54", "E55", "E79", "E86", "E87", "E93",
                "E61", "E62", "E66", "E67", "E71", "E72", "E76", "E77", "E78", "E84"),
  ROI = c(rep("fronto_central", 6),
          rep("central", 8),
          rep("central_parietal", 9),
          rep("parietal", 10))
)

mod2_df_check <- df %>%
  left_join(roi_lookup, by = "Electrode", relationship = "many-to-many") %>%
  filter(!is.na(ROI)) %>%
  count(PID, Group, Session, Stimuli, ROI)
        
mod2_df <- df %>%
  left_join(roi_lookup, by = "Electrode", relationship = "many-to-many") %>%
  filter(!is.na(ROI)) %>%
  group_by(PID, Group, Session, Stimuli, ROI) %>%
  mutate(MeanLatency = mean(Latency, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(PID, Group, Session, Stimuli, ROI, .keep_all = TRUE) %>%
  select(PID, Group, Session, Stimuli, ROI, MeanLatency)

mod2 <- aov_car(
  MeanLatency ~ Group * Session * Stimuli * ROI + 
    Error(PID/(Session * Stimuli * ROI)),
  data = mod2_df,
  factorize = FALSE
)

mod2_summary <- mod2$anova_table

# -----------------------------
# Model 3: Single Averaged ROI
# -----------------------------
pd_hotspot <- tibble(
  Electrode = c("E62", "E67", "E77", "E72", "E71", "E76"),
  ROI = "PD_hotspot"
)

mod3_df <- df %>%
  left_join(pd_hotspot, by = "Electrode", relationship = "many-to-many") %>%
  filter(!is.na(ROI)) %>%
  group_by(PID, Group, Session, Stimuli, ROI) %>%
  mutate(MeanLatency = mean(Latency, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(PID, Group, Session, Stimuli, ROI, .keep_all = TRUE) %>%
  select(PID, Group, Session, Stimuli, ROI, MeanLatency)

mod3 <- aov_car(
  MeanLatency ~ Group * Session * Stimuli + 
    Error(PID / (Session * Stimuli)),
  data = mod3_df,
  factorize = FALSE
)

mod3_summary <- mod3$anova_table
