# Load Library
library(tidyverse)
library(afex)
library(emmeans)
options(scipen = "999")

# Read the data
data <- read.csv('data.csv')

# Convert relevant columns to factors
data <- data %>%
  mutate(
    PID = as.factor(PID),
    Group = as.factor(Group),
    Session = as.factor(Session),
    Stimuli = as.factor(Stimuli)
  )

# Define the column to averaged
amp_cols <- c("Acc", "RT")

# Average sessions A and B
avg_ab <- data %>%
  filter(Session %in% c("A", "B")) %>%
  group_by(PID, Group, Stimuli) %>%
  summarise(across(all_of(amp_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "AA")

# Average sessions C and D
avg_cd <- data %>%
  filter(Session %in% c("C", "D")) %>%
  group_by(PID,  Group, Stimuli) %>%
  summarise(across(all_of(amp_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "RB")

# Average sessions HR
avg_hr <- data %>%
  filter(Session %in% c("A", "C")) %>%
  group_by(PID, Group, Stimuli) %>%
  summarise(across(all_of(amp_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "HR")

# Average Sessions RPE
avg_rpe <- data %>%
  filter(Session %in% c("B", "D")) %>%
  group_by(PID, Group, Stimuli) %>%
  summarise(across(all_of(amp_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(Session = "RPE")

# Extract baseline session
baseline <- data %>%
  filter(Session == "BL") %>%
  select(PID, Group, Stimuli, all_of(amp_cols)) %>%
  mutate(Session = "BL")

# Combine the result
df <- bind_rows(baseline, avg_ab, avg_cd)
df <- bind_rows(baseline, avg_hr, avg_rpe)

# Anova result
mod <- df %>%
  group_by(PID, Group, Session, Stimuli) %>%
  summarise(resp = mean(RT, na.rm=TRUE), .groups="drop")

mod3 <- aov_car(
  resp ~ Group * Session * Stimuli + 
    Error(PID / (Session * Stimuli)),
  data = mod,
  factorize = FALSE
)

mod3_summary <- mod3$anova_table

