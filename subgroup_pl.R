# Load Library
library(tidyverse)
library(afex)
library(emmeans)
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

# Define the column to averaged
amp_cols <- grep("^PeakLat_E", names(data), value = TRUE)

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
#new <- bind_rows(baseline, avg_ab, avg_cd)
new <- bind_rows(baseline, avg_hr, avg_rpe)

# Reshape the data: pivot electrode columns into long format
df <- new %>%
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

# Show the ANOVA result
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

roi_order <- c("fronto_central", "central", "central_parietal", "parietal")

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


mod2_df <- mod2_df %>%
  mutate(ROI = factor(ROI, levels = roi_order))

mod2 <- aov_car(
  MeanLatency ~ Group * Session * Stimuli * ROI + 
    Error(PID/(Session * Stimuli * ROI)),
  data = mod2_df,
  factorize = FALSE
)

mod2_summary <- mod2$anova_table


## Tukey-adjusted pairwise comparison
# Estimated marginal means
session_emm <- emmeans(mod2, ~ Session | ROI)

# Summarize EMMs and compute approximate SD
session_emm_summary <- summary(session_emm) %>%
  mutate(SD = SE * sqrt(df))

# Tukey-adjusted pairwise comparisons
session_pairs <- pairs(session_emm, adjust = "tukey") %>%
  summary(infer = TRUE)

# Compute effect size
session_pairs <- session_pairs %>%
  mutate(
    cohen_d = t.ratio / sqrt(df)
  )

# Split contrast into Session1 and Session2
session_pairs <- session_pairs %>%
  mutate(
    Session1 = sub(" - .*", "", contrast),
    Session2 = sub(".* - ", "", contrast)
  )

# Join with EMMs to get means and SDs
session_pairs_with_means <- session_pairs %>%
  left_join(session_emm_summary, by = c("Session1" = "Session", "ROI")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(session_emm_summary, by = c("Session2" = "Session", "ROI")) %>%
  rename(Mean2 = emmean, SD2 = SD)

session_pairs_with_means <- session_pairs_with_means %>%
  mutate(
    Session1 = factor(Session1, levels = session_order),
    Session2 = factor(Session2, levels = session_order)
  ) %>%
  mutate(
    s_low  = pmin(as.numeric(Session1), as.numeric(Session2)),
    s_high = pmax(as.numeric(Session1), as.numeric(Session2)),
    
    Session_low  = session_order[s_low],
    Session_high = session_order[s_high],
    
    contrast_clean = paste0(Session_low, " - ", Session_high)
  )

# Preserve consistent ordering
desired_order <- unique(session_pairs_with_means$contrast_clean)

session_pairs_with_means <- session_pairs_with_means %>%
  mutate(
    contrast_clean = factor(contrast_clean, levels = desired_order)
  )

# Plot
ggplot(session_pairs_with_means,
       aes(x = estimate, 
           y = contrast,
           color = contrast)) +
  
  geom_point(size = 3) +
  
  geom_errorbarh(
    aes(xmin = lower.CL.x, xmax = upper.CL.x),
    height = 0.2
  ) +
  
  facet_wrap(~ ROI) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "BL - HR" = "#D37f7f",
    "BL - RPE" = "#105446",
    "HR - RPE"  = "#874537"
    
   # "BL - HR" = "#A155B9",
  #  "BL - RPE" = "#165BAA",
   # "HR - RPE"  = "#C2185B"
  )) +
  
  scale_y_discrete(
    limits = rev(unique(session_pairs_with_means$contrast))
  ) +
  
  labs(
    x = "Estimate",
    y = "Session Comparison"
  ) +
  
  theme_classic()

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

