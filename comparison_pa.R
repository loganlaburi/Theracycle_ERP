
# Load Library
library(tidyverse)
library(afex)
library(emmeans)
options(scipen = "999")

# Read the data
data <- read.csv('pa.csv')

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
    cols = starts_with("MeanPeakAmp_E"),  # all electrode columns
    names_to = "Electrode",
    values_to = "Amplitude"
  ) %>%
  mutate(
    Electrode = gsub("MeanPeakAmp_", "", Electrode)  # clean electrode names
  )

# --------------------------------
# Model 1: Individual Electrode
# --------------------------------
mod1_electrode <- c("E6", "E11", "E55", "E62","E72", "E75", "E129")

mod1_df <- df %>% 
  filter(Electrode %in% mod1_electrode) 

mod1 <- aov_car(
  Amplitude ~ Group * Session * Stimuli * Electrode +
    Error(PID/(Session * Stimuli * Electrode)), 
  data = mod1_df, factorize = FALSE
)

# Show the ANOVA result
mod1_summary <- mod1$anova_table

## Tukey-adjusted pairwise comparison
# Groups across Stimuli and Electrode
group_emm <- emmeans(mod1, ~ Group | Stimuli * Electrode)
group_emm_summary <- summary(group_emm) %>%
  mutate(SD = SE * sqrt(df))  # approximate SD

# Get pairwise comparisons
group_pairs <- pairs(group_emm, adjust = "tukey") %>%
  summary(infer = TRUE)

group_pairs <- group_pairs %>%
  mutate(
    cohen_d = t.ratio / sqrt(df),   
  )

# Split contrast into two groups
group_pairs <- group_pairs %>%
  mutate(
    Group1 = sub(" - .*", "", contrast),
    Group2 = sub(".* - ", "", contrast)
  )

# Join with EMMs to get means/SDs for both groups
group_pairs_with_means <- group_pairs %>%
  left_join(group_emm_summary, by = c("Group1" = "Group", "Stimuli", "Electrode")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(group_emm_summary, by = c("Group2" = "Group", "Stimuli", "Electrode")) %>%
  rename(Mean2 = emmean, SD2 = SD)

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(
    mean_label = paste0(
      Group1, ": ", round(Mean1, 2), "±", round(SD1, 2),
      "\n",
      Group2, ": ", round(Mean2, 2), "±", round(SD2, 2)
    )
  )

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(
    contrast_elec = paste0(contrast, " (", Electrode, ")")
  )

mod1_electrode <- c("E6", "E11", "E55", "E62","E72", "E75", "E129")

mod1_df <- df %>% 
  filter(Electrode %in% mod1_electrode) 

mod1 <- aov_car(
  Amplitude ~ Group * Session * Stimuli * Electrode +
    Error(PID/(Session * Stimuli * Electrode)), 
  data = mod1_df, factorize = FALSE
)

# Show the ANOVA result
mod1_summary <- mod1$anova_table

## Tukey-adjusted pairwise comparison
# Groups across Stimuli and Electrode
group_emm <- emmeans(mod1, ~ Group | Stimuli * Electrode)
group_emm_summary <- summary(group_emm) %>%
  mutate(SD = SE * sqrt(df))  # approximate SD

# Get pairwise comparisons
group_pairs <- pairs(group_emm, adjust = "tukey") %>%
  summary(infer = TRUE)

group_pairs <- group_pairs %>%
  mutate(
    cohen_d = t.ratio / sqrt(df),   
  )

# Split contrast into two groups
group_pairs <- group_pairs %>%
  mutate(
    Group1 = sub(" - .*", "", contrast),
    Group2 = sub(".* - ", "", contrast)
  )

# Join with EMMs to get means/SDs for both groups
group_pairs_with_means <- group_pairs %>%
  left_join(group_emm_summary, by = c("Group1" = "Group", "Stimuli", "Electrode")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(group_emm_summary, by = c("Group2" = "Group", "Stimuli", "Electrode")) %>%
  rename(Mean2 = emmean, SD2 = SD)

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(
    mean_label = paste0(
      Group1, ": ", round(Mean1, 2), "±", round(SD1, 2),
      "\n",
      Group2, ": ", round(Mean2, 2), "±", round(SD2, 2)
    )
  )

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(Electrode = factor(Electrode, levels = mod1_electrode)) %>%
  arrange(match(Electrode, mod1_electrode))

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(
    contrast_elec = paste0(contrast, " (", Electrode, ")")
  ) %>%
  mutate(
    contrast_elec = factor(
      contrast_elec,
      levels = rev(unique(contrast_elec))  
    )
  )

# Plot 
ggplot(group_pairs_with_means,
       aes(x = estimate, 
           y = contrast_elec,
           color = Stimuli,
           group = Stimuli)) +
  
  geom_point(size = 3,
             position = position_dodge(width = 0.6)) +
  
  geom_errorbarh(
    aes(
      xmin = lower.CL.x,
      xmax = upper.CL.x
    ),
    height = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  
  scale_color_manual(values = c(
    "Bin1" = "#009dc4",
    "Bin2" = "#FDA172",
    "Bin3" = "#7A9E7E"
  )) +
  
  theme_classic()

# Sessions across Stimuli and Electrode 
session_order <- c("BL", "A", "B", "C", "D")

session_emm <- emmeans(mod1, ~ Session | Stimuli * Electrode)
session_emm_summary <- summary(session_emm) %>%
  mutate(SD = SE * sqrt(df)) 

# Get pairwise comparisons
session_pairs <- pairs(session_emm, adjust = "tukey") %>%
  summary(infer = TRUE)

session_pairs <- session_pairs %>%
  mutate(
    cohen_d = t.ratio / sqrt(df),  
  )

# Split contrast into two groups
session_pairs <- session_pairs %>%
  mutate(
    Session1 = sub(" - .*", "", contrast),
    Session2 = sub(".* - ", "", contrast)
  )

# Join with EMMs to get means/SDs for both groups
session_pairs_with_means <- session_pairs %>%
  left_join(session_emm_summary, by = c("Session1" = "Session", "Stimuli", "Electrode")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(session_emm_summary, by = c("Session2" = "Session", "Stimuli", "Electrode")) %>%
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

# This is the order that is wanted in the graph
desired_order <- c(
  "BL - A", "BL - B", "BL - C", "BL - D",
  "A - B", "A - C", "A - D",
  "B - C", "B - D",
  "C - D"
)

session_pairs_with_means <- session_pairs_with_means %>%
  mutate(
    contrast_clean = factor(contrast_clean, levels = desired_order)
  )

# Plot
ggplot(session_pairs_with_means,
       aes(x = estimate, 
           y = contrast_clean,
           color = Stimuli,
           group = Stimuli)) +
  
  geom_point(size = 2.5,
             position = position_dodge(width = 0.6)) +
  
  geom_errorbarh(
    aes(
      xmin = lower.CL.x,
      xmax = upper.CL.x
    ),
    height = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  
  facet_wrap(~ Electrode) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Bin1" = "#009dc4",
    "Bin2" = "#FDA172",
    "Bin3" = "#7A9E7E"
  )) +
  
  scale_y_discrete(limits = rev(levels(session_pairs_with_means$contrast_clean))) +  
  
  labs(
    #title = "Post Hoc Effect Sizes by Electrode and Stimulus",
    x = "Estimate",
    y = "Session Comparison",
    color = "Stimulus"
  ) +
  
  theme_classic()


write.csv(session_pairs_with_means, "posthoc.csv")

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
  mutate(MeanAmplitude = mean(Amplitude, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(PID, Group, Session, Stimuli, ROI, .keep_all = TRUE) %>%
  select(PID, Group, Session, Stimuli, ROI, MeanAmplitude)

mod2 <- aov_car(
  MeanAmplitude ~ Group * Session * Stimuli * ROI + 
    Error(PID/(Session * Stimuli * ROI)),
  data = mod2_df,
  factorize = FALSE
)

mod2_summary <- mod2$anova_table

## Tukey-adjusted pairwise comparison
# Groups across Stimuli and Electrode
session_emm <- emmeans(mod2, ~ Session | Stimuli * ROI)
session_emm_summary <- summary(session_emm) %>%
  mutate(SD = SE * sqrt(df))  # approximate SD
# Get pairwise comparisons
session_pairs <- pairs(session_emm, adjust = "tukey") %>%
  summary(infer = TRUE)
# Split contrast into two groups
session_pairs <- session_pairs %>%
  mutate(
    Session1 = sub(" - .*", "", contrast),
    Session2 = sub(".* - ", "", contrast)
  )
session_order <- c("BL", "A", "B", "C", "D")
# Join with EMMs to get means/SDs for both groups
session_pairs_with_means <- session_pairs %>%
  left_join(session_emm_summary, by = c("Session1" = "Session", "Stimuli", "ROI")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(session_emm_summary, by = c("Session2" = "Session", "Stimuli", "ROI")) %>%
  rename(Mean2 = emmean, SD2 = SD)

session_pairs_with_means <- session_pairs_with_means %>%
  mutate(
    cohen_d = t.ratio / sqrt(df)
  )

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


session_pairs_with_means <- session_pairs_with_means %>%
  mutate(
    contrast_clean = factor(contrast_clean, levels = desired_order)
  )

# Plot
ggplot(session_pairs_with_means,
       aes(x = estimate, 
           y = contrast_clean,
           color = Stimuli)) +
  
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  
  geom_errorbarh(
    aes(
      xmin = lower.CL.x,
      xmax = upper.CL.x
    ),
    position = position_dodge(width = 0.6),
    height = 0.2
  ) +
  
  facet_wrap(~ ROI) +
  
  scale_color_manual(values = c(
    "Bin1" = "#009dc4",
    "Bin2" = "#FDA172",
    "Bin3" = "#7A9E7E"
  )) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_y_discrete(limits = rev) +
  
  labs(
    x = "Estimate",
    y = "Session Comparison",
    color = "Stimulus"
  ) +
  
  theme_classic()

# Sessions across Stimuli and Electrode 
group_emm <- emmeans(mod2, ~ Group | Session * Stimuli * ROI)
group_emm_summary <- summary(group_emm) %>%
  mutate(SD = SE * sqrt(df))  
group_pairs <- pairs(group_emm, adjust = "tukey") %>%
  summary(infer = TRUE)
group_pairs <- group_pairs %>%
  mutate(
    cohen_d = t.ratio / sqrt(df),   
  )
group_pairs <- group_pairs %>%
  mutate(
    Group1 = sub(" - .*", "", contrast),
    Group2 = sub(".* - ", "", contrast)
  )
group_pairs_with_means <- group_pairs %>%
  left_join(group_emm_summary, by = c("Group1" = "Group", "Session", "Stimuli", "ROI")) %>%
  rename(Mean1 = emmean, SD1 = SD) %>%
  left_join(group_emm_summary, by = c("Group2" = "Group", "Session", "Stimuli", "ROI")) %>%
  rename(Mean2 = emmean, SD2 = SD)

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(
    contrast_roi = paste0(contrast, " (", ROI, ")")
  )
group_pairs_with_means <- group_pairs_with_means %>%
  mutate(
    label = paste0(Session, " | ", Stimuli)
  )

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(label = Session)

group_pairs_with_means <- group_pairs_with_means %>%
  mutate(
    Session = factor(Session, levels = c("BL", "A", "B", "C", "D"))
  )

# Plot
ggplot(group_pairs_with_means,
       aes(x = estimate, 
           y = Session,
           color = Stimuli,
           group = Stimuli)) +
  
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  
  geom_errorbarh(
    aes(
      xmin = lower.CL.x,
      xmax = upper.CL.x
    ),
    height = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  
  facet_wrap(~ ROI) +
  
  scale_y_discrete(
    limits = rev(levels(group_pairs_with_means$Session))
  ) +   
  
  scale_color_manual(values = c(
    "Bin1" = "#009dc4",
    "Bin2" = "#FDA172",
    "Bin3" = "#7A9E7E"
  )) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  labs(
    title = "Group Differences (PD vs HOA) by ROI",
    x = "Cohen’s d",
    y = "Session",
    color = "Stimulus"
  ) +
  
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold")
  )

write.csv(session_pairs_with_means, "posthoc.csv")

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
  mutate(MeanAmplitude = mean(Amplitude, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(PID, Group, Session, Stimuli, ROI, .keep_all = TRUE) %>%
  select(PID, Group, Session, Stimuli, ROI, MeanAmplitude)

mod3 <- aov_car(
  MeanAmplitude ~ Group * Session * Stimuli + 
    Error(PID / (Session * Stimuli)),
  data = mod3_df,
  factorize = FALSE
)

mod3_summary <- mod3$anova_table
