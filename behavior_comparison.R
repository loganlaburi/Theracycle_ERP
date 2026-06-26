
# Load Library
library(tidyverse)
library(afex)
library(emmeans)
library(mediation)
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

data$Acc <- 100*data$Acc

summary_stats <- data %>%
  group_by(Group, Session, Stimuli) %>%
  summarise(
    mean = mean(RT, na.rm = TRUE),
    sd = sd(RT, na.rm = TRUE),
  ) %>%
  mutate(
    Formatted = paste0(round(mean, 2), " ± ", round(sd, 2))
  )



mod <- aov_car(
  Acc ~ Group * Session * Stimuli + 
    Error(PID / (Session * Stimuli)),
  data = data,
  factorize = FALSE
)

mod_summary <- mod$anova_table

# ----------------- Demographic information
data <- read.csv('demo.csv')
colnames(data)

var.test(SBP_Resting ~ PD, data = data)

t.test(SBP_Resting ~ PD, data = data, alternative = "two.sided",
       var.equal = FALSE)

# --------------- Intervention Session
data <- read.csv('Intervention.csv')

# Convert relevant columns to factors
data <- data %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Group = as.factor(Group),
    Session = as.factor(Session),
  )

mod <- aov_car(
  oxygen_peak ~ Group * Session  + 
    Error(Participant_ID / Session),
  data = data,
  factorize = FALSE
)

mod_summary <- mod$anova_table


# Post-hoc for Group
group_summary <- emmeans(mod, pairwise ~ Group)
group_table <- summary(group_summary$contrasts)

# Post-hoc for Session
session_summary <- emmeans(mod, pairwise ~ Session)
session_table <- summary(session_summary$contrasts)


# Get estimated marginal means for each Group × Session
emm <- emmeans(mod, ~ Group * Session)
# Pairwise comparisons of Session within each Group
session_contrasts_by_group <- contrast(emm, method = "pairwise", by = "Group")
table <- as.data.frame(summary(session_contrasts_by_group))


session_stats <- data %>%
  group_by(Session) %>%
  summarise(
    mean = mean(oxygen_peak, na.rm = TRUE),
    sd = sd(oxygen_peak, na.rm = TRUE),
  ) %>%
  mutate(
    Formatted = paste0(round(mean, 2), " ± ", round(sd, 2))
  )


summary_stats <- data %>%
  group_by(Group, Session) %>%
  summarise(
    mean = mean(oxygen_peak, na.rm = TRUE),
    sd = sd(oxygen_peak, na.rm = TRUE),
  ) %>%
  ungroup()


#Filter data for Session A
df <- data %>% filter(Session == "A")

# Run independent t-test
var.test(HR_peak ~ Group, data = df)

t.test(HR_peak ~ Group, data = df,
       alternative = "two.sided",
       var.equal = TRUE)


