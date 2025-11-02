# ==================================================================================================================================================================
#==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  CANNASTREET - Analysis of Manuscript "Short-term effects of cannabis legalisation in Germany on driving under the influence of cannabis"
# SUPPLEMENT: Sensitivity Analysis - Difference-in-Differences Approach for the whole sample (not limited to at least monthly cannabis users)
# CODE AUTHOR:    Anna Schranz
# DATE STARTED:   281025

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================


# 0) ESSENTIALS
# ______________________________________________________________________________________________________________________

# clean workspace
rm(list=ls())

packages <- c("tidyverse", "car", "haven", "kableExtra", "boot", "emmeans", "sf", "mice", "ggthemes")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# current date:
DATE <- format(Sys.Date(), "%Y%m%d")

#write and export data? 
dataexport <- TRUE #TRUE = export data, FALSE = no export

folder_path_plots <- "Output/figures/" # folder for plots
folder_path_tables <- "Output/tables/" # folder for tables

# create folders if they do not exist
if (!dir.exists(folder_path_plots)) {
  dir.create(folder_path_plots, recursive = TRUE)
}

if (!dir.exists(folder_path_tables)) {
  dir.create(folder_path_tables, recursive = TRUE)
}

set.seed(1234)  # Seed setzen für Reproduzierbarkeit
n_boot <- 2000

## FUNCTIONS 
# ______________________________________________________________________________________________________________________
calculate_prevalence <- function(data, merkmal, weights, n_boot, conf_level = 0.95) {
  
  # Gewichtete Prävalenz berechnen
  weighted_prev <- sum(data[[merkmal]] * data[[weights]], na.rm = TRUE) / sum(data[[weights]], na.rm = TRUE)
  
  # Ungewichtete Prävalenz berechnen
  unweighted_prev <- mean(data[[merkmal]], na.rm = TRUE)
  
  # Ungewichtete und gewichtete Fallzahlen berechnen
  unweighted_n <- sum(!is.na(data[[merkmal]]))
  weighted_n <- sum(data[[weights]], na.rm = TRUE)
  
  # Bootstrap-Funktion für gewichtete Prävalenz
  boot_function_weighted <- function(data, indices) {
    sampled_data <- data[indices, ]
    sum(sampled_data[[merkmal]] * sampled_data[[weights]], na.rm = TRUE) / sum(sampled_data[[weights]], na.rm = TRUE)
  }
  
  # Bootstrap-Funktion für ungewichtete Prävalenz
  boot_function_unweighted <- function(data, indices) {
    sampled_data <- data[indices, ]
    mean(sampled_data[[merkmal]], na.rm = TRUE)
  }
  
  # Bootstrapping durchführen
  boot_results_weighted <- boot(data, statistic = boot_function_weighted, R = n_boot)
  boot_results_unweighted <- boot(data, statistic = boot_function_unweighted, R = n_boot)
  
  # Konfidenzintervalle berechnen
  ci_weighted <- boot.ci(boot_results_weighted, type = "perc", conf = conf_level)
  ci_unweighted <- boot.ci(boot_results_unweighted, type = "perc", conf = conf_level)
  
  return(tibble(
    weighted_prop = weighted_prev,
    weighted_ci_lower = ci_weighted$perc[4],
    weighted_ci_upper = ci_weighted$perc[5],
    unweighted_prop = unweighted_prev,
    unweighted_ci_lower = ci_unweighted$perc[4],
    unweighted_ci_upper = ci_unweighted$perc[5],
    unweighted_n = unweighted_n,
    weighted_n = weighted_n
  ))
}

logisticPseudoR2s <- function(LogModel) {
  dev <- LogModel$deviance
  nullDev <- LogModel$null.deviance
  modelN <- length(LogModel$fitted.values)
  R.l <- 1 - dev / nullDev
  R.cs <- 1- exp ( -(nullDev - dev) / modelN)
  R.n <- R.cs / ( 1 - ( exp (-(nullDev / modelN))))
  cat("Pseudo R^2 for logistic regression\n")
  cat("Hosmer and Lemeshow R^2 ", round(R.l, 3), "\n")
  cat("Cox and Snell R^2 ", round(R.cs, 3), "\n")
  cat("Nagelkerke R^2 ", round(R.n, 3), "\n")
}

## COLOR DEFINITIONS
# ______________________________________________________________________________________________________________________
# used for nominal variables
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors_nominal <- cbPalette
# countries
colors_country <- c("AT" = "#EF3340", "DE" = "#000000")
# sex
colors_sex <- c("weiblich" = "#0072B2", "männlich" = "#D55E00")
# agegroups
blue_gradient <- colorRampPalette(c("#ADD8E6", "#000080"), space = "Lab")
colors_age <- blue_gradient(4)
# edugroups
green_gradient <- colorRampPalette(c("#90EE90", "#006400"), space = "Lab")
colors_edu <- green_gradient(3)
# DUIC frequencies
blue_colors <- blue_gradient(4)
gray_color <- c("#808080")
freq_colors <- c(gray_color, blue_colors) 

welle_labels <- c(
  "1" = expression(t[0]),  # tiefgestellte 0
  "2" = expression(t[1])   # tiefgestellte 1
)

country_labels <- c(
  "AT" = "Austria",
  "DE" = "Germany"
)

# DATA IMPORT
# ______________________________________________________________________________________________________________________
path_GER_1 <- file.path("Data/", "cleaned/", "full_GER_wave1_cleaned2.rds")
path_GER_2 <- file.path("Data/", "cleaned/", "full_GER_wave2_cleaned2.rds")
path_AT_1 <- file.path("Data/", "cleaned/", "full_AT_wave1_cleaned2.rds")
path_AT_2 <- file.path("Data/", "cleaned/", "full_AT_wave2_cleaned2.rds")

#read in data
data_GER_1 <- readRDS(path_GER_1)
data_GER_2 <- readRDS(path_GER_2)
data_AT_1 <- readRDS(path_AT_1)
data_AT_2 <- readRDS(path_AT_2)

# DATA PREPARATION
# ______________________________________________________________________________________________________________________
#select relevant variables
data_GER_1_sel <- data_GER_1 %>%
  select(-weights) %>%
  mutate(DUIC.MIXEDUSE_old = DUIC.MIXEDUSE, weights = weights_ohneWiederbefragt, Land = "DE") %>%
  dplyr::select(ID, sex, weights, wiederbefragt, edu_group, agegroup, agegroup_basicsample, age_in_years, AUDITC2,
                GSZB2, ZS, AS1, AS2, AS3, can_freq, can_use_12M, DRIVERLICENSE, DUIC, DUIC.FREQ.01M_wave1coding, DUIC.FREQ.12M_wave1coding,
                DUIC12m_full, DUIC12m_alone, DUIC12m_mixed, DUIC30d_full, DUIC30d_alone, DUIC30d_mixed, Welle, 
                MEDICALUSE.01, MEDICALUSE.02, DUIC.MIXEDUSE_old, Land, DUIC.INTOX, DGURBA, sport_freq, sport_12M, fish_12M, SOURCE.1, gisd_5, AGS,
                gambling_freq, alcohol_freq, tobacco_freq, DUIA, use_benzos, use_heroin)#, MILEAGE, starts_with("SPEEDING")) %>%
# mutate(MILEAGE = as_factor(MILEAGE), 
#        drives = case_when(MILEAGE == "Gar nicht" ~ FALSE,
#                           MILEAGE %in% c("Bis 5.000 km", "5.001 - 20.000 km", "Mehr als 20.000 km") ~ TRUE,
#                           TRUE ~ NA),
#        #speeding as factor
#        across(starts_with("SPEEDING"), as_factor)) %>%
# select(-MILEAGE)

data_GER_2_sel <- data_GER_2 %>%
  select(-weights) %>%
  mutate(weights = weights_ohneWiederbefragt, Land = "DE") %>% 
  select(ID, sex, weights, wiederbefragt, edu_group,  age_in_years, 
         agegroup, agegroup_basicsample, GSZB2, ZS, AS1, AS2, AS3, AUDITC2,
         can_freq, can_use_12M, DUIC.MIXEDUSE, DUIC.MED, DRIVERLICENSE, CWM.PRE, DUIC, DUIC.FREQ.01M_wave1coding, DUIC.FREQ.01M_exact, DUIC.FREQ.12M_wave1coding,
         DUIC.FREQ.01M_wave2coding, DUIC.FREQ.12M_wave2coding, DUIC.FREQ.12M_exact, DUIC.MED, DUIC12m_full, DUIC12m_alone, DUIC12m_mixed, 
         DUIC30d_full, DUIC30d_alone, DUIC30d_mixed, Welle, MEDICALUSE.01, MEDICALUSE.02, Land, DUIC.INTOX, DGURBA, sport_12M, sport_freq, fish_12M, SOURCE.1, gisd_5, AGS,
         gambling_freq, alcohol_freq, tobacco_freq, DUIA, use_benzos, use_heroin) #MILEAGE, starts_with("SPEEDING")) %>%
# mutate(MILEAGE = as_factor(MILEAGE),
#        drives = case_when(MILEAGE == "Gar nicht" ~ FALSE,
#                           MILEAGE %in% c("Bis 5.000 km", "5.001 - 20.000 km", "20.001 - 30.000 km", "mehr als 30.000 km") ~ TRUE,
#                           TRUE ~ NA)) %>%
# select(-MILEAGE)

data_AT_1_sel <- data_AT_1 %>%
  mutate(DUIC.MIXEDUSE_old = DUIC.MIXEDUSE, Land = "AT") %>%
  select(ID, sex, weights, edu_group, agegroup, agegroup_basicsample,  age_in_years,  AUDITC2,
         GSZB2, ZS, AS1, AS2, AS3, can_freq, can_use_12M, DRIVERLICENSE, DUIC, 
         DUIC.FREQ.01M_wave1coding, DUIC.FREQ.12M_wave1coding, DUIC12m_full, DUIC12m_alone, DUIC12m_mixed, 
         DUIC30d_full, DUIC30d_alone, Welle, MEDICALUSE.01, MEDICALUSE.02, DUIC.MIXEDUSE_old, Land, sport_12M, sport_freq, fish_12M, DGURBA, AGS,
         gambling_freq, alcohol_freq, tobacco_freq, DUIA) #starts_with("SPEEDING"))

data_AT_2_sel <- data_AT_2 %>%
  select(ID, sex, weights, edu_group, agegroup, agegroup_basicsample, GSZB2, ZS, AS1, AS2, AS3, can_freq, age_in_years, AUDITC2,
         can_use_12M, DRIVERLICENSE, DUIC.MIXEDUSE, DUIC.MED, DUIC, DUIC.FREQ.01M_exact, DUIC.FREQ.01M_wave2coding, DUIC.FREQ.12M_wave2coding,
         DUIC.FREQ.01M_wave1coding, DUIC.FREQ.01M_exact, DUIC.FREQ.12M_wave1coding,
         DUIC.FREQ.12M_exact, DUIC.MED, DUIC12m_full, DUIC12m_alone, DUIC12m_mixed, 
         DUIC30d_full, DUIC30d_alone, Welle, MEDICALUSE.01, MEDICALUSE.02, DUIC.INTOX, sport_12M, sport_freq, fish_12M, DGURBA, AGS,
         gambling_freq, alcohol_freq, tobacco_freq, DUIA) %>%
  mutate(Land = "AT")


#combine data
data <- bind_rows(data_GER_1_sel, data_GER_2_sel, data_AT_1_sel, data_AT_2_sel)

data$Land <- relevel(factor(data$Land), ref = "AT")
data$Welle <- factor(data$Welle, levels = c(1, 2))

# prepare variables
#----------------------------------------------------------
#randomly assign sex == "divers" to "male" or "female"
divers_indices <- which(data$sex == "divers")
n_divers <- length(divers_indices)
set.seed(42)
new_values <- sample(c("weiblich", "männlich"), n_divers, replace = TRUE)
data$sex[divers_indices] <- new_values
data$sex <- factor(data$sex, levels = c("männlich", "weiblich"))

#driver license as binary variable
data <- data %>%
  mutate(DRIVERLICENSE_bin = case_when(
    DRIVERLICENSE == "Ja" ~ 1,
    DRIVERLICENSE %in% c("Nein, noch nie gehabt", "Nein, aber früher gehabt") ~ 0,
    TRUE ~ NA_real_  # NA for any other cases
  ))

data <- data %>%
  mutate(tobacco_freq = haven::as_factor(tobacco_freq),
         alcohol_freq = haven::as_factor(alcohol_freq),
         gambling_freq = haven::as_factor(gambling_freq)) %>%
  mutate(tobacco_use_12M = ifelse(tobacco_freq %in% c("(fast) täglich", "mindestens einmal pro Woche", "mindestens einmal im Monat", "seltener als einmal im Monat"), 1, 0),
         alcohol_use_12M = ifelse(alcohol_freq %in% c("(fast) täglich", "mindestens einmal pro Woche", "mindestens einmal im Monat", "seltener als einmal im Monat"), 1, 0),
         gambling_12M = ifelse(gambling_freq %in% c("(fast) täglich", "mindestens einmal pro Woche", "mindestens einmal im Monat", "seltener als einmal im Monat"), 1, 0),
         use_benzos = as_factor(use_benzos),
         use_heroin = as_factor(use_heroin),
         use_sedatives = if_else(use_benzos == "ausgewählt" | use_heroin == "ausgewählt", 1, 0))


#     alcohol_freq_num = case_when(
#       alcohol_freq_chr == "(fast) täglich" ~ 365,
#       alcohol_freq_chr == "mindestens einmal pro Woche" ~ 52,
#       alcohol_freq_chr == "mindestens einmal im Monat" ~ 12,
#       alcohol_freq_chr == "seltener als einmal im Monat" ~ 6,
#       alcohol_freq_chr == "gar nicht" ~ 0,
#       TRUE ~ NA_real_
#     ),
#     AUDITC2_num = case_when(
#       AUDITC2 == 1 ~ 1.5,
#       AUDITC2 == 2 ~ 3.5,
#       AUDITC2 == 3 ~ 5.5,
#       AUDITC2 == 4 ~ 8,
#       AUDITC2 == 5 ~ 10,
#       TRUE ~ 0
#     ),
#     drinks_per_year = case_when(
#       alcohol_freq_chr == "gar nicht" ~ 0,  # gar nicht -> 0
#       alcohol_freq_chr != "gar nicht" & is.na(AUDITC2) ~ NA_real_,  # trinken, aber C2 fehlt -> NA
#       TRUE ~ alcohol_freq_num * AUDITC2_num
#     )
#   )
# ==================================================================================================================================================================
# Sample definition (all, but without repeated cases)
# ==================================================================================================================================================================
#subsample of non-repeated cases
df_GSZB2 <- subset(data, GSZB2 == 1)

# ==================================================================================================================================================================
# Preparation of DUIC variable
# ==================================================================================================================================================================

# identifier for exclusive medical users with a prescription
df_GSZB2 <- df_GSZB2 %>%
  mutate(
    EXCL_MEDICALUSER = ifelse(
      !is.na(MEDICALUSE.01) & !is.na(MEDICALUSE.02) &
        MEDICALUSE.02 == "Ja, mir wurde medizinisches Cannabis ärztlich verschrieben (bezahlt durch Krankenkasse oder als Selbstzahler)" &
        MEDICALUSE.01 == "Ausschließlich für medizinische Zwecke",
      1,
      0
    )
  )

table(df_GSZB2$EXCL_MEDICALUSER, useNA = "ifany")

#set DUIC = 0 for EXCL_MEDICALUSER
df_GSZB2 <- df_GSZB2 %>%
  mutate(DUIC12m_full = ifelse(EXCL_MEDICALUSER == 1, 0, DUIC12m_full)) 

#set DUIC = 0 for non-cannabis users (can_use_12M == 0)
df_GSZB2 <- df_GSZB2 %>%
  mutate(DUIC12m_full = ifelse(can_use_12M == 0, 0, DUIC12m_full))


table(df_GSZB2$DUIC12m_full, df_GSZB2$can_freq, df_GSZB2$Land, useNA = "ifany")
#for 322 cases in GER wave 1 DUIC12m_full is missing (less than monthly users)

# ==================================================================================================================================================================
#IMPUTATION OF MISSING DUIC CASES IN GERMANY WAVE 1 (LESS THAN MONTHLY USERS)
# ==================================================================================================================================================================

# subset German data for imputation
df_GSZB2_DE <- df_GSZB2 %>%
  filter(Land == "DE")

# define training data: Wave 2 (less-than-monthly users, with observed DUIC)
train_imp <- df_GSZB2_DE %>%
  filter(Welle == 2, can_freq == "seltener als einmal im Monat", !is.na(DUIC12m_full))
# Fit logistic regression model on training data
fit <- glm(DUIC12m_full ~ agegroup + sex + edu_group + DGURBA + gambling_freq + alcohol_freq + sport_freq,
           data = train_imp, family = binomial)

summary(fit)
# define target data: Wave 1 (less-than-monthly users, DUIC missing)
target_imp <- df_GSZB2_DE %>%
  filter(Welle == 1, can_freq == "seltener als einmal im Monat")


# stochastic imputation of missing DUIC values
set.seed(1234)
target_imp$DUIC12m_full_imputed <- rbinom(nrow(target_imp),  # number of cases to impute 
                                          1, # brnoulli draw: each case becomes 0 or 1
                                          predict(fit, newdata = target_imp, type = "response")) # use predicted probabilities from the logistic model

# predcit returns the predicted probability that each person in target_imp has DUIC = 1, based on the logistic regression


# merge imputed data back into German dataset
df_GSZB2_DE <- df_GSZB2_DE %>%
  left_join(target_imp %>% select(ID, DUIC12m_full_imputed), by = "ID") %>%
  mutate(DUIC12m_full_final = ifelse(is.na(DUIC12m_full), DUIC12m_full_imputed, DUIC12m_full))

# merge imputed data back into main dataset
df_GSZB2 <- df_GSZB2 %>%
  left_join(df_GSZB2_DE %>% select(ID, DUIC12m_full_final), by = "ID") %>%
  mutate(DUIC12m_full_imp = coalesce(DUIC12m_full_final, DUIC12m_full))

# check DUIC prev by wave and frequency
prop.table(table(df_GSZB2[df_GSZB2$can_freq == "seltener als einmal im Monat", ]$DUIC12m_full_imp,
                 df_GSZB2[df_GSZB2$can_freq == "seltener als einmal im Monat", ]$Welle,
                 df_GSZB2[df_GSZB2$can_freq == "seltener als einmal im Monat", ]$Land), margin = c(2,3))

# ==================================================================================================================================================================
# DiD Analysis
# ==================================================================================================================================================================

#fit DiD model
# DID analysis for DUIC full prevalence

df_GSZB2$Land <- relevel(factor(df_GSZB2$Land), ref = "AT")
df_GSZB2$Welle <- factor(df_GSZB2$Welle, levels = c(1, 2))  # Welle 1 als Referenz

DiD_unadjusted <- glm(DUIC12m_full_imp ~ Welle * Land, data = df_GSZB2, family = "binomial")
summary(DiD_unadjusted)
logisticPseudoR2s(DiD_unadjusted)
DiD_null <- glm(DUIC12m_full_imp ~ 1, data = df_GSZB2, family = "binomial")
anova(DiD_null, DiD_unadjusted, test = "Chisq")
DiD_unadj_coef <- exp(cbind(OR = coef(DiD_unadjusted), confint(DiD_unadjusted)))


# Adjusted DiD analysis 1 (only covariates)
DiD_adjusted <- glm(DUIC12m_full_imp ~ agegroup + sex + edu_group + DGURBA + Welle * Land, data = df_GSZB2, family = "binomial")
summary(DiD_adjusted)
DiD_coef_adjusted <- exp(cbind(OR = coef(DiD_adjusted), confint(DiD_adjusted)))
logisticPseudoR2s(DiD_adjusted)
df_DUIC12m_full_complete <- df_GSZB2[complete.cases(df_GSZB2[, c("agegroup", "sex", "edu_group", "DGURBA")]), ]
DiD_null_complete <- glm(DUIC12m_full_imp ~ 1, data = df_DUIC12m_full_complete, family = "binomial")
anova(DiD_null_complete, DiD_adjusted, test = "Chisq")
vif(DiD_adjusted)

#save as html table
if (dataexport) {
  library(kableExtra)
  kable(DiD_coef_adjusted, format = "html", digits = 2, caption = "DiD Analysis for DUIC (Sample 1, Sensitivity Analysis)") %>%
    kable_styling("striped", full_width = F) %>%
    save_kable(file = paste0(folder_path_tables, "DiD_DUIC_Sample1_", DATE, ".html"))
}

# ==================================================================================================================================================================
# Prevalence estimates by country and wave
aggregates_DUIC <- df_GSZB2 %>%
  group_by(Land, Welle) %>%
  summarise(
    results = list(calculate_prevalence(cur_data(), "DUIC12m_full_imp", "weights", n_boot = n_boot)),
    .groups = "drop"
  ) %>%
  unnest_wider(results)

