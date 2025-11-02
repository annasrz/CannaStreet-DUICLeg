# ==================================================================================================================================================================
#==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  CANNASTREET - Analysis of Manuscript "Short-term effects of cannabis legalisation in Germany on driving under the influence of cannabis"
# SUPPLEMENT: DUIC among among at least monthly medical-only cannabis users with a prescription
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
                GSZB2, ZS, AS1, AS2, AS3, can_freq, can_use_12M, DRIVERLICENSE, DUIC, 
                DUIC12m_full, Welle, MEDICALUSE.01, MEDICALUSE.02, Land, DGURBA)


data_GER_2_sel <- data_GER_2 %>%
  select(-weights) %>%
  mutate(weights = weights_ohneWiederbefragt, Land = "DE") %>% 
  select(ID, sex, weights, wiederbefragt, edu_group,  age_in_years, 
         agegroup, agegroup_basicsample, GSZB2, ZS, AS1, AS2, AS3, 
         can_freq, can_use_12M, DUIC.MIXEDUSE, DRIVERLICENSE, DUIC12m_full, Welle, MEDICALUSE.01, MEDICALUSE.02, Land, DGURBA) 

data_AT_1_sel <- data_AT_1 %>%
  mutate(DUIC.MIXEDUSE_old = DUIC.MIXEDUSE, Land = "AT") %>%
  select(ID, sex, weights, edu_group, agegroup, agegroup_basicsample,  age_in_years,  AUDITC2,
         GSZB2, ZS, AS1, AS2, AS3, can_freq, can_use_12M, DRIVERLICENSE, DUIC, 
         DUIC.FREQ.01M_wave1coding, DUIC.FREQ.12M_wave1coding, DUIC12m_full, DUIC12m_alone, DUIC12m_mixed, 
         DUIC30d_full, DUIC30d_alone, Welle, MEDICALUSE.01, MEDICALUSE.02, DUIC.MIXEDUSE_old, Land, sport_12M, sport_freq, fish_12M, DGURBA, AGS,
         gambling_freq, alcohol_freq, tobacco_freq, DUIA) 

data_AT_2_sel <- data_AT_2 %>%
  select(ID, sex, weights, edu_group, agegroup, agegroup_basicsample, GSZB2, ZS, AS1, AS2, AS3, can_freq, age_in_years,
         can_use_12M, DRIVERLICENSE, DUIC12m_full, Welle, MEDICALUSE.01, MEDICALUSE.02) %>%
  mutate(Land = "AT")


#combine data
data <- bind_rows(data_GER_1_sel, data_GER_2_sel, data_AT_1_sel, data_AT_2_sel)

data$Land <- relevel(factor(data$Land), ref = "AT")
data$Welle <- factor(data$Welle, levels = c(1, 2))


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



# ==================================================================================================================================================================
# Sample definition (at least monthly medical-only cannabis users with a prescription, non-repeated cases)
# ==================================================================================================================================================================
#subsample of non-repeated cases
monthly_users <- subset(data, GSZB2 == 1) %>% #non-repeated cases only
  filter(can_freq %in% c("mindestens einmal im Monat", "mindestens einmal pro Woche", "(fast) täglich")) %>% #at least monthly cannabis use)
  mutate(
    EXCL_MEDICALUSER = ifelse(
      !is.na(MEDICALUSE.01) & !is.na(MEDICALUSE.02) &
        MEDICALUSE.02 == "Ja, mir wurde medizinisches Cannabis ärztlich verschrieben (bezahlt durch Krankenkasse oder als Selbstzahler)" &
        MEDICALUSE.01 == "Ausschließlich für medizinische Zwecke",
      1,
      0
    )
  )

#check 
table(monthly_users$EXCL_MEDICALUSER, df_GSZB2$MEDICALUSE.02,df_GSZB2$MEDICALUSE.01, useNA = "ifany")

# get number of medical-only users with prescription
table(df_GSZB2$Welle,df_GSZB2$Land, useNA = "ifany")
table(monthly_users$EXCL_MEDICALUSER, df_GSZB2$Welle,df_GSZB2$Land, useNA = "ifany")

medical_users <- monthly_users %>%
  filter(EXCL_MEDICALUSER == 1)

# ==================================================================================================================================================================
# DUIC prevalence among at least monthly medical-only cannabis users with a prescription
# ==================================================================================================================================================================

#calculate prevalences
aggregates_DUIC <- medical_users %>%
  group_by(Land, Welle) %>%
  summarise(
    results = list(calculate_prevalence(cur_data(), "DUIC12m_full", "weights", n_boot = n_boot)),
    .groups = "drop"
  ) %>%
  unnest_wider(results)

glm_DUIC <- glm(DUIC12m_full ~ Welle + agegroup + sex + DRIVERLICENSE_bin + edu_group, data = medical_users %>% filter(Land == "DE"), family = binomial)
summary(glm_DUIC)