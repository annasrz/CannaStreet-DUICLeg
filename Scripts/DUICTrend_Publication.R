# ==================================================================================================================================================================
#==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  CANNASTREET - Analysis of Manuscript "Short-term effects of cannabis legalisation in Germany on driving under the influence of cannabis"
# CODE AUTHOR:    Anna Schranz
# DATE STARTED:   020625

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================


# 0) ESSENTIALS
# ______________________________________________________________________________________________________________________

# clean workspace
rm(list=ls())

packages <- c("tidyverse", "car", "haven", "kableExtra", "boot", "emmeans", "sf", "rnaturalearthdata", "rnaturalearth", "ggthemes")

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


#shape file for mapping municipalities to a map 
path_degurba_mapping <- "Data/Regio/DGURBA_RG_01M_2021/DGURBA_RG_01M_2021.shp" #https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/degree-urbanisation
degurba_mapping <- sf::st_read(path_degurba_mapping) %>%
  filter(CNTR_CODE %in% c("DE", "AT")) %>% #es fehlt die Gemeindekennziffer 70370 AT, verursacht missing in AT Datensatz
  rename(AGS = LAU_ID)


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

# pepare variable DUIA12M
data <- data %>%
  mutate(DUIA = as_factor(DUIA),
         DUIA12M = case_when(
           DUIA %in% c("Ja, innerhalb der letzten 30 Tage", "Ja, innerhalb der letzten 12 Monate") ~ 1,
           DUIA %in% c("Ja, länger als 12 Monate her", "Nein, nie") ~ 0,
           is.na(DUIA) ~ 0, #as NA means no alcohol use in past 12M
           TRUE ~ NA_real_
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
# #Sample description - Sample 1 (all, but without non-repeated cases)
# ==================================================================================================================================================================
#subsample of non-repeated cases
df_GSZB2 <- subset(data, GSZB2 == 1)

table(df_GSZB2$Welle, df_GSZB2$Land)

table(df_GSZB2$Welle[df_GSZB2$Land == "DE"], df_GSZB2$gisd_5[df_GSZB2$Land == "DE"], useNA = "ifany") #GISD: 0 = lowest deprivation, 1 = highest deprivation
table(df_GSZB2$sex, df_GSZB2$Welle, df_GSZB2$Land, useNA = "ifany")
table(df_GSZB2$edu_group, df_GSZB2$Welle, df_GSZB2$Land, useNA = "ifany")
table(df_GSZB2$agegroup, df_GSZB2$Welle, df_GSZB2$Land, useNA = "ifany")
table(df_GSZB2$DGURBA, df_GSZB2$Welle, df_GSZB2$Land, useNA = "ifany")

#weighted percentages
demograf_vars <- c("sex", "edu_group", "agegroup", "DGURBA", "gisd_5")

result_list <- list()

for (var in demograf_vars) {
  # Variable inkl. NA als Faktorlevel
  var_vec <- df_GSZB2[[var]]
  var_vec <- addNA(var_vec)  # NA als eigener Faktorlevel
  tab <- with(df_GSZB2, tapply(weights, list(Welle, Land, var_vec), sum, na.rm=TRUE))
  weighted_perc <- apply(tab, c(1, 2), function(x) 100 * x / sum(x, na.rm=TRUE))
  print(paste("Weighted percentages for", var, ":"))
  print(round(weighted_perc, 1))
}

#weighted percentages for sex, pooled
tab <- with(df_GSZB2, tapply(weights, sex, sum, na.rm=TRUE))
weighted_perc <- 100 * tab / sum(tab, na.rm = TRUE)


# summary of weights by Welle and Land
weights_summary <- df_GSZB2 %>%
  group_by(Welle, Land) %>%
  summarise(
    mean_weight = mean(weights, na.rm = TRUE),
    sd_weight = sd(weights, na.rm = TRUE),
    min_weight = min(weights, na.rm = TRUE),
    max_weight = max(weights, na.rm = TRUE)
  )

#weighted mean age by Welle and Land
weighted_mean_age <- df_GSZB2 %>%
  group_by(Welle, Land) %>%
  summarise(
    weighted_mean_age = sum(age_in_years * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE),
    weighted_sd_age = sqrt(sum(weights * (age_in_years - weighted_mean_age)^2, na.rm = TRUE) / sum(weights, na.rm = TRUE))
  )

# share of source pharmacy (referenced in discussion)
source_pharmacy_share <- df_GSZB2 %>%
  filter(!is.na(SOURCE.1)) %>%
  filter(!can_freq %in% c("seltener als einmal im Monat", "gar nicht")) %>%
  mutate(SOURCE.1 = haven::as_factor(SOURCE.1)) %>%
  select(Welle, SOURCE.1, weights) %>%
  group_by(Welle) %>%
  summarise(
    n = n(),
    n_pharm = sum(SOURCE.1 == "ausgewählt", na.rm = TRUE),
    unweighted_share_pharmacy = n_pharm / n,
    weighted_share_pharmacy = sum(as.numeric(SOURCE.1 == "ausgewählt") * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
  )

# Distribution of municipalities (map)
ags_welle_land <- df_GSZB2 %>%
  distinct(AGS, Welle, Land)

gemeinden_mit_daten <- degurba_mapping %>%
  inner_join(ags_welle_land, by = "AGS") %>%
  mutate(Welle = as.character(Welle))

laendergrenzen <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin %in% c("Germany", "Austria"))  # Nur DE & AT

welle_labels_fw <- c("1" = "t0", "2" = "t1")

map <- ggplot() +
  geom_sf(data = gemeinden_mit_daten, aes(fill = Land, color = Land), alpha = 0.7, size = 0.01) +
  geom_sf(data = laendergrenzen, fill = NA, color = "black", size = 0.5) +
  facet_wrap(~ Welle, labeller = labeller(Welle = welle_labels_fw)) +
  scale_fill_manual(values = colors_country, labels = country_labels) +
  scale_color_manual(values = colors_country, labels = country_labels) +
  theme_void() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)
        #plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        #plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Save the map
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "map_sample1", DATE, ".tiff"),
    plot = map,
    width = 107,                   
    height = 90,                  
    units = "mm",                  
    dpi = 300,
    bg = "white",
    device = "tiff")
}

# 12-month prevalence of cannabis, alcohol, tobacco use, and gambling 
aggregates_canuse <- df_GSZB2 %>%
  group_by(Land, Welle) %>%
  summarise(
    results = list(calculate_prevalence(cur_data(), "can_use_12M", "weights", n_boot = n_boot)),
    .groups = "drop"
  ) %>%
  unnest_wider(results)

aggregates_tobaccouse <- df_GSZB2 %>%
  group_by(Land, Welle) %>%
  summarise(
    results = list(calculate_prevalence(cur_data(), "tobacco_use_12M", "weights", n_boot = n_boot)),
    .groups = "drop"
  ) %>%
  unnest_wider(results)

aggregates_alcoholuse <- df_GSZB2 %>%
  group_by(Land, Welle) %>%
  summarise(
    results = list(calculate_prevalence(cur_data(), "alcohol_use_12M", "weights", n_boot = n_boot)),
    .groups = "drop"
  ) %>%
  unnest_wider(results)

aggregates_gambling <- df_GSZB2 %>%
  group_by(Land, Welle) %>%
  summarise(
    results = list(calculate_prevalence(cur_data(), "gambling_12M", "weights", n_boot = n_boot)),
    .groups = "drop"
  ) %>%
  unnest_wider(results)

# ==================================================================================================================================================================
# #Sample description - Sample 2 (monthly cannabis users without medical exemption)
# ==================================================================================================================================================================
df_DUIC_full <- subset(data, AS3 == 1)

# #all speeding cols to factors
# df_DUIC_full <- df_DUIC_full %>%
#   mutate(across(starts_with("SPEEDING"), as_factor))
# 
# ,
#          across(starts_with("SPEEDING"),
#                 ~ case_when(
#                   drives == FALSE ~ "nodriving",
#                   TRUE ~ as.character(.x)
#                 )))


# flow chart: exclusions
table(df_GSZB2$Welle, df_GSZB2$Land, df_GSZB2$can_use_12M)


table(
  df_GSZB2$Welle[df_GSZB2$can_use_12M == 1],
  df_GSZB2$Land[df_GSZB2$can_use_12M == 1],
  df_GSZB2$MEDICALUSE.01[df_GSZB2$can_use_12M == 1] == "Ausschließlich für medizinische Zwecke" &
    df_GSZB2$MEDICALUSE.02[df_GSZB2$can_use_12M == 1] == "Ja, mir wurde medizinisches Cannabis ärztlich verschrieben (bezahlt durch Krankenkasse oder als Selbstzahler)",
  useNA = "ifany"
)

table(
  df_GSZB2$Land[df_GSZB2$can_use_12M == 1 & df_GSZB2$Welle == 2],
  df_GSZB2$MEDICALUSE.01[df_GSZB2$can_use_12M == 1 & df_GSZB2$Welle == 2] == "Ausschließlich für medizinische Zwecke" &
    df_GSZB2$MEDICALUSE.02[df_GSZB2$can_use_12M == 1 & df_GSZB2$Welle == 2] == "Ja, mir wurde medizinisches Cannabis ärztlich verschrieben (bezahlt durch Krankenkasse oder als Selbstzahler)",
  df_GSZB2$DUIC30d_full[df_GSZB2$can_use_12M == 1 & df_GSZB2$Welle == 2], useNA = "ifany"
)



table(df_DUIC_full$Welle, df_DUIC_full$Land)
table(df_DUIC_full$sex, df_DUIC_full$Welle, df_DUIC_full$Land, useNA = "ifany")
table(df_DUIC_full$edu_group, df_DUIC_full$Welle, df_DUIC_full$Land, useNA = "ifany")
table(df_DUIC_full$agegroup, df_DUIC_full$Welle, df_DUIC_full$Land, useNA = "ifany")
table(df_DUIC_full$DGURBA, df_DUIC_full$Welle, df_DUIC_full$Land, useNA = "ifany")
table(df_DUIC_full$can_freq, df_DUIC_full$Welle, df_DUIC_full$Land, useNA = "ifany")
table(df_DUIC_full$DRIVERLICENSE_bin, df_DUIC_full$Welle, df_DUIC_full$Land, useNA = "ifany")

# unweighted percentages
demograf_vars_DUIC <- c("sex", "edu_group", "agegroup", "DGURBA", "can_freq", "DRIVERLICENSE_bin")

result_list_DUIC <- list()
for (var in demograf_vars_DUIC) {
  var_vec <- addNA(df_DUIC_full[[var]]) 
  
  tab_abs <- table(df_DUIC_full$Welle, df_DUIC_full$Land, var_vec)
  
  tab_perc <- apply(tab_abs, c(1, 2), function(x) 100 * x / sum(x, na.rm = TRUE))
  
  cat("\nUnweighted percentages for", var, ":\n")
  print(round(tab_perc, 1))
}

summary(df_DUIC_full$age_in_years)

# ==================================================================================================================================================================
# #Sample description - Sample 3 (past-year users with DUIC in the past 30 days, t₁ only)
# ==================================================================================================================================================================


mischkonsum_weights <- c(
  "Ausschließlich Cannabis konsumiert" = 1,
  "Gelegentlich neben Cannabis auch Alkohol oder andere Drogen konsumiert" = 0.75,
  "Meistens neben Cannabis auch Alkohol oder andere Drogen konsumiert" = 0.25,
  "Immer neben Cannabis auch Alkohol oder andere Drogen konsumiert" = 0
)


DUIC30_df <- df_GSZB2 %>%
  filter(
    !(MEDICALUSE.01 == "Ausschließlich für medizinische Zwecke" &
        MEDICALUSE.02 == "Ja, mir wurde medizinisches Cannabis ärztlich verschrieben (bezahlt durch Krankenkasse oder als Selbstzahler)"),
    can_freq != "gar nicht",
    DUIC30d_full == 1,
    Welle == 2
  ) %>%
  mutate(
    DUIC_freq_30_num = case_when(
      DUIC30d_full == 0 ~ 0,
      DUIC.FREQ.01M_wave2coding == "1 mal" ~ 1,
      DUIC.FREQ.01M_wave2coding == "2-3 mal" ~ 2.5,
      DUIC.FREQ.01M_wave2coding == "4-9 mal" ~ 6.5,
      DUIC.FREQ.01M_wave2coding == "10-15 mal" ~ 12.5,
      DUIC.FREQ.01M_wave2coding == "Mehr als 15 mal" ~ DUIC.FREQ.01M_exact,
      TRUE ~ NA_real_
    )) %>%
  filter(DUIC_freq_30_num < 99) %>% 
  filter(!is.na(DUIC_freq_30_num)) %>%
  #add new variable for the approx. DUIC episodes by cannabis only und psu
  mutate(
    DUIC_only_freq_30_num = DUIC_freq_30_num * mischkonsum_weights[DUIC.MIXEDUSE],
    DUIC_psu_freq_30_num = DUIC_freq_30_num * (1 - mischkonsum_weights[DUIC.MIXEDUSE]
    )) %>%
  filter(!is.na(DUIC_only_freq_30_num), 
         !is.na(DUIC_psu_freq_30_num)) 



table(DUIC30_df$Land)
table(DUIC30_df$sex, useNA = "ifany")
table(DUIC30_df$edu_group, useNA = "ifany")
table(DUIC30_df$agegroup, useNA = "ifany")
table(DUIC30_df$DGURBA, useNA = "ifany")
table(DUIC30_df$can_freq, useNA = "ifany")
table(DUIC30_df$DRIVERLICENSE_bin, useNA = "ifany")
# unweighted percentages
prop.tab <- function(var) {
  tab <- table(DUIC30_df[[var]], useNA = "ifany")
  perc <- prop.table(tab) * 100
  return(round(perc, 1))
}

#over demograf_vars_DUIC 
lapply(demograf_vars_DUIC, prop.tab)

summary(DUIC30_df$age_in_years)

# ==================================================================================================================================================================
# #12M can use for Welle and Land - weighted DiD
# ==================================================================================================================================================================

# parallel trends before survey? 

# Daten für Erwachsene (AT: 15–64 Jahre (EUDA), GER: 18–64 Jahre (ESA))
adults_data <- data.frame(
  Land = c(rep("AT", 3), rep("DE", 5)),
  Jahr = c(2008, 2015, 2020, 2009, 2012, 2015, 2018, 2021),
  can_use_prev12M = c(3.5, 6.4, 6.3, 4.8, 4.5, 6.1, 7.1, 8.8),
  cohort = "Adults"
)

# Daten für Jugendliche (ESPAD-Daten)
adolescents_data <- data.frame(
  Land = c("AT", "DE", "AT", "DE", "AT", "DE", "AT", "DE", "AT", "DE"),
  Jahr = c(2003, 2003, 2007, 2007, 2015, 2015, 2019, 2019, 2024, 2024),
  can_use_prev12M = c(9.3, 10.9, 6.2, 5.7, 9.2, 7.7, 11.3, 10.5, 6.1, 6.8),
  cohort = "Adolescents"
)

cannabis_use_data <- rbind(adults_data, adolescents_data)

#plot trends in both countries fpr adults and adolescents
trend_canuse_pre <- ggplot(cannabis_use_data, aes(x = Jahr, y = can_use_prev12M, color = Land, linetype = cohort)) +
  geom_line(size = 0.7) +
  geom_point(size = 1.5) +
  labs(x = "Year",
       y = "Past 12-month cannabis use prevalence",
       color = "",
       linetype = ""
  ) +
  scale_y_continuous(
    breaks = seq(0, 17, by = 2),
    labels = paste0(seq(0, 17, by = 2), "%"),
    limits = c(0, 17), expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = seq(2004, 2025, by = 2),
    limits = c(2002, 2025), expand = c(NA, 0)
  ) +
  scale_color_manual(values = colors_country, labels = country_labels) +
  scale_linetype_manual(values = c("twodash", "solid")) +
  theme_gdocs(base_family = "Calibri", base_size = 9) +
  theme(
    panel.grid.minor = element_line(color = "gray80", linetype = "dotted", linewidth = 0.3),
    axis.title.y = element_text(color= "black"),
    axis.title.x = element_text(color= "black"),
    axis.text.y = element_text(color= "black"),
    axis.text.x = element_text(color= "black", angle = 45, hjust = 1),
    strip.text = element_text(color= "black", face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(color= "black"),
    legend.text = element_text(color= "black"),
    legend.box.margin = margin(t = -9, unit = "pt") #abstand zwischen legenden und plot
    # legend.margin = margin(t = 6, unit = "pt") #abstand zwischen legendeninhalt und legendenrahmen
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
# theme(
#   panel.grid.minor = element_blank(),
#   axis.title.y = element_text(size = 15),
#   axis.title.x = element_text(size = 15),
#   axis.text = element_text(size = 13),
#   strip.text = element_text(size = 14, face = "bold"),
#   legend.position = "bottom",
#   legend.title = element_text(size = 15),
#   legend.text = element_text(size = 13),
#   legend.margin = margin(t = -5, b = 0, unit = "pt"),
#   legend.box.margin = margin(t = -10, unit = "pt"),
#   plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)
# )

#save as tiff
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "trend_canuse_pre_", DATE, ".tiff"),
    plot = trend_canuse_pre,
    width = 107,                   
    height = 90,                  
    units = "mm",                  
    dpi = 400,
    device = "tiff"
  )
}



#weighted DiD analysis for can use 12M
glm_diff_in_diff <- glm(can_use_12M ~ Welle * Land, data = df_GSZB2, family = "binomial", weights = weights)
summary(glm_diff_in_diff)
DiD_coef_canuse12m <- exp(cbind(OR = coef(glm_diff_in_diff), confint(glm_diff_in_diff)))
logisticPseudoR2s(glm_diff_in_diff)
glm_diff_in_diff_canuse_null <- glm(can_use_12M ~ 1, data = df_GSZB2, family = "binomial", weights = weights)
anova(glm_diff_in_diff_canuse_null, glm_diff_in_diff, test = "Chisq")

#unweighted DiD analysis for can use 12M
glm_diff_in_diff_unweighted <- glm(can_use_12M ~ Welle * Land, data = df_GSZB2, family = "binomial")
summary(glm_diff_in_diff_unweighted)
DiD_coef_canuse12m_unweighted <- exp(cbind(OR = coef(glm_diff_in_diff_unweighted), confint(glm_diff_in_diff_unweighted)))
logisticPseudoR2s(glm_diff_in_diff_unweighted)
glm_diff_in_diff_canuse_null_unweighted <- glm(can_use_12M ~ 1, data = df_GSZB2, family = "binomial")
anova(glm_diff_in_diff_canuse_null_unweighted, glm_diff_in_diff_unweighted, test = "Chisq")


#adjusted DiD analysis for can use 12M with confounders (but no weights)
#chi square tests for categorical variables
# sex 
table_sex_can <- table(df_GSZB2[df_GSZB2$Welle == 1,]$can_use_12, df_GSZB2[df_GSZB2$Welle == 1,]$sex, useNA = "ifany")
table_sex_can
prop.table(table_sex_can, margin = 2)
chisq.test(table_sex_can) # significant

#age groups
table_age_can <- table(df_GSZB2[df_GSZB2$Welle == 1,]$can_use_12, df_GSZB2[df_GSZB2$Welle == 1,]$agegroup)
table_age_can
prop.table(table_age_can, margin = 2)
chisq.test(table_age_can) # significant

#edu groups
table_edu_can <- table(df_GSZB2[df_GSZB2$Welle == 1,]$can_use_12, df_GSZB2[df_GSZB2$Welle == 1,]$edu_group)
table_edu_can
prop.table(table_edu_can, margin = 2)
chisq.test(table_edu_can) # significant

# degree of urbanization
table_urban_can <- table(
  df_GSZB2[df_GSZB2$Welle == 1,]$can_use_12,
  df_GSZB2[df_GSZB2$Welle == 1,]$DGURBA
)
table_urban_can
prop.table(table_urban_can, margin = 2)
chisq.test(table_urban_can) # significant

# alcohol freq
table_alcohol_can <- table(df_GSZB2[df_GSZB2$Welle == 1,]$can_use_12, df_GSZB2[df_GSZB2$Welle == 1,]$alcohol_freq)
table_alcohol_can
prop.table(table_alcohol_can, margin = 2)
chisq.test(table_alcohol_can) # significant

# # sport freq
table_sport_can <- table(df_GSZB2[df_GSZB2$Welle == 1,]$can_use_12, df_GSZB2[df_GSZB2$Welle == 1,]$sport_freq, useNA = "ifany")
table_sport_can
prop.table(table_sport_can, margin = 2)
chisq.test(table_sport_can) # significant

# gambling freq
table_gambling_can <- table(df_GSZB2[df_GSZB2$Welle == 1,]$can_use_12, df_GSZB2[df_GSZB2$Welle == 1,]$gambling_freq, useNA = "ifany")
table_gambling_can
prop.table(table_gambling_can, margin = 2)
chisq.test(table_gambling_can) # significant

# Sensitivity analysis 2: adjusted DiD analysis for can use 12M with sociodemographic confounders (but no weights)
glm_diff_in_diff_adj_con <- glm(can_use_12M ~ agegroup + sex + edu_group + DGURBA + Welle * Land, data = df_GSZB2, family = "binomial")
summary(glm_diff_in_diff_adj_con)
vif(glm_diff_in_diff_adj_con)
logisticPseudoR2s(glm_diff_in_diff_adj_con)
df_GSZB2_complete <- df_GSZB2[complete.cases(df_GSZB2[, c("can_use_12M", "agegroup", "sex", "edu_group", "DGURBA")]), ]
glm_diff_in_diff_canuse_null_unweighted_complete <- glm(can_use_12M ~ 1, data = df_GSZB2_complete, family = "binomial")
anova(glm_diff_in_diff_canuse_null_unweighted_complete, glm_diff_in_diff_adj_con, test = "Chisq")
DiD_coef_canuse12m_adjust_con <- exp(cbind(OR = coef(glm_diff_in_diff_adj_con), confint(glm_diff_in_diff_adj_con)))

# Sensitivity analysis 3: adjusted DiD analysis for can use 12M with sociodemographic and behavioral confounders (but no weights)
glm_diff_in_diff_adj_con_ext <- glm(can_use_12M ~ agegroup + sex + edu_group + DGURBA + gambling_freq + alcohol_freq + sport_freq + Welle * Land, data = df_GSZB2, family = "binomial") #alcohol_freq + 
summary(glm_diff_in_diff_adj_con_ext)
vif(glm_diff_in_diff_adj_con_ext)
logisticPseudoR2s(glm_diff_in_diff_adj_con_ext)
anova(glm_diff_in_diff_canuse_null_unweighted_complete, glm_diff_in_diff_adj_con_ext, test = "Chisq")
# save as html table
# if (dataexport) {
#   library(kableExtra)
#   kable(DiD_adj_coef_canuse12m, format = "html", digits = 2, caption = "Weighted DiD Analysis for Past-Year Cannabis Use") %>%
#     kable_styling("striped", full_width = F) %>%
#     save_kable(file = paste0(folder_path_tables, "DiD_w_can_use_12M_", DATE, ".html"))
# }

# ==================================================================================================================================================================
# #12M prev of DUIC full Welle and Land - DiD
# ==================================================================================================================================================================

#descriptive statistics for DUIC full prevalence with CIs

#DUIC 12m prev
aggregates_DUIC <- df_DUIC_full %>%
  group_by(Land, Welle) %>%
  summarise(
    results = list(calculate_prevalence(cur_data(), "DUIC12m_full", "weights", n_boot = n_boot)),
    .groups = "drop"
  ) %>%
  unnest_wider(results)

##bivariate analysis to find potential confounders

df_wave1 <- subset(data, Welle == 1 & AS3 == 1) 

#chi square tests for categorical variables
# sex 
table_sex <- table(df_wave1$DUIC12m_full, df_wave1$sex, useNA = "ifany")
table_sex
prop.table(table_sex, margin = 2)
chisq.test(table_sex) # significant

#edu groups
table_edu <- table(df_wave1$DUIC12m_full, df_wave1$edu_group)
table_edu
prop.table(table_edu, margin = 2)
chisq.test(table_edu) # significant
#fisher test because of small n's in some cells
fisher.test(table_edu)

# age groups
table_age <- table(df_wave1$DUIC12m_full, df_wave1$agegroup)
table_age
prop.table(table_age, margin = 2)
chisq.test(table_age) # significant

# driver license
table_driver_license <- table(df_wave1$DUIC12m_full, df_wave1$DRIVERLICENSE_bin)
table_driver_license
prop.table(table_driver_license, margin = 2)
chisq.test(table_driver_license) # significant

# can frequency
df_wave1$can_freq <- factor(
  df_wave1$can_freq,
  levels = c("mindestens einmal im Monat",
             "mindestens einmal pro Woche",
             "(fast) täglich")  # leere levels entfernen
)
table_can_freq <- table(df_wave1$DUIC12m_full, df_wave1$can_freq)
table_can_freq
prop.table(table_can_freq, margin = 2)
chisq.test(table_can_freq) #not significant

# degree of urbanization
table_urban <- table(
  df_wave1$DUIC12m_full,
  df_wave1$DGURBA
)
table_urban
prop.table(table_urban, margin = 2)
chisq.test(table_urban) # not significant

# alcohol freq
table_alcohol <- table(df_wave1$DUIC12m_full, df_wave1$alcohol_freq)
table_alcohol
prop.table(table_alcohol, margin = 2)
chisq.test(table_alcohol) # not significant, however there is a clear asscending trend with higher frequency -> keep as confounder


#sport freq
table_sport <- table(df_wave1$DUIC12m_full, df_wave1$sport_freq, useNA = "ifany")
table_sport
prop.table(table_sport, margin = 2)
chisq.test(table_sport) # not significant and small n's -> fisher test
fisher.test(table_sport)

#gambling freq
table_gambling <- table(df_wave1$DUIC12m_full, df_wave1$gambling_freq, useNA = "ifany")
table_gambling
prop.table(table_gambling, margin = 2)
chisq.test(table_gambling) 

#DUIA 12m
table_DUIA <- table(df_wave1$DUIC12m_full, df_wave1$DUIA12M, useNA = "ifany")
table_DUIA
prop.table(table_DUIA, margin = 2)
chisq.test(table_DUIA) # significant

#sedatives
table_sedatives <- table(df_wave1$DUIC12m_full, df_wave1$use_sedatives)
table_sedatives
prop.table(table_sedatives, margin = 2)
chisq.test(table_sedatives) # significant, but no data in Austria in Welle 1 -> cannot be used as confounder
#correlation with alc freq
#table_sedatives_alc <- prop.table(table(df_wave1$use_sedatives, df_wave1$alcohol_freq), margin = 2)

# DID analysis for DUIC full prevalence

df_DUIC_full$Land <- relevel(factor(df_DUIC_full$Land), ref = "AT")
df_DUIC_full$Welle <- factor(df_DUIC_full$Welle, levels = c(1, 2))  # Welle 1 als Referenz

DiD_unadjusted <- glm(DUIC12m_full ~ Welle * Land, data = df_DUIC_full, family = "binomial")
summary(DiD_unadjusted)
logisticPseudoR2s(DiD_unadjusted)
DiD_null <- glm(DUIC12m_full ~ 1, data = df_DUIC_full, family = "binomial")
anova(DiD_null, DiD_unadjusted, test = "Chisq")
DiD_unadj_coef <- exp(cbind(OR = coef(DiD_unadjusted), confint(DiD_unadjusted)))


# Adjusted DiD analysis 1 (only covariates)
DiD_adjusted <- glm(DUIC12m_full ~ DRIVERLICENSE_bin + agegroup + sex + edu_group + Welle * Land, data = df_DUIC_full, family = "binomial")
summary(DiD_adjusted)
logisticPseudoR2s(DiD_adjusted)
df_DUIC12m_full_complete <- df_DUIC_full[complete.cases(df_DUIC_full[, c("agegroup", "sex", "edu_group")]), ]
DiD_null_complete <- glm(DUIC12m_full ~ 1, data = df_DUIC12m_full_complete, family = "binomial")
anova(DiD_null_complete, DiD_adjusted, test = "Chisq")
vif(DiD_adjusted)

# # confidence intervals
# DiD_adj_coef <- exp(cbind(OR = coef(DiD_adjusted), confint(DiD_adjusted)))
# if (dataexport) {
#   kable(DiD_adj_coef, format = "html", digits = 2, caption = "Adjusted DiD Analysis for DUIC Full Prevalence") %>%
#     kable_styling("striped", full_width = F) %>%
#     save_kable(file = paste0(folder_path_tables, "DiD_adjusted_DUIC_full_", DATE, ".html"))
# }


# Adjusted DiD analysis 2 (sociademographic + behavioral covariates)
DiD_adjusted_2 <- glm(DUIC12m_full ~ DRIVERLICENSE_bin + agegroup + sex + edu_group + alcohol_freq + gambling_freq + Welle * Land, data = df_DUIC_full, family = "binomial")
summary(DiD_adjusted_2)
vif(DiD_adjusted_2)
logisticPseudoR2s(DiD_adjusted_2)
anova(DiD_null_complete, DiD_adjusted_2, test = "Chisq")



# calculate estimated marginal means
emm <- emmeans(DiD_adjusted, ~ Welle * Land, type = "response")
emm_df <- as.data.frame(emm)

# estimated probabilities for each country and wave
p_DE_t0 <- emm_df$prob[emm_df$Welle == 1 & emm_df$Land == "DE"]
p_DE_t1 <- emm_df$prob[emm_df$Welle == 2 & emm_df$Land == "DE"]
p_AT_t0 <- emm_df$prob[emm_df$Welle == 1 & emm_df$Land == "AT"]
p_AT_t1 <- emm_df$prob[emm_df$Welle == 2 & emm_df$Land == "AT"]

# absolute differences in percentage points for each country
diff_DE <- round((p_DE_t1 - p_DE_t0) * 100, 2)
diff_AT <- round((p_AT_t1 - p_AT_t0) * 100, 2)

# DiD: absolute difference in difference in percentage points
did_abs_percent <- round((diff_DE - diff_AT), 2)

diff_DE   
diff_AT  
did_abs_percent



# ==================================================================================================================================================================
# Plot canuse and DUIC aggregates
# ==================================================================================================================================================================

aggregates_plot <- bind_rows(
  aggregates_canuse %>%
    select(Land, Welle, weighted_prop, weighted_ci_lower, weighted_ci_upper) %>%
    rename(
      prop = weighted_prop,
      lower = weighted_ci_lower,
      upper = weighted_ci_upper
    ) %>%
    mutate(variable = "a) Cannabis use"),
  
  aggregates_DUIC %>%
    select(Land, Welle, unweighted_prop, unweighted_ci_lower, unweighted_ci_upper) %>%
    rename(
      prop = unweighted_prop,
      lower = unweighted_ci_lower,
      upper = unweighted_ci_upper
    ) %>%
    mutate(variable = "b) DUIC")
)


#plot 
plot_DUIC_canuse <- aggregates_plot %>%
  ggplot(aes(x = as.factor(Welle), y = prop, color = Land, group = Land)) +
  geom_point(size = 0.7, position = position_dodge(width = 0.05)) +
  geom_line(linewidth = 0.5, position = position_dodge(width = 0.05)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.09, size = 0.4, position = position_dodge(width = 0.05)) +
  geom_label(
    aes(label = scales::percent(prop, accuracy = 0.1)),
    position = position_dodge(width = 0.81),
    vjust = -0.5,
    size = 3.1,
    show.legend = FALSE
  ) +
  facet_wrap(~ variable, ncol = 1, scales = "free_x") +
  scale_y_continuous(limits = c(0, 0.42),
                     breaks = seq(0, 0.42, by = 0.1),
                     expand = c(0,0),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    y = "12-month prevalence",
    x = "",
    color = "",
  ) +
  scale_color_manual(values = colors_country, labels = country_labels) +
  scale_x_discrete(labels = welle_labels) +
  theme_gdocs(base_family = "Calibri", base_size = 10) +
  theme(
    panel.grid.minor = element_line(color = "gray80", linetype = "dotted", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(6, "mm"),
    axis.title.y = element_text(color= "black"),
    axis.text = element_text(color= "black"),
    strip.text = element_text(color= "black", face = "bold"), 
    legend.position = "right",
    legend.title = element_text(color= "black"),
    legend.text = element_text(color= "black", hjust = 0.7),
    legend.box.margin = margin(t = -15, unit = "mm"), #abstand zwischen legenden und plot
    legend.margin = margin(t = -5, b = 0, unit = "pt")
    #  legend.box.margin = margin(t = -10, unit = "pt"),
    #plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)
    #  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

#export as tiff
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "DUIC_canuse_prevalence_", DATE, ".tiff"),
    plot = plot_DUIC_canuse,
    bg = "white",
    width = 120,                   
    height = 110,                  
    units = "mm",                  
    dpi = 400,
    device = "tiff"
  )
}


# for DSK
plot_DUIC_DSK <- aggregates_plot %>% filter (variable == "b) DUIC") %>%
  ggplot(aes(x = as.factor(Welle), y = prop, color = Land, group = Land)) +
  geom_line(size = 1, position = position_dodge(width = 0.05)) +
  geom_point(size = 3, position = position_dodge(width = 0.05)) +  # Punkte leicht versetzen damit sch die KIs nicht überlappen
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.04, alpha = 0.5, size = 0.4, position = position_dodge(width = 0.05))  + # Fehlerbalken versetzen
  geom_label(aes(label = scales::percent(prop, accuracy = 0.1)),
             size = 4,
             label.size = 0.2,
             fill = "white",
             color = "black",
             hjust = -0.2,
             vjust = -1,
             fontface = "bold",
             position = position_dodge(width = 0.1)) +  # Labels versetzen
  labs(title = str_wrap("12-Monats-Prävalenz DUIC unter mind. monatlich Cannabiskonsumierenden in DE und AT", width = 60),
       subtitle = str_wrap("Vergleich pre/post Legalisierung, exkl. ausschl. med. Konsum mit Rezept (DE, t0: n=393, t1: n=589; AT, t0: n=86, t1: n=92)", width = 75),
       y = "12-Monats-Prävalenz",
       x = "",
       color = "Land") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 0.4, by = 0.05), limits = c(0, 0.4), expand = c(0, 0)) +
  scale_x_discrete(labels = welle_labels) +
  scale_color_manual(values = colors_country) +
  theme_gdocs(base_family = "Aptos", base_size = 10) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    strip.text = element_text(color = "black"), 
    plot.title = element_text(color = "black"),
    plot.subtitle = element_text(color = "black"),
    plot.caption = element_text(color = "black"),
    panel.grid.minor = element_line(color = "gray80", linetype = "dotted", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#save
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "DUIC_canuse_prevalence_DSK_", DATE, ".tiff"),
    plot = plot_DUIC_DSK,
    width = 160,                   
    height = 110,                  
    units = "mm",                  
    dpi = 400,
    device = "tiff",
    bg = "white"
  )
}

# ==================================================================================================================================================================
# #Sensitivity analysis: 30D prev of DUIC full Welle and Land - DiD
# ==================================================================================================================================================================

DiD_DUIC30d_unadjusted <- glm(DUIC30d_full ~ Welle * Land, data = df_DUIC_full, family = "binomial")
summary(DiD_DUIC30d_unadjusted)

# confidence intervals
DiD_DUIC30d_unadj_coef <- exp(cbind(OR = coef(DiD_DUIC30d_unadjusted), confint(DiD_DUIC30d_unadjusted)))
#save as html table
# if (dataexport) {
#   library(kableExtra)
#   kable(DiD_unadj_coef, format = "html", digits = 2, caption = "Unadjusted DiD Analysis for DUIC Full Prevalence") %>%
#     kable_styling("striped", full_width = F) %>%
#     save_kable(file = paste0(folder_path_tables, "DiD_unadjusted_DUIC_full_", DATE, ".html"))
# }

# Adjusted DiD analysis
DiD_DUIC30d_adjusted <- glm(DUIC30d_full ~ DRIVERLICENSE_bin + agegroup + sex + edu_group + Welle * Land, data = df_DUIC_full, family = "binomial")
summary(DiD_DUIC30d_adjusted)
logisticPseudoR2s(DiD_DUIC30d_adjusted)
DiD_DUIC30d_null <- glm(DUIC30d_full ~ 1, data = df_DUIC12m_full_complete, family = "binomial")
anova(DiD_DUIC30d_null, DiD_DUIC30d_adjusted, test = "Chisq")

# confidence intervals
DiD_DUIC30d_adj_coef <- exp(cbind(OR = coef(DiD_DUIC30d_adjusted), confint(DiD_DUIC30d_adjusted)))
if (dataexport) {
  kable(DiD_DUIC30d_adj_coef, format = "html", digits = 2, caption = "Adjusted DiD Analysis for past 30 day DUIC Full Prevalence") %>%
    kable_styling("striped", full_width = F) %>%
    save_kable(file = paste0(folder_path_tables, "DiD_adjusted_DUIC30d_full_", DATE, ".html"))
}

#VIF adjusted model
vif(DiD_DUIC30d_adjusted)

# ==================================================================================================================================================================
# Negative Control - DiD (eating fish in past 12M) - Sample 1
# # ==================================================================================================================================================================# same for df_GSZB2

table(df_GSZB2$fish_12M, df_GSZB2$Welle, df_GSZB2$Land)
df_GSZB2 %>%
  group_by(Land, Welle, fish_12M) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Land, Welle) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  arrange(Land, Welle, desc(percent))

DiD_fish_unw_GSP <- glm(fish_12M ~ Welle * Land, data = df_GSZB2, family = "binomial")
summary(DiD_fish_unw_GSP)
# confidence intervals
DiD_fish_unw_coef_GSP <- exp(cbind(OR = coef(DiD_fish_unw_GSP), confint(DiD_fish_unw_GSP)))
# if (dataexport) {
#   kable(DiD_fish_unw_coef_GSP, format = "html", digits = 2, caption = "Unweighted DiD Analysis for eating fish in past 12M (GSZB2)") %>%
#     kable_styling("striped", full_width = F) %>%
#     save_kable(file = paste0(folder_path_tables, "DiD_unw_fish_12M_GSZB2_", DATE, ".html"))
# }

# Aweighted DiD analysis for eating fish in past 12M (GSZB2)
DiD_fish_w_GSP <- glm(fish_12M ~ Welle * Land, data = df_GSZB2, family = "binomial", weights = weights)
summary(DiD_fish_w_GSP)
# confidence intervals
DiD_fish_w_GSP_coef <- exp(cbind(OR = coef(DiD_fish_w_GSP), confint(DiD_fish_w_GSP)))
if (dataexport) {
  kable(DiD_fish_w_GSP_coef, format = "html", digits = 2, caption = "Weighted DiD Analysis for eating fish in past 12M (GSZB2)") %>%
    kable_styling("striped", full_width = F) %>%
    save_kable(file = paste0(folder_path_tables, "DiD_w_fish_12M_GSZB2_", DATE, ".html"))
}

# ==================================================================================================================================================================
# Negative Control - DiD (doing sports in past 12M) - Sample 2
# # ==================================================================================================================================================================
# DiD analysis for doing sports in past 12M
table(df_DUIC_full$sport_12M, df_DUIC_full$Welle, df_DUIC_full$Land)
DiD_sport_unadjusted <- glm(sport_12M ~ Welle * Land, data = df_DUIC_full, family = "binomial")
summary(DiD_sport_unadjusted)
# confidence intervals
DiD_sport_unadj_coef <- exp(cbind(OR = coef(DiD_sport_unadjusted), confint(DiD_sport_unadjusted)))
# if (dataexport) {
#   kable(DiD_sport_unadj_coef, format = "html", digits = 2, caption = "Unadjusted DiD Analysis for Doing Sports in Past 12M") %>%
#     kable_styling("striped", full_width = F) %>%
#     save_kable(file = paste0(folder_path_tables, "DiD_unadjusted_sport_12M_", DATE, ".html"))
# }
# Adjusted DiD analysis for doing sports in past 12M
DiD_sport_adjusted <- glm(sport_12M ~ DRIVERLICENSE_bin + agegroup + sex + edu_group + Welle * Land, data = df_DUIC_full, family = "binomial")
summary(DiD_sport_adjusted)
# confidence intervals
DiD_sport_adj_coef <- exp(cbind(OR = coef(DiD_sport_adjusted), confint(DiD_sport_adjusted)))
if (dataexport) {
  kable(DiD_sport_adj_coef, format = "html", digits = 2, caption = "Adjusted DiD Analysis for Doing Sports in Past 12M") %>%
    kable_styling("striped", full_width = F) %>%
    save_kable(file = paste0(folder_path_tables, "DiD_adjusted_sport_12M_", DATE, ".html"))
}

# ==================================================================================================================================================================
# DUIC episodes by cannabis use frequency and polysubstance use
# ==================================================================================================================================================================

DUIC30_df <- df_GSZB2 %>%
  filter(
    !(MEDICALUSE.01 == "Ausschließlich für medizinische Zwecke" &
        MEDICALUSE.02 == "Ja, mir wurde medizinisches Cannabis ärztlich verschrieben (bezahlt durch Krankenkasse oder als Selbstzahler)"),
    can_freq != "gar nicht",
    #  DUIC30d_full == 1,
    Welle == 2
  ) %>%
  mutate(
    DUIC_freq_30_num = case_when(
      DUIC30d_full == 0 ~ 0,
      DUIC.FREQ.01M_wave2coding == "1 mal" ~ 1,
      DUIC.FREQ.01M_wave2coding == "2-3 mal" ~ 2.5,
      DUIC.FREQ.01M_wave2coding == "4-9 mal" ~ 6.5,
      DUIC.FREQ.01M_wave2coding == "10-15 mal" ~ 12.5,
      DUIC.FREQ.01M_wave2coding == "Mehr als 15 mal" ~ DUIC.FREQ.01M_exact,
      TRUE ~ NA_real_
    )) %>%
  #filter(DUIC_freq_30_num < 99 | is.na(DUIC_freq_30_num)) %>% 
  #add new variable for the approx. DUIC episodes by cannabis only und psu
  mutate(
    DUIC_only_freq_30_num = DUIC_freq_30_num * mischkonsum_weights[DUIC.MIXEDUSE],
    DUIC_psu_freq_30_num = DUIC_freq_30_num * (1 - mischkonsum_weights[DUIC.MIXEDUSE]),
    analyseDUIC = ifelse(
      DUIC30d_full == 1 & !is.na(DUIC_freq_30_num) & !is.na(DUIC_only_freq_30_num)
      & !is.na(DUIC_psu_freq_30_num) & DUIC_freq_30_num < 99,
      1, 0
    ))


#check
summary(DUIC30_df %>% filter(analyseDUIC == 1) %>% pull(DUIC_only_freq_30_num))
sd(DUIC30_df %>% filter(analyseDUIC == 1) %>% pull(DUIC_only_freq_30_num))
summary(DUIC30_df %>% filter(analyseDUIC == 1) %>% pull(DUIC_psu_freq_30_num))
summary(DUIC30_df %>% filter(analyseDUIC == 1) %>% pull(DUIC_freq_30_num))
sd(DUIC30_df %>% filter(analyseDUIC == 1) %>% pull(DUIC_freq_30_num))



# calculate DUIC episodes by cannabis use frequency and polysubstance use

sum_duicep <- DUIC30_df %>%
  group_by(can_freq) %>%
  summarise(
    n_users = n(),
    n_duicep = sum(DUIC_freq_30_num[analyseDUIC == 1], na.rm = TRUE), 
    n_duicep_only = sum(DUIC_only_freq_30_num[analyseDUIC == 1], na.rm = TRUE),
    n_duic_psu = sum(DUIC_psu_freq_30_num[analyseDUIC == 1], na.rm = TRUE)
  ) %>%
  mutate(
    share_users = n_users / sum(n_users),
    share_duicep = n_duicep / sum(n_duicep),
    share_duicep_only = n_duicep_only / sum(n_duicep_only),
    share_duic_psu = n_duic_psu / sum(n_duic_psu)
  ) %>%
  ungroup() 

sum(sum_duicep$n_duicep)
sum(sum_duicep$n_duicep_only)
sum(sum_duicep$n_duic_psu)

sum(sum_duicep$n_duicep_only)/sum(sum_duicep$n_duicep)
sum(sum_duicep$n_duic_psu)/sum(sum_duicep$n_duicep)


#prepare data to have one column with "pop_of_int" with the values "users", "duicep_only" and "duic_psu"
pop_duicep_long <- sum_duicep %>%
  pivot_longer(
    cols = c(share_users, share_duicep, share_duicep_only, share_duic_psu),
    names_to = "grundgesamtheit",
    values_to = "share"
  ) %>%
  mutate(
    grundgesamtheit = dplyr::recode(grundgesamtheit,
                                    "share_users" = "Cannabis users",
                                    "share_duicep" = "DUIC episodes",
                                    "share_duicep_only" = "DUIC-only episodes",
                                    "share_duic_psu" = "DUIC-PSU episodes"),
    grundgesamtheit = factor(grundgesamtheit, levels = c("Cannabis users",
                                                         "DUIC episodes",
                                                         "DUIC-only episodes",
                                                         "DUIC-PSU episodes"),
                             labels = c("Cannabis users", 
                                        "DUIC episodes",
                                        "DUIC(–) episodes",
                                        "DUIC(+) episodes")),
    can_freq = dplyr::recode(can_freq,
                             `seltener als einmal im Monat` = "Less than monthly",
                             `mindestens einmal im Monat` = "Monthly",
                             `mindestens einmal pro Woche` = "Weekly",
                             `(fast) täglich` = "(Almost) daily"),
    can_freq = factor(can_freq, levels = c("Less than monthly", "Monthly", "Weekly", "(Almost) daily"))
  )

#plot stacked bar chart
fig1_w <- ggplot(pop_duicep_long %>% subset(grundgesamtheit != "DUIC episodes"), aes(x = grundgesamtheit, y = share, fill = can_freq)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_label(aes(x = grundgesamtheit, y = share, group = can_freq, label = scales::percent(share, accuracy = 1)),
             position = position_stack(vjust = 0.5), color = "black", fill = alpha("white", 0.9), label.size = 0, size = 3.3) +
  scale_fill_manual(values = blue_colors, labels = str_wrap(c("Less than monthly", "Monthly", "Weekly", "(Almost) daily"), width = 10)) +
  scale_x_discrete(labels = str_wrap(c("Cannabis users", "DUIC(–) episodes", "DUIC(+) episodes"), width = 10)) +
  labs(
    x = "",
    y = "Share",
    fill = "Cannabis use frequency",
    title = "",
  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), labels = scales::percent, breaks = seq(0, 1, by = 0.5)) +
  theme_gdocs(base_family = "Calibri", base_size = 10) +
  theme(
    panel.grid.minor = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black", vjust = -0.3),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    #plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  #legend in one row
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# save plot as tiff
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "pop_duiconlyep_duicpsuep_w_", DATE, ".tiff"),
    plot = fig1_w,
    width = 107,                   
    height = 110,                  
    units = "mm",                  
    dpi = 400,
    device = "tiff",
    bg = "white"
  )
}

#for DSK
#plot stacked bar chart
fig1_w_DSK <- ggplot(pop_duicep_long %>% subset(grundgesamtheit != "DUIC episodes"), aes(x = grundgesamtheit, y = share, fill = can_freq)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_label(aes(x = grundgesamtheit, y = share, group = can_freq, label = scales::percent(share, accuracy = 1)),
             position = position_stack(vjust = 0.5), color = "black", fill = alpha("white", 0.9), label.size = 0, size = 3.3) +
  scale_fill_manual(values = blue_colors, labels = str_wrap(c("< monatlich", "monatlich", "wöchentlich", "(fast) tägich"), width = 12)) +
  scale_x_discrete(labels = str_wrap(c("Cannabis- konsument:innen (n=1160)", "DUIC (–) Fahrten (n=392)", "DUIC (+) Fahrten (n=108)"), width = 10)) +
  labs(
    x = "",
    y = "Anteil",
    fill = "Konsumhäufigkeit",
    title = str_wrap("Verteilung der Cannabiskonsument:innen, DUIC(–)- und DUIC(+)-Fahrten nach Konsumhäufigkeit", width = 50),
    subtitle =  str_wrap("unter Personen mit DUIC in den letzten 30 Tagen in DE und AT (n=86), nach der Legalisierung", width = 100))+
  scale_y_continuous(expand = c(0,0), limits = c(0,1), labels = scales::percent, breaks = seq(0, 1, by = 0.5)) +
  theme_gdocs(base_family = "Aptos", base_size = 10) +
  theme(
    plot.title = element_text(color = "black"),
    plot.subtitle = element_text(color = "black"),
    panel.grid.minor = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black", vjust = -0.3),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"))
#plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)
# ) +
# #legend in one row
# guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# save plot as tiff
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "pop_duiconlyep_duicpsuep_w_DSK", DATE, ".tiff"),
    plot = fig1_w_DSK,
    width = 159,                   
    height = 110,                  
    units = "mm",                  
    dpi = 400,
    device = "tiff",
    bg = "white"
  )
}



# Sensititvity analysis: only German data

DUIC30_df_GER <- df_GSZB2 %>%
  filter(
    Land == "DE",
    !(MEDICALUSE.01 == "Ausschließlich für medizinische Zwecke" &
        MEDICALUSE.02 == "Ja, mir wurde medizinisches Cannabis ärztlich verschrieben (bezahlt durch Krankenkasse oder als Selbstzahler)"),
    can_freq != "gar nicht",
    #  DUIC30d_full == 1,
    Welle == 2
  ) %>%
  mutate(
    DUIC_freq_30_num = case_when(
      DUIC30d_full == 0 ~ 0,
      DUIC.FREQ.01M_wave2coding == "1 mal" ~ 1,
      DUIC.FREQ.01M_wave2coding == "2-3 mal" ~ 2.5,
      DUIC.FREQ.01M_wave2coding == "4-9 mal" ~ 6.5,
      DUIC.FREQ.01M_wave2coding == "10-15 mal" ~ 12.5,
      DUIC.FREQ.01M_wave2coding == "Mehr als 15 mal" ~ DUIC.FREQ.01M_exact,
      TRUE ~ NA_real_
    )) %>%
  #filter(DUIC_freq_30_num < 99 | is.na(DUIC_freq_30_num)) %>% 
  #add new variable for the approx. DUIC episodes by cannabis only und psu
  mutate(
    DUIC_only_freq_30_num = DUIC_freq_30_num * mischkonsum_weights[DUIC.MIXEDUSE],
    DUIC_psu_freq_30_num = DUIC_freq_30_num * (1 - mischkonsum_weights[DUIC.MIXEDUSE]),
    analyseDUIC = ifelse(
      DUIC30d_full == 1 & !is.na(DUIC_freq_30_num) & !is.na(DUIC_only_freq_30_num)
      & !is.na(DUIC_psu_freq_30_num) & DUIC_freq_30_num < 99,
      1, 0
    ))


#check
summary(DUIC30_df_GER %>% filter(analyseDUIC == 1) %>% pull(DUIC_only_freq_30_num))
sd(DUIC30_df_GER %>% filter(analyseDUIC == 1) %>% pull(DUIC_only_freq_30_num))
summary(DUIC30_df_GER %>% filter(analyseDUIC == 1) %>% pull(DUIC_psu_freq_30_num))
summary(DUIC30_df_GER %>% filter(analyseDUIC == 1) %>% pull(DUIC_freq_30_num))
sd(DUIC30_df_GER %>% filter(analyseDUIC == 1) %>% pull(DUIC_freq_30_num))



# calculate DUIC episodes by cannabis use frequency and polysubstance use

sum_duicep_GER <- DUIC30_df_GER %>%
  group_by(can_freq) %>%
  summarise(
    n_users = n(),
    n_duicep = sum(DUIC_freq_30_num[analyseDUIC == 1], na.rm = TRUE), 
    n_duicep_only = sum(DUIC_only_freq_30_num[analyseDUIC == 1], na.rm = TRUE),
    n_duic_psu = sum(DUIC_psu_freq_30_num[analyseDUIC == 1], na.rm = TRUE)
  ) %>%
  mutate(
    share_users = n_users / sum(n_users),
    share_duicep = n_duicep / sum(n_duicep),
    share_duicep_only = n_duicep_only / sum(n_duicep_only),
    share_duic_psu = n_duic_psu / sum(n_duic_psu)
  ) %>%
  ungroup() 

sum(sum_duicep_GER$n_duicep)
sum(sum_duicep_GER$n_duicep_only)
sum(sum_duicep_GER$n_duic_psu)

sum(sum_duicep_GER$n_duicep_only)/sum(sum_duicep$n_duicep)
sum(sum_duicep_GER$n_duic_psu)/sum(sum_duicep$n_duicep)


#prepare data to have one column with "pop_of_int" with the values "users", "duicep_only" and "duic_psu"
pop_duicep_long_GER <- sum_duicep_GER %>%
  pivot_longer(
    cols = c(share_users, share_duicep, share_duicep_only, share_duic_psu),
    names_to = "grundgesamtheit",
    values_to = "share"
  ) %>%
  mutate(
    grundgesamtheit = dplyr::recode(grundgesamtheit,
                                    "share_users" = "Cannabis users",
                                    "share_duicep" = "DUIC episodes",
                                    "share_duicep_only" = "DUIC-only episodes",
                                    "share_duic_psu" = "DUIC-PSU episodes"),
    grundgesamtheit = factor(grundgesamtheit, levels = c("Cannabis users",
                                                         "DUIC episodes",
                                                         "DUIC-only episodes",
                                                         "DUIC-PSU episodes"),
                             labels = c("Cannabis users", 
                                        "DUIC episodes",
                                        "DUIC(–) episodes",
                                        "DUIC(+) episodes")),
    can_freq = dplyr::recode(can_freq,
                             `seltener als einmal im Monat` = "Less than monthly",
                             `mindestens einmal im Monat` = "Monthly",
                             `mindestens einmal pro Woche` = "Weekly",
                             `(fast) täglich` = "(Almost) daily"),
    can_freq = factor(can_freq, levels = c("Less than monthly", "Monthly", "Weekly", "(Almost) daily"))
  )

#plot stacked bar chart
fig1_w_GER <- ggplot(pop_duicep_long_GER %>% subset(grundgesamtheit != "DUIC episodes"), aes(x = grundgesamtheit, y = share, fill = can_freq)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_label(aes(x = grundgesamtheit, y = share, group = can_freq, label = scales::percent(share, accuracy = 1)),
             position = position_stack(vjust = 0.5), color = "black", fill = alpha("white", 0.9), label.size = 0, size = 3.3) +
  scale_fill_manual(values = blue_colors, labels = str_wrap(c("Less than monthly", "Monthly", "Weekly", "(Almost) daily"), width = 10)) +
  scale_x_discrete(labels = str_wrap(c("Cannabis users", "DUIC(–) episodes", "DUIC(+) episodes"), width = 10)) +
  labs(
    x = "",
    y = "Share",
    fill = "Cannabis use frequency",
    title = "",
  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), labels = scales::percent, breaks = seq(0, 1, by = 0.5)) +
  theme_gdocs(base_family = "Calibri", base_size = 10) +
  theme(
    panel.grid.minor = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black", vjust = -0.3),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    #plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  #legend in one row
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# save plot as tiff
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "pop_duiconlyep_duicpsuep_w_GER_", DATE, ".tiff"),
    plot = fig1_w_GER,
    width = 107,                   
    height = 110,                  
    units = "mm",                  
    dpi = 400,
    device = "tiff",
    bg = "white"
  )
}

#and as svg
if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "pop_duiconlyep_duicpsuep_w_GER_", DATE, ".svg"),
    plot = fig1_w_GER,
    width = 9,
    height = 6,
    dpi = 300,
    bg = "white",
    device = "svg"
  )
}



# ==================================================================================================================================================================
# Confidence Intervals via Bootstrapping
# ==================================================================================================================================================================

# Share of DUIC-only and DUIC-PSU episodes among all DUIC episodes
boot_fun <- function(data, indices) {
  d <- data[indices, ]
  total <- sum(d$DUIC_freq_30_num, na.rm = TRUE)
  only <- sum(d$DUIC_only_freq_30_num, na.rm = TRUE)
  psu <- sum(d$DUIC_psu_freq_30_num, na.rm = TRUE)
  return(c(only / total, psu / total))
}


# include only cases with full and plausible data
boot_data <- DUIC30_df %>% 
  filter(analyseDUIC == 1)

# bootstrap
set.seed(123)
boot_results <- boot(data = boot_data, statistic = boot_fun, R = 2000)

hist(boot_results$t[, 1],
     breaks = 50,
     main = "Bootstrap-Verteilung: Anteil DUIC-only an allen DUIC",
     xlab = "Anteil",
     col = "skyblue",
     border = "white")

hist(boot_results$t[, 2],
     breaks = 50,
     main = "Bootstrap-Verteilung: Anteil DUIC-PSU an allen DUIC",
     xlab = "Anteil",
     col = "skyblue",
     border = "white")

# CIs 
boot.ci(boot_results, index = 1, type = "bca") #only/total
boot.ci(boot_results, index = 2, type = "bca") #psu/total


# Share of DUIC episodes (PSU and only) by cannabis use frequency
boot_fun_DUIC_canfreq_share <- function(data, indices) {
  d <- data[indices, ]
  
  result <- d %>%
    group_by(can_freq) %>%
    summarise(
      n_users = n(),
      n_duicep = sum(DUIC_freq_30_num[analyseDUIC == 1], na.rm = TRUE), 
      n_duicep_only = sum(DUIC_only_freq_30_num[analyseDUIC == 1], na.rm = TRUE),
      n_duic_psu = sum(DUIC_psu_freq_30_num[analyseDUIC == 1], na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(result) < length(unique(data$can_freq))) {
    return(rep(NA, 4 * length(unique(data$can_freq))))
  }
  
  result <- result %>%
    arrange(can_freq) %>%
    mutate(
      share_users = n_users / sum(n_users),
      share_duicep = n_duicep / sum(n_duicep),
      share_duicep_only = n_duicep_only / sum(n_duicep_only),
      share_duic_psu = n_duic_psu / sum(n_duic_psu)
    )
  
  vec <- c(
    result$share_users,
    result$share_duicep,
    result$share_duicep_only,
    result$share_duic_psu
  )
  
  var_names <- c("share_users", "share_duicep", "share_duicep_only", "share_duic_psu")
  group_names <- result$can_freq
  
  names(vec) <- paste0(rep(var_names, each = length(group_names)), ".", group_names)
  
  return(vec)
}

boot_result_canfreq <- boot(data = DUIC30_df, statistic = boot_fun_DUIC_canfreq_share, R = 2000)
names(boot_result_canfreq$t0)

idx <- which(names(boot_result_canfreq$t0) == "share_duic_psu.(fast) täglich") #insert name of DUICtype+canfreq here
hist(boot_result_canfreq$t[, idx],
     breaks = 50,
     main = "Bootstrap-Verteilung",
     xlab = "Anteil",
     col = "tomato",
     border = "white")
boot.ci(boot_result_canfreq, index = idx, type = "bca")

