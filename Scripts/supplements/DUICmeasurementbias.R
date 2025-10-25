# ==================================================================================================================================================================
#==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  CANNASTREET - Analysis of Manuscript "Short-term effects of cannabis legalisation in Germany on driving under the influence of cannabis"
# SUPPLEMENT: MEASUREMENT BIAS CHECKS (only GER data as AT data is limited to few variables)
# CODE AUTHOR:    Anna Schranz
# DATE STARTED:   251025

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================


# 0) ESSENTIALS
# ______________________________________________________________________________________________________________________

# clean workspace
rm(list=ls())

packages <- c("tidyverse", "ggthemes", "haven")

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


## DEFINITIONS
# ______________________________________________________________________________________________________________________
#colors
colors_country <- c("AT" = "#EF3340", "DE" = "#000000")

#labels
country_labels <- c(
  "AT" = "Austria",
  "DE" = "Germany"
)

# DATA IMPORT
# ______________________________________________________________________________________________________________________
path_GER_1 <- file.path("Data/", "cleaned/", "full_GER_wave1_cleaned2.rds")
path_GER_2 <- file.path("Data/", "cleaned/", "full_GER_wave2_cleaned2.rds")

#read in data
data_GER_1 <- readRDS(path_GER_1)
data_GER_2 <- readRDS(path_GER_2)


# DATA PREPARATION
# ______________________________________________________________________________________________________________________
#select relevant variables
data_GER_1_sel <- data_GER_1 %>%
  select(-weights) %>%
  mutate(DUIC.MIXEDUSE_old = DUIC.MIXEDUSE, weights = weights_ohneWiederbefragt, Land = "DE") %>%
  dplyr::select(ID, sex, weights, wiederbefragt, edu_group, agegroup, agegroup_basicsample, age_in_years,
                GSZB2, ZS, AS1, AS2, AS3, can_freq, can_use_12M, DRIVERLICENSE, DUIC, DUIC.FREQ.01M_wave1coding, DUIC.FREQ.12M_wave1coding,
                DUIC12m_full, DUIC12m_alone, DUIC12m_mixed, DUIC30d_full, DUIC30d_alone, DUIC30d_mixed, Welle, 
                MEDICALUSE.01, MEDICALUSE.02, DUIC.MIXEDUSE_old, Land, DUIC.INTOX, DGURBA, sport_freq, sport_12M, fish_12M, gisd_5, AGS,
                gambling_freq, alcohol_freq, tobacco_freq, DUIA, KSE_NQ_z, KSE_PQ_z)

data_GER_2_sel <- data_GER_2 %>%
  select(-weights) %>%
  mutate(weights = weights_ohneWiederbefragt, Land = "DE") %>% 
  select(ID, sex, weights, wiederbefragt, edu_group,  age_in_years, 
         agegroup, agegroup_basicsample, GSZB2, ZS, AS1, AS2, AS3,
         can_freq, can_use_12M, DUIC.MIXEDUSE, DUIC.MED, DRIVERLICENSE, CWM.PRE, DUIC, DUIC.FREQ.01M_wave1coding, DUIC.FREQ.01M_exact, DUIC.FREQ.12M_wave1coding,
         DUIC.FREQ.01M_wave2coding, DUIC.FREQ.12M_wave2coding, DUIC.FREQ.12M_exact, DUIC.MED, DUIC12m_full, DUIC12m_alone, DUIC12m_mixed, 
         DUIC30d_full, DUIC30d_alone, DUIC30d_mixed, Welle, MEDICALUSE.01, MEDICALUSE.02, Land, DUIC.INTOX, DGURBA, sport_12M, sport_freq, fish_12M, gisd_5, AGS,
         gambling_freq, alcohol_freq, tobacco_freq, DUIA, use_benzos, KSE_NQ_z, KSE_PQ_z)

#combine data
data <- bind_rows(data_GER_1_sel, data_GER_2_sel)
data$Welle <- factor(data$Welle, levels = c(1, 2))

# Sample 2
df_DUIC_full <- subset(data, AS3 == 1)

# FUNCTIONS
# ______________________________________________________________________________________________________________________
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

# MEASUREMENT BIAS CHECKS
# ______________________________________________________________________________________________________________________

#Neben den Korrelationen mit alternativen Maßen für soziale Erwünschtheit, stützen auch die Korrelationen mit den Devianzmaßen die Gültigkeit der KSE-G.
#Diese zeigen, dass Personen, die Gesetzesübertretungen wie Fahren ohne Fahrkarte oder Diebstahl zugeben,
#eher weniger dazu neigen ihre positiven Qualitäten zu übertreiben bzw. negative zu untertreiben.
#Sowohl für PQ+ als auch für NQ- finden sich geringe bis mittlere negative Zusammenhänge mit den Devianzmaßen.
#(https://zis.gesis.org/skala/Kemper-Beierlein-Bensch-Kovaleva-Rammstedt-Soziale-Erw%C3%BCnschtheit-Gamma-%28KSE-G%29?lang=de#Tabelle6a)



#correlation of DUIC and social desirability KSE_NQ_z and KSE_PQ_z before legalization (wave 1)
cor_test_NQ <- cor.test(df_DUIC_full[df_DUIC_full$Welle == 1,]$DUIC12m_full, df_DUIC_full[df_DUIC_full$Welle == 1,]$KSE_NQ_z, method = "pearson") #siginicant negative correlation -> DUIC is less likely with higher denial of negative qualities
cor_test_PQ <- cor.test(df_DUIC_full[df_DUIC_full$Welle == 1,]$DUIC12m_full, df_DUIC_full[df_DUIC_full$Welle == 1,]$KSE_PQ_z, method = "pearson") #no significant correlation with exaggeration of positive qualities

# is the influence of social desirability on DUIC different before vs after legalization?
glm_interaction_NQ <- glm(DUIC12m_full ~ KSE_NQ_z * Welle, data = df_DUIC_full, family = binomial(link = "logit"))
summary(glm_interaction_NQ) #no significant interaction -> no difference in influence of social desirability before vs after legalization

#goodness of fit
glm_null <- glm(DUIC12m_full ~ 1, data = df_DUIC_full[!is.na(df_DUIC_full$KSE_NQ_z),], family = binomial)
anova(glm_null, glm_interaction_NQ, test = "Chisq")
logisticPseudoR2s(glm_interaction_NQ)


#sample size per wave (no NAs in KSE_NQ_z)
table(df_DUIC_full$Welle[!is.na(df_DUIC_full$KSE_NQ_z)])




#odds ratios
exp(cbind(OR = coef(glm_interaction_NQ), confint(glm_interaction_NQ)))