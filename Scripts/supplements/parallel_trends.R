# ==================================================================================================================================================================
#==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  CANNASTREET - Analysis of Manuscript "Short-term effects of cannabis legalisation in Germany on driving under the influence of cannabis"
# SUPPLEMENT: CHECK OF PARALLEL TRENDS ASSUMPTION 
# CODE AUTHOR:    Anna Schranz
# DATE STARTED:   020625

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================


# 0) ESSENTIALS
# ______________________________________________________________________________________________________________________

# clean workspace
rm(list=ls())

packages <- c("tidyverse", "ggthemes")

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

# DATA
# ______________________________________________________________________________________________________________________

#CANNABIS USE 
# Adults (AT: 15–64 Jahre (EUDA/GöG), GER: 18–64 Jahre (ESA))
adults_data <- data.frame(
  Land = c(rep("AT", 3), rep("DE", 6)),
  Jahr = c(2008, 2015, 2020, 2009, 2012, 2015, 2018, 2021, 2024),
  can_use_prev12M = c(3.5, 6.4, 6.3, 4.8, 4.5, 6.1, 7.1, 8.8, 9.8),
  n = c(3000, 3000, 4650, 8030, 9061, 9204, 9101, 8986, NA),
  #dummy values for n in AT that are smaller than probabable n -> increases variance, GER: actual sample sizes
  cohort = "Adults"
)

# Daten für Jugendliche (ESPAD-Daten)
adolescents_data <- data.frame(
  Land = c("AT", "DE", "AT", "DE", "AT", "DE", "AT", "DE", "AT", "DE"),
  Jahr = c(2003, 2003, 2007, 2007, 2015, 2015, 2019, 2019, 2024, 2024),
  can_use_prev12M = c(9.3, 10.9, 6.2, 5.7, 9.2, 7.7, 11.3, 10.5, 6.1, 6.8),
  cohort = "Adolescents"
)

# MOTOR VEHICLE CRASHES

mva_filepath <- "Data/MVA/"
mva_GER <- readRDS(paste0(mva_filepath, "destatis_verkehr_cleaned.rds"))

mva_AT <- data.frame(
  Land = "AT",
  Jahr = c(2020, 2021, 2022, 2023, 2024),
  value = c(106, 101, 114, 152, 176)
)


# POPULATION DATA
pop_AT <- data.frame(
  Land = "AT",
  Jahr = c(2020, 2021, 2022, 2023, 2024),
  pop_value = c(8929910, 8932664, 8978929, 9104772, 9158750)
)

pop_DE <- data.frame(
  Land = "DE",
  Jahr = c(2019, 2020, 2021, 2022, 2023, 2024),
  pop_value = c(83019213, 83166711, 83155031, 83237124 , 83118501, 83456045)
)

pop <- rbind(pop_AT, pop_DE)

# PREPARE data
# ______________________________________________________________________________________________________________________
#yearly instead of monthly MVA data
mva_GER_yearly <- mva_GER %>%
  group_by(year) %>%
  filter(outcome == "secondaryC") %>%
  summarise(value = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  rename(Jahr = year) %>%
  mutate(Land = "DE") 

# bind data
mva_data <- rbind(mva_GER_yearly, mva_AT) %>%
  left_join(pop, by = c("Land", "Jahr")) %>%
  #calculate per 100,000 inhabitants
  mutate(value_per100k = (value / pop_value) * 100000)

cannabis_use_data <- rbind(adults_data, adolescents_data) 


# PLOTTING
# ______________________________________________________________________________________________________________________

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
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


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

trend_MVA_pre <- ggplot(mva_data, aes(x = Jahr, y = value_per100k, color = Land)) +
  geom_line(size = 0.7) +
  geom_point(size = 1.5) +
  labs(x = "Year",
       y = "Number of motor vehicle crashes involving\ndrugs (per 100,000 inhabitants)",
       color = ""
  ) +
  scale_y_continuous(
    breaks = seq(0, 3, by = 0.5),
    limits = c(0, 3), expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = seq(2020, 2025, by = 1),
    limits = c(2020, 2024)
  ) +
  #vertical line for t0 measurement Jahr = 2023
  geom_vline(xintercept = 2023, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = colors_country, labels = country_labels) +
  theme_gdocs(base_family = "Calibri", base_size = 9) +
  theme(
    panel.grid.minor.y = element_line(color = "gray80", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(color= "black"),
    axis.title.x = element_text(color= "black"),
    axis.text.y = element_text(color= "black"),
    axis.text.x = element_text(color= "black", angle = 45, hjust = 1),
    strip.text = element_text(color= "black", face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(color= "black"),
    legend.text = element_text(color= "black"),
    legend.box.margin = margin(t = -9, unit = "pt") 
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

if (dataexport) {
  ggsave(
    filename = paste0(folder_path_plots, "trend_MVA_pre_", DATE, ".tiff"),
    plot = trend_MVA_pre,
    width = 107,                   
    height = 90,                  
    units = "mm",                  
    dpi = 400,
    device = "tiff"
  )
}


# STATISTICAL TESTING of parallel trends
# ______________________________________________________________________________________________________________________

# MVAs
# poisson model to test for differences in trends pre-legalisation
mva_data_mod <- mva_data %>%
  filter(Jahr < 2024 & Jahr >= 2020) %>%
  mutate(Land = factor(Land, levels = c("DE", "AT")),
         Jahr_c = Jahr - min(Jahr)) #center year variable
         
mva_poisson <- glm(value ~ Jahr_c * Land + offset(log(pop_value)),
                            family = poisson(link = "log"),
                            data = mva_data_mod)
         
summary(mva_poisson)
# dispersion test
sum(residuals(mva_poisson, type="pearson")^2) / mva_poisson$df.residual # should be close to 1 for poisson (1.18 - acceptable)

# cannabis use 
adults_data_mod <- adults_data %>%
  filter(!is.na(n)) %>%  # 2024 entfernen
  mutate(
    can_use = round(can_use_prev12M / 100 * n),  # estimated number of cannabis users
    Jahr_c = Jahr - min(Jahr)  # center year variable
  )

can_use_mod <- glm(cbind(can_use, n - can_use) ~ Jahr_c * Land,
                     family = binomial(link = "logit"),
                     data = adults_data_mod)
summary(can_use_mod)