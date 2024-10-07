###############################################################################-
#                         Thesis: Chapter 3 Analysis                           #
#  Chapter 4: Wastewater Surveillance for AMR: A Longitudinal Study across     #
#                 Urban, Rural, and Remote Communities                         #
#                         Jade Daigle | 2024                                  #
###############################################################################-

# Load required libraries
library(tidyverse)
library(vegan)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(readxl)
library(lme4)
library(skimr)
library(lmerTest)
library(DHARMa)
library(pairwiseAdonis) 
library(xlsx)

#### Alpha Diversity Analysis ####

# Import data
qdata <- read_xlsx("chapter3_data.xlsx", sheet = "qmetaquant")

# Summary of imported data
skim(qdata)

# Catchment data filtering
catchment <- read_xlsx("chapter3_data.xlsx", sheet = "catchment") |> 
  dplyr::filter(location != "NTC")

# Join datasets
qdata <- qdata |> 
  mutate(location = str_sub(name, 1, -3)) |> 
  full_join(catchment, by = "location")

# Make data wide for diversity calculations
adiv <- qdata |> 
  dplyr::select(c(name, target, abs_abund.rrna)) |> 
  pivot_wider(names_from = "target", values_from = "abs_abund.rrna", values_fill = 0) |> 
  remove_rownames() |> 
  column_to_rownames(var = 'name')

# Create dataframe for alpha diversity measurements
alpha_div <- qdata |> 
  dplyr::select(c(name, location, catchment)) |> 
  unique() |> 
  mutate(week = str_sub(name, -1))

# Convert to factors
alpha_div$location <- factor(alpha_div$location)
alpha_div$week <- factor(alpha_div$week)

# Calculate Shannon and Simpson diversity indices
alpha_div$a_shannon <- diversity(adiv, index = 'shannon')
alpha_div$a_simpson <- diversity(adiv, index = 'invsimpson')

# Set levels for location factor
alpha_div$location <- factor(alpha_div$location, levels = c('Remote', 'Rural', 'URB-1A', 'URB-2', 'URB-1C', 'URB-1B'))

# Statistical tests for diversity indices
kruskal.test(a_shannon ~ location, data = alpha_div)
pairwise.wilcox.test(alpha_div$a_shannon, alpha_div$location, p.adjust.method = 'BH')

kruskal.test(a_simpson ~ location, data = alpha_div)
pairwise.wilcox.test(alpha_div$a_simpson, alpha_div$location, p.adjust.method = 'BH')

#### Beta Diversity Analysis ####

# NMDS ordination using Bray-Curtis distance
set.seed(2525)
nmds.arg <- metaMDS(adiv, distance = "bray", k = 3, trymax = 100)
print(paste("Stress Score:", nmds.arg$stress))

# Extract NMDS site scores and convert to dataframe
nmds.a <- as.data.frame(nmds.arg$points) |> 
  rownames_to_column(var = "name") |> 
  mutate(location = str_sub(name, 1, 3))

# Bray-Curtis distance matrix
set.seed(2525)
group_bray <- vegdist(adiv, method = "bray")

# Create metadata from rownames and merge with NMDS scores
metadata <- data.frame(name = rownames(adiv), location = str_sub(rownames(adiv), 1, -3))
nmds.a <- merge(nmds.a, metadata, by = "name")

# PERMANOVA analysis
adonis2_result <- adonis2(group_bray ~ location, data = metadata)
print(adonis2_result)

# Pairwise comparisons
pairwise_results <- pairwise.adonis(group_bray, metadata$location)
print(pairwise_results)

#### Data Modeling ####

# Load and clean data
qdata <- read_xlsx("chapter3_data.xlsx", sheet = "qmetaquant")
catchment <- read_xlsx("chapter3_data.xlsx", sheet = "catchment") |> 
  dplyr::filter(location != "NTC")
qdata <- qdata |> 
  mutate(location = str_sub(name, 1, -3)) |> 
  full_join(catchment, by = "location") |> 
  dplyr::select(c(name, location, date_sampled, target, abs_abund.rrna)) |>
  dplyr::rename("aro_name" = "target")

# CARD database
card <- read_xlsx("chapter3_data.xlsx", sheet = "card_short")

# Metadata for sequencing summary
metadata <- read_xlsx("chapter3_data.xlsx", sheet = "sequencing_summary") |>
  filter(location != "NTC") |>
  dplyr::select(c(name, days_start)) |>
  dplyr::filter(name != "Remote_4") |> 
  mutate(days_start = as.numeric(days_start))

# Merge datasets and prepare for modeling
data <- qdata |>
  full_join(card, by = "aro_name") |>
  drop_na(name) |>
  dplyr::select(name, location, new_gene_family, abs_abund.rrna) |>
  group_by(name, location, new_gene_family) |>
  dplyr::summarise(total = sum(abs_abund.rrna)) |>
  ungroup() |> 
  full_join(metadata, by = "name") |>
  mutate(days_start_scaled = scale(days_start),
         week = str_sub(name, -1)) |> 
  filter(!is.na(new_gene_family))

#### Linear Mixed-Effects Model ####

# Filter and model the data
model_data <- filtered_df

# Fit the linear mixed-effects model
jmodel <- lmer(log10(total) ~ days_start_scaled + (1 + days_start_scaled | location), data = filtered_df)

# Model summary and diagnostics
summary(jmodel)
sim_res_mixed <- simulateResiduals(jmodel)
plot(sim_res_mixed)
AIC(jmodel)
BIC(jmodel)

# Extract and process random effects
random_effects <- ranef(jmodel)
data <- as.data.frame(random_effects)
data$location <- rownames(random_effects[[1]])
