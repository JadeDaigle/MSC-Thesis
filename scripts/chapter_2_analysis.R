###############################################################################-
#                         Thesis: Chapter 2 Analysis                           #    
#     Optimizing a Wastewater DNA Extraction for AMR Detection                 #
#                         Jade Daigle | 2024-03-12                             #
###############################################################################-

# Load required libraries
library(readxl)
library(tidyverse)
library(Rmisc)
library(ggpmisc)
library(ggpubr)

# Load metadata
metadata <- read_xlsx("chapter2_data.xlsx", sheet = "metadata")

# Common theme for all plots
common_theme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 12),
  strip.text = element_text(size = 14)  
)

#### Figure 1: Pellet Mass and DNA Concentration ####

#---- Figure 1A: Pellet mass vs. Process volume
df1 <- metadata |> 
  summarySE(measurevar = "pellet_size", groupvars = "vol_ww")

p1 <- df1 |> 
  ggplot(aes(x = vol_ww, y = pellet_size)) + 
  geom_line(col = "black", size = 1) +
  geom_errorbar(aes(ymin = pellet_size - se, ymax = pellet_size + se), width = 2, size = 0.5) +
  geom_point(size = 3, col = "black") +
  scale_x_continuous(breaks = c(25, 50, 100, 200)) +
  stat_poly_eq(aes(label = after_stat(eq.label)), size = 4) + 
  stat_poly_eq(label.y = 0.9, size = 4) +
  labs(x = "Process Volume (mL)", y = "Pellet Mass (g)") +
  theme_bw() +
  common_theme

#---- Figure 1B: Qubit DNA concentration vs. Process volume
df2 <- metadata |>
  dplyr::select("ex_kit", "vol_ww", "qubit_DNA") |>
  summarySE(measurevar = "qubit_DNA", groupvars = c("ex_kit", "vol_ww"))

color <- c("#264653", "#f4a261")

p2 <- df2 |>
  ggplot(aes(x = vol_ww, y = qubit_DNA, group = ex_kit, colour = ex_kit)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = qubit_DNA - se, ymax = qubit_DNA + se), width = 2, size = 0.5) +
  geom_point(size = 3.5) +
  scale_color_manual(values = color) +
  scale_x_continuous(breaks = c(25, 50, 100, 200)) +
  coord_cartesian(ylim = c(0, 400)) +
  facet_wrap(~ ex_kit) +
  stat_poly_eq(aes(label = after_stat(eq.label)), size = 4, label.y = 0.95) +
  stat_poly_eq(label.y = 0.85, size = 4) +
  labs(x = "Process Volume (mL)", y = "DNA Concentration (ng/µL)", colour = "Extraction Kit") +
  theme_bw() +
  common_theme +
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

#---- Figure 1C: Nanodrop vs. Process volume
df3 <- metadata |>
  dplyr::select(c(ex_kit, vol_ww, nanodrop_260_280, nanodrop_260_230)) |>
  dplyr::rename(`A260:280` = nanodrop_260_280, `A260:230` = nanodrop_260_230) |>
  pivot_longer(cols = c(3:4), names_to = "variable", values_to = "value") |>
  summarySE(measurevar = "value", groupvars = c("ex_kit", "vol_ww", "variable"))

p3 <- df3 |>
  ggplot(aes(x = vol_ww, y = value, group = ex_kit, colour = ex_kit)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 2, size = 0.5) +
  geom_point(size = 3.5) +
  geom_hline(data = df3 |> filter(variable == "A260:280"), aes(yintercept = 1.7), col = "black", linetype = 5) +
  geom_hline(data = df3 |> filter(variable == "A260:230"), aes(yintercept = 2.0), col = "black", linetype = 5) +
  scale_x_continuous(breaks = c(25, 50, 100, 200)) +
  scale_color_manual(values = color) +
  facet_wrap(~ variable) +
  labs(x = "Process Volume (mL)", y = "Absorbance ratio") +
  theme_bw() +
  common_theme +
  theme(legend.position = "none")

# Arrange the plots together
combined_plot <- ggarrange(
  ggarrange(p2, ncol = 1, labels = c("A"), common.legend = TRUE),
  p3,
  ncol = 1,
  heights = c(1, 1),
  labels = c("", "B"),
  hjust = -0.5,
  vjust = 1
)

# Print and save the combined plot
print(combined_plot)
ggsave("chapter2_1.png", combined_plot, width = 10, height = 8)

#### Figure 2: qPCR Analysis ####

# Prepare qPCR data for Figure 2
qpcr_2 <- read_xlsx("chapter2_data.xlsx", sheet = "cp") |>
  dplyr::select(-linelist) |>
  pivot_longer(names_to = "target", values_to = "cp", cols = -1) |>
  full_join(metadata, by = "sample_id") |>
  dplyr::select(c(ex_kit, vol_ww, target, cp)) |>
  dplyr::filter(!grepl("10", target)) |>
  mutate(log_cp = log10(cp),
         target = str_sub(target, 1, -6),
         target = recode(target,
                         "ecoli" = "E. coli",
                         "entero" = "Entero. spp.",
                         "cpq" = "CrAssphage",
                         "kpc" = "KPC",
                         "ndm" = "NDM",
                         "rrna" = "16S rRNA")) |>
  mutate(target = fct_reorder(target, log_cp, .fun = max, .desc = TRUE, .na_rm = TRUE)) |>
  summarySE(measurevar = "log_cp", groupvars = c("ex_kit", "target", "vol_ww"))

# Plot qPCR results
p4 <- qpcr_2 |>
  ggplot(aes(x = vol_ww, y = log_cp, color = ex_kit)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = log_cp - se, ymax = log_cp + se), width = 2, size = 0.5) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(25, 50, 100, 200)) +
  scale_color_manual(values = color) +
  facet_grid(target ~ ex_kit, scales = "free") +
  labs(x = "Process Volume (mL)", y = "Concentration (copies/mL)\n", col = "Extraction Kit") +
  theme_bw() +
  common_theme +
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

# Print and save the qPCR plot
print(p4)
ggsave("chapter2_2.png", p4, width = 10, height = 8)

#### Figure 3: PMB Extraction Analysis ####

# Prepare qPCR data for Figure 3
qpcr_3 <- read_xlsx("chapter2_data.xlsx", sheet = "cp") |>
  dplyr::select(-linelist) |>
  pivot_longer(names_to = "id", values_to = "cp", cols = -1) |>
  full_join(metadata, by = "sample_id") |>
  dplyr::filter(ex_kit == "PMB") |>
  dplyr::select(c(vol_ww, id, cp)) |>
  mutate(target = str_sub(id, 1, -6),
         dilution = case_when(
           grepl("100$", id) ~ "100-fold",
           grepl("10$", id) ~ "10-fold",
           TRUE ~ "Neat"),
         target = case_when(
           grepl("entero", target) ~ "Entero spp.",
           grepl("ecoli", target) ~ "E. coli",
           grepl("cpq", target) ~ "crAssphage",
           grepl("rrna", target) ~ "16S rRNA",
           grepl("ndm", target) ~ "NDM",
           grepl("kpc", target) ~ "KPC"),
         log_cp = log10(cp)) |>
  mutate(target = fct_reorder(target, log_cp, .fun = max, .desc = TRUE, .na_rm = TRUE)) |>
  summarySE(measurevar = "log_cp", groupvars = c("target", "dilution", "vol_ww"))

# Plot qPCR results for PMB extraction
color_2 <- c('#264653', '#2a9d8f', '#f4a261')

p5 <- qpcr_3 |>
  ggplot(aes(x = vol_ww, y = log_cp, color = dilution)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = log_cp - se, ymax = log_cp + se), width = 2, size = 0.5) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(25, 50, 100, 200)) +
  scale_color_manual(values = color_2) +
  facet_wrap(~ target) +
  labs(x = "Process Volume (mL)", y = "Concentration (Log10 Copies/mL)\n", col = "Dilution") +
  theme_bw() +
  common_theme +
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

# Print and save the plot
print(p5)
ggsave("chapter2_3.png", p5, width = 10, height = 8)

#### Figure 4: Inhibition Analysis ####

# Prepare inhibition data for Figure 4
data <- read_xlsx("chapter2_data.xlsx", sheet = "inhibition") |>
  dplyr::select(c(1, 7, 24:27)) |>
  mutate('1:10' = .data[['10_avg']] / .data[['neat_avg']],
         '1:100' = .data[['100_avg']] / .data[['neat_avg']],
         pr1000 = .data[['1000_avg']] / .data[['neat_avg']],
         pr10_10 = .data[['10_avg']] / .data[['10_avg']],
         pr100_10 = .data[['100_avg']] / .data[['10_avg']],
         pr100_neat = .data[['100_avg']] / .data[['neat_avg']],
         pr1000_100 = .data[['1000_avg']] / .data[['100_avg']]) |>
  pivot_longer(names_to = 'variable', values_to = 'value', cols = c(7:13)) |>
  dplyr::filter(variable == "pr1000_100" | variable == "pr100_10" | variable == "1:10") |>
  mutate(target = recode(target, 
                         "E.coli" = "E. coli",
                         "Entero" = "Entero spp.")) |>
  mutate(variable = recode(variable,
                           "pr1000_100" = "1:100",
                           "pr100_10" = "1:10",
                           "1:10" = "Neat"))

# Set factor levels
data$variable <- factor(data$variable, levels = c("Neat", "1:10", "1:100"))
data$target <- factor(data$target, levels = c('16S rRNA', 'Entero spp.', 'CPQ', 'E. coli', 'KPC'))

# Plot inhibition analysis
p6 <- data |>
  ggplot(aes(x = target, y = value, fill = variable)) +    
  geom_bar(stat = 'identity', position = 'dodge', width = 0.5) +
  coord_cartesian(ylim = c(0.25, 4)) +
  scale_y_continuous(breaks = c(1, 2, 3, 4)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  facet_wrap(~ sample_id) +
  labs(x = "Target", y = "Fold Difference", fill = 'Dilutions') +
  scale_fill_manual(values = c('#264653', '#2a9d8f', '#f4a261')) +
  theme_bw() +
  common_theme +
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

# Print and save the inhibition plot
print(p6)
ggsave("Chapter2_4.png", p6, width = 10, height = 6)

                 
