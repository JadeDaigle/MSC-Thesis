###############################################################################-
#                         Thesis: Chapter 3 Analysis                           #
#       Development and Validation of a Wastewater-Specific Quantitative       #
#           Metagenomic Workflow that Quantifies AMR gene families             #
#                         Jade Daigle | 2024-03-12                             #
###############################################################################-

#### Figure 12: Flow and Bacteria Levels ####

# Load required libraries
library(tidyverse)
library(ggpubr)
library(patchwork)
library(readxl)
library(grid)

# Import metadata and preprocess
metadata <- read_xlsx("chapter3_data.xlsx", sheet = "sequencing_summary") |> 
  mutate(location = str_sub(name, 1, -3)) |> 
  filter(location != "NTC") |> 
  filter(name != "Remote_4") |> 
  mutate(flow = as.numeric(flow))

# Define the color palette for the plots
colors <- c("Flow (Megaliters/day)" = "#E69F00", "Total Bacteria (Log10 Copies/mL)" = "#183E35")

# Function to create plots for each location
create_plots <- function(location_name, y1_lim, y2_lim, show_x_axis_p2 = FALSE) {
  data <- metadata |> 
    mutate(`Total Bacteria (Log10 Copies/mL)` = log10(rrna_avg),
           `Flow (Megaliters/day)` = flow) |> 
    dplyr::select(date_sampled, location, `Total Bacteria (Log10 Copies/mL)`, `Flow (Megaliters/day)`) |> 
    pivot_longer(cols = c(`Total Bacteria (Log10 Copies/mL)`, `Flow (Megaliters/day)`), names_to = "variable", values_to = "value") |> 
    dplyr::filter(location == location_name)
  
  p1 <- ggplot(data, aes(x = date_sampled, y = value, color = variable, shape = variable, linetype = variable, group = interaction(location, variable))) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    facet_wrap(~ location, scales = "free_y") +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(17, 16)) +  # Triangle for Flow and Circle for Total Bacteria
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_y_continuous(name = "Measurement", breaks = scales::pretty_breaks(n = 3)) +  # Reduce y-axis ticks
    coord_cartesian(ylim = y1_lim) + 
    labs(x = NULL, color = "Variable", shape = "Variable", linetype = "Variable") +
    theme_pubr() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          axis.line.x.bottom = element_blank())
  
  p2 <- ggplot(data, aes(x = date_sampled, y = value, color = variable, shape = variable, linetype = variable, group = interaction(location, variable))) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(17, 16)) +  
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_y_continuous(name = "Measurement", breaks = scales::pretty_breaks(n = 3)) + 
    coord_cartesian(ylim = y2_lim) + 
    labs(x = NULL, color = "Variable", shape = "Variable", linetype = "Variable") +
    theme_pubr() +
    theme(axis.text.x = if (show_x_axis_p2) element_text(angle = 45, hjust = 1) else element_blank(),
          axis.ticks.x = if (show_x_axis_p2) element_line() else element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background = element_blank())
  
  return(list(p1 = p1, p2 = p2, data = data))
}

# Function to add correlation annotations
add_correlation_annotation <- function(plot, data, location_name) {
  data_filtered <- data |> dplyr::filter(location == location_name)
  data_total_bacteria <- data_filtered$value[data_filtered$variable == "Total Bacteria (Log10 Copies/mL)"]
  data_flow <- data_filtered$value[data_filtered$variable == "Flow (Megaliters/day)"]
  
  correlation <- cor.test(data_total_bacteria, data_flow, method = "pearson")
  r_value <- round(correlation$estimate, 2)
  p_value <- format(correlation$p.value, digits = 2, nsmall = 2)
  correlation_label <- paste("R =", r_value, ", p =", p_value)
  
  annotation <- annotation_custom(grob = textGrob(correlation_label, x = unit(0.1, "npc"), y = unit(0.9, "npc"), 
                                                  just = "left", gp = gpar(fontsize = 10, fontface = "italic")))
  
  plot + annotation
}

# Create and annotate plots for each location
urb_1c_plot <- create_plots("URB-1C", c(110, 160), c(8, 10), show_x_axis_p2 = TRUE)
urb_2_plot <- create_plots("URB-2", c(50, 100), c(7.5, 10), show_x_axis_p2 = TRUE)
urb_1b_plot <- create_plots("URB-1B", c(40, 60), c(7.5, 10), show_x_axis_p2 = TRUE)
urb_1a_plot <- create_plots("URB-1A", c(15, 25), c(7.5, 10))
rural_plot <- create_plots("Rural", c(15, 25), c(7.5, 10))
remote_plot <- create_plots("Remote", c(7.5, 10), c(3, 5))

urb_1c_plot$p1 <- add_correlation_annotation(urb_1c_plot$p1, urb_1c_plot$data, "URB-1C")
urb_2_plot$p1 <- add_correlation_annotation(urb_2_plot$p1, urb_2_plot$data, "URB-2")
urb_1b_plot$p1 <- add_correlation_annotation(urb_1b_plot$p1, urb_1b_plot$data, "URB-1B")
urb_1a_plot$p1 <- add_correlation_annotation(urb_1a_plot$p1, urb_1a_plot$data, "URB-1A")
rural_plot$p1 <- add_correlation_annotation(rural_plot$p1, rural_plot$data, "Rural")
remote_plot$p1 <- add_correlation_annotation(remote_plot$p1, remote_plot$data, "Remote")

# Combine the plots into a single figure
combined_plot <- ggarrange(remote_plot$p1 / remote_plot$p2, 
                           rural_plot$p1 / rural_plot$p2, 
                           urb_1a_plot$p1 / urb_1a_plot$p2, 
                           urb_1b_plot$p1 / urb_1b_plot$p2, 
                           urb_2_plot$p1 / urb_2_plot$p2, 
                           urb_1c_plot$p1 / urb_1c_plot$p2, 
                           ncol = 3, nrow = 2, heights = c(1, 1.2), common.legend = TRUE, legend = "top")

# Adjust legend appearance and annotate with axis labels
combined_plot <- combined_plot +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.key.size = unit(1, "cm"))

annotated_plot <- annotate_figure(combined_plot,
                                  bottom = text_grob("Date Sampled", size = 13),
                                  left = text_grob("Measurement", rot = 90, size = 13))

# Print and save the annotated plot
print(annotated_plot)
ggsave("chapter3_1.png", plot = annotated_plot, width = 10, height = 6)

#### Figure 13: qPCR vs wqMeta vs qMeta ####

# Load required libraries
library(tidyverse)
library(skimr)
library(readxl)
library(janitor)
library(xlsx)

# Import qPCR and wqMeta data
qpcr_cpday <- read_xlsx("chapter3_data.xlsx", sheet = "pcr_cpday")
wqmeta <- read_xlsx("chapter3_data.xlsx", sheet = "wqmeta")

# Import metadata for catchment size
catchment <- read_xlsx("chapter3_data.xlsx", sheet = "catchment") |> 
  dplyr::filter(location != "NTC")

# Prepare data for visualization
df <- wqmeta |> 
  dplyr::filter(gene != "Other") |> 
  mutate(location = str_sub(name, 1, -3)) |> 
  dplyr::select(-c(date_sampled, week, aro_term, sub_gene_family, amr_gene_family)) |> 
  group_by(name, location, gene) |> 
  dplyr::summarise(ab_qmeta = sum(qmeta),
                   ab_wqmeta = sum(wqmeta)) |>
  ungroup() |> 
  full_join(qpcr_cpday, by = c("name", "gene")) |> 
  filter(!is.na(ab_wqmeta))

# Merge with catchment data
merged <- df |> 
  dplyr::filter(gene != "MCR") |> 
  full_join(catchment, by = "location")

# Visualize data with log-transformed abundances
plot_df <- merged |>
  mutate(ab_qmeta = log10(ab_qmeta),
         ab_wqmeta = log10(ab_wqmeta),
         cpday = log10(cpday)) |>
  filter(!is.na(gene), gene != "KPC") |> 
  pivot_longer(names_to = "type", values_to = "ab_abund", cols = c(ab_wqmeta, ab_qmeta))

plot_df$gene <- factor(plot_df$gene, levels = c("CTX-M", "CMY", "mecA", "vanA", "OXA"))

# Define a manual size scale for catchment size
plot_df$catchment_size_factor <- cut(plot_df$catchment, 
                                     breaks = c(0, 5185, 51000, 86000, 220000, 404000, Inf),
                                     labels = c("5000", "50000", "100000", "200000", "400000", "400000"))

# Plot comparison between qPCR and wqMeta/qMeta
p2 <- plot_df |> 
  ggplot(aes(x = cpday)) +
  geom_point(data = subset(plot_df, type == "ab_wqmeta"), aes(y = ab_abund, color = "wqMeta", size = catchment_size_factor), alpha = 0.8) +
  geom_point(data = subset(plot_df, type == "ab_qmeta"), aes(y = ab_abund, color = "qMeta", size = catchment_size_factor), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(12, 21)) +
  facet_wrap(~ gene) +
  stat_cor(data = subset(plot_df, type == "ab_wqmeta"), aes(x = cpday, y = ab_abund), label.y = 20.5, size = 4.5, method = "pearson", color = "#264653") +
  stat_cor(data = subset(plot_df, type == "ab_qmeta"), aes(x = cpday, y = ab_abund), label.y = 19.5, size = 4.5, method = "pearson", color = "#e76f51") +
  labs(x = "\nqPCR Absolute Abundance (cp/day)",
       y = "Quantitative Metagenomics Absolute Abundance (cp/day)\n",
       size = "Catchment Size") +
  scale_color_manual(values = c("wqMeta" = "#264653", "qMeta" = "#e76f51")) +
  scale_size_manual(values = c("5000" = 2, "50000" = 3, "100000" = 3.5, "200000" = 4, "400000" = 6)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 4), title = NULL),
         size = guide_legend(order = 2, title = "Catchment Size"))

# Print and save the comparison plot
print(p2)
ggsave("chapter3_2.png", plot = p2, width = 10, height = 6)
