# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

# Read coverage data from samtools coverage output
coverage_data <- read.table("reads_coverage_summary.txt", header = TRUE, sep = "\t")

# Display the data structure
print("Coverage data:")
print(coverage_data)

# Prepare data for plotting
plot_data <- coverage_data %>%
  mutate(
    Type = case_when(
      grepl("NC_003197", X.rname) ~ "Chromosome",
      grepl("NC_003277", X.rname) ~ "Plasmid",
      TRUE ~ "Other"
    ),
    Length_Mb = endpos / 1e6,
    Coverage_Pct = coverage,
    Mean_Depth = meandepth,
    Base_Quality = meanbaseq,
    Map_Quality = meanmapq
  ) %>%
  filter(Type %in% c("Chromosome", "Plasmid"))

# Create individual plots
p1 <- ggplot(plot_data, aes(x = Type, y = Coverage_Pct, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = paste0(round(Coverage_Pct, 1), "%")), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Chromosome" = "#2E86AB", "Plasmid" = "#A23B72")) +
  labs(title = "Genome Coverage", y = "Coverage (%)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, 110)

p2 <- ggplot(plot_data, aes(x = Type, y = Mean_Depth, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = paste0(round(Mean_Depth, 1), "x")), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Chromosome" = "#2E86AB", "Plasmid" = "#A23B72")) +
  labs(title = "Sequencing Depth", y = "Mean Depth (x)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  )

p3 <- ggplot(plot_data, aes(x = Type, y = numreads, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = comma(numreads)), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Chromosome" = "#2E86AB", "Plasmid" = "#A23B72")) +
  scale_y_continuous(labels = comma) +
  labs(title = "Aligned Reads", y = "Number of Reads") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  )

p4 <- ggplot(plot_data, aes(x = Type, y = Base_Quality, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = round(Base_Quality, 1)), 
            vjust = -0.5, size = 5, fontface = "bold") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "#E63946", size = 1) +
  annotate("text", x = 1.5, y = 32, label = "Q30 threshold", 
           color = "#E63946", size = 3.5) +
  scale_fill_manual(values = c("Chromosome" = "#2E86AB", "Plasmid" = "#A23B72")) +
  labs(title = "Mean Base Quality", y = "Quality Score") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, max(plot_data$Base_Quality) + 5)

p5 <- ggplot(plot_data, aes(x = Type, y = Map_Quality, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = round(Map_Quality, 1)), 
            vjust = -0.5, size = 5, fontface = "bold") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#06A77D", size = 1) +
  annotate("text", x = 1.5, y = 52, label = "High quality", 
           color = "#06A77D", size = 3.5) +
  scale_fill_manual(values = c("Chromosome" = "#2E86AB", "Plasmid" = "#A23B72")) +
  labs(title = "Mean Mapping Quality", y = "MAPQ Score") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, max(plot_data$Map_Quality) + 10)

p6 <- ggplot(plot_data, aes(x = Type, y = Length_Mb, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = paste0(round(Length_Mb, 2), " Mb")), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Chromosome" = "#2E86AB", "Plasmid" = "#A23B72")) +
  labs(title = "Contig Size", y = "Length (Mb)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  )

# Create summary table
summary_table <- plot_data %>%
  select(Type, Contig_ID = X.rname, Length_Mb, Reads = numreads, 
         Coverage = Coverage_Pct, Mean_Depth, Base_Q = Base_Quality, 
         Map_Q = Map_Quality) %>%
  mutate(
    Length_Mb = round(Length_Mb, 2),
    Coverage = paste0(round(Coverage, 2), "%"),
    Mean_Depth = paste0(round(Mean_Depth, 1), "x"),
    Base_Q = round(Base_Q, 1),
    Map_Q = round(Map_Q, 1)
  )

# Create table plot
table_theme <- gridExtra::ttheme_default(
  core = list(
    fg_params = list(fontsize = 11),
    bg_params = list(fill = c("#2E86AB", "#A23B72"), alpha = 0.3)
  ),
  colhead = list(
    fg_params = list(fontsize = 12, fontface = "bold"),
    bg_params = list(fill = "#457B9D", alpha = 0.7)
  ),
  rowhead = list(fg_params = list(fontsize = 11))
)

table_grob <- gridExtra::tableGrob(summary_table, rows = NULL, theme = table_theme)

# Combine all plots
title_grob <- grid::textGrob(
  "Salmonella Genome Coverage Analysis\nRaw Reads vs Reference Genome",
  gp = grid::gpar(fontsize = 18, fontface = "bold")
)

subtitle_grob <- grid::textGrob(
  "Coverage Statistics Summary",
  gp = grid::gpar(fontsize = 14, fontface = "bold")
)

# Arrange plots
top_plots <- grid.arrange(p1, p2, p3, ncol = 3)
bottom_plots <- grid.arrange(p4, p5, p6, ncol = 3)

# Final combined plot
final_plot <- grid.arrange(
  title_grob,
  top_plots,
  bottom_plots,
  subtitle_grob,
  table_grob,
  ncol = 1,
  heights = c(0.8, 3, 3, 0.5, 1.5)
)

# Save the plot
ggsave("coverage_analysis.png", final_plot, width = 14, height = 10, dpi = 300)

cat("Coverage analysis plot saved as 'coverage_analysis.png'\n")