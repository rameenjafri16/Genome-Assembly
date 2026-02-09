library(tidyverse)
library(patchwork)

plasmid_id <- "NC_003277.2"
BIN <- 100L   # variant bins

# --------- Load variants ----------
variants <- read.table("variants_hq.bed", header=FALSE)
colnames(variants) <- c("chr", "start", "end", "type", "qual")
variants <- as_tibble(variants)

# --------- Load NEW fine-window coverage ----------
# Expected columns: chr, start, end, depth
coverage_small <- read.table("coverage_100bp.tsv", header=FALSE)
colnames(coverage_small) <- c("chr", "start", "end", "depth")
coverage_small <- as_tibble(coverage_small)

# --------- Variants on plasmid ----------
v <- variants %>%
  filter(chr == plasmid_id) %>%
  mutate(pos = start + 1) %>%
  mutate(type = case_when(
    type %in% c("SNP","SNV") ~ "SNP",
    type %in% c("INS","INSERTION") ~ "INS",
    type %in% c("DEL","DELETION") ~ "DEL",
    TRUE ~ as.character(type)
  ))

stopifnot(nrow(v) > 0)

# --------- Determine plasmid length ----------
plasmid_len <- max(
  coverage_small %>% filter(chr == plasmid_id) %>% pull(end),
  v %>% pull(pos),
  na.rm = TRUE
)

# --------- Variant bins (100 bp) ----------
bins <- tibble(
  bin = 0:floor((plasmid_len - 1) / BIN),
  bin_mid = (bin * BIN) + BIN/2
)

vb <- v %>%
  mutate(bin = floor((pos - 1) / BIN)) %>%
  count(bin, type, name = "n") %>%
  right_join(bins, by = "bin") %>%
  mutate(n = replace_na(n, 0))

vb_total <- vb %>%
  group_by(bin, bin_mid) %>%
  summarise(n = sum(n), .groups="drop")

# --------- Variant plot ----------
p_var <- ggplot() +
  geom_segment(
    data = vb_total,
    aes(x = bin_mid, xend = bin_mid, y = 0, yend = n),
    linewidth = 0.35, alpha = 0.20, color = "grey35"
  ) +
  geom_point(
    data = vb_total,
    aes(x = bin_mid, y = n),
    size = 0.7, alpha = 0.20, color = "grey35"
  ) +
  geom_segment(
    data = vb %>% filter(n > 0),
    aes(x = bin_mid, xend = bin_mid, y = 0, yend = n, color = type),
    linewidth = 0.55
  ) +
  geom_point(
    data = vb %>% filter(n > 0),
    aes(x = bin_mid, y = n, color = type),
    size = 1.2
  ) +
  scale_color_manual(values = c(SNP="#0066FF", INS="#FF9900", DEL="#FF3333")) +
  coord_cartesian(xlim = c(0, plasmid_len), expand = FALSE) +
  labs(
    title = paste0("Variant density (lollipop) across plasmid ", plasmid_id),
    subtitle = paste0("Variant bins: ", BIN, " bp"),
    y = "Variants per bin",
    x = NULL,
    color = "Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# --------- Coverage heatmap strip (TRUE fine windows) ----------
cov_strip <- coverage_small %>%
  filter(chr == plasmid_id) %>%
  mutate(depth_cap = pmin(depth, as.numeric(quantile(depth, 0.99, na.rm = TRUE))))

cov_strip <- cov_strip %>%
  mutate(depth_plot = ifelse(depth == 0, NA, depth_cap))

p_cov <- ggplot(cov_strip) +
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=depth_plot), color=NA) +
  scale_fill_gradientn(colors=c("white","gold","red"), na.value="grey90", name="Coverage")


p_cov <- ggplot(cov_strip) +
  geom_rect(aes(
    xmin = start, xmax = end,
    ymin = 0, ymax = 1,
    fill = depth_cap
  ), color = NA) +
  scale_fill_gradientn(colors = c("white", "gold", "red"), name = "Coverage") +
  coord_cartesian(xlim = c(0, plasmid_len), expand = FALSE) +
  labs(x = "Position on plasmid (bp)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

# --------- Stack ----------
p_var / p_cov + plot_layout(heights = c(4, 0.9))

summary(cov_strip$depth)
quantile(cov_strip$depth, c(.9, .95, .99, 1))

