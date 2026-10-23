
library(ggplot2)
library(tidyverse)
library(ragg)
library(scales)
library(grid)


data <- read.table(
  "2016-FDR 0.2.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

head(data)


if ("logPw" %in% colnames(data)) {
  color_var <- "logPw"
} else if ("logP" %in% colnames(data)) {
  color_var <- "logP"
} else {
  stop("check column names")
}


if (!"id" %in% colnames(data)) {
  data$id <- 1:nrow(data)
}


data <- data %>% arrange(id)


fc_min_show <- 0.5
fc_max_show <- 2.5
data$FC_plot <- pmin(pmax(data$FC, fc_min_show), fc_max_show)
fc_breaks <- seq(0.5, 2, by = 0.5)
inner_blank <- -0.5


label_offset <- 0.05

label_data <- data
number_of_bar <- nrow(label_data)

angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

logp_values <- data[[color_var]]
logp_values <- logp_values[is.finite(logp_values)]

sig_cut <- 1.30103

logp_min <- min(logp_values, na.rm = TRUE)
logp_cap <- as.numeric(quantile(logp_values, 0.95, na.rm = TRUE))
logp_cap <- max(logp_cap, sig_cut + 0.5)

if (abs(logp_cap - logp_min) < 1e-8) {
  logp_cap <- logp_min + 1e-6
}

raw_breaks <- c(
  logp_min,
  sig_cut,
  sig_cut + (logp_cap - sig_cut) * 0.35,
  sig_cut + (logp_cap - sig_cut) * 0.70,
  logp_cap
)

raw_colors <- c(
  "#6A91C3",
  "#A7CCE3",
  "#FFFFCC",
  "#FDBE81",
  "#DF5952"
)

keep_in_range <- raw_breaks >= logp_min & raw_breaks <= logp_cap
raw_breaks <- raw_breaks[keep_in_range]
raw_colors <- raw_colors[keep_in_range]

keep_unique <- c(TRUE, diff(raw_breaks) > 1e-8)
color_breaks <- raw_breaks[keep_unique]
color_colors <- raw_colors[keep_unique]

if (length(color_breaks) < 2) {
  color_breaks <- c(logp_min, logp_cap)
  color_colors <- c("#6A91C3", "#DF5952")
}

color_values <- rescale(color_breaks, from = c(logp_min, logp_cap))

legend_breaks <- pretty(c(logp_min, logp_cap), n = 5)
legend_breaks <- legend_breaks[legend_breaks >= logp_min & legend_breaks <= logp_cap]


p <- ggplot(data, aes(x = as.factor(id), y = FC_plot, fill = .data[[color_var]])) +
  
  geom_bar(
    stat = "identity",
    width = 1,
    color = "#525252",
    linewidth = 0.3
  ) +
  
  scale_fill_gradientn(
    colors = color_colors,
    values = color_values,
    limits = c(logp_min, logp_cap),
    oob = squish,
    breaks = legend_breaks,
    labels = round(legend_breaks, 2),
    name = expression(-log[10](italic(P)))
  ) +
  
  scale_y_continuous(
    limits = c(inner_blank, fc_max_show + label_offset),
    breaks = fc_breaks,
    labels = fc_breaks,
    expand = c(0, 0)
  ) +
  
  coord_polar(start = 0) +
  
  geom_hline(
    yintercept = fc_breaks,
    color = "grey80",
    linewidth = 0.3,
    linetype = "dashed"
  ) +

  annotate(
    "text",
    x = 1,
    y = fc_breaks,
    label = fc_breaks,
    color = "grey30",
    size = 3,
    hjust = 1
  ) +
  
  geom_text(
    data = label_data,
    aes(
      x = id,
      y = FC_plot + label_offset,
      label = individual,
      hjust = hjust
    ),
    color = "black",
    fontface = "bold",
    size = 2.2,
    angle = label_data$angle,
    inherit.aes = FALSE
  ) +
  
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = unit(rep(0.5, 4), "cm"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 9, color = "black")
  )

print(p)


ggsave(
  filename = "circular_bar_plot_2016.pdf",
  plot = p,
  width = 10,
  height = 10,
  bg = "white"
)

