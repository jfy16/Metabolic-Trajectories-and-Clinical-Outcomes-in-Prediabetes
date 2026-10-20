# Run with: Rscript circular_bar_plot_2013.R
# Alternatively, use source("/full/path/to/circular_bar_plot_2013.R").
# Keep 2013.txt in the same folder as this script.
# Input columns: id (optional), individual, FC, and logPw or logP.
# logP/logPw must already contain the values used for the -log10(P) color scale.
# Output: circular_bar_plot_2013.pdf in the same folder.

# ==============================
# Circular bar plot
# FC controls bar height; the displayed FC range is fixed at 0.5-2.5.
# logP or logPw controls bar color.
# ==============================

# Load the required packages.
library(ggplot2)
library(tidyverse)
library(ragg)
library(scales)
library(grid)

# ------------------------------
# Read the data.
# ------------------------------
# Locate this script when run with Rscript or source().
# When running selected lines interactively, set the working directory first.
get_script_dir <- function() {
  source_files <- lapply(sys.frames(), function(frame) frame$ofile)
  source_files <- Filter(
    function(path) is.character(path) && length(path) == 1L && nzchar(path),
    source_files
  )
  if (length(source_files) > 0L) {
    source_file <- source_files[[length(source_files)]]
    return(dirname(normalizePath(source_file, mustWork = TRUE)))
  }
  file_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_args) > 0L) {
    script_file <- sub("^--file=", "", file_args[[1L]])
    return(dirname(normalizePath(script_file, mustWork = TRUE)))
  }
  getwd()
}

script_dir <- get_script_dir()
input_file <- file.path(script_dir, "2013.txt")
output_file <- file.path(script_dir, "circular_bar_plot_2013.pdf")

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

data <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_columns <- c("individual", "FC")
missing_columns <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0L) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}
if (nrow(data) == 0L) {
  stop("The input file contains no data rows.")
}

head(data)

# ------------------------------
# Automatically select the color column.
# Prefer logPw; use logP if logPw is unavailable.
# ------------------------------
if ("logPw" %in% colnames(data)) {
  color_var <- "logPw"
} else if ("logP" %in% colnames(data)) {
  color_var <- "logP"
} else {
  stop("Neither logPw nor logP was found in the input data. Check the column names.")
}

# Generate sequential IDs if the id column is absent.
if (!"id" %in% colnames(data)) {
  data$id <- 1:nrow(data)
}

# Sort the data by id.
data <- data %>% arrange(id)

# ------------------------------
# Fix the displayed FC range at 0.5-2.5.
# Plot FC_plot while preserving the original FC values.
# ------------------------------
fc_min_show <- 0.5
fc_max_show <- 2.5

data$FC_plot <- pmin(pmax(data$FC, fc_min_show), fc_max_show)

# Radial FC tick marks.
fc_breaks <- seq(0.5, 2, by = 0.5)

# Inner radius setting for the empty center.
inner_blank <- -0.5

# Offset labels beyond the bar ends.
label_offset <- 0.05

# ------------------------------
# Prepare the label data.
# ------------------------------
label_data <- data
number_of_bar <- nrow(label_data)

angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

# ------------------------------
# Configure the color scale.
# ------------------------------
logp_values <- data[[color_var]]
logp_values <- logp_values[is.finite(logp_values)]

sig_cut <- 1.30103   # -log10(P) corresponding to P = 0.05.

logp_min <- min(logp_values, na.rm = TRUE)
logp_cap <- as.numeric(quantile(logp_values, 0.95, na.rm = TRUE))
logp_cap <- max(logp_cap, sig_cut + 0.5)

# Avoid a degenerate color scale when values are nearly identical.
if (abs(logp_cap - logp_min) < 1e-8) {
  logp_cap <- logp_min + 1e-6
}

# Initial color-scale breakpoints.
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

# Keep only breakpoints within the color-scale limits.
keep_in_range <- raw_breaks >= logp_min & raw_breaks <= logp_cap
raw_breaks <- raw_breaks[keep_in_range]
raw_colors <- raw_colors[keep_in_range]

# Remove duplicate breakpoints to prevent warnings.
keep_unique <- c(TRUE, diff(raw_breaks) > 1e-8)
color_breaks <- raw_breaks[keep_unique]
color_colors <- raw_colors[keep_unique]

# Retain at least two color-scale breakpoints.
if (length(color_breaks) < 2) {
  color_breaks <- c(logp_min, logp_cap)
  color_colors <- c("#6A91C3", "#DF5952")
}

color_values <- rescale(color_breaks, from = c(logp_min, logp_cap))

legend_breaks <- pretty(c(logp_min, logp_cap), n = 5)
legend_breaks <- legend_breaks[legend_breaks >= logp_min & legend_breaks <= logp_cap]

# ------------------------------
# Create the plot.
# ------------------------------
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
  
  # Add radial FC reference lines.
  geom_hline(
    yintercept = fc_breaks,
    color = "grey80",
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  
  # Add numeric FC tick labels.
  annotate(
    "text",
    x = 1,
    y = fc_breaks,
    label = fc_breaks,
    color = "grey30",
    size = 3,
    hjust = 1
  ) +
  
  # Add metabolite labels.
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

# Display the plot.
print(p)


# Export the PDF.
# ------------------------------
ggsave(
  filename = output_file,
  plot = p,
  width = 10,
  height = 10,
  bg = "white"
)

message("Saved PDF: ", output_file)
