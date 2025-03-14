#!/usr/bin/env Rscript

# command
args <- commandArgs(trailingOnly = TRUE)

# three files are needed
if (length(args) != 3) {
  stop("./spearman.plot.structure.r genotype.txt phenotype.txt PC_file.txt\n")
}

suppressMessages(library(ggplot2))
suppressMessages(library(ggExtra))

file1_path <- args[1]
file2_path <- args[2]
pc_file_path <- args[3]

path_parts <- strsplit(file2_path, "/")[[1]]
file2_name <- if (length(path_parts) > 1) path_parts[length(path_parts) - 1] else "Unknown"

hash_table <- list()

file1 <- readLines(file1_path)
for (line in file1) {
  line_split <- strsplit(line, "\t")[[1]]
  if (length(line_split) >= 2 && !is.na(as.numeric(line_split[2]))) {
    id <- line_split[1]
    num <- as.numeric(line_split[2])
    hash_table[[id]] <- num
  }
}

value1_list <- c()
value2_list <- c()
id_list <- c()

file2 <- readLines(file2_path)
for (line in file2) {
  line_split <- strsplit(line, "\t")[[1]]
  if (length(line_split) >= 2 && !is.na(as.numeric(line_split[2]))) {
    id <- line_split[1]
    num2 <- as.numeric(line_split[2])
    if (!is.null(hash_table[[id]])) {
      value1_list <- c(value1_list, hash_table[[id]])
      value2_list <- c(value2_list, num2)
      id_list <- c(id_list, id)
    }
  }
}

pc_data <- read.table(pc_file_path, header = TRUE, sep = "\t", comment.char = "")

colnames(pc_data) <- gsub("^#", "", colnames(pc_data))

pc_data$IID <- trimws(as.character(pc_data$IID))
id_list <- trimws(as.character(id_list))

pc_data <- pc_data[pc_data$IID %in% id_list, ]

pc_data <- pc_data[match(id_list, pc_data$IID), ]

df <- data.frame(
  Value1 = value1_list,
  Value2 = value2_list,
  PC1 = pc_data$PC1,
  PC2 = pc_data$PC2,
  PC3 = pc_data$PC3,
  PC4 = pc_data$PC4  #four PCs
)

if (length(value1_list) > 1) {
  suppressWarnings({
    # Spearman
    spearman_result <- cor.test(value1_list, value2_list, method = "spearman")
    spearman_r <- round(spearman_result$estimate, 4)
    spearman_p <- format(spearman_result$p.value, scientific = TRUE, digits = 4)
  })

  suppressWarnings({
    # Multiple linear regression
    model <- lm(Value2 ~ Value1 + PC1 + PC2 + PC3 + PC4, data = df)
    model_summary <- summary(model)

    beta_value1 <- round(model_summary$coefficients["Value1", "Estimate"], 4)
    beta_p <- format(model_summary$coefficients["Value1", "Pr(>|t|)"], scientific = TRUE, digits = 4)
  })

  plot_data <- data.frame(
    Value1 = value1_list,
    Value2 = value2_list,
    Color = ifelse(value1_list > 1.25, "#FF4500", "#1E90FF"),
    Category = ifelse(value1_list > 1.25, "Above", "Below")
  )

  p <- ggplot(plot_data, aes(x = Value1, y = Value2, shape = Category, color = Color)) +
    geom_point(size = 1) +
    scale_shape_manual(values = c("Above" = 4, "Below" = 3)) + 
    geom_smooth(aes(x = Value1, y = Value2), method = "lm", color = "black", se = FALSE, linewidth = 0.5, inherit.aes = FALSE) +
    labs(
      x = "Coverage Ratio (Repeat I / Repeat II)",
      y = file2_name
    ) +
    scale_color_identity() +
    theme_minimal(base_size = 6) +
    theme(
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8),
      legend.position = "none"
    ) +
    annotate(
      "text",
      x = min(value1_list),
      y = max(value2_list),
      label = paste0(
        "Spearman rho = ", spearman_r, "\nP = ", spearman_p,
        "\nMLR Beta = ", beta_value1, "\nP = ", beta_p
      ),
      hjust = 0,
      vjust = 1,
      size = 2,
      color = "black"
    )

  p_with_density <- ggExtra::ggMarginal(
    p, type = "density",
    xparams = list(fill = "#8B008B", alpha = 0.8, size = 0.1),
    yparams = list(fill = "#000080", alpha = 0.8, size = 0.1)
  )

  output_pdf <- paste0("scatterplot_with_density_", file2_name, ".pdf")
  pdf(output_pdf, width = 3.3, height = 3.3)
  print(p_with_density)
  dev.off()
  print(paste("#####", file2_path, spearman_p, beta_p, sep = " "))

  cat("Saved", output_pdf, "\n")
} else {
  cat("!!! Error\n")
}

