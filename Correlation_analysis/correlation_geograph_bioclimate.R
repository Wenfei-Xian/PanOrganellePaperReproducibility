if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

input_files <- c("1001G.plus.bio1-19.txt", "1001G_Eastern-Asia.plus.bio1-19.txt")

columns_to_analyze <- c("Longitude", "Latitude", paste0("wc2.1_2.5m_bio_", 1:19))

results <- data.frame(File = character(),
                      Variable = character(),
                      Spearman_Rho = numeric(),
                      P_Value = numeric(),
                      stringsAsFactors = FALSE)

create_plot <- function(data, x_col, y_col, x_label, y_label, annotation_text, file_name) {
  p <- ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = data[[x_col]] > 1.25, shape = data[[x_col]] > 1.25), size = 1) +
    scale_color_manual(values = c("#1E90FF", "#FF4500")) +  
    scale_shape_manual(values = c(3, 4)) +  
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
    labs(x = x_label, y = y_label) +
    theme_minimal() +
    theme(legend.position = "none") +
    annotate("text", x = Inf, y = Inf, label = annotation_text,
             hjust = 1.1, vjust = 1.1, size = 3)
  
  ggsave(file_name, plot = p, device = "pdf", width = 3.5, height = 3)
}

for (input_file in input_files) {
  data <- read.table(input_file, header = TRUE, sep = "\t")
  # colnames(data) <- c("accession", "Longitude", "Latitude", "4k", "6k", "ratio", "wc2.1_2.5m_bio_1", ..., "wc2.1_2.5m_bio_19")
  
  for (col in columns_to_analyze) {
    spearman_test <- cor.test(data$ratio, data[[col]], method = "spearman")
    
    results <- rbind(results, data.frame(
      File = input_file,
      Variable = col,
      Spearman_Rho = as.numeric(spearman_test$estimate),
      P_Value = spearman_test$p.value,
      stringsAsFactors = FALSE
    ))
    
    annotation_text <- sprintf("rho = %.2f\np = %.3g", as.numeric(spearman_test$estimate), spearman_test$p.value)
    
    output_filename <- paste0(tools::file_path_sans_ext(input_file), "_Scatter_ratio_vs_", col, ".pdf")
    
    create_plot(data, 
                x_col = "ratio", 
                y_col = col, 
                x_label = "Coverage Ratio (Repeat I / Repeat II)",
                y_label = col,
                annotation_text = annotation_text,
                file_name = output_filename)
  }
}

print(results)

write.table(results, "correlation_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("done！\n")