if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(terra)
library(dplyr)

# Define the directory containing the bio raster files
bio_data_dir <- "wc2.1_2.5m_bio"  # path

# List and load the raster files
bio_files <- list.files(bio_data_dir, pattern = "wc2\\.1_2\\.5m_bio_\\d+\\.tif$", full.names = TRUE)
bio_raster <- rast(bio_files) 

# Define the input files to process
input_files <- c("1001G.txt", "1001G_Eastern-Asia.txt")

# Loop through each input file, process and output the result
for (input_file in input_files) {
  
  # Read the input data
  data <- read.table(input_file, header = TRUE, sep = "\t")
  colnames(data) <- c("accession", "Longitude", "Latitude", "4k", "6k", "ratio")
  
  # Extract bio values using the coordinates (Longitude and Latitude)
  coords <- data[, c("Longitude", "Latitude")]
  bio_values <- terra::extract(bio_raster, coords, xy = FALSE)
  
  # Combine the original data with the extracted bio values
  data_with_bio <- bind_cols(data, as.data.frame(bio_values))
  
  # Create an output file name based on the input file name
  output_file <- paste0(tools::file_path_sans_ext(input_file), ".plus.bio1-19.txt")
  
  # Write the combined data to the output file
  write.table(data_with_bio, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Saved", output_file, "\n")
}