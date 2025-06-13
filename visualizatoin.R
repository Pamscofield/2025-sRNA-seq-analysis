# R Script for sRNA-seq Analysis Visualization

# Define the base output directory where processed files are located
aligned_dir <- '/home_local/xiaoqingpan/Downloads/github_protocols/sRNA-seq-analysis/test/output/' # Please confirm this is the correct absolute path

# --- Automatic Strain Detection ---
# This section now automatically detects strain names by listing files
# ending with '_aligned_merged_1st_nt.txt' and extracting the unique prefixes.
# This makes the script more robust if you add or remove strains.
merged_1st_nt_files <- list.files(path = aligned_dir, pattern = "_aligned_merged_1st_nt.txt$", full.names = FALSE)
strain <- unique(str_replace(merged_1st_nt_files, "_aligned_merged_1st_nt.txt", ""))

# Load necessary R packages
# `require()` is used here for brevity; `library()` is also common.
require(ggplot2)    # For creating sophisticated plots
require(ggseqlogo)  # For generating sequence logos
require(data.table) # For efficient data handling, especially `fread`
require(tidyverse)  # A collection of R packages designed for data science (includes dplyr, readr, stringr, etc.)
require(scales)     # For scale functions in ggplot2, like `comma`
require(ggthemes)   # Provides extra themes for ggplot2

# --- First Nucleotide Analysis and Visualization ---

# Initialize a list to store first nucleotide data for each strain
first_nt <- list()

# Loop through each detected strain to read and process first nucleotide data
for (s in strain) {
  # Read the tab-separated file for the first nucleotide, without column names
  tmp2 <- read_tsv(paste0(aligned_dir, s, '_aligned_merged_1st_nt.txt'),
                   col_names = FALSE)
  
  # Change 'T' to 'U' in the nucleotide column (common for RNA sequences)
  tmp2$X1 <- str_replace_all(tmp2$X1, 'T', 'U')
  
  # Store the processed data in the list, named by strain
  first_nt[[s]] <- tmp2$X1 # Extracting only the relevant column
}

# Define a custom color scheme for the sequence logo (A, U, C, G)
cs1 = make_col_scheme(chars = c('A', 'U', 'C', 'G'),
                      cols = c("#E69F00", "#0072B2", "#D55E00", "#56B4E9"))

# Generate the sequence logo plot
# `first_nt` is the list of nucleotide vectors
# `ncol = 3` arranges plots in 3 columns if multiple strains
# `method = 'prob'` uses probability to determine letter height
# `col_scheme` applies the custom colors
ggseqlogo(first_nt,
          ncol = 3,
          stack_width = 0.9,
          method = 'prob',
          col_scheme = cs1)

# Save the first nucleotide plot as a PDF
ggsave(paste0(aligned_dir, '1st_nt.pdf'))

# --- Read Length Distribution Analysis and Visualization ---

# Initialize an empty data frame to store merged length data
length_merge <- data.frame()

# Loop through each detected strain to read and process read length data
for (s in strain) {
  # Read the length file using fread for efficiency, then add a 'strain' column
  tmp <- fread(paste0(aligned_dir, s, '_aligned_merged_length.txt')) %>%
    mutate(strain = s)
  
  # Combine data from current strain with the overall merged data frame
  length_merge <- rbind(tmp, length_merge)
}

# Count occurrences of each read length for each strain
# This step is typically for inspection/summary, not directly used in the plot below.
length_merge %>%
  group_by(V1, strain) %>%
  count()

# Create the read length distribution histogram
ggplot(length_merge,
       aes(V1)) +
  geom_histogram(binwidth = 1, # Use binwidth to control bin size for integer lengths
                 fill = "#0072B2", # Optional: Set a fill color
                 color = "white") + # Optional: Set a border color for bars
  facet_wrap(~strain, ncol = 3) + # Create separate plots for each strain, arranged in 3 columns
  scale_x_continuous(limits = c(10, 30), # Set X-axis limits for read length
                     breaks = seq(10, 30, by = 2)) + # Optional: Set specific breaks for readability
  scale_y_continuous(labels = comma) + # Format Y-axis labels with commas for large numbers
  xlab('Read Length') + # X-axis label
  ylab('Count of Reads') + # Y-axis label (changed for clarity)
  theme_base() + # Apply a base theme from ggthemes
  theme(
    axis.title.y = element_text(margin = margin(r = 20), size = 18), # Adjust Y-axis title margin and size
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1), # Adjust X-axis text size and angle
    axis.text.y = element_text(size = 16), # Adjust Y-axis text size
    # panel.border = element_blank(), # Removed as theme_base() handles some of this, and explicit blanking can be tricky.
    panel.background = element_blank(), # Remove panel background
    plot.background = element_blank(), # Remove plot background
    strip.text = element_text(size = 14, face = "bold") # Adjust facet label appearance
  )

# Save the length distribution plot as a PDF
ggsave(paste0(aligned_dir, 'length_distribution.pdf'))
