# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(ggtext)
# Specify the file path
file_path <- "jan24/N50_Stool.xlsx"

# Read data from different sheets with sheet index
data_group1 <- read_excel(file_path, sheet = 1)  # Assuming the first sheet corresponds to "N50 < 1,000"
data_group2 <- read_excel(file_path, sheet = 2)  # Assuming the second sheet corresponds to "N50 = 1,000 - 10,000"
data_group3 <- read_excel(file_path, sheet = 3)  # Assuming the third sheet corresponds to "N50 > 10,000"

# Combine data from different sheets into one data frame
combined_data <- bind_rows(
  data_group1 %>% mutate(Group = "N50 < 1,000"),
  data_group2 %>% mutate(Group = "N50 = 1,000 - 10,000"),
  data_group3 %>% mutate(Group = "N50 > 10,000")
)


## Compact graph  without the facets:
scatter_plot_2<- ggplot(combined_data, aes(x = N50, y = MAG_comp, shape=Direct_seq_detection,color = Direct_seq_detection)) +
  geom_point() +
  #facet_wrap(~ Group, ncol = 1, scales = "free_x") +  # Use free_x for independent X-axis scales
  labs(
    title = "",
    x = "N50",
    y = "*Campylobacter* genome completeness [%]",
    color = "MAG *Campylobacter* species detection",
    shape = "MAG *Campylobacter* species detection"
  ) +
  scale_color_manual(
    values = c("Yes" = "black", "NO" = "red"),
    labels = c("Yes" = "*Campylobacter* species identified", "NO" = "*Campylobacter* species not identified")
  ) +  # Set legend colors and labels without italics
  scale_shape_discrete(
    labels = c("Yes" = "*Campylobacter* species identified", "NO" = "*Campylobacter* species not identified")
  ) +  # Set legend colors and labels without italics
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),  # Adjust background color
    panel.grid.major = element_line(color = "grey"),  # Adjust grid line color
    axis.title.y=element_markdown(),
    legend.text = element_markdown(),
    legend.title = element_markdown(),
  )
scatter_plot_compact = scatter_plot_2 + scale_x_log10(breaks=c(200,500,1000,2000,5000,10000,20000,50000,100000,200000),
                               labels=scales::label_comma())

ggsave("scatter_plot_compact.PNG", plot = scatter_plot_compact, width = 8, height = 5, units = "in", dpi = 300, bg = "white")


# Create a scatter plot for each group with free X-axis scales
scatter_plot <- ggplot(combined_data, aes(x = N50, y = MAG_comp, color = Direct_seq_detection)) +
  geom_point() +
  facet_wrap(~ Group, ncol = 1, scales = "free_x") +  # Use free_x for independent X-axis scales
  labs(
    title = "",
    x = "N50",
    y = "*Campylobacter* genome completeness [%]",
    color = "MAG *Campylobacter* species detection"
  ) +
  scale_color_manual(
    values = c("Yes" = "black", "NO" = "red"),
    labels = c("Yes" = "*Campylobacter* species identified",
               "NO" = "*Campylobacter* species not identified")
  ) +  # Set legend colors and labels without italics
  guides(color = guide_legend(reverse = TRUE)) +  # Reverse the order of the legend
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),  # Adjust background color
    panel.grid.major = element_line(color = "grey"),  # Adjust grid line color
    axis.title.y=element_markdown(),
    legend.text = element_markdown(),
    legend.title = element_markdown()
  )

# Print the scatter plot
scatter_plot
ggsave("scatter_plot_2.PNG", plot = scatter_plot, width = 10, height = 8, units = "in", dpi = 300, bg = "white")



######## merge horizontally
# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)

# Specify the file path
file_path <- "/users/djeghout/desktop/N50_Stool.xlsx"

# Define a function to read and preprocess data from a sheet
read_and_preprocess_data <- function(sheet_name, group_label) {
  data <- read_excel(file_path, sheet = sheet_name)
  data$Group <- group_label
  return(data)
}

# Read and preprocess data from different sheets
data_group1 <- read_and_preprocess_data(1, "N50 < 1,000")
data_group2 <- read_and_preprocess_data(2, "N50 = 1,000 - 10,000")
data_group3 <- read_and_preprocess_data(3, "N50 > 10,000")

# Combine data from different sheets into one data frame
combined_data <- bind_rows(data_group1, data_group2, data_group3)

# Create a scatter plot for each group with free X-axis scales
scatter_plot <- ggplot(combined_data, aes(x = N50, y = MAG_comp, color = Direct_seq_detection)) +
  geom_point() +
  facet_wrap(~ Group, ncol = 1, scales = "free_x") +  # Use free_x for independent X-axis scales
  labs(title = "N50, Campylobacter MAGs completeness and direct sequencing detection",
       x = "N50",
       y = "Campylobacter MAGs completeness [%]",
       color = "MAG species detection") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "grey95"),  # Adjust background color
    panel.grid.major = element_line(color = "black")  # Adjust grid line color
  )

# Print the scatter plot
print(scatter_plot)

 # Save the scatter plot as an image file
ggsave("combined_scatter_plot.PNG", plot = scatter_plot, width = 6, height = 8, units = "in", dpi = 300)


