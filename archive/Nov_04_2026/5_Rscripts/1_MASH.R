# **********************************************************************************
# Loading packages 
# **********************************************************************************
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# if(!require("pacman"))
#   install.packages("pacman")

# pacman::p_load("openxlsx", "vegan", "ggpubr", "ggrepel", "dplyr", "RColorBrewer", "plotly")

# **********************************************************************************
# Loading data
# **********************************************************************************

lines <- readLines("mash_dist.txt")  # Replace with your file path

# Extract sample names (first column)
sample_names <- sapply(strsplit(lines, "\t"), `[[`, 1)
n <- length(sample_names)

# Initialize a square matrix of zeros
mash_table <- matrix(0, nrow = n, ncol = n)
rownames(mash_table) <- colnames(mash_table) <- sample_names

# Fill the lower triangle from the file
for (i in seq_along(lines)) {
  parts <- strsplit(lines[i], "\t")[[1]]
  if (length(parts) > 1) {
    mash_table[i, 1:(length(parts) - 1)] <- as.numeric(parts[-1])
  }
}

# Make the matrix symmetric by copying the lower triangle to the upper
mash_table[upper.tri(mash_table)] <- t(mash_table)[upper.tri(mash_table)]

# Print the full symmetric distance matrix
print(mash_table)

# edit names to match with metadata
rownames(mash_table) <- gsub("\\.fq.gz", "", rownames(mash_table))
colnames(mash_table) <-rownames(mash_table)
mash_table[1:3,1:3]
dim(mash_table) #

mds.env  <- read.xlsx("MASH_ENV.xlsx",sheet="seased_env", rowNames=TRUE)
dim(mds.env) # 
mds.env <- mds.env[rownames(mds.env) %in% rownames(mash_table), ]
dim(mds.env) # 
mds.env[1:3,1:3]

mash_table <- mash_table[rownames(mash_table) %in% rownames(mds.env) , colnames(mash_table) %in% rownames(mds.env)]
dim(mash_table) # 

# **********************************************************************************
# distance ordination
# **********************************************************************************

dist_obj <- as.dist(mash_table)
# dist_ord <- metaMDS(dist_obj, k = 3, trymax = 25)
dist_ord <- metaMDS(dist_obj, k = 3, trymax = 1000)

# View NMDS stress (lower is better; < 0.2 is usually acceptable)
dist_ord$stress

ordiplot(dist_ord, display = 'sites', type = 'text')
dist_ord$eig

data.scores <- as.data.frame(scores(dist_ord, display = "sites")) 

merged_df <- merge(data.scores, mds.env, by = "row.names")
merged_df

# with sample label
ggplot(merged_df, aes(NMDS1, NMDS2, fill= Source)) +
  geom_point(shape=21, size=4) +
  coord_equal() + 
  theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
  xlab("dist_ord1") +
  ylab("dist_ord2") + 
  labs(fill = "SourceSeason") +
  geom_text_repel(aes(label = ReefName), vjust = -1, hjust = 0.5, size = 3, colour = "black") +
  stat_ellipse(aes(color = Source, group = Source), linetype = 2, level = 0.90)

# save figure as svg
x <- ggplot(merged_df, aes(NMDS1, NMDS2, fill = Source)) +
  geom_point(aes(shape = TripDescriotion), size = 4) +  
  # coord_equal() + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
  xlab("NMDS1") +
  ylab("NMDS2") + 
  labs(fill = "Source", shape = "TripDescriotion") +
  # geom_text_repel(aes(label = ReefName), vjust = -1, hjust = 0.5, size = 3, colour = "black") +
  scale_fill_manual(values = c("Seawater" = "#1f78b4", "Sediment" = "#33a02c")) +
  scale_color_manual(values = c("Seawater" = "#1f78b4", "Sediment" = "#33a02c")) +
  stat_ellipse(aes(color = Source, group = Source), linetype = 2) +
  coord_fixed(ratio = 1) +
  annotate("text", 
           label = paste0("Stress = ", round(dist_ord$stress, 3)), 
           x = Inf, y = Inf, 
           hjust = 1.2, vjust = 2.5, 
           size = 4) + 
  scale_shape_manual(values = c(21, 22, 23, 24)) + coord_fixed(ratio = 0.5)

# ggsave(file="mash_dist_NMDS.svg", plot=x, width=7.5, height=7.5)

# **********************************************************************************
# COMPARISON
# **********************************************************************************

# Convert dist object to matrix
dist_mat <- as.matrix(dist_obj)

# Get group info
group <- mds.env$Source

# Number of observations
n <- length(group)

# Create all pairwise combinations (upper triangle only)
combinations <- expand.grid(row = 1:n, col = 1:n)
combinations <- combinations[combinations$row < combinations$col, ]

# Add group info
combinations$group_row <- group[combinations$row]
combinations$group_col <- group[combinations$col]

# Extract distances
combinations$distance <- dist_mat[as.matrix(combinations[, c("row", "col")])]

# Within-group summary
within_group <- combinations %>%
  filter(group_row == group_col) %>%
  group_by(group = group_row) %>%
  summarise(
    min_distance = min(distance),
    max_distance = max(distance),
    mean_distance = mean(distance),
    sd_distance = sd(distance),
    .groups = 'drop'
  )

print(within_group)

# PERMANOVA
vegan::adonis2(dist_obj ~ Source, method="bray", sqrt.dist = TRUE, data = mds.env, permutations = 999)

dispersion<-vegan::betadisper(dist_obj, group=mds.env$Source)
dispersion

vegan::permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE)

# ANOSIM 
x <- vegan::anosim(dist_obj, mds.env$Source)
summary(x)
plot(x)

# Checking 3rd dimension 
data.scores_3d <- as.data.frame(scores(dist_ord, display = "sites"))
merged_df_d3 <- merge(data.scores_3d, mds.env, by = "row.names")

nb.cols <- length(unique(merged_df_d3$Source))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

plot_ly(x=merged_df_d3$NMDS1, y=merged_df_d3$NMDS2, z=merged_df_d3$NMDS3,
        name=merged_df_d3$Source, type="scatter3d", mode="markers",
        color=merged_df_d3$Source, colors = mycolors) %>%
  add_trace(type = 'scatter3d', size = 1,
            x = merged_df_d3$NMDS1, y = merged_df_d3$NMDS2, z = merged_df_d3$NMDS3,
            opacity=0.01,color=merged_df_d3$Source,  colors = mycolors)

