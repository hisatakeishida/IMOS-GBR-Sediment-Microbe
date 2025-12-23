# **********************************************************************************
# Loading packages 
# **********************************************************************************
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# if(!require("pacman"))
#   install.packages("pacman")
# 
# install.packages("microViz")
# remotes::install_github("cafferychen777/MicrobiomeStat")
# 
# pacman::p_load("microbiome","openxlsx", "dplyr", "knitr", "upstartr", "hrbrthemes", "RColorBrewer", "microViz", "skimr", "Maaslin2", "bit64", "ggrepel", "data.table", "ggpubr", "MicrobiomeStat", "mixOmics", "vegan", "eulerr", "adespatial", "psych", "RVAideMemoire", "ggVennDiagram")
# library("ComplexHeatmap")


# **********************************************************************************
# colors to use 
# **********************************************************************************

# c("#f2b09e", # Sampling trip 1
# "#5b5b5b", # Sampling trip 2 
# "#f3c789", # Sampling trip 3
# "slateblue")#  Sampling trip 4

# c("Seawater" = "#1f78b4", "Sediment" = "#33a02c")

# **********************************************************************************
# Loading dataset
# **********************************************************************************

#load 
sample_data <- sample_data(read.xlsx("~/OneDrive - The University of Queensland/1_IMOS/1_SingleM_ANALYSIS/IMOS_10Mseqs_SingleM_FINAL.xlsx", sheet = "seased_env_RS226", rowNames=TRUE))
otu <- read.xlsx("~/OneDrive - The University of Queensland/1_IMOS/1_SingleM_ANALYSIS/IMOS_10Mseqs_SingleM_FINAL.xlsx", sheet = "seased_cov_RS226", rowNames=TRUE)
tax_rel <- as.matrix(read.xlsx("~/OneDrive - The University of Queensland/1_IMOS/1_SingleM_ANALYSIS/IMOS_10Mseqs_SingleM_FINAL.xlsx", sheet = "seased_tax_RS226", rowNames=TRUE))
otu_rel <- sweep(otu, 2, colSums(otu), FUN = "/")
otu_rel_scaled <- otu_rel

# Initialize the new row for nonprok_scale
nonprok_scale_row <- numeric(length(colnames(otu_rel)))
names(nonprok_scale_row) <- colnames(otu_rel)

#scale relative abundance based on coverage with library size and prokaryotic fraction estimated by singlem microbial_fraction 
for (col in colnames(otu_rel)) {
# Get library size 
library_size <- sample_data$metagenome_size[rownames(sample_data) == col]
# Get prokaryotic fraction 
prok_fraction <- sample_data$read_fraction[rownames(sample_data) == col]
# Get non-prokaryotic fraction 
nonprok_fraction <- sample_data$nonprok_fraction[rownames(sample_data) == col]
# Calculate the scaling factor 
prok_scale <- (library_size * prok_fraction) / 100
nonprok_scale <- (library_size * nonprok_fraction) / 100

# Multiply the column values by the corresponding prok_scale value
otu_rel_scaled[[col]] <- otu_rel[[col]] * prok_scale
nonprok_scale_row[col] <- nonprok_scale
}

# Add the nonprok_scale_row as a new row to the otu_rel_scaled table
otu_rel_scaled <- rbind(otu_rel_scaled, nonprok_scale_row)
rownames(otu_rel_scaled)[nrow(otu_rel_scaled)] <- paste0("otu", nrow(otu_rel_scaled))

#round and convert into integer
otu_rel_scaled <- otu_rel_scaled %>% mutate_all(~ as.numeric(round(.)))

# Add the nonprok_scale_row as a new row to the tax table
new_row <- rep("", ncol(tax_rel))
names(new_row) <- colnames(tax_rel)
new_row[1] <- "SingleM non-bac/arc"
tax_rel <- rbind(tax_rel, new_row)
rownames(tax_rel)[nrow(tax_rel)] <- paste0("otu", nrow(otu_rel_scaled))

# Setting up phyloseq object
sample_data <- sample_data(sample_data)
otu_scaled_in <- otu_table(otu_rel_scaled, taxa_are_rows = TRUE)
tax_rel_in <- tax_table(tax_rel)
pseq = phyloseq(otu_table = otu_scaled_in, sample_data=sample_data, taxonomyTable =tax_rel_in)
lib.size <- apply(t((otu_table(pseq))), 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

pseq <- tax_fix(pseq, min_length = 1)
pseq <- phyloseq_validate(pseq)

# **********************************************************************************
# Composition plot with non-euks
# **********************************************************************************
######
DOMAIN<- pseq %>%
  comp_barplot("Domain", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# ggsave(file="~/OneDrive - The University of Queensland/1_IMOS/1_SingleM_ANALYSIS/singlem_10Mseqs.RS220_seased_composition_domainWithnonprok.png", plot=DOMAIN, width=7.5, height=7.5)


pseq %>%
  comp_barplot("Phylum", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

pseq %>%
  comp_barplot("Class", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

pseq %>%
  comp_barplot("Order", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

pseq %>%
  comp_barplot("Family", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

pseq %>%
  comp_barplot("Genus", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

#####

# **********************************************************************************
#  Only prokaryotic 
# **********************************************************************************

keep_bac <- c("Bacteria")
keep_arc <- c("Archaea")
keep_bacarc<- c("Bacteria", "Archaea")

pseq_bac <-  subset_taxa(pseq, Domain %in% keep_bac)
pseq_arch <-  subset_taxa(pseq, Domain %in% keep_arc)
pseq_bac_arch <-  subset_taxa(pseq, Domain %in% keep_bacarc)

# **********************************************************************************
# Composition plot
# **********************************************************************************
pseq_bac_arch %>%
  comp_barplot("Domain", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

pseq_bac_arch %>%
  comp_barplot("Phylum", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source,TripSeason), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

phylum <-pseq_bac_arch %>%
  comp_barplot("Phylum", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

class <- pseq_bac_arch %>%
  comp_barplot("Class", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))
class

class2 <- pseq_bac_arch %>%
  comp_barplot("Class", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source,TripDescriotion), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))
class2

faam <- pseq_bac_arch %>%
  comp_barplot("Family", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_wrap(vars(Source,TripDescriotion), scales = "free") + # scales = "free" is IMPORTANT!
  coord_flip() +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))
faam

# **********************************************************************************
# Ordination
# **********************************************************************************

pseq_bac_arch_family <- pseq_bac_arch %>%
  tax_transform(trans = "hellinger", rank = "Family") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    size = 6,
    fill = "Source", shape = "TripDescriotion",
    plot_taxa = 10:1,
    tax_vec_length = 0.5,
    tax_vec_style_all=
      vec_constraint(linewidth = 0, arrow = NULL),
    tax_lab_style = tax_lab_style(
      type = "label", max_angle = 90, aspect_ratio = 1,  
      size = 3, alpha = 0.9, fontface = "bold", # style the labels
      label.r = unit(0, "mm") # square corners of labels - see ?geom_label
    )
  ) +   stat_ellipse(aes(linetype = Source, colour = Source), linetype = 0.2)+ 
  theme_bw() + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), aspect_ratio = 1) +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) + 
  scale_color_manual(values = c("Seawater" = "#1f78b4", "Sediment" = "#33a02c"))  + theme(legend.position = "none")+ 
  coord_fixed(xlim=c(-0.75,0.75), ylim=c(-0.75,0.75))
# +  coord_equal()
pseq_bac_arch_family

# ggsave(file="seased_FAM_PCA.v2.svg", plot=pseq_bac_arch_family, width=7.5, height=7.5)

# **********************************************************************************
# core family
# **********************************************************************************
pseq.rel <- pseq_bac_arch  %>% tax_transform("compositional", rank = "Family") 

#By source 
source_group <- unique(as.character(meta(pseq.rel)$Source))
print(source_group)

list_core <- c() # an empty object to store information
for (n in source_group){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Source == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0/100, # 0.001 in atleast 90% samples 
                         prevalence = 0.80)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
list_core

remove_ranks <- c(
  "Domain",
  "Phylum",
  "Phylum_Class",
  "Class",
  "Order"
)

pattern <- paste0(" (", paste(remove_ranks, collapse = "|"), ")$")

list_core_clean <- lapply(list_core, function(x) {
  x[!grepl(pattern, x)]
})

common_taxa <- intersect(
  list_core_clean$Seawater,
  list_core_clean$Sediment
)

ggVennDiagram(list_core, force_upset = TRUE, order.set.by = "name")
ggVennDiagram(list_core_clean, force_upset = TRUE, order.set.by = "name")

common_taxa

mycols <- c("Seawater" = "#1f78b4", "Sediment" = "#33a02c")
plot(venn(list_core),fills = mycols)
ggVennDiagram(list_core, force_upset = TRUE, order.set.by = "name")

