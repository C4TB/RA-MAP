# Load required packages
library(tidyverse)
library(M3C)
library(ggplot2)
library(ConsensusClusterPlus)
library(SummarizedExperiment)
library(bioplotr)
library(NMF)
library(RColorBrewer)

# Summary ----
# Script to perform M3C consensus clustering on the RAMAP multi-platform data
# First perform PCA on each platform and extract first 10 principal components
# Perform M3C and plot heatmap of resulting clusters
# 
# After M3C get the cluster assignment and use bioplotr functions to look at the association between
# consensus cluster and each component for each platform (driver plots) adjusting for age, gender and baseline SDAI.
# This was to recreate a plot that Dennis originally created in JMP. 
# 
# Because of the way in which we decided to do the driver plots, taking the association between c-cluster and 
# each platform and plotting them together into a single figure, it requires a lot of awkward / ugly transposing of the 
# data and bioplotr functions wrapped up in wrapper functions. Sorry about that! 


### STEP 1 PCA platforms ----

# Read data 
full <- read_tsv("multiomics_baseline_common_104_patients.tsv")
full <- full[-1]  # remove first column (e.g. rownames)

# Get platform names and separate into list
platforms <- c("CLIN", "CD14_", "CD4", "CD8", "HLA", "MIRN", 
               "PBMC", "PROT", "SERA", "SNP", "URIN", "WHL")

get_indices <- function(pattern, anno) grep(pattern, anno)
split_data <- map(setNames(platforms, platforms), ~ full[get_indices(.x, full$ANNO), ])

# create PCA function
pca_top_components <- function(df, setname, n_pc = 10) {
  df <- as.data.frame(df)
  mat <- setNames(data.frame(t(df[,-1])), df[,1])
  pcs <- prcomp(mat)$x
  pcs <- as.data.frame(pcs)[, 1:n_pc]
  colnames(pcs) <- paste0(setname, "_PC", seq_len(n_pc))
  return(pcs)
}

# Run PCA on all platforms (except clinical data)
pca_results <- map2(split_data[names(split_data) != "CLIN"], 
                    names(split_data)[names(split_data) != "CLIN"], 
                    ~ pca_top_components(.x, .y))

# Get sample IDs
ID <- setNames(data.frame(ID = rownames(pca_results[[1]])), "ID")

# Merge all PCA data
merged_pca <- bind_cols(ID, pca_results)

# Transpose for clustering
pcs = setNames(data.frame(t(merged_pca[,-1])), merged_pca[,1])



### STEP 2 run M3C on PCAs of platforms  ----

# MC3 wont run properly using current sampleIDs because there is a leading '0' in some IDs
# so create a new vector of sample names and create a lookup for original sampleIDs 

sampleIDs <- setNames(as.data.frame(colnames(pcs)), "SampleID")
sampleIDs$var = paste0("V", 1:104) 

# change colnames of pcs df to the new IDs
colnames(pcs) = sampleIDs$var

# Run M3C
set.seed(123)
m3c_results <- M3C(pcs, maxK = 10, clusteralg = "pam", iters = 100)
# Whilst running this spits out the cluster stability plots

# Plot consensus clustering heatmap 
ccmatrix <- m3c_results$realdataresults[[2]]$consensus_matrix
annon <- m3c_results$realdataresults[[2]]$ordered_annotation
pal <- rev(colorRampPalette(brewer.pal(9, "Reds"))(256))

aheatmap(ccmatrix, annCol = annon[1], Colv = NA, Rowv = NA, 
         cexRow = 0, cexCol = 0, color = rev(pal), scale = "none")


### STEP 3 bioplotr driver plot ----

# get consensus cluster assignment and sampleIDs. 
# sampleIDs are carried over using rownames. Create column instead.
annon$var = rownames(annon)
annon = left_join(annon, sampleIDs) # map back to sampleIDs using lookup above


# get clinical data, transpose and transfer SampleIDs into a column  
clin_s = as.data.frame(split_data$CLIN)
clin_s = setNames(data.frame(t(clin_s[,-1])), clin_s[,1]) # transpose
clin_s$SampleID = row.names(clin_s)

# add the specific clinical variables to consensus cluster assignment
annon = left_join(annon, select(clin_s, SampleID, CLIN_AGE, CLIN_GENDER, CLIN_SDAI0M))


generate_driver_plot_from_split <- function(platform_data, platform_name, annon, select_clin_vars = c("consensuscluster", "CLIN_AGE", "CLIN_GENDER", "CLIN_SDAI0M")) {
  
  # Transpose, using ANNO / platform names as colnames
  expr_mat <- setNames(data.frame(t(platform_data[,-1])), platform_data[,1])
  expr_mat$ID <- rownames(expr_mat)

  # Merge with annotation (which holds sample metadata and consensuscluster)
  merged <- merge(annon, expr_mat, by.x = "SampleID", by.y = "ID")

  # Split expression and phenotype
  pdata <- merged[1:6]
  assay_mat <- setNames(data.frame(t(merged[,-c(1:6)])), pdata$var)
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment(assay_mat, colData = pdata)

  cnts <- assay(se)
  clins <- as_tibble(colData(se)) %>% select(all_of(select_clin_vars))

  # Run plot_drivers
  plot_out <- plot_drivers(cnts, clins, bivariate = FALSE, n_pc = 10, p_adj = "BH")
  plot_out$Platform <- platform_name
  
  return(plot_out)
}



# Apply to each platform
plot_list <- setNames(
  lapply(platforms, function(name) {
    generate_driver_plot_from_split(split_data[[name]], name, annon)
  }),
  platforms
)

# the warnings here just refer to a deprecated scale argument

# Combine each platform association results from bioplotr into a df
plot_data <- do.call(rbind, lapply(names(plot_list), function(name) {
  df <- plot_list[[name]]$data
  df$Platform <- name
  df
}))


# Keep only "consensuscluster" rows as these are the only ones we wanted to plot
plot_data <- as.data.frame(plot_data[plot_data$Feature == "consensuscluster" & plot_data$Platform != "CLIN", ])


# Define platform and PC order
level_order <- paste0("PC", 1:10)
level_order2 <- c( "CD4", "CD14_","CD8", "HLA", "MIRN", 
                  "PBMC", "PROT", "SERA", "SNP", "URIN", "WHL")


# Build plot
p <- ggplot(plot_data, aes(x = factor(PC, level = level_order),
                           y = factor(Platform, level = level_order2),
                           fill = Association,
                           text = Association,
                           color = Significant)) +
  geom_tile(size = 1L, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("grey90", "black")) +
  guides(color = "none") +
  scale_fill_gradient2(midpoint = 15, low = "white", mid = "red", high = "brown", name = expression(~-log[10]~italic(p))) +
  labs(title = "Variation By Feature",
       x = "Principal Component",
       y = "Platform") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

