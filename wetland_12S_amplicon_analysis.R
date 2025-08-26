################################################################################
####~~~~~~~~~~~~~~~~~~ Maine wetland 12S metabarcoding 2023 ~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~ Code for analyzing Illumina data ~~~~~~~~~~~~~~~~~~~####
################################################################################

# 0. RUN CUSTOM FUNCTIONS SCRIPT ----------------------
source("custom_functions_settings.R")

# 1. LOAD PACKAGES INTO R -----------------------------

# Define CRAN and Bioconductor packages
cran_packages <- c(
  "readxl", "lubridate", "reshape2", "readr", "vegan", "ape", "lme4", "lmerTest", "emmeans",
  "permute", "randomForest", "ggplot2", "ggpubr", "corrplot", "stringr", "tidyr", "tidyverse", 
  "tibble", "RColorBrewer", "corrplot", "viridis","microbiomeutilities"
)

bioc_packages <- c(
  "dada2", "phyloseq", "phylotools", "DESeq2", "DECIPHER", "ade4", "ecodist", "microbiome",
  "ranacapa"
)

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing CRAN packages
cran_missing <- cran_packages[!(cran_packages %in% rownames(installed.packages()))]
if (length(cran_missing) > 0) {
  install.packages(cran_missing)
}

# Install missing Bioconductor packages
bioc_missing <- bioc_packages[!(bioc_packages %in% rownames(installed.packages()))]
if (length(bioc_missing) > 0) {
  BiocManager::install(bioc_missing, ask = FALSE, update = TRUE)
}

# Load all packages
all_packages <- c(cran_packages, bioc_packages)
invisible(lapply(all_packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE))
}))

# 2. LOAD DATA FROM ASV PROCESSING SCRIPT----------------------- 
seqtab.decontam <- readRDS('data/results/seqtab_decontam.rds')
ps.dat.samps <- readRDS('data/results/phyloseq_data_samples_wetlands_2023.rsd')

# 3. OMIT NON-TARGET TAXONOMY AND UNASSIGNED ASVs ------------------------
# make a new object to store phyloseq sample data
ps.tax.samps <- ps.dat.samps 
ps.tax.samps <- subset_taxa(ps.dat.samps,
                            (Phylum == "Chordata" & 
                               ((Order != "Primates" | is.na(Order)) &
                                  (Family != "Bovidae" | is.na(Family)) &
                                  (Species %notin% c("Felis catus",
                                                     "Gallus gallus", 
                                                     "Gallus gallus gallus",
                                                     "Canis lupus familiaris") | is.na(Species)))))

# fix the reptile class
tax.tab <- data.frame(tax_table(ps.tax.samps))
tax.tab$Class[tax.tab$Class == "character(0)"] <- "Reptilia"
# fix the moronid fishes
tax.tab$Order[tax.tab$Order == "character(0)"] <- "Acanthuriformes"
# reassign corrected tax table to ps object
tax.tab <- as.matrix(tax.tab)
tax_table(ps.tax.samps) <- tax.tab

### Remove all OTUs that are 0s across all samples
# function to prune missing OTUs from list of ps objects
prune_missing_otus <- function(ps) {
  # Get the taxa sums, i.e., the total count for each OTU across all samples
  otu_sums <- taxa_sums(ps)
  # Identify OTUs with no positive reads (sums equal to 0)
  taxa_to_keep <- otu_sums > 0
  # Prune the phyloseq object by keeping only the OTUs with positive reads
  ps_pruned <- prune_taxa(taxa_to_keep, ps)
  return(ps_pruned)
}

ps.dat.samps <- prune_missing_otus(ps.dat.samps)
ps.tax.samps <- prune_missing_otus(ps.tax.samps)

# summarize total ASVs and taxa at different taxonomic scales
# total reads, ASVs, and taxonomic levels for each assay (omitting unknown taxa)
tax.sum.tab <- data.frame(Samples = NA,
                          Reads = NA,
                          ASVs = NA,
                          Phylum = NA,
                          Class = NA,
                          Order = NA,
                          Family = NA,
                          Genus = NA,
                          Species = NA)
tax.dat <- psmelt(ps.tax.samps) %>% dplyr::filter(Abundance > 0)
# sample size
tax.sum.tab$Samples <- nrow(otu_table(ps.tax.samps))
# reads
tax.sum.tab$Reads <- sum(otu_table(ps.tax.samps))
# ASVs
tax.sum.tab$ASVs <- length(unique(tax.dat$OTU))
# taxonomic levels
for(level in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  tax.sum.tab[,level] <- length(na.omit(unique(tax.dat[,level])))
}

# print summary table
tax.sum.tab
#Samples   Reads ASVs Phylum Class Order Family Genus Species
#140     3084553  228      1     5    22     49    74      84


# 4. AGGLOMERATE TAXA BY GROUPS -------
# agglomerate reads by different taxonomic scales
split_species = function(string, n = 2) {
  splits = str_split(string, "/", n + 1)
  res = map_if(splits, ~length(.x) > 2, ~.x[1:n]) %>%
    map_chr(str_c, collapse = "/")
  return(res)
}
add_taxonomy_column = function(physeq, num_species = 2) {
  tax_df = as.data.frame(tax_table(physeq)@.Data) %>%
    rownames_to_column("OTU") %>%
    mutate(Species = split_species(Species, n = num_species)) %>%
    mutate(Taxonomy =
             case_when(
               is.na(Class)  ~ Phylum,
               is.na(Order)  ~ Class,
               is.na(Family)  ~ Order,
               is.na(Genus)   ~ Family,
               is.na(Species) ~ Genus,
               TRUE ~ Species)
    )
  tax = as.matrix(tax_df[, -1])
  rownames(tax) = tax_df$OTU
  tax_table(physeq) = tax_table(tax)
  return(physeq)
} 

# custom groupings (lowest taxonomic group)
ps.tax.samps.cust <- add_taxonomy_column(ps.tax.samps)
ps.tax.samps.cust <- tax_glom(ps.tax.samps.cust, "Taxonomy", NArm = TRUE)

# species-level groupings
ps.tax.samps.spp <- tax_glom(ps.tax.samps, "Species", NArm = TRUE)

# 5. READ-BASED RAREFACTION CURVES -----
## 5a. For only taxonomically-verified metabarcoding data ----
nstep <- 100

# prune data to have total read counts above the rarefaction step size
rar.dat <- prune_samples(sample_sums(ps.tax.samps) > nstep, ps.tax.samps)

# get rarefaction curves in ggplot with ggrare from the ranacapa package
ggrar.test.verts <- ggrare(rar.dat, step = nstep, se = FALSE) +
  labs(x = "Sample Size", y = "ASVs") +
  ggtitle("Only vertebrates") +
  theme_bw() + mythemes

ggsave("figures/sequence.based.rarefaction.verts.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

## 5b. For all metabarcoding data -----
# prune data to have total read counts above the rarefaction step size
rar.dat <- prune_samples(sample_sums(ps.dat.samps) > nstep, ps.dat.samps)

# get rarefaction curves in ggplot
ggrar.test.all <- ggrare(rar.dat, step = nstep, se = FALSE) +
  labs(x = "Sample Size", y = "ASVs") +
  ggtitle("All ASVs") +
  theme_bw() + mythemes

ggrar.test.all
ggsave("figures/sequence.based.rarefaction.all.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

ggarrange(ggrar.test.verts, ggrar.test.all,
          nrow = 2, labels = c("A", "B"))

ggsave("figures/sequence.based.rarefaction.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

# one sample is really quite an outlier in terms of sequencing depth

## 5c. Rarefy samples to even sampling depth ----

# rarefy without replacement 
# (depth of 0.9 times min after omitting samples with few [i.e., <1000] reads)
# rarefy BEFORE using tax_glom
set.seed(1)

# rarefaction
ps.rar.samps <- 
  rarefy_even_depth(ps.tax.samps, rngseed=1, 
                    sample.size=0.9*min(sample_sums(ps.tax.samps)[
                      sample_sums(ps.tax.samps) > 1000]), replace=F)

# custom groupings (lowest taxonomic group)
ps.rar.samps.cust <- add_taxonomy_column(ps.rar.samps)
ps.rar.samps.cust <- tax_glom(ps.rar.samps.cust, "Taxonomy", NArm = TRUE)

# species-level groupings
ps.rar.samps.spp <- tax_glom(ps.rar.samps, "Species", NArm = TRUE)

# 6. GENERAL OVERVIEW OF BIODIVERSITY-----------------------

# some general figures
tax.level <- "Class"
rar.glom <- tax_glom(ps.rar.samps, taxrank = tax.level, NArm=FALSE)
rar.melt <- psmelt(rar.glom) %>%
  group_by(Site, Wetland_type2, !!sym(tax.level)) %>%
  dplyr::summarize(Reads = sum(Abundance)) %>%
  filter(Reads > 0) %>%
  mutate(Site = factor(Site)) %>%
  ungroup()

# community composition by site
ggplot(rar.melt, aes(Site, Reads, fill = !!sym(tax.level))) +
  geom_bar(stat = "identity", position = "fill", col = "black") +
  labs(x = "Site", y = "Prop. reads") +
  scale_x_discrete(drop=FALSE) +
  theme_bw() + mythemes + 
  scale_fill_manual(values = brewer.pal(n = 6, name = "Accent")) +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "top") +
  guides(fill=guide_legend(ncol=2))

# community composition by site and wetland type
ggplot(rar.melt, aes(Wetland_type2, Reads, fill = !!sym(tax.level))) +
  geom_bar(stat = "identity", position = "fill", col = "black") +
  labs(x = "Wetland type", y = "Prop. reads") +
  scale_x_discrete(drop=FALSE) +
  theme_bw() + mythemes + 
  coord_flip() +
  scale_fill_manual(values = brewer.pal(n = 6, name = "Accent")) +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.key.size = unit(0.4, "cm")) +
  guides(fill=guide_legend(ncol=2))

# calculate diversity metrics
Diversity <- estimate_richness(ps.rar.samps.spp, measures=c("Observed", "Shannon"))
Diversity$Evenness <- Diversity$Shannon/log(Diversity$Observed)
Diversity$ReadFile_Name <- rownames(Diversity)
Diversity <- left_join(Diversity, metadata.runs, by = "ReadFile_Name")

# plot diversity by individual sites
ggplot(Diversity, aes(reorder(Site, Shannon), Shannon, fill = Wetland_type2)) +
  geom_boxplot() + 
  coord_flip() +
  labs(x = "Site", y = "Shannon diversity") +
  theme_bw() + mythemes

# plot diversity by wetland type
Diversity_melt_all <- melt(Diversity, 
                           id.vars = c("ReadFile_Name", "Site", "Wetland_type2"), 
                           measure.vars = c("Observed", "Shannon", "Evenness"),
                           variable.name = "Index")

ggplot(Diversity_melt_all, aes(Wetland_type2, value)) +
  geom_violin(trim=TRUE, fill = "grey", col = "grey") + 
  geom_boxplot(width = 0.1, aes(group = Wetland_type2)) + 
  theme_bw() + mythemes + labs(x = "Wetland type", y = "Index", title = "All vertebrates") +
  facet_wrap(~Index, scales = "free") +
  theme(legend.position = "none") 

# plot diversity by wetland type (omitting fish)
ps.rar.samps.spp.nofish <- subset_taxa(ps.rar.samps.spp, Class != "Actinopteri")

# calculate diversity metrics
Diversity <- estimate_richness(ps.rar.samps.spp.nofish, measures=c("Observed", "Shannon"))
Diversity$Evenness <- Diversity$Shannon/log(Diversity$Observed)
Diversity$ReadFile_Name <- rownames(Diversity)
Diversity <- left_join(Diversity, metadata.runs, by = "ReadFile_Name")

# plot diversity by wetland type
Diversity_melt_nofish <- melt(Diversity, 
                              id.vars = c("ReadFile_Name", "Site", "Wetland_type2"), 
                              measure.vars = c("Observed", "Shannon", "Evenness"),
                              variable.name = "Index")

ggplot(Diversity_melt_nofish, aes(Wetland_type2, value)) +
  geom_violin(trim=TRUE, fill = "grey", col = "grey") + 
  geom_boxplot(width = 0.1, aes(group = Wetland_type2)) + 
  theme_bw() + mythemes + labs(x = "Wetland type", y = "Index", title = "No fish") +
  facet_wrap(~Index, scales = "free") +
  theme(legend.position = "none") 

# plot both together
library(ggh4x)
Diversity_melt_all$Data <- "All vertebrates"
Diversity_melt_nofish$Data <- "No fish"
Diversity_melt_combined <- rbind(Diversity_melt_all, Diversity_melt_nofish)
ggplot(Diversity_melt_combined, aes(Wetland_type2, value)) +
  geom_violin(trim=TRUE, fill = "grey", col = "grey") + 
  geom_boxplot(width = 0.1, aes(group = Wetland_type2)) + 
  theme_bw() + mythemes + labs(x = "Wetland type", y = "Index") +
  facet_grid2(Data ~ Index, scales = "free", independent = "all") +
  theme(legend.position = "none") 

## IDEAS FOR ANALYSIS / MANUSCRIPTS
### 1. eDNA vs relative abundance with traditional monitoring methods (egg mass counts, trapping, cameras)
#### 1a. eDNA vs egg-mass counts
# We only have data for spotted salamander eggs (n = 15 sites)
amma_data <- 
  subset_samples(ps.tax.samps.spp, Egg_mass_survey == "Yes") %>%
  psmelt() %>% 
  mutate(eDNA_day = paste(month(Date, label = T, abbr = T), day(Date))) %>%
  filter(Species == "Ambystoma maculatum") %>%
  group_by(Site, eDNA_day) %>% 
  summarize(Reads = mean(Abundance)/Area_sq_m, Eggs = mean(AMMA_eggs)/Area_sq_m) %>%
  ungroup()

ggplot(amma_data, aes(Eggs, Reads)) +
  geom_point() +
  geom_smooth(method = "lm")



#### cluster analysis to compare communities detected from either method
#### paired bar plots comparing species reads vs number of detections (something quantitative and z-scored)
#### compare simple summary of richness from either method
#### regression relating reads vs camera detections

### 2. Comparing eDNA biodiversity between pools and other forested wetlands

### 3. Misc. species summaries
# amphibian occurence table
amphibian_data <- subset_taxa(ps.tax.samps.spp, Class == "Amphibia") %>%
  psmelt() %>% group_by(Site, Species) %>% summarize(Reads = mean(Abundance)) %>%
  ungroup()

amphibian_site_table <- dcast(amphibian_data, Site ~ Species, value = "Reads")

# summary of data for private landowners
reidy_sites_melt <- subset_samples(ps.tax.samps.cust, Site %in% c("CR-1", "CR-2", "CR-Beaver", "CR-Beaver-CT")) %>% psmelt()
reidy_tax_summary <- reidy_sites_melt %>% group_by(Site, Class) %>% filter(Abundance > 0) %>% summarize(taxa = unique(Taxonomy))

jamison_sites_melt <- subset_samples(ps.tax.samps.cust, Site %in% c("GB-Willow", "GB-Essex")) %>% psmelt()
jamison_tax_summary <- jamison_sites_melt %>% group_by(Site, Class) %>% filter(Abundance > 0) %>% summarize(taxa = unique(Taxonomy))

larson_sites_melt <- subset_samples(ps.tax.samps.cust, Site %in% c("Larson")) %>% psmelt()
larson_tax_summary <- larson_sites_melt %>% group_by(Site, Class) %>% filter(Abundance > 0) %>% summarize(taxa = unique(Taxonomy))

