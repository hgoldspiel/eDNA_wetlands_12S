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
                                  (Species %notin% c("Gallus gallus", 
                                                     "Gallus gallus gallus") | is.na(Species)))))
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
#1     140 3085888  232      1     5    22     49    75      86

write.csv(tax.sum.tab, "data/assigned_taxa_sums.csv", row.names = FALSE)


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
  group_by(Site, Wetland_type, !!sym(tax.level)) %>%
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
  guides(fill=guide_legend(nrow=1))

ggsave("figures/community.composition.site.jpg", 
       width = 8, height = 5, units = "in", dpi = 300)

# community composition by site and wetland type
ggplot(rar.melt, aes(factor(Wetland_type, levels = rev(c("Pool", "Semipermanent", "Permanent"))), 
                            Reads, fill = !!sym(tax.level))) +
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
ggplot(Diversity, aes(reorder(Site, Shannon), Shannon, 
                      fill = factor(Wetland_type, levels = c("Pool", "Semipermanent", "Permanent")))) +
  geom_boxplot() + 
  coord_flip() +
  labs(x = "Site", y = "Shannon diversity") +
  theme_bw() + mythemes +
  guides(fill=guide_legend(title = "Wetland type"))

ggsave("figures/ASV.diversity.shannons.site.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

# plot diversity by wetland type
Diversity_melt_all <- reshape2::melt(Diversity, 
                                     id.vars = c("ReadFile_Name", "Site", "Wetland_type"), 
                                     measure.vars = c("Observed", "Shannon", "Evenness"),
                                     variable.name = "Index")

Diversity_melt_all$Wetland_type <- factor(Diversity_melt_all$Wetland_type,
                                          levels = c("Pool", "Semipermanent", "Permanent"))
  
ggplot(Diversity_melt_all, aes(Wetland_type, value)) +
  geom_violin(trim=TRUE, fill = "grey", col = "grey") + 
  geom_boxplot(width = 0.1, aes(group = Wetland_type)) + 
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
Diversity_melt_nofish <- reshape2::melt(Diversity, 
                                        id.vars = c("ReadFile_Name", "Site", "Wetland_type"), 
                                        measure.vars = c("Observed", "Shannon", "Evenness"),
                                        variable.name = "Index")

Diversity_melt_nofish$Wetland_type <- factor(Diversity_melt_nofish$Wetland_type,
                                          levels = c("Pool", "Semipermanent", "Permanent"))

ggplot(Diversity_melt_nofish, aes(Wetland_type, value)) +
  geom_violin(trim=TRUE, fill = "grey", col = "grey") + 
  geom_boxplot(width = 0.1, aes(group = Wetland_type)) + 
  theme_bw() + mythemes + labs(x = "Wetland type", y = "Index", title = "No fish") +
  facet_wrap(~Index, scales = "free") +
  theme(legend.position = "none") 

# plot both together
library(ggh4x)
Diversity_melt_all$Data <- "All vertebrates"
Diversity_melt_nofish$Data <- "No fish"
Diversity_melt_combined <- rbind(Diversity_melt_all, Diversity_melt_nofish)
ggplot(Diversity_melt_combined, aes(Wetland_type, value)) +
  geom_violin(trim=TRUE, fill = "grey", col = "grey") + 
  geom_boxplot(width = 0.1, aes(group = Wetland_type)) + 
  theme_bw() + mythemes + labs(x = "Wetland type", y = "Index") +
  facet_grid2(Data ~ Index, scales = "free", independent = "all") +
  theme(legend.position = "none") 

# 7. ANALYSIS -------
## 7a. eDNA validation with conventional surveys (egg mass counts, trapping, cameras) ----
### 7a-1. eDNA vs egg-mass counts ----

# We only have data for spotted salamander eggs (n = 15 sites) 
# summarized read count data per survey
# normalize reads by total reads from sample?
raw.read.counts <- data.frame(ReadFile_Name = sample.names.shortF,
                              Raw_reads = filtoutput[,"reads.in"])

amma_data <- 
  subset_samples(ps.tax.samps.spp, Egg_mass_survey == "Yes") %>%
  psmelt() %>% 
  left_join(raw.read.counts, by = "ReadFile_Name") %>%
  mutate(eDNA_day = paste(month(Date, label = T, abbr = T), day(Date)),
         DOY = yday(Date),
         Hatch = ifelse(Date < "2023-06-01", "No", "Yes")) %>%
  filter(Species == "Ambystoma maculatum" & !is.na(AMMA_eggs)) %>%
  group_by(Site, DOY, eDNA_day, Hatch) %>% 
  dplyr::summarize(Reads = mean(Abundance),
            Reads_norm = mean(Abundance/Raw_reads),
            Reads_per_ml = mean(Abundance/Sample_Volume_ml),
            Reads_per_ml_norm = mean((Abundance/Raw_reads)/Sample_Volume_ml),
            Temp = mean(Temp..C),
            pH = mean(pH)) %>%
  left_join(metadata.sites, by = c("Site" = "SiteID")) %>%
  ungroup()

# add volume data and adjust other variables relative to volume
# we will treat wetlands as paraboloids (bowls) for a rough vol. calculation
amma_data$Wetland_vol_m3 <- 0.5*(amma_data$Spring_depth_cm/100)*amma_data$Area_sq_m
# correct egg and read counts
amma_data$AMMA_eggs_per_m2 <- amma_data$AMMA_eggs / amma_data$Area_sq_m
amma_data$Reads_per_ml_m3 <- amma_data$Reads_per_ml / amma_data$Wetland_vol_m3
amma_data$Reads_per_ml_m3_norm <- amma_data$Reads_per_ml_norm / amma_data$Wetland_vol_m3

# relationship, uncorrected
(raw.amma.egg.cor <- ggplot(amma_data, aes(AMMA_eggs, Reads)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() + mythemes +
  labs(x = "Egg masses", y = "eDNA read count"))

# relationship, corrected for vol
ggplot(amma_data, aes(AMMA_eggs_per_m2, Reads_per_ml_m3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() + mythemes

# relationship, corrected for vol, log transformed data (check that its kosher)
ggplot(amma_data, aes(log(AMMA_eggs_per_m2), log(Reads_per_ml_m3+1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() + mythemes

### better, but lots of noise still

# relationship with time of eDNA sample (before or after egg hatching)
ggplot(amma_data, aes(AMMA_eggs_per_m2, Reads_per_ml_m3, color = Hatch)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() + mythemes

# log transformed reads
ggplot(amma_data, aes(log(AMMA_eggs_per_m2), log(Reads_per_ml_m3+1), color = Hatch)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() + mythemes

# eDNA vs date
(raw.amma.date.cor <- ggplot(amma_data, aes(DOY, Reads)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() + mythemes +
  labs(x = "DOY", y = "eDNA read count"))

# gamm to account for nonlinear effects of time of eDNA sample
library(mgcv)
amma_data$Site <- as.factor(amma_data$Site)
eggXdate_eDNA_gamm <- 
  gam(Reads_per_ml_m3 ~ 
        te(AMMA_eggs_per_m2, DOY, bs = c("tp", "tp"), m = 2) + 
        s(Site, bs = "re"),
      nthreads = 2, 
      data = amma_data)

egg_eDNA_gamm <- 
  gam(Reads_per_ml_m3 ~ 
        s(AMMA_eggs_per_m2, bs = "tp") +
        s(Site, bs = "re"),
      nthreads = 2, 
      data = amma_data)

date_eDNA_gamm <- 
  gam(Reads_per_ml_m3 ~ 
        s(DOY, bs = "tp") +
        s(Site, bs = "re"),
      nthreads = 2, 
      data = amma_data)

null_eDNA_gamm <- 
  gam(Reads_per_ml_m3 ~ 1 +
        s(Site, bs = "re"),
      nthreads = 2, 
      data = amma_data)

library(bbmle)
AICtab(eggXdate_eDNA_gamm, egg_eDNA_gamm, date_eDNA_gamm, null_eDNA_gamm)
# looks like DOY is stronger influence on read counts than egg mass counts

summary(eggXdate_eDNA_gamm) # barely significant interaction term

# pred vs observed reads
qplot(amma_data$Reads_per_ml_m3, predict(eggXdate_eDNA_gamm)) +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed") +
  theme_classic() + mythemes

## let's create prediction surfaces 
pred.grads <- expand.grid(AMMA_eggs_per_m2 = seq(min(amma_data$AMMA_eggs_per_m2, na.rm = T),
                                                 max(amma_data$AMMA_eggs_per_m2, na.rm = T)*1.01, 
                                                 length.out = 30),
                          DOY = seq(min(amma_data$DOY), 
                                    max(amma_data$DOY)*1.01, 
                                    1),
                          Site = "RE")

egg.date.pred <- as.data.frame(predict(eggXdate_eDNA_gamm, exclude = c("s(Site)"),
                                       newdata = pred.grads, se.fit = TRUE))

egg.date.pred <- cbind(pred.grads, egg.date.pred)

library(scico)
ggplot(egg.date.pred, aes(x = AMMA_eggs_per_m2, y = DOY, z = fit)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  geom_contour(colour = "white") + 
  geom_point(data = amma_data, inherit.aes = FALSE, aes(x = AMMA_eggs_per_m2, y = DOY), shape = "+", 
             color = "grey", size = 4) +
  scale_fill_scico(palette = 'davos') +
  theme_bw() + mythemes +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar(title = "Reads/ml/m3"))

# the gam seems to be overfitting

# use linear model 

# parametric model to assess significance with lmer
library(lme4)

# inspect data distribution
hist(amma_data$Reads_per_ml_m3, prob = T); 
lines(density(amma_data$Reads_per_ml_m3), col = "blue", lwd = 2)

hist(sqrt(amma_data$Reads_per_ml_m3), prob = T); 
lines(density(sqrt(amma_data$Reads_per_ml_m3)), col = "blue", lwd = 2)

hist(log(amma_data$Reads_per_ml_m3+1), prob = T); 
lines(density(log(amma_data$Reads_per_ml_m3 + 1)), col = "blue", lwd = 2)

# seems like some transformation would be helpful
# this paper used log transformation on reads for validating amphibian eDNA surveys 
# https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/1365-2656.13468
lmer.full <- lmer(sqrt(Reads_per_ml_m3) ~ 
                    AMMA_eggs_per_m2 * DOY + (1|Site), 
                  data = amma_data)
lmer.add <- lmer(sqrt(Reads_per_ml_m3) ~ 
                   AMMA_eggs_per_m2 + DOY + (1|Site), 
                 data = amma_data)
lmer.eggs <- lmer(sqrt(Reads_per_ml_m3) ~ 
                    AMMA_eggs_per_m2 + (1|Site), 
                  data = amma_data)
lmer.DOY <- lmer(sqrt(Reads_per_ml_m3) ~ 
                    DOY + (1|Site), 
                  data = amma_data)
lmer.null <- lmer(sqrt(Reads_per_ml_m3) ~ 1 + (1|Site), 
                  data = amma_data)

library(performance)
library(MuMIn)
check_normality(lmer.full) # no good!
check_heteroscedasticity(lmer.full) # good!
AICtab(lmer.full, lmer.add, lmer.eggs, lmer.DOY, lmer.null) 
anova(lmer.add, lmer.full) # additive model is better
summary(lmer.add)
performance::r2(lmer.add)

## let's create prediction surfaces 
egg.date.pred <- as.data.frame(predict(lmer.add, re.form= NA,
                                       newdata = pred.grads, se.fit = TRUE))

egg.date.pred <- cbind(pred.grads, egg.date.pred)

library(scico)
lmer.pred.eggs.eDNA <- ggplot(egg.date.pred, aes(x = AMMA_eggs_per_m2, y = DOY, z = fit)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  geom_contour(colour = "white") + 
  geom_point(data = amma_data, inherit.aes = FALSE, 
             aes(x = AMMA_eggs_per_m2, y = DOY), shape = "+", 
             color = "grey", size = 4) +
  scale_fill_scico(palette = 'davos') +
  theme_bw() + mythemes +
  labs(x = "Salamander egg masses / m2") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar(title = "sqrt(Reads/ml/m3)"))

lmer.pred.eggs.eDNA

# glmer with lognormal distribution
glmer.full <- glmer(Reads_per_ml_m3 ~ 
                     scale(AMMA_eggs_per_m2) * scale(DOY) + (1|Site), 
                   family = gaussian(link = "log"),
                   data = amma_data,
                   mustart=pmax(amma_data$AMMA_eggs,1e-3))
glmer.add <- glmer(Reads_per_ml_m3 ~ 
                      scale(AMMA_eggs_per_m2) + scale(DOY) + (1|Site), 
                    family = gaussian(link = "log"),
                    data = amma_data,
                   mustart=pmax(amma_data$AMMA_eggs,1e-3))
glmer.eggs <- glmer(Reads_per_ml_m3 ~ 
                     scale(AMMA_eggs_per_m2) + (1|Site), 
                   family = gaussian(link = "log"),
                   data = amma_data,
                   mustart=pmax(amma_data$AMMA_eggs,1e-3))
glmer.DOY <- glmer(Reads_per_ml_m3 ~ 
                    scale(DOY) + (1|Site), 
                  family = gaussian(link = "log"),
                  data = amma_data,
                  mustart=pmax(amma_data$AMMA_eggs,1e-3))
glmer.null <- glmer(Reads_per_ml_m3 ~ 1 + (1|Site), 
                   family = gaussian(link = "log"),
                   data = amma_data,
                   mustart=pmax(amma_data$AMMA_eggs,1e-3))

bbmle::AICctab(glmer.full, glmer.add, glmer.eggs, glmer.DOY, glmer.null, weights = TRUE)
anova(glmer.full, glmer.add)
summary(glmer.add)
performance::r2(glmer.add)

## let's create prediction surfaces 
egg.date.pred <- as.data.frame(predict(glmer.add, re.form = NA, type = "response",
                                       newdata = pred.grads, se.fit = TRUE))

egg.date.pred <- cbind(pred.grads, egg.date.pred)

egg.date.pred.avg <- as.data.frame(modavgPred(list(glmer.full, glmer.add, glmer.eggs, glmer.DOY, glmer.null), 
                                          modnames = c("Full", "Additive", "Eggs", "DOY", "Null"),
                                          re.form = NA, type = "response", 
                                          newdata = pred.grads, se.fit = TRUE))

egg.date.pred.avg <- cbind(pred.grads, egg.date.pred.avg)

# for one model
glmer.pred <- ggplot(egg.date.pred, aes(x = AMMA_eggs_per_m2, y = DOY, z = fit)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  geom_contour(colour = "white") + 
  geom_point(data = amma_data, inherit.aes = FALSE, 
             aes(x = AMMA_eggs_per_m2, y = DOY), shape = "+", 
             color = "grey", size = 4) +
  scale_fill_scico(palette = 'davos') +
  theme_bw() + mythemes +
  labs(x = expression(paste("Egg masses / ", m^{2}))) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar(title = expression(paste("Reads / ml / ", m^{3}))))

glmer.pred

# model averaged estimates
glmer.avg.pred <- ggplot(egg.date.pred.avg, aes(x = AMMA_eggs_per_m2, y = DOY, z = mod.avg.pred)) +
  geom_raster(aes(fill = mod.avg.pred), interpolate = TRUE) +
  geom_contour(colour = "white") + 
  geom_point(data = amma_data, inherit.aes = FALSE, 
             aes(x = AMMA_eggs_per_m2, y = DOY), shape = "+", 
             color = "grey", size = 4) +
  scale_fill_scico(palette = 'davos') +
  theme_bw() + mythemes +
  labs(x = "Salamander egg masses / m2") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar(title = "Reads/ml/m3"))

glmer.avg.pred

# combine some plots for sharing
ggarrange(ggarrange(raw.amma.egg.cor + ylab("eDNA abundance") + xlab("Egg masses") +
                      theme(plot.margin=unit(c(0,1,0,1), 'cm')), 
                    raw.amma.date.cor + ylab("eDNA abundance") +
                      theme(plot.margin=unit(c(0,1,0,1), 'cm')), labels = c("A", "B")),
          glmer.pred, nrow = 2, labels = c("", "C"))

ggsave("figures/egg.mass.validation.jpg", 
       width = 9, height = 7, units = "in", dpi = 300)

### 7a-2. eDNA vs trapping ----
#### load trapping data
traps <- read.csv("data/trapping_data.csv")
traps$Date <- as.Date(traps$Date, format = "%m/%d/%Y")

# create binary, daily detection histories for amphibian species of interest
# convert empty traps to NAs for species
traps$Species[traps$Species %in% c("", "N/A")] <- NA

# fix other species spellings
traps$Species[traps$Species == "AMMA "] <- "AMMA"
traps$Species[traps$Species == "LISY "] <- "LISY"
traps$Species[traps$Species == "AMBY "] <- "AMBY"

# scratch code
# convert variable formats
## traps
factors <- c("Trap_ID", "Site", "Glow", "Metal", "Species")
traps[factors] <- lapply(traps[factors], factor)
traps$Date_f <- as.factor(traps$Date)
traps$Total <- as.numeric(traps$Total)

# add variables for trap night and overall day in the study
traps$TrapNight <- 1 
traps <- traps %>%
  filter(Date != "2023-04-15") %>%
  group_by(Site) %>%
  mutate(TrapNight = case_when(Date == dplyr::first(Date) ~ 1,
                               Date == nth(unique(Date), 2) ~ 2,
                               TRUE ~ 3)) %>%
  ungroup() %>%
  mutate(TrapNight = as.factor(TrapNight)) 

# extract a data frame with all trap x species x site x year combinations
newdat <- list()
alllevs <- do.call(expand.grid, 
                   lapply(traps[c("TrapNight", "Species", "Site")], levels))
for(i in 1:nrow(alllevs)) {
  if(alllevs$TrapNight[i] %in% traps$TrapNight[traps$Site == alllevs$Site[i]]) {
    newdat <- rbind(newdat, alllevs[i,])
  }
  else(next)
}

# append data frame with 0s 
# remove rows where species were trapped for a given site and day
zerodata <- list()
for(i in unique(traps$Site)) {
  for(j in unique(traps$TrapNight[traps$Site == i])) {
    sitezeroes <- newdat[newdat$Site == i & newdat$TrapNight == j,]
    sitezeroes$Total <- NA
    zerodata <- bind_rows(zerodata, 
                          sitezeroes[sitezeroes$Species %notin% 
                                       traps$Species[traps$Site == i &
                                                       traps$TrapNight == j],])
  }
}

# combine with original trap dataset
newtraps <- bind_rows(traps, zerodata) %>%
  arrange(Date, Site, TrapNight) %>%
  filter(!is.na(Species)) %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# fill in missing dates
for(i in 1:nrow(newtraps)) {
  newtraps$Date[i] <- ifelse(is.na(newtraps$Date[i]), 
                             unique(traps$Date[traps$Site == newtraps$Site[i] & 
                                                 traps$TrapNight == newtraps$TrapNight[i]]),
                             newtraps$Date[i])
}

newtraps$Method <- "Trapping"

# format species-specific capture histories for unmarked
library(auk)
trap_occu_data <- list()
for(species in c("AMBY", "AMLA", "AMMA", "LISY")) {
  # aggregate site-level species-specific trapping data
  spp_data <- newtraps[newtraps$Species == species & !is.na(newtraps$Species),] %>%
    group_by(Site, Date, Method) %>%
    summarize(Total = sum(Total)) %>%
    ungroup() %>%
    mutate(Site = if_else(Site == "RW-S", "RW-Stream", Site))
  # reformat for unmarked
  trap_occu_data[[species]] <- format_unmarked_occu(spp_data, 
                                                    site_id = "Site", 
                                                    response = "Total", 
                                                     site_covs = "Method", 
                                                     obs_covs = "Date") %>%
    mutate(across(c(y.1, y.2, y.3), makebinary))
}

#### format eDNA data in same way for occupancy analysis
eDNA_occu_data <- list()
for(species in c("Ambystoma barbouri", "Ambystoma laterale", 
                 "Ambystoma maculatum", "Lithobates sylvaticus")) {
  # aggregate site-level species-specific eDNA data at trapping sites (and locations)
  spp_data <- ps.tax.samps.spp %>%
    psmelt() %>%
    mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
    filter(Species == species &
             Date < "2023-05-05" &
             Trapping == "Yes") %>%
    mutate(Method = "eDNA") %>%
    arrange(Date)
  # reformat for unmarked
  eDNA_occu_data[[species]] <- format_unmarked_occu(spp_data,
                                                    site_id = "Site",
                                                    response = "Abundance",
                                                    site_covs = "Method",
                                                    obs_cov = "Date") %>%
    mutate(across(c(y.1, y.2, y.3, y.4), makebinary))
}

# combine trapping and eDNA data for analysis
combined_occu_data <- list() 
for(spp in 1:length(trap_occu_data)) {
  combined_occu_data[[names(trap_occu_data[spp])]] <-
    bind_rows(trap_occu_data[[spp]], eDNA_occu_data[[spp]])
}

combined_occu_data_wide <- list() 
for(spp in 1:length(trap_occu_data)) {
  combined_occu_data_wide[[names(trap_occu_data[spp])]] <-
    left_join(trap_occu_data[[spp]], eDNA_occu_data[[spp]], by = "site")
}

# run occupancy analysis 
library(ubms)
library(unmarked)
library(AICcmodavg)
occu_mod_outputs <- list()
for(spp in names(combined_occu_data_wide)) {
  # extract species-level occu data and reformat date as DOY
  dat <- combined_occu_data_wide[[spp]] %>%
    mutate(across(starts_with("Date"), yday)) %>%
    mutate(Method.1.x = "Trapping", Method.2.x = "Trapping", Method.3.x = "Trapping",
           Method.1.y = "eDNA", Method.2.y = "eDNA", Method.3.y = "eDNA", Method.4 = "eDNA") %>%
    dplyr::select(-c("Method.x", "Method.y"))
  # reformat data for umbs
  umf_dat <- unmarkedFrameOccu(y = dplyr::select(dat, starts_with("y")),
                               siteCovs = data.frame(Site = dat$site),
                               obsCovs = list(Date = dplyr::select(dat, starts_with("Date")),
                                              Method = dplyr::select(dat, starts_with("Method"))))
  # fit model in unmarked
  mod.list <- list()
  mod.list[["Method * Date"]] <- occu(~Method * Date ~ 1, data = umf_dat)
  mod.list[["Method + Date"]] <- occu(~Method + Date ~ 1, data = umf_dat)
  mod.list[["Method"]] <- occu(~Method ~ 1, data = umf_dat)
  mod.list[["Date"]] <- occu(~Date ~1, data = umf_dat)
  mod.list[["Null"]] <- occu(~1 ~ 1, data = umf_dat)
  
  # save list of models
  occu_mod_outputs[[spp]][["mods"]] <- mod.list
  
  # goodness of fit test 
  occu_mod_outputs[[spp]][["GoF"]] <- mb.gof.test(mod.list[["Method * Date"]], 
                                                  plot.hist = TRUE, nsim = 1000)
  
  # aictab
  occu_mod_outputs[[spp]][["aictab"]] <- aictab(cand.set = mod.list,
                                                modnames = c("Method * Date",
                                                             "Method + Date",
                                                             "Method",
                                                             "Date",
                                                             "Null"),
                                                sort = TRUE,
                                                digits = 4)
}

occu_aic_tables <- cbind(data.frame(Species = rep(names(occu_mod_outputs), 
                                            each = length(occu_mod_outputs[[1]]$mods))),
                         rbind(occu_mod_outputs$AMBY$aictab,
                               occu_mod_outputs$AMLA$aictab,
                               occu_mod_outputs$AMMA$aictab,
                               occu_mod_outputs$LISY$aictab))

# GoF from global models
for(spp in names(occu_mod_outputs)) {
  print(spp)
  print(occu_mod_outputs[[spp]]$GoF)
}

                         
# likelihood ratio tests between additive and interactive
## Unisexual
anovaOD(occu_mod_outputs$AMBY$mods[[2]], occu_mod_outputs$AMBY$mods[[1]]) # additive better
summary(occu_mod_outputs$AMBY$mods$`Method + Date`) # both method and date significant
## Blue-spot
anovaOD(occu_mod_outputs$AMLA$mods[[2]], occu_mod_outputs$AMLA$mods[[1]]) # interactive better
summary(occu_mod_outputs$AMLA$mods$`Method * Date`) # converges but huge errors, no significance
summary(occu_mod_outputs$AMLA$mods$`Method + Date`) # weak effect of date
## Spotted
anovaOD(occu_mod_outputs$AMMA$mods[[2]], occu_mod_outputs$AMMA$mods[[1]]) # additive better
summary(occu_mod_outputs$AMMA$mods$`Method + Date`) # both method and date significant
## Wood frog
anovaOD(occu_mod_outputs$LISY$mods[[2]], occu_mod_outputs$LISY$mods[[1]]) # additive better
summary(occu_mod_outputs$LISY$mods$`Method + Date`) # both method and date significant

# plot effects of top model (or model-averaged predictions)
modavgdat.spp <- NULL
for(spp in names(occu_mod_outputs)) {
  mod.list <- occu_mod_outputs[[spp]][["mods"]]
  modavgdat <- data.frame(expand.grid(c("eDNA", "Trapping"), 
                                      108:123))
  colnames(modavgdat) <- c("Method", "Date")
  modavgdat$Method <- factor(modavgdat$Method, levels = c("Trapping", "eDNA"))
  modavg <- modavgPred(cand.set = mod.list, newdata = modavgdat, parm.type = "detect")
  modavg <- cbind(modavgdat, modavg)
  modavg$Species <- spp
  modavgdat.spp <- bind_rows(modavgdat.spp, modavg)
}

modavgdat.spp$Date.mdy <- as.Date("2023-01-01") + (modavgdat.spp$Date - 1)
ggplot(modavgdat.spp, aes(Date.mdy, mod.avg.pred, color = Method, fill = Method, lty = Method)) +
  geom_line(lwd = 1.1) +
  geom_ribbon(aes(ymin = lower.CL,
                  ymax = upper.CL), alpha = 0.3) +
  facet_wrap(~Species) +
  scale_x_date(limits = c(as.Date("2023-04-18"), as.Date("2023-05-03")),
               breaks = c(as.Date("2023-04-18"), as.Date("2023-04-25"), as.Date("2023-05-02")),
               date_labels = "%b %d") +
  labs(x = "Date", y = "Prob. of detection") +
  scale_color_manual(values = col2[c(1,2)]) +
  scale_fill_manual(values = col2[c(1,2)]) +
  theme_bw() + mythemes

# top model effects, bivariate predictions
fig.labs <- c("Wood Frog", "Spotted Salamander", 
              "Unisexual Salamander", "Blue-spotted Salamander")
fig.lets <- LETTERS[1:4]
names(fig.lets) <- c("LISY", "AMMA", "AMBY", "AMLA")
fig.lets.labs <- paste(fig.lets, "  ", fig.labs)
names(fig.lets.labs) <- c("LISY", "AMMA", "AMBY", "AMLA")

pred.spp <- NULL
for(spp in names(occu_mod_outputs)) {
  newdat <- data.frame(expand.grid(c("eDNA", "Trapping"), 
                                      108:123))
  colnames(newdat) <- c("Method", "Date")
  newdat$Method <- factor(newdat$Method, levels = c("Trapping", "eDNA"))
  pred <- predict(occu_mod_outputs[[spp]]$mods$`Method + Date`,
                    newdata = newdat, type = "det")
  pred <- cbind(newdat, pred)
  pred$Species <- spp
  pred.spp <- bind_rows(pred.spp, pred)
}

pred.spp$Date.mdy <- as.Date("2023-01-01") + (pred.spp$Date - 1)
pred.spp$Species <- factor(pred.spp$Species, levels = c("LISY", "AMMA", "AMBY", "AMLA"))

ggplot(pred.spp, aes(Date.mdy, Predicted, color = Method, fill = Method, lty = Method)) +
  geom_line(lwd = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  facet_wrap(~Species, labeller = as_labeller(fig.lets.labs)) +
  scale_x_date(limits = c(as.Date("2023-04-18"), as.Date("2023-05-03")),
               breaks = c(as.Date("2023-04-18"), as.Date("2023-04-25"), as.Date("2023-05-02")),
               date_labels = "%b %d") +
  labs(x = "Date", y = "Prob. of detection") +
  scale_color_manual(values = col2[c(1,2)]) +
  scale_fill_manual(values = col2[c(1,2)]) +
  theme_bw() + mythemes

# top model effects, univariate predictions
occu.plots.univar <- list()
for(variable in c("Method", "Date")) {
  for(spp in names(occu_mod_outputs)) {
    occu.plots.univar[[variable]][[spp]] <- 
      plotEffectsData(occu_mod_outputs[[spp]]$mods$`Method + Date`,
                      type = "det", 
                      covariate = variable) %>%
      mutate(Species = spp)
  }
  # combine
  occu.plots.univar[[variable]][["Data"]] <- bind_rows(occu.plots.univar[[variable]])
  occu.plots.univar[[variable]][["Data"]]$Species <- 
    factor(occu.plots.univar[[variable]][["Data"]]$Species, 
           levels = c("LISY", "AMMA", "AMBY", "AMLA"))
  # make plots
  if(variable == "Method") {
    occu.plots.univar[[variable]][["Data"]]$Method <- 
      factor(occu.plots.univar[[variable]][["Data"]]$covariateValue, 
             levels = c("Trapping", "eDNA"))
    occu.plots.univar[[variable]][["Plot"]] <-
      ggplot(occu.plots.univar[[variable]][["Data"]], 
             aes(x = Method, y = Predicted, shape = Method, fill = Method)) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
      geom_point(size = 3) +
      scale_fill_manual(values = col2[1:2]) +
      scale_shape_manual(values = c(21,24)) +
      facet_wrap(~Species, labeller = as_labeller(fig.lets.labs)) +
      labs(x = variable, y = "Probability of detection") +
      theme_bw() + mythemes
  } else {
    occu.plots.univar[[variable]][["Plot"]] <-
      ggplot(occu.plots.univar[[variable]][["Data"]], aes(x = covariateValue, y = Predicted)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      facet_wrap(~Species, labeller = as_labeller(fig.lets.labs)) +
      labs(x = variable, y = "Probability of detection") +
      theme_bw() + mythemes
  }
}

# print plots
occu.plots.univar$Method$Plot
occu.plots.univar$Date$Plot

# get new labels for species for plotting
spp_names <- list(
  'AMMA'="Spotted salamander",
  'AMLA'="Blue-spotted salamander",
  'AMBY'="Unisexual salamander",
  'LISY'="Wood frog"
)

# make function for overriding labels
spp_labeller <- function(variable,value){
  return(spp_names[value])
}

# update figure without editing data
amph.method.det <- occu.plots.univar$Method$Plot + facet_wrap(~Species, labeller=spp_labeller)
amph.method.det <- occu.plots.univar$Method$Plot
amph.method.det
ggsave("figures/trapping.validation.jpg", 
       width = 9, height = 7, units = "in", dpi = 300)


# build cumulative detection functions from additive models

# From a probability of detection, the rarefaction curve is just
# a cumulative probability from a series of coin flips
#		y = 1 - (1-p)^x, 
#		where y is the cumulative probability of detection
#					x is the number of trials
#					p is the probability of detection on any one trial

# take average probabilities for each species across sites to build species-level rarefaction curves
## NOT ACCOUNTING FOR EFFECT OF DATE!
curve_data <- data.frame()
for(spp in names(occu_mod_outputs)) {
  x_vals <- 1:10
  for(method in c("Trapping", "eDNA")) {
    y_vals <- 1 - (1 - occu.plots.univar$Method$Data[
      occu.plots.univar$Method$Data$Method == method &
      occu.plots.univar$Method$Data$Species == spp,]$Predicted)^x_vals
    ymin_vals <- 1 - (1 - occu.plots.univar$Method$Data[
      occu.plots.univar$Method$Data$Method == method &
        occu.plots.univar$Method$Data$Species == spp,]$lower)^x_vals
    ymax_vals <- 1 - (1 - occu.plots.univar$Method$Data[
      occu.plots.univar$Method$Data$Method == method &
        occu.plots.univar$Method$Data$Species == spp,]$upper)^x_vals
    temp_data <- data.frame(
      Method = method,
      Species = spp,
      x = x_vals,
      y = y_vals,
      ymin = ymin_vals,
      ymax = ymax_vals
    )
    curve_data <- rbind(curve_data, temp_data)
  }
}

curve_data$Method <- factor(curve_data$Method, levels = c("Trapping", "eDNA"))
curve_data$Species <- factor(curve_data$Species, levels = c("LISY", "AMMA", "AMBY", "AMLA"))

(amph.method.det.curve <- ggplot(curve_data, aes(x = x, y = y, shape = Method)) +
  geom_line(aes(color = Method)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = Method), width = 0, position = position_dodge(0.4)) +
  geom_point(aes(color = Method, fill = Method), size = 2, position = position_dodge(0.4)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey") +
  scale_color_manual(values = col2[1:2]) +
  scale_fill_manual(values = col2[1:2]) +
  scale_shape_manual(values = c(21,24)) +
  scale_x_continuous(breaks = c(1:10)) +
  facet_wrap(~Species, nrow = 4) +
  labs(x = "Sampling units", y = "Cumulative detection prob.") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        panel.background = element_blank(),
        strip.background = element_blank()) + mythemes)

amph.method.det.curve + facet_wrap(~Species, labeller = as_labeller(fig.lets.labs))
ggsave("figures/trapping.validation.curve.jpg", 
       width = 9, height = 7, units = "in", dpi = 300)

fig.lets.labs.2col <- c("E    Wood frog", "F    Spotted Salamander", 
                        "G    Unisexual Salamander", "H    Blue-spotted Salamander")
names(fig.lets.labs.2col) <- names(fig.lets)
ggarrange(amph.method.det + 
            facet_wrap(~Species, labeller = as_labeller(fig.lets.labs), nrow = 4) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  strip.background = element_blank()),
          amph.method.det.curve + facet_wrap(~Species, labeller = as_labeller(fig.lets.labs.2col), nrow = 4),
          common.legend = T, legend = "right")
ggsave("figures/trapping.validation.curve.composite.jpg", 
       width = 9, height = 7, units = "in", dpi = 300)


### 7a-3. eDNA vs cameras ----
#### load and combine eDNA and camera datasets ----
# import camera data
cameras <- read.csv("data/camera_data_2025-11-12.csv")
cameras$Species.binom <- paste(cameras$Genus, cameras$Species.1)

# reformat and correct some species
camera.spp <- unique(cameras$Species.Common.Name)
camera.spp <- camera.spp[camera.spp %notin% c("", "Unknown", "Unknown bird", "Unknown mammal", 
                                              "Unknown amphibian", "Unknown reptile", 
                                              "Other (comment)")]

library(data.table)
# subset eDNA from camera sites, get average reads for each taxa frome each site on each occasion
cameras.edna <- subset_samples(ps.tax.samps.spp, Cameras == "Yes") %>% psmelt()
cameras.edna.agg <- cameras.edna %>% 
  mutate(Site = ifelse(Site == "RW-Stream", "RW-S", Site),
         Date_simple = format(Date, "%Y-%m-%d")) %>%
  group_by(Site, Date_simple, OTU) %>%
  dplyr::summarize(Mean_reads = mean(Abundance),
            DateTimeSurvey = max(Date)) %>% ungroup() %>%
  dplyr::group_by(Site) %>%
  mutate(Occasion = frank(Date_simple, ties.method = "dense")) %>% ungroup() %>%
  left_join(dplyr::select(cameras.edna, OTU, Class, Order, Family, Genus, Species),
            by = "OTU", multiple = "first")

# subset camera data that overlaps with eDNA surveys
cameras.validation <- cameras %>%
  filter(Date.Time <= max(cameras.edna$Date))

# extract bird and mammal data from eDNA 
cameras.edna.birds.mammals <- 
  filter(cameras.edna.agg, Class %in% c("Aves", "Mammalia"))

# make some post-hoc species assignments to both datasets
cameras.validation <- cameras.validation %>%
  mutate(Species.binom = case_when(Species.binom == "Canis familiaris" ~ "Canis lupus familiaris",
                                   Species.binom == "Felix catus" ~ "Felis catus",
                                   TRUE ~ Species.binom))

cameras.edna.birds.mammals <- cameras.edna.birds.mammals %>%
  mutate(Species = case_when(Species == "Canis lupus" ~ "Canis lupus familiaris",
                             TRUE ~ Species))

# combine monthly datasets of reads and camera trap detections ----
## assign site-specific occasion variable to camera dataset
camera.site.deployments <- as.POSIXct(c("2023-04-06 01:00:00",
                                        "2023-04-05 01:00:00",
                                        "2023-04-13 01:00:00",
                                        "2023-03-31 01:00:00",
                                        "2023-04-13 01:00:00",
                                        "2023-04-05 01:00:00",
                                        "2023-04-13 01:00:00",
                                        "2023-04-05 01:00:00"),
                                      format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
names(camera.site.deployments) <- c("FS", "SB", "HP", "RW", "RW-S", "UFTS1", "UFTS5", "UFTS42")

cameras.validation$Occasion <- 1
cameras.validation$Occasion.dur <- 1
for(site in unique(cameras.edna.agg$Site)) {
  occs <- seq(min(cameras.edna.agg$Occasion[cameras.edna.agg$Site == site]),
              max(cameras.edna.agg$Occasion[cameras.edna.agg$Site == site]), 1)
  dates <- unique(cameras.edna.agg$DateTimeSurvey[cameras.edna.agg$Site == site])
  occ.dat <- data.frame(occ = occs,
                        date1 = dates,
                        date0 = lag(dates))
  occ.dat$date0[is.na(occ.dat$date0)] <- camera.site.deployments[site]
  occ.dat$dur <- as.integer(occ.dat$date1 - occ.dat$date0)
  occ.find <- function(x) occ.dat$occ[x >= occ.dat$date0 & x < occ.dat$date1]
  cameras.validation$Occasion[cameras.validation$Site.Name == site] <-
    unlist(sapply(cameras.validation$Date.Time[cameras.validation$Site.Name == site], FUN = occ.find))
  for(occ in 1:4) {
    cameras.validation$Occasion.dur[cameras.validation$Site.Name == site &
                                      cameras.validation$Occasion == occ] <-
      occ.dat$dur[occ.dat$occ == occ]
  }
}

## calculate relative abundance of detections for species in each site x occasion
cameras.validation.agg <- cameras.validation %>%
  filter(Family != "" & Attention.Needed.Flag == "false") %>% # omit empty or ambiguous images
  mutate(Date = date(ymd_hms(Date.Time))) %>% # create day variable for daily counts
  group_by(Site.Name, Occasion, Class, Order, Family, Genus, Species.binom) %>%
  dplyr::summarize(Total.images = n()/unique(Occasion.dur),
            Total.images.water = sum(Contacting.water == "Yes")/unique(Occasion.dur),
            Total.indiv = sum(Sighting.Quantity)/unique(Occasion.dur),
            Total.indiv.water = sum(Sighting.Quantity[Contacting.water == "Yes"])/unique(Occasion.dur),
            Total.days = length(unique(Date))/unique(Occasion.dur),
            Total.days.water = length(unique(Date[Contacting.water == "Yes"]))/unique(Occasion.dur),
            .groups = "rowwise") %>%
  ungroup() %>%
  dplyr::rename(Site = Site.Name,
                Species = Species.binom)

## combine monthly eDNA and camera data into new data frame
cameras.edna.combined <- full_join(cameras.edna.agg, cameras.validation.agg,
                                   by = c("Site", "Occasion", "Class",
                                          "Order", "Family", "Genus", "Species")) %>%
  # fill in zeroes for either method missing data
  mutate(Mean_reads = ifelse(is.na(Mean_reads), 0, Mean_reads),
         Total.images = ifelse(is.na(Total.images), 0, Total.images),
         Total.images.water = ifelse(is.na(Total.images.water), 0, Total.images.water),
         Total.indiv = ifelse(is.na(Total.indiv), 0, Total.indiv),
         Total.indiv.water = ifelse(is.na(Total.indiv.water), 0, Total.indiv.water),
         Total.days = ifelse(is.na(Total.days), 0, Total.days),
         Total.days.water = ifelse(is.na(Total.days.water), 0, Total.days.water)) %>%
  # filter only birds and mammals with detections from either method 
  dplyr::filter(Class %in% c("Aves", "Mammalia") &
                  (Mean_reads > 0 | Total.images > 0)) %>%
  # remove non-native species detected from eDNA
  dplyr::filter(Species %notin% c("Cairina moschata", "Aplonis panayensis", "Carpodacus erythrinus",
                                  "Chloris sinica", "Mus spretus", "Rattus norvegicus")) %>%
  # remove native species without reference sequences for eDNA
  dplyr::filter(Species %notin% c("Astur cooperii", "Baeolophus bicolor","Empidonac traillii")) %>%
  # rename / group species to be shared across methods
  mutate(Species = ifelse(Species == "Leuconotopicus villosus", "Picoides villosus", Species),
         Species = ifelse(Species == "Canis lupus", "Canis lupus familiaris", Species),
         Species = ifelse(Species %in% c("Peromyscus maniculatus", 
                                         "Peromyscus  maniculatus", 
                                         "Peromyscus leucopus"), 
                          "Peromyscus sp.", Species))

# create table of eDNA species detections by number of sites detected
cameras.edna.site.tab <- cameras.edna.birds.mammals %>%
  filter(Mean_reads > 0 & 
           Species %notin% c("Cairina moschata", "Aplonis panayensis", "Carpodacus erythrinus",
                             "Chloris sinica", "Mus spretus", "Rattus norvegicus",
                             "Astur cooperii", "Baeolophus bicolor","Empidonac traillii")) %>%
  group_by(Class, Species) %>%
  dplyr::summarize(n.sites.eDNA = length(unique(Site)),
            prop.sites.eDNA = n.sites.eDNA / 8) %>% ungroup()

# create table of camera species detections by number of sites detected
cameras.site.tab <- cameras.validation %>%
  filter(Attention.Needed.Flag == "false" &
           Class %in% c("Aves", "Mammalia")) %>%
  group_by(Class, Species.binom) %>%
  dplyr::summarize(n.sites.cameras = length(unique(Site.Name)),
            prop.sites.cam = n.sites.cameras / 8) %>% ungroup() %>%
  dplyr::rename(Species = Species.binom) %>%
  mutate(Species = ifelse(Species == "Leuconotopicus villosus", "Picoides villosus", Species))
  
# combine tables into one
cameras.contr.tab <- cameras.site.tab %>%
  full_join(cameras.edna.site.tab, by = c("Class", "Species")) %>%
  replace(is.na(.), 0) %>%
  arrange(Class, Species)

# quick correlation of eDNA / camera species prevalence (total mess!)
ggplot(cameras.contr.tab, aes(prop.sites.cam, prop.sites.eDNA, color = Class)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  ylim(c(0,1)) +
  theme_classic() + mythemes

#### NMDS and PERMANOVA to compare community data structure between methods ----
#### Mantel tests to compare presence/absence data at the camera sites ----

# monthly composition comparison

# overall (site-level) composition comparison

#### Regression relating reads vs camera detections ----
# simple correlation between reads and camera detection frequencies
## standardize eDNA reads and images for each species

## by class
ggplot(cameras.edna.combined, aes(scale(Total.images), scale(Mean_reads), color = Class)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  theme_classic() + mythemes +
  labs(x = "images/day (z-score)", y = "eDNA reads (z-score)")

## by order
library(ggh4x)
ggplot(cameras.edna.combined, aes(scale(Total.images),  scale(Mean_reads))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", col = "grey") +
  facet_wrap(~Class + Order) +
  theme_bw() + mythemes +
  theme(strip.text.x = element_text(size = 11)) +
  labs(x = "images/day (z-score)", y = "eDNA reads (z-score)")

# remove species with <3 records
library(scales)
cameras.edna.combined %>%
  group_by(Species) %>%
  filter(n() > 3) %>%
  ggplot(aes(Total.images,  scale(Mean_reads))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", col = "grey") +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))) +
  facet_wrap(~Class + Order, scales = "free_x") +
  theme_bw() + mythemes +
  theme(strip.text.x = element_text(size = 11)) +
  labs(x = "images/day", y = "eDNA reads (z-score)")

## with logistic curves (probability of detecting species with eDNA)
cameras.edna.combined %>%
  group_by(Species) %>%
  filter(n() > 3) %>%
  ggplot(aes(Total.images,  makebinary(Mean_reads))) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  facet_wrap(~Class + Order, scales = "free_x") +
  theme_bw() + mythemes +
  theme(strip.text.x = element_text(size = 11)) +
  labs(x = "images/day", y = "eDNA detection")

## only including observations with image data
cameras.edna.combined %>%
  filter(Total.images > 0) %>%
  ggplot(aes(scale(Total.images),  scale(Mean_reads))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", col = "grey") +
  facet_wrap(~Class + Order, scales = "free") +
  theme_bw() + mythemes +
  theme(strip.text.x = element_text(size = 11)) +
  labs(x = "images/day (z-score)", y = "eDNA reads (z-score)")


### 7b. Landscape eDNA analysis

# Comparing eDNA biodiversity in different wetlands with different hydroperiods
# Drivers of community structure, include other environmental covariates
# Ambystoma complex co-occurrence analysis
# Community occupancy analysis?
## a la https://esj-journals.onlinelibrary.wiley.com/doi/full/10.1002/1438-390X.12219

# Misc. species summaries
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

#### 7b. Overall diversity comparison ----
# 7b-1. Venn diagrams and community composition ----
# create empty data frame to populate
## get total number of amphibian species detected with either method
methods.venn <- data.frame(Class = c("Amphibian", "Mammal", "Bird"),
                           Trad = NA,
                           eDNA = NA,
                           Both = NA)
amphibian.trad <- c("Lithobates sylvaticus", "Ambystoma barbouri", "Ambystoma maculatum", 
                        "Ambystoma laterale", "Pseudacris crucifer", "Lithobates septentrionalis",
                        "Lithobates clamitans", "Notophthalmus viridescens")
amphibian.eDNA.melt <- ps.tax.samps.spp %>%
  psmelt() %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(Class == "Amphibia" &
           # remove non-native species
           Species %notin% c("Ambystoma opacum", "Ambystoma jeffersonianum") &
           Date < "2023-05-05" &
           Trapping == "Yes") %>%
  mutate(Species = ifelse(Species == "Lithobates catesbeianus", "Lithobates clamitans", Species))
amphibian.eDNA <- unique(amphibian.eDNA.melt$Species)

## populate in empty summary data frame for amphibians
methods.venn[methods.venn$Class == "Amphibian", "Trad"] <- 
  sum(amphibian.trad %notin% amphibian.eDNA)
methods.venn[methods.venn$Class == "Amphibian", "eDNA"] <- 
  sum(amphibian.eDNA %notin% amphibian.trad)
methods.venn[methods.venn$Class == "Amphibian", "Both"] <- 
  sum(amphibian.trad %in% amphibian.eDNA)

## populate in empty summary data frame for mammal
cameras.edna.combined.venn <- full_join(cameras.edna.agg, cameras.validation.agg,
                                        by = c("Site", "Occasion", "Class",
                                               "Order", "Family", "Genus", "Species")) %>%
  # fill in zeroes for either method missing data
  mutate(Mean_reads = ifelse(is.na(Mean_reads), 0, Mean_reads),
         Total.images = ifelse(is.na(Total.images), 0, Total.images),
         Total.images.water = ifelse(is.na(Total.images.water), 0, Total.images.water),
         Total.indiv = ifelse(is.na(Total.indiv), 0, Total.indiv),
         Total.indiv.water = ifelse(is.na(Total.indiv.water), 0, Total.indiv.water),
         Total.days = ifelse(is.na(Total.days), 0, Total.days),
         Total.days.water = ifelse(is.na(Total.days.water), 0, Total.days.water)) %>%
  # filter only birds and mammals with detections from either method 
  dplyr::filter(Class %in% c("Aves", "Mammalia") &
                  (Mean_reads > 0 | Total.images > 0)) %>%
  # remove non-native species detected from eDNA
  dplyr::filter(Species %notin% c("Cairina moschata", "Aplonis panayensis", "Carpodacus erythrinus",
                                  "Chloris sinica", "Mus spretus", "Rattus norvegicus")) %>%
  mutate(Species = ifelse(Species == "Leuconotopicus villosus", "Picoides villosus", Species),
         Species = ifelse(Species == "Canis lupus", "Canis lupus familiaris", Species),
         Species = ifelse(Species %in% c("Peromyscus maniculatus", "Peromyscus  maniculatus", "Peromyscus leucopus"), "Peromyscus sp.", Species))

mammal.trad <- unique(cameras.edna.combined.venn[cameras.edna.combined.venn$Class == "Mammalia" & 
                                              cameras.edna.combined.venn$Total.images > 0,]$Species)
mammal.eDNA <- unique(cameras.edna.combined.venn[cameras.edna.combined.venn$Class == "Mammalia" & 
                                              cameras.edna.combined.venn$Mean_reads > 0,]$Species)
methods.venn[methods.venn$Class == "Mammal", "Trad"] <- 
  sum(mammal.trad %notin% mammal.eDNA)
methods.venn[methods.venn$Class == "Mammal", "eDNA"] <- 
  sum(mammal.eDNA %notin% mammal.trad)
methods.venn[methods.venn$Class == "Mammal", "Both"] <- 
  sum(mammal.trad %in% mammal.eDNA)

## populate in empty summary data frame for bird
bird.trad <- unique(cameras.edna.combined.venn[cameras.edna.combined.venn$Class == "Aves" & 
                                              cameras.edna.combined.venn$Total.images > 0,]$Species)
bird.eDNA <- unique(cameras.edna.combined.venn[cameras.edna.combined.venn$Class == "Aves" & 
                                              cameras.edna.combined.venn$Mean_reads > 0,]$Species)
methods.venn[methods.venn$Class == "Bird", "Trad"] <- 
  sum(bird.trad %notin% bird.eDNA)
methods.venn[methods.venn$Class == "Bird", "eDNA"] <- 
  sum(bird.eDNA %notin% bird.trad)
methods.venn[methods.venn$Class == "Bird", "Both"] <- 
  sum(bird.trad %in% bird.eDNA)

## add a total column
methods.venn$Total <- methods.venn$Trad + methods.venn$eDNA + methods.venn$Both

## venn diagrams
library(ggvenn)
amph <- list(Traditional = amphibian.trad, eDNA = amphibian.eDNA)
bird <- list(Traditional = bird.trad, eDNA = bird.eDNA)
mammal <- list(Traditional = mammal.trad, eDNA = mammal.eDNA)

# crete separate legend to attach to composite fig
dummy_legend <- 
  ggplot(data.frame(Method = factor(c("Traditional", "eDNA"),
                                    levels = c("Traditional", "eDNA")),
                    y = c(1, 1)),
         aes(x = Method, y = y, fill = Method)) +
  geom_col() +
  scale_fill_manual(
    values = c("dodgerblue", "purple"),
    name = "Method") +
  theme_void() +
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.box.spacing = unit(0, "pt"))+
  guides(fill=guide_legend(ncol= 1))

legend.grob <- cowplot::get_legend(dummy_legend)

amph.venn <- ggvenn(amph, auto_scale = TRUE, fill_color = c("dodgerblue", "purple"), 
                    fill_alpha = 0.4, show_percentage = FALSE, set_name_color = "transparent") +
  ggtitle("Amphibians \n (S = 12)") +
  theme(plot.title.position = 'panel', 
        plot.title = element_text(hjust = 0.5, vjust = -0.4, size = 16),
        plot.margin=unit(c(0,0,0,0), "pt"))+
  coord_fixed(expand = FALSE)

bird.venn <- ggvenn(bird, auto_scale = TRUE, fill_color = c("dodgerblue", "purple"), 
                    fill_alpha = 0.4, show_percentage = FALSE, set_name_color = "transparent") +
  ggtitle("Birds \n (S = 31)") +
  theme(plot.title.position = 'panel', 
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.margin=unit(c(0,0,0,0), "pt"))+
  coord_fixed(expand = FALSE)

mammal.venn <- ggvenn(mammal, auto_scale = TRUE, fill_color = c("dodgerblue", "purple"), 
                      fill_alpha = 0.4, show_percentage = FALSE, set_name_color = "transparent") +
  ggtitle("Mammals \n (S = 25)") +
  theme(plot.title.position = 'panel', 
        plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 16),
        plot.margin=unit(c(0,0,0,0), "pt"))+
  coord_fixed(expand = FALSE)

venn.comp.1r <- ggarrange(amph.venn, bird.venn, mammal.venn, legend.grob,
                       ncol = 4, labels = c("A", "B", "C", ""),
                       label.y = 0.88)
venn.comp.1r

venn.comp <- ggarrange(amph.venn, bird.venn, NULL, mammal.venn, legend.grob,
                       nrow = 2, ncol = 2, labels = c("A", "B", "C", ""))
venn.comp

ggsave("figures/trapping.validation.venn.jpg", 
       width = 6, height = 6, units = "in", dpi = 300)

# combine venn diagrams with 12S composition in seasonal versus semi/permanent wetlands
rar.melt$Wetland_type2 <- rar.melt$Wetland_type
rar.melt$Wetland_type2[rar.melt$Wetland_type2 == "Pool"] <- "Seasonal"
rar.melt$Wetland_type2 <- factor(rar.melt$Wetland_type2, 
                                 levels = c("Seasonal", "Semipermanent", "Permanent"))
rar.melt <- rar.melt[!is.na(rar.melt$Class),]

class.comp.type <- rar.melt %>% group_by(Wetland_type2, Class) %>%
  dplyr::summarize(Reads = sum(Reads)) %>%
  ggplot(aes(factor(Wetland_type2, levels = rev(c("Seasonal", "Semipermanent", "Permanent"))), 
             Reads, fill = Class)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.75) +
  labs(x = NULL, y = "Prop. reads") +
  scale_x_discrete(drop=FALSE) +
  coord_flip() +
  scale_fill_manual(values = brewer.pal(n = 6, name = "Accent")) +
  theme_void() + mythemes +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        axis.text.y = element_text(hjust = 1),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "top") +
  guides(fill=guide_legend(nrow = 1, reverse = TRUE))

class.comp.type

ggarrange(venn.comp.1r,
          class.comp.type, 
          nrow = 2, heights = c(0.75,0.5),
          labels = c("", "D"))

ggsave("figures/trapping.validation.eDNA_composition.1row.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

# 7b-2. Relating naive occupancy between methods ----

# amphibian data
## create validation dataset
amphibian.eDNA.val <- 
  ps.tax.samps.spp %>%
  psmelt() %>%
  filter(Class == "Amphibia" &
           # remove non-native species
           Species %notin% c("Ambystoma opacum", "Ambystoma jeffersonianum") &
           Trapping == "Yes" | Egg_mass_survey == "Yes") %>%
  mutate(Species = ifelse(Species == "Lithobates catesbeianus", "Lithobates clamitans", Species),
         Site = ifelse(Site == "RW-Stream", "RW-S", Site))

# fix spelling of RW-Stream to be same as trapping dataset
metadata.sites2 <- metadata.sites %>%
  mutate(Site = ifelse(SiteID == "RW-Stream", "RW-S", SiteID))

## populate table for naive occupancy
amph.naive.occ <- data.frame(Species = unique(c(amph$Traditional, amph$eDNA)),
                             nsites.trad = NA, nsites.eDNA = NA, nsites.both = NA, 
                             nsites.tot = NA, nsites.surv = NA)
amph.naive.occ$Code <- c("LISY", "AMBY", "AMMA", "AMLA", "PSCR", "LISE", "LICL", "NOVI",
                         "PLCI", "HESC", "LIPI", "EUBI")
## loop through amphibian species and datasets
for(sp in amph.naive.occ$Code) {
  if(sp == "AMMA") {
    # list of validation sites
    sites <- metadata.sites2$Site[metadata.sites2$Trapping == "Yes" | 
                                    metadata.sites2$Egg_mass_survey == "Yes"]
    amph.naive.occ$nsites.surv[amph.naive.occ$Code == sp] <- length(sites)
    
    # sites where AMMA found with traditional methods (trapping or egg mass counts)
    sites.trad <- na.omit(unique(c(as.character(traps$Site[traps$Species == sp]), 
                                   as.character(metadata.sites2$Site[metadata.sites2$AMMA_eggs > 0]))))
    amph.naive.occ$nsites.trad[amph.naive.occ$Code == sp] <- length(sites.trad)
    
    # validation sites where AMMA found with eDNA
    sites.eDNA <- unique(amphibian.eDNA.val$Site[
      amphibian.eDNA.val$Species == amph.naive.occ$Species[amph.naive.occ$Code == sp] &
        amphibian.eDNA.val$Abundance > 0])
    
    amph.naive.occ$nsites.eDNA[amph.naive.occ$Code == sp] <- length(sites.eDNA)
    
    # sites where AMMA found with both methods
    sites.both <- sites.trad[sites.trad %in% sites.eDNA]
    amph.naive.occ$nsites.both[amph.naive.occ$Code == sp] <- sum(sites.trad %in% sites.eDNA)
    
    # sites where AMMA found with either method
    amph.naive.occ$nsites.tot[amph.naive.occ$Code == sp] <- length(unique(c(sites.trad, sites.eDNA)))
  } else if(sp == "LISY") {
    # list of validation sites
    sites <- metadata.sites2$Site[metadata.sites2$Trapping == "Yes"]
    amph.naive.occ$nsites.surv[amph.naive.occ$Code == sp] <- length(sites)
    
    # sites where LISY found with traditional methods (trapping or egg mass counts)
    sites.trad <- na.omit(unique(c(as.character(traps$Site[traps$Species == sp]), 
                                   as.character(metadata.sites2$Site[metadata.sites2$LISY_eggs > 0]))))
    amph.naive.occ$nsites.trad[amph.naive.occ$Code == sp] <- length(sites.trad)
    
    # validation sites where AMMA found with eDNA
    sites.eDNA <- unique(amphibian.eDNA.val$Site[
      amphibian.eDNA.val$Species == amph.naive.occ$Species[amph.naive.occ$Code == sp] &
        amphibian.eDNA.val$Abundance > 0 & amphibian.eDNA.val$Trapping == "Yes"])
      
      amph.naive.occ$nsites.eDNA[amph.naive.occ$Code == sp] <- length(sites.eDNA)
      
      # sites where LISY found with both methods
      sites.both <- sites.trad[sites.trad %in% sites.eDNA]
      amph.naive.occ$nsites.both[amph.naive.occ$Code == sp] <- sum(sites.trad %in% sites.eDNA)
      
      # sites where LISY found with either method
      amph.naive.occ$nsites.tot[amph.naive.occ$Code == sp] <- length(unique(c(sites.trad, sites.eDNA)))
  } else {
    # for all other species
    # list of validation sites
    sites <- metadata.sites2$Site[metadata.sites2$Trapping == "Yes"]
    amph.naive.occ$nsites.surv[amph.naive.occ$Code == sp] <- length(sites)
    
    # sites where sp found with traditional methods (trapping or egg mass counts)
    sites.trad <- na.omit(unique(as.character(traps$Site[traps$Species == sp])))
    amph.naive.occ$nsites.trad[amph.naive.occ$Code == sp] <- length(sites.trad)
    
    # validation sites where sp found with eDNA
    sites.eDNA <- unique(amphibian.eDNA.val$Site[
      amphibian.eDNA.val$Species == amph.naive.occ$Species[amph.naive.occ$Code == sp] &
        amphibian.eDNA.val$Abundance > 0 & amphibian.eDNA.val$Trapping == "Yes"])
    
    amph.naive.occ$nsites.eDNA[amph.naive.occ$Code == sp] <- length(sites.eDNA)
    
    # sites where sp found with both methods
    sites.both <- sites.trad[sites.trad %in% sites.eDNA]
    amph.naive.occ$nsites.both[amph.naive.occ$Code == sp] <- sum(sites.trad %in% sites.eDNA)
    
    # sites where sp found with either method
    amph.naive.occ$nsites.tot[amph.naive.occ$Code == sp] <- length(unique(c(sites.trad, sites.eDNA)))
  }
}

# bird and mammal data
## populate table for naive occupancy
birdmamm.naive.occ <- cameras.edna.combined.venn %>%
  dplyr::group_by(Species, Site) %>%
  dplyr::summarize(eDNA = makebinary(mean(Mean_reads)),
                   cameras = makebinary(mean(Total.images))) %>%
  ungroup() %>%
  group_by(Species) %>%
  dplyr::summarize(nsites.trad = length(unique(Site[cameras == 1])),
            nsites.eDNA = length(unique(Site[eDNA == 1])),
            nsites.both = length(unique(Site[cameras == 1 & eDNA == 1])),
            nsites.tot = length(unique(Site[cameras == 1 | eDNA == 1])),
            nsites.surv = 8) %>%
  ungroup()

# combine with bird and mammal data
naive.occ.data <- data.frame(
  Species = c(amph.naive.occ$Species, cameras.contr.tab$Species),
  nsites_trad = c(amph.naive.occ$nsites.trad, cameras.contr.tab$n.sites.cameras),
  nsites_eDNA = c(amph.naive.occ$nsites.eDNA, cameras.contr.tab$n.sites.eDNA),
  nsites_surv = c(amph.naive.occ$nsites.surv, rep(8, nrow(cameras.contr.tab))),
  Class = c(rep("Amphibia", nrow(amph.naive.occ)),
            cameras.contr.tab$Class))
naive.occ.data$nsites_trad_prop = naive.occ.data$nsites_trad/naive.occ.data$nsites_surv
naive.occ.data$nsites_eDNA_prop = naive.occ.data$nsites_eDNA/naive.occ.data$nsites_surv

# remove species not found in both surveys
naive.occ.data.both <- naive.occ.data[naive.occ.data$nsites_trad > 0 & 
                                        naive.occ.data$nsites_eDNA > 0,]

occ.mod.list <- list()
occ.mod.tab <- data.frame(Class = c("Amphibia", "Aves", "Mammalia"),
                             Intercept = NA,
                             slope = NA,
                             SE = NA,
                             p_val = NA,
                             R2 = NA)
# run models
for(group in c("Amphibia", "Aves", "Mammalia")) {
  mod <- glm(nsites_eDNA ~ nsites_trad,
            naive.occ.data.both[naive.occ.data.both$Class == group,],
            family = poisson)
  # check for overdispersion
  overdisp <- 1 - pchisq(deviance(mod), df.residual(mod))
  print(paste(group, "overdispersion p =", round(overdisp,3), "(Pearson)"))
  # if overdispersion, use poisson
  if(overdisp < 0.05) {
    mod <- glm(nsites_eDNA ~ nsites_trad,
               naive.occ.data.both[naive.occ.data.both$Class == group,],
               family = quasipoisson)
  }
  occ.mod.list[[group]] <- mod
  occ.mod.tab[occ.mod.tab$Class == group,]$Intercept <- coef(mod)[[1]]
  occ.mod.tab[occ.mod.tab$Class == group,]$slope <- coef(mod)[[2]]
  occ.mod.tab[occ.mod.tab$Class == group,]$SE <- summary(mod)$coefficients[2,2]
  occ.mod.tab[occ.mod.tab$Class == group,]$p_val <- summary(mod)$coefficients[2,4]
  occ.mod.tab[occ.mod.tab$Class == group,]$R2 <- as.numeric(r2(mod))
}

# predict models
occ.plot.data <- NULL
vert.classes <- c("Amphibians" = "Amphibia", "Birds" = "Aves", "Mammals" = "Mammalia")
for(group in c("Amphibia", "Aves", "Mammalia")) {
  newdat <- data.frame(
    nsites_trad = seq(1,max(naive.occ.data.both$nsites_trad[naive.occ.data.both$Class == group]),
                      length.out = 30))
  # predict on link scale to get appropriate CIs
  preds <- as.data.frame(predict(mod.list[[group]],
                                 newdata = newdat,
                                 type = "link", 
                                 se.fit = TRUE))
  critval <- 1.96 ## approx 95% CI
  upr <- occ.mod.list[[group]]$family$linkinv(preds$fit + (critval * preds$se.fit))
  lwr <- occ.mod.list[[group]]$family$linkinv(preds$fit - (critval * preds$se.fit))
  fit <- occ.mod.list[[group]]$family$linkinv(preds$fit) 
  preds <- cbind(newdat,fit,lwr,upr)
  preds$Group = names(vert.classes[vert.classes == group])
  occ.plot.data <- bind_rows(occ.plot.data, preds)
}
  
# plot relationship
fig.lets.labs.class <- c("Amphibians" = "A      Amphibians", 
                         "Birds" = "B      Birds", 
                         "Mammals" = "C      Mammals")
naive.occ.data.both <- naive.occ.data.both %>%
  mutate(Group = ifelse(Class == "Amphibia", "Amphibians", 
                        ifelse(Class == "Aves", "Birds", "Mammals")))
(naive.occ.plot <- ggplot(occ.plot.data, aes(nsites_trad, fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_point(data = naive.occ.data.both, aes(nsites_trad, nsites_eDNA),size = 3, fill = "black") +
    facet_wrap(~Group, nrow = 1, scale = "free", labeller = as_labeller(fig.lets.labs.class)) +
    theme_bw() + mythemes +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          strip.background = element_blank()) +
  labs(x = "n sites (traditional)", y = "n sites (eDNA)"))

ggsave("figures/trapping.validation.naive.occ.jpg", 
       width = 8, height = 5, units = "in", dpi = 300)

# 7b-3. Comparing species detection (naive richness) with either method in poisson GLM ----
## get number of species detected from sites with trad vs eDNA (only from concurrent surveys!)
### amphibians
#### trapping
amph.rich.trad <- traps %>%
  filter(!is.na(Species %notin% c("LITH", "FISH", "CRAY"))) %>%
  group_by(Site) %>%
  # get number of species detected in each site with traditional methods
  dplyr::summarize(Trapping = length(unique(Species))) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup()
#### eDNA
amph.rich.edna <- ps.tax.samps.spp %>%
  psmelt() %>%
  # only include eDNA data from trapping sites taken concurrently with trapping
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
         Site = ifelse(Site == "RW-Stream", "RW-S", Site)) %>%
  filter(Class == "Amphibia" &
           Date < "2023-05-05" &
           Trapping == "Yes" &
           Abundance > 0) %>%
  group_by(Site) %>%
  # get number of species detected in each site with eDNA
  dplyr::summarize(eDNA = length(na.omit(unique(Species)))) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  ungroup()
##### combine together
amph.rich.both <- left_join(amph.rich.trad, amph.rich.edna, by = "Site")
amph.rich.both$Class <- "Amphibians"

### birds
bird.rich.trad <- cameras.edna.combined %>%
  filter(Total.images > 0 & Class == "Aves") %>%
  group_by(Site) %>%
  dplyr::summarize(Trapping = length(na.omit(unique(Species))))%>%
  ungroup()
bird.rich.edna <- cameras.edna.combined %>%
  filter(Mean_reads > 0 & Class == "Aves") %>%
  group_by(Site) %>%
  dplyr::summarize(eDNA = length(na.omit(unique(Species))))%>%
  ungroup()
bird.rich.both <- left_join(bird.rich.trad, bird.rich.edna, by = "Site")
bird.rich.both$Class <- "Birds"

### mammals
mamm.rich.trad <- cameras.edna.combined %>%
  filter(Total.images > 0 & Class == "Mammalia") %>%
  group_by(Site) %>%
  dplyr::summarize(Trapping = length(na.omit(unique(Species)))) %>%
  ungroup()
mamm.rich.edna <- cameras.edna.combined %>%
  filter(Mean_reads > 0 & Class == "Mammalia") %>%
  group_by(Site) %>%
  dplyr::summarize(eDNA = length(na.omit(unique(Species)))) %>%
  ungroup()
mamm.rich.both <- left_join(mamm.rich.trad, mamm.rich.edna, by = "Site")
mamm.rich.both$Class <- "Mammals"

## combine all data
verts.rich <- bind_rows(amph.rich.both,
                        bird.rich.both,
                        mamm.rich.both) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  reshape2::melt(id.vars = c("Site", "Class"), variable.name = "Method", value.name = "S") %>%
  mutate(Method = factor(Method, levels = c("Trapping", "eDNA")))
  

## run analysis (poisson GLM)
rich.mod.list <- list()
rich.mod.tab <- data.frame(Class = c("Amphibians", "Birds", "Mammals"),
                           Intercept = NA,
                           eDNA = NA,
                           SE = NA,
                           p_val = NA,
                           R2 = NA)
# run models
for(group in c("Amphibians", "Birds", "Mammals")) {
  mod <- glm(S ~ Method,
             verts.rich[verts.rich$Class == group,],
             family = poisson)
  # check for overdispersion
  overdisp <- 1 - pchisq(deviance(mod), df.residual(mod))
  print(paste(group, "overdispersion p =", round(overdisp,3), "(Pearson)"))
  # if overdispersion, use poisson
  if(overdisp < 0.05) {
    mod <- glm(S ~ Method,
               verts.rich[verts.rich$Class == group,],
               family = quasipoisson)
  }
  rich.mod.list[[group]] <- mod
  rich.mod.tab[rich.mod.tab$Class == group,]$Intercept <- coef(mod)[[1]]
  rich.mod.tab[rich.mod.tab$Class == group,]$eDNA <- coef(mod)[[2]]
  rich.mod.tab[rich.mod.tab$Class == group,]$SE <- summary(mod)$coefficients[2,2]
  rich.mod.tab[rich.mod.tab$Class == group,]$p_val <- summary(mod)$coefficients[2,4]
  rich.mod.tab[rich.mod.tab$Class == group,]$R2 <- as.numeric(r2(mod))
}
  
# predict models
rich.plot.data <- NULL
for(group in c("Amphibians", "Birds", "Mammals")) {
  newdat <- data.frame(
    Method = c("Trapping", "eDNA")
  )
  # predict on link scale to get appropriate CIs
  preds <- as.data.frame(predict(rich.mod.list[[group]],
                                 newdata = newdat,
                                 type = "link", 
                                 se.fit = TRUE))
  critval <- 1.96 ## approx 95% CI
  upr <- rich.mod.list[[group]]$family$linkinv(preds$fit + (critval * preds$se.fit))
  lwr <- rich.mod.list[[group]]$family$linkinv(preds$fit - (critval * preds$se.fit))
  fit <- rich.mod.list[[group]]$family$linkinv(preds$fit) 
  preds <- cbind(newdat,fit,lwr,upr)
  preds$Class <- group
  rich.plot.data <- bind_rows(rich.plot.data, preds)
}

fig.lets.labs.class.2col <- c("Amphibians" = "D     Amphibians", 
                         "Birds" = "E     Birds", 
                         "Mammals" = "F     Mammals")
# plot relationship
rich.plot.data$Method <- factor(rich.plot.data$Method, levels = c("Trapping", "eDNA"))
(naive.rich.plot <- ggplot(rich.plot.data, aes(Method, fit, shape = Method, fill = Method)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0) +
  geom_point(size = 3) +
  scale_fill_manual(values = col2[1:2]) +
  scale_shape_manual(values = c(21,24)) +
  facet_wrap(~Class, nrow = 1, scale = "free", labeller = as_labeller(fig.lets.labs.class.2col)) +
  theme_bw() + mythemes +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +
  labs(x = "Method", y = "Richness"))

ggsave("figures/trapping.validation.naive.rich.jpg", 
       width = 8, height = 5, units = "in", dpi = 300)

# combine naive occ and rich plots together
ggarrange(naive.occ.plot + theme(plot.margin = unit(c(0,1,0,1), 'lines')),
          NULL,
          naive.rich.plot + theme(plot.margin = unit(c(0,1,0,1), 'lines')),
          nrow = 3, align = "v", heights = c(0.8,0.05,1))

ggsave("figures/trapping.validation.naive.combined.jpg", 
       width = 9, height = 6, units = "in", dpi = 300)
