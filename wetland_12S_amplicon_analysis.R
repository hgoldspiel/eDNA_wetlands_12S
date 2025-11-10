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
  guides(fill=guide_legend(ncol=2))

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

# plot diversity by wetland type
Diversity_melt_all <- melt(Diversity, 
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
Diversity_melt_nofish <- melt(Diversity, 
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
         Hatch = ifelse(Date < "2023-06-01", "No", "Yes"),
         Wetland_vol = 0.5*Area_sq_m*Spring_depth_cm) %>%
  filter(Species == "Ambystoma maculatum" & !is.na(AMMA_eggs)) %>%
  group_by(Site, DOY, eDNA_day, Hatch) %>% 
  summarize(Reads = mean(Abundance),
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
ggplot(amma_data, aes(AMMA_eggs, Reads)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() + mythemes

# relationship, corrected for vol
ggplot(amma_data, aes(AMMA_eggs_per_m2, Reads_per_ml_m3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() + mythemes

# relationship, corrected for vol, log transformed data (check that its kosher)
ggplot(amma_data, aes(log(AMMA_eggs_per_m2+1), log(Reads_per_ml_m3+1))) +
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
ggplot(amma_data, aes(log(AMMA_eggs_per_m2+1), log(Reads_per_ml_m3+1), color = Hatch)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() + mythemes

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
                                                 length.out = 100),
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

# parametric model to assess significance with lmer
library(lme4)

# inspect data distribution
hist(amma_data$Reads_per_ml_m3, prob = T); 
lines(density(amma_data$Reads_per_ml_m3), col = "blue", lwd = 2)

hist(log(amma_data$Reads_per_ml_m3+1), prob = T); 
lines(density(log(amma_data$Reads_per_ml_m3 + 1)), col = "blue", lwd = 2)

# seems like a gamma or lognormal distribution, could also log Reads_per_ml_m3 or sqrt transform reads
# this paper used log transformation on reads for validating amphibian eDNA surveys 
# https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/1365-2656.13468
lmer.full <- lmer(log(Reads_per_ml_m3+1) ~ 
                    scale(log(AMMA_eggs_per_m2+1)) * scale(DOY) + (1|Site), 
                  data = amma_data)
lmer.add <- lmer(log(Reads_per_ml_m3+1) ~ 
                   scale(log(AMMA_eggs_per_m2+1)) + scale(DOY) + (1|Site), 
                 data = amma_data)
lmer.eggs <- lmer(log(Reads_per_ml_m3+1) ~ 
                    scale(log(AMMA_eggs_per_m2+1)) + (1|Site), 
                  data = amma_data)
lmer.DOY <- lmer(log(Reads_per_ml_m3+1) ~ 
                    scale(DOY) + (1|Site), 
                  data = amma_data)
lmer.null <- lmer(log(Reads_per_ml_m3+1) ~ 1 + (1|Site), 
                  data = amma_data)

library(performance)
check_normality(lmer.full) # no good!
check_heteroscedasticity(lmer.full) # good!
AICtab(lmer.full, lmer.add, lmer.eggs, lmer.DOY, lmer.null) 
anova(lmer.add, lmer.full) # additive model is better
summary(lmer.add)

## let's create prediction surfaces 
egg.date.pred <- as.data.frame(predict(lmer.add, re.form= NA,
                                       newdata = pred.grads, se.fit = TRUE))

egg.date.pred <- cbind(pred.grads, egg.date.pred)

library(scico)
ggplot(egg.date.pred, aes(x = AMMA_eggs_per_m2, y = DOY, z = fit)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  geom_contour(colour = "white") + 
  geom_point(data = amma_data, inherit.aes = FALSE, 
             aes(x = AMMA_eggs_per_m2, y = DOY), shape = "+", 
             color = "grey", size = 4) +
  scale_fill_scico(palette = 'davos') +
  theme_bw() + mythemes +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar(title = "log(Reads/ml/m3)"))

# Extract start value estimates from the null lmer
fixef_start <- fixef(lmer.null)
theta_start <- as.numeric(getME(lmer.null, "theta"))

# glmer with lognormal distribution (can't find starting values!)
# running with Gamma function causes RStudio to crash...
glmer.full <- glmer(Reads_per_ml_m3 ~ 
                     scale(AMMA_eggs_per_m2) * scale(DOY) + (1|Site), 
                   family = gaussian(link = "log"),
                   start = list(fixef = fixef_start, theta = theta_start),
                   data = amma_data)
glmer.add <- glmer(Reads_per_ml_m3 ~ 
                      scale(AMMA_eggs_per_m2) + scale(DOY) + (1|Site), 
                    family = gaussian(link = "log"),
                    start = list(fixef = fixef_start, theta = theta_start),
                    data = amma_data)
glmer.eggs <- glmer(Reads_per_ml_m3 ~ 
                     scale(AMMA_eggs_per_m2) + (1|Site), 
                   family = gaussian(link = "log"),
                   start = list(fixef = fixef_start, theta = theta_start),
                   data = amma_data)
glmer.DOY <- glmer(Reads_per_ml_m3 ~ 
                    scale(DOY) + (1|Site), 
                  family = gaussian(link = "log"),
                  start = list(fixef = fixef_start, theta = theta_start),
                  data = amma_data)
glmer.null <- glmer(Reads_per_ml_m3 ~ 1 + (1|Site), 
                   family = gaussian(link = "log"),
                   start = list(fixef = fixef_start, theta = theta_start),
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)),
                   data = amma_data)

AICtab(glmer.full, glmer.add, glmer.eggs, glmer.DOY, glmer.null)


### 7a-2. eDNA vs trapping ----
#### load trapping data ----
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

#### format eDNA data in same way for occupancy analysis ----
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

# combine trapping and eDNA data for analysis----
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
ggplot(pred.spp, aes(Date.mdy, Predicted, color = Method, fill = Method, lty = Method)) +
  geom_line(lwd = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  facet_wrap(~Species) +
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
  # make plots
  if(variable == "Method") {
    occu.plots.univar[[variable]][["Data"]]$Method <- 
      factor(occu.plots.univar[[variable]][["Data"]]$covariateValue, 
             levels = c("Trapping", "eDNA"))
    occu.plots.univar[[variable]][["Plot"]] <-
      ggplot(occu.plots.univar[[variable]][["Data"]], 
             aes(x = Method, y = Predicted, shape = Method, color = Method)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
      scale_color_manual(values = col2[1:2]) +
      facet_wrap(~Species) +
      labs(x = variable, y = "Probability of detection") +
      theme_bw() + mythemes
  } else {
    occu.plots.univar[[variable]][["Plot"]] <-
      ggplot(occu.plots.univar[[variable]][["Data"]], aes(x = covariateValue, y = Predicted)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      facet_wrap(~Species) +
      labs(x = variable, y = "Probability of detection") +
      theme_bw() + mythemes
  }
}

# print plots
occu.plots.univar$Method$Plot
occu.plots.univar$Date$Plot

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

ggplot(curve_data, aes(x = x, y = y, color = Method)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3) +
  geom_point(size = 2) +
  geom_line() +
  scale_color_manual(values = col2[1:2]) +
  scale_x_continuous(breaks = c(1:10)) +
  facet_wrap(~Species) +
  labs(x = "Sampling units (eDNA samples, trapping occasions)", y = "Cumulative detection prob.") +
  theme_bw() + mythemes

### 7a-3. eDNA vs cameras ----

# import camera data
cameras <- read.csv("data/camera_data_2025-10-27.csv")

# add latin names
library(taxize)
library(taxizedb)
# reformat and correct some species
cameras$Species.Common.Name[cameras$Species.Common.Name == " Beaver"] <- "Beaver"
cameras$Species.Common.Name[cameras$Species.Common.Name == "Opossum"] <- "Virginia opossum"
camera.spp <- unique(cameras$Species.Common.Name)
camera.spp <- camera.spp[camera.spp %notin% c("", "Unknown", "Unknown bird", "Unknown mammal", 
                                              "Unknown amphibian", "Unknown reptile", "Other (comment)")]
camera.scinames <- comm2sci(camera.spp)
camera.tax.names <- data.frame(Common.Name = NULL, Sci.Name = NULL)
for(i in names(camera.scinames)) {
  temp <- data.frame(Com.Name = i,
                     Sci.Name = ifelse(length(camera.scinames[[i]]) == 0, NA, 
                                       camera.scinames[[i]]))
  camera.tax.names <- bind_rows(camera.tax.names, temp)
}

# fill in some gaps for species
camera.tax.names <- camera.tax.names %>%
  mutate(Sci.Name = case_when(
    Com.Name == "Spotted sandpiper" ~ "Actitus macularius",
    Com.Name == "Bicknell's thrush" ~ "Catharus bicknelli" ,
    Com.Name == "Barred owl" ~ "Strix varia",
    Com.Name == "Grey fox" ~ "Urocyon cineroargenteus",
    Com.Name == "Porcupine" ~ "Erethizon dorsatum",
    Com.Name == "Eastern coyote" ~ "Canis latrans",
    Com.Name == "Snowshoe hare" ~ "Lepus americanus",
    Com.Name == "Deer mouse" ~ "Peromyscus sp.",
    Com.Name == "Woodthrush" ~ "Hylocichla mustelina",
    Com.Name == "Brown thrasher" ~ "Toxostoma rufum",
    Com.Name == "Vesper's sparrow" ~ "Pooecetes gramineus",
    Com.Name == "Winter wren" ~ 
    TRUE ~ as.character(Sci.Name)))
# add other taxonomy fields
# join everything to full camera dataset

tax <- name2taxid(spp, out_type = "summary")
blt$Kingdom <- NA
blt$Phylum <- NA
blt$Class <- NA
blt$Order <- NA
blt$Family <- NA
blt$Genus <- NA
blt$sscinames <- spp
# for each unique blasted species, populate extra taxonomic classes from taxize
for(i in 1:nrow(tax)) {
  classes <- classification(tax$id[i])
  blt[blt$sscinames == tax$name[i],"Kingdom"] <-
    ifelse(length(classes[[1]]$name[classes[[1]]$rank == "kingdom"]) > 0,
           classes[[1]]$name[classes[[1]]$rank == "kingdom"], NA)
  blt[blt$sscinames == tax$name[i],"Phylum"] <-
    ifelse(length(classes[[1]]$name[classes[[1]]$rank == "phylum"]) > 0,
           classes[[1]]$name[classes[[1]]$rank == "phylum"], NA)
  blt[blt$sscinames == tax$name[i],"Class"] <-
    ifelse(length(classes[[1]]$name[classes[[1]]$rank == "class"]) > 0,
           classes[[1]]$name[classes[[1]]$rank == "class"], NA)
  blt[blt$sscinames == tax$name[i],"Order"] <-
    ifelse(length(classes[[1]]$name[classes[[1]]$rank == "order"]) > 0,
           classes[[1]]$name[classes[[1]]$rank == "order"], NA)
  blt[blt$sscinames == tax$name[i],"Family"] <-
    ifelse(length(classes[[1]]$name[classes[[1]]$rank == "family"]) > 0,
           classes[[1]]$name[classes[[1]]$rank == "family"], NA)
  blt[blt$sscinames == tax$name[i],"Genus"] <-
    ifelse(length(classes[[1]]$name[classes[[1]]$rank == "genus"]) > 0,
           classes[[1]]$name[classes[[1]]$rank == "genus"], NA)
}


# Cluster analysis to compare communities detected from either method
# NMDS and PERMANOVA to compare community data structure
# Mantel tests to compare monthly presence/absence data at the camera sites
# Regression relating reads vs camera detections
# Paired bar plots comparing species reads vs number of detections (something quantitative and z-scored)
# Compare simple summary of richness from either method


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

