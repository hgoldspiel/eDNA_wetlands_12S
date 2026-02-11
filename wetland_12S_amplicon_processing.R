################################################################################
####~~~~~~~~~~~~~~~~~~ Maine wetland 12S metabarcoding 2023 ~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~ Code for processing Illumina data ~~~~~~~~~~~~~~~~~~####
################################################################################

# 0. RUN CUSTOM FUNCTIONS SCRIPT ----------------------
source("custom_functions_settings.R")

# 1. LOAD PACKAGES INTO R -----------------------------

# Define CRAN and Bioconductor packages
cran_packages <- c(
  "readxl", "lubridate", "reshape2", "rBLAST", "archive",
  "readr", "vegan", "ape", "lme4", "lmerTest", "emmeans",
  "permute", "randomForest", "ggpubr", "corrplot", "stringr",
  "tidyr", "tibble"
)

bioc_packages <- c(
  "dada2", "GenomeInfoDbData", "ShortRead", "phyloseq",
  "phylotools", "biomartr", "DESeq2", "DECIPHER"
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

# 2. LOAD RAW DATA ---------------------------
## 1.  Specify folder path to all raw (unfiltered) sequence files
path_raw <- paste0(getwd(), "/data/raw_seqs/")

## note that we already ran our cutadapt script in Git Bash to splice out the primers and Nextera tails

## 2.  Specify file paths for forward (R1) and reverse (R2) reads
fnsF <- list.files(path_raw, recursive = TRUE, pattern = "_R1_001_trimmed")
fastqsF <- fnsF[grepl('.gz$', fnsF)] 

fnsR <- list.files(path_raw, recursive = TRUE, pattern = "_R2_001_trimmed")
fastqsR <- fnsR[grepl('.gz$', fnsR)] 

## 3.  Fully specify the path, and file names (fns) 
fnsF <- file.path(path_raw, fastqsF)
fnsR <- file.path(path_raw, fastqsR)

## 4. Pull the sample names (the run names, not the eDNA sample ID)
sample.namesF <- sapply(strsplit(basename(fastqsF), '_R1'), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(fastqsR), '_R2'), `[`, 1)

## 5. Pull the shortened sample names
sample.names.shortF <- substr(sample.namesF, 1, 5)
sample.names.shortR <- substr(sample.namesR, 1, 5)

## 6.  Assign the sample names to the list of fastq files
names(fnsF) <- sample.names.shortF 
names(fnsR) <- sample.names.shortR

## 7. Load metadata and group the sequence files based on primer set
metadata.runs <- read_excel("metadata/Metadata for Bioinformatics.xlsx", 1)
metadata.samples <- read_excel("metadata/Metadata for Bioinformatics.xlsx", 2,
                               col_types = c(rep("text", 6), "date", "date", rep("numeric", 3),
                                             "text", "text", "date", "text", "text", "numeric",
                                             "text", "text", "text", rep("numeric", 14)))
metadata.sites <- read_excel("metadata/Site_data_2023.xlsx",
                             col_types = c(rep("text",3), rep("numeric", 4), "text", "numeric",
                                           rep("text", 4), "date", rep("numeric", 4), 
                                           rep("text", 4)))
# make wetland type variable two classes (ephemeral, semi/permanent)
metadata.sites$Wetland_type2 <- metadata.sites$Wetland_type
metadata.sites$Wetland_type2[metadata.sites$Wetland_type2 %in% 
                               c("Semipermanent", "Permanent")] <- "Semi/permanent"

# impute missing pH value for one sample
metadata.samples$pH[metadata.samples$Tube_ID == "UFTS27-S-1"] <- 
  mean(na.omit(metadata.samples$pH[
    metadata.samples$Site == "UFTS27"]))

# classify samples and runs
metadata.samples$sample_type <- "Sample"
metadata.samples$sample_type[metadata.samples$Location == "B"] <- "Field blank"
metadata.samples$sample_type[metadata.samples$Tube_ID == "Qiagen-Blank-Pool"] <- "Qiagen blank"
metadata.samples$sample_type[metadata.samples$Tube_ID %in% c("NTC-1-NTC-3", "NTC2")] <- "PCR blank"
metadata.samples$sample_type[metadata.samples$Tube_ID %in% c("Riaz-Mock-1", "Riaz-Mock-2", "Riaz-Mock-3")] <- "Mock"

# join up the modified sample metadata with the run metadata
metadata.runs <- left_join(metadata.runs, metadata.samples, by = "Tube_ID")
# join site-level data to samples
metadata.runs <- left_join(metadata.runs, metadata.sites, by = c("Site" = "SiteID"))

# 3. CHECK SEQUENCE QUALITY ---------------------------
## 1. Review a few random samples for sequencing quality.   
# Quality stores >30 are good, < 20 are poor and should be trimmed
# The green line is the average quality score for all sequences in that sample, 
# the black/grey are relative abundance of that score at each base.

plotQualityProfile(fnsF[1]) # pretty great until 110
plotQualityProfile(fnsR[1]) # about the same, great until 110

## 2.  Look at the whole dataset to aggregate all samples onto one graph. 
### forward reads
qual.plot.f <- plotQualityProfile(fnsF, aggregate=T) +
  geom_hline(aes(yintercept = 30), lty = "dashed") +
  geom_hline(aes(yintercept = 20), lty = "dashed") +
  ggtitle("Forward")
reads.sum.f <- max(qual.plot.f$data$Cycle)

### reverse reads
qual.plot.r <- plotQualityProfile(fnsR, aggregate=T) +
  geom_hline(aes(yintercept = 30), lty = "dashed") +
  geom_hline(aes(yintercept = 20), lty = "dashed") +
  ggtitle("Reverse")
reads.sum.r <- max(qual.plot.r$data$Cycle)

# save figure
ggarrange(qual.plot.f, qual.plot.r,
          nrow = 2, ncol = 1,
          labels = c("A", "B"))
ggsave("figures/read.quality.jpg", width = 8, height = 12, units = "in", dpi = 300) 

### Summary of quality 

# This dataset contains 196 samples, with 8,513,187 total reads

# The R1 (Forward) avg quality goes below 30 at cycle (the nucleotide position) 150
# The R2 (Reverse) avg quality goes below 30 at cycle 150

# 4. FILTER/TRIM RAW SEQUENCES ---------------------------
## 1. Create folders on your computer to place your filtered data into.
dir.create("data/filtered") 
dir.create("data/filtered_forward") 
dir.create("data/filtered_reverse") 

## 2. Specify the folder path to place filtered files for reads 1 and 2.
path_filt_F <- "data/filtered_forward"
path_filt_R <- "data/filtered_reverse"

### Filter and trim Forward/Reverse
# specify where to send the filtered files, and what to name them
filtFs <- file.path(path_filt_F, paste0(sample.namesF, "_F_filt.fastq.gz"))
filtRs <- file.path(path_filt_R, paste0(sample.namesR, "_R_filt.fastq.gz"))

# associate the same names with the filtered files
names(filtFs) <- sample.names.shortF
names(filtRs) <- sample.names.shortR

# filter
filtoutput <- filterAndTrim( 
    fwd = fnsF, filt = filtFs, 
    rev = fnsR, filt.rev = filtRs, 
    trimRight = c(10,10),
    minLen = 50,
    maxEE = c(2,2), ## The max number of expected errors allowed in a read
    maxN = 0, # max number of ambiguous bases (N) allowed, DADA2 doesn't allow any
    rm.phix = TRUE, # remove any PhiX DNA (used as positive control)
    verbose = 0)


# save the filtered variable for sample totals
saveRDS(filtoutput, "data/filtered/filtered_output_wetlands_2023.rds") 
filtoutput <- readRDS("data/filtered/filtered_output_wetlands_2023.rds")

# take a look at the counts
print(filtoutput)

# order the info by the first column (reads.in)
filtoutput[order(filtoutput[,1], decreasing = FALSE),]
# order the info by the second column (reads.out) which are filtered reads
filtoutput[order(filtoutput[,2], decreasing = FALSE),]

# sum the reads
sum(filtoutput[,1]) # raw
sum(filtoutput[,2]) # filtered

# Summary of filtering
## Which parameters did you use and why? 
## Trim left = 0, Trim right = 10 (approx. threshold of qualscore < 30)

## What was the largest and smallest number of raw reads? 
## Largest was HG176 (a blank, UFTS27-B-2) with 329531
## Smallest was HG186 (a blank, FB-C-1) with 45

## What was the largest and smallest number of filtered reads? largest = 524,038, smallest = 0
## Largest was HG176 (a blank, UFTS27-B-2) with 312503
## Smallest was HG083 (a blank, FS-S-B-4) with 10

## How many total raw reads?
## 8513187

## How many total filtered reads?
## 8003563

# 5. DADA2 LEARN ERROR RATES ----------------------------
# Best to learn error rates per unique run and prior to dereplication

set.seed(1234) 

# Paths to filtered files
path_filt_F <- "data/filtered_forward"
path_filt_R <- "data/filtered_reverse"

# Get names for filtered files
## forward
fnsfiltsF <- list.files(path_filt_F, full.names = TRUE)
filtsF <- fnsfiltsF[grepl('.gz$', fnsfiltsF)]
sample.namesF <- sapply(strsplit(basename(filtsF), '_filt'), `[`, 1) 
names(filtsF) <- sample.namesF 
sample.names.shortF <- substr(sample.namesF, 1, 5)

# reverse
fnsfiltsR <- list.files(path_filt_R, full.names = TRUE)
filtsR <- fnsfiltsR[grepl('.gz$', fnsfiltsR)] 
sample.namesR <- sapply(strsplit(basename(filtsR), '_'), `[`, 1) 
names(filtsR) <- sample.namesR
sample.names.shortR <- substr(sample.namesR, 1, 5)

# simpler sample names
names(filtsF) <- sample.names.shortF
names(filtsR) <- sample.names.shortR

# 4. loop to Learn error rates from each sequencing run
errF <- learnErrors(filtsF, nbases = 1e8, multithread=TRUE, randomize=TRUE, verbose = 0)
errR <- learnErrors(filtsR, nbases = 1e8, multithread=TRUE, randomize=TRUE, verbose = 0)

# 5. Save the error profiles to output files
saveRDS(errF, "data/results/error_profile_F.RDS") 
saveRDS(errR, 'data/results/error_profile_R.RDS') 

# 6. Optional, look at error plots to view the prediction rates of different sequencing errors.
# From the function info page: "This function plots the observed frequency of each transition (eg. A->C) as a function of the associated quality score. It also plots the final estimated error rates (if they exist). The initial input rates and the expected error rates under the nominal definition of quality scores can also be shown."
# Thread on plot interpretation: https://github.com/benjjneb/dada2/issues/975 

# forward reads
err.plots.f <- plotErrors(errF, nominalQ=TRUE) + ggtitle("Forward")
err.plots.f
# save plot
ggsave(filename = "figures/errF.jpg",  width = 8, height = 6, units = "in", dpi = 300)

# reverse reads
err.plots.r <- plotErrors(errR, nominalQ=TRUE) + ggtitle("Reverse")
err.plots.r
# save plot
ggsave(filename = "figures/errR.jpg", width = 8, height = 6, units = "in", dpi = 300)

# 6. DEREPLICATE AND INFER ASVs ----------------------
derep.fw <- derepFastq(filtsF, verbose=FALSE)
names(derep.fw) <- metadata.runs$ReadFile_Name
derep.rv <- derepFastq(filtsR, verbose=FALSE)
names(derep.rv) <- metadata.runs$ReadFile_Name

# Inferring ASVs from sequences using error profiles
dadaFs <- dada(derep.fw, err = errF, multithread = TRUE, verbose = 0, pool = TRUE)
dadaRs <- dada(derep.rv, err = errR, multithread = TRUE, verbose = 0, pool = TRUE)

# 2. Construct a sequence table. DADA2 may drop samples that didn't have any sequences which passed the quality control steps.
# merge pairs
mergers <- mergePairs(dadaFs, derep.fw, dadaRs, derep.rv, minOverlap = 20, verbose = 0)
# make sequence table
seqtab <- makeSequenceTable(mergers)

# 3. Get the dimensions of your table, which is a great way to count how many samples are and SVs included.
dim(seqtab)
View(seqtab)

## 196 samples and 6722 sequences

# 4. Save the seqtab before you quit!
saveRDS(seqtab, 'data/results/seqtab_wetlands.rds')


# 7. DADA2 REMOVE CHIMERAS ----------------------------	
# Chimeras are accidentally created sequences during lab protocols. Remove them.	

# 1. Remove chimeras. Leave the method as consensus. 
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose=TRUE) 

## 1720 bimeras

# 2. Check dimensions of your cleaned sequence tables.
dim(seqtab.nochim)

## 196 samples and 5002 sequences

# 3. Save your chimera-free sequence table.
saveRDS(seqtab.nochim, 'data/results/seqtab_nochim.rds')


# 8. DECONTAM / WORKFLOW VERIFICATION --------------------------

# Add blanks with no reads back in the seqtab.nochim as zero-ed rows
seqtab.nochim <- readRDS('data/results/seqtab_nochim.rds') # reload the seqtab

# 1. Track reads through the analysis here to see how many were lost at each QC step, in case you were too stringent. Only need to count F reads, even if you did paired F/R
getN <- function(x) sum(getUniques(x)) #set function to count the unique reads

# reformat files for joining metadata to ASV data
filtoutput.df <- as.data.frame(filtoutput)
filtoutput.df$sample <- substr(rownames(filtoutput.df),1,5)
filtoutput.df.nochim <- filtoutput.df[filtoutput.df$sample %in% substr(rownames(seqtab.nochim), 1,5),] 

# 2. Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor(s), all into a new variable
track <- left_join(filtoutput.df, metadata.runs, by = join_by("sample" == "ReadFile_Name"))
seqtab.ASVs.df <- as.data.frame(rowSums(seqtab > 0))
seqtab.reads.df <- as.data.frame(rowSums(seqtab))
seqtab.nochim.reads.df <- as.data.frame(rowSums(seqtab.nochim))
seqtab.nochim.ASVs.df <- as.data.frame(rowSums(seqtab.nochim > 0))
colnames(seqtab.ASVs.df) <- "ASVs.out"
colnames(seqtab.reads.df) <- "reads.dada"
colnames(seqtab.nochim.reads.df) <- "reads.nonchimeras"
colnames(seqtab.nochim.ASVs.df) <- "ASVs.nonchimeras"
seqtab.ASVs.df$sample <- rownames(seqtab.ASVs.df)
seqtab.reads.df$sample <- rownames(seqtab.reads.df)
seqtab.nochim.reads.df$sample <- rownames(seqtab.nochim.reads.df)
seqtab.nochim.ASVs.df$sample <- rownames(seqtab.nochim.ASVs.df)
track <- left_join(track, seqtab.ASVs.df, by = "sample")
track <- left_join(track, seqtab.reads.df, by = "sample")
track <- left_join(track, seqtab.nochim.reads.df, by = "sample")
track <- left_join(track, seqtab.nochim.ASVs.df, by = "sample")
track <- track[,c("Tube_ID", "sample", 
                  "reads.in", "reads.out", "reads.dada", "reads.nonchimeras", 
                  "ASVs.out", "ASVs.nonchimeras", "sample_type")]
track$sample_type <- factor(track$sample_type, 
                               levels = c("Sample", "Mock", "Field blank", "Qiagen blank", "PCR blank"))
# force NAs to 0
track[is.na(track)] <- 0

track.melt <- melt(track, id = c("Tube_ID", "sample", "sample_type")) %>%
  dplyr::rename(dada = variable)

# 3. Save the tracking variable and data as an R file
saveRDS(track, "data/results/tracked_seqs_pools_2021.RDS")

# 4. Plot the number of reads by sample type and quality control stage
track.melt %>%
  filter(dada %notin% c("ASVs.out", "ASVs.nonchimeras")) %>%
  ggplot(aes(x = sample_type, y = value, color = dada)) + 
  geom_point(position = position_dodge(0.4), alpha = 0.2) + 
  geom_boxplot(width = 0.3, position = position_dodge(0.4)) +
  ylab("Reads") + xlab("Sample type") +
  scale_color_manual(values = col1[1:4]) + 
  theme_bw() + theme(legend.position = "bottom", 
                     axis.text.x = element_text(angle = 25, hjust = 1)) +
  ggtitle("Total reads by sample type") + mythemes

ggsave("figures/reads.dada.nochim.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

# It looks like we have some serious contamination (mostly human DNA)

# 5. Remove contamination from blanks

# Add blanks with no reads back in the seqtab.nochim as zero-ed rows
seqtab.nochim <- readRDS('data/results/seqtab_nochim.rds') # reload the seqtab

# Loop through sampling stages to remove reads from paired blanks
mdat <- metadata.runs
dat <- seqtab.nochim # new standalone matrix
# PCR blank IDs
b_PCRs <- mdat$ReadFile_Name[mdat$sample_type == "PCR blank"]
# Extraction blank ID
b_Qiagen <- mdat$ReadFile_Name[mdat$sample_type == "Qiagen blank"]
# 1. Remove the PCR contamination
for(i in 1:nrow(dat)) {
  readID <- rownames(dat)[i]
  b_PCR <- mdat$PCR_blank[mdat$ReadFile_Name == readID]
  b_PCR <- mdat$ReadFile_Name[mdat$Tube_ID == b_PCR]
  ## 1a. If the sample is a PCR blank, keep original read counts
  if(readID %in% b_PCRs) {
    dat[i,] <- dat[i,]
  }
  ## 1c. Otherwise, subtract the number of reads from the matching PCR blank
  else {
    reads <- dat[i,]  
    bp_reads <- dat[which(rownames(dat) == b_PCR),]
    dat[i,] <- reads - bp_reads
  }
} 
dat[dat < 0] <- 0
# 2. Remove the Qiagen extraction contamination
for(i in 1:nrow(dat)) {
  readID <- rownames(dat)[i]
  ## 2a. If the sample is the Zymo or a PCR blank, keep original read counts
  if(b_Qiagen == readID | readID %in% b_PCRs) {
    dat[i,] <- dat[i,] 
  }
  ## 2b. Subtract the number of reads from the Zymo blank
  else {
    reads <- dat[i,] 
    bq_reads <- dat[which(rownames(dat) == b_Qiagen),]
    dat[i,] <- reads - bq_reads
  }
} 
# Force negative values to be 0
dat[dat < 0] <- 0
# 3. Remove the field contamination
for(i in 1:nrow(dat)) {
  readID <- rownames(dat)[i]
  b_field <- mdat$Field_blank[mdat$ReadFile_Name == readID]
  b_field <- 
    ifelse(!is.na(b_field), mdat$ReadFile_Name[mdat$Tube_ID == b_field],
           ifelse(is.na(b_field) & mdat$sample_type[mdat$ReadFile_Name == readID] == "Field blank",
                  readID, NA))
  ## 3a. If the sample has no field blank or is the field blank, keep reads
  if(is.na(b_field) | b_field == readID) {
    dat[i,] <- dat[i,] 
  } 
  ## 3b. Subtract the number of reads for that ASV found in the field blank
  else {
    reads <- dat[i,] 
    bf_reads <- dat[which(rownames(dat) == b_field),]
    dat[i,] <- reads - bf_reads
  }
} 
# Force negative values to be 0
dat[dat < 0] <- 0

# 4. Remove containation from "tag-jumping" based on our mock community reads
## identify mock community sequences and calculate threshold for tag jumping
col_max_value <- apply(dat, 2, max)  # Column-wise maximum values
mocks <- mdat$ReadFile_Name[mdat$sample_type == "Mock"]
mock.ASVs <- numeric()
for (col in 1:ncol(dat)) {
  # Check if maximum value in column is found at specific rows
  if (all(max(dat[mocks, col]) == col_max_value[col])) {
    mock.ASVs <- c(mock.ASVs, col)
  }
}

# find the maximum reads of mock comm spp in other samples
tagjump.reads <- max(dat[rownames(dat) %notin% mocks, mock.ASVs]) 

# populate the list of thresholds to refer to later, if needed
tagjump.thresholds <- tagjump.reads

# now subtract these reads from all reads in the matrix
dat <- dat - tagjump.reads

# Force negative values to be 0
dat[dat < 0] <- 0

# Convert numeric matrix to integer matrix
storage.mode(dat) <- "integer"
seqtab.decontam <- dat


# summary of ASVs and reads by assay
dim(seqtab.decontam)
sum(seqtab.decontam)

# Save your contamination-free sequence table.
saveRDS(seqtab.decontam, 'data/results/seqtab_decontam.rds')  

# Load the new decontaminated seq table, if needed
seqtab.decontam <- readRDS('data/results/seqtab_decontam.rds')

# 6. Append reads after decontam to tracking df and plot reads through dada workflow
seqtab.decontam.reads.df <- as.data.frame(rowSums(seqtab.decontam))
seqtab.decontam.ASVs.df <- as.data.frame(rowSums(seqtab.decontam > 0))
colnames(seqtab.decontam.reads.df) <- "reads.decontam"
colnames(seqtab.decontam.ASVs.df) <- "ASVs.decontam"
seqtab.decontam.reads.df$sample <- rownames(seqtab.decontam.reads.df)
seqtab.decontam.ASVs.df$sample <- rownames(seqtab.decontam.ASVs.df)
track <- left_join(track, seqtab.decontam.reads.df, by = "sample")
track <- left_join(track, seqtab.decontam.ASVs.df, by = "sample")
track$sample_type <- factor(track$sample_type, 
                               levels = c("Sample", "Mock", "Field blank", 
                                          "PCR blank", "Qiagen blank"))
track[is.na(track)] <- 0

# 7. Plot reads through workflow
track.melt <- melt(track, id = c("Tube_ID", "sample", "sample_type")) %>%
  dplyr::rename(dada = variable)

## Fig 1. Reads by sample type and decontamination stage
track.melt %>%
  filter(dada %notin% c("ASVs.out", "ASVs.nonchimeras", "ASVs.decontam") &
           sample_type %in% c("Sample", "Field blank", "Mock", "Qiagen blank")) %>%
  ggplot(aes(x = sample_type, y = value, color = dada)) + 
  geom_point(position = position_dodge(0.7), alpha = 0.2) + 
  geom_boxplot(width = 0.5, position = position_dodge(0.7)) +
  ylab("Reads") + xlab("Sample type") +
  scale_color_manual(values = col1[1:5]) + 
  theme_bw() + theme(legend.position = "bottom", 
                     axis.text.x = element_text(angle = 25, hjust = 1)) +
  ggtitle("Total reads by sample type") + mythemes

ggsave("figures/reads.sample.dada.decontam.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

## Fig 2. Total reads in samples by decontamination stage
track.melt %>%
  filter(sample_type == "Sample" &
           dada %notin% c("ASVs.out", "ASVs.nonchimeras", "ASVs.decontam")) %>%
  mutate(step = case_when(dada == "reads.in" ~ "raw",
                          dada == "reads.out" ~ "filtered",
                          dada == "reads.dada" ~ "dada",
                          dada == "reads.nonchimeras" ~ "nochim",
                          dada == "reads.decontam" ~ "decontam"),
         step = factor(step, levels = c("raw", "filtered", "dada", "nochim", "decontam"))) %>%
  ggplot(aes(step, value, fill = step)) +
  geom_violin() +
  geom_point(alpha = 0.1, position = position_jitter(0.05)) +
  scale_fill_manual(values = col1[1:5]) + 
  ylab("Reads") + xlab("dada stage") +
  theme_bw() + theme(legend.position = "bottom", 
                     axis.text.x = element_text(angle = 25, hjust = 1)) +
  ggtitle("Total reads over bioinformatic pipeline") + mythemes

ggsave("figures/reads.dada.decontam.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

## Fig 3. ASVs by sample type and decontamination stage
track.melt %>%
  filter(dada %in% c("ASVs.out", "ASVs.nonchimeras", "ASVs.decontam") &
           sample_type %in% c("Sample", "Field blank", "Mock", "Qiagen blank")) %>%
  ggplot(aes(x = sample_type, y = value, color = dada)) + 
  geom_point(position = position_dodge(0.7), alpha = 0.2) + 
  geom_boxplot(width = 0.5, position = position_dodge(0.7)) +
  ylab("ASVs") + xlab("Sample type") +
  scale_color_manual(values = col1[1:4]) + 
  theme_bw() + theme(legend.position = "bottom", 
                     axis.text.x = element_text(angle = 25, hjust = 1)) +
  ggtitle("Total ASVs by sample type") + mythemes

ggsave("figures/ASVs.sample.dada.decontam.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

## Fig 4. Total ASVs by decontamination stage
track.melt %>%
  filter(sample_type == "Sample" &
           dada %in% c("ASVs.out", "ASVs.nonchimeras", "ASVs.decontam")) %>%
  ggplot(aes(dada, value, fill = dada)) +
  geom_violin() +
  geom_point(alpha = 0.1, position = position_jitter(0.05)) +
  scale_fill_manual(values = col1[1:4]) + 
  ylab("ASVs") + xlab("dada stage") +
  theme_bw() + theme(legend.position = "bottom", 
                     axis.text.x = element_text(angle = 25, hjust = 1)) +
  ggtitle("Total ASVs over dada workflow") + mythemes

ggsave("figures/ASVs.dada.decontam.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

# 9. DADA2 ASSIGNING TAXONOMY FROM REF DB FILE ---------------------------------

# summary of sequences (total reads and number of sites present)
seqsums <- apply(seqtab.decontam, 2, sum) # sequence-wise read counts
seqsamp <- seqtab.decontam
seqsamp[seqsamp > 0] <- 1 # convert reads to binary presence/absence
seqsamp <- apply(seqsamp, 2, sum) # site-wise presences
seqsums <- data.frame(ASV = names(seqsums), 
                      total.reads = seqsums, 
                      total.samples = seqsamp,
                      row.names = NULL)

# ref libs
## Custom Riaz ref lib for the northeast (developed by Erin Grey, UMaine)
ref.libs.custom.edited <- "ref_libs/REFDB_12S_2025-02-20.fasta"
tax.levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# assign taxonomy
set.seed(123)
all.taxa.custom.edited <- 
  assignTaxonomy(seqs = seqtab.decontam, 
                 refFasta = ref.libs.custom.edited,
                 tryRC = TRUE, 
                 minBoot = 70,
                 taxLevels = tax.levels,
                 outputBootstraps = TRUE,
                 multithread = TRUE)

saveRDS(all.taxa.custom.edited, 'data/results/species_customrefs_edits_wetlands.rsd') 

# 10. BLAST ALL UNASSIGNED SEQUENCES ------------------------------
# create FASTAs for unassigned ASVs
ASV.list <- list()
# collate all ASVs and create unique sequence IDs
ASV.list$all <- data.frame(seq = rownames(all.taxa.custom.edited$tax),
                           ID = paste0("Sequence", seq(nrow(all.taxa.custom.edited$tax))))
# identify ASVs with missing assignments to Species
ASVs.unknown <- rownames(all.taxa.custom.edited$tax)[
      is.na(all.taxa.custom.edited$tax[,"Species"])]
# create FASTA for missing ASVs for each assay
ASV.list$unk <- 
  data.frame(seq.name = ASV.list$all$ID, 
             seq.text = ASV.list$all$seq) %>%
  filter(seq.text %in% ASVs.unknown)
# export FASTA
dat2fasta(ASV.list$unk, "data/results/ASV_sequences.fasta")

# download the entire NCBI nt database (2024-03-24)
# options(timeout = 300000)
# biomartr::download.database.all(db = "nt_euk", path = "data/nt_euk")

# unzip the nt references (2024-03-24)
ref.seqs <- list.files("data/nt_euk/") ## NOTE: THIS MUST BE MOVED TO THE LOCAL DATA FOLDER FOR ANY GIVEN ANALYSIS
#for(ii in ref.seqs) {
# archive_extract(paste0("data/nt_euk/",ii), dir = "data/nt_euk/")
#}

blast.f6 <-
  ## The fields you want back from blast
  c('qseqid', 'sseqid', 'scomnames', 'sscinames', 'sskingdoms', 'qcovs', 'pident')
# get path to BLAST software
blastn <- "C:/Program Files/NCBI/blast-2.15.0+/bin/blastn.exe"
# get path to local nt database
ntdb <- "nt_euk"
# get path to taxdb data
taxdb <- "C:/Program Files/NCBI/blast-2.15.0+/bin/taxdb"

# blast away!
local.blast <- function() {
  # change working directory to where the nt_euk sequences are
  cur <- getwd()
  on.exit(setwd(cur))
  setwd("data/nt_euk")
  system2(command = blastn,
          args = c(
            '-db', ntdb,
            '-num_threads', '10', # only works with local blasting
            '-outfmt', sprintf('"6 %s"', paste(collapse = ' ', blast.f6)),
            '-perc_identity','97',
            '-max_target_seqs', '6',
            '-query', "C:/Users/harrison.goldspiel/Desktop/eDNA_wetlands_12S/data/results/ASV_sequences.fasta",
            '-out', 'C:/Users/harrison.goldspiel/Desktop/eDNA_wetlands_12S/data/results/blast.txt'),
  wait = TRUE,
  stdout = TRUE  # capture the output in a character vector
  )
}

local.blast()

# read in the blast results
blast.list <- read.table(paste0(getwd(),"/data/results/blast.txt"),
               quote = "", sep = "\t", col.names = blast.f6)

# get additional taxonomy data from BLAST hits
library(taxize)
library(taxizedb)
blt <- blast.list
spp <- ifelse(grepl("\\;", blast.list$sscinames),
              substr(blast.list$sscinames, 1, 
                     regexpr(";", blast.list$sscinames, fixed=TRUE)-1),
              blast.list$sscinames)
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
# keep only one match (the highest match with most overlap) for each seq
blt <- blt %>%
  distinct(.) %>%
  group_by(qseqid) %>%
  arrange(desc(qcovs), desc(pident)) %>%
  dplyr::slice(1) %>%
  ungroup()

# append the new, reliable BLASTed taxonomy data to the existing custom list
# make new copy without the NCBI taxonomy codes
assign.tax.list <- all.taxa.custom.edited$tax
blt.valid <- blt # new copy to populate only with accurate and relevant matches
names(assign.tax.list) <- substring(names(assign.tax.list), 1, nchar(names(assign.tax.list))-4)
# remove the tax ID from the original list
assign.tax.list <- apply(assign.tax.list,
                         FUN = function(x) gsub('_.*','',x), 
                         MARGIN = 2)

# select thresholds of matches to keep (qcov > 0.97, pident > 0.99)
blt.valid <- blt.valid %>%
  filter(qcovs >= 97 & pident >= 99)

# connect the BLT to the ASV sequence ID
unknown.tax.list <- left_join(ASV.list$unk, 
                              blt.valid,
                              by = c("seq.name" = "qseqid"))
rownames(unknown.tax.list) <- unknown.tax.list$seq.text
unknown.tax.list <- unknown.tax.list %>%
  dplyr::select(Kingdom, 
                Phylum, 
                Class, 
                Order, 
                Family, 
                Genus, 
                Species = sscinames)

# add blasted taxonomy to the original assigned lists
final.tax.list <- merge(assign.tax.list,
                        unknown.tax.list,
                        by = "row.names",
                        all.x = TRUE)

row.names(final.tax.list) <- final.tax.list$Row.names
final.tax.list <- subset(final.tax.list, 
                         select = -c(Row.names))

# do some final cleaning to keep one set of taxonomic assignments
final.tax.list$Phylum <- ifelse(is.na(final.tax.list$Phylum.y),
                                final.tax.list$Phylum.x,
                                final.tax.list$Phylum.y)
final.tax.list$Class <- ifelse(is.na(final.tax.list$Class.y),
                               final.tax.list$Class.x,
                               final.tax.list$Class.y)
final.tax.list$Order <- ifelse(is.na(final.tax.list$Order.y),
                               final.tax.list$Order.x,
                               final.tax.list$Order.y)
final.tax.list$Family <- ifelse(is.na(final.tax.list$Family.y),
                                final.tax.list$Family.x,
                                final.tax.list$Family.y)
final.tax.list$Genus <- ifelse(is.na(final.tax.list$Genus.y),
                               final.tax.list$Genus.x,
                               final.tax.list$Genus.y)
final.tax.list$Species <- ifelse(is.na(final.tax.list$Species.y),
                                 final.tax.list$Species.x,
                                 final.tax.list$Species.y)
# save the final assignments
final.tax.list <- final.tax.list[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

# rename species that are misassigned to native taxa containing identical DNA barcodes
final.tax.list$Species[final.tax.list$Species %in% "Turdus cardis"] <- "Turdus migratorius"
final.tax.list$Species[final.tax.list$Species %in% "Cyanopica cyanus koreensis"] <- "Cyanocitta cristata"

saveRDS(final.tax.list, 'data/results/spp_assign_final_wetlands_2023.rsd') 

# Examine species data prior to and after decontamination
seqtab.nochim <- readRDS('data/results/seqtab_nochim.rds')
spptab.nochim <- seqtab.nochim
colnames(spptab.nochim) <- final.tax.list[,"Species"][match(colnames(spptab.nochim), 
                                                              rownames(final.tax.list))]
seqtab.decontam <- readRDS('data/results/seqtab_decontam.rds')
spptab.decontam <- seqtab.decontam
colnames(spptab.decontam) <- final.tax.list[,"Species"][match(colnames(spptab.decontam), 
                                                              rownames(final.tax.list))]

# 11. CREATE PHYLOSEQ OBJECTS ------------------------------
final.tax.list <- readRDS('data/results/spp_assign_final_wetlands_2023.rsd')
# create phyloseq objects
samples <- sample_data(as.data.frame(metadata.runs))
row.names(samples) <- metadata.runs$ReadFile_Name

ps.dat <- phyloseq(otu_table(as.matrix(seqtab.decontam), taxa_are_rows=FALSE), 
                   sample_data(samples),
                   tax_table(as.matrix(final.tax.list)))
# change taxa names
taxa_names(ps.dat) <- paste0("Seq", seq(ntaxa(ps.dat)))

## create a phyloseq object with only field samples
ps.dat.samps <- subset_samples(ps.dat, sample_type == "Sample")

## create a principle component for water quality variables
wq.vars <- c("Depth_cm", "Cond_uS_cm", "ODO_mg_L", "Sal_psu", "SpCond_uS_cm", "Temp..C", "TDS_mg_L", "pH")
pairs(data.frame(sample_data(ps.dat.samps)[,wq.vars.redced]))
wq.vars.reduced <- c("Depth_cm", "Cond_uS_cm", "ODO_mg_L", "Temp..C","pH")
wetland.wq.pc <- prcomp(sample_data(ps.dat.samps)[,wq.vars.reduced][
  complete.cases(sample_data(ps.dat.samps)[,wq.vars.reduced])], center = TRUE, scale. = TRUE)
summary(wetland.wq.pc) # first 2 PCs are 0.64 cumulative variance
wetland.wq.pc # first PC is depth (+), conductivity (+), and ODO (+), # second PC is temp (-) and pH (-)

# examine raw phyloseq objects
ps.dat
ps.dat.samps

saveRDS(ps.dat.samps, 'data/results/phyloseq_data_samples_wetlands_2023.rsd') 

# 12. SUMMARY STATS OF ASVs AND READS ---------
# get second dataset for just taxonomically-verified ASVs from each assay
final.tax.list.valid <- final.tax.list
final.tax.list.valid <- final.tax.list.valid[final.tax.list.valid$Phylum == "Chordata" &
                                    !is.na(final.tax.list.valid$Phylum),]

# Create table to populate proportions of ASVs assigned to higher and lower taxa
tax.assign.props <- data.frame(Kingdom = NA,
                               Phylum = NA,
                               Class = NA,
                               Order = NA,
                               Family = NA, 
                               Genus = NA,
                               Species = NA)

# make a copy of table to get proportions of validated taxa from each assay
tax.assign.props.valid <- tax.assign.props

# populate table through nested for loop
for(level in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  tax.assign.props[,level] = 
    sum(!is.na(final.tax.list[,level]))/nrow(final.tax.list)
  tax.assign.props.valid[,level] = 
    sum(!is.na(final.tax.list.valid[,level]))/nrow(final.tax.list)
}


# melt table for plotting summary of assignment success for ASVs
## all sequences
tax.assign.melt <- reshape2::melt(tax.assign.props, 
                                  variable.name = "Level", 
                                  value.name = "Prop_ASVs") 

## verified sequences
tax.assign.melt.valid <- reshape2::melt(tax.assign.props.valid, 
                                        variable.name = "Level", 
                                        value.name = "Prop_ASVs_valid")

tax.assign.melt$Prop_ASVs_valid <- tax.assign.melt.valid$Prop_ASVs_valid 
tax.assign.melt <- tax.assign.melt[-1,]

# plot taxonomic assignment success 
ASV.props <- ggplot(tax.assign.melt, aes(Level, Prop_ASVs)) +
  geom_col(fill = "grey") +
  geom_col(aes(Level, Prop_ASVs_valid), fill = "black", width = 0.25) +
  labs(x = "Taxonomic level", y = "Proportion of ASVs") +
  ylim(c(0,1)) +
  theme_bw() + mythemes +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4),
        legend.position = "none")

ASV.props

# save plot
ggsave("figures/ASVs.prop.assigned.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

# let's do the same as above, but for total reads instead of ASVs
## create output tables for proportions of reads
reads.assign.props <- data.frame(Kingdom = NA,
                                 Phylum = NA,
                                 Class = NA,
                                 Order = NA,
                                 Family = NA, 
                                 Genus = NA,
                                 Species = NA)

reads.assign.props.valid <- reads.assign.props

# populate proportion tables
ps.dat.melt <- psmelt(ps.dat)
ps.dat.melt.valid <- ps.dat.melt
# validation by assay-specific taxonomy
ps.dat.melt.valid <- ps.dat.melt.valid[ps.dat.melt.valid$Phylum == "Chordata",]
# calculate proportion of reads for each taxonomic level
for(level in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      reads.assign.props[,level] = 
        sum(ps.dat.melt$Abundance[!is.na(ps.dat.melt[,level])])/sum(ps.dat.melt$Abundance)
      reads.assign.props.valid[,level] = 
        sum(ps.dat.melt.valid$Abundance[!is.na(ps.dat.melt.valid[,level])])/sum(ps.dat.melt$Abundance)
}

# melt for plotting
reads.assign.melt <- reshape2::melt(reads.assign.props, 
                                    variable.name = "Level", 
                                    value.name = "Prop_reads")

reads.assign.melt.valid <- reshape2::melt(reads.assign.props.valid, 
                                          variable.name = "Level", 
                                          value.name = "Prop_reads_valid") 

reads.assign.melt$Prop_reads_valid <- reads.assign.melt.valid$Prop_reads_valid
reads.assign.melt <- reads.assign.melt[-1,]

# plot taxonomic assignment success (by reads) 
reads.props <- ggplot(reads.assign.melt, aes(Level, Prop_reads)) +
  geom_col(fill = "grey") +
  geom_col(aes(Level, Prop_reads_valid), fill = "black", width = 0.25) +
  labs(x = "Taxonomic level", y = "Proportion of reads") +
  ylim(c(0,1)) +
  theme_bw() + mythemes +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4),
        legend.position = "none")

reads.props

# save plot
ggsave("figures/Reads.prop.assigned.jpg", 
       width = 8, height = 6, units = "in", dpi = 300)

# combine ASV and reads plot
ggarrange(ASV.props + theme(axis.text.x = element_blank(), axis.title.x = element_blank()), 
          reads.props,
          nrow = 2, labels = c("A", "B"), heights = c(0.45,0.55))

# save plot
ggsave("figures/ASV.reads.prop.assigned.jpg", 
       width = 8, height = 10, units = "in", dpi = 300)

# how many sequences have no taxonomic information at all
final.tax.list %>%
    filter(if_all(everything(), is.na)) %>%
    nrow() / nrow(final.tax.list)

# how many reads have no taxonomic information at all
psmelt(ps.dat) %>%
      filter(if_all(c(Kingdom, Phylum, Class, Order, Family, Species), is.na)) %>%
      summarize(prop.unknown = sum(Abundance) / sum(psmelt(ps.dat)$Abundance)) 

# prop of assigned ASVs that were assigned to target taxa
tax.assign.correct <- as.matrix(tax.assign.props.valid) / as.matrix(tax.assign.props)
min(tax.assign.correct, na.rm = T)

# prop of assigned reads that were assigned to target taxa
reads.assign.correct <- as.matrix(reads.assign.props.valid) / as.matrix(reads.assign.props)
min(reads.assign.correct, na.rm = T)

# Contamination review
Riaz.blanks <- metadata.runs$ReadFile_Name[metadata.runs$sample_type %in% c("Field blank", 
                                                                            "Qiagen blank", 
                                                                            "PCR blank")]
View(spptab.nochim[Riaz.blanks,])

## prop of reads in raw data that were cow
cows.12s <- which(tax_table(ps.dat)[,"Genus"] %in% "Bos")
sum(seqtab.nochim[,cows.12s])/sum(seqtab.nochim)

## prop of reads in raw data that are human
human.12s <- which(tax_table(ps.dat)[,"Genus"] %in% "Homo")
sum(seqtab.nochim[,human.12s])/sum(seqtab.nochim)

# prop of reads in raw data that were "contaminants" and filtered out
(sum(seqtab.nochim) - sum(seqtab.decontam)) / (sum(seqtab.nochim))

# order cleaned samples by sequencing depth
sort(rowSums(otu_table(ps.dat.samps)))
## HG089 (FS-S-N-4) has over 200,000 cycles, >100K greater than all other samples




