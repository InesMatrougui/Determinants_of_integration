library(dplyr)
library(ggplot2)
library(tidyverse)

# ============================================================================
# DEPTH BUSCO ANALYSIS
# ============================================================================
# Calculate sequencing depth on host genome compared to BUSCO genes

# Host BUSCO coordinates were identified using BUSCO v. 5.7.1 with the insect_odb10 database
# Command unix: busco -i ${host}.fasta -c 20 -l insect_odb10 -m geno -o BUSCO_${host}.fasta -f
# Mean depth for each BUSCO gene was calculated using CovWindows software https://github.com/HeloiseMuller/CovWindows
# Command unix: CovWindows/build/CovWindows -w BUSCO_sharded_${samples}.txt -f -c ${samples}_trimmed_vs_${host}_coverage_positions -d -o ${samples}_BUSCO_coverage_positions.txt
# Coverage positions file generated from WorkflowBowBlast (bowtie2=TRUE, coverage=TRUE)

setwd("dir")

# Define sample names
samples = c("BfCsIn1","BfCsIn2", "BfCsIn3", "BfCsCo1", "BfCsCo2", "BfCsCo3", "In1", "In2", "In3", "Co1","Co2", "Co3", "a1",
            "FAW4", "FAW5", "FAW6", "Cp4", "Cp5", "Cp6_b", "SfCi1", "SfCi2")

# Function to calculate mean log-transformed depth for host BUSCO genes
mean_log <- function(ind) {
  ind <- read.table(paste0(ind, "_BUSCO_coverage_positions.txt"), header=FALSE, sep=' ', quote = "") %>% 
    select(-V2, -V4, -V6) %>% 
    # Replace zero values with small pseudocount to avoid log(0)
    mutate(V8=if_else(V8==0, 0.00001, V8)) %>% 
    mutate(log_mean=log(V8)) %>% 
    summarise(log.host_BUSCO_depth=mean(log_mean),
              host_BUSCO_depth=mean(V8)) %>% 
    mutate(sample=ind)
  return(ind)
}

# Apply function to all samples and combine results
BUSCO_host <- bind_rows(lapply(samples, mean_log))

# ============================================================================
# WASP BUSCO DEPTH CALCULATION
# ============================================================================
# Calculate sequencing depth on wasp genome compared to BUSCO genes

# Wasp BUSCO coordinates were identified using BUSCO v. 5.7.1 with the insect_odb10 database
# Command unix: busco -i ${wasp}.fasta -c 20 -l insect_odb10 -m geno -o BUSCO_${wasp}.fasta -f
# Mean depth for each BUSCO gene was calculated using CovWindows software https://github.com/HeloiseMuller/CovWindows
# Command unix: CovWindows/build/CovWindows -w BUSCO_sharded_${wasp}.txt -f -c ${samples}_trimmed_vs_${wasp}_coverage_positions -d -o ${samples}_BUSCO_coverage_positions.txt


# Function to calculate mean depth for wasp BUSCO genes
mean <- function(ind) {
  dt <- read.table(paste0(ind, "_BUSCO_wasp_coverage_positions.txt"), header=FALSE, sep=' ', quote = "") %>% 
    select(-V2, -V4, -V6)
  mean_BUSCO <- data.frame(wasp_BUSCO_depth=sum(dt$V8)/length(dt$V8), 
                           sample=ind)
  return(mean_BUSCO)
}

# Apply function to all samples and combine results
BUSCO_wasp <- bind_rows(lapply(samples, mean))

# Merge host and wasp depth data
host_wasp_depth <- left_join(BUSCO_host, BUSCO_wasp, by ='sample')

# ============================================================================
# DATA INTEGRATION
# ============================================================================
# Merge all data: integration data and depth data for Bayesian analyses

# Define segments containing Heritable Integrated Mutations (HIM)
SEG_HIM <- data.frame(Segment=c('Segment_1', 'Segment_4', 'Segment_7', 
                                'Segment_10', 'Segment_11', 'Segment_12', 
                                'Segment_14', 'Segment_16', 'Segment_17', 
                                'Segment_18', 'Segment_24', 'Segment_26', 
                                'Segment_27', 'Segment_28', 'Segment_32', 
                                'Segment_35'), HIM='y')

# Define segments and their corresponding replication units
SEG_RU <- data.frame(Segment=c('Segment_1', 'Segment_4', 'Segment_7', 
                                'Segment_10', 'Segment_17','Segment_11', 'Segment_12', 
                                'Segment_14', 'Segment_26', 
                                'Segment_27', 'Segment_28', 'Segment_32', 
                                'Segment_35', 'Segment_18', 'Segment_24','Segment_16', 'Segment_15', 
                                'Segment_20/33', 'Segment_2', 'Segment_36', 'Segment_13', 'Segment_9', 'Segment_31', 'Segment_22',
                                'Segment_19', 'Segment_25', 'Segment_30', 'Segment_23', 'Segment_6', 'Segment_5', 
                                'Segment_21', 'Segment_8', 'Segment_37'),
                      Replication_Unit=c('RU5', 'RU7', 'RU4', rep('RU3', 2), rep('RU6.1', 2), 'RU6.2', 'RU8', 
                                         rep('RU2.3',8), rep('RU2.1', 7), rep('RU1', 6), rep('RU9',2), 'unknown'))

# Read integration and depth data
ChimD <- read.table("chimera_depth.txt", header = TRUE, sep = "\t")

# Add metadata to dataset
ChimD <- mutate_all(ChimD, ~replace(., is.na(.), 0)) %>% 
  # Standardize segment names
  mutate(Segment= if_else(Segment=="Segment_33" | Segment=="Segment_20", "Segment_20/33", Segment)) %>% 
  # Add HIM presence/absence information
  merge(SEG_HIM, by = "Segment", all = TRUE) %>% 
  # Add host susceptibility information (specific wasp-host pairs known to be susceptible)
  mutate(host_susceptibility= if_else((wasp=="cotesia_sesamiae_mombasa"& host=="busseola_fusca")|
                                        (wasp=="cotesia_flavipes"& host=="spodoptera_frugiperda"), "n", "y")) %>% 
  # Replace remaining NA values with 'n'
  mutate_all( ~replace(., is.na(.), 'n')) %>%
  # Add replication unit information
  merge(SEG_RU, by = "Segment", all = TRUE) 

# Convert numeric columns to numeric type
ChimD$int  <- as.numeric(ChimD$int)
ChimD$depth  <- as.numeric(ChimD$depth)

# Merge integration data with BUSCO depth data for Bayesian analysis
ChimD_Bayes <- left_join(ChimD, host_wasp_depth , by = c("sample")) %>%
    select(-host_BUSCO_depth)
    
# Export dataset for Bayesian analysis
write.table(ChimD_Bayes, "chimera_depth_Bayes.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

# ============================================================================
# DATA STANDARDIZATION FOR VISUALIZATION
# ============================================================================

ChimD_norm <- left_join(ChimD, host_wasp_depth , by = c("sample")) %>% 
  select(-log.host_BUSCO_depth) %>% 
  # Normalize depth by dividing by host BUSCO depth (sequencing effort) and subtracting wasp BUSCO depth (contamination)
  mutate(depth=((depth*1)/host_BUSCO_depth)-((wasp_BUSCO_depth*1)/host_BUSCO_depth)) %>% 
  # Ensure depth values are non-negative
  mutate(depth=if_else(depth<0, 0, depth)) %>% 
  # Normalize integration counts by host BUSCO depth (sequencing effort)
  mutate(int=(int*1)/host_BUSCO_depth)%>% 
  # Keep only segments with HIM
  filter(HIM=='y') %>% 
  # Create combined wasp-host category column
  unite(wasp_h, wasp, host, sep = " : ", remove = FALSE)

# Reshape data to long format for visualization
Chim_point <- ChimD_norm %>% 
  pivot_longer(cols=c(int, depth), names_to = "Mean", values_to = "values")%>% 
  ungroup() %>% 
  # Replace non-numeric values with 0
  mutate(values=if_else(values=="NaN"|values=="Inf",0,values)) %>%
  select(-HIM, -host_susceptibility, -wasp_BUSCO_depth, -host_BUSCO_depth)

Chim_Int <- Chim_point

# Export formatted data
write.table(Chim_Int, "Chim_Int.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
