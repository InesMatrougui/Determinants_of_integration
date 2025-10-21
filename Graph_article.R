# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
# Libraries for data manipulation and visualization
library(dplyr)
library(ggplot2)
library(gridExtra)      # For arranging multiple plots
library(ggpubr)         # For publication-ready figures
library(ggpirate)       # For pirate plots
library(ape)            # For phylogenetic tree manipulation

setwd("dir")

# ============================================================================
# DATA LOADING AND PREPARATION
# ============================================================================

# Load processed integration and depth data
Chim_Int <- read.table("Chim_Int.txt", header = TRUE, sep = "\t")

# ============================================================================
# FIGURE 2: INTEGRATION EFFICIENCY BY CIRCLE AND WASP-HOST PAIR
# ============================================================================
# Convert segment names to factor and rename to "Circle" format

# ============================================================================
# REORDER SEGMENT FACTORS FOR COMPREHENSIVE VISUALIZATION
# ============================================================================
# Convert segment names to factor and rename to "Circle" format

Chim_Int$Segment <- factor(Chim_Int$Segment, 
                           levels = c("Segment_1", "Segment_2","Segment_4", "Segment_5", "Segment_6", "Segment_7", "Segment_8", "Segment_9", "Segment_10", 
                                      "Segment_11", "Segment_12", "Segment_13", "Segment_14", "Segment_15", "Segment_16","Segment_17", "Segment_18", 
                                      "Segment_19","Segment_20/33", "Segment_21","Segment_22", "Segment_23","Segment_24", "Segment_25", "Segment_26", 
                                      "Segment_27", "Segment_28", "Segment_30","Segment_31","Segment_32", "Segment_35", "Segment_36", "Segment_37"), 
                           labels = c("Circle_1", "Circle_2","Circle_4", "Circle_5", "Circle_6", "Circle_7", "Circle_8", "Circle_9", "Circle_10", 
                                      "Circle_11", "Circle_12", "Circle_13", "Circle_14", "Circle_15", "Circle_16","Circle_17", "Circle_18", 
                                      "Circle_19","Circle_20/33", "Circle_21","Circle_22", "Circle_23","Circle_24", "Circle_25", "Circle_26", 
                                      "Circle_27", "Circle_28", "Circle_30", "Circle_31","Circle_32", "Circle_35", "Circle_36", "Circle_37"))

# Create pirate plot showing integration counts and depth values
# Stratified by wasp-host pairs and displayed separately for each circle with HIM
figure_bilan <- ggplot(Chim_Int) +
  aes(x=wasp_h, y=values, fill=wasp, color=Mean)  +
  # Pirate plot: combines points, lines (means), and bars with confidence intervals
  geom_pirate(points = TRUE, lines = TRUE, bars = TRUE, cis = FALSE, violins = FALSE, 
              points_params = list(shape = 16, size = 1.5), jitter_width = 0) +
  scale_color_manual(values = c("lightgray", "black")) +
  scale_fill_manual(values = c("#D95F02","#E6AB02", "#E7298A","#66A61E", "#7570B3")) +
  ylab("Number of integrations and depth(x)") +
  theme_bw()+
  coord_flip() +
  # Separate panels for each circle
  facet_grid(cols = vars(Segment)) +
  theme(strip.text.x = element_text(size = 7), axis.text.x = element_text(size = 6))+
  # Set specific order for wasp-host pairs on x-axis for consistency
  scale_x_discrete(limits=c("cotesia_icipe : spodoptera_frugiperda", "cotesia_flavipes : chilo_partellus", "cotesia_flavipes : spodoptera_frugiperda",
                            "cotesia_typhae : sesamia_nonagrioides", "cotesia_sesamiae_kitale : busseola_fusca", "cotesia_sesamiae_kitale : sesamia_calamistis",
                            "cotesia_sesamiae_mombasa : busseola_fusca", "cotesia_sesamiae_mombasa : sesamia_calamistis"))

# Export comprehensive integration efficiency figure
ggsave(filename = "figure_bilan.pdf", 
       plot = figure_bilan, 
       height = 9, 
       width = 15, 
       units = "in")

# ============================================================================
# PHYLOGENETIC TREE VISUALIZATION
# ============================================================================
# Read phylogenetic tree in Newick format with branch lengths and bootstrap values

myTree <- ape::read.tree(text = "(Cotesia_icipe:0.0424299990,(Cotesia_flavipes:0.0095996158,((Cotesia_sesamiae_kitale:0.0025423157,Cotesia_sesamiae_mombasa:0.0032904776)100/100:0.0038223243,Cotesia_typhae:0.0061466487)100/100:0.0219293140)100/100:0.0131139270);")

# Display and save phylogenetic tree with scale bar
pdf(file = "figure_bilan_phylo.pdf") 
add.scale.bar(plot(myTree))
dev.off()

# ============================================================================
# SAMPLE NAME MAPPING
# ============================================================================
# Create lookup table to map long sample identifiers to short names for plotting
# Format: "HostSpecies + WaspSpeciesAbbreviation"
# BfCsm = Busseola fusca + Cotesia sesamiae mombasa
# BfCsk = Busseola fusca + Cotesia sesamiae kitale
# ScCsm = Sesamia calamistis + Cotesia sesamiae mombasa
# ScCsk = Sesamia calamistis + Cotesia sesamiae kitale
# SfCf = Spodoptera frugiperda + Cotesia flavipes
# CpCf = Chilo partellus + Cotesia flavipes
# SfCi = Spodoptera frugiperda + Cotesia icipe
# SnCt = Sesamia nonagrioides + Cotesia typhae

names_samples <- data.frame(
  name = c("BfCsm1", "BfCsm2", "BfCsm3", "BfCsk1", "BfCsk2", "BfCsk3", "ScCsm1", "ScCsm2", "ScCsm3",
           "ScCsk1", "ScCsk2", "ScCsk3", "SfCf1", "SfCf2", "SfCf3", "CpCf1", "CpCf2", "CpCf3", "SfCi1", "SfCi2", "SnCt"), 
  sample = c("BfCsCo1", "BfCsCo2", "BfCsCo3", "BfCsIn1", "BfCsIn2", "BfCsIn3", "Co1", "Co2", "Co3", 
             "In1", "In2", "In3", "FAW4", "FAW5", "FAW6", "Cp4", "Cp5", "Cp6_b", "SfCi1", "SfCi2", "a1")
)

# Join short names to main dataset
Chim_Int <- left_join(Chim_Int, names_samples, by="sample")

# ============================================================================
# STATISTICAL SUMMARIES
# ============================================================================

# Summary 1: Total number of integrations per sample (HIM circles only)
nb_int_sample <- filter(Chim_Int, HIM=='y') %>% 
  group_by(name, wasp_h) %>% 
  summarize(nb_int = sum(int))

# Summary 2: Range and mean of integration counts per wasp-host pair
nb_int_mean_paires <- group_by(nb_int_sample, wasp_h) %>% 
  summarize(range_int = max(nb_int) - min(nb_int), 
            mean_nb_int = mean(nb_int))

# Summary 3: Mean depth stratified by HIM presence/absence
depth_HIM <- group_by(Chim_Int, HIM) %>% 
  summarize(mean_depth = mean(depth))

# Summary 4: Total sequencing depth per sample and wasp-host pair
depth_sample <- group_by(Chim_Int, name, wasp_h) %>% 
  summarize(sum_depth = sum(depth))

# Summary 5: Range and mean of summed depth per wasp-host pair
depth_mean_paires <- group_by(depth_sample, wasp_h) %>% 
  summarize(range_depth = max(sum_depth) - min(sum_depth), 
            mean_death = mean(sum_depth)) 

# Summary 6: Average number of integrations per circle and wasp-host pair (HIM circles only)
nb_int_segment <- filter(Chim_Int, HIM=='y') %>% 
  group_by(Segment, wasp_h) %>% 
  summarize(nb_int = mean(int))

# Summary 7: Average depth per circle and wasp-host pair
depth_segment <- group_by(Chim_Int, Segment, wasp_h) %>% 
  summarize(sum_depth = mean(depth))

# ============================================================================
# FIGURE S6A: DEPTH BY SAMPLE
# ============================================================================
# Bar plot showing total sequencing depth per sample, colored by wasp-host pair

p <- ggplot(left_join(depth_sample, depth_mean_paires, by=c("wasp_h")), aes(x=name, y=sum_depth, fill=wasp_h)) + 
  geom_col(color="black", linewidth=.5, alpha=.75) + 
  xlab("Samples") + ylab("depth(x) per hhg") + labs(title="A") + theme_bw() + 
  # Set sample order for consistent visualization
  scale_x_discrete(limits = c("BfCsm1", "BfCsm2", "BfCsm3",  "BfCsk1", "BfCsk2", "BfCsk3", "ScCsm1", "ScCsm2", "ScCsm3", 
                              "ScCsk1", "ScCsk2", "ScCsk3", "SfCf1", "SfCf2", "SfCf3", "CpCf1", "CpCf2", "CpCf3", "SnCt","SfCi1", "SfCi2")) +
  # Custom color palette for each wasp-host pair
  scale_fill_manual(values = c("#02818a","#80cdc1","gold1","#762a83","#c2a5cf","#a6dba0","#1b7837","#ef6548")) +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        legend.position = "none")

# ============================================================================
# FIGURE S3A: DEPTH DISTRIBUTION BY CIRCLE
# ============================================================================
# Box plot showing depth distribution across all circles, stratified by HIM presence

p1 <- ggplot(Chim_Int, aes(x=Segment, y=depth, color=HIM)) + 
  geom_boxplot() + 
  geom_point(alpha=0.5) +
  theme_bw() + xlab("Circles") + ylab("depth(x) per hhg") + labs(title="A") +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        legend.position = "none") +
  scale_color_manual(labels=c("circle without HIM", "circle with HIM"), values = c("#b2dfee", "#cd919e"))

# ============================================================================
# FIGURE S6B: INTEGRATION COUNTS BY SAMPLE
# ============================================================================
# Bar plot showing total integration counts per sample, colored by wasp-host pair

p2 <- ggplot(left_join(nb_int_sample, nb_int_mean_paires, by=c("wasp_h")), aes(x=name, y=nb_int, fill=wasp_h)) + 
  geom_col(color="black", linewidth=0.5, alpha=.75) + 
  xlab("Samples") + ylab("Number of integrations per hhg") + labs(title="B", fill = "host wasp pair") + theme_bw() + 
  # Set sample order for consistent visualization
  scale_x_discrete(limits = c("BfCsm1", "BfCsm2", "BfCsm3", "BfCsk1", "BfCsk2", "BfCsk3", "ScCsm1", "ScCsm2", "ScCsm3",
                              "ScCsk1", "ScCsk2", "ScCsk3", "SfCf1", "SfCf2", "SfCf3", "CpCf1", "CpCf2", "CpCf3", "SnCt","SfCi1", "SfCi2")) +
  # Custom colors and scientific names for legend
  scale_fill_manual(values = c("#02818a","#80cdc1","gold1","#762a83","#c2a5cf","#a6dba0","#1b7837","#ef6548"), 
                    labels=c(expression(italic("Cotesia flavipes - Chilo partellus")),
                             expression(italic("Cotesia flavipes - Spodoptera frugiperda")),
                             expression(italic("Cotesia icipe - Spodoptera frugiperda")),
                             expression(italic("Cotesia sesamiae kitale - Busseola fusca")),
                             expression(italic("Cotesia sesamiae kitale - Sesamia calamistis")),
                             expression(italic("Cotesia sesamiae mombasa - Busseola fusca")),
                             expression(italic("Cotesia sesamiae mombasa - Sesamia calamistis")),
                             expression(italic("Cotesia typhae - Sesamia nonagrioides"))))+
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

# ============================================================================
# FIGURE S3B: INTEGRATION COUNTS DISTRIBUTION BY CIRCLE
# ============================================================================
# Box plot showing integration distribution across circles, stratified by HIM presence
# Filter out zero values for better visualization

p3 <- ggplot(mutate(Chim_Int, int=if_else(int==0, NA, int)), aes(x=Segment, y=int, color=HIM)) + 
  geom_boxplot() + 
  geom_point(alpha=.5) +
  theme_bw() + xlab("Circles") + ylab("Number of integrations per hhg") + labs(title="B") +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_color_manual(labels=c("circle without HIM", "circle with HIM"), values = c("#b2dfee", "#cd919e"))

# ============================================================================
# SAVE COMBINED FIGURE: DEPTH AND INTEGRATIONS BY CIRCLE
# ============================================================================

ggsave("figure_depth_segment.pdf", 
       grid.arrange(p1, p3, ncol=2, nrow=1, widths = c(1.63, 2)),
       height = 6, 
       width = 15)

# ============================================================================
# SAVE COMBINED FIGURE: DEPTH AND INTEGRATIONS BY SAMPLE
# ============================================================================

ggsave("figure_depth_samples.pdf", 
       grid.arrange(p, p2, ncol=2, nrow=1, widths = c(1.2, 2)),
       height = 6, 
       width = 15)

# ============================================================================
# FUNCTION FOR CREATING INDIVIDUAL INTEGRATION PLOTS
# ============================================================================
# Function to generate stacked bar plots of integration numbers by circle
# for each wasp-host pair separately
# 
# Parameters:
#   ind = wasp-host pair identifier (e.g., "cotesia_flavipes:chilo_partellus")
#   col_h = vector of colors for samples
#   titre = title of the plot (scientific names)
#   ABC = panel label (A, B, C, etc.)

figure_ind <- function(ind, col_h, titre, ABC) {
  ind <- ggplot(filter(Chim_Int, wasp_h == ind & HIM == "y"), aes(x=Segment, y=int, fill=name)) + 
    geom_col() + 
    xlab("Circles") + ylab("Number of integrations per hhg") + 
    labs(title = ABC, subtitle = titre, fill = "Samples") + 
    theme_classic() + 
    scale_fill_manual(values = col_h) + 
    ylim(0, 8.5) +
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10), 
          plot.subtitle = element_text(face = "italic", size=12))
  
  return(ind)
}

# ============================================================================
# GENERATE INDIVIDUAL PLOTS FOR EACH WASP-HOST PAIR
# ============================================================================

# Panel A: Cotesia flavipes - Chilo partellus
p_Cf_Cp <- figure_ind(ind = "cotesia_flavipes:chilo_partellus", 
                      col_h = c("darkorange1", "orange", "#fee391"), 
                      titre = "Cotesia flavipes - Chilo partellus", 
                      ABC = "A")

# Panel B: Cotesia flavipes - Spodoptera frugiperda
p_Cf_Sf <- figure_ind(ind = "cotesia_flavipes:spodoptera_frugiperda", 
                      col_h = c("darkorange1", "orange", "#fee391"), 
                      titre = "Cotesia flavipes - Spodoptera frugiperda", 
                      ABC = "B")

# Panel C: Cotesia sesamiae kitale - Sesamia calamistis
p_Csk_Sc <- figure_ind(ind = "cotesia_sesamiae_kitale:sesamia_calamistis", 
                       col_h = c("darkgreen", "chartreuse3", "darkseagreen2"), 
                       titre = "Cotesia sesamiae kitale - Sesamia calamistis", 
                       ABC = "C")

# Panel D: Cotesia sesamiae kitale - Busseola fusca
p_Csk_Bf <- figure_ind(ind = "cotesia_sesamiae_kitale:busseola_fusca", 
                       col_h = c("darkgreen", "chartreuse3", "darkseagreen2"), 
                       titre = "Cotesia sesamiae kitale - Busseola fusca", 
                       ABC = "D")

# Panel E: Cotesia sesamiae mombasa - Sesamia calamistis
p_Csm_Sc <- figure_ind(ind = "cotesia_sesamiae_mombasa:sesamia_calamistis", 
                       col_h = c("palevioletred3", "hotpink", "lightpink"), 
                       titre = "Cotesia sesamiae mombasa - Sesamia calamistis", 
                       ABC = "E")

# Panel F: Cotesia sesamiae mombasa - Busseola fusca
p_Csm_Bf <- figure_ind(ind = "cotesia_sesamiae_mombasa:busseola_fusca", 
                       col_h = c("palevioletred3", "hotpink", "lightpink"), 
                       titre = "Cotesia sesamiae mombasa - Busseola fusca", 
                       ABC = "F")

# Panel G: Cotesia icipe - Spodoptera frugiperda
p_Ci <- figure_ind(ind = "cotesia_icipe:spodoptera_frugiperda", 
                   col_h = c("yellow","gold1"), 
                   titre = "Cotesia icipe - Spodoptera frugiperda", 
                   ABC = "G")

# Panel H: Cotesia typhae - Sesamia nonagrioides
p_Ct <- figure_ind(ind = "cotesia_typhae:sesamia_nonagrioides", 
                   col_h = c("#7570B3"), 
                   titre = "Cotesia typhae - Sesamia nonagrioides", 
                   ABC = "H")

# ============================================================================
# SAVE COMBINED FIGURE S4: ALL WASP-HOST PAIRS
# ============================================================================
# Arrange all 8 panels in a 4x2 grid

ggsave("figure_nb_int_all_samples.pdf", 
       grid.arrange(p_Cf_Cp, p_Cf_Sf, p_Csk_Sc, p_Csk_Bf, p_Csm_Sc, p_Csm_Bf, p_Ci, p_Ct, ncol=2, nrow=4),
       height = 12, 
       width = 9)

# ============================================================================
# FIGURE 3A: CORRELATION BETWEEN DEPTH AND INTEGRATION NUMBERS
# ============================================================================
# Scatter plot with linear regression line showing the relationship between
# sequencing depth and number of integrations for HIM circles only

p4 <- ggplot(filter(Chim_Int, HIM == "y"), aes(x = depth, y = int)) +
  geom_point() + 
  xlab("depth(x) per hhg") + ylab("Number of integrations per hhg") +
  theme_bw() +
  # Add linear regression line without confidence interval
  geom_smooth(method = lm, se = FALSE, color = "black") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24))

# Export correlation plot
ggsave("figure_correlation.pdf", 
       p4,
       height = 8, 
       width = 10)
