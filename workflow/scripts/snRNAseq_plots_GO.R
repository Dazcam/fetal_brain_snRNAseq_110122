# -------------------------------------------------------------------------------------
#
#    snRNAseq GO plot  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figure 6

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(janitor)
library(readxl)

##  Info  -----------------------------------------------------------------------------

# 1. Table has 2 sets of col values - https://stackoverflow.com/questions/62613535
# 2. Y = cell-types, X = GO terms
# 3. Circle size = fold-change enrichment, colour FDR / P-value 
# Changing Fold enrichment dot radius:
# https://ggplot2-book.org/scale-other.html / https://stackoverflow.com/questions/56098080 

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/resources/sheets/"
FIG_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/figures/"
CELL_TYPES <- c("FC_ExN_2", "FC_ExN_3", "FC_ExN_4", 
                "GE_InN_2", "Hipp_ExN_4", "Hipp_ExN_6")


##  Load data  ------------------------------------------------------------------------
GO_DATA <- read_excel(paste0(DATA_DIR, "FinalGOterms_to_plot_for_figure6.xlsx"), skip = 2, col_names = F) %>%  
  mutate(across(everything(), ~replace_na(.x, 0))) # Get rid of NAs

colnames(GO_DATA) <- c("Term", "FC_ExN_2-FE", "FC_ExN_2-FDR", "FC_ExN_3-FE", "FC_ExN_3-FDR", 
                       "FC_ExN_4-FE", "FC_ExN_4-FDR", "GE_InN_2-FE", "GE_InN_2-FDR", "Hipp_ExN_3-FE", 
                       "Hipp_ExN_3-FDR", "Hipp_ExN_5-FE", "Hipp_ExN_5-FDR")

# Need to add abbr. for long terms
GO_DATA$Term <- c("GO:0007399~nervous system development", "GO:0048666~neuron development", 
                  "GO:0022008~neurogenesis", "GO:0021872~forebrain generation of neurons", 
                  "GO:0030182~neuron differentiation", "GO:0001764~neuron migration", 
                  "GO:0048667~cell morphogenesis involved in neuron differentiation", 
                  "GO:0007156~homophilic cell adhesion via PAM", 
                  "GO:0031175~neuron projection development", "GO:0061564~axon development", 
                  "GO:0030516~regulation of axon extension", "GO:0050803~regulation of synapse structure or activity", 
                  "GO:0050808~synapse organization", "GO:0007416~synapse assembly", 
                  "GO:0099536~synaptic signaling", "GO:0050804~modulation of synaptic transmission", 
                  "GO:0007186~G-protein coupled receptor signaling pathway", "GO:0042391~regulation of membrane potential", 
                  "GO:0043269~regulation of ion transport", "GO:0006836~neurotransmitter transport", 
                  "GO:0034762~regulation of transmembrane transport", "GO:0007610~behavior", 
                  "GO:0007613~memory")


# Create factor for y-axis order
ORDERED_LIST <- GO_DATA %>% select(Term) %>%  
  #  mutate(Term = gsub(".*~", "", Term)) %>%  # Remove GO numbers
  #  mutate(Term = R.utils::capitalize(Term)) %>%
  pull() %>%
  as_factor()


# Prep data
GO_DATA <- GO_DATA %>% 
  pivot_longer(-Term) %>% 
  separate(name, into = c("cell_type", "Score"), '-') %>% 
  pivot_wider(names_from = Score, values_from = value) %>%
  mutate_at('cell_type', str_replace_all, "_", "-") 
  # mutate(bio_term = gsub(".*~", "", Term)) %>%
  # mutate(bio_term = R.utils::capitalize(bio_term)) %>%

colnames(GO_DATA) <- c("Term", "cell_type", "Fold Enrichment", "FDR")

# Plot
GO_plot <- ggplot(data = GO_DATA, aes(y = factor(Term, level = rev(ORDERED_LIST)), x = cell_type, 
                                      color = -log10(FDR), size = `Fold Enrichment`)) +
  geom_point() +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        #      panel.grid.major = element_blank(), 
        #    panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.text=element_text(size = 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(colour = "#000000")) +
  ylab("") + 
  xlab("") +
  scale_radius(limits = c(1, 6), range = c(1,10))

# Add to axis x if needed: , angle = 45, vjust = 1, hjust = 1

# Save plot
tiff(paste0(FIG_DIR, "Fig_6.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
print(GO_plot)
dev.off()

# jpeg
jpeg(paste0(FIG_DIR, "Fig_6.jpg"), width = 1200, height = 720, 
     units = "px", pointsize = 12, quality = 150)
print(GO_plot)
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------