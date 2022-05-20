# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA GO plots 
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figure 7

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(reshape2) # For melt

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/magma_GO/"
FIG_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/figures/'
SIG_CELLS <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'GE_InN_2', 'Hipp_ExN_3', 'Hipp_ExN_5')

GO_DF <- read.table(paste0(DATA_DIR, 'magma_GO.gsa.out'), header = FALSE) %>%
  janitor::row_to_names(row_number = 1)

# Need to recode GO Terms as Magma truncated them - https://stackoverflow.com/questions/65178820
# COMPLETE_TERMS <- c("GO:0007399~nervous" = "GO:0007399~nervous system development", 
#                     "GO:0048666~neuron" = "GO:0048666~neuron development", 
#                     "GO:0022008~neurogenesis" = "GO:0022008~neurogenesis",
#                     "GO:0021872~forebrain" = "GO:0021872~forebrain generation of neurons", 
#                     "GO:0030182~neuron" = "GO:0030182~neuron differentiation", 
#                     "GO:0001764~neuron" = "GO:0001764~neuron migration", 
#                     "GO:0048667~cell" = "GO:0048667~cell morphogenesis involved in neuron differentiation", 
#                     "GO:0007156~homophilic" = "GO:0007156~homophilic cell adhesion via PAM", 
#                     "GO:0031175~neuron" = "GO:0031175~neuron projection development", 
#                     "GO:0061564~axon" = "GO:0061564~axon development", 
#                     "GO:0030516~regulation" = "GO:0030516~regulation of axon extension", 
#                     "GO:0050803~regulation" = "GO:0050803~regulation of synapse structure or activity", 
#                     "GO:0050808~synapse" = "GO:0050808~synapse organization", 
#                     "GO:0007416~synapse" = "GO:0007416~synapse assembly", 
#                     "GO:0099536~synaptic" = "GO:0099536~synaptic signaling", 
#                     "GO:0050804~modulation" = "GO:0050804~modulation of synaptic transmission", 
#                     "GO:0007186~G-protein" = "GO:0007186~G-protein coupled receptor signaling pathway", 
#                     "GO:0042391~regulation" = "GO:0042391~regulation of membrane potential", 
#                     "GO:0043269~regulation" = "GO:0043269~regulation of ion transport", 
#                     "GO:0006836~neurotransmitter" = "GO:0006836~neurotransmitter transport", 
#                     "GO:0034762~regulation" = "GO:0034762~regulation of transmembrane transport", 
#                     "GO:0007610~behavior" = "GO:0007610~behavior", 
#                     "GO:0007613~memory" = "GO:0007613~memory")

COMPLETE_TERMS <- c("GO:0007399~nervous" = "Nervous system development", 
                    "GO:0048666~neuron" = "Neuron development", 
                    "GO:0022008~neurogenesis" = "Neurogenesis",
                    "GO:0021872~forebrain" = "Forebrain generation of neurons", 
                    "GO:0030182~neuron" = "Neuron differentiation", 
                    "GO:0001764~neuron" = "Neuron migration", 
                    "GO:0048667~cell" = "Cell morphogenesis involved in neuron differentiation", 
                    "GO:0007156~homophilic" = "Homophilic cell adhesion via PAM", 
                    "GO:0031175~neuron" = "Neuron projection development", 
                    "GO:0061564~axon" = "Axon development", 
                    "GO:0030516~regulation" = "Regulation of axon extension", 
                    "GO:0050803~regulation" = "Regulation of synapse structure or activity", 
                    "GO:0050808~synapse" = "Synapse organization", 
                    "GO:0007416~synapse" = "Synapse assembly", 
                    "GO:0099536~synaptic" = "Synaptic signaling", 
                    "GO:0050804~modulation" = "Modulation of synaptic transmission", 
                    "GO:0007186~G-protein" = "G-protein coupled receptor signaling pathway", 
                    "GO:0042391~regulation" = "Regulation of membrane potential", 
                    "GO:0043269~regulation" = "Regulation of ion transport", 
                    "GO:0006836~neurotransmitter" = "Neurotransmitter transport", 
                    "GO:0034762~regulation" = "Regulation of transmembrane transport", 
                    "GO:0007610~behavior" = "Behavior", 
                    "GO:0007613~memory" = "Memory")

for (CELL_TYPE in SIG_CELLS) {

GO_DF_FILT <- filter(GO_DF, grepl(CELL_TYPE, FULL_NAME)) %>%
  mutate(Term = gsub(paste0("^[^", CELL_TYPE, "_]*", CELL_TYPE,"_"), "", FULL_NAME)) %>%
  select(Term, P) %>%
  mutate(Full_Term = COMPLETE_TERMS[as.character(Term)]) %>%
  mutate(cell_type = rep(CELL_TYPE, length(Full_Term))) %>%
  mutate_at('cell_type', str_replace_all, "_", "-") 

assign(paste0(CELL_TYPE, '_df'), GO_DF_FILT)

}

# Combine filtered dfs
GO_DF_GROUP <- rbind(FC_ExN_2_df, FC_ExN_3_df, FC_ExN_4_df, 
      GE_InN_2_df, Hipp_ExN_3_df, Hipp_ExN_5_df) 


# Create factor for y-axis order
ORDERED_LIST <- GO_DF_GROUP %>% select(Full_Term) %>%  
  #  mutate(Term = gsub(".*~", "", Term)) %>%  # Remove GO numbers
  #  mutate(Term = R.utils::capitalize(Term)) %>%
  pull() %>% 
  unname() %>%
  as_factor()

# levels(ORDERED_LIST) <- c("GO:0007399~nervous system development", "GO:0048666~neuron development", 
#                           "GO:0022008~neurogenesis", "GO:0021872~forebrain generation of neurons", 
#                           "GO:0030182~neuron differentiation", "GO:0001764~neuron migration", 
#                           "GO:0048667~cell morphogenesis involved in neuron differentiation", 
#                           "GO:0007156~homophilic cell adhesion via PAM", 
#                           "GO:0031175~neuron projection development", "GO:0061564~axon development", 
#                           "GO:0030516~regulation of axon extension", "GO:0050803~regulation of synapse structure or activity", 
#                           "GO:0050808~synapse organization", "GO:0007416~synapse assembly", 
#                           "GO:0099536~synaptic signaling", "GO:0050804~modulation of synaptic transmission", 
#                           "GO:0007186~G-protein coupled receptor signaling pathway", "GO:0042391~regulation of membrane potential", 
#                           "GO:0043269~regulation of ion transport", "GO:0006836~neurotransmitter transport", 
#                           "GO:0034762~regulation of transmembrane transport", "GO:0007610~behavior", 
#                           "GO:0007613~memory")

levels(ORDERED_LIST) <- c("Nervous system development", "Neuron development", 
                          "Neurogenesis", "Forebrain generation of neurons", 
                          "Neuron differentiation", "Neuron migration", 
                          "Cell morphogenesis involved in neuron differentiation", 
                          "Homophilic cell adhesion via PAM", 
                          "Neuron projection development", "Axon development", 
                          "Regulation of axon extension", "Regulation of synapse structure or activity", 
                          "Synapse organization", "Synapse assembly", 
                          "Synaptic signaling", "Modulation of synaptic transmission", 
                          "G-protein coupled receptor signaling pathway", "Regulation of membrane potential", 
                          "Regulation of ion transport", "Neurotransmitter transport", 
                          "Regulation of transmembrane transport", "Behavior", 
                          "Memory")


GO_PLOT <- ggplot(GO_DF_GROUP, aes(x = -log10(as.numeric(P)), y = Full_Term)) +
  geom_bar(stat = "identity", color = 'black', fill = '#F8766D') +
  facet_grid(~cell_type) +
  geom_vline(xintercept=-log10(0.05/121), linetype = "dashed", color = "black") +
  #  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 16),
        axis.title.y = element_text(colour = "#000000", size = 16),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, colour = "black")) +
  xlab(expression(-log[10](P))) +
  ylab('Term') +
  xlim(0, 10.5) +
  scale_x_continuous(breaks= scales::pretty_breaks()) +
  scale_y_discrete(limits = rev(levels(ORDERED_LIST)))
  

# Tiff
tiff(paste0(FIG_DIR, "Fig_7.tiff"), height = 20, width = 40, units='cm', 
     compression = "lzw", res = 300)
GO_PLOT
dev.off()

# Jpeg
jpeg(paste0(FIG_DIR, "Fig_7.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
GO_PLOT
dev.off()



# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------