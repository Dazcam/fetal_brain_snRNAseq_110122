# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA conditional analyses and gene overlap plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 5 and 6 and 10????? 
#  https://stackoverflow.com/questions/55855426

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(janitor)

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/"
MAGMA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/magma_conditional/"
FIG_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/figures/"
SKENE_CELL_TYPES <- c("skene_CA1", "skene_InN", "skene_MSN", "skene_SS")
BRYOIS_CELL_TYPES <- c("bryois_exCA1", "bryois_exCA3", "bryois_exDG", "bryois_exPFC1")
SIG_CELLS <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'GE_InN_2', 'Hipp_ExN_3', 'Hipp_ExN_5')
SIG_CELLS_EDIT <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'GE-InN-2', 'Hipp-ExN-3', 'Hipp-ExN-5')

##  Load data -------------------------------------------------------------------------
fetal_matrix <- readRDS(paste0(DATA_DIR, 'fetal_overlap_SCZ_magma_P_0.05_matrix.rds')) 
skene_matrix <- readRDS(paste0(DATA_DIR, 'skene_overlap_SCZ_magma_P_0.05_matrix.rds'))
bryois_matrix <- readRDS(paste0(DATA_DIR, 'bryois_overlap_SCZ_magma_P_0.05_matrix.rds'))

colnames(skene_matrix) <- c("FC_ExN_2 (538)", "FC_ExN_3 (476)", "FC_ExN_4 (480)", 
                            "GE_InN_2 (439)", "Hipp_ExN_3 (432)", "Hipp_ExN_5 (398)")
rownames(skene_matrix) <- c('Adult_mus_SS (547)', 'Adult_mus_InN (519)', 'Adult_mus_MSN (561)', 'Adult_mus_CA1 (543)')
colnames(bryois_matrix) <- c("FC_ExN_2 (538)", "FC_ExN_3 (476)", "FC_ExN_4 (480)", 
                            "GE_InN_2 (439)", "Hipp_ExN_3 (432)", "Hipp_ExN_5 (398)")
rownames(bryois_matrix) <- c('Adult_hum_exCA1 (489)', 'Adult_hum_exCA3 (511)', 'Adult_hum_exDG (504)', 'Adult_hum_exPFC1 (516)')

##  Plot ------------------------------------------------------------------------------
# Figure 5A - Adult human vs. fetal gene overlap grid  --------------------------------
fig_5A_data <- reshape2::melt(bryois_matrix) %>%
  arrange(Var1) %>%
  group_by(Var1) %>%
  mutate(Var1 = R.utils::capitalize(Var1))

fig_5A <- fig_5A_data %>%  
  ggplot(aes(x=Var1, y=Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_text(aes(label = value, size = 13)) +
  theme(axis.text.x = element_text(angle = -90, hjust = TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = "#000000", size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "none",
        panel.grid = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  scale_x_discrete(labels = as.factor(gsub("_", "-", as.vector(unique(fig_5A_data$Var1))))) +
  scale_y_discrete(limits = rev(levels(fig_5A_data$Var2)), labels = rev(gsub("_", "-", as.vector(unique(fig_5A_data$Var2))))) +
  xlab("") + 
  ylab("") +
  coord_equal(ratio = 1) 


for (CELL_TYPE in BRYOIS_CELL_TYPES) {
  
  CELL_TYPE_EDIT <- gsub("bryois", "adult_hum", CELL_TYPE)
  
  ##  Load Data  ----------------------------------------------------------------------
  BRYOIS_DATA <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_', CELL_TYPE, '.gsa.out'), header = FALSE) %>%
    row_to_names(row_number = 1) %>% 
    filter(!str_detect(VARIABLE, 'FC_ExN_5|skene|FC_InN_1|GE_InN_1|Thal|bryois')) %>%
    mutate(VARIABLE = paste0(VARIABLE, " no ", CELL_TYPE_EDIT)) %>%
    mutate(VARIABLE = gsub("_", "-", VARIABLE))
  
  assign(CELL_TYPE, BRYOIS_DATA)
  
}

fetal_vs_adult_human <- rbind(bryois_exCA1, bryois_exCA3, bryois_exDG, bryois_exPFC1) 
fetal_vs_adult_human$VARIABLE <- as.factor(fetal_vs_adult_human$VARIABLE)

fig_5B <- fetal_vs_adult_human %>%
  ggplot(aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = '#F8766D', color = 'black') +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  #  ggtitle(toupper(TITLE)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 15),
        axis.title.y = element_text(colour = "#000000", size = 15),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13)) +
  xlab(expression(-log[10](P))) +
  ylab('Cell type')



# Figure 6A - Adult mouse vs. fetal gene overlap grid  --------------------------------
fig_6A_data <- reshape2::melt(skene_matrix) %>%
  arrange(Var1) %>%
  group_by(Var1) %>%
  mutate(Var1 = R.utils::capitalize(Var1))

fig_6A <- fig_6A_data %>%  
  ggplot(aes(x=Var1, y=Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_text(aes(label = value, size = 13)) +
  theme(axis.text.x = element_text(angle = -90, hjust = TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = "#000000", size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "none",
        panel.grid = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  scale_x_discrete(labels = as.factor(gsub("_", "-", as.vector(unique(fig_6A_data$Var1))))) +
  scale_y_discrete(limits = rev(levels(fig_6A_data$Var2)), labels = rev(gsub("_", "-", as.vector(unique(fig_6A_data$Var2))))) +
  xlab("") + 
  ylab("") +
  coord_equal(ratio = 1) 


for (CELL_TYPE in SKENE_CELL_TYPES) {
  
  CELL_TYPE_EDIT <- gsub("skene", "adult_mus", CELL_TYPE)
  
  ##  Load Data  ----------------------------------------------------------------------
  SKENE_DATA <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_', CELL_TYPE, '.gsa.out'), header = FALSE) %>%
    row_to_names(row_number = 1) %>% 
    filter(!str_detect(VARIABLE, 'FC_ExN_5|skene|FC_InN_1|GE_InN_1|Thal|bryois')) %>%
    mutate(VARIABLE = paste0(VARIABLE, " no ", CELL_TYPE_EDIT)) %>%
    mutate(VARIABLE = gsub("_", "-", VARIABLE))
  
  assign(CELL_TYPE, SKENE_DATA)
  
}

fetal_vs_adult_mouse <- rbind(skene_CA1, skene_InN, skene_MSN, skene_SS) 
fetal_vs_adult_mouse$VARIABLE <- as.factor(fetal_vs_adult_mouse$VARIABLE)

fig_6B <- fetal_vs_adult_mouse %>%
  ggplot(aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = '#F8766D', color = 'black') +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  #  ggtitle(toupper(TITLE)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 15),
        axis.title.y = element_text(colour = "#000000", size = 15),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13)) +
  xlab(expression(-log[10](P))) +
  ylab('Cell type')

# Figure 4A - Fetal cell type gene overlap plot  --------------------------------------
fig_10A_data <- reshape2::melt(fetal_matrix) %>%
  arrange(Var1) %>%
  group_by(Var1) %>%
  filter(row_number() >= which(Var1 == Var2)) %>%
  mutate(bold = case_when(Var1 == Var2 ~ TRUE,                  # Lines to omit
                          Var1 != Var2 ~ FALSE))                # if bold boxes
frames = fig_10A_data[fig_10A_data$bold, c("Var1", "Var2")]       # not
frames$Var1 = as.integer(frames$Var1)                           # required
frames$Var2 = rev(as.integer(frames$Var2))                      # and geom_rect

fig_10A <- fig_10A_data %>%  ggplot(aes(x = Var1, y = Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_rect(data = frames, size = 1, fill = NA, colour = "black",
            aes(xmin = Var1 - 0.5, xmax = Var1 + 0.5, ymin = Var2 - 0.5, ymax = Var2 + 0.5)) +
  geom_text(aes(label = value, size = 20)) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = "#000000", size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        axis.title.x = element_text(colour = "#000000", size = 15),
        axis.title.y = element_text(colour = "#000000", size = 15),
        legend.position = "none",
        panel.grid = element_blank()) +
  scale_y_discrete(limits = rev(levels(fig_10A_data$Var2))) +
  xlab("") + 
  ylab("") +
  coord_fixed() 

# Conditional analysis bar chart conditioning fetal cell types on FC-ExN-2
fig_10B_data <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_FC_ExN_2.gsa.out'), header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>% 
  filter(!str_detect(VARIABLE, "FC_ExN_2|FC_ExN_5|skene|FC_InN_1|GE_InN_1|Thal|bryois")) %>%
  mutate(VARIABLE = R.utils::capitalize(VARIABLE)) %>%
  mutate(VARIABLE = gsub("_", "-", VARIABLE))

fig_10B <- ggplot(data = fig_10B_data, aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = '#F8766D', color = 'black') +
  #  geom_vline(xintercept=-log10(0.05/14), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 15),
        axis.title.y = element_text(colour = "#000000", size = 15),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13)) +
  xlab(expression(-log[10](P))) +
  ylab('Cell type') +
  xlim(0, 5)

# Save plots
# Fig 5
tiff(paste0(FIG_DIR, "Fig_5.tiff"), height = 20, width = 40, units='cm', 
     compression = "lzw", res = 300)
plot_grid(fig_5A, fig_5B, align = 'h', labels = 'AUTO', label_size = 18, axis = 'tb')
dev.off()

# Fig 6 
tiff(paste0(FIG_DIR, "Fig_6.tiff"), height = 20, width = 40, units='cm', 
     compression = "lzw", res = 300)
plot_grid(fig_6A, fig_6B, align = 'h', labels = 'AUTO', label_size = 18, axis = 'tb')
dev.off()

# Fig 10 
tiff(paste0(FIG_DIR, "Fig_10.tiff"), height = 20, width = 40, units='cm', 
     compression = "lzw", res = 300)
plot_grid(fig_10A, fig_10B, align = 'h', labels = 'AUTO', label_size = 18, axis = 'tb')
dev.off()

# Jpegs
# Fig 5
jpeg(paste0(FIG_DIR, "Fig_5.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
plot_grid(fig_5A, fig_5B, align = 'h', labels = 'AUTO', label_size = 18, rel_heights = c(1, 0.1), axis = 'tb')
dev.off()

# Fig 6
jpeg(paste0(FIG_DIR, "Fig_6.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
plot_grid(fig_6A, fig_6B, align = 'h', labels = 'AUTO', label_size = 18, rel_heights = c(1, 0.1), axis = 'tb')
dev.off()

# Fig 10 
jpeg(paste0(FIG_DIR, "Fig_10.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
plot_grid(fig_10A, fig_10B, align = 'h', labels = 'AUTO', label_size = 18, rel_heights = c(1, 0.1), axis = 'tb')
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------