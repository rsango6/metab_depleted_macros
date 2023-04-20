library(reshape2)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(ggrepel)
library(pheatmap)
library(plyr)
library(ggsignif)


##### for stats #####
library(labelled)   # labeling data
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(GGally)     # advanced plot
library(car)        # useful for anova/wald test
library(Epi)        # easy getting CI for model coef/pred
library(lme4)       # linear mixed-effects models
library(lmerTest)   # test for linear mixed-effects models
library(emmeans)    # marginal means
library(multcomp)   # CI for linear combinations of model coef
library(geepack)    # generalized estimating equations
library(ggeffects)  # marginal effects, adjusted predictions
library(gt) 
####################


my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=12),
                  title = element_text(size=12, face="bold", vjust=2),
                  panel.background = element_rect(fill = 'gray99',
                                                  colour = "black", 
                                                  linewidth=0.5),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=1.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line(),
                  panel.grid.major = element_line(colour = "gray40", linetype="dotted"),
                  panel.grid.minor = element_line(colour = "gray40", linetype="dashed"),
                  legend.justification=c(1,1),
                  legend.position=c(1,1),
                  legend.title=element_blank(),
                  legend.text = element_text(size = 14))

dir = '/Users/rokosango/PhD/SFB_Andrea/'

df = read.csv(paste0(dir, '20230301_VogelAndrea copy.csv'))


df = df[df$X != "sample19",] #remove sample19


dfNum = df %>%
  dplyr::select(-X, -mouse, -Group, -Day) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::select(where(is.numeric))

dfNum = dfNum[ , colSums(is.na(dfNum))==0] #remove cols with NA's.

logm = log(dfNum, 2)

logm = scale(logm, center = T, scale = T)
rownames(logm) = df$X
logm = as.data.frame(logm)

logm = logm %>%
  dplyr::mutate(Group = df$Group) %>%
  dplyr::mutate(Day = df$Day) %>%
  dplyr::mutate(Sample = df$X) %>%
  dplyr::mutate(Mouse = df$mouse) %>%
  relocate(Group, .before = Spermidine) %>%
  relocate(Day, .before = Group) %>%
  relocate(Mouse, .before = Day) %>%
  relocate(Sample, .before = Mouse) %>%
  dplyr::arrange(Mouse)

logm$Day = factor(logm$Day, levels = c("d4", "d11", "d21 (EP)"))
logm$Group = as.factor(logm$Group)

MetsLME = list()
LME_Results = list()
Anova_res = list()
SignificantMets = list()

mets_names = names(logm)

for (i in 5:ncol(logm)) {
  
  MetsLME[[i]] = cbind(logm$Sample, logm$Mouse,
                      logm$Group, logm$Day,
                     logm[i])
  names(MetsLME)[i] = mets_names[i]
  names(MetsLME[[i]]) = c("Sample", "Mouse", "Group", "Day", "ZScore")
  

}

MetsLME[1:4] = NULL

#We can test whether the mean response is constant over time and groups by testing the null hypothesis 
#that all the regression coefficients used to model Day and Group are simultaneously equal to zero
# H0 : beta 1 = beta 2 = beta 3 = beta 4 = 0 


for (i in 1:length(MetsLME)) {
  
  LME_Results[[i]] =  lmer(ZScore ~ Group * Day + (1 | Mouse), data = MetsLME[[i]])
  names(LME_Results)[i] = names(MetsLME)[i]
  
  Anova_res[[i]] = Anova(LME_Results[[i]], type = 3)
  names(Anova_res)[i] = names(LME_Results)[i]
  
  
}

for (i in names(Anova_res)) {
    if (any(Anova_res[[i]]$`Pr(>Chisq)` < 0.05)) {
      SignificantMets[[i]] = Anova_res[[i]]
    }
}

SigMetsUnlisted = data.frame(values = unlist(SignificantMets))
SigMetsUnlisted$Metabolites = rownames(SigMetsUnlisted)

SigMetsUnlisted = SigMetsUnlisted %>% #weird: need to do filter and mutate here separately else error
  dplyr::filter(str_detect(Metabolites, ">Chisq")) %>%
  dplyr::mutate(Interaction = rep(rownames(SignificantMets[["Glucose"]]),
                                  times = nrow(SigMetsUnlisted) / nrow(SignificantMets[["Glucose"]]))) %>%
  dplyr::arrange(values, ascending = T) %>%

GroupDayOnly = SigMetsUnlisted %>%
  dplyr::filter(Interaction == "Group:Day") %>%
  dplyr::filter(values < 0.05) %>%
  dplyr::mutate(Signif = case_when(
    values < 0.001 ~ "***",
    values > 0.001 & values < 0.01 ~ "**",
    values > 0.01 & values < 0.05 ~ "*",
    values > 0.05 ~ "ns"))

plotlist = list()

setwd("~/PhD/SFB_Andrea/plots")

for (i in unique(GroupDayOnly$Metabolites)){
  
  
  plotlist[[i]] = ggplot(MetsLME[[i]], aes(Group, ZScore, fill = Day)) +
    geom_boxplot() +
     annotate(
       "text", label = GroupDayOnly$Signif[GroupDayOnly$Metabolites == i],
       x = 1.5, y = max(MetsLME[[i]]$ZScore) + 0.5, size = 6, colour = "black", fontface = "bold") +
    labs(y = "Met Z Score", title = i) +
    theme_bw()
  #ggsave(filename = paste("Boxplot", i, ".pdf", sep = "_"), plot = plotlist[[i]],
  #       device = "pdf", width = 3, height = 3, units = "in")

}

Pairwise_comparisons_list = list()
Pairwise_p_vals = list()

for (i in unique(GroupDayOnly$Metabolites)) {
  for (j in 1:3) {
  
  Pairwise_comparisons_list[[i]] = emmeans(LME_Results[[i]], pairwise ~ Group | Day)
  Pairwise_p_vals[[i]][j] = tidy(Pairwise_comparisons_list[[i]][["contrasts"]][j])[8]
  
  }
}


PairwiseMetsUnlisted = data.frame(values = unlist(Pairwise_p_vals))
PairwiseMetsUnlisted$Metabolites = rownames(PairwiseMetsUnlisted)
PairwiseMetsUnlisted$Contrast = rep(c("d4", "d11", "d21 (EP)"), 
                                    times = nrow(PairwiseMetsUnlisted) / length(Pairwise_p_vals[["stearoylcarnitine"]]))
  
PairwiseMetsUnlisted = PairwiseMetsUnlisted %>% 
  dplyr::filter(values < 0.05) %>%
  dplyr::arrange(Metabolites, ascending = T) %>%
  dplyr::mutate(Metabolites = gsub("[[:digit:]]+", "", Metabolites)) %>% #remove a digit at the end of alphanumeric string
  dplyr::mutate(Signif = case_when(
    values < 0.001 ~ "***",
    values > 0.001 & values < 0.01 ~ "**",
    values > 0.01 & values < 0.05 ~ "*",
    values > 0.05 ~ "ns"))

row.names(PairwiseMetsUnlisted) = 1:nrow(PairwiseMetsUnlisted)

unique(PairwiseMetsUnlisted$Metabolites) #32

setwd("~/PhD/SFB_Andrea/plots")

for (i in unique(PairwiseMetsUnlisted$Metabolites)){
  
  annotation <- data.frame(
    x = c(0.7,0.7,0.7),
    y = c(max(MetsLME[[i]]$ZScore) + 0.55,
          max(MetsLME[[i]]$ZScore) + 0.35,
          max(MetsLME[[i]]$ZScore) + 0.15),
    label = c(paste0("d4:", 
                     subset(PairwiseMetsUnlisted, Contrast == 'd4' & Metabolites == i)[4]), 
              paste0("d11:", 
                     subset(PairwiseMetsUnlisted, Contrast == 'd11' & Metabolites == i)[4]), 
              paste0("d21: ", 
                     subset(PairwiseMetsUnlisted, Contrast == 'd21 (EP)' & Metabolites == i)[4])),
    color = factor(c("#F8766D", "#00BA38", "#619CFF"), levels = c("#F8766D", "#00BA38", "#619CFF")),
    size= 2, fontface="italic",
    Group = c("d4", "d11", "d21 (EP)"))
  
  rownames(annotation) = c("d4", "d11", "d21 (EP)")
  
  if (nrow(subset(PairwiseMetsUnlisted, Contrast == 'd4' & Metabolites == i)[4]) == 0) {
    annotation = annotation[!row.names(annotation) %in% 'd4',]
  }
  
  if (nrow(subset(PairwiseMetsUnlisted, Contrast == 'd11' & Metabolites == i)[4]) == 0) {
    annotation = annotation[!row.names(annotation) %in% 'd11',]
  }
  
  if (nrow(subset(PairwiseMetsUnlisted, Contrast == 'd21 (EP)' & Metabolites == i)[4]) == 0) {
    annotation = annotation[!row.names(annotation) %in% 'd21 (EP)',]
  }
  
  plotlist[[i]] = ggplot(MetsLME[[i]], aes(Group, ZScore, fill = Day)) +
    geom_boxplot() +
    geom_text(data = annotation, aes(x = x, y = y, label = label, 
                                     size = size, fontface = fontface), 
              inherit.aes = FALSE, show.legend = FALSE) + #inherit.aes = FALSE else previous ggplot aes is inherited!!!
    annotate(
      "text", label = paste0("Group*Day: ", GroupDayOnly$Signif[GroupDayOnly$Metabolites == i]),
      x = 2.2, y = max(MetsLME[[i]]$ZScore) + 0.5, size = 3.5, colour = "black", fontface = "bold") +
    labs(y = "Met Z Score", title = i) +
    theme_bw()
  ggsave(filename = paste("Boxplot", i, ".pdf", sep = "_"), plot = plotlist[[i]],
        device = "pdf", width = 4.5, height = 4, units = "in")
  
}

setwd("~/PhD/SFB_Andrea/")
write.csv(GroupDayOnly, 'GroupDayOnlyPvals.csv')
write.csv(PairwiseMetsUnlisted, 'PairwiseMetsUnlisted.csv')

#####################################


fun4everything = function(df, day, title) {
  
  suppressWarnings({

  dfNum = df %>%
    dplyr::filter(Day == day) %>%
    dplyr::select(-X, -mouse, -Group, -Day, -X.1) %>%
    dplyr::mutate_if(is.character, as.numeric) %>%
    dplyr::select(where(is.numeric))
  
  dfNum = dfNum[ , colSums(is.na(dfNum))==0] #remove cols with NA's.
  
  logm = log(dfNum, 2)
  
  logm = scale(logm, center = T, scale = T)
  rownames(logm) = df$X[df$Day == day]
  t.logm = t(logm)
  
  t.logm = as.data.frame(t.logm)
  
  pca <- prcomp(logm, center = F, scale = F)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
  
  pcaresults <- summary(pca)
  
  scree.data <- as.data.frame(pcaresults$importance) # eigenvalues / standard deviations
  score.data <- as.data.frame(pcaresults$x) # coordinates of the samples (i.e., observations) 
  loading.data <- as.data.frame(pcaresults$rotation) # loading scores for individual variables
  
  score.data <- score.data[, 1:2]
  score.data$Group <- df$Group[df$Day == day]
  rownames(score.data) = df$X[df$Day == day]
  perc.var = 100 * pcaresults[["importance"]][2, ]
  
  
  scores_plot <- ggplot(score.data, aes(x = PC1, y = PC2)) +
    geom_point(aes(shape = Group, colour = Group)) + 
    geom_label_repel(label = rownames(score.data)) +
    stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = Group)) +
    xlab(paste("PC1:", perc.var[1], "%")) + 
    ylab(paste("PC2:", perc.var[2], "%")) +
    my.theme +
    ggtitle(title)
  
  ######################################
  
  pvals <- apply(t.logm, 1, function(x) {
    t.test(x[1:4], x[5:8])$p.value
  })
  
  pvals = pvals[order(pvals)]
  #pvals = pvals[pvals < 0.05]
  pvals_adj = p.adjust(pvals, method = "BH")
  
  if (day == "d21 (EP)") {
    df_filtered <- cbind(t.logm[, 1:3], t.logm[, 4:7])
    Group <- rep(c("Iso","Depletion"), times = c(3,4))
  } else {
    df_filtered <- cbind(t.logm[, 1:4], t.logm[, 5:8])
    Group <- rep(c("Iso","Depletion"), each = 4)
  }

  df_filtered <- df_filtered[names(pvals_adj), ]
  
  df_filtered = as.data.frame(df_filtered)
  
  annot <- data.frame(Group = Group)
  rownames(annot) <- colnames(t.logm)
  
  
  heatmap = pheatmap(t.logm,
            color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
            border_color = "gray40",
            scale = "row",
            angle_col = 45,
            fontsize_row = 5,
            cluster_cols = F,
            cellwidth = 15,
            annotation_col =  annot,
            annotation_names_col = F,
            main = title
           )
  
  heatmap_filtered = pheatmap(df_filtered,
                     color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
                     border_color = "gray40",
                     scale = "row",
                     angle_col = 45,
                     fontsize_row = 10,
                     cluster_cols = F,
                     cellwidth = 15,
                     annotation_col = annot,
                     annotation_names_col = F,
                     main = title
  )
  
  
  return_list = list(scores_plot, loadings_plot, loadings.sig, t.logm, heatmap, df_filtered, 
                     heatmap_filtered, pvals, pvals_adj)
  names(return_list) = c("scores_plot", "loading_plot", "loadings_df", "scaled_log_dataframe", "heatmap_all_mets", 
                         "filtered4significance_df", "filtered_heatmap", "pvals", "adj_p_vals")
  return(return_list)
  
  })
  
}

day_4_results = fun4everything(df, "d4", "Day 4")
day_11_results = fun4everything(df, "d11", "Day 11")
day_21_results = fun4everything(df, "d21 (EP)", "Day 21")


MSEADay4 = read.csv(paste0(dir, "MSEA_Day4.csv"))
MSEADay11 = read.csv(paste0(dir, "MSEA_Day11.csv"))
MSEADay21 = read.csv(paste0(dir, "MSEA_Day21_new.csv"))


MSEAfun = function(df, title, list) {

  df = df %>%
    mutate(BH = p.adjust(P.value, method = "BH")) %>%
    filter(BH < 0.1) %>%
    mutate(EnrichmentRatio = Hits / Total)
  
  
  g <- ggplot(df, aes(x = EnrichmentRatio, y = reorder(str_wrap(Metabolite.Set,
                                                               width = 20), -BH),
                   size = Hits, fill = BH)) +
    geom_point(shape=21, stroke = 1) +
    scale_size(range = c(2, 6)) + 
    scale_fill_gradient(
      low = "#D6604D",
      high = "#4393C3",
      aesthetics = "fill")  + 
    labs(y = "", title = title) + 
    theme(axis.text = element_text(face = 'bold',
                                   size = 10))
  
  
  
  list[[length(list) + 1]] = g #append to the list with double brackets, not single!
  names(list)[length(list)] = "MSEAplot"
  
  return(list)
    
}


day_4_results = MSEAfun(MSEADay4, "MSEA Day 4", day_4_results)
day_11_results = MSEAfun(MSEADay11, "MSEA Day 11", day_11_results)
day_21_results = MSEAfun(MSEADay21, "MSEA Day 21", day_21_results)

