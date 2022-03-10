### ===========================================================================================
### R code to plot the mitochondrial genome phylogeny (Figure 1b) in
### "A plea for a conservative approach to the report of new species records from eDNA surveys"
### by Oscar Puebla and Martin Helmkampf (2022)
### ===========================================================================================

library(ggtree)
library(treeio)
library(tidyverse)
library(ggtext)

setwd("/Users/martin/Documents/Projects/2_Other/egem/3_phylo")


### Set color default
clr_neutral <- rgb(0.2, 0.2, 0.2)


### Read in and root tree
tree_file <- read.tree("raxml/egem_mtg2_f_gtr_500.raxml.support")
tree_rooted <- root(phy = tree_file, outgroup = "PL17_160floflo")


### Shorten root branch (= last branch)
tree_rooted$edge.length[1] <- tree_rooted$edge.length[1] * 0.5


### Convert to ggtree object, prepare labels and add support categories
tree_data <- ggtree(tree_rooted) %>%
  .$data %>%
  mutate(
    spec = ifelse(isTip, str_sub(label, -6, -4), "ungrouped"),
    loc = ifelse(isTip, str_sub(label, -3, -1), "ungrouped"),
    loc = replace(loc, label == "Ref_puemtg", "pan"),
    support = as.numeric(label),
    support_class = cut(support, breaks = c(0, 50, 70, 90, 100),
                                 labels = c("(0,50]", "(50,70]", "(70,90]", "(90,100]")),
    branch_type = if_else(node %in% 234, "broken", "whole")
  )


### Initialize tree
tree <- ggtree(tree_data)


### Plot tree
(tree <-
  ggtree(tree_data, layout = "rectangular", aes(color = spec, linetype = branch_type)) +
  geom_tippoint(aes(
    color = spec,
    shape = loc,
    fill = after_scale(color)
  )) +
  geom_nodepoint(
    data = tree %>% filter(!isTip, support_class != "(0,50]"),   # Apply to nodes with support >50 only
    aes(
      fill = support_class,
      size = support_class
    ),
    shape = 21,
    color = clr_neutral
  ) +
  # Add scale bar
  geom_treescale(
    width = 0.001,
    x = 0, y = 30,
    offset = -3, fontsize = 3,
    color = "gray70"
  ) +
  # Set colors and shapes
  scale_color_manual(
    values = c(
      pue = "#E17366",
      nig = "#000000",
      uni = "#CCCCCC",
      gem = "dodgerblue",
      gum = "#E3A258",
      may = "#7EA7C2",
      flo = "#A39CCF"
    ),
    labels = c("*H. puella*", "*H. nigricans*", "*H. unicolor*", "*H. gemma*", "*H. gummigutta*", "*H. maya*", "*H. floridae*")
  ) +
  scale_shape_manual(
    values = c(
      bel = 21,
      flo = 24,
      hon = 22,
      pan = 23,
      mtg = 04
    ),
    labels = c("Belize", "Florida", "Honduras", "Panama", "unknown")
  ) +
  scale_fill_manual(
    values = c(
      `(50,70]`  = "white",
      `(70,90]`  = "gray",
      `(90,100]` = "black"
    ),
    drop = FALSE
  ) +
  scale_size_manual(
    values = c(
      `(50,70]`  = 1.5,
      `(70,90]`  = 1.5,
      `(90,100]` = 1.5
    ),
    na.value = 0,
    drop = FALSE
  ) +
  scale_linetype_manual(values = c(whole = 1, broken = 3), guide = "none") +
  # Add legend
  guides(
    color = guide_legend(title = "Species", title.position = "top", nrow = 4, order = 1, override.aes = list(shape = NA, size = 0.6)),
    shape = guide_legend(title = "Site", title.position = "top", nrow = 5, order = 2),
    fill = guide_legend(title = "Support", title.position = "top", nrow = 4),
    size = guide_legend(title = "Support", title.position = "top", nrow = 4)
  ) +
  # Add asterisks
  annotate(geom = "text",
           label = "*", 
           x = c(tree_data$x[tree_data$label == "ENA_gemmtg"] + 0.0005, tree_data$x[tree_data$label == "Ref_puemtg"] + 0.0005),
           y = c(0.5, 100.5), color = "gray20"
    ) +
  # Set theme
  theme_void() +
  theme(
    legend.position = "bottom", legend.justification = "left", legend.title.align = 0,
    legend.text = element_markdown(color = "gray20", margin = margin(r = 30, unit = "pt")),   # space between legends
    legend.title = element_text(color = "gray20"),
    legend.key = element_rect(size = 1, color = "white")
  )
)


### Print tip labels and node numbers
# (tree <- tree + geom_tiplab(size = 1.5, hjust = -0.3))
#
# (tree <- tree + geom_text2(aes(subset =! isTip, label = node), size = 2.5, hjust = -0.3))


### Turn into ggplot object, position tree and legend
(gt <- ggplot() +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = 0
  ) +
  annotation_custom(
    grob = ggplotGrob(tree + theme(legend.position = "none")),
    ymin = 0.05, ymax = 1.23,
    xmin = 0.02, xmax = 1
  ) +
  annotation_custom(
    grob = cowplot::get_legend(tree),
    ymin = -0.22, ymax = -0.02,
    xmin = 0.065, xmax = 1
  ) +
  theme_void()
)


### Save to file
ggsave(
  plot = gt,
  filename = "tree_plots/egem_mtg2_f_gtr_500_v1.pdf",
  width = 6,
  height = 9,
  device = cairo_pdf,
  bg = "transparent"
)


### *** Note ***
### Plotting haplotype network object creates a list that could not be converted into a ggplot object
### Haplotype network and phylogenetic tree plots were combined manually in Inkscape
