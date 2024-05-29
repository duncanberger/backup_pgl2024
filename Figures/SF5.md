# Supplementary Figure 5
```{r}
# Load necessary libraries
pacman::p_load(readxl, dplyr, tidyverse, openxlsx, ggtree, phytools, ggplot2, ape, TreeSummarizedExperiment)

# Import table of per loci allele missingness
missing_list <- read.xlsx("missing.xlsx")

# Import table of per isolate clonal complex assignment
cc_list <- read.xlsx("Supp_data_4.xlsx") %>% 
  select(c("PubMLST.ID", "Clonal.complex.(MLST)"))

# Merge tables on common identifier
missing_merge <- merge(cc_list, missing_list, by.x = c("PubMLST.ID"), by.y = c("id"))

# Rename the clonal complex column for clarity
missing_merge <- missing_merge %>% 
  dplyr::rename(CC = "Clonal.complex.(MLST)")

# Import phylogenetic tree and root it at the midpoint
tree3 <- midpoint.root(read.newick("19k_tree.nwk"))

# Plot a simple phylogeny using ggtree
tree4 <- ggtree(tree3, layout = "rectangular", size = 0.1, color = "grey30") 

# Import tip labels
tiplabel <- read.xlsx("tiplabel.xlsx")

# Find the most recent common ancestor (MRCA) for a specific group of tips
PS <- grep("ps", tiplabel$tip_label, value = TRUE)
PS_MRCA <- findMRCA(tree3, tips = PS, type = c("node"))

# Filter data for specific clonal complexes
CC344_label <- missing_merge %>% filter(CC == "344")
CC448_label <- missing_merge %>% filter(CC == "448")
CC344_448 <- rbind(CC344_label, CC448_label)
CC344_448_label <- merge(tiplabel, CC344_448, by.y = "PubMLST.ID", by.x = "tip_label") %>% 
  select("tip_label")
CC344_448_values <- CC344_448_label$tip_label

# Find the MRCA for the combined clonal complexes
CC_MRCA <- findMRCA(tree3, tips = CC344_448_values, type = c("node"))

# Find descendant nodes from the MRCA for both groups
PS_tips <- findDescendant(tree = tree3, node = PS_MRCA, only.leaf = FALSE, self.include = TRUE)
CC_tips <- findDescendant(tree = tree3, node = CC_MRCA, only.leaf = FALSE, self.include = TRUE)

# Extract specific nodes
CC_tips <- CC_tips$NODE_2951
PS_tips <- PS_tips$NODE_1999

# Highlight specific nodes in the phylogeny
tree4 <- tree4 +
  geom_tree(aes(color = ifelse(node %in% PS_tips, "pseudopneumoniae", 
                               ifelse(node %in% CC_tips, "CC", "other")))) +
  scale_color_manual(values = c("pseudopneumoniae" = "red", "CC" = "blue", "other" = "black"))

# Display the phylogeny
tree4

# Plot the phylogeny with an additional panel showing missing alleles data
facet_plot(tree4, panel = "missing alleles", 
           data = missing_merge, 
           geom_segment, 
           mapping = aes(x = 0, xend = na_row, y = y, yend = y, 
                         color = ifelse(CC == 344 | CC == 448, "CC344", "other"))) + 
  scale_color_manual(values = c("CC344" = "blue", "other" = "black", 
                                "pseudopneumoniae" = "red", "CC" = "blue")) +
  theme_tree2() +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  labs(x = "", y = "", colour = "")

```
