# Supplementary Figure 4
```{r}
# Load libraries
pacman::p_load(readxl,dplyr,tidyverse,openxlsx,ggrepel)

# Import table of per loci allele missingness
missing <- read.xlsx("missing.xlsx")

# Import table of per isolate clonal complex assignment
cc_list <- read.xlsx("Supp_data_4.xlsx") %>% select(c("PubMLST.ID", "Clonal.complex.(MLST)"))

# Merge tables
missing_merge <- merge(cc_list, missing, by.x=c("PubMLST.ID"), by.y=c("id"))

# Rename column for clarity
missing_merge <- missing_merge %>% rename(CC="Clonal.complex.(MLST)")

# Get mean missing counts per clonal complex
missing_merge_counts <- missing_merge %>% 
  group_by(CC) %>% 
  select(CC, na_row) %>% 
  summarise(count=n(), mean=mean(na_row), sstdev=sd(na_row), min=min(na_row), max=max(na_row))

# Filter clonal complexes with counts less than 10
missing_merge_counts_10p <- missing_merge_counts %>% filter(count>=10)

# Drop rows with missing values
missing_merge_counts_10p2 <- missing_merge_counts_10p %>% drop_na()

# Plotting

# Order data by mean count of missing alleles
missing_merge_counts_10p2 <- missing_merge_counts_10p2[order(missing_merge_counts_10p2$mean), ]

# Create a label column based on mean count condition
missing_merge_counts_10p2$label <- NA
missing_merge_counts_10p2$label <- ifelse(missing_merge_counts_10p2$mean > 25, missing_merge_counts_10p2$CC, missing_merge_counts_10p2$label)

# Order levels of CC for plotting
order <- missing_merge_counts_10p2$CC
missing_merge_counts_10p2$CC <- ordered(missing_merge_counts_10p2$CC, levels = order)

# Adjust mean by standard deviation
missing_merge_counts_10p2$m_sd <- missing_merge_counts_10p2$mean - missing_merge_counts_10p2$sstdev
missing_merge_counts_10p2$m_sd <- ifelse(missing_merge_counts_10p2$m_sd < 0, 0, missing_merge_counts_10p2$m_sd)



# Create first plot object
one <- ggplot(data=missing_merge_counts_10p2) + 
  geom_pointrange(aes(x=mean, y=CC, xmin=m_sd, xmax=mean+sstdev)) + 
  theme_bw() +
  xlab("Mean number of missing alleles") +
  ylab("Clonal complex") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(), legend.position = "none",
        axis.text.x=element_text(color="black", size=7),
        axis.text.y=element_text(color="black", size=3),
        axis.title.y = element_text(face="bold", color="black", size=8),
        axis.title.x = element_text(face="bold", color="black", size=8))

# Create second plot object with labels
two <- ggplot(data=missing_merge_counts_10p2 %>% filter(is.na(label) == FALSE)) + 
  geom_pointrange(aes(x=mean, y=CC, xmin=m_sd, xmax=mean+sstdev)) + 
  theme_bw() +
  xlab("Mean number of missing alleles") +
  ylab("Clonal complex") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(), legend.position = "none",
        axis.text.x=element_text(color="black", size=7),
        axis.text.y=element_text(color="black", size=3),
        axis.title.y = element_text(face="bold", color="black", size=8),
        axis.title.x = element_text(face="bold", color="black", size=8)) +
  geom_label(aes(label = label, y=CC, x=mean+sstdev), size = 3, position = position_dodge(0.3)) +
  scale_x_continuous(breaks = c(0,25,50,75,100)) +
  theme(axis.text.y=element_blank())

# Create final visualization combining two plots
ggplot(data=missing_merge_counts_10p2) + 
  geom_pointrange(aes(x=mean, y=CC, xmin=m_sd, xmax=mean+sstdev)) + 
  theme_bw() +
  xlab("Mean number of missing alleles") +
  ylab("Clonal complex") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(), legend.position = "none",
        axis.text.x=element_text(color="black", size=7),
        axis.text.y=element_text(color="black", size=3),
        axis.title.y = element_text(face="bold", color="black", size=8),
        axis.title.x = element_text(face="bold", color="black", size=8)) +
  theme(panel.border = element_blank()) +  
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y=element_blank()) +
  annotation_custom(
    ggplotGrob(two), 
    xmin = 30, xmax = 70, ymin = 50, ymax =200) +
  scale_x_continuous(expand=c(0,0), breaks = c(0,25,50,75,100)) +
  scale_y_discrete(expand = c(0,0)) +
  annotate("rect", xmin=28, xmax=55, ymin=297, ymax=305, fill = "grey60", alpha=0.4, col = "black", lty=1) +
  geom_segment(x=45, yend=200, y=297, arrow = arrow(length = unit(0.5, "cm")))

```
