library(tidyverse)

# Barplot of number of genes in orthogroups
# load the data
gc <- read_tsv("Orthogroups.GeneCount.tsv") %>%
  select(-Total) %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Genes") %>%
  mutate(Species = recode(Species,
                          Osativa_indica_WRKY = "Rice",
                          Sbicolor_WRKY = "Sorghum",
                          Zmays_WRKY = "Maize"),
         Species = factor(Species, levels = c("Maize", "Sorghum", "Rice")),  # Rice ends up on top
         Group = case_when(Genes == 0 ~ "0 (absent)",
                           Genes >= 4 ~ "4 or more",
                           TRUE ~ as.character(Genes)),
         Group = factor(Group, levels = c("4 or more", "3", "2", "1", "0 (absent)"))) %>%
  count(Species, Group)

# define colours
pal <- c("0 (absent)" = "#E5E5E5", "1" = "#C6DBEF", "2" = "#6BAED6",
         "3" = "#2171B5", "4 or more" = "#08306B")

# create plot
ggplot(gc, aes(x = n, y = Species, fill = Group)) +
  geom_col(colour = "white", linewidth = 0.6, width = 0.65) +
  geom_text(aes(label = n, group = Group,
                colour = Group %in% c("3", "4 or more")),
            position = position_stack(vjust = 0.5),
            size = 3.6, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey20")) +
  scale_fill_manual(values = pal, breaks = rev(levels(gc$Group)),
                    name = "Genes in the orthogroup") +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 91),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 91)) +
  labs(x = "Number of orthogroups (91 in total)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.line.x = element_line(colour = "grey40"),
        axis.ticks.x = element_line(colour = "grey40"))

ggsave("orthogroup_size_by_species.png", width = 8, height = 3, dpi = 300)


# Boxplot for MAFFFT results
# load data
cs <- read_csv("WRKY_conservation_summary(Conservation Summary).csv")

# define pairs
pairs <- cs %>%
  select(Orthogroup,
         `Rice vs maize`     = Rice_Maize_MeanIdentity,
         `Rice vs sorghum`   = Rice_Sorghum_MeanIdentity,
         `Maize vs sorghum`  = Maize_Sorghum_MeanIdentity) %>%
  pivot_longer(-Orthogroup, names_to = "Pair", values_to = "Identity") %>%
  drop_na(Identity) %>%                      # removes the NA orthogroups
  mutate(Pair = factor(Pair, levels = c("Rice vs maize", "Rice vs sorghum", "Maize vs sorghum")))

pairs %>%
  group_by(Pair) %>%
  summarise(n = n(), mean = mean(Identity), median = median(Identity),
            sd = sd(Identity), min = min(Identity), max = max(Identity))

# create plot
ggplot(pairs, aes(x = Pair, y = Identity, fill = Pair)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey25", alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.5, colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", colour = "black") +
  scale_fill_manual(values = c("Rice vs maize" = "#9ECAE1",
                               "Rice vs sorghum" = "#6BAED6",
                               "Maize vs sorghum" = "#2171B5"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = NULL, y = "Mean pairwise amino acid identity per orthogroup (%)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey40"),
        axis.ticks = element_line(colour = "grey40"))

ggsave("identity_by_species_pair.png", width = 6, height = 4.5, dpi = 300)


# Venn diagram of WRKY orthogroups shared among rice, sorghum and maize
library(ggVennDiagram)

og <- read_tsv("Orthogroups.tsv", show_col_types = FALSE)

has_genes <- function(col) og$Orthogroup[!is.na(col) & col != ""]

sets <- list(
  Rice    = has_genes(og$Osativa_indica_WRKY),
  Sorghum = has_genes(og$Sbicolor_WRKY),
  Maize   = has_genes(og$Zmays_WRKY)
)

stopifnot(
  lengths(sets) == c(Rice = 73, Sorghum = 83, Maize = 77),
  length(Reduce(intersect, sets)) == 59
)

cols <- c(Rice = "#0072B2", Sorghum = "#E69F00", Maize = "#009E73")

p <- ggVennDiagram(
  sets,
  label       = "count",
  label_alpha = 0,
  label_size  = 5,
  set_size    = 5.5,
  set_color   = cols,
  edge_size   = 1
) +
  scale_fill_gradient(low = "white", high = "grey80", guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  theme(plot.margin = margin(10, 10, 10, 10))

p

ggsave("orthogroup_venn.png", p, width = 6, height = 5, dpi = 300)
