# plots to check ages
# Uses dat3 object from reformat_covs_for_limix.R
library(tidyr)
library(ggplot2)
library(cowplot)
library(yarrr)
library(ggsignif)

pal <- piratepal(palette = "xmen") %>% unname

ageplot <- ggplot(filter(dat3, !is.na(phase) & !is.na(batch)), 
                  aes(y = ageDif, x = as.factor(phase), fill = as.factor(batch))) +   
  geom_jitter(pch = 21, width = 0.3, height = 0.01, size = 2) +
  scale_fill_manual(values = pal, name = "Batch") +
  labs(y = "Years after initial appointment", x = "INTERVAL 'phase' used for RNA")

ageBoxPlot <- ggplot(filter(dat3, !is.na(phase) & !is.na(batch)), 
                     aes(y = age_RNA, x = as.factor(phase), fill = as.factor(batch))) +
  geom_boxplot() +
  geom_point(pch = 20, size = 1, position = position_jitterdodge()) +
  scale_fill_manual(values = pal, name = "Batch") +
  labs(y = "Age", x = "INTERVAL 'phase' used for RNA") +
  geom_signif(comparisons = list(c("24m","48m"), c("24m","p3"), c("48m","p3")),
              map_signif_level = T,
              y_position = c(80,90,85))

ageDiffBoxPlot <- ggplot(filter(dat3, !is.na(phase) & !is.na(batch)), 
                         aes(y = ageDif, x = as.factor(phase), fill = as.factor(batch))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 20, size = 1, position = position_jitterdodge(jitter.height = 0.01)) +
  scale_fill_manual(values = pal, name = "Batch") +
  labs(y = "Years since initial appointment", x = "INTERVAL 'phase' used for RNA") +
  geom_signif(comparisons = list(c("24m","48m"), c("24m","p3"), c("48m","p3")),
              map_signif_level = T,
              y_position = c(4,4.4,4.2))


bothPlot <- plot_grid(ageBoxPlot, ageDiffBoxPlot)

ggsave(plot = ageplot, filename = "scripts/RNAseq/ageplot.png")
ggsave(plot = ageBoxPlot, filename = "scripts/RNAseq/ageboxplot.png")
ggsave(plot = bothPlot, filename = "scripts/RNAseq/bothplot.png")

datTall <- dat2 %>% 
  select(sample_id, agePulse, age1, age24m, age48m, agep3) %>%
  gather(age1:agep3, key = "phase", value = "age", -sample_id)

ggplot(datTall, aes(x = agePulse, y = age, colour = phase)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) 
