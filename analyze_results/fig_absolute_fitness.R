# Code for reproducing Fig. S21 in manuscript's appendix
# Fig is lifetime fitness as a function of \widehat{s} and z
# to examine effects of longevity on sensitivity to z

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

w.by.s.by.z = expand.grid(
  # Phenotype
  z = (1:300)/100,
  # \widehat{s}
  s.hat = (1:99)/100,
  # Fixing max. intrinsic fitness at 3
  W.max = 3
) %>%
  # W(z), lifetime fitness
  mutate(W = W.max * (1 - s.hat)*exp(-z^2/2) / (1 - s.hat*exp(-z^2/2)))

w.by.s.by.z %>%
  ggplot(aes(x = z, y = W, colour = s.hat)) +
  annotate('segment', x = 0, xend = 3, y = 1, yend = 1, linetype = 2, colour = 'gray88') +
  geom_line(aes(group = s.hat), linewidth = 1.2) +
  geom_point(
    data = w.by.s.by.z %>% filter(s.hat %in% c(0.1, 0.5, 0.9)),
    aes(fill = factor(s.hat, labels = c('low', 'medium', 'high'))),
    size = 3, shape = 21, stroke = 0.125
  ) +
  scale_colour_viridis_c(breaks = (1:4)/5, name = expression(hat(s))) +
  scale_fill_manual(
    values = c("#999999", "#56B4E9", "#E69F00"), 
    breaks = c('high', 'medium', 'low'),
    'simulation\nlongevity'
  ) +
  labs(x = 'z', y = 'W(z)') +
  theme(
    panel.background = element_rect(fill = NA),
    legend.position = 'top',
    legend.title = element_text(vjust = 0.75, hjust = 0.5),
    text = element_text(family = 'ArialMT')
  )

# Export figure 
ggsave('analyze_results/figs_out/figS1_lifetime_fitness.pdf', width = 8, height = 5)
  