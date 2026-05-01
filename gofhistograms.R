rm(list=ls())

library(ggplot2)
library(dplyr)
library(gridExtra)

gofall <- read.csv("gofall.csv")

custom_breaks <- seq(from = 0, to = 1, length.out = 21)

gofplots <- list()

bplpplot <- ggplot(gofall, aes(x = BPLP, fill = Taxa)) +
  geom_histogram(breaks = custom_breaks, alpha = 0.8) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(0, 1)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Bounded Power Law G.O.F. test p-values", y = "Number of datasets") +
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = c(0.85, 0.8))+
   labs(tag = "A")

lognormpplot <- ggplot(gofall, aes(x = lognormP, fill = Taxa)) +
  geom_histogram(breaks = custom_breaks,  alpha = 0.8) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(0, 1)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Log-normal G.O.F. test p-values", y = "Number of datasets") +
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = c(0.85, 0.8)) +
  labs(tag = "B")

# corals only for minus sampled data
msbplpplot <- ggplot(gofall, aes(x = MSBPLP)) +
  geom_histogram(breaks = custom_breaks, alpha = .4) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(0, 1)) +
  labs(x = "Corals: Minus-sampled Bounded Power Law G.O.F. test p-values", y = "Number of datasets") +
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(tag = "C")

mslnpplot <- ggplot(gofall, aes(x = MSlognormP)) +
  geom_histogram(breaks = custom_breaks, alpha = .4) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(0, 1)) +
  labs(x = "Corals: Minus-sampled Log-normal G.O.F. test p-values", y = "Number of datasets") +
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(tag = "D")

gofplots <- c(bplpplot,lognormpplot,msbplpplot,mslnpplot)

ggsave(
  filename = "gofplots.pdf", 
  plot = marrangeGrob(grobs= gofplots, nrow=2, ncol=2, list(top=NULL),
                      layout_matrix = rbind(c(1,2), c(3,4))), 
  width = 17, height = 15
)
