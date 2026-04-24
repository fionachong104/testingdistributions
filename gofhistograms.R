rm(list=ls())

library(ggplot2)
library(dplyr)

gofall <- read.csv("gofall.csv")

custom_breaks <- seq(from = 0, to = 1, length.out = 21)


ggplot(gofall, aes(x = lognormP, fill = Taxa)) +
  geom_histogram(breaks = custom_breaks,  alpha = 0.8) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(-0.2, 1.25)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Log-normal GOF test P-value", y = "Number of observations") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 ggplot(gofall, aes(x = BPLP, fill = Taxa)) +
  geom_histogram(breaks = custom_breaks, alpha = 0.8) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(-0.2, 1.25)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Bounded Power Law GOF test P-value", y = "Number of observations") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# corals only for minus sampled data
ggplot(gofall, aes(x = MSBPLP)) +
  geom_histogram(breaks = custom_breaks, alpha = .4) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(-0.2, 1.25)) +
  labs(x = "Minus-sampled Bounded Power Law GOF test P-value", y = "Number of observations") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(gofall, aes(x = MSlognormP)) +
  geom_histogram(breaks = custom_breaks, alpha = .4) +
  stat_bin(
    aes(label = after_stat(if_else(count > 0, as.character(count), "")), group = 1),
    geom = "text",
    breaks = custom_breaks,
    vjust = -1
  ) +
  scale_x_continuous(breaks = custom_breaks, limits = c(-0.2, 1.25)) +
  labs(x = "Minus-sampled Log-normal GOF test P-value", y = "Number of observations") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
