#=======================================================================
# Loading Packages
#=======================================================================
library(tidyverse)
library(TSclust)
library(ape)
#=======================================================================
# Loading data
#=======================================================================
load("fitted_all.RData")
#=======================================================================
# Plots
#=======================================================================
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

fitted_values$country <- as.factor(fitted_values$country)
levels(fitted_values$country)[c(25,29,38,40,54,58,67,145,162,
                                164,166,195,198,199,200,201,202)] <-
  c("British Virgin Isl","BSES","IC Japan","CAE","DRC",
    "Dominican Rep","Falkland Isl","NMI","SKN","SVG","STP","TC Islands",
    "UAE","UK","Tanzania","USA","US Virgin Isl")

#BSES = Bonaire, Saint Eustatius and Saba
#IC Japan = Cases on an international conveyance - Japan
#CAE = Central African Republic
#DRC = Democratic Republic of the Congo
#NMI = Northern Mariana Islands
#SKN = Saint Kitts and Nevis
#SVG = Saint Vincent and the Grenadines
#STP = Sao Tome and Principe
#TC Islands = Turks and Caicos Islands
#UAE = United Arab Emirates

last_60 <- fitted_values %>%
  filter(as.numeric(day) - 18262 > 79) %>%
  select(day, country, ar) %>%
  pivot_wider(id_cols = 1:2,
              names_from = "country",
              values_from = "ar") %>%
  select(-day)

tsdist <- diss(t(last_60), "DTWARP")
names(tsdist) <- colnames(last_60)
hc <- hclust(tsdist, "ward.D2")

cols <- gg_color_hue(10)
cols5 <- cols[5]
cols[5] <- cols[2]
cols[2] <- cols5
clus10 <- cutree(hc, 10)

country_query <- c("Brazil","USA")
country_fonts <- rep(1, nrow(tsdist))
country_fonts[names(last_60) %in% country_query] <- 2
country_cex <- rep(.7, nrow(tsdist))
country_cex[names(last_60) %in% country_query] <- 1

pdf("cluster_last60.pdf", w = 10, h = 10)
plot(as.phylo(hc), type = "fan", tip.color = cols[clus10],
     cex = country_cex,
     font = country_fonts,
     edge.width = .7)
dev.off()
#=======================================================================
