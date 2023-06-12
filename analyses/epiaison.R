rm(list = ls())
setwd("~/Stage/analyses")
library(ggplot2)

getwd()

load("../donnees/bac")

ggplot(bac , aes(x = epiaison)) + geom_histogram()


ggplot(bac , aes(x = epiaison , fill = semis)) + geom_histogram()
