###
### Load required packages
###

library(readr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(foreach)
library(zoo)
library(broom)
library(RColorBrewer)
library(grid)
library(dplyr)

###
### Define color palettes
###

color.srd = '#E74C3C'
color.geu = '#3498DB'
color.dgn = '#BDC3C7'

color.overexpression = "#DCBD23"
color.underexpression = "#0365C0"

color.outliers = '#C0392B'
color.nonoutliers = '#BDC3C7'

###
### Define ggplot2 theme for paper
###

theme_srd = function(fontsize=14) {
	base_size = fontsize
	border_color = 'black'
	tick_color = 'black'
	linewidth = .5
	theme_bw(base_size=14) +
		theme(panel.background = element_blank(), 
					panel.border = element_blank(),
					panel.grid.major = element_blank(), #line(colour='#ecf0f1'),
					panel.grid.minor = element_blank(),
					axis.text.x = element_text(vjust = 1, size=base_size*.8, color='black'), 
					axis.text.y = element_text(hjust = 1, size=base_size*.8, color='black'), 
					axis.text = element_text(size=base_size),
					legend.title = element_blank(), 
					legend.key = element_blank(),
					legend.text=element_text(size=base_size*.8),
					legend.background=element_rect(fill='NA'),
					plot.title=element_text(size=base_size*.9,hjust=.5),
					axis.line = element_line())
}

###
### Shared functions
###

# Returns a numeric vector of labels for x after binning by quantile (default: by decile)
label_by_bin = function(x, probs=seq(0,1,.1)) {
	cut(x, breaks=quantile(x, probs=probs), include.lowest=T, labels=1:(length(probs)-1))
}