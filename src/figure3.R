source('src/helpers.R')

###
### Panel 1 - Shared expression outlier example
###

outliers = read_tsv('data/outliers/top.outliers.tsv')

# Load the population expression data for the example gene: ENSG00000187994
target = read_tsv('data/outliers/rare_qtl_models/ENSG00000187994-chr19-39368871.csv')
target.outlier = filter(outliers, simple.id == 'ENSG00000187994')

# Draw the histogram
data = data.frame(z=target$expression)
histogram = ggplot(data, aes(x=z)) + 
	geom_histogram(color='black', fill='#ecf0f1', bins=75) + 
	theme_srd() +
	labs(x='',y='') + 
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.border=element_blank(), axis.line.y = element_blank())

# Annotate the histogram
xl = .025
ggdraw() + 
	draw_plot(histogram , .2, -.1, .8, .9) + 
	draw_grob(linesGrob(x=c(.32,.35),y=c(.55,.3), arrow=arrow(ends="last", length=unit(10,'pt')), gp=gpar(lwd=3, col='#c0392b'))) +
	draw_grob(linesGrob(x=c(.39,.37),y=c(.43,.31), arrow=arrow(ends="last", length=unit(10,'pt')), gp=gpar(lwd=3, col='#c0392b'))) +
	draw_grob(linesGrob(x=c(.7975,.77),y=c(.72,.54), arrow=arrow(ends="last", length=unit(10,'pt')), gp=gpar(lwd=3, col='#7f8c8d'))) +
	draw_label('Gene expression (z-score)',x=.6,y=.925, hjust=0.5, size=20) +
	draw_label('Outlier parent',x=.32,y=.58, hjust=0.5, vjust=0, size=16, colour = '#c0392b') +
	draw_label('Outlier child',x=.39,y=.46, hjust=0.5, vjust=0, size=16, colour = '#c0392b') +
	draw_label('Non-outlier parent',x=.7975,y=.75, hjust=0.5, vjust=0, size=16) +
	draw_grob(linesGrob(x=c(0.075+xl,0.15+xl),y=c(.75,.75), gp=gpar(lwd=2, color='#34495e'))) +
	draw_grob(linesGrob(x=c(0.1125+xl,0.1125+xl),y=c(.75,.45), gp=gpar(lwd=2, color='#34495e'))) +
	draw_grob(rectGrob(width=unit(19,'pt'),height=unit(19,'pt'),x=0.075+xl,y=0.75, hjust=0.5, vjust=0.5, gp=gpar(color='black',fill='#c0392b', lwd=2))) +
	draw_grob(circleGrob(r=unit(10,'pt'), x=.15+xl,y=.75, gp=gpar(fill='#ecf0f1', color='', lwd=2))) +
	draw_grob(circleGrob(r=unit(10,'pt'), x=.1125+xl,y=.45, gp=gpar(fill='#c0392b', color='', lwd=2))) +
	draw_label(round(target.outlier$parent.residual,2),x=0.075+xl,y=.615, col='#c0392b', hjust=.6, size=16)+
	draw_label(round(target.outlier$nonparent.residual,2),x=0.15+xl,y=.615, col='black', hjust=.5, size=16)+
	draw_label(round(target.outlier$child.residual,2),x=0.1125+xl,y=.31, col='#c0392b', hjust=.6, size=16) +
	draw_label("Discovery Trio", x=.1125+xl,y=.2,hjust=0.5, size=16) 

ggsave('figs/fig3_panel1.png', width=12,height=2)

###
### Panel 2 - Outlier scatterplot
###

ggplot(subset(outliers, fdr.group == 3), aes(x=child.residual, y=parent.residual)) + 
	geom_point(color='#bdc3c7', alpha=.5) + 
	theme_srd(fontsize=20) + 
	geom_point(aes(x=child.residual, y=parent.residual, color=factor(fdr.group)), data=subset(outliers, fdr.group != 3), alpha=.5, size=3) +
	geom_hline(yintercept=0, color='gray', lty=3) + 
	geom_vline(xintercept=0, color='gray', lty=3) +
	labs(x='Child expression (z-score)', y='Parent expression (z-score)') +
	theme(legend.position = c(0,1), legend.justification = c(0,1)) +
	scale_color_manual(name='',values=c('#c0392b','#FFB31B'), labels=c('Outlier 5% FDR','Outlier 10% FDR','Top Outliers')) +
	scale_x_continuous(limits=c(-10,10)) + 
	scale_y_continuous(limits=c(-10,10))

ggsave('figs/fig3_panel2.png')

###
### Panel 3 - Allelic imbalance
###

# Load ASE data
ase = read_tsv('data/outliers/ase_data.tsv')
significant_outliers = filter(outliers, fdr.group < 3)
outlier.keys = c(paste(significant_outliers$child.id, significant_outliers$gene.id, sep=':'),
								 paste(significant_outliers$parent.id, significant_outliers$gene.id, sep=':'))
up_regulated_genes = subset(significant_outliers, child.residual > 0)$gene.id

# Separate ASE data sets by outliers and non-outliers
ase$is.outlier = ifelse(ase$set.key %in% outlier.keys, 'Outliers', 'Non-outliers')
ase.outliers = subset(ase, set.key %in% outlier.keys)
ase.nonoutliers = subset(ase, gene.id %in% significant_outliers$gene.id & !(set.key %in% outlier.keys))

# Re-arrange data and plot
plot.data = bind_rows(ase.outliers, ase.nonoutliers)
plot.data$effect.direction = ifelse(plot.data$gene.id %in% up_regulated_genes, 'Over-expression', 'Under-expression')

plot.data$xscale = 0
plot.data[plot.data$effect.direction == 'Over-expression' & plot.data$is.outlier == 'Non-outliers',]$xscale = 1
plot.data[plot.data$effect.direction == 'Over-expression' & plot.data$is.outlier == 'Outliers',]$xscale = 2
plot.data[plot.data$effect.direction == 'Under-expression' & plot.data$is.outlier == 'Non-outliers',]$xscale = 3
plot.data[plot.data$effect.direction == 'Under-expression' & plot.data$is.outlier == 'Outliers',]$xscale = 4

labels = paste(c('Non-outliers', 'Outliers'), '\nn = ',c(sum(plot.data$xscale == 1),sum(plot.data$xscale == 2),sum(plot.data$xscale == 3),sum(plot.data$xscale == 4)), sep='')

plot.data %>%
	ggplot(aes(x=as.factor(xscale), fill=is.outlier, y=abs(.5-ref.count/(ref.count+alt.count)))) + 
	geom_boxplot(position=position_dodge(.75), outlier.shape=NA, width=.9, lwd=.8) + 
	theme_srd()  + 
	labs(x='', y='Allelic imbalance (AI)') + 
	geom_vline(xintercept=2.5, color='black', lty=3) + 
	scale_x_discrete(limits=1:4, labels=labels) +
	annotate('text',x=1.5,y=.475,label='Over-expression',hjust=.5,size=5) +
	annotate('text',x=3.5,y=.475,label='Under-expression',hjust=.5,size=5) + 
	scale_fill_manual(limits=c('Non-outliers','Outliers'), values=c(color.nonoutliers,color.outliers)) +
	theme(legend.position = 0)

ggsave('figs/fig3_panel3.png')

###
### Panel 4 - Intra-family correlations
###

# Calculate expression correlation
expression_correlation_matrix = cor(filter(outliers, fdr.group < 3) %>% select(nonparent.residual, parent.residual, child.residual))
colnames(expression_correlation_matrix) = c('Non-outlier\nParent', 'Outlier\nParent', 'Outlier\nChild')
rownames(expression_correlation_matrix) = c('Non-outlier\nParent', 'Outlier\nParent', 'Outlier\nChild')
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
expression_correlation_plot = reshape2::melt(expression_correlation_matrix) %>%
	ggplot(aes(x=Var1, y=Var2, fill=value)) + 
	geom_tile() + 
	theme_srd() + labs(x='', y='') + 
	scale_fill_gradientn(guide=F, name='', colours=c("#c0392b", "#e74c3c", "#FFFFFF", "#3498db", "#2980b9"), limits=c(-1,1)) + 
	geom_text(aes(x=Var1, y=Var2, label=round(value,2)), size=6, color='white') +
	theme(panel.border = element_blank(), axis.text.x= element_text(angle=45, hjust=1, color='white'), axis.ticks.x=element_line(color='white')) +
	scale_x_discrete(limits=c('Outlier\nChild','Outlier\nParent','Non-outlier\nParent')) +
	scale_y_discrete(limits=rev(c('Outlier\nChild','Outlier\nParent','Non-outlier\nParent'))) +
	ggtitle('Expression') + 
	theme(axis.line.x = element_blank(), axis.line.y = element_blank())

# Calculate allele-specific expression correlation
parents = unique(outliers$parent.id) # 122
children = unique(outliers$child.id) # 61
require(plyr)
outlier_gene_ase_data = subset(ase, gene.id %in% significant_outliers$gene.id)
melted_ase = ddply(significant_outliers, .(gene.id), function(x) {
	child = x$child.id[1]
	parent = x$parent.id[1]
	other_parent = ifelse(x$parent.sex[1] == 'mother', x$father.id[1], x$mother.id[1])
	
	ase_data = subset(outlier_gene_ase_data, gene.id == x$gene.id[1] & sample_id %in% c(child, parent, other_parent))
	
	grouped_ase = ddply(ase_data, .(stop), function(y) {
		if(nrow(y) == 3) {
			child_data = subset(y, sample_id %in% children & is.outlier == 'Outliers')
			oparent_data = subset(y, sample_id %in% parents & is.outlier == 'Outliers')
			nparent_data = subset(y, sample_id %in% parents & is.outlier != 'Outliers')
			child_ase = abs(.5-child_data$ref.count/(child_data$ref.count+child_data$alt.count))
			parent_ase = abs(.5-oparent_data$ref.count/(oparent_data$ref.count+oparent_data$alt.count))
			nonparent_ase = abs(.5-nparent_data$ref.count/(nparent_data$ref.count+nparent_data$alt.count))
			child_ase = child_data$ref.count/(child_data$depth)
			parent_ase = oparent_data$ref.count/(oparent_data$depth)
			nonparent_ase = nparent_data$ref.count/(nparent_data$depth)
			data.frame(c=child_ase,p=parent_ase,n=nonparent_ase)
		}
	})
	grouped_ase
})

ase_correlation_matrix = cor(as.matrix(melted_ase[,c(15,13,14)]))
colnames(ase_correlation_matrix) = c('Non-outlier\nParent', 'Outlier\nChild', 'Outlier\nParent')
rownames(ase_correlation_matrix) =c('Non-outlier\nParent', 'Outlier\nChild', 'Outlier\nParent')
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
ase_correlation_plot = ggplot(reshape2::melt(ase_correlation_matrix), aes(x=Var1, y=Var2, fill=value)) +  
	geom_tile() + 
	theme_srd() + labs(x='', y='') + 
	scale_fill_gradientn(guide=F, name='', colours=c("#c0392b", "#e74c3c", "#FFFFFF", "#3498db", "#2980b9"), limits=c(-1,1)) + 
	geom_text(aes(x=Var1, y=Var2, label=round(value,2)), size=6, color='white') +
	theme(panel.border = element_blank(), axis.text.x= element_text(angle=45, hjust=1, color='white'), axis.ticks.x=element_line(color='white')) +
	scale_x_discrete(limits=c('Outlier\nChild','Outlier\nParent','Non-outlier\nParent')) +
	scale_y_discrete(limits=rev(c('Outlier\nChild','Outlier\nParent','Non-outlier\nParent'))) +
	ggtitle('Allelic imbalance') + 
	theme(axis.line.x = element_blank(), axis.line.y = element_blank())

ggsave('figs/fig3_panel4_expression.png', expression_correlation_plot)
ggsave('figs/fig3_panel4_ase.png', ase_correlation_plot)

###
### Panel 5 - Relationship between allelic imbalance and expression
###

# Calculate mean/std. error of allelic imbalance in bins
library(purrr)
extract_ai = function(data) { abs(.5 - data$ref.count / (data$ref.count + data$alt.count))}
merged.plot.boxes = plot.data %>%
	left_join(significant_outliers, by='gene.id') %>%
	mutate(effect.bin = as.numeric(cut(abs(average.effect), c(1,2,3,4,5,6,7,8,Inf)))) %>%
	group_by(effect.bin, is.outlier) %>%
	nest() %>%
	mutate(ai = map(data, extract_ai)) %>%
	unnest() %>%
	group_by(effect.bin, is.outlier) %>%
	dplyr::summarize(mean = mean(ai),
						n = length(ai),
						sd = sd(ai),
						std = sd/sqrt(length(ai))) %>%
	filter(n > 10)

# Generate figure	
dodge = position_dodge(.5)
ggplot(merged.plot.boxes, aes(x=effect.bin, color=is.outlier, y=mean)) + 
	geom_point(position=dodge, size=3) + 
	geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position=dodge, width=.3, alpha=.3, lwd=1) +
	scale_x_continuous(limits=c(1.5,7.5)) +
	theme_srd() + 
	scale_y_continuous(limits=c(0,.5)) + 
	labs(y='Allelic imbalance (AI)', x='Average expression (z-score)') +
	theme(legend.position=c(0,1),legend.justification=c(0,1), legend.background=element_rect(fill='transparent',color=NA)) +
	scale_color_manual(limits=c('Non-outliers','Outliers'), values=c(color.nonoutliers,color.outliers))

ggsave('figs/fig3_panel5.png')

# Test correlation between ASE and expression
outlier_ase_data = plot.data %>%
	left_join(significant_outliers, by='gene.id') %>%
	filter(is.outlier == 'Outliers') %>%
	mutate(ai = abs(.5 - (ref.count / ref.count + alt.count)))
cor.test(abs(outlier_ase_data$ref.count/outlier_ase_data$depth - .5),abs(outlier_ase_data$average.effect),method='spearman', exact=F) # rho = .338
