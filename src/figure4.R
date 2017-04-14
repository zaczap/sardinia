source('src/helpers.R')

###
### Panel 1 - TSS/TES enrichments
###

outliers = read_tsv('data/outliers/top.outliers.tsv')
variants = read_tsv('data/outliers/variant_data.tsv')

# How many shared rare variants in outliers vs. non-outlier lineages
variants %>%
	group_by(n_outlier_lineages) %>%
	summarize(N=n())
N.outlier = sum(variants$n_outlier_lineages)
N.nonoutlier = sum(variants$n_nonoutlier_lineages & variants$n_outlier_lineages == 0)

# Gene body +/- 250 KB
variants %>%
	filter(dist_to_tss >= -250000 & dist_to_tes <= 250000) %>%
	summarize(n_outlier = sum(n_outlier_lineages),
						n_nonoutlier = sum(n_nonoutlier_lineages > 0 & n_outlier_lineages == 0)) %>%
	mutate(relative_enrichment = n_outlier/(n_nonoutlier/121))

# Around TSS
variants %>%
	filter(abs(dist_to_tss) <= 250000) %>%
	summarize(n_outlier = sum(n_outlier_lineages),
						n_nonoutlier = sum(n_nonoutlier_lineages > 0 & n_outlier_lineages == 0)) %>%
	mutate(relative_enrichment = n_outlier/(n_nonoutlier/121))

# Around TES
variants %>%
	filter(abs(dist_to_tes) <= 250000) %>%
	summarize(n_outlier = sum(n_outlier_lineages),
						n_nonoutlier = sum(n_nonoutlier_lineages > 0 & n_outlier_lineages == 0)) %>%
	mutate(relative_enrichment = n_outlier/(n_nonoutlier/121))

# Calculate enrichment points in sliding window for TSS
point.data.tss = foreach(w = seq(-250000,250000,1000), .combine=rbind) %do% {
	x = subset(variants, w <= dist_to_tss & dist_to_tss <= w+1000)
	a = sum(x$n_outlier_lineages)
	b = nrow(x) - a
	data.frame(mid=(2*w + 1000)/2, y=a/((b+1)/121))
}

# Calculate enrichment points in sliding window for TES
point.data.tes = foreach(w = seq(-250000,250000,1000), .combine=rbind) %do% {
	x = subset(variants, w <= dist_to_tes & dist_to_tes <= w+1000)
	a = sum(x$n_outlier_lineages)
	b = nrow(x) - a
	data.frame(mid=(2*w + 1000)/2, y=a/((b+1)/121))
}

# Calculate expanding enrichment line for TSS
ov = filter(variants, n_outlier_lineages > 0)
nv = filter(variants, n_outlier_lineages == 0)
tss.data = foreach(w = seq(5000,250000,5000), .combine=rbind) %do% {
	a = nrow(ov[ov$dist_to_tss <= w & ov$dist_to_tss >= (-w),])
	b = nrow(nv[nv$dist_to_tss <= w & nv$dist_to_tss >= (-w),])
	enrich = a/(b/121)
	data.frame(n_rv=a, n_bg_rv=b, E=enrich, c_rv=a/(2*w), c_bg_rv=b/(2*w), type='TSS')
}

# Calculate expanding enrichment for TES
tes.data = foreach(w = seq(5000,250000,5000), .combine=rbind) %do% {
	a = nrow(ov[ov$dist_to_tes <= w & ov$dist_to_tes >= (-w),])
	b = nrow(nv[nv$dist_to_tes <= w & nv$dist_to_tes >= (-w),])
	enrich = a/(b/121)
	data.frame(n_rv=a, n_bg_rv=b, E=enrich, c_rv=a/(2*w), c_bg_rv=b/(2*w), type='TES')
}

# Calculate smoothed line for TSS and TES
line.data = plyr::ddply(rbind(tss.data, tes.data), plyr::.(type), function(x) {
	E = c(rev(x$E)[1:(nrow(x))], x$E[2:nrow(x)])
	x = as.vector(quantile(-250000:250000, probs=seq(0, 1, 1/(length(E)-1))))
	data.frame(bin=x, enrichment=E)
})

# Generate smoothed point data w/ moving average
smoothed.point.data.tss = with(point.data.tss, data.frame(mid = rollmean(mid, 4), y = rollmean(y, 4), type='TSS'))
smoothed.point.data.tes = with(point.data.tes, data.frame(mid = rollmean(mid, 4), y = rollmean(y, 4), type='TES'))
point.data = bind_rows(smoothed.point.data.tss, smoothed.point.data.tes)
	
# Adjust x-coordinates of points and lines
buffer = 550000
point.data[point.data$type == 'TES',]$mid = point.data[point.data$type == 'TES',]$mid + buffer
line.data[line.data$type == 'TES',]$bin = line.data[line.data$type == 'TES',]$bin + buffer

# Generate plot
ggplot(subset(line.data, type=='TSS'), aes(x=bin, y=enrichment)) + geom_hline(yintercept=1,color='gray',lty=3) +
	geom_point(aes(x=mid, y=y, color=y), data=point.data) + 
	scale_colour_gradient2(low="gray", high='#c0392b') + 
	geom_line(lwd=1, color='#c0392b') +
	geom_line(aes(x=bin, y=enrichment), lwd=1, color='#c0392b', data=subset(line.data, type=='TES')) +
	theme_srd(fontsize=20) + theme(legend.position=0) +
	scale_x_continuous(limits=c(-250000,250000+buffer), breaks=c(-200000,0,200000,buffer-200000,buffer,buffer+200000), labels=c('-200kb','TSS','200kb','-200kb','TES','200kb')) +
	labs(x='', y='Enrichment (relative risk)') + scale_y_continuous(limits=c(0,5))

ggsave('figs/fig4_panel1.png')

###
### Panel 2 - Chromatin enrichments
###

# Subset to a specific window of variants and collapse annotations
variants_in_region = filter(variants, abs(dist_to_tss) <= 50000)
enhancer.elements = c("6_EnhG", "7_Enh", "12_EnhBiv")
promoter.elements = c("1_TssA", "2_TssAFlnk", "10_TssBiv")
transcribed.elements = c("3_TxFlnk","4_Tx", "5_TxWk" , "11_BivFlnk"  )
null.elements = c("9_Het","15_Quies","8_ZNF/Rpts" )
repressive.elements = c("14_ReprPCWk","13_ReprPC" )
mapping = data.frame(chromatin=c(enhancer.elements,promoter.elements,transcribed.elements,null.elements,repressive.elements),collapsed=rep(c('enhancer','promoter','transcribed','other','repressive'),c(3,3,4,3,2)))
variants_in_region = merge(variants_in_region, mapping, by='chromatin') %>% tbl_df()

# Perform Fisher's exact test on chromatin states
library(purrr)
enrichmentTest = function(x) {
	n.outlier = sum(x$n_outlier_lineages)
	n.nonoutlier = nrow(x) - n.outlier
	t = matrix(c(n.outlier, n.nonoutlier,N.outlier-n.outlier,N.nonoutlier-n.nonoutlier), nrow=2)
	tidy(fisher.test(t))
}

stateData = variants_in_region %>% 
	group_by(collapsed) %>% 
	nest() %>% 
	mutate(model = map(data, enrichmentTest)) %>%
	unnest(model) %>%
	select(-data)
	
# Perform Fisher's exact test on splicing annotation
x = filter(variants_in_region, splicing == 'splicing')
n.outlier = sum(x$n_outlier_lineages)
n.nonoutlier = nrow(x) - n.outlier
t = matrix(c(n.outlier, n.nonoutlier,N.outlier-n.outlier,N.nonoutlier-n.nonoutlier), nrow=2)
spliceData = tidy(fisher.test(t))
spliceData = cbind(collapsed='splicing', spliceData)

# Combine chromatin and splicing data
enrichment_data = bind_rows(stateData, spliceData)
enrichment_data$estimate = log(enrichment_data$estimate)
enrichment_data$conf.low = log(enrichment_data$conf.low)
enrichment_data$conf.high = log(enrichment_data$conf.high)

ggplot(enrichment_data, aes(x=collapsed, y=estimate)) + 
	geom_hline(yintercept=0, color='black', lty=3) + 
	geom_errorbar(aes(ymin=conf.low,ymax=conf.high), width=0, lwd=1, color='#34495e') +
	geom_point(size=4, color='#34495e') + 
	theme_srd() +
	labs(x='', y='ln (odds ratio)') +
	scale_x_discrete(limits=c('splicing','promoter','enhancer','transcribed','other','repressive')) +
	theme(legend.position=0) + 
	coord_flip() +
	scale_y_continuous(limits=c(-1,6), breaks=-1:6)

ggsave('figs/fig4_panel2.png')

###
### Panel 3 - Rare eQTL replication
###

population_rare_eqtl = read_tsv('data/outliers/rare_eqtl_results.tsv')
variants$merge.key = paste(variants$ensembl_gene_id, variants$snp, sep='-')
plot.data = merge(population_rare_eqtl, variants, by='merge.key')

top = plyr::ddply(subset(plot.data, abs(dist_to_tss) < 250000), plyr::.(ensembl_gene_id), function(x) x[which.max(x$CADD),])

ggplot(top, aes(x=dist_to_tss, y=beta, size=-log10(bh), color=beta > 0)) + 
	geom_hline(yintercept=0) + 
	geom_smooth(aes(fill=beta > 0)) + 
	labs(x='', y=expression('eQTL Effect '~ (beta))) + 
	geom_point(alpha=.5) + 
	theme_srd() + 
	theme(legend.position='bottom') + 
	scale_color_manual(values=c(color.underexpression,color.overexpression)) + 
	scale_fill_manual(values=c(color.underexpression,color.overexpression),label=c('Under-expression','Over-expression')) +
	scale_y_continuous(limits=c(-6.5,6.5)) +
	scale_x_continuous(breaks=c(-200000,-100000,0,100000,200000),labels=c('-200KB','-100KB','TSS','+100KB','+200KB')) +
	guides(size=FALSE,color=FALSE,fill = guide_legend(override.aes = list(linetype = 0))) +
	theme(axis.title.x=element_blank())

ggsave('figs/fig4_panel3.png')

###
### Panel 4 - Functional constraint metrics
###

library(plyr)

plot.data = merge(population_rare_eqtl, variants, by='merge.key')

plot.data$gerp_bin = label_by_bin(plot.data$gerp, probs=c(0,.25,.50,.75,.90,.99,1))
plot.data$CADD_bin = label_by_bin(plot.data$CADD, probs=c(0,.25,.50,.75,.90,.99,1))
plot.data$phylop_bin = label_by_bin(plot.data$phylop, probs=c(0,.25,.50,.75,.90,.99,1))
plot.data$fitCons_bin = label_by_bin(plot.data$fitCons, probs=c(0,.25,.50,.75,.90,.99,1))

gerp_data = ddply(plot.data, .(gerp_bin), function(x) {
	m = mean(-log10(x$pvalue))
	s = sd(-log10(x$pvalue))/sqrt(nrow(x))
	l = m - s
	u = m + s
	data.frame(mean=m, conf.low=l, conf.high=u, type='GERP', n=nrow(x), s=s)
})

CADD_data = ddply(plot.data, .(CADD_bin), function(x) {
	m = mean(-log10(x$pvalue))
	s = sd(-log10(x$pvalue))/sqrt(nrow(x))
	l = m - s
	u = m + s
	data.frame(mean=m, conf.low=l, conf.high=u, type='CADD', n=nrow(x), s=s)
})

phylop_data = ddply(plot.data, .(phylop_bin), function(x) {
	m = mean(-log10(x$pvalue))
	s = sd(-log10(x$pvalue))/sqrt(nrow(x))
	l = m - s
	u = m + s
	data.frame(mean=m, conf.low=l, conf.high=u, type='PhyloP', n=nrow(x), s=s)
})

fitCons_data = ddply(plot.data, .(fitCons_bin), function(x) {
	m = mean(-log10(x$pvalue))
	s = sd(-log10(x$pvalue))/sqrt(nrow(x))
	l = m - s
	u = m + s
	data.frame(mean=m, conf.low=l, conf.high=u, type='FitCons', n=nrow(x), s=s)
})

plot.data = data.table::rbindlist(list(gerp_data, CADD_data, phylop_data, fitCons_data))
colnames(plot.data) <-  c('bin','mean','conf.low','conf.high','type','N','sd')

ggplot(plot.data, aes(x=bin, y=mean, group=type, color=type)) + 
	geom_errorbar(aes(ymin=conf.low, ymax=conf.high, fill=type), alpha=.3, position=position_dodge(.5), width=0, lwd=2) +
	geom_point(size=3, position=position_dodge(.5)) +
	theme_srd(fontsize=20) + 
	theme(legend.position=c(0,1), legend.justification=c(0,1)) + 
	labs(x='Quantile (%)', y='-log10 eQTL p-value') +
	scale_x_discrete(labels=c('0-25','25-50','50-75','75-90','90-99','99-100')) +
	scale_color_manual(name='',breaks=c('GERP','CADD','PhyloP','FitCons'),values=c('#3498db','#e67e22','#2ecc71','#8e44ad')) +
	scale_fill_manual(name='',breaks=c('GERP','CADD','PhyloP','FitCons'),values=c('#3498db','#e67e22','#2ecc71','#8e44ad'))

ggsave('figs/fig4_panel4.png')
