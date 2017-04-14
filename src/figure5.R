source('src/helpers.R')

###
### Panel 1 - Rare splicing variants
###

library(annotables)
ensembl2gene = data.frame(ensgene=grch37$ensgene, symbol=grch37$symbol)

variants = read_tsv('data/outliers/variant_data.tsv')

rare_splicing_events = subset(variants, n_outlier_lineages > 0 & splicing != 'not_splicing')

rare_splicing_data = plyr::ddply(rare_splicing_events, plyr::.(ensembl_gene_id), function(x) {
	dir = 'data/outliers/rare_qtl_models/'
	filename = paste0(x$ensembl_gene_id,'-','chr',x$chrom,'-',x$snp,'.csv')
	cbind(read.table(paste0(dir, filename)), snp=paste0(x$chrom,':',x$snp), hgvs=paste0('chr',x$chrom,':g.',x$snp,x$ref,'>',x$alt))
})

rare_splicing_data = merge(rare_splicing_data, ensembl2gene, by.x='ensembl_gene_id', by.y='ensgene')
rare_splicing_data$label = paste0(rare_splicing_data$symbol)

ggplot(rare_splicing_data, aes(x=label, y=expression, color=factor(genotype))) + 
	geom_hline(yintercept=0, color='black', lty=3) + 
	geom_point(position=position_jitterdodge(.2,dodge.width=.5)) + 
	theme_srd(fontsize=20) + 
	scale_color_manual(labels=c('Ref / Ref', 'Ref / Alt'), limits=c(0,1),values=c('gray','red')) +
	theme(legend.position=c(0,1), legend.direction='horizontal', legend.justification=c(0,1), legend.background=element_rect(fill='transparent',color=NA)) +
	labs(x='',y='Expression (z-score)') + 
	theme(axis.text.x=element_text(angle=0, hjust=.5, vjust=0.8))

ggsave('figs/fig5_panel1.png')
