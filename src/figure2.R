source('src/helpers.R')

###
### Panel 1 - Allele frequency comparison
###

deltaAF = read_tsv('data/deltaAF/deltaAF_data.tsv', col_types = 'ccccdicdddicidd') %>%
	rename(DeltaAF=DeltaAf,EUR_AF=EURUF503_AF,SRD_AF=SRDUF691_AF)
deltaAF$bin = as.numeric(cut(abs(deltaAF$AbsDeltaAf), c(-Inf,.05,.1,.15,.2, Inf)))
deltaAF$bin = ifelse(deltaAF$IsDeltaAfOutlier == 1, 6, deltaAF$bin)
deltaAF$bin = as.factor(deltaAF$bin)

reds = brewer.pal(5, "Reds")
cutoffs = c('0.00','0.05','0.10','0.15','0.20')
color.labels = paste0(c(rep('|ΔAF| > ',5),'Top 1%'), c(cutoffs, ''))

deltaAF %>% 
	arrange(bin) %>%
	ggplot(aes(x=EUR_AF,SRD_AF,color=bin)) +
		geom_point() +
		geom_hline(yintercept=.5, lty=2, color='gray') +
		geom_vline(xintercept=.5, lty=2, color='gray') +
	theme_srd() +
	scale_color_manual(limits=c(1,2,3,4,5,6),values=c(reds,'#3498db'),labels=c(color.labels,'Top 1%')) +
	labs(x='European AF',y='Sardinian AF') +
	theme(legend.position=c(1,0),legend.justification = c(1,0))

ggsave('figs/fig2_panel1.png')
 
###
### Panel 2 - LD decay comparison
###

ld = read.table("data/LD/ld_decay.16Jan2015.unique.txt", header=T)
ldnull = read.table("data/LD/ld_decay_null_snps.txt", header=T)
p1 = read.table("data/LD/MerlinEqtl.PopGen.22Columns.Rfriendly.tsv", header=T)
q1 = merge(p1,ld,by.x="SnpIdQuery",by.y="SNP")
p2 = read.table("data/LD/eQTLs.Fdr0.05.WithConditional.supercededSoUnique.tsv", header=T)
q = merge(p2,q1,by.x="SNP",by.y="SnpIdQuery")

frequency = c(0.00,0.05,0.10,0.15,0.20)
distance =  seq(0,297000,by=3000)
rcol = brewer.pal(length(frequency), "Reds")
bcol = brewer.pal(length(frequency), "Blues")

rsq.data = plyr::ldply(frequency,  function(x) {
	maf_threshold=x
	qs = q[(q$SRDUF691_AF-q$EURUF503_AF) > x & !is.na(q$EURUF503_AF) & q$EURUF503_AF != 0 & !is.na(q$SRDUF691_AF), ]
	x = seq(0,297000,by=3000)
	y = as.vector(smooth(as.numeric(lapply(qs[,32:ncol(qs)], function(x) mean(x, na.rm=T)))))
	data.frame(distance=x, mean.r2 = y, AF=maf_threshold, grouping = paste('delta',maf_threshold), label = paste0("> ", format(round(maf_threshold, 2), nsmall = 2), " (n = ", nrow(qs),")") , stringsAsFactors = F)
})

rsq2.data = plyr::ldply(frequency,  function(x) {
	maf_threshold=x
	qs = q[q$SRDUF691_AF > x & !is.na(q$EURUF503_AF) & q$EURUF503_AF != 0 & abs(q$SRDUF691_AF-q$EURUF503_AF) < 0.025 & !is.na(q$SRDUF691_AF),]
	x = seq(0,297000,by=3000)
	y = as.vector(smooth(as.numeric(lapply(qs[,32:ncol(qs)], function(x) mean(x, na.rm=T)))))
	data.frame(distance=x, mean.r2 = y, AF=maf_threshold, grouping = paste('srd',maf_threshold), label = paste0(" > ", format(round(maf_threshold, 2), nsmall = 2), " (n = ", nrow(qs),")") , stringsAsFactors = F)
})

null.data = data.frame(distance=seq(0,297000,by=3000), mean.r2=as.vector(smooth(as.numeric(lapply(ldnull[,2:ncol(ldnull)], function(x) mean(x, na.rm=T))))), AF=0, grouping='null', label='> 0.2% (n = 417)', stringsAsFactors = F)

color.labels = unique(rsq.data$label)
color.limits = paste('delta', frequency)
color.values = c(bcol)

color.labels = c(paste0('|ΔAF| ', color.labels), ' ', paste0('AF ', unique(rsq2.data$label)), '  ','Non-eQTL |ΔAF| > 0.20','(n = 417)')
color.limits = c(color.limits, 'b', paste('srd', frequency), 'c', 'null','d')
color.values = c(rcol, 'white', bcol, 'white', '#000000','white')

sardinia_scale = scale_color_manual(limits=color.limits, values=color.values, labels=color.labels)

plot.data = rbind(rsq.data,rsq2.data)
plot.data %>%
	ggplot(aes(x=distance/1000,  y=mean.r2, group=grouping, color=grouping)) + 
	geom_line(lwd=1) + 
	theme_srd() + 
	labs(x='Distance (in kb)', y=expression('Mean r'^2~' value')) + 
	scale_y_continuous(limits = c(0,.2)) + 
	theme(legend.position=c(1,1), legend.justification=c(1,1)) + 
	scale_x_continuous(breaks=seq(0,300,50)) +
	sardinia_scale + 
	geom_line(data=null.data, aes(x=distance/1000,  y=mean.r2, group=grouping, color=grouping), lwd=1) + 
	guides(color= guide_legend(ncol=1))

ggsave('figs/fig2_panel2.png')

###
### Panel 3 - MS and malaria resistance
###

delta = deltaAF$DeltaAF
ms_sub = subset(deltaAF, IsMsGwas == 1)$DeltaAF
malaria_sub = subset(deltaAF, IsMalariaGene == 1)$DeltaAF
all_sub = delta

ms_percents = 100*table(as.numeric(cut(abs(ms_sub), c(-Inf,.05,.1,.15,.2, Inf))))/length(ms_sub)
malaria_percents = 100*table(as.numeric(cut(abs(malaria_sub), c(-Inf,.05,.1,.15,.2, Inf))))/length(malaria_sub)
all_percents = 100*table(as.numeric(cut(abs(delta), c(-Inf,.05,.1,.15,.2, Inf))))/length(delta)

reds = brewer.pal(5, "Reds")
cutoffs = c('0.00','0.05','0.10','0.15','0.20')
color.labels = paste0('|ΔAF| > ', cutoffs)

plot.data = data.frame(group=rep(c('MS','Malaria','All'),c(5,5,5)), percents=c(rev(ms_percents),rev(malaria_percents),rev(all_percents)), color=rev(c(1:5,1:5,1:5)))

n_vector = c(length(ms_sub), length(malaria_sub), length(all_sub))
ggplot(plot.data, aes(x=group,y=percents,fill=factor(color))) + 
	geom_bar(stat='identity',position='stack') + 
	theme_srd() +
	scale_x_discrete(limits=c('MS','Malaria','All'),labels=c('MS\neQTLs','Malaria\neQTLs','All\neQTLs')) +
	scale_fill_manual(limits=1:5,labels=color.labels, values=reds) +
	labs(x='', y='Percentage of hits')

ggsave('figs/fig2_panel3.png')

###
### Panel 4 - deltaAF/Fst enrichments
###

# Path to deltaAF hits
daf_gwas = "data/disease_enrichment/DeltaAf.GwasEnrichment.FisherEnrichment.SortedPvalue.tsv"
daf_immunobase = "data/disease_enrichment/DeltaAf.ImmunoBaseEnrichment.FisherEnrichment.SortedPvalue.tsv"    
daf_malaria = "data/disease_enrichment/DeltaAf.PhenopediaEnrichment.ByDisease_GenesCountsGe0.FisherEnrichment.Malaria.tsv" 

# Path to Fst hits
fst_gwas = "data/disease_enrichment/FstWBackgroundSnps.GwasEnrichment.FisherEnrichment.SortedPvalue.tsv"
fst_immunobase = "data/disease_enrichment/FstWBackgroundSnps.ImmunoBaseEnrichment.FisherEnrichment.SortedPvalue.tsv"
fst_malaria = "data/disease_enrichment/FstWBackgroundSnps.PhenopediaEnrichment.ByDisease_GenesCountsGe0.FisherEnrichment.Malaria.tsv"

# Specify header format
header = c("trait", "R1C1", "R2C1", "R1C2", "R2C2", "odds_ratio", "pvalue")

# Read deltaAF data
daf_gwas = read_tsv(daf_gwas, col_names=header) %>% mutate(signal="deltaAF", db="GWAS")
daf_immunobase = read_tsv(daf_immunobase, col_names=header) %>% mutate(signal="deltaAF", db="ImmunoBase")
daf_malaria = read_tsv(daf_malaria, col_names=header, skip=1) %>% mutate(signal="deltaAF", db="Phenopedia", trait="malaria")

# Read Fst data
fst_gwas = read_tsv(fst_gwas, col_names=header, skip=1) %>% mutate(signal="Fst", db="GWAS")
fst_immunobase = read_tsv(fst_immunobase, col_names=header,skip=1) %>% mutate(signal="Fst", db="ImmunoBase")
fst_malaria = read_tsv(fst_malaria, col_names=header, skip=1) %>% mutate(signal="Fst", db="Phenopedia", trait="malaria")

# Combine data into a single data frame for deltaAF / Fst, and combined
daf_data = bind_rows(daf_gwas, daf_immunobase, daf_malaria)
fst_data = bind_rows(fst_gwas, fst_immunobase, fst_malaria)
disease_data = bind_rows(daf_data, fst_data)

# Add the number of 'tests' performed:
test_data = disease_data %>% 
	group_by(db, signal) %>% 
	summarize(n_tests = n()) %>% 
	right_join(disease_data) %>% 
	ungroup 
test_data[test_data$db == 'GWAS',]$n_tests = 354 # some GWAS traits with missing records

# Perform Bonferroni correction
test_data = test_data %>%
	mutate(bonf_pvalue = pmin(n_tests * pvalue, 1), log10pvalue = -log10(bonf_pvalue)) %>%
	arrange(desc(log10pvalue)) %>%
	mutate(trait=as.factor(trait), db=factor(db), idx=1:nrow(.))

# These are the tests w/ enough data to consider:
supplemental_table_data = test_data %>% filter(R1C1 + R2C1 >= 10)

# Select the data we want to plot
plot_data = filter(supplemental_table_data, trait %in% unique(supplemental_table_data$trait)[1:10])
plot_data$trait = factor(plot_data$trait, levels=rev(unique(plot_data$trait)))

trait_labels = rev(c('Multiple\nsclerosis (MS)',
										 'Malaria', 
										 'Height',
										 'Type 2 diabetes',
										 'Rheumatoid arthritis',
										 'HDL cholesterol', 
										 'Mean platelet\nvolume', 
										 'Bipolar disorder', 
										 'Blood pressure', 
										 'Bone mineral\ndensity'))

# Manually specify an offset for points that were hard to see:
plot_data$offset = ifelse(plot_data$signal == 'deltaAF', .1, -.1)
plot_data$offset = ifelse(plot_data$trait %in% c('multiple_sclerosis','malaria'), 0, plot_data$offset)

# Point to significance threshold
cutoff = -log10(0.05)

plot_data %>%
	ggplot(aes(x=as.numeric(trait)+offset, y=log10pvalue, color=signal, shape=db)) + 
	geom_point(size=5) + 
	scale_x_continuous(breaks=1:10,labels=trait_labels) + 
	coord_flip() + 
	theme_srd() + 
	geom_hline(yintercept=cutoff, lty=2) + 
	labs(y=expression("-log"[10]~"adjusted p-value"), x='') + 
	scale_shape(name='', breaks=c('GWAS','ImmunoBase','Phenopedia'), labels=c('GWAS (n = 354)','ImmunoBase (n = 19)', 'Phenopedia (n = 1)')) + 
	theme(legend.box.just = "left") + 
	theme(legend.position=c(1,0), legend.justification=c(1,0)) + 
	scale_color_manual(name='', breaks=c('deltaAF','Fst'), values=c('#e74c3c','#2980b9'), labels=c('Top 1% ΔAF','Top 1% Fst')) 

ggsave('figs/fig2_panel4.png')
