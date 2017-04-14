source('src/helpers.R')

###
### Panel 1 - eQTL comparison
###

# Read in effect size comparison data
data = read_tsv('data/effect_sizes/effect_size_comparisons.tsv')

# Select top K eQTLs
limits = tibble(K=c(100, 500, 1000, 3000))
topK_eqtls = crossing(limits, data) %>% 
	filter(rhoRank <= K)

# Generate figure
scale_fill = scale_fill_manual(name='Population', labels=c('Sardinia','Geuvadis (Europe)','DGN (Europe)'), limits=c('SRD188','GEU188','DGN188'), values=c(color.srd, color.geu, color.dgn))
topK_eqtls %>%
	ggplot(aes(x=factor(K), y=abs(rho), fill=factor(group))) +
	geom_boxplot(outlier.shape=NA) +
	scale_fill + 
	theme_srd() +
	theme(legend.position=c(0,0),legend.justification=c(0,0)) +
	labs(x='Top N eQTLs', y='Correlation of genotype and expression') +
	scale_y_continuous(limits=c(.3, 1), breaks=seq(.3,1,.1))

ggsave('figs/fig1_panel1.png')

###
### Panel 2 - aseQTL comparison
###

# Load aseQTL data
aseqtls.srd = read_tsv('data/aseqtls/aseQTLs.SRD.sub.txt')
aseqtls.geu = read_tsv('data/aseqtls/aseQTLs.GEU.sub.txt')
aseqtls.dgn = read_tsv('data/aseqtls/aseQTLs.DGN.sub.txt')

# Parse aseQTL data
parse_aseqtls = function(aseq, study) {
	aseq$ASECHR = unlist(lapply(strsplit(aseq$ASE, '_'), function(x) {as.numeric(x[2])}))
	aseq$ASEPOS = unlist(lapply(strsplit(aseq$ASE, '_'), function(x) {as.numeric(x[3])}))
	aseq$SNPCHR = unlist(lapply(strsplit(aseq$SNP, '_'), function(x) {as.numeric(x[2])}))
	aseq$SNPPOS = unlist(lapply(strsplit(aseq$SNP, '_'), function(x) {as.numeric(x[3])}))
	aseq$BH = p.adjust(aseq$BF, method = 'BH')
	aseq$STUDY = study
	aseq$RHO = as.numeric(aseq$RHO)
	aseq$CIS_MAF = as.numeric(aseq$CIS_MAF)
	aseq$ASE_MAF = as.numeric(aseq$ASE_MAF)
	aseq = aseq[order(aseq$BH, -abs(aseq$RHO)), ]
	aseq$ID = paste(aseq$ASECHR, aseq$ASEPOS, sep = "_")
	aseq$DIST = aseq$ASEPOS - aseq$SNPPOS
	aseq$CIS_MAF = ifelse(aseq$CIS_MAF > .5, 1 - aseq$CIS_MAF, aseq$CIS_MAF)
	aseq$ASE_MAF = ifelse(aseq$ASE_MAF > .5, 1 - aseq$ASE_MAF, aseq$ASE_MAF)
	aseq$rhoRank = rank(-abs(aseq$RHO))
	return(aseq)
}

aseqtls.srd = parse_aseqtls(aseqtls.srd, "SRD188")
aseqtls.geu = parse_aseqtls(aseqtls.geu, "GEU188")
aseqtls.dgn = parse_aseqtls(aseqtls.dgn, "DGN188")
aseqtls = bind_rows(aseqtls.srd, aseqtls.geu, aseqtls.dgn)

# Select top K aseQTLs
limits = tibble(K=c(100, 250, 500, 1000))
topK_aseqtls = crossing(limits, aseqtls) %>% 
	filter(rhoRank <= K)

# Generate figure
topK_aseqtls %>%
	ggplot(aes(x=factor(K), y=abs(RHO), fill=factor(STUDY))) +
	geom_boxplot(outlier.shape=NA) +
	scale_fill + 
	theme_srd() +
	theme(legend.position=0) +
	labs(x='Top N aseQTLs', y='Correlation of genotype and ASE') + 
	scale_y_continuous(limits=c(.5, 1), breaks=seq(.5,1,.1))

ggsave('figs/fig1_panel2.png')
