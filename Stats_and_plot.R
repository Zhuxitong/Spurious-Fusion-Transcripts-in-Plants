library(tidyverse)
library(magrittr)
library(scales)
library(rstatix)
library(ggpubr)
library(colorspace)
library(tidyplots)
library(ggview)
library(gg.gap)
library(ggVennDiagram)
library(GenomicFeatures)
library(rtracklayer)
library(pheatmap)

# HC LC number -------------------------------------------------------

all_num <- read.table('HC.txt', header = T)
p <- all_num %>% 
  group_by(data) %>% mutate(freq=num/sum(num)*100) %>% 
  ungroup() %>% 
  mutate_at(.vars = 'data', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>% 
  ggplot(aes(x=data, y = num, fill = type)) + geom_bar(stat = 'identity', position = position_dodge2(), width = 0.8) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + scale_y_break(c(150,650), scales = 0.4, ticklabels = c(650,700)) + scale_y_break(c(700,3800), scales = 0.8, ticklabels = c(3850,3950)) + geom_text(mapping = aes(x = data, y = num, label = num), position = position_dodge(width = 0.9)) + scale_fill_discrete_sequential(palette='OrYel', rev=F)
ggsave(plot = p, filename = 'Number_fusion_events.pdf', width = 4, height = 4)


# Venn for sequencing tech ------------------------------------------------

isoseq_venn <- read.table('Venn/isoseq.venn.id', comment.char = '$') %>% pull(V1)
drs_venn <- read.table('Venn/drs.venn.id', comment.char = '$') %>% pull(V1)
cDNA_venn <- read.table('Venn/cDNA.venn.id', comment.char = '$') %>% pull(V1)
tech_list <- list(IsoSeq = isoseq_venn,
                  dRNA = drs_venn,
                  cDNA = cDNA_venn)
p <- ggVennDiagram(x = tech_list, shape_id = '301', label = 'count') + scale_fill_continuous_sequential('PurpOr')
p + canvas(5,3)
ggsave(plot = p, filename = 'Venn/Sequencing_tech_overlap.pdf', width = 5, height = 3)

venn <- Venn(tech_list)
venn_data <- process_data(venn, shape_id = '301')
plot_venn(venn_data, label = 'count', edge_size = 0.75)
genes <- paste0("gene",1:1000)
set.seed(20210302)
gene_list <- list(A = sample(genes,100),
                  B = sample(genes,200))
p <- ggVennDiagram(gene_list, shape_id = '201', label = 'count')
ggsave(plot = p, filename = 'Venn/dre_rep_overlap.pdf', width = 4, height = 4)
p <- ggVennDiagram(gene_list, shape_id = '201f', label = 'count')
ggsave(plot = p, filename = 'Venn/cDNA_rep_overlap.pdf', width = 4, height = 4)


# fusion type in cDNA both replicates -------------------------------------

p <- data.frame(
  stringsAsFactors = FALSE,
                type = c("HC", "LC", "NS", "PRT"),
                num = c(1L, 4L, 14L, 0L)
) %>% ggplot(aes(x=type, y = num)) + geom_bar(stat = 'identity', width = 0.6, fill='#6AADC1') + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + geom_text(aes(label = num))
ggsave(plot = p, filename = 'Venn/cDNA_rep_overlap_type.pdf', width = 3, height = 4)
'#8577AE'


# fusion type in RNA-seq replicates ----------------------------------
p <- data.frame(
  stringsAsFactors = FALSE,
                type = c("Read-through","Duplication",
                       "Inversion","Translocation"),
                num = c(55L, 2L, 4L, 10L),
                order = 1:4
) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('Read-through','Duplication','Inversion','Translocation'))) %>% 
  ggplot(aes(x=type, y = num, fill=order)) + geom_bar(stat = 'identity', width = 0.6) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'right', axis.text.x = element_text(angle = 30, vjust=0.6)) + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + geom_text(aes(label = num)) + scale_fill_gradient(low = '#A6C9D4', high = '#6AADC1')
ggsave(plot = p, filename = 'Venn/RNAseq_rep_overlap_type.pdf', width = 3, height = 4)


# TPM cor -----------------------------------------------------------------
exp_isoseq <- read.table('01exp_readcounts/mh63_isoseq_salmon_readCount.tsv', header = T)
exp_drs <- read.table('01exp_readcounts/mh63_drs_salmon_readCount.tsv', header = T)
exp_cDNA <- read.table('01exp_readcounts/mh63_cDNA_salmon_readCount.tsv', header = T)
pdf(file = '01exp_readcounts/readCount_cor.pdf', width = 5, height = 4)
merge(x = exp_isoseq %>% select(1,2), y = exp_drs %>% select(c(1:3)), by = 'Name') %>% merge(x=., y=exp_cDNA %>% select(1:3), by='Name') %>% column_to_rownames(var = 'Name') %>% cor(method = 's') %>% pheatmap(color = sequential_hcl(palette = 'RedOr', n = 100, rev = T))
dev.off()

merge(x = exp_isoseq %>% select(1,2), y = exp_drs %>% select(c(1:3)), by = 'Name') %>% merge(x=., y=exp_cDNA %>% select(1:3), by='Name') %>% column_to_rownames(var = 'Name') %>% cor(method = 's')


# Tech reads num ----------------------------------------------------------
p <- read.table('02data_quality/Read_num.txt', header = T) %>% 
  mutate(n=num/1e+06) %>%
  mutate_at(.vars = 'type', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA','RNAseq') %>% rev())) %>% 
  filter(type!="RNAseq") %>% 
  ggplot(aes(x=type, y = n, shape=rep, color=rep)) + geom_point(size=4) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'right') + xlab('') + ylab(expression(paste('Number of reads (', 10^6,')'))) + scale_color_manual(values = c('#C55B5A', '#435A8F')) + scale_y_continuous(breaks = breaks_extended(n = 8))+ coord_flip()
ggsave(plot = p, filename = '02data_quality/Read_num_1.pdf', width = 4, height = 3)


# Reads length and N50 ----------------------------------------------------

reads_len <- read.table('02data_quality/reads.len', header = F) %>% set_colnames(c('sample','len')) %>% mutate(n=len/1000)
reads_len %>% group_by(sample) %>% summarise(n=mean(len))
p <- reads_len %>% 
  mutate_at(.vars = 'sample', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>%
  ggbarplot(data = ., x = 'sample', y = 'len', add = c('mean_sd'), width = 0.3, fill = 'sample') + theme(axis.text.y = element_text(angle=90, hjust = 0.5),) + xlab('') + ylab('Mean reads length') + theme(legend.position = 'none') + scale_fill_manual(values = c('#81B3D6','#6D60BB','#F57F73')) + scale_y_continuous(breaks = seq(0,4500,500))
ggsave(plot = p, filename = '02data_quality/reads.len.pdf', width = 3, height = 4)

p <- data.frame(
  stringsAsFactors = FALSE,
                sample = c("IsoSeq", "dRNA", "cDNA"),
                len = c(2635L, 1074L, 1022L)
) %>% mutate_at(.vars = 'sample', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>%
  ggbarplot(data = ., x = 'sample', y = 'len', width = 0.3, fill = 'sample') + theme(axis.text.y = element_text(angle=90, hjust = 0.5),) + xlab('') + ylab('Mean reads length') + theme(legend.position = 'none') + scale_fill_manual(values = c('#B3D2E7','#BEBBD8','#F7B2AB'))
ggsave(plot = p, filename = '02data_quality/reads.len.N50.pdf', width = 3, height = 4)


# Transcript novelty ------------------------------------------------------

trans_novelty <- read.table('02data_quality/transcripts_type.txt', header = T)
trans_novelty %>% 
  group_by(type) %>% mutate(freq=num/sum(num)*100) %>% pivot_wider(id_cols = type,names_from = cate, values_from = freq)
p <- trans_novelty %>% group_by(type) %>% mutate(freq=num/sum(num)*100) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>% 
  ggplot(aes(x=type, y = freq, fill=cate)) + geom_bar(stat = 'identity', width = 0.5, ) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + xlab('') + ylab('Transcript novelty (%)') + scale_y_continuous(expand = c(0,0)) +  scale_fill_manual(values = c('#7F80A3','#FDBF70','#7F7F7F','#EB786B','#B2AAD3','#70B0E0','#7FBDB0'))
ggsave(plot = p, filename = '02data_quality/transcripts_type.pdf', width = 3, height = 4)


# Sequencing tech frag_mapped  --------------------------------------------
frag_mapped <- read.table('02data_quality/frag_mapped.txt', header = F) %>% set_colnames(c('type','pct'))
p <- frag_mapped %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>% 
  filter(pct!='NA') %>% 
  ggplot(aes(x=pct, fill = type)) + geom_histogram(position = position_dodge2(), bins = 30) + scale_x_continuous(limits = c(0.5,1)) + scale_y_break(c(2000,8000)) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'NA') + xlab('') + ylab('Number of mapped reads (10^3)')
ggsave(plot = p, filename = '02data_quality/frag_mapped.pdf', width = 4, height = 4)

frag_mapped %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>% 
  filter(pct>=0.9) %>% 
  group_by(type) %>% count()

# Distribution of partner genes -------------------------------------------

partner_genes_chr <- read.table('02data_quality/Circos/Distribution_of_partner_genes_on_chr.txt', header = T)

cor.test(partner_genes_chr$IsoSeq, partner_genes_chr$genome)
cor.test(partner_genes_chr$dRNAseq, partner_genes_chr$genome)
cor.test(partner_genes_chr$cDNAseq, partner_genes_chr$genome)

partner_genes_chr_all <- partner_genes_chr %>% mutate(n=IsoSeq+dRNAseq+cDNAseq)
cor.test(partner_genes_chr_all$n, partner_genes_chr_all$genome)


# SHS length distribution  ------------------------------------------------

break_type_cumsum <- break_type %>%
  filter(breakpoint=="SHS") %>%
  select(len) %>% mutate(a=if_else(len>=16, 16, len)) %>% 
  group_by(a) %>% count() %>% ungroup() %>%
  mutate(cum_n=cumsum(n)) %>% mutate(freq=cum_n/sum(n)*100)
  
p <- break_type %>%
  filter(breakpoint=="SHS") %>%
  select(len) %>% mutate(a=if_else(len>=16, 16, len)) %>% 
  group_by(a) %>% count() %>% ungroup %>% 
  ggplot(aes(x=a,y=n)) + geom_bar(stat = 'identity', fill='#BAC9EA') + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black')) + scale_y_continuous(limits = c(0,500),expand = c(0,0), sec.axis = sec_axis(~./100,breaks = seq(0,5,1), labels = seq(0,1,0.2))) + scale_x_continuous(breaks = seq(3,15,2), labels = seq(3,15,2)) + xlab('SHS length (bp)') + ylab('Number of fusions') + geom_line(data = break_type_cumsum, aes(x=a, y = freq*5), color='#5F8BEE', size=1) + geom_point(data = break_type_cumsum, aes(x=a, y = freq*5), color='#5F8BEE', size=2) + geom_point(data = shs_len_random_sum, aes(x=len,y = m), shape=18, color='#4C4C4C', size=2.5) + geom_errorbar(data = shs_len_random_sum, aes(x=len,ymin = m-se, ymax = m+se), width = 0.5, color='#4C4C4C')
ggsave(plot = p, filename = '03shs/SHS_length_distribution.pdf', width = 4, height = 3)

shs_len_random <- read.table('03shs/random_shs.count', header = F) %>% set_colnames(c('cishu','len','n'))
shs_len_random_sum <- shs_len_random %>% 
  filter(len<20 & len!=3) %>% 
  mutate(num=round(n*(3261/725),digits = 0)) %>% 
  group_by(len) %>% summarise(m=median(num),me=mean(num),sd=sd(num), n=n(), se=sd/sqrt(n)) %>% ungroup() %>% add_row(len=15,m=0,me=0,sd=0,n=0,se=0) %>% add_row(len=16,m=0,me=0,sd=0,n=0,se=0)
shs_len_random %>% 
  filter(len<20 & len!=3) %>% 
  mutate(num=round(n*(3261/725),digits = 0)) %>% 
  group_by(cishu) %>% summarise(great4=sum(num)) %>% 
  get_summary_stats()


# SHS sequence freq  ------------------------------------------------------

shs_seq <- read.table('03shs/all_shs.tab', header = F) %>% set_colnames(c('id','len','seq'))
shs_seq %>% select(len, seq) %>% 
  filter(len>=4) %>% 
  group_by(len, seq) %>% count() %>% 
  ungroup %>% group_by(len) %>% 
  arrange(len, desc(n), .by_group = T) %>%
  slice_head(n=50) %>% write_tsv(file = '03shs/all_shs.txt')

# SHS seqLogo -------------------------------------------------------------
library(DiffLogo)
library(seqLogo)
logo_shs_4 <- getPwmFromFastaFile(filename = '03shs/logo_shs_4bp.fa', alphabet = )
logo_rand_4 <- getPwmFromFastaFile('03shs/logo_rand_4bp.fa')
seqLogo::seqLogo(logo_shs_4, ic.scale = F)
seqLogo::seqLogo(logo_rand_4, ic.scale = F)
DiffLogo::seqLogo(pwm = logo_shs_4, stackHeight = sumProbabilities, baseDistribution = probabilities)
diffLogoFromPwm(pwm1 = logo_shs_4, pwm2 = logo_rand_4, ymin = -0.1, ymax = 0.1, stackHeight = shannonDivergence)
diffLogoFromPwm(pwm1 = logo_shs_4, pwm2 = logo_rand_4, ymin = -0.1, ymax = 0.1, stackHeight = sumOfAbsICDifferences)
stackHeight = sumOfAbsProbabilityDifferences
stackHeight = sumProbabilities
pwmDivergence(pwm_left = logo_shs_4, pwm_right = logo_rand_4,divergence = shannonDivergence)

diffLogoObj <- createDiffLogoObject(pwm1 = logo_shs_4, pwm2 = logo_rand_4, stackHeight = shannonDivergence)
enrichDiffLogoObjectWithPvalues(diffLogoObj = diffLogoObj, n1 = 824, n2 = 824, numberOfPermutations = 100)

logo_shs_4 <- getPwmFromFastaFile(filename = '03shs/logo_shs_4bp.fa', alphabet = RNA)
logo_rand_4 <- getPwmFromFastaFile('03shs/logo_rand_4bp.fa', alphabet = RNA)
diffLogoObj <- createDiffLogoObject(pwm1 = logo_rand_4, pwm2 = logo_shs_4, stackHeight = shannonDivergence, alphabet = RNA, baseDistribution = normalizedDifferenceOfProbabilities)
diffLogoObj <- enrichDiffLogoObjectWithPvalues(diffLogoObj = diffLogoObj, n1 = 412, n2 = 412, numberOfPermutations = 100)
diffLogoObj$pvals
pwmDivergence(pwm_left = logo_shs_4, pwm_right = logo_rand_4,divergence = shannonDivergence)
pdf(file = '03shs/logo_4bp.pdf', width = 5, height = 5)
diffLogo(diffLogoObj = diffLogoObj, ymin = -0.01, ymax = 0.01, sparse = T)
dev.off()

logo_shs_5 <- getPwmFromFastaFile(filename = '03shs/logo_shs_5bp.fa', alphabet = RNA)
logo_rand_5 <- getPwmFromFastaFile('03shs/logo_rand_5bp.fa', alphabet = RNA)
diffLogoObj <- createDiffLogoObject(pwm1 = logo_rand_5, pwm2 = logo_shs_5, stackHeight = shannonDivergence, alphabet = RNA, baseDistribution = normalizedDifferenceOfProbabilities)
diffLogoObj <- enrichDiffLogoObjectWithPvalues(diffLogoObj = diffLogoObj, n1 = 370, n2 = 370, numberOfPermutations = 100)
diffLogoObj$pvals
pwmDivergence(pwm_left = logo_shs_5, pwm_right = logo_rand_5,divergence = shannonDivergence)
pdf(file = '03shs/logo_5bp.pdf', width = 5, height = 5)
diffLogo(diffLogoObj = diffLogoObj, ymin = -0.01, ymax = 0.01,sparse = T)
dev.off()

logo_shs_6 <- getPwmFromFastaFile(filename = '03shs/logo_shs_6bp.fa', alphabet = RNA)
logo_rand_6 <- getPwmFromFastaFile('03shs/logo_rand_6bp.fa', alphabet = RNA)
diffLogoObj <- createDiffLogoObject(pwm1 = logo_rand_6, pwm2 = logo_shs_6, stackHeight = shannonDivergence, alphabet = RNA, baseDistribution = normalizedDifferenceOfProbabilities)
diffLogoObj <- enrichDiffLogoObjectWithPvalues(diffLogoObj = diffLogoObj, n1 = 350, n2 = 350, numberOfPermutations = 100)
diffLogoObj$pvals
pwmDivergence(pwm_left = logo_shs_6, pwm_right = logo_rand_6,divergence = shannonDivergence)
pdf(file = '03shs/logo_6bp.pdf', width = 5, height = 5)
diffLogo(diffLogoObj = diffLogoObj, ymin = -0.01, ymax = 0.01,sparse = T)
dev.off()

# expression for all breakpoint types -------------------------------------
exp_breakpoints_type <- read.table('03shs/exp_all_types.txt', header = F) %>% set_colnames(c('geneID', 'type','tpm'))
exp_breakpoints_type %>% 
  group_by(type) %>% get_summary_stats()

tmp_median <- exp_breakpoints_type %>% 
  mutate(ltpm=log10(tpm+0.1)) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned','Random'))) %>% 
  ggplot(aes(x=type, y = ltpm)) + geom_boxplot()
dat <- ggplot_build(tmp_median)$data[[1]]

p <- exp_breakpoints_type %>% 
  mutate(ltpm=log10(tpm+0.1)) %>% 
  ggviolin(x = 'type', y = 'ltpm', width = 0.5, add = c('boxplot'), draw_quantiles = F, fill='type',color='NA', add.params = list(fill='white', color='NA')) + xlab('') + ylab(expression(paste(Log[10](TPM+0.1)))) + theme(legend.position = 'none') + scale_fill_manual(values = c('#73C594','#A291B2','#91AEDB','#EF8A47','#8A8A8A')) + geom_segment(data = dat, aes(x=xmin+0.25, xend=xmax-0.25, y=middle, yend=middle), color='#4C4C4C') + stat_compare_means(comparisons = list(c('SHS','Random'),c('Adjacent','Random'),c('Unknown','Random'), c('Misaligned','Random')), method = "wilcox.test", label = 'p.signif')
ggsave(plot = p, filename = '03shs/exp_all_types.pdf', width = 5, height = 6)


# SHS tissue-specific expression  -----------------------------------------
all_gene_exp <- read.table('03shs/gene_tpm_4tissue.matrix', header = F) %>% set_colnames(c('geneID','flagleaf','leaf','panicle','root'))
all_gene_exp$tau <- all_gene_exp %>% select(-geneID) %>%
  apply(1,FUN = function(x){sum(1-x/max(x))/(length(x)-1)}) %>% round(digits = 2)
all_gene_exp %>% 
  filter(tau>=0.8) %>% nrow()

library(ggalluvial)
merge(exp_breakpoints_type, all_gene_exp, by = 'geneID') %>%
  mutate_at(.vars = 'tau', ~ifelse(is.nan(.x), 0, .x)) %>% 
  select(type, tau) %>% 
  mutate(sp=if_else(tau>=0.7, 'Y','N')) %>%
  group_by(type, sp) %>% count() %>% 
  ungroup() %>% group_by(type) %>% mutate(freq=n/sum(n)*100) %>% 
  write_tsv('03shs/shs_tissueSpecific.txt')

read.table('03shs/shs_tissueSpecific.txt', header = T) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned','Random'))) %>% 
  ggplot(aes(x=type, y = freq, fill=sp, stratum = sp, alluvium = sp)) + geom_stratum(color='NA') + geom_alluvium(curve_type = 'linear') + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none') + scale_y_continuous(expand = c(0,0)) + xlab('') + ylab('Percentage of fuison events') + coord_flip()

  read.table('03shs/shs_tissueSpecific.txt', header = T) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned','Random') %>% rev())) %>% 
  mutate(new=if_else(sp=='Y', 0-freq, freq)) %>% 
  ggplot(aes(x=type, y = new, fill=sp, stratum = sp, alluvium = sp)) + geom_stratum(color='NA') + geom_alluvium(curve_type = 'linear') + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none') + scale_y_continuous(expand = c(0,0)) + xlab('') + ylab('Percentage of fuison events') + scale_fill_manual(values = c('#8B7BB8', '#EE754C') %>% rev)
p + canvas(3,4)
ggsave(plot = p, filename = '03shs/shs_tissue_specific.pdf', width = 4, height = 4)


# Unknown sequences length ------------------------------------------------
break_type_new <- read.table('03shs/breakpoint_newType.txt', header = F) %>% set_colnames(c('tech','reads','class','start','end','len','bp','bp_new'))
p <- break_type_new %>% 
  select(tech, len, bp, bp_new) %>% 
  filter(bp_new=="Unknown") %>%
  ggplot(aes(x=len, y = after_stat(scaled))) + geom_density(color="#5F8BEE", fill="#B9C9E9", linewidth=1) +  theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black')) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + xlab('Length of unknown sequences') + geom_vline(xintercept = 26, size=1, linetype = 'dashed', color='#4C4C4C') + ylab('Density')
p + canvas(5,4)
ggsave(plot = p, filename = '03shs/length_unknown_seq.pdf', width = 5, height = 4)

# shs length include >15 --------------------------------------------------
break_type_new %>% 
  select(tech, len, bp_new) %>% 
  mutate(bp_len=if_else(bp_new=='SHS' & len>15, 'SHS(>15)', bp_new)) %>%
  group_by(tech, bp_len) %>% count() %>%
  ungroup() %>% group_by(bp_len) %>% summarise(total=sum(n))


# unknown sequences GC ----------------------------------------------------
unknown_gc <- read.table('03shs/pacbio_flnc.unknown.all.gc', header = F) %>% set_colnames(c('id','seq','at','gc'))
unknown_gc_random <- read.table('03shs/pacbio_flnc.unknown.random.gc', header = F) %>% set_colnames(c('gc')) %>% add_column(type='random')
p <- unknown_gc %>% select(gc) %>% add_column(type='unknown') %>% bind_rows(., unknown_gc_random) %>% ggboxplot(x = 'type', y = 'gc', width = 0.3, notchwidth = 0.6, fill='#ABA6E6', color='#ABA6E6', notch=T) + stat_compare_means(method = 'wilcox.test') + xlab('') + ylab('GC (%)') #+ canvas(4,4)
ggsave(plot = p, filename = '03shs/GC_unknown_seq.pdf', width = 4, height = 4)

unknown_gc %>% select(gc) %>% add_column(type='unknown') %>% bind_rows(., unknown_gc_random) %>% 
  group_by(type) %>% get_summary_stats()

# Unknown seq source ------------------------------------------------------

generate_random_dna <- function(M) {
  bases <- c("A", "C", "G", "T")
  dna_sequence <- sample(bases, M, replace = TRUE, prob = c(0.2306,0.2694,0.2694,0.2306))
  dna_sequence <- paste(dna_sequence, collapse = "")
  return(dna_sequence)
}

M <- 5
out <- data.frame()
for (i in 1:422){
  s1 <- generate_random_dna(M)
  fa <- data.frame(a=paste0("random_",i), b=paste0(s1))
  out <- bind_rows(out, fa)
}
write_tsv(x = out, file = '03shs/unknown_random_5bp.txt', col_names = F)

# Adjacent sequence -------------------------------------------------------

break_type_new %>% 
  filter(bp_new=='Adjacent') %>%
  filter(class=='LC' | class=='NS') %>% nrow


# breakpoints type reads length -------------------------------------------

bp_reads_len_cDNA <- read.table('03shs/reads_length.cDNA', header = F) %>% set_colnames(c('id','type','len'))
bp_reads_len_cDNA.rand <- read.table('03shs/reads_length.cDNA.random', header = F) %>% set_colnames(c('tech','len')) %>% add_column(type='Random')

bind_rows(bp_reads_len_cDNA %>% select(type,len), bp_reads_len_cDNA.rand %>% select(type,len)) %>% 
  filter(type!="SHS(3bp)") %>% 
  group_by(type) %>% get_summary_stats()

tmp_median <- bind_rows(bp_reads_len_cDNA %>% select(type,len), bp_reads_len_cDNA.rand %>% select(type,len)) %>% 
  filter(type!="SHS(3bp)") %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned','Random'))) %>% 
  ggplot(aes(x=type,y=len)) + geom_boxplot()
dat <- ggplot_build(tmp_median)$data[[1]]

p <- bind_rows(bp_reads_len_cDNA %>% select(type,len), bp_reads_len_cDNA.rand %>% select(type,len)) %>% 
  filter(type!="SHS(3bp)") %>% 
  filter(len<4000) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned','Random'))) %>% 
  ggviolin(x = 'type', y = 'len', width = 0.6, add = c('mean'),draw_quantiles = F, fill='type',color='NA') + xlab('') + ylab('Reads length (kb)') + theme(legend.position = 'none') + scale_fill_manual(values = c('#73C594','#A291B2','#91AEDB','#EF8A47','#8A8A8A')) + stat_compare_means(ref.group = 'Random', method = "wilcox.test", label = 'p.signif') + geom_segment(data = dat, aes(x=xmin+0.25, xend=xmax-0.25, y=middle, yend=middle), color='#4C4C4C') + scale_y_continuous(limits = c(0,4000))
ggsave(plot = p, filename = '03shs/reads_length.pdf', width = 5, height = 6)


# trans-splicing and helitron ---------------------------------------------
helitron_random <- read.table('04error_Alignment/random.helitron.freq', header = F) %>% set_colnames(c('cishu','num'))
p <- helitron_random %>% 
  ggplot(aes(x=num, y=after_stat(scaled))) + geom_density(color = '#6276DD', fill='#D1D4F4', size=1.2) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none') + scale_y_continuous(expand = c(0,0)) + xlab('Number of Helitron-related alignment') + ylab('Density') + geom_vline(xintercept = c(94,122), size=1.1, linetype=2)
ggsave(plot = p, filename = '04error_Alignment/random.helitron.pdf', width = 4, height = 4)


# adjacent types ----------------------------------------------------------
adjacent_types <- read.table('04error_Alignment/adjacent_types.txt', header = F) %>% set_colnames(c('main_type','type'))
adjacent_types %>% 
  group_by(main_type, type) %>% count()


# adjacent types alignment and coverage -----------------------------------
library(ggExtra)

adjacent_len_cov <- read.table('04error_Alignment/alignment_len_coverage_plot.txt', header = F) %>% set_colnames(c('len','cov','main_type','type'))

adjacent_len_cov_plot <- adjacent_len_cov %>% 
  mutate(new_type=if_else(type=='ShortAlign' | type=='specific-SV','Short','Others')) %>% filter(cov<75)
sp <- ggscatter(adjacent_len_cov_plot, x = 'len', y = 'cov', color = 'new_type', size=2, alpha=.8) + theme(legend.position = 'none') + scale_color_manual(values = c('#F57F73','#81B3D6') %>% rev) + xlab('Length of minimal alignment') + ylab('Coverage of minimal alignment')
p <- ggMarginal(sp, type="density", groupColour = T, groupFill = T)
ggsave(plot = p, filename = '04error_Alignment/alignment_len_coverage_plot.pdf', width = 4, height = 3)

# Misaligned alignment score ----------------------------------------------
alignment_minQuality <- read.table('04error_Alignment/Misaligned_alignmentScore.txt', header = T, comment.char = '$')
alignment_minQuality %>% 
  select(type, score) %>% 
  group_by(type) %>% get_summary_stats()

p <- ggbarplot(data = alignment_minQuality %>% mutate_at(.vars = 'type', ~factor(.x, levels=c('shs','join','unknown','error'))), x = 'type', y = 'score', add = c("mean_se"), width = 0.3, fill = 'type') + scale_y_continuous(expand = c(0,0)) + xlab('') + ylab('Normalized alignment score') + theme(legend.position = 'none') + scale_fill_manual(values = c('#73C594','#A291B2','#91AEDB','#EF8A47'))
p + canvas(4,4)
ggsave(plot = p, filename = '04error_Alignment/Misaligned_alignmentScore.pdf', width = 4, height = 4)
alignment_minQuality %>% wilcox_test(score~type)

# Misaligned alignment mapped length --------------------------------------

alignment_lenMin <- read.table('04error_Alignment/Misaligned_alignmentLen.txt',header = T, comment.char = '$')
alignment_lenMin %>% 
  group_by(type) %>% get_summary_stats()

tmp_median <- alignment_lenMin %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('shs','join','unknown','error'))) %>%
  ggplot(aes(x=type, y=len)) + geom_boxplot()
dat <- ggplot_build(tmp_median)$data[[1]]
p <- alignment_lenMin %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('shs','join','unknown','error'))) %>%
  ggplot(aes(x=type, y=len)) + geom_boxplot(width=0.4, fill='#1C7DC1', color='#1C7DC1',notch = T, notchwidth = 0.6) + scale_y_continuous(oob = scales::oob_squish, breaks = scales::breaks_extended(n = 6)) + theme_pubr() + theme(axis.title.x = element_blank()) + ylab('') + geom_segment(data = dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color='white')
ggsave(plot = p, filename = '04error_Alignment/Misaligned_alignmentLen.pdf', width = 4, height = 4)
alignment_lenMin %>% wilcox_test(len~type)


# error alignment parental cov --------------------------------------------

error_align_parental_cov <- read.table('04error_Alignment/error_alignment.txt', header = T) %>% set_colnames(c('g1','g2','tech','type'))

error_align_parental_cov %>% 
  group_by(type) %>% get_summary_stats()

p <- error_align_parental_cov %>% 
  filter(g1<=100 & g2<=100) %>% 
  select(type, g1,g2) %>% pivot_longer(!type, names_to = 'gene', values_to = 'cov') %>% filter(type!="SHS(3bp)") %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned'))) %>% 
  ggplot(aes(x = type, y = cov, fill = gene)) + geom_boxplot(width=0.6, outlier.colour = 'NA', notch = T, notchwidth = 0.7) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none') + xlab('') + ylab('Parental gene coverage (%)') + scale_fill_manual(values = c(c('#91AEDB','#EF8A47')))
ggsave(plot = p, filename = '04error_Alignment/error_alignment.pdf', width = 4, height = 4)

error_align_parental_cov %>% 
  filter(g1<=100 & g2<=100) %>% 
  select(type, g1,g2) %>% 
  filter(type!="SHS(3bp)") %>% 
  pull(type) %>% table()


# error alignment and similarity ------------------------------------------
adj_xtab <- as.table(rbind(c(4,40-4),c(375,4718-375)))
dimnames(adj_xtab) <- list(group=c('adj','total'), si=c('Y','N'))
fisher.test(adj_xtab)

shs_xtab <- as.table(rbind(c(3,40-3),c(3481,4718-3481)))
dimnames(shs_xtab) <- list(group=c('shs','total'), si=c('Y','N'))
fisher.test(shs_xtab)
fisher_test(shs_xtab,detailed = T)

mis_xtab <- as.table(rbind(c(32,40-32),c(440,4718-440)))
dimnames(mis_xtab) <- list(group=c('mis','total'), si=c('Y','N'))
fisher.test(mis_xtab)
fisher_test(mis_xtab,detailed = T)

unknown_xtab <- as.table(rbind(c(1,40-1),c(422,4718-422)))
dimnames(unknown_xtab) <- list(group=c('unknown','total'), si=c('Y','N'))
fisher.test(unknown_xtab)


p <- data.frame(
  stringsAsFactors = FALSE,
                type = c("SHS", "Adjacent", "Unknown", "Misaligned"),
                enrich = c(-5.116128, 0.3637183, -1.937482, 5.280208),
                p = c(1.86e-18, 0.5562, 0.2569, 1.06e-25)
) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned') %>% rev())) %>% 
  ggplot(aes(x=type, y = enrich, fill=type)) + geom_bar(stat = 'identity', color=NA, width = 0.3) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none') + coord_flip() + geom_hline(yintercept = 0) + ylab('Log2(Enrichment)') + xlab('')
ggsave(plot = p, filename = '04error_Alignment/error_alignment_repeatGenes.pdf', width = 4, height = 4)


# breakpoints and network -------------------------------------------------

break_type_new %>% 
    select(tech, bp_new) %>% 
    group_by(tech, bp_new) %>% count() %>% 
    ungroup %>% group_by(tech) %>% mutate(freq=round(n/sum(n)*100, digits = 2)) %>% 
    pivot_wider(id_cols = tech, names_from = bp_new, values_from = freq, values_fill = 0)

break_type_new %>% 
  group_by(bp_new) %>% select(bp_new) %>% count()

break_type_new %>% 
  filter(grepl('OsMH_11G0087600', reads)) %>% 
  select(tech, bp_new) %>% 
  group_by(tech, bp_new) %>% count()


# net degree --------------------------------------------------------------

gene_degree <- read.table('04error_Alignment/degree.txt', header = F) %>% set_colnames(c('tech','id','degree','exp'))
gene_degree %>% 
  ggplot(aes(x=degree, color=tech)) + stat_ecdf() + scale_x_log10()

gene_degree %>%
  group_by(tech) %>% 
  summarise(
    mean_degree = mean(degree),
    median_degree = median(degree),
    total_edges = sum(degree)/2,
    top_genes_30perc = sum(head(sort(degree, decreasing = TRUE), n = 15)) / sum(degree)
  )

library(poweRlaw)
degree_dist <- displ$new(gene_degree$degree)
est <- estimate_xmin(degree_dist)
degree_dist$setXmin(est)

bs_p <- bootstrap_p(degree_dist, threads = 4)
print(bs_p)

ggplot(gene_degree, aes(x = degree, fill=tech)) +
  geom_density(alpha = 0.5) + labs(title = "Gene Degree Distribution", x = "Gene degree", y = "Density") + scale_x_log10() +theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'right') + scale_y_continuous(expand = c(0,0))

p <- gene_degree %>% 
  mutate(new_degree=if_else(degree>=10, 10, degree)) %>% 
  group_by(tech, new_degree) %>% count() %>% 
  mutate_at(.vars = 'tech', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>%
  ggplot(aes(x=new_degree, y=n, fill = tech)) + geom_bar(stat = 'identity', position = position_dodge2()) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + scale_y_continuous(expand = c(0,0)) + xlab('Gene degree') + ylab('Gene number') + scale_fill_manual(values = c('#81B3D6','#6D60BB','#F57F73')) + scale_y_break(breaks = c(800,2000), scales = 0.5) + scale_x_continuous(breaks = breaks_extended(n = 10))
p + canvas(5,4)
ggsave(plot = p, filename = '04error_Alignment/gene_degree.pdf', width = 5, height = 4)

gene_degree %>% 
  mutate(new_degree=if_else(degree>=10, 10, if_else(degree!=1, 2, degree))) %>%
  select(tech, new_degree) %>% 
  group_by(tech, new_degree) %>% count()

degree_freq <- as.data.frame(table(gene_degree$degree)) %>%
  mutate(Degree = as.numeric(Var1), 
         Cumulative = cumsum(Freq)/sum(Freq))

ggplot(degree_freq, aes(x = Degree, y = 1-Cumulative)) +
  geom_point() +
  scale_x_log10() + 
  scale_y_log10() +
  labs(title = "Log-Log Degree Distribution", 
       x = "Degree (log10)", y = "P(X ≥ x) (log10)") +
  theme_minimal()

gene_ranked <- gene_degree %>%
  group_by(tech) %>% 
  arrange(desc(degree)) %>%
  mutate(
    cum_perc = cumsum(degree)/sum(degree),
    gene_label = ifelse(row_number() <= 10, id, "")
  )

gene_ranked_100 <- bind_rows(gene_ranked %>% filter(tech=='IsoSeq') %>% slice(1:100) %>% add_column(x=sprintf('%03d', 1:100)), gene_ranked %>% filter(tech=='dRNA') %>% slice(1:100) %>% add_column(x=sprintf('%03d', 1:100)), gene_ranked %>% filter(tech=='cDNA') %>% slice(1:100) %>% add_column(x=sprintf('%03d', 1:100)))

p <- gene_ranked_100 %>% 
  mutate_at(.vars = 'tech', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>% 
  ggplot(aes(x = x, y = degree)) +
  geom_bar(aes(fill=tech), stat = "identity") +
  geom_line(aes(y = cum_perc * 919, color=tech, group = tech)) +
  scale_y_continuous(
    name = "Degree",
    expand = c(0,0),
    sec.axis = sec_axis(~./max(gene_ranked$degree), 
                        name = "Cumulative Percentage")) +
  labs(x = "Gene (Ordered by Degree)") +
  theme(axis.text.x = element_text(angle = 90, size=4)) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top', axis.text.x = element_text(angle = 90, size=5)) + scale_fill_manual(values = c('#81B3D6','#6D60BB','#F57F73')) + scale_color_manual(values = c('#81B3D6','#6D60BB','#F57F73'))
p+canvas(7,4)
ggsave(plot = p, filename = '04error_Alignment/gene_degree_cum.pdf', width = 7, height = 4)

gene_ranked_10 <- bind_rows(gene_ranked %>% filter(tech=='IsoSeq') %>% slice(1:10) %>% add_column(x=sprintf('%03d', 1:10)), gene_ranked %>% filter(tech=='dRNA') %>% slice(1:10) %>% add_column(x=sprintf('%03d', 1:10)), gene_ranked %>% filter(tech=='cDNA') %>% slice(1:10) %>% add_column(x=sprintf('%03d', 1:10)))

gene_ranked_10 %>% write_tsv('04error_Alignment/top_10_gene_degree.txt')

p <- read.table('04error_Alignment/top_10_gene_degree.txt', header = T) %>% 
  mutate_at(.vars = 'tech', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>%
  mutate(category_sorted=reorder_within(id, degree, tech)) %>%
  ggplot(aes(x=category_sorted, y=degree)) + geom_bar(aes(fill=log10(exp)), stat = "identity", width = 0.8) + facet_grid(tech~., scales = 'free_y') + coord_flip() + xlab('') + ylab('Gene degree') + geom_text(aes(label=degree), hjust=-0.5, size=3) + theme_pubr(base_size = 8) + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + scale_fill_continuous_sequential(palette='RedOr')
p+canvas(4,7)
ggsave(plot = p, filename = '04error_Alignment/top_10_gene_degree.pdf', width = 4, height = 7)


p <- gene_degree %>% 
  ggplot(aes(x=degree %>% log2(), y = log10(exp+0.1))) + geom_point(color='#AAA6AD') + geom_smooth(method = 'lm', formula = y~x, fill="#FEEAA5", color='#EFC000FF') + scale_x_continuous(breaks = breaks_extended(n = 6), labels = 2^(seq(0,10,2))) + xlab('Gene degree') + ylab('Log10(TPM+0.1)') + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'))
p+canvas(4,4)
ggsave(plot = p, filename = '04error_Alignment/gene_degree_vs_exp.pdf', width = 4, height = 4)


gene_degree %>% 
  filter(degree>1) %>% 
  mutate(ldegree=log2(degree),lexp=log10(exp+0.1)) %>% 
  cor_test(ldegree, lexp, method = 's')

# chr10 OsMH_10G0000800 gene -----------------------------------------------

p <- read.table('05network/OsMH_10G0000800_type.txt', header = F) %>% 
  mutate_at(.vars = 'V1', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>%
  mutate_at(.vars = 'V2', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned','SHS(3bp)'))) %>% 
  ggplot(aes(x=V2, y=V3, fill=V2)) + geom_bar(stat = 'identity', width = 0.5) + facet_grid(.~V1, scales = 'free', space = 'free') + scale_y_break(breaks = c(40,170)) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text.x = element_blank(), legend.position = 'none', axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_blank())+ ylab('Number') + geom_text(aes(label = V3, vjust=-0.5)) + scale_fill_manual(values = c('#73C594','#A291B2','#91AEDB','#EF8A47','#8A8A8A'))
p + canvas(7,2)
ggsave(plot = p, filename = '05network/OsMH_10G0000800_type.pdf', width = 7, height = 2)


# correct and uncorrected DRS reads ---------------------------------------

p <- data.frame(
  stringsAsFactors = FALSE,
              type = c("Adjacent", "Misaligned", "SHS", "Unknown"),
           correct = c(38L, 77L, 20L, 11L),
         uncorrect = c(819L, 3046L, 314L, 633L)
) %>% pivot_longer(-type, names_to = 'correct', values_to = 'n') %>% 
  group_by(correct) %>% mutate(freq=round(n/sum(n)*100, digits = 2)) %>%
  mutate_at(.vars = 'correct', ~factor(.x, levels=c('uncorrect','correct'))) %>%
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned') %>% rev)) %>% 
  ggplot(aes(x=correct, y = freq, fill=type, stratum = type, alluvium = type)) + geom_stratum(color='NA', width = 0.3) + geom_alluvium(curve_type = 'linear', width = 0.3) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'bottom', axis.title.y = element_blank()) + scale_y_continuous(expand = c(0,0), position = 'right') + ylab('Percentage of fuison events') + geom_text(aes(label = n), position = position_stack(vjust = 0.5))  + coord_flip()+ scale_fill_manual(values = c('#DF6F4F','#7C6EA1','#8EA9D0','#E59323'))
p + canvas(5,3)
ggsave(plot = p, filename = '04error_Alignment/correct_vs_uncorrected.pdf', width = 5, height = 3)

data.frame(
  stringsAsFactors = FALSE,
  type = c("Adjacent", "Misaligned", "SHS", "Unknown"),
  correct = c(38L, 77L, 20L, 11L),
  uncorrect = c(819L, 3046L, 314L, 633L)
) %>% column_to_rownames('type') %>% 
  row_wise_prop_test()


# RNA-seq types -----------------------------------------------------------

rnaseq_all_num <- read.table('06RNAseq_results/all_types.txt', header = T)
p <- rnaseq_all_num %>% 
  pivot_longer(-type, names_to = 'tissue', values_to = 'num') %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('Readthrough','Inversion','Translocation') %>% rev)) %>%
  mutate_at(.vars = 'tissue', ~factor(.x, levels=c('flagleaf','root','leaf','panicle') %>% rev())) %>%
  ggplot(aes(x=tissue, y = num, fill=type)) + geom_bar(stat = 'identity', width = 0.6, position = position_dodge2()) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'bottom', axis.title.y = element_blank()) + ylab('Number of fusion events') + scale_fill_discrete_sequential(palette='Reds', nmax=9, order=4:6) + scale_y_continuous(position = 'right', expand = c(0,0))+ coord_flip()
p + canvas(4,3)
ggsave(plot = p, filename = '06RNAseq_results/rnaseq_type.pdf', width = 4.5, height = 3)


# Inter or intra fusion events --------------------------------------------
p <- read.table('03shs/Inter_intra.txt', header = T) %>% 
  group_by(type) %>% mutate(freq=round(num/sum(num)*100, digits = 2)) %>% 
  ggplot(aes(x=type, y=num, fill=chr)) + geom_bar(stat = 'identity', position = position_dodge2(), width = 0.6) + geom_text(aes(label = freq),position = position_dodge(width = 0.7)) + scale_y_break(c(1000,3500), scales = 0.6) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'right')  + xlab('') + ylab('Number of fusion events') + scale_fill_manual(values = c('#62cac2','#90d4c9'))
ggsave(plot = p, filename = '03shs/Inter_intra.pdf', width = 5, height = 5)

# TPM comparison ----------------------------------------------------------
## cDNA
cDNA_tpm <- read.table('03shs/tpm_cDNA.txt', header = F) %>% set_colnames(c('id1','id2','id1_exp','id2_exp','ave'))
cDNA_tpm %>% ggplot(aes(x=ave)) + geom_density() + scale_x_log10()
cDNA_tpm_random_all <- read.table('03shs/tpm_cDNA_random_allGene.txt', header = F) %>% set_colnames(c('id1','id2','id1_exp','id2_exp','ave'))
cDNA_tpm_random_all %>% ggplot(aes(x=ave)) + geom_density() + scale_x_log10()
cDNA_tpm_random_exp <- read.table('03shs/tpm_cDNA_random_expressedGene.txt', header = F) %>% set_colnames(c('id1','id2','id1_exp','id2_exp','ave'))
cDNA_tpm_random_exp %>% ggplot(aes(x=ave)) + geom_density() + scale_x_log10()

isoseq_tpm <- read.table('03shs/tpm_isoseq.txt', header = F) %>% set_colnames(c('id1','id2','id1_exp','id2_exp','ave'))
isoseq_tpm %>% ggplot(aes(x=ave)) + geom_density() + scale_x_log10()
isoseq_tpm_random_exp <- read.table('03shs/tpm_isoseq_random_expressedGene.txt', header = F) %>% set_colnames(c('id1','id2','id1_exp','id2_exp','ave'))
isoseq_tpm_random_exp %>% ggplot(aes(x=ave)) + geom_density() + scale_x_log10()

drs_tpm <- read.table('03shs/tpm_drs.txt', header = F) %>% set_colnames(c('id1','id2','id1_exp','id2_exp','ave'))
drs_tpm %>% ggplot(aes(x=ave)) + geom_density() + scale_x_log10()
drs_tpm_random_exp <- read.table('03shs/tpm_drs_random_expressedGene.txt', header = F) %>% set_colnames(c('id1','id2','id1_exp','id2_exp','ave'))
drs_tpm_random_exp %>% ggplot(aes(x=ave)) + geom_density() + scale_x_log10()


bind_rows(drs_tpm %>% add_column(type='real'), drs_tpm_random_exp %>% add_column(type='random')) %>% 
  ggplot(aes(x=ave, fill=type)) + geom_density() + scale_x_log10() + ggtitle('dRNA')
bind_rows(cDNA_tpm %>% add_column(type='real'), cDNA_tpm_random_exp %>% add_column(type='random')) %>% 
  ggplot(aes(x=ave, fill=type)) + geom_density() + scale_x_log10() + ggtitle('cDNA')
bind_rows(isoseq_tpm %>% add_column(type='real'), isoseq_tpm_random_exp %>% add_column(type='random')) %>% 
  ggplot(aes(x=ave, fill=type)) + geom_density() + scale_x_log10() + ggtitle('isoseq')

fusion_exp_obs <- bind_rows(isoseq_tpm, drs_tpm, cDNA_tpm) %>% add_column(type='Observed')
fusion_exp_random <- bind_rows(isoseq_tpm_random_exp, drs_tpm_random_exp, cDNA_tpm_random_exp) %>% add_column(type='random')


p <- bind_rows(fusion_exp_obs, fusion_exp_random) %>% 
  filter(ave!=0) %>%
  ggplot(aes(x=ave, fill=type)) + geom_density(color=NA) + scale_x_log10(label=label_number_auto()) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none') + scale_y_continuous(breaks = breaks_extended(), expand = c(0,0)) + scale_fill_manual(values = c('#B3D2E7','#81B3D6') %>% rev) + xlab('Average expression (TPM)') + ylab('Density') + geom_vline(xintercept = c(median(fusion_exp_obs$ave), median(fusion_exp_random$ave)), linetype = 2, linewidth = 0.8, colour = '#7A7A7A')
ggsave(plot = p, filename = '03shs/fusion_tpm_comparison.pdf', width = 5, height = 5)


# Overlap with other tools ------------------------------------------------
p <- data.frame(
  stringsAsFactors = FALSE,
  type = c("HC", "LC", "NS", "PRT"),
  num = c(1L, 3L, 100L, 0L)
) %>% 
  ggplot(aes(x=type, y = num)) + geom_bar(stat = 'identity', width = 0.6, fill='#6AADC1') + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + geom_text(aes(label = num))
ggsave(plot = p, filename = 'Venn/jaffal_vs_longGF_overlap.pdf', width = 3, height = 4)

# interaction upsetplot ---------------------------------------------------

interaction_matrix <- read.table('02data_quality/obversed.matrix', header = F, comment.char = '$') %>% set_colnames(c('me', 'pair','type','tech'))
a <- list(dna2dna_H3K4me3 = interaction_matrix %>% filter(me=='dna2dna_H3K4me3') %>% pull('pair'),
     dna2dna_RNAPII = interaction_matrix %>% filter(me=='dna2dna_H3K4me3') %>% pull('pair'),
     dna2rna_H3K4me3 = interaction_matrix %>% filter(me=='dna2rna_H3K4me3') %>% pull('pair'),
     dna2rna_H3K9ac = interaction_matrix %>% filter(me=='dna2rna_H3K9ac') %>% pull('pair'),
     rna2rna = interaction_matrix %>% filter(me=='rna2rna') %>% pull('pair'))
pdf(file = '02data_quality/obversed.pdf', width = 5, height = 4)
upset(fromList(a), sets = c('dna2dna_H3K4me3','dna2dna_RNAPII','dna2rna_H3K4me3','dna2rna_H3K9ac','rna2rna') %>% rev(), keep.order = T,order.by = 'degree', decreasing = F, main.bar.color = '#ABA6E6', mainbar.y.label = 'Number of fusion events', point.size = 2.5, line.size = 0.8, shade.color = '#898989', shade.alpha=0.2, matrix.dot.alpha = 0.99, matrix.color='#6FA2D7', text.scale=c(1.5,1,0.5,0.5,1.5,1))
dev.off()

# p-value
library(ggridges)
interaction_pvalue <- read.table('02data_quality/shanfeng.txt', header = F) %>% set_colnames(c('sample', 'num','type'))
interaction_num <- data.frame(
  stringsAsFactors = FALSE,
  type = c("dna2dna_H3K4me3","dna2dna_RNAPII","dna2rna_H3K4me3","dna2rna_H3K9ac","rna2rna"),
  num = c(17L, 9L, 1L, 4L, 12L))

p <- interaction_pvalue %>% 
  mutate(lnum=log2(num)) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('dna2dna_H3K4me3','dna2dna_RNAPII','dna2rna_H3K4me3','dna2rna_H3K9ac','rna2rna') %>% rev())) %>% 
  ggplot(aes(x=lnum, y=type, fill=type)) + geom_density_ridges(color='NA', scale=1, rel_min_height=0.01) + scale_y_discrete(expand = c(0,0)) + scale_fill_discrete_sequential(palette='Sunset', rev=F) + theme_pubr() + theme(legend.position = 'none') + xlab('Obversed') + ylab('') + scale_x_continuous(breaks = seq(0,8,2), labels = 2^seq(0,8,2)) + geom_text(data = interaction_num, aes(x=log2(num), y=type, label = num)) + geom_vline(data = interaction_num, aes(xintercept =log2(num)))
ggsave(plot = p, filename = '02data_quality/obversed.pvalue.pdf', width = 5, height = 2)

# adjust p-value
interaction_adjust_pvalue <- read.table('02data_quality/Interaction_adjust_p.txt', header = T)
interaction_adjust_pvalue %>% adjust_pvalue(p.col = 'p', method = 'BH')


# fused_split_com  --------------------------------------------------------

a <- read.table(text="MH63_split	great_2	1748
MH63_fused	great_2	1535
MSU7_split	great_2	772
MSU7_fused	great_2	366
MH63_split	great_3	299
MH63_fused	great_3	93
MSU7_split	great_3	33
MSU7_fused	great_3	19", header=F)
a$V1 <- factor(a$V1, levels = c('MH63_split','MH63_fused','MSU7_split','MSU7_fused'))
a$V2 <- factor(a$V2, levels = c('great_3','great_2'))
p <- a %>% ggplot(aes(x=V1, y=V3, fill=V2)) + geom_bar(stat = 'identity', width = 0.4) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), axis.text.x = element_text(angle = 45, vjust = 0.6)) + ylab('Gene number') + scale_fill_manual(values = c('#6FA2D7','#B3D2E7') %>% rev()) + geom_text(aes(label = V3))
p + canvas(4,5)
ggsave(filename = '06RNAseq_results/fusion_split.pdf', plot = p, width = 3, height = 5)

## in mean within annotated genes
b <- read.table(text = 'MH63	gene	loci	2532
MH63	isoform	transcript	2609
MH63	gene_in	loci	0
MH63	isoform_in	transcript	4180
MH63_man	gene	loci	4582
MH63_man	isoform	transcript	5631
MH63_man	gene_in	loci	0
MH63_man	isoform_in	transcript	5486
MSU7	gene	loci	11172
MSU7	isoform	transcript	11202
MSU7	gene_in	loci	0
MSU7	isoform_in	transcript	6267
IRGSP1	gene	loci	9729
IRGSP1	isoform	transcript	10125
IRGSP1	gene_in	loci	0
IRGSP1	isoform_in	transcript	3945', header=F)
p <- b %>% mutate_at(.vars = 'V3', ~factor(.x, levels=c('transcript','loci'))) %>% 
  mutate_at(.vars = 'V1', ~factor(.x, levels=c('MH63','MH63_man','MSU7','IRGSP1'))) %>%
  ggplot(aes(x=V3,y=V4, fill=V2)) + geom_bar(stat = 'identity', width = 0.5) + theme_pubr() + facet_grid(V1~., scales = 'free') + coord_flip() + theme(axis.title = element_text(colour = 'black'), axis.title.y = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black')) + ylab('Number of genes/transcripts') + scale_fill_manual(values = c("#F57F73", "#4C72B0","#B3D2E7","#81B3D6", )) + scale_y_continuous(expand = c(0,0))
ggsave(plot = p, filename = '06RNAseq_results/new_gene_isoform.pdf', width = 4, height = 4)


# Overlap fused genes -----------------------------------------------------
p <- data.frame(
  stringsAsFactors = FALSE,
       check.names = FALSE,
         tech = c("IsoSeq", "dRNA", "cDNA"),
         read1 = c(243L, 340L, 0L),
         read2 = c(259L, 680L, 629L)
) %>% 
  mutate_at(.vars = 'tech', ~factor(.x, levels=c('IsoSeq','dRNA','cDNA'))) %>%
  pivot_longer(-tech, names_to = 'reads', values_to = 'num') %>% 
  ggplot(aes(x=tech, y = num, fill = reads)) + geom_bar(stat = 'identity', position = position_dodge2(padding = 0.1), width = 0.6) + geom_text(aes(label = num), position = position_dodge2(width = 0.6)) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none', axis.title.x = element_blank()) + scale_y_continuous(expand = c(0,0)) + ylab('Fused genes') + scale_fill_manual(values = c('#60CABE','#B6E3E6'))
p + canvas(2,2)
ggsave(plot = p, filename = '06RNAseq_results/fused_supported.pdf', width = 4, height = 4)


# SY63 control sample cor -------------------------------------------------

pdf(file = '07control/sample_cor.pdf', width = 5, height = 5)
read.table('07control/sample_cor.txt', header = T, row.names = 1) %>% 
  pheatmap(border_color = 'white', display_numbers = F, number_color = 'white', show_rownames = F)
dev.off()


# trans mate pairs --------------------------------------------------------

p <- data.frame(
  stringsAsFactors = FALSE,
  sample = c("MH63_leaf", "MH63_root", "SY63_leaf", "SY63_root"),
  Intergenic = c(4532L, 4044L, 3808L, 4058L),
  Genic = c(7052L, 7994L, 7200L, 8336L)
) %>% pivot_longer(-sample, names_to = 'type', values_to = 'num') %>% 
  group_by(sample) %>% mutate(freq=round(num/sum(num)*100, digits = 2)) %>% 
  ggplot(aes(x=sample, y = num/1000, fill=type)) + geom_bar(stat = 'identity', position = position_dodge2(), width = 0.6) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none', axis.title.x = element_blank()) + scale_y_continuous(expand = c(0,0)) + geom_text(aes(label=num), position = position_dodge2(0.6))
p + canvas(5,5)
ggsave(plot = p, filename = '07control/trans_reads_num.pdf', width = 5, height = 5)


# trans gene num ----------------------------------------------------------

p <- data.frame(
  stringsAsFactors = FALSE,
            sample = c("Leaf", "Root", "Leaf", "Root"),
              type = c("Same", "Same", "Diff", "Diff"),
               num = c(2316L, 2418L, 1488L, 1264L)
) %>%
  mutate_at(.vars = 'type', ~factor(.x, levels=c('Same','Diff'))) %>%
  ggplot(aes(x=sample, y=num/1000, fill=type)) + geom_bar(stat = 'identity', position = position_dodge2(), width = 0.4) + geom_text(aes(label = num), position = position_dodge2(width = 0.4), vjust=-0.3) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none', axis.title.x = element_blank()) + scale_y_continuous(expand = c(0,0)) + ylab('Number of genes (Kb)') + scale_fill_manual(values = c('#8AD0C6', '#BAE2D9'))
p + canvas(5,5)
ggsave(plot = p, filename = '07control/trans_gene_num.pdf', width = 5, height = 5)


# cor leaf and isoseq -----------------------------------------------------

a <- read.table('07control/cor.leaf_and_isoseq', header = F) %>% set_colnames(c('id', 'leaf','isoseq'))
cor.test(a$V2, a$V3, method = 'spearman')
a %>%
  mutate(x=log10(leaf+0.1), y=log10(isoseq+0.1)) %>% 
  cor_test(x,y, method = 's')

b <- read.table('07control/cor.flagleaf_and_isoseq', header = F) %>% set_colnames(c('id', 'flagleaf','isoseq'))
b %>%
  mutate(x=log10(flagleaf+0.1), y=log10(isoseq+0.1)) %>% 
  cor_test(x,y, method = 's')
a %>% nrow()
b %>% nrow()
merge(a,b, by = 'id') %>% 
  cor_test(leaf, flagleaf)

b %>% cor_test(flagleaf, isoseq, method = 'p')


# Diff type num -----------------------------------------------------------
diff_type <- read.table('07control/Diff_type.txt', header = T)
diff_type %>% 
  group_by(tissue, type) %>% count() %>% 
  ungroup() %>% group_by(tissue) %>% mutate(freq=round(n/sum(n)*100,digits = 2)) %>% 
  write_tsv('07control/Diff_type_summary.txt')

# Same type num -----------------------------------------------------------
same_type <- read.table('07control/Same_type.txt', header = T)
p <- same_type %>% 
  group_by(tissue, type) %>% count() %>% 
  ungroup() %>% group_by(tissue) %>% mutate(freq=sprintf("%.3f",n/sum(n)*100)) %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('TS','Singleton','Control-only','Low-support','manual-review'))) %>%
  ggplot(aes(x=type, y = n, fill=type)) + geom_bar(stat = 'identity', width = 0.5) + facet_grid(.~tissue) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none', axis.title.y = element_blank()) + scale_y_continuous(expand = c(0,0)) + coord_flip() + ylab('Number of fusions') + scale_fill_manual(values = c('#7F82C7', '#B989DA','#DD9ACD','#F3B5B6','#FAD29F') %>% rev) + geom_text(aes(label=n), hjust=-0.1)
p + canvas(6,4)
ggsave(plot = p, filename = '07control/Same_type.pdf', width = 5, height = 3)


# nipp HC type number -----------------------------------------------------
library(janitor)
nipp_drs_type <- read.table('08nipp_niNanJie/nipp_drs_type.txt', header = T)
p <- nipp_drs_type %>% 
  mutate_at(.vars = 'tissue', ~factor(.x, levels=c('embryo','leaf','pistil','root','stamen','stem','seedling','floral_buds'))) %>%
  ggplot(aes(x=type, y = num, fill=type))+ geom_bar(stat = 'identity') + facet_grid(.~tissue) + scale_fill_manual(values = c('#9486A9','#8DA9FD','#F9AB85','#F0947F')) +theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top', axis.title.x = element_blank()) + scale_y_continuous(expand = c(0,0))+ ylab('Number of fusion events') + geom_text(aes(label = num), vjust=-0.3)
ggsave(plot = p, filename = '08nipp_niNanJie/nipp_drs_type.pdf', width = 11, height = 3.5)


# ZS97 fusions on ZS97 genome ---------------------------------------------
p <- read.table('10fusions_excel/ZS97_fusions.txt', header = T) %>% 
  group_by(type, fusion) %>% count() %>% 
  ungroup %>% group_by(type) %>% mutate(freq=round(n/sum(n)*100)) %>% 
  ggplot(aes(x=type, y = freq, fill=fusion)) + geom_bar(stat = 'identity', width = 0.5) + theme_pubr() + theme(axis.title = element_text(colour = 'black'), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'bottom', axis.title.x = element_blank()) + ylab('Number of fusions') + scale_fill_manual(values = c('#7F82C7','#FAD29F'))
p + canvas(4,5)
ggsave(plot = p, filename = '10fusions_excel/ZS97_fusions.pdf', width = 5, height = 6)
#

# Sanger Sequncing --------------------------------------------------------
library(sangerseqR)
sanger_559 <- readsangerseq('10fusions_excel/0003_31425052200559_(SHS02-F1R2)_[SHS02-F1].ab1')
pdf(file = '10fusions_excel/test.pdf')
chromatogram(sanger_559)
dev.off()


# Fluorescent expression --------------------------------------------------

f_exp <- read.table('f_exp.txt', header = T)
f_exp %>% tidyplot(x=type, y=exp, fill=gene) %>% 
  add_mean_bar() %>% 
  add_sem_errorbar()
p <- f_exp %>% 
  ggbarplot(x = 'type', y = 'exp', fill='gene', add = 'mean_se', position = position_dodge(0.8), color = 'NA', add.params = list(color='black', width=0.3, size=0.2)) + scale_fill_manual(values = c('#D58B94','#EFAE69','#709CC9','#786B9A')) + scale_y_continuous(expand = c(0,0), breaks = breaks_extended(n=8)) + theme(axis.title = element_text(colour = 'black'), legend.position = 'none', axis.text = element_text(color = 'black'), axis.title.x = element_blank()) + scale_y_break(c(10,20), scales = 1)
p + canvas(7,2)
ggsave(plot = p, filename = 'f_exp.pdf', width = 5, height = 2)


# Revised cDNA replicates -------------------------------------------------

revised_cDNA_all <- data.frame(
  stringsAsFactors = FALSE,
              type = c("HC-a", "HC-b", "LC", "PRT"),
          original = c(2L, 129L, 4031L, 2L),
        downsample = c(2L, 101L, 3303L, 2L))
revised_cDNA_all %>% pivot_longer(cols = -type, names_to = 'exp', values_to = 'n')
p <- revised_cDNA_all %>% 
  pivot_longer(cols = -type, names_to = 'exp', values_to = 'n') %>% 
  mutate_at(.vars = 'exp', ~factor(.x, levels=c('original','downsample'))) %>%
  ggplot(aes(x=type, y = n, fill=exp)) + geom_bar(stat = 'identity', position = position_dodge2(width = 0.8)) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + scale_y_break(c(200,3000), scales = 1) + geom_text(mapping = aes(x = type, y = n, label = n), position = position_dodge(width = 0.9)) + scale_fill_manual(values = c('#EB786B','#70B0E0'))
ggsave(plot = p, filename = 'Revised_cDNA_number_fusion_events.pdf', width = 4, height = 4)


revised_cDNA_breakpoints <- data.frame(
  stringsAsFactors = FALSE,
  type = c("Adjacent","Misaligned","SHS","Unknown"),
          original = c(207L, 352L, 2914L, 317L),
        downsample = c(208L, 304L, 2485L, 411L))
p <- revised_cDNA_breakpoints %>% 
  pivot_longer(cols = -type, names_to = 'exp', values_to = 'n') %>% 
  mutate_at(.vars = 'exp', ~factor(.x, levels=c('original','downsample'))) %>% 
  ggplot(aes(x=type, y = n, fill=exp)) + geom_bar(stat = 'identity', position = position_dodge2(width = 0.8)) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + scale_y_break(c(400,2000), scales = 1) + geom_text(mapping = aes(x = type, y = n, label = n), position = position_dodge(width = 0.9)) + scale_fill_manual(values = c('#EB786B','#70B0E0'))
ggsave(plot = p, filename = 'Revised_cDNA_number_breakpoints.pdf', width = 4, height = 4)

# Revised iso-seq correlation ---------------------------------------------

revised_isoseq_cor <- read.table('toMH63_matrix.heatmap.tab', row.names = 1, header = T)
row_anno <- data.frame( a=rownames(revised_isoseq_cor), b=rownames(revised_isoseq_cor) %>% str_split(pattern = '_', n = 2) %>% map_chr(~.x[[2]])) %>% column_to_rownames('a')
pdf(file = 'toMH63_matrix.heatmap.pdf', width = 6, height = 5)
revised_isoseq_cor %>% 
  pheatmap(border_color = 'white', display_numbers = T, number_color = 'white', show_rownames = F, annotation_row = row_anno, annotation_colors = list(b=c(flagLeaf='#9DBE47',panicle='#8B7BB8', root='#EE754C')))
dev.off()


# Revised MH63 flagleaf isoseq SMRT celles correlation --------------------

revised_isoseq_cells <- read.table('isoseq_SMRT_cells.tab', row.names = 1, header = T)
row_anno <- data.frame( a=rownames(revised_isoseq_cells), b=rownames(revised_isoseq_cells) %>% str_split(pattern = '_', n = 4) %>% map_chr(~.x[[2]])) %>% column_to_rownames('a')
pdf(file = 'isoseq_SMRT_cells.pdf', width = 6, height = 5)
revised_isoseq_cells %>% 
  pheatmap(border_color = 'white', display_numbers = T, number_color = 'white', show_rownames = F, annotation_row = row_anno, annotation_colors = list(b=c('1-2kb'='#ED665D','2-4kb'='#A8786E', '4kb+'='#729ECE')))
dev.off()
  

# Revised isoseq subreads -------------------------------------------------

revised_isoseq_subreads_fusion <- data.frame(
  stringsAsFactors = FALSE,
  type = c("HC-a","HC-b","LC","PRT"),
  num = c(8L, 546L, 3097L, 1L))
p <- revised_isoseq_subreads_fusion %>% 
  ggplot(aes(x=type, y = num, fill=type)) + geom_bar(stat = 'identity', width = 0.5) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + scale_y_break(c(800,2500)) + geom_text(mapping = aes(x = type, y = num, label = num)) + scale_fill_discrete_sequential(palette='OrYel', rev=F)
p + canvas(4,4)
ggsave(plot = p, filename = 'Revised_isoseq_subreads_number.pdf', width = 4, height = 4)


revised_isoseq_subreads_breakpoints <- data.frame(
  stringsAsFactors = FALSE,
  type = c("Adjacent","Misaligned","SHS","Unknown"),
  num = c(213L,2398L,874L,167L))
p <- revised_isoseq_subreads_breakpoints %>% 
  mutate_at(.vars = 'type', ~factor(.x, levels=c('SHS','Adjacent','Unknown','Misaligned'))) %>% 
  ggplot(aes(x=type, y = num, fill=type)) + geom_bar(stat = 'identity', width = 0.5) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + scale_y_break(c(900,2000)) + geom_text(mapping = aes(x = type, y = num, label = num)) + scale_fill_manual(values = c('#A0A1C3','#94BBC6','#90D3C1','#BCE8A8'))
ggsave(plot = p, filename = 'Revised_isoseq_subreads_breakpoints.pdf', width = 4, height = 4)
#

interaction_matrix <- read.table('02data_quality/obversed.matrix', header = F, comment.char = '$') %>% set_colnames(c('me', 'pair','type','tech'))
interaction_matrix %>% 
  group_by(tech) %>% count()


interaction_source <- read.table('obversed.plot', header = T, comment.char = '@')
interaction_source %>% 
  group_by(tech, type) %>% count()

p <- data.frame(
  stringsAsFactors = FALSE,
  tech = c("IsoSeq","IsoSeq","IsoSeq","IsoSeq","cDNAseq","cDNAseq","cDNAseq","cDNAseq","dRNAseq","dRNAseq","dRNAseq","dRNAseq"),
  type = c("HC-a","HC-b","LC","PRT","HC-a","HC-b","LC","PRT","HC-a","HC-b","LC","PRT"),
  n = c(0L, 0L, 32L, 16L, 0L, 3L, 39L, 1L, 0L, 3L, 8L, 0L))%>% 
  mutate_at(.vars = 'tech', ~factor(.x, levels=c('IsoSeq','dRNAseq','cDNAseq'))) %>% 
  ggplot(aes(x=tech, y = n, fill = type)) + geom_bar(stat = 'identity', position = position_dodge2(width = 0.8)) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Number of fusion events') + scale_y_continuous(expand = c(0,0)) + geom_text(mapping = aes(x = tech, y = n, label = n), position = position_dodge(width = 0.9), vjust = -0.3) + scale_fill_discrete_sequential(palette='OrYel', rev=F)
ggsave(plot = p, filename = 'obversed.plot.pdf', width = 4, height = 4)
#



# Revised propration of unknown junction sequences ------------------------
p <- data.frame(
  stringsAsFactors = FALSE,
              all  = c('Alignment','Alignment'),
              type = c("Partners", "Non-partners"),
             total = c(339L, 339L),
                 n = c(64L, 275L)
) %>% 
  mutate(pct=n/total*100) %>% 
  ggplot(aes(x=all, y=pct, fill=type)) + geom_bar(stat = 'identity', width = 0.3) + theme_pubr()  + theme(axis.title = element_text(colour = 'black'), axis.title.x = element_blank(), legend.title = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'top') + ylab('Proportion of junctions') + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = c('#E7988E','#9EBCD6'))
ggsave(plot = p, filename = '11_NBT_revised/unknown_to_partners.pdf', width = 3, height = 5)