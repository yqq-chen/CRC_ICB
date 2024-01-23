#---- Figure 4 ----#


# Library and data loading

library(data.table)
require(dplyr)
require(Seurat)
require(ggpubr)
library(ggplot2)

load('obj/ttr_obj.rda')
load('obj/tcr_obj.rda')
sample_info = data.frame(fread('obj/clinical_info_09.csv'))


## UMAP

ttr_obj@meta.data %>% 
ggplot() + 
    geom_point(aes(x=UMAP1, y=UMAP2, color=Ttr_sub), size=2, shape=16, stroke=0) + 
    scale_color_manual(values = tex_col, name='Annotation') + 
    # facet_grid(cols=vars(factor(Stage, levels=c('Pre','On','Post'))),
    #            rows=vars(Response)) +
	theme_pubr() +
	theme(plot.title = element_text(size=22),
	      strip.text.x = element_text(size = 18),
	      strip.background =element_rect(fill='white', color = 'black', size = 1), 
	      axis.text.x = element_text(size=18, angle = 0, hjust = .5, vjust=0),
	      axis.text.y = element_text(size=15),
	      axis.title = element_text(size = 18),
	      axis.line = element_line(size = .8),
	      legend.position = 'bottom',
	      legend.key.size = unit(1.4,'line'),
	      legend.text = element_text(size = 15),
	      legend.title = element_text(size = 15, face = 'bold'),
	      plot.margin = margin(b=20, l=10, t=10, r=10)) +
	guides(colour = guide_legend(override.aes = list(size=5, alpha = 1), ncol=2)) 



## Tumor-reactive CD8 T cell subtypes over Tumor-reactive CD8 T cell population
### Calculating abundance matrix

tmp1 = ttr_obj@meta.data %>% 
	group_by(Patient, Stage) %>% 
	summarise(n1=n())
tmp2 = ttr_obj@meta.data %>% 
	group_by(Patient, Stage, Ttr_sub) %>% 
	summarise(n2=n())
ttr_freq_df = tmp2 %>% 
	left_join(tmp1, by = intersect(colnames(tmp1),colnames(tmp2))) %>%
    mutate(freq=n2/n1) %>% 
    left_join(sample_info, by='Patient')

### Plotting
tmp %>%
	filter(Response != 'SD') %>%
	ggplot(aes(x=factor(Stage, levels=c('Pre','On','Post')), y=freq, 
	           fill=factor(Stage, levels=c('Pre','On','Post')))) + 
	geom_boxplot(aes(col=factor(Stage, levels=c('Pre','On','Post'))),
	             alpha=0.3, position=position_dodge(width=.8), width=.7, lwd=.8, outlier.color = NA) +
	geom_point(aes(col=factor(Stage, levels=c('Pre','On','Post'))),
	           position = position_jitterdodge(jitter.width =.3, dodge.width = .8)) +
	scale_fill_manual(values = c('#BB5C63','#DEB3C3','#66527E'), name='Stage') +
	scale_color_manual(values = c('#BB5C63','#DEB3C3','#66527E'), name='') + 
	facet_grid(cols = vars(Ttr_sub), rows = vars(Response), scales = 'free') +
	xlab('') + ylab('Proportion (% Ttr cells)') +
	# ylim(0,1) +
	theme_pubr() +
	theme(strip.text = element_text(size = 13, face='bold'),
	      strip.background =element_rect(fill=NA, color = NA, size = 1), 
	      axis.text.x = element_text(size=13, angle = 0, hjust = .5, vjust=.5),
	      axis.text.y = element_text(size=13),
	      axis.title = element_text(size = 18),
	      axis.line = element_line(size = .6),
	      legend.position = 'right',
	      legend.key.size = unit(1.4,'line'),
	      legend.text = element_text(size = 15),
	      legend.title = element_text(size = 15, face = 'bold'),
	      plot.margin = margin(b=0, l=10, t=10, r=10))



## Tumor-reactive CD8 T cell subtypes over T cell population in tumor/blood
### Calculating abundance matrix

ttr_global_df = tcr_obj@meta.data %>% 
    group_by(Tissue, Patient, Stage) %>% 
    summarise(n1 = n()) 

ttr_sub_df = tcr_obj@meta.data %>% 
    group_by(Tissue, Ttr_sub, Patient, Stage) %>%
    summarise(n = n()) %>% 
    replace_na(list(n=0))

ttr_freq_df = tcr_obj@meta.data %>% 
	mutate(batch = paste(Tissue, Patient, Stage, sep = '-')) %>% 
    tidyr::expand(batch, Ttr_sub) %>% 
    separate(batch, into = c('Tissue','Patient','Stage'), sep = '-') %>%
    left_join(ttr_sub_df, by = c('Tissue', 'Patient', 'Stage','Ttr_sub')) %>% 
    left_join(ttr_global_df, by = c('Tissue', 'Patient', 'Stage')) %>% 
    replace_na(list(n=0)) %>%
    mutate(freq = n/n1) %>% 
    left_join(sample_info[,c('Patient','Response','MSI','SIZE_CHANGE')], by = 'Patient')

ttr_freq_df$Stage = factor(ttr_freq_df$Stage, levels = c('Pre','On','Post'))

tmp = ttr_freq_df %>% 
	group_by(Tissue, Ttr_sub, Patient) %>% 
	summarise(freq_dnmc = freq[which(Stage=='Post')] - freq[which(Stage=='Pre')])
ttr_freq_df = ttr_freq_df %>% left_join(tmp, by = c('Tissue','Ttr_sub','Patient'))


### Plotting

ttr_freq_df %>% 
	filter(Tissue == 'Tumor') %>%
	filter(Response %in% c('CR','PR')) %>%
	ggplot(aes(x=factor(Stage, levels=c('Pre','On','Post')), y=freq, 
	           fill=factor(Stage, levels=c('Pre','On','Post')))) + 
	geom_boxplot(aes(col=factor(Stage, levels=c('Pre','On','Post'))),
	             alpha=0.3, position=position_dodge(width=.8), width=.7, lwd=.8, outlier.color = NA) +
	geom_point(aes(col=factor(Stage, levels=c('Pre','On','Post'))),
	           position = position_jitterdodge(jitter.width =.3, dodge.width = .8)) +
	scale_fill_manual(values = c('#BB5C63','#DEB3C3','#66527E'), name='Stage') +
	scale_color_manual(values = c('#BB5C63','#DEB3C3','#66527E'), name='') + 
	facet_grid(cols = vars(Ttr_sub), rows = vars(Response), scales = 'free') +
	xlab('') + ylab('Proportion (%T cells)') +
	theme_pubr() +
	theme(strip.text = element_text(size = 13, face='bold'),
	      strip.background =element_rect(fill=NA, color = NA, size = 1), 
	      axis.text.x = element_text(size=13, angle = 0, hjust = .5, vjust=.5),
	      axis.text.y = element_text(size=13),
	      axis.title = element_text(size = 18),
	      axis.line = element_line(size = .6),
	      legend.position = 'right',
	      legend.key.size = unit(1.4,'line'),
	      legend.text = element_text(size = 15),
	      legend.title = element_text(size = 15, face = 'bold'),
	      plot.margin = margin(b=0, l=10, t=10, r=10))

## Pie chart
### Ttr subtypes at single clone level

clone_df = ttr_obj@meta.data %>% 
    group_by(clone.id, Treatment, Ttr_sub, Patient, Response) %>% 
    summarise(n=n()) %>%
    pivot_wider(names_from = 'Ttr_sub', values_from = c('n')) %>%
    arrange(factor(clone.id, levels=clone_order)) %>% 
    ungroup() %>%
    group_by(clone.id) %>% 
    mutate(n=n()) %>% 
    filter(n>1)

input = clone_df %>% ungroup() %>%
    pivot_longer(cols = contains('_'), names_to = 'Ttr_sub', values_to = 'cell_counts') %>%
    replace_na(list(cell_counts=0)) %>%
    group_by(clone.id, Patient, Response) %>% 
    mutate(n_clone=sum(cell_counts)) %>% 
    ungroup() %>%
    arrange(desc(Ttr_sub)) %>%
    group_by(clone.id, Treatment, Patient, Response) %>% 
    mutate(n_treatment=sum(cell_counts),
    	prop = cell_counts / n_treatment *100) %>% 
    mutate(ypos = cumsum(prop)- 0.5*prop ) %>% 
    ungroup()
    
clone_order = input %>% 
	arrange(desc(n_clone)) %>% 
	pull(clone.id) %>% 
	unique()

### Plotting

input %>% 
    filter(clone.id == 'clone:31') %>%
    ggplot(aes(x=1, y=prop, fill=Ttr_sub)) + 
    geom_bar(stat="identity", width=1, color="white", position='fill') +
    coord_polar("y", start=0) +
    facet_grid(cols=vars(factor(Treatment, levels=c('I','II','III','IV'))),
               rows=vars(factor(clone.id, levels=clone_order))) +
    geom_text(aes(label=Patient), x=1, y=0) +
    theme_void() + 
    theme(legend.position="right",
          strip.text = element_text(size = 13,face='bold'),
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13, face = 'bold')) +
    scale_fill_manual(values = tex_col , name='Ttr sub-types')