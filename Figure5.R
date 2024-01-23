#---- Figure 5 ----#


# Library and data loading

require(dplyr)
require(pheatmap)
require(Seurat)
require(ggpubr)
library(ggplot2)

load('obj/ttr_obj.rda')
load('obj/cd8_tcr_obj_11.rda')
load('/raid1/cyq/tissue_col.rda')

## Cell or clone distributions
### UMAP

ttr_obj@meta.data %>%
ggplot() + 
    geom_point(aes(x=UMAP1, y=UMAP2, color=Tissue), size=.5) + 
    scale_color_manual(values = c(get_palette('lancet',15)[c(2)], 
                                  '#F0A73A', get_palette('lancet',15)[c(7)]), name='') + 
    # facet_wrap(~Tissue, nrow = 1) +
    theme_nothing() +
    theme(plot.title = element_text(size=22),
          strip.text = element_text(size = 18),
          strip.background =element_rect(fill=NA, color = NA), 
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=10, r=10)) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1), ncol=1)) 

input = ttr_obj@meta.data %>% filter(Response != 'SD') 
input %>%
ggplot() + 
    geom_point(aes(x=UMAP1, y=UMAP2, color=blood_share), size=.8) + 
    # facet_grid(cols=vars(factor(Stage, levels=c('Pre','On','Post'))),
    #            rows=vars(Response)) +
    scale_color_manual(values = c('#EEB50C','#C72D22','lightgrey', '#616FAB','#EAD8C3'), name='Blood_share') + 
    theme_pubr() + 
    theme(plot.title = element_text(size=22),
          strip.text = element_text(size = 18, face='bold'),
          strip.background =element_rect(fill=NA, color = NA, size = 1), 
          axis.text.x = element_text(size=13, angle = 0, hjust = .5, vjust=0),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .8),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=10, r=10)) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 

### histogram

selected_clone = ttr_obj@meta.data %>% 
    filter(Stage %in%c('Post')) %>%
    pull(clone.id) %>% 
    unique()
input = ttr_obj@meta.data %>% 
	filter(clone.id %in% selected_clone, Tissue=='Blood') %>%
    group_by(Patient, Response) %>% 
    mutate(n_blood=n()) %>% 
    ungroup() %>%
    group_by(Patient, Response, n_blood, clone_origin) %>% 
    summarise(n_origin=n()) %>% 
    mutate(freq_origin = n_origin/n_blood) %>% ungroup()

input %>%
    ggplot(aes(x=Patient, y=freq_origin, 
               fill=factor(clone_origin, levels=c('Pre','On','Post')))) +
    geom_bar(stat='identity', position='stack') +
    scale_fill_manual(values = c('#BB5C63','#DEB3C3','#66527E'), name='Clone origin in blood') +
    xlab('') +
    ylab('Frequency') +
    ggtitle(paste0('Post-treatment blood-asso Ttr-like clones')) +
    theme_pubr() +
    theme(plot.title = element_text(size=18),
          strip.text = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = 1.5), 
          axis.text.x = element_text(size=15, angle = 60, hjust = 1, vjust=1),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .5),
          legend.position = 'bottom',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=0, r=0))


## Expressions
### Heatmap

blood_t_clone = cd8_tcr_obj@meta.data %>% 
	filter(Tissue=='Blood') %>% 
	pull(clone.id) %>% 
	unique()

cd8_tcr_obj@meta.data = cd8_tcr_obj@meta.data %>% 
    mutate(tissue_share = case_when(clone.id %in% blood_t_clone ~ 'Blood_shared',
                                    TRUE ~ 'Others'))
# table(cd8_tcr_obj$tissue_share)
cd8_tcr_obj@meta.data = cd8_tcr_obj@meta.data %>% 
	mutate(tissue_ttr = case_when(Tissue=='Blood' & Ttr=='Non-Ttr' ~ 'Blood_Non-Ttr',
		Tissue=='Blood' & Ttr=='Ttr' ~ 'Blood_Ttr',
		Ttr=='Non-Ttr' ~ 'Tissue_Non-Ttr',
		Ttr=='Ttr' & tissue_share=='Blood_shared' ~ 'Tissue_Ttr_blood_asso',
		Ttr=='Ttr' & tissue_share!='Blood_shared' ~ 'Tissue_Ttr_blood_unasso',
		TRUE ~ 'Others'))
# table(cd8_tcr_obj@meta.data$tissue_ttr, cd8_tcr_obj@meta.data$Ttr)

genes=c('CXCL13','HAVCR2','ENTPD1','ITGAE','LAYN','PDCD1','IFNG','IL2RA','TNFRSF9', 
        'HLA-DQB1','HLA-DRA','HLA-DQA1','CD74','CCL4','CCL4L2','GZMK')
orders = c('Blood_Non-Ttr','Blood_Ttr','Tissue_Ttr_blood_asso','Tissue_Ttr_blood_unasso')
input = data.frame(
    tissue = cd8_tcr_obj$Tissue,
    tissue_share = cd8_tcr_obj$tissue_ttr,
    t(as.matrix(cd8_tcr_obj@assays$RNA@data[genes,]))
) %>% 
    filter(tissue != 'Normal') %>%
    select(-tissue) %>%
    filter(tissue_share %in% orders) %>%
    group_by(tissue_share) %>% 
    summarise_all(list(mean)) %>% 
    tibble::column_to_rownames('tissue_share') %>%
    scale()
input = input[orders,]
# head(input)

pheatmap(t(input), cluster_cols = F, cluster_rows = F, display_numbers = F,
#          filename = 'figures/fig5_Ttr_tissue_share_heatmap.pdf', 
         height = 7, width = 3 ,fontsize = 12, scale = 'none',
        color = rev(brewer.pal(11, 'RdYlBu')), 
         border_color = NA)

## Proportion

input0 = ttr_obj@meta.data %>% 
	filter(Response != 'SD', Tissue!='Blood') %>%
	mutate(pt_rsp = paste0(Patient, '-', Response))

input1 = input0 %>% 
    group_by(Patient, Response, Stage, tissue_share) %>% 
    summarise(n=n()) %>% 
    ungroup() %>%
    group_by(Patient, Stage, Response) %>% 
    mutate(n0 = sum(n))

tmp = input0 %>% 
	tidyr::expand(pt_rsp, Stage, tissue_share) %>% 
	separate(pt_rsp, into = c('Patient','Response'), sep='-') %>% 
	distinct() %>% 
    left_join(input1, by=c('Patient','Response', 'Stage','tissue_share')) %>%
    mutate(prop = n / n0) %>%
    replace_na(list(prop=0))

### boxplot

input = tmp %>% filter(tissue_share=='Blood_associated')

input %>% 
ggplot(aes(col=Response, y=prop, x=Response)) + 
    geom_boxplot(width=.5) +
    scale_shape_manual(values = c(17,15,16)) +
    geom_jitter(aes(shape=Stage), width=.25, size=1.5) +
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
          plot.margin = margin(b=0, l=10, t=10, r=10)) +
    scale_color_manual(values = c('#616FAB','#EEB50C') , name='Response')+ 
    guides(shape = guide_legend(override.aes = list(size=5, alpha = 1)))


### Barplot

input = tmp %>% 
    select(-n,-n0) %>%
    group_by(Stage, Response,tissue_share) %>%
    summarise(prop_mean = mean(prop), prop_median=median(prop)) %>%
    select(-prop_median) %>%
    pivot_wider(names_from = 'Response', values_from = 'prop_mean') %>%
    mutate(log2FC = log2((CR+0.001)/(PR+0.001)))
# head(input)

input %>% ggplot(aes(y=log2FC, x=tissue_share)) + 
    stat_summary(fun = mean, geom = 'bar') +
    facet_wrap(~factor(Stage, levels=c('Pre','On','Post')), nrow=1) +
    geom_hline(yintercept = 0, linetype='dashed') +
    xlab('') + ylab('log2FC(CR/PR)') +
    theme_pubr() + 
    theme(strip.text = element_text(size = 13, face='bold'),
          strip.background =element_rect(fill=NA, color = NA, size = 1), 
          axis.text.x = element_text(size=13, angle = 60, hjust = 1, vjust=1),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .6),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=10, r=10))