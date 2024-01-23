#---- Figure 1 ----#


## Library and data loading

library(data.table)
library(dplyr)
library(tidyr)
library(Seurat)
library(ggpubr)
library(ggplot2)

sample_info = data.frame(fread('obj/clinical_info_09.csv'))
df = read.table('obj/df_0618.txt')
load('/raid1/cyq/major_celltype_col.rda')
load('/raid1/cyq/subtype_col.rda')
# dim(df)

## UMAP
### Cell subtypes

df %>% 
ggplot(aes(x=UMAP1, y=UMAP2, color=annotation)) + 
    geom_point(size=.1, shape=16, stroke=0) + 
    theme_pubr() + 
    scale_color_manual(values = subtype_col, name='') + 
    theme(legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          legend.position = 'none')+
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))

### Major cell types

df %>% 
ggplot(aes(x=UMAP1, y=UMAP2, color=MajorCellType)) + 
    geom_point(size=.1, shape=16, stroke=0) + 
    theme_pubr() + 
    scale_color_manual(values = major_cell_colors, name='') + 
    theme(legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          legend.position = 'right')+
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


## Proportion of major celltypes over all immune cells in tumor
### Calculating abundance matrix

imm_cls = c('T','B','Mye','ILC')
meta = df %>% filter(Tissue == 'Tumor')
anno_df = data.frame()
global_df = meta %>% 
	filter(MajorCellType %in% imm_cls) %>%
    group_by(Tissue,Patient, Stage) %>% 
    summarise(n = n())

sub_df = meta %>% 
	filter(MajorCellType %in% imm_cls) %>%
    group_by(MajorCellType, Tissue, Patient, Stage) %>% 
    summarise(n1 = n())

anno_df= meta %>% 
	filter(MajorCellType %in% imm_cls) %>%
    tidyr::expand(MajorCellType, Batch, Tissue) %>%
    separate(Batch, into=c('Patient', 'Stage'), sep='-') %>%
    select(MajorCellType, Patient, Stage, Tissue) %>% 
    left_join(sub_df, by=c('MajorCellType', 'Patient', 'Stage','Tissue')) %>% 
    left_join(global_df, by = c('Patient', 'Stage','Tissue')) %>% 
    replace_na(list(n1=0,n=0)) %>%
    mutate(freq = n1/n) 

tmp = anno_df %>% 
	group_by(MajorCellType, Patient) %>% 
	summarise(freq_dnmc = freq[which(Stage=='Post')] - freq[which(Stage=='Pre')])
anno_df = anno_df %>% left_join(tmp, by = c('MajorCellType','Patient'))
# head(anno_df)

### Plotting
#### Pre-treatment proportion

input = anno_df %>% 
    select(-c(n1,n)) %>% 
    distinct() %>% 
    filter(Stage == 'Pre') 

input %>%
ggplot(aes(x=Response, y=freq)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response, shape=MSI), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    scale_y_continuous(limits=c(0,1)) +
    scale_color_manual(values = c('#5778A2','#DCB44B','#8A281E'), name='Response') +
    facet_wrap(~MajorCellType, ncol=4, scales = 'free_y') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10))

#### Proportion dynamics

input = anno_df %>% 
    select(-c(n1,n)) %>% 
    distinct() %>% 
    filter(Stage == 'Pre') 

input %>%
ggplot(aes(x=Response, y=freq_dnmc)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = 1, width=.5, fill=NA) +
    geom_point(aes(color=Response, shape=MSI), size=1.5, stroke=1, position = position_jitter(width=.15)) +
    theme_pubr() + 
    scale_color_manual(values = c('#5778A2','#DCB44B','#8A281E'), name='Response') +
    facet_wrap(~MajorCellType, ncol=4, scales = 'free_y') +
    theme(strip.text.x = element_text(size = 15),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 15),
          axis.line = element_line(size = .5),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=30, l=10)) 



