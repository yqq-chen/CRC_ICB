#---- Figure 3 ----#


# Library and data loading
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
require(Seurat)
require(ggpubr)

sample_info = data.frame(fread('obj/clinical_info_09.csv'))
load('/raid1/cyq/tumor_t_obj.rda')
load('/raid1/cyq/t_obj.rda')

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

## UMAP

tumor_t_obj@meta.data %>% 
ggplot() + 
    geom_point(aes(x=UMAP1, y=UMAP2, color=annotation), size=.01, alpha=1) + 
    scale_color_manual(values = col_yy, name='T sub-types') + 
#     geom_text(aes(x=m1, y=m2, label=anno_brief), data = text_anno, size=5) +
    theme_pubr() +
    theme(plot.title = element_text(size=22),
          strip.text = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = .1), 
          axis.text.x = element_text(size=18, angle = 0, hjust = .5, vjust=0),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .8),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=10, r=10)) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


tumor_t_obj@meta.data %>% 
ggplot(aes(x=UMAP1, y=UMAP2)) + 
    facet_grid(cols=vars(factor(Stage, levels=c('Pre','On','Post'))),
               rows=vars(Response)) +
    stat_density_2d(aes(fill = stat(piece)),
                  geom = "polygon", 
                  n = 30, 
                  bins = 10,
                  contour = T,
                  contour_var = "ndensity", adjust = 1/4
                   ) + 
    scale_fill_viridis(option = "A") +
    theme_pubr() +
    theme(plot.title = element_text(size=22),
          strip.text = element_text(size = 18),
          strip.background =element_rect(fill=NA, color = NA, size = .1), 
          axis.text = element_blank(),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .8),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=10, r=10))

## Proportion

global_df = t_obj@meta.data %>% 
            group_by(Tissue, Patient, Stage) %>% summarise(n1 = n()) 

sub_df = t_obj@meta.data %>% 
         group_by(Tissue, annotation, Patient, Stage) %>% 
         summarise(n = n())


sub_freq_df =
    t_obj@meta.data %>% 
    mutate(sample_stage = paste(Patient, Stage, Tissue, sep = '-')) %>% 
    tidyr::expand(sample_stage,annotation) %>%
    separate(sample_stage, into = c('Patient','Stage','Tissue'), sep = '-') %>%
    left_join(sub_df, by = c('Tissue', 'Patient', 'Stage','annotation')) %>% 
    left_join(global_df, by = c('Tissue', 'Patient', 'Stage')) %>%
    replace_na(list(n=0, n1=0)) %>%
    mutate(freq = n/n1) %>% 
    left_join(sample_info, by = 'Patient')

sub_freq_df$Stage = factor(sub_freq_df$Stage, levels = c('Pre','On','Post'))

# head(sub_freq_df)
# dim(sub_freq_df)

### Scatter plot

cls = 'c23_CD8_Tex_LAYN'
input = sub_freq_df %>% 
    filter(Tissue == 'Tumor') %>%
    select(-c(n1,n)) %>% distinct() %>% 
    pivot_wider(names_from = Stage, values_from = freq) %>%
    replace_na(list(On=0, Post=0, Pre=0)) %>%
    filter(annotation %in% cls) 

lm = lm(Pre ~ SIZE_CHANGE, data = input)
print(paste0('R^2 = ',format(summary(lm)$r.squared, digits = 2)))
print(paste0('p.value = ',format(lmp(lm), digits = 2)))

input %>%
    ggplot(aes(x = Pre, y = SIZE_CHANGE)) +
    geom_smooth(method='lm', col='black', lwd=.5) +
    geom_point(aes(col = Response, shape = MSI), size = 2, stroke=2) +
    scale_color_manual(values = c('#5778A2','#DCB44B','#8A281E')) +
    scale_shape_manual(values = c(16,17,1,4),name='MSI') +
    theme_pubr() + labs(title = cls) +
    geom_text(x = .05, y = 1.1, label = paste0('R = ',format(sqrt(summary(lm)$r.squared), digits = 2))) +
    geom_text(x = .045, y = 1.05, label = paste0('P = ',format(lmp(lm), digits = 2))) +
    geom_text(aes(x=Pre, y=SIZE_CHANGE-0.08, label=Patient), data=input) +
    theme(legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          axis.title = element_text(size=15),
          plot.margin = margin(r = 10, t=10)) +
    xlab('Frequency at baseline') +
    ylab('Change of tumor size (%)')


### Line chart

input = sub_freq_df %>% 
	filter(Tissue == 'Tumor') %>%
	filter(annotation==cls) 
input %>%
	ggplot(aes(x=Stage, y=freq, color=Response, shape=MSI, group=Patient)) + 
	geom_point(aes(shape=MSI), size=2, stroke=2) +
	geom_line(aes(linewidth=MSI), alpha=.5) + 
	scale_color_manual(values = get_palette('nejm', 3)[c(2,3,1)], name='Response') +
	scale_shape_manual(values = c(16,17,1,4), name='MSI state') +
	scale_linewidth_discrete(range = c(.5,1.5), name='MSI state') +
	xlab('Sampling stage') + ylab(paste0('Proportion')) +
	facet_wrap(~Response) +
	theme_pubr() +
	theme(strip.text.x = element_text(size = 18),
	      strip.background =element_rect(fill=NA, color = NA, size = 1), 
	      axis.text.x = element_text(size=15, angle = 0, hjust = .5),
	      axis.text.y = element_text(size=15),
	      axis.title = element_text(size = 18),
	      axis.line = element_line(size = .8),
	      legend.position = 'right',
	      legend.key.size = unit(1.4,'line'),
	      legend.text = element_text(size = 13),
	      legend.title = element_text(size = 15, face = 'bold'),
	      plot.margin = margin(b = 10, t=10, l=10, r=10))


## Expansion 

meta_pre = t_obj@meta.data %>% filter(Tissue=='Tumor', Stage=='Pre')
in.dat = meta_pre[,c('cell_name', 'clone.id', 'clone.status', 'Patient', 'annotation', 'Tissue')]
colnames(in.dat) = c('Cell_Name', 'clone.id', 'clone.status', 'patient', 'majorCluster', 'loc')
out_pre = Startrac.run(in.dat, proj="CRC",verbose=F) %>% suppressWarnings()
expa_pre = out_pre@cluster.data %>% filter(aid != 'CRC') %>% select('aid','majorCluster','expa') %>%
               rename(Patient=aid, annotation=majorCluster, exp_pre=expa) %>%
               mutate_all(~replace(., is.nan(.), 0))


expa_df = expa_pre %>% full_join(expa_on, by=c('Patient','annotation')) %>% 
                       full_join(expa_post, by=c('Patient','annotation')) %>%
                       left_join(sample_info, by='Patient')
# head(expa_df)

### scatter plot

input = expa_df %>% rename(Pre=exp_pre, On=exp_on, Post=exp_post) %>%
    filter(annotation == 'c23_CD8_Tex_LAYN')
lm = lm(Pre ~ SIZE_CHANGE, data = input)
# print(paste0('R^2 = ',format(summary(lm)$r.squared, digits = 2)))
# print(paste0('p.value = ',format(lmp(lm), digits = 2)))

input %>%
ggplot(aes(x = Pre, y = SIZE_CHANGE)) +
    geom_smooth(method='lm', col='#525252') +
    geom_point(aes(col = Response, shape = MSI), size = 2, stroke = 2) +
    scale_color_manual(values = c('#5778A2','#DCB44B','#8A281E')) +
    scale_shape_manual(values = c(16,17,1,4), name='MSI') +
    geom_text(x = .032, y = 1.3, label = paste0('R = ',format(sqrt(summary(lm)$r.squared), digits = 2))) +
    geom_text(x = .036, y = 1.2, label = paste0('P = ',format(lmp(lm), digits = 2))) +
    geom_text(aes(x=Pre, y=SIZE_CHANGE-0.08, label=Patient), data=input) +
    theme_pubr() + 
    theme(legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          axis.title = element_text(size=15),
          plot.margin = margin(l = 10, b=10, t=10)) +
    labs(title = paste0('c23_CD8_Tex_LAYN at baseline')) +
    ylab('Tumor Regression') +
    xlab('STARTRAC-expa')


### line chart

input=expa_df %>% 
rename(Pre=exp_pre, On=exp_on, Post=exp_post) %>%
pivot_longer(cols = c('Pre','On','Post'),
    names_to = "Stage",
    values_to = "exp",
    values_drop_na = TRUE) %>% 
filter(annotation=='c23_CD8_Tex_LAYN') 

input %>%
ggplot(aes(x=factor(Stage, levels=c('Pre','On','Post')), y=exp,
           group=Patient, col=Response)) + 
geom_point(aes(shape=MSI), size=2, stroke=2) +
geom_line(aes(linewidth=MSI), alpha=.5) + 
scale_color_manual(values = get_palette('nejm', 3)[c(2,3,1)], name='Response') +
scale_shape_manual(values = c(16,17,1,4), name='MSI state') +
scale_linewidth_discrete(range = c(.5,1.5), name='MSI state') +
xlab('') + ylab(paste0('STARTRAC-expa')) +
facet_wrap(~Response) +
theme_pubr() +
theme(strip.text.x = element_text(size = 18),
      strip.background =element_rect(fill=NA, color = NA, size = 1), 
      axis.text.x = element_text(size=15, angle = 0, hjust = .5),
      axis.text.y = element_text(size=15),
      axis.title = element_text(size = 18),
      axis.line = element_line(size = .8),
      legend.position = 'right',
      legend.key.size = unit(1.4,'line'),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15, face = 'bold'),
      plot.margin = margin(b = 10, t=10, l=10, r=10)) +
guides(colour = guide_legend(override.aes = list(size=3, alpha = 1),ncol=1),
       shape=guide_legend(ncol=1))
