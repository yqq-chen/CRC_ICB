#---- Figure 2 ----#


# Library and data loading

library(tidyr)
library(data.table)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
require(ggpubr)
library(Seurat)
library(NMF)
library(ggradar)
library(corrr)

load('/raid1/cyq/subtype_col.rda')
load('/raid1/cyq/proportion_subtype.rda')
load('raid1/cyq/tumor_obj.rda')
load('raid1/cyq/score_module.rda')
sample_info = data.frame(fread('obj/clinical_info_09.csv'))

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

lm_reg_coeff = function(data, group_name, y, x){
    data %>% 
        tidyr::nest(-!!sym(group_name)) %>%
        mutate(fit = purrr::map(data, ~lm(!!sym(y) ~ !!sym(x), data=.))) %>%
        mutate(glanced = purrr::map(fit, broom::glance), 
               tided = purrr::map(fit, broom::tidy)) %>%
        select(!!sym(group_name), glanced, tided) %>%
        tidyr::unnest(glanced) %>% 
        select(!!sym(group_name), r.squared, tided) %>%
        tidyr::unnest(tided) %>% 
        filter(term == x) %>%
        select(!!sym(group_name), r.squared, estimate, p.value) %>%
        mutate(Ti = sign(estimate) * r.squared)
}

## Scatter plot

df_reg_dnmc = anno_df_t %>% 
    select(annotation, SIZE_CHANGE, freq_dnmc) %>%
    distinct() %>% 
    lm_reg_coeff('annotation', 'SIZE_CHANGE','freq_dnmc') %>% 
    as.data.frame() %>% 
    select(annotation, Ti, p.value)%>% 
    rename(pval_Ti = p.value)
# head(df_reg_dnmc)

Roe = with(df %>% filter(Tissue %in% c('Normal', 'Tumor')),
    chisq.test(annotation, Tissue))
roe_pearson_residual_df = Roe$residuals %>% 
    as.data.frame %>% 
    pivot_wider(names_from = "Tissue", values_from = "Freq")
roe_ti_input = df_reg_dnmc %>% 
    left_join(roe_pearson_residual_df, by='annotation') %>% 
    rename(SubCellType=annotation)

lm_input = lm(Ti ~ Tumor, data=roe_ti_input)
inproe_ti_inputut %>%
    ggplot(aes(x=Tumor, y=Ti)) +
    geom_smooth(method='lm', col='black', lwd=.5) +
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    geom_point(aes(col=SubCellType)) + #, size=-log10(fdr_ti))) +
    scale_color_manual(values = subtype_col, name = 'Cell subtype') +
    # scale_size(range = c(1, 8), name='-log10(FDR_Ti)') +
    geom_text(x = 50, y = .45, label = paste0('R = ',format(sqrt(summary(lm_input)$r.squared), digits = 2)), size=4) +
    geom_text(x = 50, y = .4, label = paste0('P = ',format(lmp(lm_input), digits = 2)), size=4) +
    labs(x="Tumor enrichment", y="Therapeutic Index") +
    theme_pubr() +
    theme(plot.title = element_text(face = "bold", size=16),
          axis.text.x = element_text(size=15, angle = 0, hjust = .5),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .8),
          legend.position = 'right',
          legend.key.size = unit(2,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t = 0, l=0) 
         ) +
    guides(color = 'none') 



## NMF program
# Constructing cellular abundance matrix

meta = tumor_obj@meta.data %>% 
  mutate(pt_stg = paste(Patient, Stage, sep='-')) 

sub_freq_df = data.frame()
for (i in MajorCellTypes) {
    global_df = meta %>% filter(MajorCellType == i) %>% 
        group_by(pt_stg) %>% summarise(n = n())

    sub_df = meta %>% filter(MajorCellType == i) %>% 
        group_by(annotation, pt_stg) %>% summarise(n1 = n())

    tmp = meta %>% filter(MajorCellType == i) %>% 
        tidyr::expand(annotation, pt_stg) %>%
        left_join(sub_df, by=c('annotation', 'pt_stg')) %>% 
        left_join(global_df, by = c('pt_stg')) %>% replace_na(list(n1=0,n=0)) %>%
        mutate(freq = n1/n) 
    sub_freq_df = rbind(sub_freq_df, tmp)
}

nmf_input = sub_freq_df %>% 
	  mutate(pt_stg_2 = pt_stg) %>% separate(pt_stg_2, into = c('Patient','Stage'), sep = '-') %>% 
	  left_join(sample_info[, c('Patient','Response','MSI','SIZE_CHANGE')], by = 'Patient') %>%
	  mutate(rsp_stg = paste(Response, Stage, sep = '-')) %>%
    tidyr::pivot_wider(
        id_cols = c('pt_stg','Response','Stage','MSI','rsp_stg'),
        names_from = annotation, 
        values_from = freq) %>% 
	  arrange(Response, Stage) 


# NMF analysis

nrank=5
mtx = nmf_input[, 6:ncol(nmf_input)]
nmf.res = nmf(t(nmf_input[, 6:ncol(nmf_input)]), rank=nrank, nrun=20, seed=123456)
w = basis(nmf.res)
rownames(w) = colnames(mtx)
colnames(w) = paste0('NMF', seq(1,nrank))

w_ = w %>% t %>% scale
p=ComplexHeatmap::Heatmap(w_, width = 11, height = 4, column_km = nrank, 
    cluster_rows = F, clustering_method_columns = 'single') 
heatmap = draw(p)


# Activity dynamics

h = coef(nmf.res) %>% t
colnames(h) = paste0('NMF', seq(1,nrank))

h_df = h %>% as.data.frame() %>% 
    mutate(rsp_stg = nmf_input$rsp_stg, 
    response = nmf_input$Response,
    pt_stg = nmf_input$pt_stg)

input = h_df %>% separate(pt_stg, into = c('Patient','Stage'), sep = '-') %>%
    select(-rsp_stg) %>%
    pivot_longer(names_to = 'NMF', values_to = 'value', cols = contains('NMF')) %>%
    pivot_wider(names_from = 'Stage', values_from = 'value') %>%
    mutate(dnmc=Post-Pre)

input %>%
ggplot(aes(x=response, y=dnmc)) +
    geom_boxplot(aes(col=response),alpha=.5, width=.6, lwd=.8, outlier.color = NA) +
    geom_point(aes(col=response), position = position_jitter(width = .15)) +
    geom_hline(yintercept = 0, linetype='dashed') +
    facet_wrap(~NMF, nrow=1) +
    scale_color_manual(values = c('#5778A2','#DCB44B','#8A281E'), name='Response') +
    theme_pubr() + ylab('NMF score dynamic') + xlab('') +
    theme(plot.title = element_text(face = "bold", size=16),
          strip.text = element_text(size = 15, color='black'),
          strip.background =element_rect(fill=NA, color = NA, size = .5), 
          axis.text.x = element_text(size=15, angle = 0, hjust = .5, vjust=.5),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .5),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t = 0, l=0))

# Radar

radar_info = data.frame()
for (r in c('CR','PR','SD')) {
    if (r=='SD'){
        for (s in c('Pre','Post')) {
            radar_info_ = score_module_df %>% #t() %>% as.data.frame() %>% 
#                 tibble::rownames_to_column('sample') %>%
                mutate(pt_stg = sample) %>% separate(pt_stg, into = c('Patient','Stage'), sep='-') %>%
                left_join(sample_info[,c('Patient','Response')], by='Patient') %>%
                filter(Stage==s,Response==r) 
            radar_info_label_ = radar_info_ %>% select('sample','Patient','Stage','Response')
            radar_info_ = radar_info_ %>% select(-c('sample','Patient','Stage','Response'))
            radar_info_ = cbind(radar_info_label_, ((radar_info_-min(radar_info_))/(max(radar_info_)-min(radar_info_))))
            radar_info = rbind(radar_info, radar_info_)
        }
    }
    else {
        for (s in c('Pre','On','Post')) {
            radar_info_ = score_module_df %>% #t() %>% as.data.frame() %>% 
#                 tibble::rownames_to_column('sample') %>%
                mutate(pt_stg = sample) %>% separate(pt_stg, into = c('Patient','Stage'), sep='-') %>%
                left_join(sample_info[,c('Patient','Response')], by='Patient') %>%
                filter(Stage==s,Response==r) 
            radar_info_label_ = radar_info_ %>% select('sample','Patient','Stage','Response')
            radar_info_ = radar_info_ %>% select(-c('sample','Patient','Stage','Response'))
            radar_info_ = cbind(radar_info_label_, ((radar_info_-min(radar_info_))/(max(radar_info_)-min(radar_info_))))
            radar_info = rbind(radar_info, radar_info_)
        }
    }
}

label_cols = radar_info %>% 
    mutate(col=case_when(Response=='CR' ~ '#5778A2',
                         Response=='PR' ~ '#DCB44B',
                         Response=='SD' ~ '#8A281E')) %>%
    pull(col) #%>% tibble::column_to_rownames('sample')
names(label_cols) = radar_info$sample
# label_cols

stg='Post'
rsp='CR'
radar_input = radar_info %>%
filter(Stage==stg, Response==rsp) %>% 
select(-c('Patient','Stage','Response')) %>% distinct() 

# similarity_df = radar_input %>%
# pivot_longer(cols=-1) %>%
# pivot_wider(names_from = sample, values_from = value) %>%
# corrr::correlate() %>%
# mutate(across(where(is.numeric), .fns = ~replace_na(.,1))) %>%
# arrange(desc(`CRC01-Post`))

pdf(paste0('figres/radar_',stg,'_CR.pdf'), height = 12, width = 8)
radar_input %>%
# mutate(sample = factor(sample, levels=similarity_df$term)) %>%
ggradar(background.circle.colour = "white",
        axis.line.colour = "black",
        gridline.min.colour = "black",
        gridline.mid.colour = "black",
        gridline.max.colour = "black",
        group.colours = label_cols,
        group.point.size=0,
        fill=F,
#         fill.alpha=.25,
        plot.title = stg
       )
