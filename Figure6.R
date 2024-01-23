#---- Figure 6 ----#


# Library and data loading
library(ggplot2)
library(caret)
library(pROC)
require(Seurat)
require(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)

load('raid1/cyq/cd8_blood_obj.rda')
load('raid1/cyq/cd8_tcr_obj.rda')
sample_info = data.frame(fread('obj/clinical_info_09.csv'))

## Clone diversity

diversity_df = cd8_blood_obj@meta.data %>% 
	ungroup() %>% 
    select(Patient, Response, Stage, clone.id) %>% 
    distinct(clone.id, .keep_all = TRUE) %>% 
    ungroup() %>%
    group_by(Patient, Response, Stage) %>% 
    summarise(clone_n = n())

diversity_df %>%
	mutate(Response2 = case_when(Response == 'SD' ~ 'Non-responder', TRUE ~ 'Responder')) %>%
    filter(Stage=='Pre') %>%
ggplot(aes(x=factor(Response2, levels=c('Responder','Non-responder')), y=clone_n)) + 
    geom_boxplot(aes(color=Response2), outlier.color = NA, lwd = 1, width=.5, alpha=.3) +
    geom_point(aes(color=Response2), size=1.5, stroke=1, position = position_jitter(width=.2)) +
    theme_pubr() + 
    scale_color_manual(values = c('#8A281E','#5778A2','#DCB44B'), name='Response') +
    scale_shape_manual(values = c(16,1,17,4), name='MSI state') +
    xlab('') + 
    ylab('Clone Diversity') + 
    ggtitle('Circulating CD8 T cells') +
    theme(strip.text.x = element_text(size = 18),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .5),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(l=10,t=10))


## Scatter plot

input = cd8_tcr_obj@meta.data %>% 
    filter(Tissue=='Blood') %>%
    filter(Clone_type %in% c('New','Expanded')) %>%
    select(Patient, clone.id) %>% distinct(clone.id, .keep_all = TRUE) %>%
    group_by(Patient) %>% summarise(n_exp = n(), log2_n_exp = log2(n_exp)) %>%
    left_join(sample_info, by = 'Patient')

lm = lm(SIZE_CHANGE ~ n_exp, data = input)

# print(paste0('R^2 = ',format(summary(lm)$r.squared, digits = 2)))
# print(paste0('P = ',format(lmp(lm), digits = 2)))

p = input %>%
    ggplot(aes(y = SIZE_CHANGE, x = n_exp)) +
    geom_smooth(method='lm', col='black', lwd=.7) +
    geom_point(aes(col = Response, shape = MSI2), size = 2, stroke=2) +
    scale_color_manual(values = c('#5778A2','#DCB44B','#8A281E')) +
    scale_shape_manual(values = c(16,17,1,4)) +
    labs(title = 'Expanded/new CD8 T clones in Blood') +
    ylab('Tumor regression ratio') + 
    xlab('Clone counts') +
    geom_text(x = 500, y = 1.4, label = paste0('R = ',format(sqrt(summary(lm)$r.squared), digits = 2))) +
    geom_text(x = 500, y = 1.2, label = paste0('P = ',format(lmp(lm), digits = 2))) +
    theme_pubr() + 
    theme(strip.text.x = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = 1), 
          axis.text.x = element_text(size=15, angle = 0, hjust = 1),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .8),
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t = 10)) +
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) 


## UMAP

cd8_blood_obj@meta.data %>% 
ggplot() + 
    geom_point(aes(x=UMAP1, y=UMAP2, col=annotation), size=.1, alpha=.5) + 
    scale_color_manual(values = get_palette('lancet', 9), name='') + 
    theme_nothing() +
    theme(plot.title = element_text(size=22),
          strip.text = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = .1), 
          legend.position = 'bottom',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=50, t=10, r=60)) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1), ncol=2)) 

cd8_blood_obj@meta.data %>% 
ggplot() + 
    geom_point(aes(x=UMAP1, y=UMAP2, col=Ttr), size=.5, alpha=.5) + 
    scale_color_manual(values = c('lightgrey','#8A281E'), name='') + 
    theme_nothing() +
    theme(plot.title = element_text(size=22),
          strip.text = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = .1), 
          legend.position = 'right',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=10, r=10)) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 


## Density plot

cd8_blood_obj@meta.data %>% 
    group_by(clone.id,Ttr) %>% 
    summarise(clone_size = n()) %>%
ggplot(aes(y = Ttr, x = log10(clone_size), col=Ttr, fill=Ttr)) +
    geom_density_ridges2(alpha = 0.5, 
                        jittered_points = F, 
                        scale = 1.5, 
                        position = position_points_jitter(height = 0)) +
    scale_y_discrete(expand = c(.05,.05)) +
    scale_color_manual(values = get_palette('nejm',2)[c(2,1)]) +
    scale_fill_manual(values =get_palette('nejm',2)[c(2,1)]) +
    theme_pubr() + 
    ggtitle(paste0('Circulating CD8 T cell clone')) + 
    ylab('') + 
    xlab('Log10 clone size') +
    theme(strip.text = element_text(size = 15, face='bold'),
          strip.background =element_rect(fill=NA, color = 'black', size = 1), 
          axis.text.x = element_text(size=13, angle = 0, hjust = .5, vjust=.5),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .6),
          legend.position = 'none',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(b=0, l=10, t=10, r=10))


## Differential expression analysis

obj = cd8_blood_obj
Idents(obj) = obj$Ttr
de_blood_cd8_Ttr_nonTtr =FindMarkers(obj, `ident.1` = 'Ttr', `ident.2` = 'Non-Ttr',logfc.threshold = 0, test.use = 't')

### Plot

de_blood_cd8_Ttr_nonTtr %>% 
ggplot() + 
	theme_pubr() +
	xlab('Log2(fold change)') + 
	ylab('-Log10 (adj. P value)') +
	geom_point(aes(x=avg_log2FC, y = -log10(p_val_adj), 
		color = abs(avg_log2FC) > 0 & p_val_adj <=0.05), size=1) + 
	scale_color_manual(values = c('grey', '#BD4739'), name='') + 
	geom_text_repel(aes(x=avg_log2FC, y = -log10(p_val_adj), label=rowname), 
		size=7, max.overlaps=20, fontface = "italic") +
    theme(strip.text.x = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = 1), 
          axis.text.x = element_text(size=20, angle = 0, hjust = .5),
          axis.text.y = element_text(size=20),
          axis.title = element_text(size = 20),
          axis.line = element_line(size = .8),
          legend.position = 'none',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=10, l=10,r=10)) +
	guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))

## Geneset enrichment

de_gene = de_blood_cd8_Ttr_nonTtr %>% arrange(desc(avg_log2FC))

entrez_gene = AnnotationDbi::select(hs, 
	keys = rownames(de_gene),
	columns = c("ENTREZID", "SYMBOL"),
	keytype = "SYMBOL") %>% 
    filter(!is.na(ENTREZID)) %>% 
	distinct(SYMBOL, .keep_all = T) %>% 
	filter(!is.na(ENTREZID))

original_gene_list = de_gene %>% 
	filter(rownames(.) %in% entrez_gene$SYMBOL) %>% 
	pull(avg_log2FC)
names(original_gene_list) = entrez_gene$ENTREZID
gene_list = na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

gse = gseGO(geneList=gene_list, 
             ont ='ALL',
             keyType = "ENTREZID", 
             minGSSize = 10, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = hs,
             pAdjustMethod = "BH"
            )

### Plot

gse %>% 
	as.data.frame() %>% 
	group_by(sign(NES)) %>% 
	slice_max(abs(NES), n=50) %>%
	arrange(desc(NES)) %>%
	ggplot(aes(x=NES, y=fct_reorder(Description, NES), fill=factor(sign(NES), levels=c(-1,1)))) + #-log10(qvalues)
	geom_bar(stat='identity') +
	geom_vline(xintercept = 0) + 
	scale_fill_manual(values = c('#86C0C4','#377DA6'), name='') +
    theme_pubr() + 
	labs(title = 'Pathway enrichment in blood CD8 Tem Tprolif Ttr - Non-Ttr') + 
    theme(strip.text.x = element_text(size = 18),
    	strip.background =element_rect(fill='white', color = 'black', size = 1), 
        axis.text.x = element_text(size=15, angle = 0, hjust = .5),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 15),
        axis.line.y = element_blank(),
        legend.position = 'none',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t=0, l=0)) +
    xlab('NES') +
    ylab('Pathway') 


## Lasso regression

blood_obj = cd8_blood_obj
lasso_input = data.frame(
    t(data.frame(blood_obj@assays$RNA@data[mhc_ii_protein_complex,])),
    Patient = blood_obj$Patient,
    Response = blood_obj$Response,
    Stage = blood_obj$Stage
) %>% group_by(Patient, Response, Stage) %>% 
    summarise_all(list(mean)) %>%
    filter(Response %in% c('CR','PR')) %>%
    filter(Stage=='Pre') %>% 
    ungroup() %>%
    select(-Stage,-Patient) %>% 
    mutate(Response = case_when(Response=='CR' ~ 1, Response=='PR' ~ 0))

lasso_input$Response = as.factor(lasso_input$Response)
head(lasso_input)


ctrlspecs = trainControl(method = 'LOOCV',  savePredictions = 'all')
lambda_vector = 10^seq(0, -3, length=500)

set.seed(123)
model_pre = train(Response ~ .,
                  data = lasso_input_pre,
                  preProcess = c('center','scale'),
                  method = 'glmnet',
                  tuneGrid = expand.grid(alpha=1, lambda=lambda_vector),
                  trControl = ctrlspecs,
                  na.action = na.omit,
                  intercept = F,
                  family = "binomial"
             ) %>% suppressWarnings()

# head(model_pre$bestTune$lambda)
round(coef(model_pre$finalModel, model_pre$bestTune$lambda), 3)


### predictive score

predict_score = c('CD74','HLA-DQA2','HLA-DQB1','HLA-DRA')
blood_obj = AddModuleScore(blood_obj, features = list(predict_score), name = 'predict_score')
HNSCC_obj = AddModuleScore(HNSCC_obj, features = list(predict_score), name = 'predict_score')
cd8_predict_df = blood_obj@meta.data %>%
    filter(Stage %in% c('Pre'), Response != 'SD') %>%
    group_by(Patient, Response) %>% 
    summarise(score = mean(predict_score1)) %>%
    mutate(rsp=case_when(Response == 'CR' ~ 1, Response == 'PR' ~0))

#### Baseline prediction

cd8_predict_df %>%
ggplot(aes(x=Response, y=score, group=Response)) +
    geom_point(aes(color=Response), size=1, stroke=1, position = position_jitter(width=.2)) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = .8, width=.5, fill=NA) +
    theme_pubr() + 
    scale_color_manual(values = c('#5778A2','#DCB44B'), name='Response') + 
    xlab('') + ylab('Four-gene score')+
    theme(strip.text = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = 1), 
          axis.text.x = element_text(size=15, angle = 0, hjust = .5, vjust=.5),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .8),
          legend.position = 'none',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=10, l=10, r=50))
# ROC curve

res.roc = roc(rsp ~ score, data=cd8_predict_df)  
plot.roc(res.roc, print.auc = TRUE)

#### Score dynamics

cd8_predict_dnmc_df = blood_obj@meta.data %>%
    filter(Response != 'SD') %>%
    group_by(Patient, Stage, Response) %>% summarise(score = mean(predict_score1))

cd8_predict_dnmc_df %>%
filter(Stage %in% c('Pre','Post')) %>%
ggplot(aes(x=factor(Stage, levels=c('Pre','Post')), y=score)) +
    geom_line(aes(group=Patient), lwd=.3) +
    geom_point(aes(color=Response), size=1, stroke=1) +
    geom_boxplot(aes(color=Response), outlier.color = NA, lwd = .8, width=.5, fill=NA) +
    theme_pubr() + 
    scale_color_manual(values = c('#5778A2','#DCB44B'), name='Response') +
    facet_wrap(~Response) +
    xlab('') + ylab('MHC II 4-gene score')+
    theme(strip.text = element_text(size = 18),
          strip.background =element_rect(fill='white', color = 'black', size = 1), 
          axis.text.x = element_text(size=15, angle = 0, hjust = .5, vjust=.5),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 18),
          axis.line = element_line(size = .8),
          legend.position = 'none',
          legend.key.size = unit(1.4,'line'),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15, face = 'bold'),
          plot.margin = margin(t=10, l=10, r=50))
