library(data.table)
library(MOFA2)

# Preparing data by centring and scaling

intra.lip = read.csv('data/IntracellularLipidomicsRaw.csv',row.names = 1)
rownames(intra.lip)=intra.lip$LipidIon
intra.lip = intra.lip[,4:ncol(intra.lip)]
# Per-total normalisation
intra.lip = as.data.frame(t(t(intra.lip)/colSums(intra.lip)))
# Take average of replicates for each cell line
intra.lip=data.frame('Control.A'=rowMeans(intra.lip[,1:3]),
                     'Control.B'=rowMeans(intra.lip[,4:6]),
                     'Control.C'=rowMeans(intra.lip[,7:9]),
                     'Control.A.treatment'=rowMeans(intra.lip[,10:12]),
                     'Control.B.treatment'=rowMeans(intra.lip[,13:15]),
                     'Control.C.treatment'=rowMeans(intra.lip[,16:18]),
                     'PMS.A'=rowMeans(intra.lip[,19:21]),
                     'PMS.B'=rowMeans(intra.lip[,22:24]),
                     'PMS.C'=rowMeans(intra.lip[,25:27]),
                     'PMS.D'=rowMeans(intra.lip[,28:30]),
                     'PMS.A.treatment'=rowMeans(intra.lip[,31:33]),
                     'PMS.B.treatment'=rowMeans(intra.lip[,34:36]),
                     'PMS.C.treatment'=rowMeans(intra.lip[,37:39]),
                     'PMS.D.treatment'=rowMeans(intra.lip[,40:42]))
intra.lip[intra.lip==0]=NA
# Centre and scale
intra.lip=t(scale(t(intra.lip)))


secret = read.csv('data/SecretomeRaw.csv',row.names = 1)
rownames(secret)=paste0(secret$Gene,'_',secret$Protein.ID)
secret = secret[,1:(ncol(secret)-1)]
secret = secret[,9:ncol(secret)]
secret = secret[,!(colnames(secret)%like%'NC')]
# Per-total normalisation
secret = as.data.frame(t(t(secret)/colSums(secret)))
# Take average of replicates for each cell line
secret=data.frame('Control.A'=rowMeans(secret[,1:3],na.rm=T),
                  'Control.B'=rowMeans(secret[,4:6],na.rm=T),
                  'Control.C'=rowMeans(secret[,7:9],na.rm=T),
                  'Control.A.treatment'=rowMeans(secret[,10:12],na.rm=T),
                  'Control.B.treatment'=rowMeans(secret[,13:15],na.rm=T),
                  'Control.C.treatment'=rowMeans(secret[,16:18],na.rm=T),
                  'PMS.A'=rowMeans(secret[,19:21],na.rm=T),
                  'PMS.B'=rowMeans(secret[,22:24],na.rm=T),
                  'PMS.C'=rowMeans(secret[,25:27],na.rm=T),
                  'PMS.D'=rowMeans(secret[,28:30],na.rm=T),
                  'PMS.A.treatment'=rowMeans(secret[,31:33],na.rm=T),
                  'PMS.B.treatment'=rowMeans(secret[,34:36],na.rm=T),
                  'PMS.C.treatment'=rowMeans(secret[,37:39],na.rm=T),
                  'PMS.D.treatment'=rowMeans(secret[,40:42],na.rm=T))
# Filter out proteins only expressed in 20% of samples
secret = secret[rowSums(secret!=0)>0.2*ncol(secret),]
secret[secret==0]=NA
# Centre and scale
secret=t(scale(t(secret)))

# Create MOFA object
MOFAobject <- create_mofa(list('IntracellularLipidomics'=as.matrix(intra.lip),
                               'Secretomics'=as.matrix(secret)))

# Set options to default except for 3 factors
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors = 3
head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)


# Prepare MOFA
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
outfile = file.path(getwd(),"treatmentmodel_minimalfiltering.hdf5")
reticulate::use_python('/usr/bin/python3')
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

print(plot_data_overview(MOFAobject.trained))

model = MOFAobject.trained
model <- impute(model)
Nsamples = sum(model@dimensions$N)

# Sample metadata
sample_metadata <- data.frame(
  sample = samples_names(model)[[1]],
  condition = c(rep('Control',6),rep('PMS',8)),
  treatment = c(rep(c('None','Treated'),each=3),rep(c('None','Treated'),each=4)),
  condition_treatment = paste0(c(rep('Control',6),rep('PMS',8)),'_',
                               c(rep(c('None','Treated'),each=3),
                                 rep(c('None','Treated'),each=4))),
  patient = c('Control_A','Control_B','Control_C','Control_A_treatment',
              'Control_B_treatment','Control_C_treatment',
              'PMS_A','PMS_B','PMS_C','PMS_D' ,'PMS_A_treatment',
              'PMS_B_treatment','PMS_C_treatment','PMS_D_treatment')
)

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)

# Total variance explained per view and group
head(model@cache$variance_explained$r2_total[[1]]) # group 1
# Variance explained for every factor in per view and group
head(model@cache$variance_explained$r2_per_factor[[1]]) # group 1
plot_variance_explained(model, x="factor", y="view")
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
# Plot out sample distribution in each factor
plot_factor(model, 
            factor = 1:3,
            color_by = "treatment",
            shape_by = "condition"
)
plot_factor(model, 
            factor = 1:3,
            color_by = "condition",
            shape_by = "treatment"
)
plot_factor(model, 
            factor = 1:3,
            color_by = "condition_treatment",
            shape_by = "condition_treatment"
)
plot_factor(model, 
            factor = 1:3,
            color_by = "patient",
            shape_by = "condition_treatment"
)

library(ggplot2)
p <- plot_factor(model, 
                 factors = c(1,2,3),
                 color_by = "condition",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c("Control"="black", "PMS"="red")) +
  scale_fill_manual(values=c("Control"="black", "PMS"="red"))

print(p)

p <- plot_factor(model, 
                 factors = c(1,2,3),
                 color_by = "treatment",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

p <- p + 
  scale_color_manual(values=c("None"="black", "Treated"="red")) +
  scale_fill_manual(values=c("None"="black", "Treated"="red"))

print(p)

p <- plot_factor(model, 
                 factors = c(1,2,3),
                 color_by = "condition_treatment",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

print(p)

plot_factors(model, 
             factors = 1:3,
             color_by = "condition"
)
plot_factors(model, 
             factors = 1:3,
             color_by = "treatment"
)
plot_factors(model, 
             factors = 1:3,
             color_by = "condition_treatment"
)

# Take top features in each factor and plot them out
for (modality in c('IntracellularLipidomics','Secretomics')){
  print(plot_weights(model,
                     view = modality,
                     factor = 1,
                     nfeatures = 10,     # Number of features to highlight
                     scale = T,          # Scale weights from -1 to 1
                     abs = F,text_size = 3           # Take the absolute value?
  ))
  print(plot_weights(model,
                     view = modality,
                     factor = 2,
                     nfeatures = 10,     # Number of features to highlight
                     scale = T,          # Scale weights from -1 to 1
                     abs = F,text_size = 3           # Take the absolute value?
  ))
}

for (modality in c('IntracellularLipidomics','Secretomics')){
  print(plot_top_weights(model,
                         view = modality,
                         factor = 1,
                         nfeatures = 30,
  ))
  print(plot_top_weights(model,
                         view = modality,
                         factor = 2,
                         nfeatures = 30
  ))
  if (modality!='Secretomics'){
    print(plot_data_heatmap(model,
                            view = modality,         # view of interest
                            factor = 1,             # factor of interest
                            features = 30,   # number of features to plot (they are selected by weight)
                            
                            # extra arguments that are passed to the `pheatmap` function
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE
    ))
    print(plot_data_heatmap(model,
                            view = modality,         # view of interest
                            factor = 2,             # factor of interest
                            features = 30,          # number of features to plot (they are selected by weight)
                            # extra arguments that are passed to the `pheatmap` function
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE
    ))
  }
  
  print(plot_data_scatter(model,
                          view = modality,         # view of interest
                          factor = 1,             # factor of interest
                          features = 5,           # number of features to plot (they are selected by weight)
                          add_lm = TRUE,          # add linear regression
                          color_by = "condition_treatment"
  ))
  print(plot_data_scatter(model,
                          view = modality,         # view of interest
                          factor = 2,             # factor of interest
                          features = 5,           # number of features to plot (they are selected by weight)
                          add_lm = TRUE,          # add linear regression
                          color_by = "condition_treatment"
  ))
  
  if (modality=='Secretomics'){
    print(plot_data_heatmap(model,
                            view = modality,         # view of interest
                            factor = 1,             # factor of interest
                            features = 20,          # number of features to plot (they are selected by weight)
                            # extra arguments that are passed to the `pheatmap` function
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE
    ))
    
    print(plot_data_heatmap(model,
                            view = modality,         # view of interest
                            factor = 2,             # factor of interest
                            features = 20,          # number of features to plot (they are selected by weight)
                            # extra arguments that are passed to the `pheatmap` function
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE
    ))
    print(plot_data_heatmap(model,
                            view = modality,         # view of interest
                            factor = 1,             # factor of interest
                            features = 30,          # number of features to plot (they are selected by weight)
                            imputed=T,
                            # extra arguments that are passed to the `pheatmap` function
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE
    ))
    
    print(plot_data_heatmap(model,
                            view = modality,         # view of interest
                            factor = 2,             # factor of interest
                            features = 30,          # number of features to plot (they are selected by weight)
                            imputed=T,
                            # extra arguments that are passed to the `pheatmap` function
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE
    ))
    print(plot_data_scatter(model,
                            view = modality,         # view of interest
                            factor = 1,             # factor of interest
                            features = 5,           # number of features to plot (they are selected by weight)
                            add_lm = TRUE,          # add linear regression
                            color_by = "condition_treatment",imputed = T
    ))
    print(plot_data_scatter(model,
                            view = modality,         # view of interest
                            factor = 2,             # factor of interest
                            features = 5,           # number of features to plot (they are selected by weight)
                            add_lm = TRUE,          # add linear regression
                            color_by = "condition_treatment",imputed = T
    ))
    
  }
  
  
}


a=correlate_factors_with_covariates(model,covariates=c('condition','treatment'),plot='r',method='color',abs=T,transpose=T,return_data=T)
a = data.frame(a)
a$metadata = rownames(a)
a = tidyr::pivot_longer(a,cols=1:3)
a$metadata = ifelse(a$metadata=='condition','disease',a$metadata)
a$metadata = factor(a$metadata,levels=c('treatment','disease'))
ggplot(a,aes(x=name,y=metadata,fill=value))+
  geom_tile()+ 
  scale_fill_gradient2(low = "white", high = "darkred",limits=c(0,1))+
  theme_classic()+xlab('')+ylab('')+labs(fill='corr')


