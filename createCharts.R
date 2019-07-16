# Analyse X versus Y nuclei in WT, Yqdel and shSLY.
# This script should be run via RScript, and will generate
# all figure templates and results tables

# .libPaths(c("Y:/bms41/pkg/3.5.3",.libPaths()))
# setwd("Y:/bms41/Mice/Yqdel warping/Pre-FISH/2019-02-21_16-40-56/")
source("functions.R") # Generic functions

load.generic.functions()

# Packages that are needed for the full analysis
install.missing(packages = c( "devtools", "Rcpp", "glue", "tibble", "colorspace", "tidyverse", "dplyr", "tidyr", "lmtest",
                             "stringr", "cowplot", "class", "mvtnorm", "scales",
                             "Matrix", "survey", "ggplot2", "svglite",
                             "cluster", "factoextra", "mclust", "robustbase", "dendextend", "PerformanceAnalytics", 
                             "relaimpo", "MASS", "coin", "DescTools",  "rcompanion", "Rtsne", "raster", 
                             "betareg", "pscl", "broom", "umap"))

# Quick saving & split violin plots
load.gg.extensions()
source("analysis_functions.R") # Other analyses
seed = 42

# Read in the raw data, filter, distinguish X and Y, and make plotting factors
# The data is loaded into the global environment
read.and.process.data = function(){
  
  cat("Reading sample data\n")
  
  sample.orders = c("WT P1", "WT P2", "WT I1", "WT I2", "WT I3",
                    "Yqdel P1", "Yqdel P2", "Yqdel I1", "Yqdel I2",
                    "Yqdel I3", "Yqdel I4", "Yqdel I5", 
                    "WT 1", "WT 2",
                    "Yqdel 1", "Yqdel 2", "Yqdel 3", "Yqdel 4",
                    "shSLY 1", "shSLY 2", "shSLY 3", "shSLY 4")
  
  read.data = function(file) read.csv(file, header = T, stringsAsFactors = F, sep="\t")
  
  read.C57.data = function(){
    C57.data = read.data("Raw_data/C57_XY_stats.txt")
    
    genotype = c("Cajou 125" = "WT", "Cajou 181" = "WT","Cajou 184" = "WT","Cajou 126" = "WT",
                 "Cajou 120" = "shSLY", "Cajou 180" = "shSLY","Cajou 183" = "shSLY","Cajou 185" = "shSLY",
                 "Mac 86" = "Yqdel", "Mac 88" = "Yqdel", "Mac 89" = "Yqdel", "Mac 92" = "Yqdel")
        
    C57.filt = C57.data %>% mutate(    Sample = str_extract(Dataset, "(Cajou|Mac)\\s\\d+"),
                                       Population = str_extract(Dataset, "(Left|Right)"),
                                       Bkg    = str_extract(Dataset, "(C57|MF1)"),
                                       Fluor  = str_extract(Dataset, "(TxR|FITC)"),
                                       Chr    = str_extract(Dataset, "(X|Y)$")) %>%
      dplyr::mutate(Type   = genotype[Sample]) %>%
      dplyr::mutate(Strain = paste(Bkg, Type)) %>%
      dplyr::filter(Difference_from_median<35) %>% # Remove last potatoes
      dplyr::ungroup() %>% 
      dplyr::mutate(Sample = case_when(Sample=="Cajou 184" ~ "WT 1",
                                       Sample=="Cajou 181" ~ "WT 2",
                                       Sample=="Mac 86" ~ "Yqdel 1",
                                       Sample=="Mac 88" ~ "Yqdel 2",
                                       Sample=="Mac 89" ~ "Yqdel 3",
                                       Sample=="Mac 92" ~ "Yqdel 4",
                                       Sample=="Cajou 120" ~ "shSLY 1",
                                       Sample=="Cajou 180" ~ "shSLY 2",
                                       Sample=="Cajou 183" ~ "shSLY 3",
                                       Sample=="Cajou 185" ~ "shSLY 4"))
    C57.filt$Sample  = factor(C57.filt$Sample, levels = sample.orders)
    C57.filt
  }
  
  C57.filt <<- read.C57.data()
  
  read.MF1.data = function(){
    
    MF1.data = read.data("Raw_data/MF1_XY_stats.txt") # excludes cells with abnormal swelling
    
    # These come from the following wrk files:
    # J = \\cgs-fs5\filespace\BBSRC & Leverhulme work
    # J:\Images\Sperm morphology and NO\2015-04-14 Kent\DAPI\2019-03-02_12-36-35\2015-04-14 Analysis.wrk
    # J:\Images\Sperm morphology and NO\Yq-del\2015-11-16 Analysis.wrk
    # J:\Images\Sperm morphology and NO\2017-05-19 Emma Capture-recapture\Emma's analysis.wrk
    # Y:\bms41\Mice\Yqdel warping\Pre-FISH\2019-02-21_16-40-56\2016-07-28 Analysis.wrk
    # Only the mapped X and Y sperm were exported
        
    MF1.data = MF1.data %>%
      dplyr::mutate(Type = str_extract(Dataset, "(WT|Yq)"),
                    Sample = str_extract(Dataset, "(Yq(del)?|WT)\\s?([R|r]ep|Ind|[A|B]) \\d?"),
                    Population = str_extract(Dataset, "(Left|Right)"),
                    Fluor  = str_extract(Dataset, "(Tx?R|FITC)"),
                    Chr    = str_extract(Dataset, "(X|Y)$"),
                    Bkg    = "MF1" ) %>%
      dplyr::mutate(Type  = gsub("Yq", "Yqdel", Type),
                    Fluor = gsub("TR", "TxR", Fluor)) %>%
      dplyr::mutate(Strain = paste(Bkg, Type)) %>%
      dplyr::filter(Difference_from_median<35)  %>%
      # dplyr::filter(Difference_from_median<35,
                    # Sample!="Yq Rep 4")  %>%  # This sample has too few nuclei to be useful
      dplyr::ungroup() %>%
      dplyr::mutate(Sample = case_when(Sample=="WT rep 1" ~ "WT P1",
                              Sample=="WT Rep 2" ~ "WT P2",
                              Sample=="WT Ind 1" ~ "WT I1",
                              Sample=="WTA "      ~ "WT I2",
                              Sample=="WTB "      ~ "WT I3",
                              Sample=="Yq Rep 5" ~ "Yqdel P1",
                              Sample=="Yqdel Rep 5" ~ "Yqdel P1",
                              Sample=="Yq Rep 6" ~ "Yqdel P2",
                              Sample=="Yq Ind 4" ~ "Yqdel I1",
                              Sample=="Yq Ind 5" ~ "Yqdel I2",
                              Sample=="Yq Ind 6" ~ "Yqdel I3",
                              Sample=="YqA "     ~ "Yqdel I4",
                              Sample=="YqB "     ~ "Yqdel I5"))
    MF1.data$Sample  = factor(MF1.data$Sample, levels = sample.orders)
    MF1.data
  }
  
  MF1.filt <<- read.MF1.data() 
  
  full.data  = dplyr::bind_rows(C57.filt, MF1.filt)
  chr.order = c("X", "Y")
  full.data$Chr = factor(full.data$Chr, levels = chr.order)
  genotype.order = c("MF1 WT", "MF1 Yqdel", "C57 WT", "C57 Yqdel", "C57 shSLY")
  
  full.data$Sample = factor(full.data$Sample, levels = sample.orders)
  full.data$Strain = factor(full.data$Strain, levels = genotype.order)
  full.data
}

# Create plots for the grouped data and the individual sample data
make.sample.plots = function(){
  cat("Creating sample plots\n")
  # Create the group and individual plots for each background
  make.plot = function(plot.function, file.name, file.width){
    vars      = c("Area_square_microns", "Difference_from_median")
    yaxts     = list( expression( paste("Area (", mu, m^2, ")") ), "Variability" )
    plot.list = mapply(plot.function, vars, yaxts, SIMPLIFY = F)
    do.call(plot_grid, c(plot.list, align='v', nrow = length(vars), ncol = 1))
    save.function = switch(file.width, "single" = ggsave.single,
                           "double" = ggsave.double,
                           "row"    = ggsave.row,
                           ggsave.double)
    save.function(file.name)
    dev.off()
  }
  
  # Area and variability plots for all data
  make.group.plot = function(yvar, yaxt){
    hline.data = full.data %>% dplyr::select_("Strain", "Sample", "Chr", yvar) %>% 
      dplyr::group_by_("Strain") %>% 
      dplyr::summarise_at(yvar, median) %>% mutate_(Median = yvar)
    hline.data$Chr = "X"
    hline.Y = hline.data
    hline.Y$Chr = "Y"
    hline.data = rbind(hline.data, hline.Y)
    
    n.data = full.data %>% dplyr::select_("Strain", "Sample", "Chr", yvar) %>% 
      dplyr::group_by(Strain, Chr) %>% dplyr::summarise(n=n())
    
    n.data.y = max(full.data[[yvar]])*1.1 # go 10% above the max for n labels
    
    ggplot(full.data, aes(x=Chr, y=full.data[[yvar]], fill = Chr)) +
      geom_hline(aes(yintercept=Median), hline.data) +
      geom_violin(trim=TRUE, lwd=0.9) +
      geom_boxplot(width=0.1, lwd=0.9, outlier.shape = NA) +
      geom_text(data=n.data, aes(x=Chr, y=n.data.y, label=paste0("n=",n)), size=2)+
      facet_grid(.~Strain)+
      ylab( yaxt )+
      theme_classic() +
      theme( legend.position="none", axis.title.x = element_blank())
  }
  make.plot(make.group.plot, "All data strain area", "single")
  
  # Area and variability plots for individual samples
  make.individual.plot = function(sample.data, yvar, yaxt){
    
    # Create median column and difference column
    # Set the y value based on the variable being drawn
    # Variability should be absolute, others should be difference
    # to sample median.
    diff.data = sample.data %>% 
      dplyr::select_("Sample", "Chr", yvar) %>% 
      dplyr::group_by(Sample) %>% 
      dplyr::mutate_(Value = yvar) %>% 
      dplyr::mutate(Median = median(Value)) %>%
      dplyr::mutate(Diff = Value - Median) %>%
      dplyr::mutate(hlineValue = case_when(yvar=="Difference_from_median" ~ Median,
                                           TRUE ~ 0),
                    yValue     = case_when(yvar=="Difference_from_median" ~ Value,
                                           TRUE ~ Diff))
    
    n.data = diff.data %>% 
      dplyr::group_by(Sample, Chr, Median, hlineValue) %>% 
      dplyr::summarise(n=n())
    n.data.y = max((diff.data$yValue))
    
    ggplot(diff.data, aes(x=Chr, y=yValue)) +
      geom_hline(data=n.data, aes(yintercept=hlineValue)) +
      geom_violin(trim=TRUE, lwd=0.9, aes(fill=Chr)) +
      geom_boxplot(width=0.1, lwd=0.9, outlier.shape = NA, aes(fill=Chr)) +
      geom_text(data=n.data, aes(x=Chr, y=n.data.y, label=paste0("n=",n)), size=2)+
      facet_grid(.~Sample)+
      ylab(yaxt )+
      theme_classic() +
      theme( legend.position="none", axis.title.x = element_blank())
  }
  
  # Individual MF1 sample plots
  make.MF1.individual.plot = function(yvar, yaxt) make.individual.plot(MF1.filt, yvar, yaxt)
  make.plot(make.MF1.individual.plot, "MF1.individual", "row")
  
  # Individual C57 sample plots
  make.C57.individual.plot = function(yvar, yaxt) make.individual.plot(C57.filt, yvar, yaxt)
  make.plot(make.C57.individual.plot, "C57.individual", "row")
}

# tSNE tests revealed 100x1000 are reasonable parameters.
# Run these on the full data, and use clustering to divide the
# tSNE matrix.
run.tsne.clustering = function(sample.data, sample.name, n.clusters, seed){
  cat("Running tSNE and clustering on", sample.name, "with", n.clusters, "clusters from", nrow(sample.data), "samples using seed", seed,"\n")
  
  # Read saved files where possible since tSNE and agnes take time
  tSNE.file = paste0("Rdata/",sample.name, " tSNE values ",seed,".Rdata")
  
  create.tSNE.values = function(){
    if(!file.exists(tSNE.file)){
      cat("No tSNE values file\n")
      input = sample.data %>% ungroup() %>% dplyr::select(one_of(paste0("Angle_profile_", seq(0,99,1))))
      sample.data$Strain = as.factor(sample.data$Strain)
      rownames(input) = paste0(sample.data$Strain, "_", sample.data$CellID)
      cat("Running tSNE on", sample.name, "\n")
      # run Rtsne with optimised parameters = perp=100, iter=1000
      set.seed(seed)
      rtsne_out = Rtsne(as.matrix(input), perplexity=100, max_iter=1000)
      values    = as.data.frame(rtsne_out$Y)
      saveRDS(values, file=tSNE.file)
      return(values)
    } else {
      cat("Reading tSNE values file\n")
      return(readRDS(tSNE.file))
    }
  }
  values = create.tSNE.values()
  
  
  umap.file = paste0("Rdata/",sample.name, " UMAP values ",seed,".Rdata")
  create.UMAP.values = function(){
    if(!file.exists(umap.file)){
      cat("No UMAP values file\n")
      input = sample.data %>% ungroup() %>% dplyr::select(one_of(paste0("Angle_profile_", seq(0,99,1))))
      sample.data$Strain = as.factor(sample.data$Strain)
      rownames(input) = paste0(sample.data$Strain, "_", sample.data$CellID)
      cat("Running UMAP on", sample.name, "\n")
      # run UMAP
      set.seed(seed)
      
      umap_out = umap(input, random_state=seed)
      values   = as.data.frame(umap_out$layout)
      saveRDS(values, file=umap.file)
      return(values)
    } else {
      cat("Reading UMAP values file\n")
      return(readRDS(umap.file))
    }
  }
  umap.values = create.UMAP.values()
  
  dend.file = paste0("Rdata/",sample.name, " tSNE agnes dendrogram ",seed,".Rdata")
  
  create.dendrogram = function(){
    if(!file.exists(dend.file)){
      cat("No tSNE dendrogram file; clustering\n")
      # Cluster the tSNE results with agnes
      hc = cluster::agnes(values, method = "ward")
      cat("Clustering using ward, AC=", format(hc$ac, digits=2), "\n")
      dend = as.dendrogram(hc)
      saveRDS(dend, file=dend.file)
      return(dend)
    } else {
      cat("Reading tSNE dendrogram file\n")
      return(readRDS(dend.file))
    }
  }
  dend = create.dendrogram()
  
  # Map the clusters back to the input data
  cat("Cutting dendrogram\n")
  clusters = dendextend::cutree(dend, n.clusters)
  sample.data$agnes  = clusters
  
  # Set the grades by cluster - note, this is by manual assigment on seed 42
  # A = normal
  # B, B/C, C = increasing grades of acrosomal curvature and hook shortening
  # D = mild tail socket malformation with dorsal angle accentuation
  # E = severe abnormality, with sperm shortening and basal widening
  # F = severe abnormality with hook effacement and basal compression.
  cat("Adding factor labels\n")
  mut.data = sample.data %>% dplyr::mutate(agnes = case_when(agnes == 1 & Bkg=="C57" ~ "A3",
                                                             agnes == 2 & Bkg=="C57" ~ "S",
                                                             agnes == 3 & Bkg=="C57" ~ "B2",
                                                             agnes == 4 & Bkg=="C57" ~ "B1",
                                                             agnes == 5 & Bkg=="C57" ~ "A2",
                                                           agnes == 6 & Bkg=="C57" ~ "A1",
                                                             agnes == 7 & Bkg=="C57" ~ "N",
                                                             agnes == 1 & Bkg=="MF1" & n.clusters==4 ~ "N2",
                                                             agnes == 2 & Bkg=="MF1" & n.clusters==4 ~ "N1",
                                                             agnes == 3 & Bkg=="MF1" & n.clusters==4 ~ "A1",
                                                             agnes == 4 & Bkg=="MF1" & n.clusters==4 ~ "A2",
                                                             TRUE ~ as.character(agnes))) # default
  
  # Set the order of clusters
  cluster.orders = switch(paste0(sample.name, n.clusters), "C57 data7" = c("N", "A1", "A2", "A3", "B1", "B2", "S"),
                          "MF1 data4" = c("N1", "N2", "A1", "A2"))
  
  mut.data$agnes = factor(mut.data$agnes, levels=cluster.orders)
  
  values$agnes   = factor(mut.data$agnes, levels=cluster.orders)
  values$agnes  = as.factor(mut.data$agnes)
  values$Strain = mut.data$Strain
  values$Sample = mut.data$Sample
  values$Chr    = mut.data$Chr
  
  umap.values$agnes   = factor(mut.data$agnes, levels=cluster.orders)
  umap.values$agnes  = as.factor(mut.data$agnes)
  umap.values$Strain = mut.data$Strain
  umap.values$Sample = mut.data$Sample
  umap.values$Chr    = mut.data$Chr
  
  plot.agnes.dendrogram = function(){
    cat("Plotting agnes dendrogram\n")
    # colours = hue_pal()(n.clusters)
    colours = as.numeric(mut.data$agnes)
    colours = colours[order.dendrogram(dend)]
    labels_colors(dend) = colours

    # Draw the tree
    png(filename = paste0("Report/",sample.name, " ", n.clusters, "-cluster tSNE dendrogram seed ",seed, ".png"), res=300, width=170*1.25, height=170*1.25, units="mm")
    plot(dend, main = paste0(sample.name, " tSNE clustering"))
    dev.off()
  }
  plot.agnes.dendrogram()
  
  cat("Plotting tSNE\n")
  # Plot the tSNE results with no colouring
  black.plot = ggplot(values, aes(x=V1, y=V2))+
    geom_point(size=0.5)+
    theme_void()+theme(legend.position = "none")

  # Plot the cluster results along with the actual division
  plot.tsne = function(colour.column){
    ggplot(values, aes(x=V1, y=V2, col=values[[colour.column]]))+
      geom_point(size=0.5)+
      ggtitle(colour.column)+
      theme_void()+theme(legend.position = "none")
  }
  
  # Plot the UMAP data with the tSNE cluster results
  plot.umap = function(colour.column){
    ggplot(umap.values, aes(x=V1, y=V2, col=values[[colour.column]]))+
      geom_point(size=0.5)+
      ggtitle(colour.column)+
      theme_void()+theme(legend.position = "none")
  }

  three.cols = hue_pal()(3)
  orig.plot  = plot.tsne("Strain")
  ggsave.single(paste0(sample.name, " ", n.clusters, "-cluster tSNE strain seed ",seed))
  
  clust.plot = plot.tsne("agnes")+scale_colour_brewer(palette="Dark2")
  ggsave.single(paste0(sample.name, " ", n.clusters, "-cluster tSNE agnes seed ",seed))
  sm.plot    = plot.tsne("Sample")
  ggsave.single(paste0(sample.name, " ", n.clusters, "-cluster tSNE sample seed ",seed))
  xy.plot    = plot.tsne("Chr")
  ggsave.single(paste0(sample.name, " ", n.clusters, "-cluster tSNE chr seed ",seed))
  
  umap.orig.plot  = plot.umap("Strain")
  umap.clust.plot = plot.umap("agnes")+scale_colour_brewer(palette="Dark2")
  ggsave.single(paste0(sample.name, " ", n.clusters, "-cluster UMAP agnes only seed ",seed))
  umap.sm.plot    = plot.umap("Sample")
  umap.xy.plot    = plot.umap("Chr")

  make.strain.tSNE.plot = function(strain){
    ggplot(values %>% dplyr::filter(Strain==strain), aes(x=V1, y=V2, col=Chr))+
      geom_point(size=0.5)+
      ggtitle(strain)+
      theme_void()+theme(legend.position = "top")
  }
  strain.plots = mapply(make.strain.tSNE.plot, unique(values$Strain), SIMPLIFY = FALSE)
  plot.row.count = ceiling(length(unique(values$Strain))/2)
  do.call(plot_grid, c(strain.plots, align='hv', nrow = plot.row.count, ncol = 2))
  ggsave.custom(paste0(sample.name, " ", n.clusters, "-cluster tSNE strains seed ",seed), width=170, height=plot.row.count*85)

  do.call(plot_grid, c( list(orig.plot, sm.plot, xy.plot, clust.plot), align='hv', nrow = 2, ncol = 2))
  ggsave.custom(paste0(sample.name, " ", n.clusters, "-cluster tSNE seed ",seed), width=170, height=170)
  
  do.call(plot_grid, c( list(umap.orig.plot, umap.sm.plot, umap.xy.plot, umap.clust.plot), align='hv', nrow = 2, ncol = 2))
  ggsave.custom(paste0(sample.name, " ", n.clusters, "-cluster UMAP seed ",seed), width=170, height=170)

  # Export the clusters as map files
  # Syntax is CellID<\t>cluster number
  map.data = sample.data %>% dplyr::select(CellID, agnes)
  write.table(map.data, file = paste0("Report/",sample.name, " ", n.clusters, "-cluster tSNE map seed ", seed, ".tsv"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  
  make.strain.umap.plot = function(strain){
    ggplot(umap.values %>% dplyr::filter(Strain==strain), aes(x=V1, y=V2, col=Chr))+
      geom_point(size=0.5)+
      ggtitle(strain)+
      theme_void()+theme(legend.position = "top")
  }
  umap.strain.plots = mapply(make.strain.umap.plot, unique(umap.values$Strain), SIMPLIFY = FALSE)
  plot.row.count = ceiling(length(unique(umap.values$Strain))/2)
  do.call(plot_grid, c(umap.strain.plots, align='hv', nrow = plot.row.count, ncol = 2))
  ggsave.custom(paste0(sample.name, " ", n.clusters, "-cluster UMAP strains seed ",seed), width=170, height=plot.row.count*85)
  
  mut.data
}

calculate.cluster.stats = function(mut.data, sample.name, n.clusters, seed){

  cat("Generating cluster summaries\n")
  
  # Calcuate proportions for agnes cluster and export
  summ = mut.data %>% 
    dplyr::select(Strain, agnes, Chr) %>%
    dplyr::group_by(Strain, Chr) %>% 
    dplyr::mutate(Strain_total = n()) %>%
    dplyr::group_by(Strain, agnes, Chr, Strain_total) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::mutate( PctStrain = n/Strain_total)
  
  write.headed.table(summ, file = paste0(sample.name, " ", n.clusters, "-cluster tSNE strain proportions seed ", seed,  ".tsv"))

  summ.with.error = summ %>% dplyr::select(-Strain_total, -PctStrain) %>%
    tidyr::spread(Chr, n) %>%
    dplyr::mutate(FrX = X/(X+Y),
                  Sep = SEP(FrX, X+Y)) %>%
    dplyr::mutate(PctX = FrX*100, SepX = Sep*100)
  
  write.headed.table(summ.with.error, file = paste0(sample.name, " ", n.clusters, "-cluster tSNE proportions and SEP seed ", seed,".tsv"))
  
  cat("Calculating X:Y ratios\n")
  
  # Plot the pct of total strain within each cluster
  summ.by.cluster = mut.data %>% 
    dplyr::group_by(Type) %>%
    dplyr::mutate(Total = n()) %>%
    dplyr::group_by(Type, agnes, Chr, Total) %>% 
    dplyr::summarise(Number=n()) %>%
    tidyr::spread(Chr, Number) %>%
    dplyr::mutate(PctX = X/(Total), PctY = Y/(Total)) %>%
    dplyr::mutate(SEM = SEP(PctX, Total)) %>%
    dplyr::select(-X, -Y, -Total) %>%
    tidyr::gather("Chr", "Pct",-Type, -agnes, -SEM) %>%
    dplyr::mutate(Chr = gsub("Pct", "", Chr), 
                  Comb = paste(Type, Chr), # variable for grouping in plot
                  Pct = Pct*100, 
                  SEM = SEM*100) 
  summ.by.cluster$Type = factor(summ.by.cluster$Type, levels = c("WT", "Yqdel", "shSLY"))
  write.headed.table(summ.by.cluster, file=paste0(sample.name, " ", n.clusters, "-cluster XvY by cluster seed ",seed, ".tsv"))
  
  ggplot(summ.by.cluster, aes(x=agnes, y=Pct, group=Comb, col=Type)) +
    geom_point(position= position_dodge(width=0.4), size=3, aes(shape=Chr))+
    geom_line(aes(linetype=Type), position=position_dodge(width=0.4))+
    geom_errorbar(aes(ymin=Pct-SEM, ymax=Pct+SEM), width=0.2, position=position_dodge(width=0.4), col="black")+
    ylab("Percent of strain")+
    ylim(0, 40)+
    scale_colour_manual(values = c("#057451", "#2852af", "#d74810"))+
    scale_linetype_manual(values=c("solid", "longdash", "dashed"))+
    scale_shape_manual(values = c(16, 17))+
    theme(legend.position = "top", axis.title.x = element_blank())
  ggsave.custom(paste0(sample.name, " ", n.clusters, "-cluster XvY by cluster seed",seed), width=85, height=70)

  # What is the natural X:Y ratio for each genotype based on our sampling? i.e. our 50:50 baseline?
  natural.ratios = mut.data %>% dplyr::group_by(Strain, Chr) %>%
    dplyr::summarise(n=n()) %>%
    tidyr::spread(Chr, n) %>%
    dplyr::mutate(X_ratio = X/(X+Y))
  write.headed.table(natural.ratios, file = paste0(sample.name, " ", n.clusters, "-cluster natural XY ratios seed ", seed,".tsv"))
  
  # Summarise the data by strain, and calculate the proportion of sperm in each cluster
  summ.by.strain = summ %>% dplyr::select( -n) %>% tidyr::spread(agnes, PctStrain) %>% ungroup()
  write.headed.table(summ.by.strain,paste0(sample.name, " ", n.clusters, "-cluster strain proportions seed ", seed, ".tsv"))
  
  # Create profile charts for each of the clusters
  make.profile.charts = function(){
    # Generate profile charts for each cluster
    profile.data = mut.data %>% ungroup() %>% dplyr::select(Dataset, agnes, Chr, matches("Angle_profile_\\d+$")) %>%
      tidyr::gather(Profile_position, value, -Dataset, -agnes, -Chr ) %>%
      dplyr::group_by(Chr, agnes, Profile_position) %>%
      dplyr::mutate(Position = as.numeric(gsub("Angle_profile_", "", Profile_position)),
                    Median = median(value),
                    Q25 = quantile(value, 0.25),
                    Q75 = quantile(value, 0.75)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-value, -Profile_position) %>%
      dplyr::distinct()
    
    ggplot(profile.data, aes(x=Position, y=Median, fill=agnes, group=Chr))+
      geom_hline(yintercept=180)+
      geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 1)+
      geom_line(aes(col=Chr))+
      xlab("Position")+ylab("Angle")+labs(fill = "Cluster")+
      ylim(50,250)+
      scale_fill_brewer(palette="Dark2")+
      theme_classic() +
      facet_wrap(~agnes)+
      theme(legend.position = "top")
    
    plot.row.count = ifelse(n.clusters<4, 1, 2)
    ggsave.custom(paste0(sample.name, " ", n.clusters, "-cluster profiles seed ", seed), 170, 85*plot.row.count)
    profile.data
  }
  cat("Making profile charts\n")
  profile.data = make.profile.charts()
  
  cat("Making violin charts\n")
  # Make area and variability charts for the clustered data
  # Area and variability plots for all data
  make.clusters.plot = function(yvar, yaxt, input.data){
    hline.data = input.data %>% dplyr::select_("Strain", "agnes", yvar) %>% 
      dplyr::group_by_("agnes") %>% 
      dplyr::summarise_at(yvar, median) %>% dplyr::mutate_(Median = yvar)
    
    n.data = input.data %>% dplyr::select_("agnes", yvar) %>% 
      dplyr::group_by(agnes) %>% dplyr::summarise(n=n())
    
    n.data.y = max(input.data[[yvar]])*1.1 # go 10% above the max for n labels
    
    ggplot(input.data, aes(x=agnes, y=input.data[[yvar]], fill = agnes)) +
      geom_hline(aes(yintercept=median(input.data[[yvar]]))) +
      geom_violin(trim=TRUE, lwd=0.9) +
      geom_boxplot(width=0.1, lwd=0.9, outlier.shape = NA) +
      geom_text(data=n.data, aes(x=agnes, y=n.data.y, label=paste0("n=",n)), size=2)+
      ylab( yaxt )+
      scale_fill_brewer(palette="Dark2")+
      theme_classic() +
      theme( legend.position="none", axis.title.x = element_blank())
  }
  vars      = c("Area_square_microns", "Difference_from_median")
  yaxts     = list( expression( paste("Area (", mu, m^2, ")") ), "Variability" ) 
  plot.list = mapply(make.clusters.plot, vars, yaxts, list(mut.data, mut.data), SIMPLIFY = F)
  do.call(plot_grid, c(plot.list, align='v', nrow = length(vars), ncol = 1))
  ggsave.row(paste0(sample.name, " ", n.clusters, "-cluster violin plots seed ", seed))
}

# Create the figure
create.area.acrosome.figure = function(){
  # Make a mean and sem error bar plot for the given variable
  make.sem.plot = function(variable, ylabel){
    point.size = 1.5
    dodge.width = 0.3
    agg.data = full.data %>% ungroup() %>% 
      dplyr::select("Strain", "Type", "Chr", variable) %>%
      dplyr::mutate_(Variable = variable) %>%
      dplyr::group_by(Strain, Type, Chr) %>%
      dplyr::summarise(N = n(),
                        SD = sd(Variable),
                        SEM = sem(Variable),
                        CIL = MeanCI(Variable, conf.level = 0.99)[2],
                        CIU = MeanCI(Variable, conf.level = 0.99)[3],
                        Mean = mean(Variable)) 
    agg.data$Type = factor(agg.data$Type, levels = c("WT", "Yqdel", "shSLY"))
    
    agg.data$Type = factor(agg.data$Type, levels = c("WT", "Yqdel", "shSLY"))
    ggplot(agg.data, aes(x=Strain, y=Mean, col=Type, shape=Chr))+
      geom_point(size=point.size, position=position_dodge(width=dodge.width))+
      geom_errorbar(aes(ymin=CIL, ymax=CIU), width=0.2, position=position_dodge(width=dodge.width))+
      scale_colour_manual(values = create.strain.colours())+
      scale_alpha_manual(values = c(1, 0.5))+
      scale_x_discrete(labels = c("WT", "Yqdel", "WT", "Yqdel", "shSLY"))+
      guides(alpha=FALSE)+
      ylab(ylabel)+
      theme_classic()+
      theme(axis.title.x = element_blank(), legend.position = "none")
  }
  
  create.acro.violin = function(){
    acrosome.data = full.data %>% ungroup() %>% 
      dplyr::select(Strain, Type, Chr, matches("Angle_profile_\\d+$")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Acrosome = Angle_profile_85) %>%
      dplyr::select(-contains("Angle_profile")) 
    
    acrosome.data$Type = factor(acrosome.data$Type, levels = c("WT", "Yqdel", "shSLY"))
    
    agg.data = acrosome.data %>%
      dplyr::group_by(Strain, Type, Chr) %>%
      dplyr::summarise(Acrosome = mean(Acrosome),
                       SEM = sem(Acrosome))
    agg.data$Type = factor(agg.data$Type, levels = c("WT", "Yqdel", "shSLY"))

    ggplot(acrosome.data, aes(x=Strain, y=Acrosome, fill=Type, alpha=Chr))+
      geom_split_violin(draw_quantiles=c(.50))+
      scale_fill_manual(values = create.strain.colours())+
      scale_alpha_manual(values = c(1, 0.5))+
      scale_x_discrete(labels = c("WT", "Yqdel", "WT", "Yqdel", "shSLY"))+
      guides(alpha=FALSE)+
      ylab("Index 85 angle")+
      theme_classic()+
      theme(axis.title.x = element_blank(), legend.position = "none")
  }
  acro.violin = create.acro.violin()

  acro.sem = make.sem.plot("Angle_profile_85", "Index 85 angle")
  area.sem = make.sem.plot("Area_square_microns", expression(paste("Area (", mu, m^2, ")")))
  
  # Make the area XvY plot using a split violin
  create.area.violin = function(){
    full.data$Type = factor(full.data$Type, levels = c("WT", "Yqdel", "shSLY"))
    ggplot(full.data, aes(x=Strain, y=Area_square_microns, fill=Type, alpha=Chr))+
      geom_split_violin(draw_quantiles=c(.50))+
      scale_fill_manual(values = create.strain.colours())+
      scale_alpha_manual(values = c(1, 0.5))+
      guides(alpha=FALSE)+
      ylab( expression(paste("Area (", mu, m^2, ")")))+
      scale_x_discrete(labels = c("WT", "Yqdel", "WT", "Yqdel", "shSLY"))+
      theme_classic()+
      theme(axis.title.x = element_blank(), legend.position = "none")
    # ggsave.custom("Area XvY", width=85, height=45)
  }
  area.violin = create.area.violin()

  
  # Create figure plots
  # Acrosome XvY and Y/Y
  # Area XvY and Y/X
  create.area.dotplot = function(){
    
    all.summ =  full.data %>% 
      dplyr::group_by(Type, Strain, Sample, Chr) %>% 
      dplyr::select(Strain, Sample, Chr, Type, Area_square_microns) %>% 
      dplyr::summarise(Mean = mean(Area_square_microns), 
                       SEM  = sem(Area_square_microns)) %>%
      dplyr::group_by(Type, Strain, Sample, Chr, Type) %>%
      tidyr::gather("Function", "Value", -Type, -Strain, -Sample, -Chr) %>%
      dplyr::mutate(ChrFunc = paste0(Function, Chr)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-Chr, -Function) %>%
      tidyr::spread(ChrFunc, Value) %>%
      dplyr::mutate(Diff = MeanY/MeanX,
                    SemXNorm = SEMX/MeanX, SemYNorm = SEMY/MeanY) %>%
      dplyr::mutate(SemBoth = sqrt((SemXNorm^2+SemYNorm^2))) %>%
      dplyr::mutate(SemUpper = Diff+SemBoth, SemLower = Diff-SemBoth) %>%
      dplyr::mutate(Type = factor(Type, levels = c("WT", "Yqdel", "shSLY")))
    
    
    wt.mean = mean(all.summ[all.summ$Type=="WT",]$Diff)
    
    ggplot(all.summ, aes(x=Strain, y=Diff, col=Type, group=Strain))+
      geom_hline(yintercept = wt.mean)+
      geom_pointrange(data=all.summ,aes(x=Strain,ymin=SemLower, ymax=SemUpper),
                      size=0.5, 
                      position = position_jitterdodge(seed=1, jitter.width=0.9, dodge.width=0.1))+
      ylab("Y/X mean area ratio")+
      scale_colour_manual(values = create.strain.colours())+
      scale_x_discrete(labels = c("WT", "Yqdel", "WT", "Yqdel", "shSLY"))+
      theme_classic() +
      theme( legend.position="none", axis.title.x = element_blank())
    
  }
  area.dots = create.area.dotplot()

  
  print.to.png = function(){
    # Arrange for a 5 panel layout
    # Move to a new page
    png(filename = "Report/Figure xxxx - acro and area combined.png", width = 170, height = 85, units = "mm", res = 300)
    grid.newpage()
    # Create layout : nrow = 2, ncol = 2
    pushViewport(viewport(layout = grid.layout(2, 3)))
    # A helper function to define a region on the layout
    define_region <- function(row, col){
      viewport(layout.pos.row = row, layout.pos.col = col)
    } 
    # Arrange the plots
    print(acro.violin, vp=define_region(1, 1))
    print(area.violin, vp = define_region(1, 2))
    print(acro.sem, vp = define_region(2, 1))
    print(area.sem, vp=define_region(2, 2))
    print(area.dots, vp=define_region(1:2, 3))
    dev.off()
  }
  print.to.png()
  
  print.to.svg = function(){
    # Arrange for a 5 panel layout
    # Move to a new page
    svg(filename = "Report/Figure xxxx - acro and area combined.svg", width = 6.6929, height = 3.3464)
    grid.newpage()
    # Create layout : nrow = 2, ncol = 2
    pushViewport(viewport(layout = grid.layout(2, 3)))
    # A helper function to define a region on the layout
    define_region <- function(row, col){
      viewport(layout.pos.row = row, layout.pos.col = col)
    } 
    # Arrange the plots
    print(acro.violin, vp=define_region(1, 1))
    print(area.violin, vp = define_region(1, 2))
    print(acro.sem, vp = define_region(2, 1))
    print(area.sem, vp=define_region(2, 2))
    print(area.dots, vp=define_region(1:2, 3))
    dev.off()
  }
  print.to.svg()

}

# The text has various values comparing X and Y.
# Invent them here.
define.percentages.for.text = function(){
  cat("Defining percentages\n")
  all.summ =  full.data %>% 
    dplyr::group_by(Type, Strain, Sample, Chr) %>% 
    dplyr::select(Strain, Sample, Chr, Type, Area_square_microns) %>% 
    dplyr::summarise(Mean = mean(Area_square_microns), 
                     SEM  = sem(Area_square_microns)) %>%
    dplyr::group_by(Type, Strain, Sample, Chr, Type) %>%
    tidyr::gather("Function", "Value", -Type, -Strain, -Sample, -Chr) %>%
    dplyr::mutate(ChrFunc = paste0(Function, Chr)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Chr, -Function) %>%
    tidyr::spread(ChrFunc, Value) %>%
    dplyr::mutate(Diff = MeanY/MeanX,
                  SemXNorm = SEMX/MeanX, SemYNorm = SEMY/MeanY) %>%
    dplyr::mutate(SemBoth = sqrt((SemXNorm^2+SemYNorm^2))) %>%
    dplyr::mutate(SemUpper = Diff+SemBoth, SemLower = Diff-SemBoth) %>%
    dplyr::mutate(Type = factor(Type, levels = c("WT", "Yqdel", "shSLY"))) %>% na.omit
  
  # Not needed unless referees ask for it?
  print.strain.values = function(strain){
    values = all.summ[all.summ$Strain==strain,]
    diff = (1-mean(values$Diff))
    vol = 1-((1-diff))^(3/2)
    
    sems = sqrt(sum(values$SemBoth^2))
    vsems = 1-((1-sems))^(3/2)
    
    cat(strain, ": Y is on average", format(diff*100, digits=2), " +/-", format(sems*100, digits=2), "% smaller than X by area\n",
        file = paste0("Report/XY average difference.txt"), append=T)
    cat(strain, ": Y is on average", format(vol*100, digits=2), " +/-", format(vsems*100, digits=2), "% smaller than X by volume\n",
        file = paste0("Report/XY average difference.txt"), append=T)
  }
  invisible(lapply(levels(full.data$Strain), print.strain.values))
  
  print.type.values = function(type){
    values = all.summ[all.summ$Type==type,]
    diff = (1-mean(values$Diff))
    vol = 1-((1-diff))^(3/2)
    
    sems = sqrt(sum(values$SemBoth^2))
    vsems = 1-((1-sems))^(3/2)
    
    cat(type, ": Y is on average", format(diff*100, digits=2), " +/-", format(sems*100, digits=2), "% smaller than X by area\n",
        file = paste0("Report/XY average difference.txt"), append=T)
    cat(type, ": Y is on average", format(vol*100, digits=2), " +/-", format(vsems*100, digits=2), "% smaller than X by volume\n",
        file = paste0("Report/XY average difference.txt"), append=T)
  }
  invisible(lapply(unique(full.data$Type), print.type.values))
  
  strain.values = full.data %>% 
    dplyr::group_by(Bkg, Type, Strain, Chr) %>% 
    dplyr::select(Bkg, Strain,  Chr, Type, Area_square_microns) %>% 
    dplyr::summarise(N=n(),
                     Mean = mean(Area_square_microns), 
                     SEM  = sem(Area_square_microns)) %>%
    dplyr::group_by(Bkg, Type, Strain, Chr, Type) %>%
    tidyr::gather("Function", "Value", -Bkg, -Type, -Strain, -Chr) %>%
    dplyr::mutate(ChrFunc = paste0(Function, Chr)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Chr, -Function) %>%
    tidyr::spread(ChrFunc, Value) %>%
    dplyr::mutate(Ratio = MeanY/MeanX,
                  SemXNorm = SEMX/MeanX, SemYNorm = SEMY/MeanY) %>%
    dplyr::mutate(SemBoth = sqrt((SemXNorm^2+SemYNorm^2))) %>%
    dplyr::mutate(SemUpper = Ratio+SemBoth, SemLower = Ratio-SemBoth) %>%
    dplyr::mutate(Type = factor(Type, levels = c("WT", "Yqdel", "shSLY")),
                  Bkg = factor(Bkg, levels = c("MF1", "C57"))) %>% na.omit %>%
    dplyr::mutate(YFctDiff = (1-Ratio),
                  YVolSEM = 1-((1-SemBoth))^(3/2)) %>%
    dplyr::mutate(YVolDiff = 1-((1-YFctDiff))^(3/2)) %>%
    dplyr::select(-SEMX, -SEMY, -SemXNorm, -SemYNorm, -SemUpper, -SemLower, -Strain) %>%
    dplyr::mutate(Sample = "Aggregate") %>%
    dplyr::arrange(Bkg, Type)
  
  write.headed.table(values, file="Table xxxx - YvX.tsv")
  
  sample.values = full.data %>%
    dplyr::group_by(Bkg, Type, Strain, Chr, Sample) %>% 
    dplyr::select(Bkg, Strain,  Chr, Type, Sample, Area_square_microns) %>% 
    dplyr::summarise(N=n(),
                     Mean = mean(Area_square_microns), 
                     SEM  = sem(Area_square_microns)) %>%
    dplyr::group_by(Bkg, Type, Strain, Chr, Type, Sample) %>%
    tidyr::gather("Function", "Value", -Bkg, -Type, -Strain, -Chr, -Sample) %>%
    dplyr::mutate(ChrFunc = paste0(Function, Chr)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Chr, -Function) %>%
    tidyr::spread(ChrFunc, Value) %>%
    dplyr::mutate(Ratio = MeanY/MeanX,
                  SemXNorm = SEMX/MeanX, SemYNorm = SEMY/MeanY) %>%
    dplyr::mutate(SemBoth = sqrt((SemXNorm^2+SemYNorm^2))) %>%
    dplyr::mutate(SemUpper = Ratio+SemBoth, SemLower = Ratio-SemBoth) %>%
    dplyr::mutate(Type = factor(Type, levels = c("WT", "Yqdel", "shSLY")),
                  Bkg = factor(Bkg, levels = c("MF1", "C57"))) %>% na.omit %>%
    dplyr::mutate(YFctDiff = (1-Ratio),
                  YVolSEM = 1-((1-SemBoth))^(3/2)) %>%
    dplyr::mutate(YVolDiff = 1-((1-YFctDiff))^(3/2)) %>%
    dplyr::select(-SEMX, -SEMY, -SemXNorm, -SemYNorm, -SemUpper, -SemLower, -Strain) %>%
    dplyr::arrange(Bkg, Type)
  
  values = rbind(strain.values, sample.values)
  
  write.headed.table(values, file="Table xxxx - YvX per-sample.tsv")
}

# Run a regression or other regression across the Yqdel phenotype abnormalities only
# i.e for MF1: (N1+N2) vs A1 vs A2
# for C57: N vs A1 - A3
# We are looking for differences in the proportion of X and Y
# Use a chi square test, then Fisher's exact test
run.yqdel.explicit.x.v.y.tests = function(mut.data, sample.name){
  cat("Testing XvY ratios\n")
  sub.data = mut.data %>% dplyr::mutate(sub.agnes = case_when(agnes=="N1" ~ "N",
                                                   agnes=="N2" ~ "N",
                                                   TRUE~as.character(agnes))) %>%
    dplyr::filter(sub.agnes %in% c("N", "A1", "A2", "A3")) %>%
    dplyr::mutate(sub.agnes = factor(sub.agnes, levels=c("N", "A1", "A2", "A3")))
  
  orig.summ = mut.data %>% 
    dplyr::select(Strain, agnes, Chr) %>%
    dplyr::group_by(Strain, Chr) %>% 
    dplyr::mutate(Strain_total = n()) %>%
    dplyr::group_by(Strain, agnes, Chr, Strain_total) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::mutate( PctStrain = n/Strain_total) %>%
    dplyr::mutate(SEP = SEP(PctStrain, n))
  
  summ = sub.data %>% 
    dplyr::select(Strain, sub.agnes, Chr) %>%
    dplyr::group_by(Strain, Chr) %>% 
    dplyr::mutate(Strain_total = n()) %>%
    dplyr::group_by(Strain, sub.agnes, Chr, Strain_total) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::mutate( PctStrain = n/Strain_total) %>%
    dplyr::mutate(SEP = SEP(PctStrain, n))
  
  summ.by.chr = sub.data %>% 
    dplyr::select(Strain, sub.agnes, Chr) %>%
    dplyr::group_by(sub.agnes, Chr) %>% 
    dplyr::summarise(n=n()) %>%
    tidyr::spread(Chr, n) %>%
    dplyr::mutate( PctX = X/(X+Y)) %>%
    dplyr::mutate(SEP = SEP(PctX, (X+Y)))
  
  run.prop.test = function(cluster){
    values = summ.by.chr %>% dplyr::ungroup() %>% 
      dplyr::filter(sub.agnes==cluster)
    result = prop.test(values$X, values$X+values$Y, p=0.5)
    adj.p = p.adjust(result$p.value, method = "BH", n=length(levels(droplevels(sub.data$sub.agnes))))
    cat(cluster, ": Z-test p=", adj.p, "\n", file=paste0("Report/", sample.name, " XvY Z-tests.txt"), append=T)
  }
  cat("Testing X:Y ratios Z-test against 50:50\n")
  tryCatch(invisible(lapply(levels(droplevels(sub.data$sub.agnes)), run.prop.test)),
           error=function(e) cat("Error calculating Z-test values for strain\n") )
  
  
  ggplot(summ.by.chr, aes(x=sub.agnes, y=PctX))+
    geom_hline(yintercept = 0.5)+
    ylab("Fraction X")+
    geom_point()+
    geom_errorbar(aes(ymin=PctX-SEP, ymax=PctX+SEP))+
    theme_classic()+
    theme(axis.title.x = element_blank())
  ggsave.single(paste0(sample.name, " X ratio by cluster"))
  
  chi.test.strain = function(strain){
    values = summ %>% dplyr::ungroup() %>%
      dplyr::filter(Strain==strain)
    
    x.values = values %>% dplyr::filter(Chr=="X")
    y.values = values %>% dplyr::filter(Chr=="Y")
    chi.result = chisq.test(x.values$n, y.values$n, simulate.p.value = TRUE)
    adj.p = p.adjust(chi.result$p.value, method = "BH", n=length(unique(sub.data$Strain)))
    cat(strain, ": Chi test XVY p=", adj.p, "\n", file=paste0("Report/", sample.name, " XvY chi-tests.txt"), append=T)
  }

  cat("Testing X:Y ratios are different between strains by chi-square\n")
  tryCatch(invisible(lapply(unique(sub.data$Strain), chi.test.strain)),
           error=function(e) cat("Error calculating chi-square values for strain\n") )

  fishers.exact.test = function(strain){
    # make a contingency table of counts - XYvsAgnes
    cluster.strain = summ %>% dplyr::ungroup() %>%
      dplyr::filter(Strain==strain) %>%
      dplyr::select(-Strain,-PctStrain, -Strain_total, -SEP) %>%
      tidyr::spread(Chr, n) %>%
      dplyr::select(-sub.agnes) %>%
      dplyr::mutate(X = as.numeric(X), Y=as.numeric(Y)) %>%
      as.matrix %>% t

    fisher.results = fisher.test(cluster.strain,alternative="two.sided", simulate.p.value=TRUE)
    adj.p = p.adjust(fisher.results$p.value, method = "BH", n=length(levels(droplevels(sub.data$sub.agnes))))
    cat(strain, ": Fisher's Exact Test for count data (X/Y vs cluster): p =", adj.p, "\n", file=paste0("Report/", sample.name, " XvY sub-cluster Fisher tests.txt"), append=T)
   
    cluster.original = orig.summ %>% dplyr::ungroup() %>%
      dplyr::filter(Strain==strain) %>%
      dplyr::select(-Strain,-PctStrain, -Strain_total, -SEP) %>%
      tidyr::spread(Chr, n) %>%
      dplyr::select(-agnes) %>%
      dplyr::mutate(X = as.numeric(X), Y=as.numeric(Y)) %>%
      as.matrix %>% t
    
    fisher.results = fisher.test(cluster.original,alternative="two.sided", simulate.p.value=TRUE)
    adj.p = p.adjust(fisher.results$p.value, method = "BH", n=length(levels(droplevels(mut.data$agnes))))
    cat(strain, ": Fisher's Exact Test for count data (X/Y vs cluster): p =", adj.p, "\n", file=paste0("Report/", sample.name, " XvY full Fisher tests.txt"), append=T)
  }

  cat("Testing X:Y ratios are different between strains by Fisher's exact test\n")
  tryCatch(invisible(lapply(unique(sub.data$Strain), fishers.exact.test)),
           error=function(e) cat("Error calculating Fisher values for cluster\n"))
  
  
  # run a logistic regession on X/Y counts versus cluster and strain
  cat("Building a logistic regression model for cluster and strain\n")
  # regression building a model across all strains
  summ.by.chr.strain = sub.data %>% 
    dplyr::select(Strain, sub.agnes, Chr) %>%
    dplyr::group_by(Strain, sub.agnes, Chr) %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Chr = as.factor(Chr),
                  Strain = as.factor(Strain),
                  sub.agnes = as.factor(sub.agnes)) %>%
    dplyr::mutate(int.agnes = unclass(sub.agnes))

  # per strain modelling
  run.logit.on.strain = function(strain) {
    values = summ.by.chr.strain %>% dplyr::ungroup() %>%
      dplyr::filter(Strain==strain)

    cat(strain, "logistic regression\n")
    
    # Treat cluster as interval data
    int.model = glm(Chr~int.agnes, family=binomial(link = "logit"),data=values, weights = n)
    
    # Treat cluster as categorical data
    full.model = glm(Chr~sub.agnes, family=binomial(link = "logit"),data=values, weights = n)
    null.model = glm(Chr~1, family=binomial(link = "logit"),data=values, weights = n)
    write.headed.table(tidy(int.model), paste0(strain, " logistic regresssion interval model.tsv"))
    write.headed.table(tidy(full.model), paste0(strain, " logistic regresssion categorical model.tsv"))

    chi = anova(full.model, test="Chisq")
    write.headed.table(tidy(chi), paste0(strain, " logistic regresssion model anova.tsv"))
    
    int.chi = anova(full.model, int.model, test="LRT")
    write.headed.table(tidy(int.chi), paste0(strain, " logistic regresssion model categorical to interval model.tsv"))
    
    lr.data = lrtest(int.model, full.model)
    write.headed.table(tidy(lr.data), paste0(strain, " logistic regresssion model categorical to interval LR test.tsv"))
    
    cat(strain, "int AIC:", AIC(int.model), "\n", file = "Report/Logistic regression AIC.txt", append=T)
    cat(strain, "full AIC:", AIC(full.model), "\n", file = "Report/Logistic regression AIC.txt", append=T)
    cat(strain, "null AIC:", AIC(null.model), "\n", file = "Report/Logistic regression AIC.txt", append=T)
  }
  
  tryCatch(invisible(lapply(unique(sub.data$Strain), run.logit.on.strain)),
           error=function(e) cat("Error calculating logistic regression values for strain\n"))
}

run.pipeline = function(seed){
  
  # Clear old results
  unlink("Report/*")
  
  # Store the session info
  writeLines(capture.output(devtools::session_info()), "Report/Session info.txt")

  full.data <<- read.and.process.data()
  # make.sample.plots()
  MF1.mut <<- run.tsne.clustering(MF1.filt,"MF1 data", 4, seed)
  C57.mut <<- run.tsne.clustering(C57.filt,"C57 data", 7, seed) 
  
  calculate.cluster.stats(MF1.mut,"MF1 data", 4, seed)
  calculate.cluster.stats(C57.mut,"C57 data", 7, seed)

  # run.statistical.tests()
  
  # calculateRelativeImportanceOfFactors(MF1.mut)
  # calculateRelativeImportanceOfFactors(MF1.mut)
  
  # calculateMagnitudeDifferencesOfParameters(full.data, "All data", "mean")

  # createProfileDecilePlots(C57.filt, "C57 data")
  # createProfileDecilePlots(MF1.filt, "MF1 data")
  
  # create.area.acrosome.figure()
  
  # define.percentages.for.text()
  
  # run.yqdel.explicit.x.v.y.tests(MF1.mut, "MF1" )
  # run.yqdel.explicit.x.v.y.tests(C57.mut, "C57")
}

run.pipeline(seed)

