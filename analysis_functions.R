# Misc analysis functions

check.hook.length = function(){
  filt.dat = full.data %>% 
    dplyr::filter(Width_of_body_microns+Length_of_hook_microns >= Bounding_width_microns &
                    Length_of_hook_microns>0 & Width_of_body_microns>0)
  
  filt.dat = filt.dat[filt.dat$Length_of_hook_microns<quantile(filt.dat$Length_of_hook_microns, 0.99) &
                        filt.dat$Length_of_hook_microns>quantile(filt.dat$Length_of_hook_microns, 0.01),]
  
  filt.dat = filt.dat %>%  dplyr::group_by(Type, Strain, Sample, Chr) %>%
    dplyr::select(Strain,  Chr, Type, Sample, Bounding_width_microns) %>%
    dplyr::summarise(Mean = mean(Bounding_width_microns),
                     SEM  = sem(Bounding_width_microns)) %>%
    dplyr::group_by(Type, Strain, Chr, Type) %>%
    tidyr::gather("Function", "Value", -Type, -Strain, -Chr, -Sample) %>%
    dplyr::mutate(ChrFunc = paste0(Function, Chr)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Chr, -Function) %>%
    tidyr::spread(ChrFunc, Value) %>%
    dplyr::mutate(Diff = MeanY/MeanX,
                  SemXNorm = SEMX/MeanX, SemYNorm = SEMY/MeanY) %>%
    dplyr::mutate(SemBoth = sqrt((SemXNorm^2+SemYNorm^2))) %>%
    dplyr::mutate(SemUpper = Diff+SemBoth, SemLower = Diff-SemBoth) %>%
    dplyr::mutate(Type = factor(Type, levels = c("WT", "Yqdel", "shSLY"))) %>% na.omit
  
  
  ggplot(filt.dat, aes(x=Strain, y=Diff, group=Sample, col=Type))+
    geom_hline(yintercept = 1)+
    geom_point(position = position_dodge(width=0.3))+
    geom_errorbar(aes(ymin=SemLower, ymax=SemUpper), width=0.2, position = position_dodge(width=0.3))+
    scale_colour_manual(values = create.strain.colours())+
    ylab("Y/X hook length ratio")+
    theme_classic()+
    theme(axis.title.x = element_blank(), legend.position = "none")
  ggsave.single("Hook length")
  
  ggplot(filt.dat, aes(x=Strain, y=Length_of_hook_microns, fill=Type, alpha=Chr))+
    geom_split_violin(draw_quantiles=c(.50))+
    scale_fill_manual(values = create.strain.colours())+
    scale_alpha_manual(values = c(1, 0.5))+
    scale_x_discrete(labels = c("WT", "Yqdel", "WT", "Yqdel", "shSLY"))+
    guides(alpha=FALSE)+
    ylab("Hook length")+
    theme_classic()+
    theme(axis.title.x = element_blank(), legend.position = "none")
}

# What is the coefficient of variation in wild type samples?
check.cov = function(){
  filt.data = full.data %>% dplyr::filter(Type=="WT") %>%
    dplyr::group_by(Bkg, Sample) %>%
    dplyr::select(Bkg, Sample, Chr, Area_square_microns) %>%
    dplyr::summarise(Mean = mean(Area_square_microns),
                     CV = raster::cv(Area_square_microns))
}

createProfileDecilePlots = function(sample.data, sample.name){
  
  quantile.data = sample.data %>% ungroup() %>% 
    dplyr::select(Strain, Chr, matches("Angle_profile_\\d+$")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(AcrosomeMean = mean(Angle_profile_75 : Angle_profile_85)) %>%
    dplyr::select(-contains("Angle_profile")) %>%
    dplyr::group_by(Strain) %>%
    dplyr::mutate(AcrosomeQuantile = findInterval(AcrosomeMean,sort(AcrosomeMean))/length(AcrosomeMean)) %>%
    dplyr::mutate(AcrosomeBin = cut(AcrosomeQuantile,10, include.lowest=T, labels=F))
  
  
  profile.data = quantile.data %>%
    dplyr::group_by(Strain, AcrosomeBin) %>%
    dplyr::mutate(AcrosomeMeanMean = mean(AcrosomeMean)) %>%
    dplyr::group_by(Strain, Chr, AcrosomeBin, AcrosomeMeanMean) %>%
    dplyr::summarise(N = n()) %>% #
    tidyr::spread(Chr, N) %>%
    dplyr::mutate(FractX = X/(X+Y))
  
  p0 =  ggplot(quantile.data, aes(x=AcrosomeMean, y=AcrosomeQuantile, col=Strain))+
    geom_line()+xlab("Mean of indexes 75-85")+ylab("Quantile of angle")+theme(legend.position = "top")
  
  p3 =  ggplot(quantile.data, aes(x=AcrosomeQuantile, y=AcrosomeBin, col=Strain))+
    geom_line()+ylab("Bin of quantile")+xlab("Quantile of angle")+theme(legend.position = "top")+scale_y_continuous(breaks=c(1:10), labels=c(1:10))
  
  p1 = ggplot(profile.data, aes(x=AcrosomeBin, y=FractX, col=Strain))+
    geom_hline(yintercept = 0.5)+xlab("Bin of quantile")+ylab("Fraction of X")+
    geom_point()+geom_line()+theme(legend.position = "top")+scale_x_continuous(breaks=c(1:10), labels=c(1:10))
  
  p2 = ggplot(profile.data, aes(x=AcrosomeBin, y=AcrosomeMeanMean, col=Strain))+
    geom_point()+geom_line()+xlab("Bin of quantile")+ylab("Mean of index means in bin")+
    theme(legend.position = "top")+scale_x_continuous(breaks=c(1:10), labels=c(1:10))
  
  do.call(plot_grid, c(list(p0, p3, p1, p2), align="v", nrow=2, ncol=2))
  
  ggsave.double(paste0(sample.name, " Quantiles and bins"))
}

# Identify significant factors influencing nuclear area
calculateRelativeImportanceOfFactors = function(sample.data){
  
  # We know that nuclear shape and size is driven by a number of factors, each with a small contribution.
  # Which of these factors has the greatest impact?
  q1 = sample.data %>% dplyr::ungroup() %>% 
    dplyr::select(Strain, Sample, Chr, starts_with("Length_seg_"), Area_square_microns) %>%
    dplyr::mutate(Strain = as.numeric(Strain), Sample = as.numeric(Sample), Chr = as.numeric(Chr))
  
  inter.formula = formula("Area_square_microns ~ Chr + Strain +Sample")
  
  fit1 = lm(inter.formula, data=q1)
  # step = stepAIC(fit1, direction="both")
  # step$anova # display results 
  # eval(step$call)
  # 
  # boot = boot.relimp(fit1, b = 1000, type = c("lmg","last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)
  # booteval.relimp(boot) # print result
  # png("segments.png", width=400, height=200, res=300, units="mm")
  # plot(booteval.relimp(boot,sort=TRUE)) # plot result 
  # dev.off()
}

# Create the profile plots for Figure 2
make.profile.plots = function(){
  cat("Creating profile plots\n")
  profile.data = full.data %>% ungroup() %>% dplyr::select(Bkg, Type, Strain, matches("Angle_profile_\\d+$")) %>%
    tidyr::gather(Profile_position, value, -Strain, -Bkg, -Type) %>%
    dplyr::group_by(Strain, Profile_position) %>%
    dplyr::mutate(Position = as.numeric(gsub("Angle_profile_", "", Profile_position)),
                  Median = median(value),
                  Q25 = quantile(value, 0.25),
                  Q75 = quantile(value, 0.75)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-value, -Profile_position) %>%
    dplyr::distinct()
  
  fill.colours = create.strain.colours()
  
  profile.data$Type = with(profile.data, factor(Type, levels = c("WT", "Yqdel", "shSLY")))
  profile.data$Bkg = with(profile.data, factor(Bkg, levels = c("MF1", "C57")))
  
  ggplot(profile.data, aes(x=Position, y=Median, fill=Type))+
    geom_rect(xmin=7, xmax=18, ymin=-Inf, ymax=Inf, fill="grey")+
    geom_rect(xmin=80, xmax=90, ymin=-Inf, ymax=Inf, fill="grey")+
    geom_hline(yintercept=180)+
    geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.5)+
    geom_line(aes(col=Type))+
    xlab("Position")+ylab("Angle")+labs(fill = "Type")+
    ylim(50,250)+
    scale_fill_manual(values=fill.colours)+
    scale_colour_manual(values=fill.colours)+
    facet_wrap(~Bkg)+
    theme(legend.position = "top")
  ggsave.custom("Figure 3B - profiles ", width=145,  height=70)
  
}


# What power do we have to find differences in the datasets?
calculate.powers = function(){
  library(pwr)
  # What is the effect size for two populations with differing proportions of 1 and 0.95?
  # h = ES.h(1, 0.95)
  # Given our sample sizes and a 1% significance level, what size XvY effect can we detect with 95% power?
  # Compute power of test, or determine parameters to obtain target power
  # pwr.2p2n.test(n1=1392, n2=1384, h = ES.h(p1 = 1, p2 = 0.99), sig.level = 0.01)
  
  # Generate possible proportions of difference between X and Y
  make.proportions = function(n1, n2){
    props = seq(0.95, 1.0, 0.0005)
    # Make h values from the proportions
    combs = ES.h(p1=1, p2=props)
    test.result = pwr.2p2n.test(n1=n1, n2=n2, h = combs, sig.level = 0.01)$power
    bound.data = as.data.frame(cbind(props, test.result, n1, n2))
    colnames(bound.data) = c("Proportion", "Power", "n1", "n2")
    bound.data
  }
  
  # Generate ranges of sample sizes to compare
  nrange = seq(50, 2000, 50)
  all.n = expand.grid(nrange, nrange)
  all.results = do.call(rbind, mapply(make.proportions, all.n$Var1, all.n$Var2, SIMPLIFY = F))
  
  # Plot the sample size matrixes for the desired proportions
  make.plot = function(proportion){
    filt = all.results %>% dplyr::filter(Proportion==proportion)
    ggplot(filt, aes(x=n1, y=n2, fill=Power>0.8))+
      geom_raster()+
      xlab("Sample 1 size")+ylab("Sample 2 size")+
      ggtitle(paste("Difference of\n", proportion,"versus 1 at p=0.01"))+
      theme(legend.position = "none")
  }
  props = c(0.95, 0.99, 0.995)
  plot.list = lapply(props, make.plot)
  do.call(plot_grid, c(plot.list, align='v', nrow=1, ncol = 3))
  ggsave.row(paste0("Proportion power"))
  
  # Make the sample size power plots for our actual datasets
  calculate.actual.sample.power = function(s){
    cat("Sample ", s, "\n")
    filt = full.data %>% filter(Sample==s)
    nX = nrow(filt %>% filter(Chr=="X"))
    nY = nrow(filt %>% filter(Chr=="Y"))
    props = seq(0.90, 1.0, 0.0005)
    # Make h values from the proportions
    combs = ES.h(p1=1, p2=props)
    test.result = pwr.2p2n.test(n1=nX, n2=nY, h = combs, sig.level = 0.01)$power
    bound.data = as.data.frame(cbind(props, test.result, nX, nY, s, levels(droplevels(unique(filt$Strain)))), stringsAsFactors = F)
    colnames(bound.data) = c("Proportion", "Power", "nX", "nY", "Sample", "Strain")
    bound.data$Power = as.numeric(bound.data$Power)
    bound.data$Proportion = as.numeric(bound.data$Proportion)
    bound.data$nX = as.numeric(bound.data$nX)
    bound.data$nY = as.numeric(bound.data$nY)
    bound.data
  }
  
  actual.sample.results = do.call(rbind, mapply(calculate.actual.sample.power, unique(full.data$Sample), SIMPLIFY = F))
  
  # Plot the actual power of our datasets to detect changes at p=0.01
  ggplot(actual.sample.results, aes(x=Proportion, y=Power, group=Sample, col=Strain))+
    geom_hline(yintercept = 0.8)+
    geom_line()+
    xlab("Y area as proportion of X area")+
    ylab("Power of test at p=0.01")+
    facet_grid(Strain~.)+
    theme_classic()+
    theme(legend.position = "none")
  
  ggsave.double("Actual sample power")
  
  # Calculate the power by strain for our actual datasets
  calculate.actual.strain.power = function(s){
    cat("Strain ", s, "\n")
    filt = full.data %>% filter(Strain==s)
    nX = nrow(filt %>% filter(Chr=="X"))
    nY = nrow(filt %>% filter(Chr=="Y"))
    props = seq(0.90, 1.0, 0.0005)
    # Make h values from the proportions
    combs = ES.h(p1=1, p2=props)
    test.result = pwr.2p2n.test(n1=nX, n2=nY, h = combs, sig.level = 0.01)$power
    bound.data = as.data.frame(cbind(props, test.result, nX, nY, s), stringsAsFactors = F)
    colnames(bound.data) = c("Proportion", "Power", "nX", "nY", "Strain")
    bound.data$Power = as.numeric(bound.data$Power)
    bound.data$Proportion = as.numeric(bound.data$Proportion)
    bound.data$nX = as.numeric(bound.data$nX)
    bound.data$nY = as.numeric(bound.data$nY)
    bound.data
  }
  
  actual.strain.results = do.call(rbind, mapply(calculate.actual.strain.power, levels(unique(full.data$Strain)), SIMPLIFY = F))
  
  ggplot(actual.strain.results, aes(x=Proportion, y=Power, col=Strain))+
    geom_hline(yintercept = 0.9)+
    geom_vline(xintercept = 0.995)+
    geom_line()+
    xlab("Y area as proportion of X area")+
    ylab("Power of test at p=0.01")+
    theme_classic()
  
  ggsave.single("Actual strain power")
  
}
# calculate.powers()

# How consistent are clustering methods?
# What clustering method should we use for the full cluster?
cluster.subset = function(sample.data, sample.name){
  cat("Testing clustering on", sample.name, "\n")
  # Find the best clustering method for the data
  choose.method = function(){
    # We need a subset of the input, otherwise it will be way too big to plot
    set.seed(42)
    steps = seq(0,99,5) # choose every 5th angle
    input.data = sample.data %>% dplyr::select(one_of(paste0("Angle_profile_", steps)))
    subset.data = input.data[sample(nrow(input.data), size=2000),]
    
    # Agglomerative
    # Using agnes over stats::hclust because we get the agglomerative coefficient
    # Values near 1 = strong cluster structure
    
    
    
    # which method gives the best clustering?
    cat("Testing cluster methods\n")
    m = c( "average", "single", "complete", "ward")
    names(m) = m
    ac = function(x)  cluster::agnes(subset.data, method = x)$ac
    vals = map_dbl(m, ac) 
    best.method = names(which.max(vals))
    cat("Selected method", best.method, "\n")
    
    # Cluster with the best method
    hc.a = cluster::agnes(subset.data, method = best.method)
    dend.agg = as.dendrogram(hc.a)
    
    
    # Cluster with divisive method
    cat("Clustering with diana\n")
    hc.d = cluster::diana(subset.data)
    dend.diana = as.dendrogram(hc.d)
    
    # Create a tanglegram to compare the clusters
    cat("Making chart\n")
    png(filename = paste0(sample.name, " agglomerative vs divisive.png"), res=300, width=170*1.25, height=170*1.25, units="mm")
    tanglegram(dend.agg, dend.diana, sort=T, lwd=1.5, edge.lwd=1.5, main_left=paste0(best.method,", ac=", format(hc.a$ac, digits=2)), main_right=paste0("Diana, dc=", format(hc.d$dc, digits=2)))
    dev.off()
    best.method
  }
  cluster.method = choose.method()
  
  # Find the best profile offset for 5% sampling
  choose.offset = function(){
    # We need a subset of the input, otherwise it will be way too big to plot
    set.seed(42)
    offsets = seq(0,4,1)
    
    cluster.offset = function(offset){
      steps = seq(offset,99,5) # choose every 5th angle
      input.data = sample.data %>% ungroup() %>% dplyr::select(one_of(paste0("Angle_profile_", steps)))
      subset.data = input.data[sample(nrow(input.data), size=500),]
      hc = cluster::agnes(subset.data, method = cluster.method)
      cat("Offset", offset, ": AC=", hc$ac, "\n")
      hc
    }
    
    clusters = lapply(offsets, cluster.offset)
    
    offset.zero = as.dendrogram(clusters[[1]])
    
    offsets = seq(2,5,1)
    make.comparison.plot = function(offset){
      dend = as.dendrogram(clusters[[offset]])
      png(filename = paste0(sample.name, " cluster test offset ", offset,".png"), res=300, width=170*1.25, height=170*1.25, units="mm")
      tanglegram(dend, offset.zero, lwd=1.5, edge.lwd=1.5, main_left="Offset 1", main_right=paste0("Offset ", offset),
                 common_subtrees_color_branches=T)
      dev.off()
    }
    
    
    lapply(offsets, make.comparison.plot)
    
    # Make a clustering from the full profile
    steps = seq(0,99,1) # choose every 5th angle
    input.data = sample.data %>% ungroup() %>% dplyr::select(one_of(paste0("Angle_profile_", steps)))
    subset.data = input.data[sample(nrow(input.data), size=500),]
    full.dend = cluster::agnes(subset.data, method = cluster.method) %>% as.dendrogram()
    
    
    make.bk.plots = function(){
      cat("Comparing to full profile\n")
      make.full.bk.plot = function(offset){
        dend = as.dendrogram(clusters[[offset]])
        corr = cor_bakers_gamma(full.dend, dend)
        cat("Offset", offset, ", Baker's gamma = ", corr, "\n")
        corr = cor_cophenetic(full.dend, dend)  
        cat("Offset", offset, ", Cophenetic = ", corr, "\n")
        Bk_plot(full.dend, dend, main=paste0("Offset ", offset), k=c(1:50), xlim=c(0,50), col_line_Bk = "blue", col_line_asymptotic = "green")
      }
      png(filename = paste0(sample.name, " Bk plots.png"), res=300, width=170*1.25, height=170*1.25, units="mm")
      par(mfrow=c(2,3))
      lapply(seq(1,5,1), make.full.plot)
      dev.off()
      
      cat("Comparing to Offset 1\n")
      make.offsets.bk.plot = function(offset){
        dend = as.dendrogram(clusters[[offset]])
        corr = cor_bakers_gamma(offset.zero, dend)
        cat("Offset", offset, ", Baker's gamma = ", corr, "\n")
        corr = cor_cophenetic(offset.zero, dend)  
        cat("Offset", offset, ", Cophenetic = ", corr, "\n")
        Bk_plot(offset.zero, dend, main=paste0("Offset ", offset), k=c(1:50), xlim=c(0,50), col_line_Bk = "blue", col_line_asymptotic = "green")
      }
      png(filename = paste0(sample.name, " Bk offset plots.png"), res=300, width=170*1.25, height=170*1.25, units="mm")
      par(mfrow=c(2,2))
      lapply(seq(2,5,1), make.offsets.bk.plot)
      dev.off()
    }
    make.bk.plots()
  }
  choose.offset()
}


# Given the results of the subset clustering, we want to 
# compare genotype proportion to cluster using the full
# profile sampling
cluster.all = function(sample.data, sample.name){
  cluster.method = "ward"
  make.cluster = function(){
    
    input.data = sample.data %>% dplyr::select(one_of(paste0("Angle_profile_", seq(0,99,5))))
    hc = cluster::agnes(input.data,method = cluster.method)
    cat("Full data using", cluster.method, "AC=", format(hc$ac, digits=2), "\n")
    hc
  }
  hc = make.cluster()
  dend = as.dendrogram(hc)
  colours = as.numeric(as.factor(sample.data$Strain))
  colours = colours[order.dendrogram(dend)]
  labels_colors(dend) = colours
  
  png(filename = paste0(sample.name, " full clustering.png"), res=300, width=170*1.25, height=170*1.25, units="mm")
  plot(dend, main = paste0(sample.name, " full clustering, AC = ", format(hc$ac, digits=2)))
  dev.off()
}
# cluster.all(MF1.filt, "MF1 data")

# test the tSNE paramaters
run.tSNE.parameter.tests = function(sample.data, sample.name){
  library(Rtsne)
  # set.seed(42)
  sample.data = take.sample(sample.data, 1000)
  input = sample.data %>% ungroup() %>% dplyr::select(one_of(paste0("Angle_profile_", seq(0,99,1))))
  sample.data$Strain = as.factor(sample.data$Strain)
  rownames(input) = paste0(sample.data$Strain, "_", sample.data$CellID)
  
  run.tsne.perp = function(perplexity, iterations){
    cat("Running tSNE with perplexity", perplexity, "for", iterations, "iterations\n")
    # run Rtsne with default parameters
    set.seed(42)
    rtsne_out = Rtsne(as.matrix(input), perplexity=perplexity, max_iter=iterations)
    values    = as.data.frame(rtsne_out$Y)
    
    ggplot(values, aes(x=V1, y=V2, col=sample.data$Strain))+
      geom_point()+
      ggtitle(paste0(perplexity, "x", iterations))+
      labs(col="Strain")+
      theme_void()+theme(legend.position = "none")
  }
  
  perps = seq(50, 300, 50)
  iters = seq(1000, 5000, 500)
  vars = expand.grid(perps, iters)
  
  plot.list = mapply(run.tsne.perp, vars$Var1, vars$Var2, SIMPLIFY = F)
  do.call(plot_grid, c(plot.list, align='hv', nrow = 9, ncol = 6))
  ggsave(filename=paste0(sample.name, " profile tSNE perplexity test.png"), width = 170*1.25, height=250*1.25, units ="mm" )
}

# Statistical testing of X versus Y distributions
run.statistical.tests = function(){
  cat("Running statistical tests\n")
  # Are the data normally distributed?
  normality.test = function(){
    
    # visual test
    ggqqplot(full.data$Area_square_microns)
    ggsave.single("Full data area qq plot")
    
    # Shapiro-Wilk test of normality
    sw.data = with(full.data, tapply(Area_square_microns, Dataset, shapiro.test))
    
    get.p.val = function(sw) sw$p.value
    pvals = unlist(lapply(sw.data, get.p.val))
    # if(min(pvals)<0.05) cat("Normality cannot be assumed\n")
    # cat("P values - min:", min(pvals), ", max", max(pvals), "\n")
  }
  normality.test()
  
  # Since the data are not normal, use a permutation test
  # to control for multiple factors
  permutation.test = function(){
    # There is a strong effect of strain, influenced by fluor
    # Suggests that there is a clear deviation from independence
    ind.result = independence_test(Area_square_microns ~ Strain | Chr,
                                   data = full.data, teststat="quad")
    statistic(ind.result, type="standardized")
    qperm(ind.result, 0.95)
    coin::pvalue(ind.result, method = "single-step")
  }
  permutation.test()
  
  # Non-parametric equivalent to a two-way ANOVA - Scheirer-Ray-Hare test
  # Is there an interaction between chromosome and strain to determine area?
  srh.test = function(){
    scheirerRayHare(Area_square_microns ~ Chr*Strain,
                    data = full.data)
    
    scheirerRayHare(Area_square_microns ~ Chr*Length_seg_0_microns,
                    data = full.data)
  }
  srh.test()
}

calculateMagnitudeDifferencesOfParameters = function(sample.data, sample.name, avg.method){
  
  # Calculate the global WT average and sem
  wt.global = sample.data %>% dplyr::group_by(Bkg, Type) %>%
    dplyr::filter(Type=="WT") %>%
    dplyr::select(-contains("pixels"), -contains("profile"), -ends_with("_start"), -ends_with("_end")) %>% # get just the columns of interest
    dplyr::summarise_if(is.numeric, .funs=c("average" = avg.method, "sem" = sem)) %>%
    tidyr::gather("Variable", "Value", -Type, -Bkg) %>% # Gather the variables
    dplyr::mutate(Func = str_extract(Variable, "(average|sem)"), # Record which function generated the row value
                  Variable = gsub("(_average|_sem)", "", Variable)) %>% # Rename the variables to remove trailing function names
    tidyr::unite(Temp, Type, Func) %>% # combine the type and function 
    tidyr::spread(Temp, Value) %>% # and spread values by type and function
    tidyr::gather(Type, Value, -Variable, -Bkg) %>% # gather back, retaining the WT values as columns
    dplyr::mutate(Func = str_extract(Type, "(average|sem)"),
                  Type = gsub("(_average|_sem)", "", Type)) %>%
    tidyr::spread(Func, Value) %>%
    dplyr::mutate(WT_global_avg = average, WT_global_sem = sem) %>%
    dplyr::select(-average, -sem, -Type)
  
  # Calulate average of each parameter by strain, and find the 
  # difference to wild type
  median.data = sample.data %>% dplyr::group_by(Bkg, Type, Chr) %>%
    dplyr::select(-contains("pixels"), 
                  -contains("profile"),
                  -contains("Width_of_body"),
                  -contains("Length_of_hook"),
                  -ends_with("_start"), 
                  -ends_with("_end")) %>% # get just the columns of interest
    dplyr::summarise_if(is.numeric, .funs=c("average" = avg.method, "sem" = sem)) %>% # calculate the average of the groups
    dplyr::group_by(Bkg, Type, Chr) %>% # make the groups we want for strain averaging
    tidyr::gather("Variable", "Value", -Type, -Chr, -Bkg) %>% # Gather the variables
    dplyr::mutate(Func = str_extract(Variable, "(average|sem)"), # Record which function generated the row value
                  Variable = gsub("(_average|_sem)", "", Variable)) %>% # Rename the variables to remove trailing function names
    tidyr::unite(Temp, Type, Func) %>% # combine the type and function 
    tidyr::spread(Temp, Value) %>% # and spread values by type and function
    tidyr::gather(Type, Value, -Chr, -Variable, -Bkg) %>% # gather back, retaining the WT values as columns
    dplyr::mutate(Func = str_extract(Type, "(average|sem)"),
                  Type = gsub("(_average|_sem)", "", Type)) %>% # Recreate the type and function columns
    merge(., wt.global, by=c("Variable", "Bkg")) %>%
    tidyr::spread(Func, Value) %>% # we now have all the strains aligned to WT
    dplyr::mutate(FractAvg = average/WT_global_avg) %>% # calculate the fractional mean
    dplyr::mutate(SemWTNorm = WT_global_sem/WT_global_avg, SemNorm = sem/average) %>% # calculate the normalised sem
    dplyr::mutate(SemBoth = sqrt((SemWTNorm^2+SemNorm^2))) %>% # calcluate the combined sem
    dplyr::mutate(SemUpper = FractAvg+SemBoth, SemLower = FractAvg-SemBoth) %>% # make upper and lower bounds for plotting sem
    na.omit # remove any rows with missing values - ie MF1 shSLY does not exist
  
  # Now test if each parameter XvY is significant by Mann-Whitney U test
  # and bind the results to the mean values
  run.strain.test = function(strain){
    bkg  = (sample.data %>% dplyr::filter(Strain==strain) %>% dplyr::distinct(Bkg))$Bkg
    type = (sample.data %>% dplyr::filter(Strain==strain) %>% dplyr::distinct(Type))$Type
    run.mann.whitney.test = function(variable){
      x.data = sample.data %>% dplyr::filter(Strain==strain, Chr=="X")
      y.data = sample.data %>% dplyr::filter(Strain==strain, Chr=="Y")
      wilcox.test(x.data[[variable]], y.data[[variable]], paired=FALSE)$p.value
    }
    pvalues = do.call(rbind, lapply(unique(median.data$Variable), run.mann.whitney.test))
    res = cbind("Variable" = unique(median.data$Variable), "pval"= pvalues)
    res = as.data.frame(res)
    res$Type = type
    res$Bkg = bkg
    res
  }
  comb = do.call(rbind, lapply(unique(sample.data$Strain), run.strain.test))
  comb$V2 = as.numeric(as.character(comb$V2))
  median.data = merge(median.data, comb, by = c("Variable", "Type", "Bkg"))
  
  # Adjust the pvalues for Bonferroni correction and determine if significant
  median.data$Significant = median.data$V2*(nrow(median.data)/2)<0.01
  
  # Clean variable names
  median.data = median.data %>% dplyr::mutate(Variable = gsub("_((square_)?microns|degrees)", "", Variable),
                                              Variable = gsub("_", " ", Variable)) %>% # polish the variable names
    dplyr::filter(Variable!="Elongation", # Something weird here.
                  Variable!="Angle between reference points") %>%  # Not really useful, and is a long name
    dplyr::filter(Variable %in% c("Area", "Difference from median", "Bounding height", "Bounding width"))
  
  # Create a combined variable to use for shape determination
  median.data$Shape = paste(median.data$Chr, median.data$Significant)
  
  # Set variable order for plotting
  # We want the size variables separate from the size-independent variables
  variables = as.data.frame(unique(median.data$Variable))
  colnames(variables) = c("Variable")
  variables = variables %>%
    dplyr::mutate(Group = case_when( grepl("([L|l]ength|[W|w]idth|height)", Variable) ~ "Size-dependent",
                                     grepl("diameter|feret", Variable) ~ "Size-dependent",
                                     grepl("Area|Perimeter", Variable) ~ "Size-dependent",
                                     TRUE~"Size-independent"))
  
  variable.order = variables %>% dplyr::arrange(desc(Group), desc(Variable)) %>% dplyr::select(Variable)
  median.data$Variable = with(median.data, factor(Variable, levels = variable.order$Variable))
  
  median.data$Type = with(median.data, factor(Type, levels = c("WT", "Yqdel", "shSLY")))
  median.data$Bkg = with(median.data, factor(Bkg, levels = c("MF1", "C57")))
  
  # Set colours - want to use Set2, but change the order
  fill.colours = RColorBrewer::brewer.pal(length(unique(sample.data$Type)), "Set2")
  temp = fill.colours[2]
  fill.colours[2] = fill.colours[3]
  fill.colours[3] = temp
  
  # Plot
  ggplot(median.data, aes(x=Variable, y=FractAvg, col=Type, shape=Shape))+
    geom_hline(yintercept = 1)+
    geom_rect(xmin=5.2, xmax=5.8, ymin=-Inf, ymax=Inf, fill="white", col=NA)+ # break between size and shape variables
    geom_pointrange(aes(x=Variable,ymin=SemLower, ymax=SemUpper),
                    size=0.7, position = position_jitterdodge(seed=1, jitter.width=0.5, dodge.width=0.05))+
    coord_flip()+
    scale_y_continuous(breaks=seq(0.5,2,0.1))+
    scale_shape_manual(values=c(1, 16, 2, 17))+ # filled and hollow shapes
    scale_colour_manual(values=fill.colours)+
    facet_grid(.~Bkg)+
    ylab(paste("Strain chromosome",avg.method,"/ WT global",avg.method))+
    labs(shape="X/Y", col="Strain")+
    guides(shape=F)+
    theme(axis.text.y = element_text(size=10), 
          axis.text.x = element_text(size=10),
          axis.title.y = element_blank(),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey"),
          legend.position = "top")
  
  
  ggsave.custom(paste0(sample.name, " ", avg.method , " magnitude differences to WT"), width=170, height=90)
  
}

