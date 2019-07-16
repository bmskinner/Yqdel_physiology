# Generic functions
# Loads these into the global environment
load.generic.functions = function(){
  
  is.installed <<- function(pkg) pkg %in% rownames(installed.packages())
  
  install.missing <<- function(packages, biopackages=c(), github=c(), repos = "https://cran.ma.imperial.ac.uk/" ) {

    install.or.load = function(pkg){
      if(!is.installed(pkg)){
        install.packages(pkg, dependencies = TRUE, repos = repos)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
    }

    install.or.load.bioconductor = function(pkg){
      install.or.load("BiocManager") # Bootstrap installer
      if(!is.installed(pkg)) {
        BiocManager::install(pkg)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
    }
    
    install.or.load.github = function(pkg){
      install.or.load("devtools") # Bootstrap installer
      if(!is.installed(pkg)) {
        devtools::install_github(pkg)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
    }

    tryCatch({
      lapply(packages, install.or.load)
      lapply(github, install.or.load.github)
      lapply(biopackages, install.or.load.bioconductor)
    }, error=function(e){
      cat("Error loading required libraries\n")
      message(e)
    })
  }
  
  # Utility functions for saving images
  ggsave.single <<- function(file) {
    ggsave(paste0("Report/",file, ".svg"), width = 85*1.25, height=85*1.25, units ="mm" )
    ggsave(paste0("Report/",file, ".png"), width = 85*1.25, height=85*1.25, units ="mm" )
  }
  ggsave.double <<- function(file){
    ggsave(paste0("Report/",file, ".svg"), width = 170*1.25, height=170*1.25, units ="mm" )
    ggsave(paste0("Report/",file, ".png"), width = 170*1.25, height=170*1.25, units ="mm" )
  }
  ggsave.row <<- function(file){
    ggsave(paste0("Report/",file, ".svg"), width = 170*1.25, height=85*1.25, units ="mm" )
    ggsave(paste0("Report/",file, ".png"), width = 170*1.25, height=85*1.25, units ="mm" )
  }
  ggsave.custom <<- function(file, width, height){
    ggsave(paste0("Report/",file, ".svg"), width = width*1.25, height=height*1.25, units ="mm" )
    ggsave(paste0("Report/",file, ".png"), width = width*1.25, height=height*1.25, units ="mm" )
  }
  
  write.headed.table <<- function(x, file) write.table(x, file = paste0("Report/",file), sep = "\t", row.names = F, col.names = T, quote = F)
  
  sem <<- function(x) sqrt(var(x)/length(x)) 
  
  SEP <<- function(p, n) sqrt((p*(1-p)) / n)
  
  create.strain.colours <<- function(){
    # Set colours - want to use Set2, but change the order
    # fill.colours = RColorBrewer::brewer.pal(3, "Set2")
    # temp = fill.colours[2]
    # fill.colours[2] = fill.colours[3]
    # fill.colours[3] = temp
    # fill.colours
    
    c("#057451", "#2852af", "#d74810")
  }
}

# Extensions for ggplot to allow
# e.g. split violin charts
load.gg.extensions = function(){
  GeomSplitViolin <<- ggproto("GeomSplitViolin", GeomViolin, 
                              draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                grp <- data[1, "group"]
                                newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                                
                                if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                  stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                            1))
                                  quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                  aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                  aesthetics$alpha <- rep(1, nrow(quantiles))
                                  both <- cbind(quantiles, aesthetics)
                                  quantile_grob <- GeomPath$draw_panel(both, ...)
                                  ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                }
                                else {
                                  ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                                }
                              })
  
  geom_split_violin <<- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                                 draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                                 show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
          params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
  }
}