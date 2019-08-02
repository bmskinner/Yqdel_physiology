# Datasets and scripts for Rathje et al - Differential sperm motility mediates the sex ratio drive shaping mouse sex chromosome evolution

The measurements of pre-FISH nuclei exported from Nuclear Morphology Analysis are in Raw_data.
They should be unzipped before running the analysis.

The R script ```createCharts.R``` performs the analyses and generates charts. The other R files contain helper functions.

Supplementary tables containing the embryo experiments and breeding data are also in the Raw_data folder:

```Breeding_data.csv```

(1) Animals from the Cambridge colony (MF1 background) were used for IVF, embryo recovery, timed mating and sperm morphology experiments  
(2) Animals from the Kent colony (MF1 background) were used for sperm morphology and sperm motility experiments  
(3) Animals from the France colony (C57Bl6/N background) were used for sperm morphology experiments  
(4) IVF experiments tested whether the presence of cumulus cells affected the sex ratio from IVF using XYRIIIqdel fathers  
(5) Embryo recovery experiments tested whether sperm aging affected the sex ratio from natural mating in XYRIII and/or XYRIIIqdel fathers  
(6) Timed mating experiments tested whether maternal genotype affected the sex ratio from natural mating in XYRIII and/or XYRIIIqdel fathers  

```XgYqdel_IVF.csv```  
XgfpYqdel fathers only. IVF with/without cumulus removal using hyaluronidase, counted at blastocyst stage. IVF performed with/without hyaluronidase, embryos counted at blastocyst stage. Females were CBAB6/F1 hybrids, 8 per IVF batch.

```Two_cell_recovery.csv```  
CBAB6/F1 hybrid females superovulated (2 per male), mated for in vivo fertilisation, embryos flushed at 2 cell stage and cultured for GFP counting at blastocyst stage, and the same with MF1 females.


```Embryo_data.csv```  
XgfpYqdel and XgfpY control fathers mated to WT females vs Yqdel daughters. Males and females housed in pairs for timed mating, and females boxed out when plugged, i.e. in vivo fertilisation without superovulation. Dissection and embryo counting performed mid gestation.