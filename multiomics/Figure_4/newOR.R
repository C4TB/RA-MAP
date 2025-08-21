# RAMap: new attempt at OR analy using generalised clusterprofiler
# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html -> ch12 generalised 
# additional parameters TERM2GENE and (opt) TERM2NAME. 
# TERM2GENE is a data.frame with first column of term ID and second column of corresponding mapped gene
# can give universe, but might be easier to set up TERM2GENE with only desired genes, as must make it anyway

# july 2023 - with msig collections
# jun 2024 - chaussabel mods
library(clusterProfiler)
library(dplyr)
library(readr)
library(ggplot2)
library(enrichplot)
#library(pathview)
# limma for strsplit2
library(limma)
library(annotate)
library(org.Hs.eg.db)
library(viridis)
#library(AnnotationDbi)
#library(KEGG.db)

#library(ReactomePA)


# prev was biomart call here

# new - read set of files like univ_file="../../annot_data/WHL_sel_ann_comp.txt"
# univ_annot_dir set in filenames.R 
source("filenames.R")


####
#from Gostats.R 

# main prog
# usage: GOStats.R scaled_metabolomics  sc_met_gostats
# runs on all files *.node in scaled_metabolomics/

# commandA=commandArgs(trailingOnly = TRUE)
# debug

based="/data/ngrs2/RA-MAP/RAMap_Feb21/all_platforms_104_patients_0m/analy/"
based="/data/ngrs2/RA-MAP/RAMap_Feb21/all_platforms_104_patients_6m/analy/"
inclusdir="clusters_p1e-10_0.50"
clustype = "highMI"
inclusdir="neighbours_of_clinical/clusters_all/"
clustype = "NoC"
maxcl = 6
h = 5; w = 8; 
if (FALSE) {
inclusdir="clusters_p1e-10_0.50"
clustype = "highMI"
maxcl = 25
h=10; w=14
}

commandA=c(file.path(based, inclusdir), "clustest", "0m", "Reactome")
commandA=c(file.path(based, inclusdir), "clustest", "0m", "ChaussabelModules")
commandA=c(file.path(based, inclusdir), "clustest", "6m", "ChaussabelModules")
# cat ("running in DEBUG MODE (fixed input and output dirs)/n")
run_direc=commandA[1]
outdir=commandA[2]
timept=commandA[3]
msigtype=commandA[4]

dir.create(outdir)

# skip if fewer than this nodes in file
min_lines = 12

t2gf = paste0("../../annot_data/", msigtype, "_multiomic_term_list.tsv")
t2g = read_tsv(t2gf)

files = list.files(path=run_direc, pattern="*.node")
#networks = c("CD14.log.wide.filtered", "CD14.DAS28CRP", "CD14.SDAI",
#             "CD4.log.wide.filtered", "CD4.DAS28CRP", "CD4.SDAI",
#             "CD8.log.wide.filtered", "CD8.DAS28CRP", "CD8.SDAI",
#            "PBMCs.log.wide.filtered", "PBMCs.DAS28CRP", "PBMCs.SDAI",
#            "WHL.log.wide.filtered", "WHL.DAS28CRP", "WHL.SDAI")

if (FALSE) {
  # run tests on each cluster for all networks
  # test get_enriched_terms("c_1.node")
  ct = 0
  xl = list()
  for(i in files) {
    #  if (! file.exists(file.path(outdir, paste0(i, ".hyperGTest.GO.summary")))){
    ng = read.table(file.path(run_direc,i), stringsAsFactors = FALSE)
    
    if (nrow(ng) < min_lines) {
      cat("network",i,"too short, not analysing\n")
    } else {
      cat("doing",i,"containing",nrow(ng),"\n")
      x <- enricher(ng$V1, TERM2GENE = t2g[,c(1,2)])
      xl[[i]] <- x
      ct=ct+1
      # got here - for debugging
      #         if (ct>3)
      #      	 stop()
    }
    #  }
  }
  cat("did", ct, "out of", length(files), "networks\n") 
  library(viridis) # #... + scale_color_viridis() 
  # tried  + scale_color_brewer(type = "qual", palette = 1) but gives error - see help - 
  # for contin scale use "distiller" version
  
  p <- dotplot(xl[["call.node"]], showCategory=30) + 
    ggtitle(paste("dotplot for c_all")) 
  # scale_color_viridis(option = "H", direction = -1)
  show(p)
}

# ch14 of book: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html

gl = list()
for(i in 1:maxcl) {
  #  if (! file.exists(file.path(outdir, paste0(i, ".hyperGTest.GO.summary")))){
  key = paste0("c_", i)
  file = paste0("c_", i, ".node")
  ng = read.table(file.path(run_direc,file), stringsAsFactors = FALSE)
  
  if (nrow(ng) < min_lines) {
    cat("network",i,"too short, not analysing\n")
  } else {
    gl[[key]] = ng$V1

  }
  #  }
}
# remove HALLMARK_ at the start, etc 
str_clip = paste0(msigtype, "_")
t2g_to_use = t2g[,c(1,2)] %>% mutate(gs_name = gsub(str_clip, "", gs_name, ignore.case = TRUE))


# come back to this: 
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html?q=compareCluster#formula-interface-of-comparecluster

ck <- compareCluster(geneCluster = gl, fun = function(x) {enricher(x, TERM2GENE = t2g_to_use)} )
# can add values = logscale but not really helpful logscale = scales::rescale(log(seq(0.1,exp(1),0.1)))
# label_format is for string wrapping. Note can be function - see FAQ
title_str = paste(timept, "dotplot for c1:c",maxcl, ", ", msigtype, "gene sets")
# , showCategory = 20 OK for some but makes GOBP too busy
p <- dotplot(ck, label_format = 70) + ggtitle(title_str) + 
#  scale_color_viridis(option = "magma") + geom_point(shape = 1, stroke = 1)
  scale_color_distiller(palette = "Blues", direction = -1) + geom_point(shape = 1, stroke = 1, color = 1)
show(p)
outf <- file.path(outdir, paste0(clustype, "_", msigtype, "_", timept,".png"))
ggsave(outf, p, w=w, h=h)
table(as.data.frame(ck)$Cluster)

# shape types to consider 1 or 16 or 21 
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html note example has cowplot grid example 
  # see ch 14 - biological theme comparison - for comparing clusters

# or made a plot list by looping, then cowplot::plot_grid(plotlist=pl, nrow = 3)
