
library(ComplexHeatmap)
library(RColorBrewer)
library(gridtext)

# run on 2 lists of paths-to-cluster-files
# Rscript cluster_jaccard.R l1.txt l2.txt row_title col_title output_png
# l1.txt applied to rows, l2 to cols 
# l1.txt contains lines like ./high_SDAI_clusters/clusters_nocut_p1e-10_0.50/c_1.node

# for ((i=1;i<=30;i++)); do echo high_SDAI_clusters/clusters_nocut_p1e-10_0.50/c_$i.node >> t1.txt; done
# for ((i=1;i<=30;i++)); do echo low_SDAI_clusters/clusters_nocut_p1e-10_0.50/c_$i.node >> t2.txt; done



argv=commandArgs(trailingOnly = TRUE)
##argv=c("t1.txt","t2.txt", "high_SDAI_clusters", "low_SDAI_clusters","jaccard_matrix.png")

flist1=scan(argv[1],what=character())
flist2=scan(argv[2],what=character())
row_title=argv[3]
column_title=argv[4]
outfile=argv[5]

quiet=TRUE
nl1=list()
for (inf in flist1) {
    cid=sub(".node","",basename(inf))
    nodes=scan(inf,what=character(),quiet=quiet)
    nl1[[cid]]=nodes
    if (! quiet)
       cat("nl1: added",length(nodes),"nodes to data for",cid,"\n")
}
nl2=list()
for (inf in flist2) {
    cid=sub(".node","",basename(inf))
    nodes=scan(inf,what=character(),quiet=quiet)
    nl2[[cid]]=nodes
    if (! quiet)
       cat("nl2: added",length(nodes),"nodes to data for",cid,"\n")
}

nr=length(nl1)
nc=length(nl2)
dm=matrix(rep(0.0,nr*nc),nrow=nr,ncol=nc)
rownames(dm)=names(nl1)
colnames(dm)=names(nl2)
print(dim(dm))
for (i in 1:nr){
    for (j in 1:nc){
    # length(intersect(nl1[["c_1"]],nl2[["c_2"]]))/length(union(nl1[["c_1"]],nl2[["c_2"]]))
    # even with numerical index need double brackets to unlist it 
    dm[i,j]=length(intersect(nl1[[i]],nl2[[j]]))/length(union(nl1[[i]],nl2[[j]]))
    }
}
# rev(heat.colors(40)),
mypal = c("#FFFFFF", brewer.pal(9, "YlOrBr")[4:9])
hm = Heatmap(dm, cluster_rows = FALSE, cluster_columns = FALSE,
             row_title_gp = gpar(fontsize = 18), column_title_gp = gpar(fontsize = 18), column_title_side = "bottom",
             row_names_side = "left", column_names_rot = 45, 
	     row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16),
             name = "jaccsim",
	     heatmap_legend_param = list(title = gt_render("Jaccard\nsimilarity", padding = unit(c(4,0,4,0), "pt")),
	                                 title_gp = gpar(fontsize = 16, fontface = "bold"), 
	                                 legend_height = unit(2.5, "cm"), labels_gp = gpar(fontsize = 16)),
	     col = mypal, 
	     rect_gp = gpar(col = "black", lwd = 1),
 	     row_title = row_title, column_title = column_title)
show(hm)
png(outfile,h=20*nr,w=20*nc)
show(hm)
dev.off()
cat("made",outfile,"\n")
