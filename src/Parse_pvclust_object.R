library("pvclust")

basedir <- dirname(sys.frame(1)$ofile)
basedir = dirname(basedir)
print(basedir)

coexpression_dir = file.path(basedir, "data", "AC049_TCGA_KNN_imputed_min_ecdna_3_matrix", "co_expression")
dir.create(coexpression_dir, showWarnings = FALSE)

pvclust_fn = file.path(basedir, "data", "pvclust", "pvclust.wardD2.TCGA_expression.ecDNA3.RData")
load(pvclust_fn)

edge_info = pvclust.TCGA.genes$edges

au95_clusters = pvpick(pvclust.TCGA.genes, alpha=0.95, pv="au", type="geq", max.only=TRUE)
au_95_df_fn = file.path(coexpression_dir, "au_95.txt")
au_95_plot_fn = file.path(coexpression_dir, "au_95.pdf")

edge_list = c()
edge_au_list = c()
edge_si_list = c()
cluster_list = c() 
for (i in 1:length(au95_clusters$edges)) {
  edge = au95_clusters$edges[i]
  edge_au = edge_info[edge, "au"]
  edge_si = edge_info[edge, "si"]
  cluster_members = au95_clusters$clusters[i][[1]]
  for (j in 1:length(cluster_members)) {
    gene = cluster_members[j]
    
    edge_list <- c(edge_list, edge)
    cluster_list <- c(cluster_list, gene)
    edge_au_list <- c(edge_au_list, edge_au)
    edge_si_list <- c(edge_si_list, edge_si)
    
  }
}

au_95_df = data.frame("gene" = cluster_list,
                      "edge" = edge_list,
                      "edge_au" = edge_au_list,
                      "edge_si" = edge_si_list)

write.table(au_95_df, file = au_95_df_fn, row.names=FALSE, col.names=TRUE, sep="\t", quote = FALSE)

pdf(file=file.path(au_95_plot_fn), width=100, height=50)
pvclust_plot <- plot(pvclust.TCGA.genes,  print.pv=c("au"), hang = -1, cex = 0.8)
pvrect(pvclust.TCGA.genes, alpha = 0.95)
dev.off()
