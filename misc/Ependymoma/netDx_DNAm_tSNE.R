# plot tSNE of DNAm of ependymoma
require(Rtsne)

# -----------------------------------------------
# helper plotting function
plot_cluster=function(data, var_cluster, palette)  
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
} 
# -----------------------------------------------

mData <- "/home/shraddhapai/BaderLab/2018_Epen_DNAm/input/GSE90496_EPN_PFAB_beta.txt.gz"

orig_dat <- read.delim(mData,sep="\t",h=T,as.is=T)

# tSNE
cat("computing tSNE\n")
set.seed(9)
t0 <- Sys.time()
dtsne <- Rtsne(as.matrix(orig_dat),check_duplicates=FALSE,pca=TRUE,
	perplexity=30,theta=0.5,dims=2)
t1 <- Sys.time()
cat(sprintf("tSNE computation took %i seconds\n", t1-t0))
dat <- as.data.frame(dtsne$Y)

# kmeans
cat("kmeans clustering\n")
clust_km <- kmeans(scale(dat),2)
cl_kmeans <- factor(clust_km)

#cat("hclustering\n")
#fit_km_hclust <- hclust(dist(scale(dat))
#cl_hclust <- factor(cutree(fit_km_hclust,k=2))

pdf("tsne_DNAm.pdf")
plot_cluster(dat,"kmeans","Accent")
