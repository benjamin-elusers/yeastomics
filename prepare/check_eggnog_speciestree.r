# Trim the alignment > threshold % gaps
wexac_dir = "/media/WEXAC/EGGNOG/4751_Fungi/4751_Fungi_179sp/"

library(ape)
library(phytools)
fudir = here::here('data','eggnog','4751_Fungi')
ncbidir =  here::here('data','ncbi')
fungi_class = treeio::read.tree(file.path(fudir,"test.fungi.tree"))
fungi_class$edge.length=rep(1,length(fungi_class$edge.length))
fungi_class_ultra=chronos(fungi_class, lambda=0) |> multi2di()


fungi_concat = treeio::read.tree(file.path(fudir,"basic_sptree/cog_100-alg_concat_default-raxml_default/fungi179-19orthogroups-seqs-sorted.fa.gz.final_tree.nwx"))
fungi_concat_ultra = force.ultrametric(fungi_concat) |> chronos(lambda = 0)


#association <- cbind(sort(fungi_class$tip.label), sort(fungi_concat$tip.label))
library(dendextend)
entanglement(fungi_class_ultra,fungi_concat_ultra)
tanglegram(rank_branches(fungi_class_ultra |> as.dendrogram()), rank_branches(fungi_concat_ultra |> as.dendrogram()),
           lab.cex = 1, edge.lwd = 1,
           margin_inner = 3.5, type = "t", center = TRUE,
           dLeaf = -0.1, xlim = c(20,0), columns_width = c(5, 1, 5))



library(phytools)
tmp=cophylo(tr1 = fungi_class_ultra, tr2 =fungi_concat_ultra, rotate=T)
edge.col=1 + (tmp$trees[[1]]$tip.label == tmp$trees[[2]]$tip.label)

pdf(here::here("plots","4751_Fungi-species_tree-cophyloplot.pdf"),width = 10,height=20)
plot(tmp,link.type="curved")
legend('topleft',legend = "NCBI\ntaxonomy", bty = 'n',text.font = 2)
legend('topright',legend = "Eggnog\n(19 orthogroups)", bty = 'n',text.font = 2)
dev.off()

comparePhylo(tmp$trees[[1]],tmp$trees[[2]],plot = T)
TreeDist::TreeDistance(tmp$trees[[1]],tmp$trees[[2]])
phytools::multiRF(tmp$trees)
library(Quartet)
library(TreeDist)
x=Quartet::QuartetStatus(tmp$trees)
d=Quartet::SimilarityMetrics(x, similarity = T)

pdf(here::here("plots","4751_Fungi-species_tree-ncbi_vs_eggnog.pdf"))
Quartet::VisualizeQuartets(tmp$trees[[1]],tmp$trees[[2]], style='pie')
dev.off()
#TreeDist::VisualizeMatching(NyeSimilarity, tmp$trees[[1]],tmp$trees[[2]],TreeDist::MapTrees())
