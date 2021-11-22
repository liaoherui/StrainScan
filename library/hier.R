x <- read.table("distance_matrix_rebuild.txt",header=T, row.names=1)
d <- as.dist(as(x,"matrix"))
#hc <- hclust(d,method="single")
hc <- hclust(d,method="complete")
res <- sort(cutree(hc,h=0.01))
res

