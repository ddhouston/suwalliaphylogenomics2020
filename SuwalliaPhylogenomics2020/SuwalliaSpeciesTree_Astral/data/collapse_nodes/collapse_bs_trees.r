library(ape)
library(ips)

args <- commandArgs(trailingOnly=TRUE)

#collapse nodes with bootstrap under 10
tr <- read.tree(args[1])
tr.bs <- lapply(tr, collapseUnsupportedEdges, value="node.label", 10)
class(tr.bs)<-"multiPhylo"

tr.sp <- strsplit(args[1], ".tre")
new.st <- paste(tr.sp, "_bs10.tre", sep="")
write.tree(tr.bs, new.st)


