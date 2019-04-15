options(stringsAsFactors = F)
rm(list=ls())
library(pacman)
p_load(cmapR)
p_load(dplyr)


## import gct files
gct.str <- dir('data', full.names = T, pattern='\\.gct$')
gct <- lapply(gct.str, parse.gctx)
names(gct) <- c('5_acK', '1_CNA', '4_pSTY', '3_Protein', '2_RNAseq')


## rdescr
gct.rdesc <- lapply(gct, function(x)x@rdesc)

## expression
gct.expr <- lapply(gct, function(x) x@mat)

## column ids
gct.cid <-lapply(gct, function(x) x@cid)

## row ids
gct.rid <-lapply(gct, function(x) x@rid)


## common sample names
samp.common <- Reduce(f = intersect, gct.cid)
gct.expr <- lapply(gct.expr, function(x)x[, samp.common])

###################
## rdesc 
gct.gene <- lapply(gct.rdesc, function(x) try(x[ ,grep('geneSymbol', colnames(x), value=T, ignore.case = T)[1]], silent=T))
#gct.gene <- lapply(gct.gene.id, )
if(sum(sapply(gct.gene, length) < 2) > 0){
  idx <- which( sapply(gct.gene, length) < 2 )
  for(i in idx)
    gct.gene[[i]] <- gct.rid[[i]]
}

## data type
gct.data.type <- lapply(names(gct.gene), function(x) rep(x, nrow(gct.expr[[x]])))
names(gct.data.type) <- names(gct.gene)

## ids
gct.id <- lapply(gct.rdesc, function(x)x[, 'id'])
#gct.id <- lapply(gct.rdesc, function(x) try(grep('id', colnames(x), value=T), silent=T))

if(sum(sapply(gct.id, length) < 2) > 0){
  idx <- which( sapply(gct.id, length) < 2 )
  for(i in idx)
    gct.id[[i]] <- gct.rid[[i]]
}

## combine
gct.rdesc <- lapply(names(gct.gene), function(x) data.frame(ID=gct.id[[x]], geneSymbol=gct.gene[[x]], DataType=gct.data.type[[x]]))

rdesc <- Reduce(f = rbind, gct.rdesc)

###############################
## cdesc from proteome
cdesc <- gct[['3_Protein']]@cdesc

## add nmf
nmf <- read.delim('data/clin_anno_nmf.txt')
nmf <- nmf[, c('Sample.ID', 'NMF.consensus')]
nmf$NMF.consensus <- paste0('C', nmf$NMF.consensus)
cdesc <- left_join(cdesc, nmf, 'Sample.ID')
rownames(cdesc) <- cdesc$Sample.ID

#################################
## single data frames
column.anno <- cdesc[samp.common, ]
row.anno <- rdesc
tab.expr.all <- Reduce(f=rbind, gct.expr)
rownames(tab.expr.all) <- make.unique(rownames(tab.expr.all)) 
rownames(row.anno) <- rownames(tab.expr.all)

## reorder to PAM50
ord.idx <- order(column.anno$PAM50)
column.anno <- column.anno[ord.idx,]
tab.expr.all <- tab.expr.all[, ord.idx]

## reorder rows
row.idx <- with(row.anno, order(geneSymbol, DataType))
tab.expr.all <- tab.expr.all[row.idx, ]
row.anno <- row.anno[row.idx, ]

## export
save(column.anno, row.anno, tab.expr.all, file = 'data/data_bc2019.RData')

