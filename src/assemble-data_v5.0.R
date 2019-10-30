options(stringsAsFactors = F)
rm(list=ls())
library(pacman)
p_load(cmapR)
p_load(dplyr)
p_load(glue)

setwd('c:/Users/karsten/Dropbox/Devel/CPTAC-BRCA2019/')

## import gct files
data.dir <- 'g:/My Drive/CPTAC/prospBRCA/data/data_freeze_v5.0/'
label <- 'v5.0'

ord.column <- 'PAM50'

gct.str <- c(glue("{data.dir}prosp-brca-v5.0-acetylome-ratio-norm-NArm.gct"),
             glue("{data.dir}prosp-brca-v5.0-gene-level-cnv-gistic2-all_data_by_genes.gct"),
             glue("{data.dir}prosp-brca-v5.0-phosphoproteome-ratio-norm-NArm.gct"),
             glue("{data.dir}prosp-brca-v5.0-proteome-ratio-norm-NArm.gct"),
             glue("{data.dir}prosp-brca-v5.0-rnaseq-fpkm-log2-row-norm-median-mad.gct")
             )

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


## remove genes without HGNC symbol (CNA data)
## e.g. "7SK|ENSG00000232512.2"
gct.expr[['1_CNA']] <- gct.expr[['1_CNA']][ -grep('\\|', gct.rid[['1_CNA']]),  ]
gct.rdesc[['1_CNA']] <- gct.rdesc[['1_CNA']][ -grep('\\|', gct.rid[['1_CNA']]),  ]
gct.rid[['1_CNA']] <- gct.rid[['1_CNA']][ -grep('\\|', gct.rid[['1_CNA']]) ]


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

#################################
## single data frames
column.anno <- cdesc[samp.common, ]
row.anno <- rdesc
tab.expr.all <- Reduce(f=rbind, gct.expr)
rownames(tab.expr.all) <- make.unique(rownames(tab.expr.all)) 
rownames(row.anno) <- rownames(tab.expr.all)

## reorder columns
ord.idx <- order(column.anno[, ord.column])
column.anno <- column.anno[ord.idx,]
tab.expr.all <- tab.expr.all[, ord.idx]

## reorder rows
row.idx <- with(row.anno, order(geneSymbol, DataType))
tab.expr.all <- tab.expr.all[row.idx, ]
row.anno <- row.anno[row.idx, ]

## export
save(column.anno, row.anno, tab.expr.all, file = glue('data/data-brca2019-{label}.Rdata'))

