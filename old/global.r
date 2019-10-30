#################################################################
## Filename: global.r
## Created: April 3, 2019
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC2.0 prospective Breast
##          Cancer cohort
## This file imports the underlying data and contains global functions
## and variables used by 'ui.R' and 'server.R'
#################################################################

#source('pheatmap.r')
library(pacman)
p_load(scales)
p_load(gtable)
p_load(ComplexHeatmap)
p_load(RColorBrewer)
p_load(circlize)

source('pheatmap.r')


## import the data
load('data/data_bc2019.RData')

## global parameters
GENE.COLUMN <<- 'geneSymbol' 
DATATYPE.COLUMN <<- 'DataType'

GENESSTART <<- c('TP53', 'ERBB2', 'PIK3CA', 'GATA3', 'ESR1', 'PGR')

GENEMAX <<- 20
TITLESTRING <<- 'CPTAC prospective BRCA (v5.0 freeze)'
WINDOWTITLE <<- 'CPTAC-BC2019'
GAPSIZEROW <<- 20
FILENAMESTRING <<- 'CPTAC2.0_prospBRCA'

cellwidth <<- 8
cellheight <<- 10

p_load(RColorBrewer)
p_load(gplots)
p_load(WriteXLS)
p_load(grid)
p_load(bcrypt)

anno.all <- c('PAM50'='PAM50', 
              'NMF.cluster'='NMF.cluster', 
              'ER'='ER', 
              'PR'='PR',
              'ERBB2.Proteogenomic.Status'='ERBB2.Proteogenomic.Status', 
              'HER2.Amplified'='HER2.Amplified', 
              'PAM50.Her2.HER2.status'='PAM50.Her2.HER2.status', 
              'TP53.mutation.status'='TP53.mutation.status', 
              'PIK3CA.mutation.status'='PIK3CA.mutation.status')

##################################################################
## 21060613 bcrypt
authenticateUser <- function(passphrase){
  if(nchar(as.character(passphrase)) > 0){
    return(checkpw(as.character(passphrase), "$2a$12$XXsIDMNZMaOgRcaq8JX8muN3gsG93YGltGcqfiqkGGgDkJKV4cwri"))
  } else {
    return(FALSE)
  }
}



###################################################################
## heatmap using ComplexHeatmap package
MyComplexHeatmap <- function(m, rdesc, cdesc, cdesc.color, max.val){
  
  if(is.null(max.val)){
    m.max <- ceiling(max(abs(m), na.rm=T))
  } else {
    m.max <- max.val
    m[m > m.max] <- m.max
    m[m < -m.max] <- -m.max
  }
  ## #####################################
  ## complexheatmap
  col.hm <- colorRamp2(seq(-m.max, m.max, length.out=11), rev(brewer.pal (11, "RdBu")))
  cdesc.ha <- HeatmapAnnotation(df=cdesc, col=cdesc.color,
                                
                                show_legend = T, show_annotation_name = T, 
                                annotation_name_side = 'left',
                               
                                 annotation_legend_param=list(
                                  direction='horizontal',
                                  vt_gap = unit(1, 'cm')
                                  
                                  #title_position = "leftcenter"
                                )
  )
  hm <- Heatmap(m, col=col.hm,
                cluster_columns = F,
                cluster_rows = F,
                top_annotation = cdesc.ha,
                split = rdesc$geneSymbol,
                #row_dend_side = 'right',
                name='relative abundance',
                show_row_names = F,
                show_column_names = F,
                use_raster = FALSE,
                  )
  ## plot
  draw(hm, annotation_legend_side='bottom')
}


##################################################################
## function to extract gene names from a string of
## character
extractGenes <- function(genes.char){

    gene.max=GENEMAX

    #cat('TEST:',genes.char, '\n')
    if(is.null(genes.char))
      return(NULL)
    #if( nchar(genes.char) == 0 ){
    if( length(genes.char) == 0 ){
      
        return(NULL)
    }
    ## extract genes
    genes.vec= unlist(strsplit(genes.char, ','))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ' '))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ';'))

    ## unique gene names
    genes.vec <- unique(genes.vec)

    ## limit to 'gene.max' genes
    if(length(genes.vec) > gene.max){
        warning(paste('more than', gene.max,'gene ids submitted! Showing results for the first 20 genes in the list.\n'))
        genes.vec <- genes.vec[1:gene.max]
    }
    return(genes.vec)
}

##################################################################
## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n.entries, n.genes){
  
  height = (n.entries+2)*11 + (n.genes-1)*GAPSIZEROW + 140
  
  return(height)
}


#######################################################
## find all entries with associated gene name in
## the dataset. returns vector of indices.
findGenesInDataset <- function(gene, show.sites){
  
  ## remove spaces
  gene <- gsub(' ', '', gene )
  gene <- unique(gene)
  ## remove emtpy strings
  gene.nchar=which(nchar(gene) == 0)
  if(length(gene.nchar) > 0)
    gene <- gene[-gene.nchar]
  
  if(length(gene) == 0) return()
  
  ## check whether the genes are present in the dataset
  gene.idx <- grep( paste(paste('(^|,)', gene, '($|,)', sep=''), collapse='|'), gsub(' ', '', row.anno[, GENE.COLUMN]) )
  if( length(gene.idx) == 0 ){
    stop('None of the gene ids you have entered could be found in the dataset!\n')
  }

  ## use row names
  gene.idx <- rownames(tab.expr.all)[gene.idx]
  
  ## exract data and remove empty rows
  data.tmp <- tab.expr.all[gene.idx, ]
  row.anno.tmp <- row.anno[gene.idx, ]
  
  ## most variable site
  if(show.sites=='most variable'){
    
    ## extract MS phospho
    vm.idx <- grep('pSTY|acK', row.anno.tmp[, DATATYPE.COLUMN])
    if( length(vm.idx) > 0 ){
      
      vm.sd <- apply(data.tmp[ vm.idx, ], 1, sd, na.rm=T)
      
      rm.idx <- tapply(vm.sd, row.anno.tmp[vm.idx, GENE.COLUMN],  function(x) names(x)[which.max(x)])
      rm.idx <- setdiff( names(vm.sd), unlist(rm.idx) )
      
      gene.idx <- setdiff(gene.idx, rm.idx)
      
      data.tmp <- data.tmp[gene.idx, ]
      row.anno.tmp <- row.anno.tmp[gene.idx, ]
      
    }
  }

  return(gene.idx)
}


#################################################################
## draw the actual heatmap
##
#################################################################
makeHM <- function(gene, filename=NA, expr=tab.expr.all, 
                   column.anno=column.anno, row.anno=row.anno, zscore="none", 
                   anno.class='PAM50', sort.dir, 
                   show.sites='all', min.val=-3, max.val=3, ...){

    n.bins=12


    ## ################################################
    ## reorder
    #ord.idx <- order(column.anno[, anno.class], decreasing = ifelse(sort.dir == 'descending', T, F))
    ord.idx <- with(column.anno, order( eval(parse( text=anno.class ))))
    column.anno <- column.anno[ord.idx,]
    expr <- expr[, ord.idx]
    
   
    gene.idx <- findGenesInDataset(gene, show.sites) 
    
    
    #################################
    ## extract sample ids
    sampleIDs <- colnames(expr)

    #################################
    ## extract genes of interest
    expr.select <- expr[gene.idx, ]
    row.anno.select <- row.anno[gene.idx, ]
    
    ## order
    
    #####################################################
    ## labels for the rows in the heatmap
    featureIDs.anno.select <- paste(row.anno.select[ , 'geneSymbol'],
                                    gsub( '5_acK', 'acK',
                                    gsub( '4_pSTY', 'pSTY', 
                                          gsub('3_Protein', 'Protein', 
                                               gsub('2_RNAseq', 'RNA-Seq', 
                                                    gsub('1_CNA', 'CNA', row.anno.select[ , 'DataType'])
                                                    ) 
                                               )
                                          )
                                    )
                                    )

    #####################################################
    ## add phosphosite annotation
    #####################################################
    ## RPPA
    #rppa.psty.idx <- grep('RPPA pSTY', featureIDs.anno.select)
    #if(length(rppa.psty.idx)>0){
    #    featureIDs.anno.select[ rppa.psty.idx ] <- paste( sub(' pSTY', '', featureIDs.anno.select[ rppa.psty.idx ]), sub('.*_(.*?)\\..*', '\\1', row.anno.select[ rppa.psty.idx, 'ID2']))
    #}
    ## MS
    ms.psty.idx <- grep('pSTY', featureIDs.anno.select)
    if(length(ms.psty.idx)>0){
        featureIDs.anno.select[ ms.psty.idx ] <- paste( sub('pSTY', '', 
                                                            featureIDs.anno.select[ ms.psty.idx ]), 
                                                        paste('p', sub('.*_([S|T|Y][0-9]*)[s|t|y].*', '\\1', row.anno.select[ ms.psty.idx, 'ID']), sep=''), sep='' )
    }
    ms.ack.idx <- grep('acK', featureIDs.anno.select)
    if(length(ms.ack.idx)>0){
      featureIDs.anno.select[ ms.ack.idx ] <- paste( sub('acK', '', 
                                                         featureIDs.anno.select[ ms.ack.idx ]), 
                                                     paste('ac', sub('.*_([K][0-9]*)[k].*', '\\1', row.anno.select[ ms.ack.idx, 'ID']), sep=''), sep='' )
    }
    
    
    #################################
    ## apply zscore
    rownames(expr.select) <- featureIDs.anno.select

    if(zscore == "row"){
      
        ## exclude CNA data from Z-scoreing
        expr.select.zscore.tmp <- lapply( rownames(expr.select), function(xx){x=expr.select[xx,];
            if( length(grep( 'CNA', xx)) == 0)return((x-mean(x, na.rm=T))/sd(x, na.rm=T));
            if( length(grep( 'CNA', xx)) > 0)return(x);
        })
        expr.select.zscore <- matrix(unlist(expr.select.zscore.tmp), ncol=ncol(expr.select), byrow=T, dimnames=list(rownames(expr.select), colnames(expr.select)))
    } else {
        expr.select.zscore <- expr.select
    }

    ## cap at -3/3
    expr.select.zscore[which(expr.select.zscore < min.val)] <- min.val
    expr.select.zscore[which(expr.select.zscore > max.val)] <- max.val

    ##############################
    ## column annotation
    #column.anno.fig <- column.anno[, c('PAM50', 'ER', 'PR', 'HER2')]
    column.anno.fig <- column.anno[, anno.all]
    colnames(column.anno.fig) <- names(anno.all)

    ##############################
    ## colors for column annotation
    column.anno.col <- list(
        PAM50=c(Basal='red', Her2='violet', LumA='blue', LumB='cyan', 'Normal-like'='grey'),
        ER=c(positive='black', negative='white', unknown='grey'),
        PR=c(positive='black', negative='white', unknown='grey'),
        ERBB2.Proteogenomic.Status=c(positive='black', negative='white', unknown='grey'),
        NMF.cluster=c(C1='yellow', C2='red', C3='blue', C4='violet'),
        TP53.mutation.status=c('1'='darkblue', '0'='white'),
        PIK3CA.mutation.status=c('1'='darkblue', '0'='white'),
        HER2.Amplified=c('0'='white', '1'='darkgreen'),
        PAM50.Her2.HER2.status=c(positive='darkgreen', negative='lightgreen')
    )

    ################################
    ## gaps
    ################################
    if(sort.dir == 'descending')
      gaps.column=cumsum(c( rev(table(column.anno[, anno.class])) ))
    if(sort.dir == 'ascending')
      gaps.column=cumsum(c(  table(column.anno[, anno.class]) ))
    gaps.row=cumsum(table(sub(' .*', '', featureIDs.anno.select)))


    ################################
    ## colors misc
    color.breaks = seq(min.val, max.val, length.out=n.bins)
    color.hm =  colorRampPalette( c('blue', 'grey', 'red'))(length(color.breaks))
    color.border = 'white'

    legend_breaks=seq(-3, 3, 1)
    legend_labels=c('-3               ', '-2', '-1', ' 0', '+1', '+2' ,'+3')


    ###############################
    ## heatmap
    #View(row.anno.select)
    #MyComplexHeatmap <- function(m, cdesc, cdesc.color, class.variable, variable.other, max.val)
   # MyComplexHeatmap(expr.select.zscore, row.anno.select, column.anno.fig, column.anno.col, max.val=NULL )
    
    pheatmap(expr.select.zscore, cluster_row=F, cluster_col=F,  annotation_col=column.anno.fig, annotation_colors=column.anno.col,  
             scale = "none", labels_row=featureIDs.anno.select, border_color=color.border, gaps_col=gaps.column, 
             gaps_row=gaps.row, color=color.hm, cellwidth=cellwidth, cellheight=cellheight, 
             labels_col=sampleIDs, breaks=color.breaks, legend_breaks=legend_breaks, legend_labels=legend_labels, na_col='white', filename=filename, ...)

  
    #########################################################################################
    ##
    ## - return part of table that is shown in the heatmap and that can be downloaded
    ## - change the CNA values (-3, -1, 0, 1, 3) back to the orignial values (-1, -.3, 0, .3, 1)
    ##
    #########################################################################################
    cna.idx <- grep('CNA', rownames(expr.select))
    if(length(cna.idx) > 0){
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == -3 ] <- -1
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == -1 ] <- -.3
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == 1 ] <- .3
        expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == 3 ] <- 1
    }
    ## add row annotation
    mat.row <- as.data.frame( cbind( row.anno[gene.idx, ], expr.select, deparse.level=0 ), stringsAsFactors=F )
    rownames(mat.row) <- make.names(featureIDs.anno.select, unique = T)

    ## column annotation
    mat.col <- as.data.frame( cbind(matrix('', nrow=ncol(column.anno.fig), ncol=ncol(row.anno)), t(column.anno.fig), deparse.level=0), stringsAsFactors=F)
    colnames(mat.col) <- colnames(mat.row)

    ## put everything together
    mat <- rbind( mat.col, mat.row)

    return(mat)
}