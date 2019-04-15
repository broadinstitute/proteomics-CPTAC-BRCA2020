#################################################################
## Filename: server.R
## Created: April 3, 2019
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC2.0 prospective Breast
##          Cancer cohort
## This defines the server logic of the Shiny-app.
#################################################################
library(shiny)


########################################################
## Define server logic
########################################################
shinyServer( function(input, output, session) {

    global <- reactiveValues(
      auth=T,
      all.genes=unique(row.anno[, GENE.COLUMN])
    )

    ##############################
    ## text field for passphrase
    output$auth.user <- renderUI({
      
      if(global$auth) return()
      list(
        passwordInput('passphrase',label='Enter passphrase', width=120),
        actionButton('authbutton', 'GO')
      )
    })
    
    ##############################
    ## check passphrase
    observeEvent(input$authbutton, {
      global$auth <- authenticateUser(input$passphrase)
    })
    
    ##############################
    ## ui
    ##############################
    output$ui.input <- renderUI({
      
      if(!global$auth) return()
      #if(is.null(input$genes)) return()
        
      list(
        ## text input
        #textInput('genes', label=paste('Paste a list of gene names (max. ', GENEMAX,')', sep=''), value=GENESSTART),
        selectizeInput('genes', label=paste('Enter your genes of interest (max. ', GENEMAX,')', sep=''), 
                       choices=global$all.genes,selected=GENESSTART, multiple=T),
        
        HTML('<br><br>'),

        fluidRow(
          ##column(2, submitButton('GO')),
          column(3, radioButtons('zscore', label='Z-score', choices=c('row', 'none'), selected='none')),
          column(3, radioButtons('allsites', label='pSTY sites', choices=c('most variable', 'all'), selected='most variable')),
          column(3, textInput('min.val', label='min', value=-3, width='80%')),
          column(3, textInput('max.val', label='max', value=3, width='80%'))
        ),
        fluidRow(
            column(6, selectizeInput('sort.after', 'Sort by', 
                                     choices=c('PAM50', 'NMF', 'HER2', 'HER2pam50.HER2amp','ERBB2.CNA', 'HER2pos.HER2amp', 'TP53', 'PIK3CA'), 
                                     selected='PAM50', multiple=TRUE)),
            column(6, radioButtons('sort.dir', '', choices=c('ascending', 'descending'), selected='ascending'))
            
        ),
          
        HTML('<br><br>'),
          
        ## download buttons
        # fluidRow(
        #     column(6, downloadButton('downloadHM', 'Download pdf')),
        #     column(6, downloadButton('downloadTab', 'Download Excel'))
        # ),
        #   
        HTML('<br><br>'),
        HTML('<p><b>Getting started</b></p>'),
        helpText('Enter your gene names of interest (official gene symbols, e.g. PIK3R2) into the text field. You can enter up to 20 genes.')#,
        #HTML('<p>For more details please see our publication <a href="http://cancerres.aacrjournals.org/content/78/10/2732" target="_blank_">Mundt <i>et al.</i> Cancer Research. 2018</a></p>')
        
        )
    })
    
    ##############################
    ## update list of input genes
    observeEvent(input$genes, {
      
      if(is.null(input$genes)) return()
      if(!global$auth) return()
      
      #cat('TEST2', input$genes, '\n')
      
      global$genes.input <- extractGenes(input$genes)
    })
    
    
    
    ##############################
    ## generate the heatmap
    output$plot <- renderPlot({
      
      if(!global$auth) return()
      if(is.null(input$genes)) return()
      
      #cat('TTT', input$genes)
      
      genes.vec <- extractGenes( input$genes )
      
      if(length(genes.vec)==0) return()
      
      hm=makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, zscore=input$zscore, 
                show.sites=input$allsites, min.val=input$min.val, max.val=input$max.val,
                anno.class=input$sort.after, sort.dir=input$sort.dir)
      global$expr.select <- hm
    },
    width = function(){ width=1400},
    height= function(){ height=ifelse( global$auth, 
                                       dynamicHeightHM(length( findGenesInDataset(extractGenes( input$genes ), input$allsites) ), 
                                                       length(unique(extractGenes( input$genes ))) ), 0 )}
    )
    
    #############################
    ## download heatmap
    output$downloadHM <- downloadHandler(

        ##filename = paste( FILENAMESTRING, ifelse(input$zscore, 'Zscore', ''),'.pdf', sep=''),
        filename = paste( FILENAMESTRING,'.pdf', sep=''),
        content = function(file){
            genes.vec <- extractGenes( input$genes )
            if(length(genes.vec)==0) return()
            hm=makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, filename=file, main=TITLESTRING, height=ifelse(length(genes.vec) < 4, 6.5, NA), as.logical(input$zscore))
            }
    )
    #############################
    ## download Excel
    output$downloadTab <- downloadHandler(
        filename = function(){paste(FILENAMESTRING, '.xlsx', sep='')},
        content = function(file){
            tab=as.data.frame(global$expr.select)
            WriteXLS('tab', ExcelFileName=file, SheetNames=FILENAMESTRING, FreezeCol=6, FreezeRow=5, row.names=T)
            }
    )
})


