###########################################################################
## Filename: server.R
## Created: Oct 24, 2019
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC prospective BRCA cohort
## This defines the server logic of the Shiny-app.
##########################################################################
library(shiny)


########################################################
## Define server logic
########################################################
shinyServer( function(input, output, session) {

    global <- reactiveValues(
      zscore='row',
      sort.dir='ascending',
      auth=T,
      init=F,
      all.genes=unique(row.anno[, GENE.COLUMN])
      
    )

    ##############################
    ## text field for passphrase
    output$auth.user <- renderUI({
      
      if(global$auth) return()
      list(
        passwordInput('passphrase',label='Enter password', width=120, placeholder = 'Password'),
        actionButton('authbutton', 'GO', )
      )
    })
    
    ##############################
    ## check passphrase
    observeEvent(input$authbutton, {
      global$auth <- authenticateUser(input$passphrase)
      global$init <- F
    })

    ##############################
    ## ui
    ##############################
    output$ui.input <- renderUI({
      
      if(global$init) return()
      validate(need(expr=(global$auth), message=HTML('wrong password!'), label = 'authenticate'))
        
      list(
        ## text input
        selectizeInput('genes', label=paste('Enter your genes of interest (max. ', GENEMAX,')', sep=''), 
                       choices=global$all.genes, selected=GENESSTART, multiple=T),
        
        HTML('<br><br>'),

        fluidRow(
          column(3, radioButtons('zscore', label='Z-score', choices=c('row', 'none'), selected='row')),
         # column(3),
          column(3, radioButtons('allsites', label='phospho/acetyl sites', choices=c('most variable', 'all'), selected='most variable')),
          column(3, textInput('min.val', label='min', value=-2, width='80%')),
          column(3, textInput('max.val', label='max', value=2, width='80%'))
        ),
        fluidRow(
            column(12, selectizeInput('sort.after', 'Sort by', 
                                     choices=names(columns.to.sort), 
                                     selected=names(columns.to.sort)[1], multiple=FALSE))#,
            #column(6, radioButtons('sort.dir', '', choices=c('ascending', 'descending'), selected='ascending'))
            
        ),
          
        HTML('<br><br>'),
          
        ## download buttons
        fluidRow(
             column(6, downloadButton('downloadHM', 'Download PDF')),
             column(6, downloadButton('downloadTab', 'Download Excel'))
        ),
           
        HTML('<br><br>'),
        HTML('<p><b>Getting started</b></p>'),
        helpText(glue("Enter your gene names of interest (official gene symbols, e.g. {GENESSTART[sample(1:length(GENESSTART), 1)]}) into the text field. You can enter up to 20 genes.")),
        HTML('<p><b>Dataset</b></p>'),
        HTML('<p>Copy number aberrations are relative to matching normal blood sample and are on log2(CNA)-1 scale. For other data types the heatmap depicts abundances relative to common tumor reference.</p>'),
        HTML(glue( "<table border-spacing:5px><tr><th>Type</th><th># features</th><th># samples</th></tr>\n
                   <tr><td>CNA</td><td>{sum(grepl('_CNA', row.anno$DataType))}</td><td>{ncol(tab.expr.all)}</td></tr>\n
                   <tr><td>mRNA</td><td>{sum(grepl('_RNAseq', row.anno$DataType))}</td><td>{ncol(tab.expr.all)}</td></tr>\n
                   <tr><td>Protein</td><td>{sum(grepl('_Protein', row.anno$DataType))}</td><td>{ncol(tab.expr.all)}</td></tr>\n
                   <tr><td>phosphosites </td><td>{sum(grepl('_pSTY', row.anno$DataType))}</td><td>{ncol(tab.expr.all)}</td></tr>\n
                   <tr><td>acetylsites</td><td>{sum(grepl('_acK', row.anno$DataType))}</td><td>{ncol(tab.expr.all)}</td></tr>\n
                   </table>")),
        
        HTML('<br><p>For more details please see our publication <a href="https://www.cell.com/cell/fulltext/S0092-8674(20)31400-8" target="_blank_">Krug <i>et al.</i> Cell, 2020</a></p>')
        
        )
    })
    
    
    
    ##############################
    ## update list of input genes
    observeEvent(input$genes, {
      
      if(is.null(input$genes)) return()
      if(!global$auth) return()
      
      global$genes.input <- extractGenes(input$genes)
    })
    
    
    
    ##############################
    ## generate the heatmap
    output$plot <- renderPlot({
      
      validate(need(global$auth, message = '' ,label = 'auth_plot'))
      validate(need(!is.null(input$genes), message = 'Nothing to show' ,label = 'no-genes-entered'))
      
      
      if(!global$auth) return()
      if(is.null(input$genes)) return()
      
      genes.vec <- extractGenes( input$genes )
      
      if(length(genes.vec)==0) return()
      
     hm=makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, 
               row.anno=row.anno, zscore=input$zscore, show.sites=input$allsites, 
               min.val=as.numeric(input$min.val), 
               max.val=as.numeric(input$max.val), 
               anno.class=columns.to.sort[input$sort.after], 
               sort.dir=global$sort.dir)
     
     global$expr.select <- hm
    },
    width = function(){ width=1500},
    height= function(){ height=ifelse( global$auth, 
                                       dynamicHeightHM(length( findGenesInDataset(extractGenes( input$genes ), input$allsites) ), 
                                                       length(unique(extractGenes( input$genes ))) ), 100 )}
    )
    
    #############################
    ## download heatmap
    output$downloadHM <- downloadHandler(
        filename =  function(){paste(FILENAMESTRING,  '-',  gsub(' |\\:','-', Sys.time()), '.pdf', sep='')},
        content = function(file){
          genes.vec <- extractGenes( input$genes )
          if(length(genes.vec)==0) return()

          pdf(file, width=1500/72, height=dynamicHeightHM(length( findGenesInDataset(extractGenes( input$genes ), input$allsites) ), 
                                                          length(genes.vec))/72 )
          hm=try(makeHM(genes.vec, expr=tab.expr.all, column.anno=column.anno, row.anno=row.anno, 
                        zscore=input$zscore, show.sites=input$allsites, 
                        min.val=as.numeric(input$min.val), max.val=as.numeric(input$max.val), 
                        anno.class=columns.to.sort[input$sort.after], sort.dir=global$sort.dir, filename = file, main = TITLESTRING.HEATMAP))
          dev.off()
          }
     )      
    
    #############################
    ## download Excel
    output$downloadTab <- downloadHandler(
        filename = function(){paste( FILENAMESTRING, '-',  gsub(' |\\:','-', Sys.time()) ,'.xlsx', sep='')},
        content = function(file){
            tab=as.data.frame(global$expr.select)
            WriteXLS('tab', ExcelFileName=file, SheetNames=FILENAMESTRING, FreezeCol=4, FreezeRow=10, row.names=T)
            }
    )
})


