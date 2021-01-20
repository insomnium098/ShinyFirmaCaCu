## Only run examples in interactive R sessions
library(shiny)
#library(sas7bdat)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(matrixStats)
library(ComplexHeatmap)
options(repos = BiocManager::repositories())



ui <-fluidPage(htmlOutput("text"),navlistPanel(
  tabPanel("Inicio",textOutput("welcome_message"),
           textOutput("welcome_message_2"),
           helpText("Para poder usar este software es necesario usar y
                    llenar la siguiente plantilla con los datos de expresion."),
           fluidRow(downloadButton('download_plantilla', 'Descargar Plantilla')),
           img(src="INCAN.jpg", height = 400, width = 400),
           img(src="Logo_Lab_Genomica.jpg", height = 300, width = 600)),
  
  tabPanel("1: Ingresar datos de expresion", fileInput("file1", "Elige la plantilla con los datos de expresion",
                                    accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
           inputPanel(tableOutput("contents")),
           helpText("La plantilla debera estar guardada en el formato original .csv")),

  
  #tabPanel("Plots"),
  
  tabPanel("2: Resultados",fluidRow(downloadButton('download_results', 'Descargar Tabla de Resultados')),
           dataTableOutput("value"),
           ),
           
  
  tabPanel("3: Heatmap",plotOutput("heatmap",
                                width = 1000,
                                height = 1000))
  
  
  
  
  
  
))

server <-function(input, output) {
  
  # You can access the value of the widget with input$file, e.g.
  
  output$welcome_message <- renderText("Firma Molecular CaCu V1.0")
  output$welcome_message_2 <- renderText("")
  
  output$contents <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    if(file.exists("a.csv")){
      file.remove("a.csv")
    }
    if(file.exists("pronostico_final.csv")){
      file.remove("pronostico_final.csv")
    }
    if(file.exists("DF_EXPRESSION_Y_CLASIFICADOR.csv")){
      file.remove("DF_EXPRESSION_Y_CLASIFICADOR.csv")
    }
    if(file.exists("plot.pdf")){
      file.remove("plot.pdf")
    }
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    a <- read.csv(inFile$datapath, sep = ",")
    write.csv(a,"a.csv")
    read.csv("a.csv")
    
  })
  
  output$value <- renderDataTable({
    #iris
    #inFile <- input$file1
    
    #if (is.null(inFile))
     # return(NULL)
    
    #a <- read.csv("a.csv")#read.csv(inFile$datapath, sep = ",")
    
    #b <- a[,1:3]
    datos_clasificador <- read.csv("CLASIFICADOR_FINAL.csv",
                                   row.names = 1, header = TRUE)
    datos_exp <- read.csv("a.csv",
                          row.names = 2, header = TRUE)
    datos_exp <- datos_exp[,2:length(colnames(datos_exp))]
    
    expresion_merge <- merge(datos_clasificador, datos_exp,
                             by.x="gene",
                             by.y = "row.names")
    pacientes <- colnames(datos_exp)
    
    for(paciente in pacientes){
      
      filt_1 <- expresion_merge %>% select(gene,clas_nr,
                                           clas_cr,paciente)
      
      cor_mala_respuesta <- rcorr(filt_1[,4],filt_1[,2], type="pearson")
      r_mala_respuesta <- cor_mala_respuesta$r[1,2]
      
      cor_buena_respuesta <- rcorr(filt_1[,4],filt_1[,3], type="pearson")
      r_buena_respuesta <- cor_buena_respuesta$r[1,2]
      
      score <- r_mala_respuesta - r_buena_respuesta
      score <- score + 0.1
      
      if(score > 0){
        pronostico <- "mal_pronostico"
      } else if (score < 0){
        pronostico <- "buen_pronostico"
      } else{
        pronostico <- "INTERMEDIO"
      }
      
      df_fin <- as.data.frame(paciente)
      df_fin$pronostico <- pronostico
      df_fin$score <- score
      
      #Si el datasset no existe, se crea
      if (!exists("df_pronostico_final")){
        df_pronostico_final <- df_fin
      }
      #Si el dataset existe, se une 
      if (exists("df_pronostico_final")){
        temp <-df_fin
        df_pronostico_final<- rbind(df_pronostico_final,df_fin)
        rm(temp)
      }
      
    }
    
    df_pronostico_final <- df_pronostico_final[2:length(rownames(df_pronostico_final)),]
    write.csv(df_pronostico_final,"pronostico_final.csv")
    
    
    ###CALCULO ZSCORE
    
    df_zscore <- data.frame(datos_exp, stringsAsFactors=FALSE) 
    names_datos <- rownames(df_zscore)
    
    #df_zscore <- as.data.frame(sapply(df_zscore, as.character))
    df_zscore <- as.data.frame(sapply(df_zscore, as.numeric))
    
    
    
    #Calculo de zscore
    #Convertir en matrix
    matriz_selected <- as.matrix((df_zscore))
    #restar matriz menos el promedio
    matriz_menos_ROWMEANS <- matriz_selected - rowMeans(matriz_selected)
    #CONVERTIR EN AS.MATRIX
    matriz_menos_ROWMEANS <- as.matrix(matriz_menos_ROWMEANS)
    #dividir matriz menos la desviacion estandar
    matriz_desvest <- matriz_menos_ROWMEANS / rowSds(matriz_selected)
    #renombrar la matriz ya con la resta y division de matriz
    matriz_zscore <- matriz_desvest
    
    rownames(matriz_zscore) <- names_datos
    
    #####
    
    
    ###UNIR EXPRESION Y CLASIFICADOR
    
    matriz_zscore <- as.data.frame(t(matriz_zscore))
    rownames(df_pronostico_final) <- df_pronostico_final$paciente
    
    df_cool_final <- cbind(df_pronostico_final, matriz_zscore)
    write.csv(df_cool_final,"DF_EXPRESSION_Y_CLASIFICADOR.csv")
    colnames(df_cool_final)[1:3] <- c("Paciente","Pronostico","Score")
    df_cool_final[,1:3]
    
  })
  
  output$selection <- reactive({return(input$file1[,c(1,2)])})
  
  #output$selected <- renderTable({selection()})
  
  test <- reactive({input$file1()})
  
  output$heatmap <- renderPlot({     # Refers  to putputs with output$<id>
  #  ggplot2::ggplot(data = iris, aes(y = Sepal.Length)) + geom_boxplot() # Refers to inputs with input$<id>
    datos <- read.csv("DF_EXPRESSION_Y_CLASIFICADOR.csv",
                      header = TRUE,
                      row.names = 1,
                      stringsAsFactors = FALSE)
    datos[is.na(datos)] <- 0
    
    
    ###cambiar x de pacientes a P
    
    pacientes <- datos$paciente#gsub("X","P", datos$paciente, fixed = TRUE)
    datos$paciente <- pacientes
    rownames(datos) <- pacientes
    
    #datos_heat <- as.data.frame(t(datos))
    datos_exp <- as.matrix(datos[,4:length(colnames(datos))])
    
    
    Heatmap(datos_exp)
    
    
    
    ha <- HeatmapAnnotation(foo = anno_empty(border = TRUE, width = unit(60, "native")),
                            which = "row",
                            gap= unit(2,"points"))
    
    
    ###Se ordenara el heatmap de acuerdo al score
    orden_rows <- order(datos$score)
    
    
    heatmap_1 <- Heatmap(datos_exp, cluster_columns = TRUE,
                         cluster_rows = FALSE,#row_names_gp = gpar(fontsize = 6),
                         show_row_names = TRUE,#top_annotation = top,col = colores,#colorRamp2(c(0, 15, 30), c("darkblue", "white", "firebrick1")),
                         name = "Expression",
                         row_order = orden_rows,
                         column_names_rot = 45 ,left_annotation = ha,
                         row_names_gp = gpar(fontsize = 10),
                         width = unit(10, "inches"), height = unit(10, "inches"))#,
    #row_names_gp = gpar(fontsize = 8,col = ifelse(vals> 0, "red", ifelse(vals== 0, "black", "blue"))))
    
    
    
    
    ###vals index son los vals ya ordenados
    ###Obtener order de los rows despues de hacer el clustering
    r_order <- row_order(heatmap_1)
    vals_index <- datos$score[r_order]
    
    ##vals son el score
    
    vals <- datos$score
    
    ###respuesta ordenada de acuerdo al score
    
    respuesta <- datos$pronostico[r_order]
    
    
    ####Plot para descargar
    ht = draw(heatmap_1)
    
    
    decorate_annotation("foo", {
      # value on x-axis is always 1:ncol(mat)
      x = 1:length(vals)#1:10
      # while values on y-axis is the value after column reordering
      value = vals_index#value[co]
      value_respuesta = respuesta
      pushViewport(viewport(yscale = c(0.4, (length(vals))+0.5), xscale = c(-0.7, 0.7)))
      grid.lines(c(0, 0), c(0.4, (length(vals))+0.5), gp = gpar(lty = 1, col = "red"),
                 default.units = "native")
      grid.points(value, rev(x), pch = 21, size = unit(2, "mm"),#16
                  gp = gpar(col = ifelse(value_respuesta == "mal_pronostico", "black", ifelse(value == "NR", "black", "black")),
                            fill = ifelse(value_respuesta == "buen_pronostico", "white", ifelse(value == "NR", "black", "black"))
                  ), default.units = "native")
      grid.xaxis(at = c(-0.5, 0, 0.5), gp = gpar(fontsize = 5), main=TRUE)
      grid.text("Score", unit(0.5, "npc"),unit(1.02, "npc") ,just = "center",
                gp = gpar(fontsize = 10))
      popViewport()
      lgd = Legend(labels = c("Bueno","Malo"), type = "points", pch = c(21,21),title = "Pronostico",
                   legend_gp = gpar(fill = c("white","black")))
      draw(lgd, x = unit(14.5, "npc"), y = unit(0.3, "npc"), just = c("right", "bottom"))
      popViewport()
      #####
      #Plot para graficar

      
      
      
    })
  })
  
  
  ############DESCARGAS
  output$download_results <- downloadHandler(
    filename = function(){
      ###Este es el nombre con que se guarda
      paste("pronostico_final","csv",sep=".")
    },
    content = function(con){
      
      ##Este es el nombre del archivo en el servidor
      file.copy("pronostico_final.csv", con)
    })
  
  output$download_plantilla <- downloadHandler(
    filename = function(){
      ###Este es el nombre con que se guarda
      paste("plantilla","csv",sep=".")
    },
    content = function(con){
      
      ##Este es el nombre del archivo en el servidor
      file.copy("PLANTILLA.csv", con)
    })
  
  
  
  
}

shinyApp(ui, server)
#}