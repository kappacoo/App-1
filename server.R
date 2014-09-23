library(shiny)
source("helper.R")
stderr <- function(x) sqrt(var(x, na.rm=TRUE)/sum(!is.na(x)))
colSE<- function(x) apply(x, 2, stderr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
		
		# Show example table in Tab 2 ==========
		
		isolate({
		output$exampleTable<- renderPrint({
			print(exampleTable)
			})	
	
		# loading in the file ======================
		dataInput <- reactive({
			inFile <- input$file
		
			if (!is.null(inFile)){
				dataFile<-read.csv(inFile$datapath, stringsAsFactors=F)
				return(dataFile)
				} else 
			return(exampleTable)
		
			})
		# Metabolites to list ======================
				
		m_options <- reactive({
			number_meta <- which(meta[,1] %in% names(dataInput()))
			m_options1<- NULL
			for (i in number_meta) {
				m_options1 <- append(m_options1, meta[i,1])
				}
		
			setNames(m_options1, m_options1)
		
			})
		
		# Groups to List ==========================
		g_options <- reactive({
			g_options1<- unique(dataInput()[["Group"]])
			g_options1<-as.character(g_options1)
			setNames(g_options1, g_options1)
		})
		
		})
		
		
		# Output of metabolite list ====================
		
		output$text <- renderPrint({
			if(!is.null(g_options())){
				str(g_options())
   			g_options()[[1]]
					}
 			 })
		
		
		# Select metabolite UI ======================
		
		output$variables <- renderUI({
			if(!is.null(m_options())){
			
     		selectInput("Name", 'Metabolite',
     			m_options(), 
     			selectize=F)
     			
     			}
			})
			
		# Collected metabolite isotopes into one table =====
		
		metaboliteTable <- reactive({
			df<-dataInput()
			names(df)<- sub("[[:punct:]]([13c]{3}|[M]{1})[[:punct:]]?", "_M", names(df), ignore.case=T)
			if(!is.null(input$Name)){
		
				metabolite<-data.frame("M0"= df[[input$Name]])
				for (i in seq(1:9)){
					if(!is.null(df[[sprintf('%s_M%s',input$Name, i)]])){
						metabolite[[sprintf("M%s", i)]] <- df[[sprintf('%s_M%s', input$Name, i)]]
					}
				}
			return(metabolite)
			}
		})
		
		# Isoform UI ====================================
		
		output$isotope<-renderUI({
			  checkboxGroupInput("i_options", label = "Plot specific isotopes:", 
    choices = names(metaboliteTable()),
    selected = names(metaboliteTable()), inline=TRUE)

		})
		
		# Metabolite list with input from isoform UI =========
		
		metaboliteNew<- reactive({
			metabolite<- metaboliteTable()[,which(names(metaboliteTable()) %in% input$i_options)]
			
			return(metabolite)
		})
		
			
		# Select groups UI ============================
		
		output$variable2 <- renderUI({
			if(!is.null(g_options())){

			selectInput("selected2",
			 'Choose groups to plot; these will be ploted in the order thay are listed.', 
			 g_options(), selected=g_options(), multiple = TRUE)
			} else (return("Edit data table to include column with Group header"))
			})


		# Pool of all isotope ============================
		
		pool <- reactive({
			apply(metaboliteTable(), 1, function(x) sum(x, na.rm=TRUE))
		})
		
		# Pool Avg by group ===========================
		poolAvg1 <- reactive({
			if(!is.null(input[["Name"]])&!is.null(input[["selected2"]])){
			df<-dataInput()
			df$Group<-factor(df$Group, 
				levels = input$selected2)
			
			s<-split(pool(),df$Group)
			avg<-sapply(s, function(x) mean(x, na.rm=T))
			return(avg)
			}
		})
		
		
		poolAvg <- reactive({
			if(!is.null(input[["Name"]])&!is.null(input[["selected2"]])){
			df<-dataInput()
			df$Group<-factor(df$Group, 
				levels = input$selected2)
			
			s<-split(pool(),df$Group)
			if(input$type==1){
				return(poolAvg1())
				} else{
					perc<-NULL
				for(i in seq(1,length(input$selected2))){
					perc[[i]]<-s[[i]]/poolAvg1()[i]*100
				}
				
				avg<-sapply(perc, function(x) mean(x, na.rm=T))
				
				return(avg)
				}		
			}
		})
		
		# Pool SEM by group ===========================
		
		poolSEM <- reactive({
						if(!is.null(input[["Name"]])&!is.null(input[["selected2"]])){
			df<-dataInput()
			df$Group<-factor(df$Group, 
				levels = input$selected2)
			
			s<-split(pool(),df$Group)
			SEM<-sapply(s, function(x) stderr(x))
			if(input$type==1){
				return(SEM)
				} else{
					perc<-NULL
				for(i in seq(1,length(input$selected2))){
					perc[[i]]<-s[[i]]/poolAvg1()[i]*100
				}
				
				SEM<-sapply(perc, function(x) stderr(x))
				
				return(SEM)
				}		
			}
		})
		
		# Total percent plot ==========================
		
		plotingAvg <- reactive({
			if(input$merge_m == F | is.null(input$merge_m)){
			if(!is.null(input[["Name"]])&!is.null(input[["selected2"]])){
			df<-dataInput()
			df$Group<-factor(df$Group, 
				levels = input$selected2)
			
			s<-split(metaboliteNew(),df$Group)
			if (length(input$i_options)==1){
				avg<-sapply(s, function(x) mean(x, na.rm=T))
				} else{
				avg<-sapply(s, function(x) colMeans(x, na.rm=T))
				}
			if (input$norm == T){
				x<-input$selected2[1]
				normValue<-t(avg[,1])
				metaboliteTable1<-metaboliteNew()/normValue[rep(1,nrow(metaboliteTable()))]
				
				s<-split(metaboliteTable1,df$Group)
				if (length(input$i_options)==1){
					avg<-sapply(s, function(x) mean(x, na.rm=T))
					} else
					avg<-sapply(s, function(x) colMeans(x, na.rm=T))
			
				} else {
				if(input$type==2){
					perc<-NULL
					for(i in seq(1,length(input$selected2))){
						perc[[i]]<-s[[i]]/poolAvg1()[i]*100
						}
					if (length(input$i_options)==1){
						avg<-sapply(perc, function(x) mean(x, na.rm=T))
						
						} else {
						avg<-sapply(perc, function(x) colMeans(x, na.rm=T))
						colnames(avg)<-levels(df$Group)
						}
					}	
				}
			return(avg)
			}
			} else {return()}
			
		})
		
		# SEM =============
		
		plotingSEM<- reactive({
			if(input$merge_m == F | is.null(input$merge_m)){
			if(!is.null(input[["Name"]])&!is.null(input[["selected2"]])){
			df<-dataInput()
			df$Group<-factor(df$Group, 
				levels = input$selected2)
			
			s<-split(metaboliteNew(),df$Group)
			if (length(input$i_options)==1){
				avg<-sapply(perc, function(x) mean(x, na.rm=T))
				SEM<-sapply(s, function(x) stderr(x))
				} else {
			avg<-sapply(s, function(x) colMeans(x, na.rm=T))
			SEM<-sapply(s, function(x) colSE(x))
			}
			if (input$norm == T){
				x<-input$selected2[1]
				normValue<-t(avg[,1])
				metaboliteTable1<-metaboliteNew()/normValue[rep(1,nrow(metaboliteTable()))]
				
				s<-split(metaboliteTable1,df$Group)
				if (length(input$i_options)==1){
				avg<-sapply(perc, function(x) mean(x, na.rm=T))
				SEM<-sapply(s, function(x) stderr(x))
				} else {
				avg<-sapply(s, function(x) colMeans(x, na.rm=T))
				SEM<-sapply(s, function(x) colSE(x))
				}
			} else{
				if(input$type==2){
			perc<-NULL
				for(i in seq(1,length(input$selected2))){
					perc[[i]]<-s[[i]]/poolAvg1()[i]*100
				}
				if (length(input$i_options)==1){
				SEM<-sapply(perc, function(x) stderr(x))
				} else {
				SEM<-sapply(perc, function(x) colSE(x))
				}
				}	
			}
			return(SEM)
			}
			} else {return()}
			

			})
			
		# If merge is TRUE ===========================
		merge <- reactive({
			if(input$merge_m ==T){
				if(class(metaboliteNew())=="data.frame"){
				rowSums(metaboliteNew())
				} else{
					metaboliteNew()
				}
				
			}
		})
			
		
		mergeAvg <- reactive({
			if(input$merge_m ==T){
			df<-dataInput()
			df$Group<-factor(df$Group, 
				levels = input$selected2)
			
			s<-split(merge(),df$Group)
			avg<-sapply(s, function(x) mean(x, na.rm=T))
			setNames(avg,levels(df$Group))
			if(input$type==1){
				return(avg)
				} else{
					perc<-NULL
				for(i in seq(1,length(input$selected2))){
					perc[[i]]<-s[[i]]/poolAvg1()[i]*100
				}
				
				avg<-sapply(perc, function(x) mean(x, na.rm=T))
				
				return(avg)
				}		
			}
		})
		
		# Pool SEM by group ===========================
		
		mergeSEM <- reactive({
			if(input$merge_m==T){
						if(!is.null(input[["Name"]])&!is.null(input[["selected2"]])){
			df<-dataInput()
			df$Group<-factor(df$Group, 
				levels = input$selected2)
			
			s<-split(merge(),df$Group)
			SEM<-sapply(s, function(x) stderr(x))
			if(input$type==1){
				return(SEM)
				} else{
					perc<-NULL
				for(i in seq(1,length(input$selected2))){
					perc[[i]]<-s[[i]]/poolAvg1()[i]*100
				}
				
				SEM<-sapply(perc, function(x) stderr(x))
				
				return(SEM)
				}		
			}
			}
		})


	   
	  
	   # test ===========================================
	   #=================================================
	   # need to write verbatimTextOutput("test") in ui.R
	   output$test <- renderPrint({ 
									
				return(mergeAvg())
			
		})	   
		
		output$test2 <- renderPrint({
			
			return(plotingSEM())
		})	
	   #==================================================
	   # Plot ==============

			# Printing plot =========================
			
 			output$plot1 <- renderPlot({
 				print(plotInput2())
 			})
 			
 			# non-reactive function for printing plot ===========
 	plotInput2 <- function(){ 
			if(is.null(input$selected2)){return()}
			if(input$merge_m == F | is.null(input$merge_m)){
				avg<-plotingAvg()
				SEM<-plotingSEM()
			} else {
				avg<-mergeAvg()
				SEM<-mergeSEM()
			}			
			
			if(is.null(SEM)){
				SEM<-avg*0
			}
			
			
			ylabel <- if(input$type==1){
				"Relative Peak Area"
			} else {
				"Percent of total"
			}
			max <- if(input$type2==1){
				max(avg)+max(SEM)
				} else {
				max(poolAvg())+max(poolSEM())
				}
			
			if(input$type2==2){
				barplot(poolAvg(), ylim=c(0,max), ylab="", main="", xaxt="n", yaxt="n", col="white")
				par(new=T)
				}	
			plotbeside<- if(input$type2==1){
				TRUE} else {FALSE}
		
			points<-barplot(avg, names=input$selected2, beside=plotbeside, ylab=ylabel, main=input$Name,ylim=c(0,max))
			
			errlength <- 1/length(points)
			
			
			if(input$type2==2){
				for (i in seq(1,length(points))) {
					if(poolSEM()[i]>0){
		arrows(points[i], poolAvg()[i], points[i], poolAvg()[i]+poolSEM()[i], angle=90,length=errlength)
		}
					}
				if(input$merge_m==T){
					for (i in seq(1,length(points))) {
						if(SEM[i]>0){
		arrows(points[i], avg[i], points[i], avg[i]+SEM[i], angle=90,length=errlength)
		arrows(points[i], avg[i], points[i], avg[i]-SEM[i], angle=90,length=errlength)
						}
					}
				}
			} else {
			if(class(points)=="matrix"& ncol(points)>1){
			for(i in seq(1,ncol(points))){
				for (k in seq(1,nrow(points))){
					if(SEM[k,i]>0){
		arrows(points[k,i], avg[k,i], points[k,i], avg[k,i]+SEM[k,i], angle=90,length=errlength)
		arrows(points[k,i], avg[k,i], points[k,i], avg[k,i]-SEM[k,i], angle=90,length=errlength)
		}
				}
				}
				} else{
					for (i in seq(1,length(points))) {
						if(SEM[i]>0){
		arrows(points[i], avg[i], points[i], avg[i]+SEM[i], angle=90,length=errlength)
		arrows(points[i], avg[i], points[i], avg[i]-SEM[i], angle=90,length=errlength)
					}
					}
			}	
			}
			
			}
		
 	# Download Plot ===============================
 	
 	output$downloadData <- downloadHandler(
 
    filename = function() { paste(input$file,input$Name, '.pdf', sep='') },
    content = function(file) {
    	pdf(file)
    	plotInput2()
	dev.off()
    }
  )		
			

})
