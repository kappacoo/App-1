library(shiny)

shinyUI(navbarPage("Fluxomics Analysis",
	tabPanel("About",
		fluidRow(
			column(6,
				wellPanel(
					h4("Fluxomics Analysis Tool"),
					p("This tool is designed for interactive plotting of fluxomics data. Use it to visualize data from isotopic labeling experiments."),
					h4("How to use"),
					p("Loading in your data table. Make sure things are correctly labeled and formated to have metabolites values in columns and individual samples in rows."),
					p("Head to Barplots. Pick a metabolite from the dropdown menu. Specify the plotting order of your groups and whether you want to look at isotope percentages vs. total metabolite pool sizes."),
					p("Like the plot you made? Download it by clicking the button.")
					)
				)
			)
		),
	tabPanel("Load Data",
		fluidRow(
			column(6,
				fileInput("file", label=h3("File input")),
    			p("Enter a .csv file with metabolites in columns and individual samples in rows. Metabolite isotope distribution should be identified in one of the following manners:", tags$ul(
    			tags$li("Metabolite, Metabolite_13C1, Metabolite_13C2"),
    			tags$li("Metabolite, Metabolite M+1, Metabolite M+2"),
    			tags$li("Metabolite, Metabolite 13C-1, Metabolite 13C-2")),
    			 "**Include a column with group assignments (label as 'Group'). See below for an example table.")
    			)),
    	fluidRow(
    		column(8, h5("Example Table:"),
    		verbatimTextOutput("exampleTable")
    		))
 		),

	tabPanel("Barplots",
		fluidRow(
			column(3, 
				wellPanel(		
					uiOutput("variables"),
					br(),
					uiOutput("variable2"),
					checkboxInput("norm", 
					label="Normalize data to first group"),
					radioButtons("type", "Plot as", choices= list("Relative Peak area" = 1, "Percent of total"=2)),
					radioButtons("type2","", choices= list("Side-by-side" = 1, "Stacked" = 2), inline=T)
					
					)
				),
			column(6,
						
			uiOutput("isotope"),
			checkboxInput("merge_m", "Merge values"),
			plotOutput("plot1"),
			downloadButton('downloadData', 'Download Plot')

			)	
			)
	
		)
	)
)