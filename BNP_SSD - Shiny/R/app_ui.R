app_ui <- function() {
  tagList(
    tags$head(tags$title("BNP_SSD")),
    shinyjs::useShinyjs(),
    waiter::useWaiter(),
    add_external_resources(),
    tags$style(
      type = "text/css",
      ".shiny-output-error { visibility: hidden; }",
      ".shiny-output-error:before { visibility: hidden; }"
    ),
    navbarPage(
      title = uiOutput("ui_navtitle"), windowTitle = "BNP_SSD", id="inTabset",

      # data tab ----------------------------------------------------------------
      tabPanel(
        title = span(tagList(icon("table"), inline(uiOutput("ui_nav1")))), value="nav1",
        fluidRow(
          br(),
          column(
            4,
            uiOutput("ui_1choose"),
            # demo data
            inline(uiOutput("ui_1data1")),
            inline(actionLink("infoDemo1", icon = icon("info-circle"), label = NULL)),
            shinyjs::hidden(div(
              id = "infoDemoText1",
              uiOutput("ui_1datahelp1")
            )),
            br(),
            inline(uiOutput("ui_1data2")),
            inline(actionLink("infoDemo2", icon = icon("info-circle"), label = NULL)),
            shinyjs::hidden(div(
              id = "infoDemoText2",
              uiOutput("ui_1datahelp2")
            )),
            # upload csv with data
            br(),
            inline(uiOutput("ui_1csv")),
            inline(actionLink("infoUpload", icon = icon("info-circle"), label = NULL)),
            shinyjs::hidden(div(
              id = "infoUploadText",
              uiOutput("ui_1csvhelp")
            )),
            uiOutput("ui_1csvupload"),

            # input data in DataTable
            inline(uiOutput("ui_1table1")),
            inline(actionLink("infoHands", icon = icon("info-circle"), label = NULL)),
            shinyjs::hidden(div(
              id = "infoHandsText",
              uiOutput("ui_1tablehelp")
            )),
            rhandsontable::rHandsontableOutput("hot")
          ),
          column(
            8,
            uiOutput("ui_1preview"),
            uiOutput("ui_viewupload"),
            # br(), br(),
            # conditionalPanel(
            #   condition = "! (0%in%dim(output.read_data()))",
            #   actionButton('jumpToP2', 'Go to Fit'))
          )
        ),
        div(
          id = "note",
          uiOutput("ui_1note1")
        )
      ),
      # fit tab -----------------------------------------------------------------
      tabPanel(
        title = span(tagList(icon("stats", lib = "glyphicon"), inline(uiOutput("ui_nav2")))), value="nav2",
        fluidRow(
          column(
            4,
            br(),
            wellPanel(
              uiOutput("ui_censor_type"),
              uiOutput("ui_2Conc"),
              uiOutput("ui_2log"),
              uiOutput("ui_2center"),
              br(),
              uiOutput("selectSpecies"),
              uiOutput("ui_2selectNit"),
              uiOutput("ui_2png"),
              div(
                id = "divFormatFit",
                br(),
                inline(uiOutput("ui_2width")),
                inline(uiOutput("ui_2height")),
                inline(uiOutput("ui_2dpi"))
              )
            )
          ),
          column(
            8,
            br(),
            conditionalPanel(
              condition = "output.checkfit",
              htmlOutput("hintFi")
            ),
            conditionalPanel(
              condition = "!output.checkfit",
              actionButton("ui_makeFig", HTML("<b>Fit the model</b>"), icon("paper-plane"), 
                           style="color: #000000; background-color: #f5f5f5; border-color: #d0cdcd")
            ),
            conditionalPanel(
              condition = "output.distPlot1",
              splitLayout(cellWidths = c("70%", "30%"), uiOutput("ui_2plot"), uiOutput("ui_2plot_cpo"))
            ),
            inline(conditionalPanel(
              condition = "output.distPlot1",
              uiOutput("ui_2dlplot")
            )),
            br(), br(),
            splitLayout(cellWidths = c("70%", "30%"), plotOutput("distPlot1"), plotOutput("cpoPlot1")),
            br(),
            conditionalPanel(
              condition = "output.distPlot1",
              textOutput("ui_2textGOF")
            ),
            br(),
            plotOutput("GofPlot"),
            # br(), br(),
            # conditionalPanel(
            #   condition = "output.distPlot1",
            #   actionButton('jumpToP3', 'Go to HC estimation'))
          )
        )
      ),
      # quantile tab --------------------------------------------------------------
      tabPanel(
        title = span(tagList(icon("calculator"), inline(uiOutput("ui_nav3")))), value="nav3",
        fluidRow(
          column(
            4,
            br(),
            wellPanel(
              uiOutput("ui_3selectQ")
            )
          ),
          column(
            8,
            br(),
            conditionalPanel(
              condition = "output.checkpred",
              htmlOutput("hintPred")
            ),
            conditionalPanel(
              condition = "output.QQPlot",
              uiOutput("ui_3plot")
            ),
            inline(conditionalPanel(
              condition = "output.QQPlot",
              uiOutput("ui_3dlplot")
            )),
            br(), br(),
            plotOutput("QQPlot"),
            br(),
            textOutput("ui_3HC"),
            br(),
            DT::dataTableOutput("QQTable"),
            # br(), br(),
            # conditionalPanel(
            #   condition = "output.QQPlot",
            # actionButton('jumpToP4', 'Go to clustering'))
          )
        )
      ),
      # custer tab --------------------------------------------------------------
      tabPanel(
        title = span(tagList(icon("circle-nodes", library = "font-awesome"), inline(uiOutput("ui_nav4")))), value="nav4",
        fluidRow(
          column(
            4,
            br(),
            wellPanel(
              uiOutput("selectGroup"),
            )
          ),
          column(
            8,
            br(),
            conditionalPanel(
              condition = "output.checkclust",
              htmlOutput("hintCl")
            ),
            conditionalPanel(
              condition = "output.ClustPlot",
              uiOutput("ui_4plot")
            ),
            inline(conditionalPanel(
              condition = "output.ClustPlot",
              uiOutput("ui_4dlplot")
            )),
            br(), br(),
            plotOutput("ClustPlot"),
            br(), 
            inline(conditionalPanel(
              condition = "output.ClustPlot",
              textOutput("ui_4clust"),
            )),
            br(),
            DT::dataTableOutput("clustTable")
          )
        )
      )
    )
  )
}

add_external_resources <- function() {
  tagList(tags$link(rel = "stylesheet", type = "text/css", href = "www/style.css"))
}