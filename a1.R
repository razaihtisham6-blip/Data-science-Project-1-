library(shiny)
library(ggplot2)
library(MASS)
library(DT)  # <-- added for the Compare Models table

ui <- fluidPage(
  titlePanel("Least Median of Squares Regression"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      # Explanatory variables (can pick multiple)
      selectizeInput("x_vars", "X Variables", choices = NULL, multiple = TRUE),
      
      # Response variable (pick exactly one)
      selectInput("y_var", "Y Variable", choices = NULL, multiple = FALSE),
      
      sliderInput("sample_size", "Sample size for LMS (if < n)",
                  min = 1, max = 100, value = 50),
      sliderInput("hist_bins", "Histogram bins", min = 5, max = 50, value = 10),
      actionButton("run", "Run Analysis")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions",
                 h4("How to use this app"),
                 p("1. Upload your CSV file"),
                 p("2. Select X and Y variables"),
                 p("3. Adjust sample size for LMS (smaller is faster)"),
                 p("4. Click 'Run Analysis'")
        ),
        tabPanel("Plot",
                 plotOutput("regression_plot", height = 450),
                 conditionalPanel(
                   condition = "input.x_vars && input.x_vars.length == 1",
                   tags$hr(),
                   h4("Single Predictor Fit (LMS vs OLS)"),
                   plotOutput("single_fit_plot", height = 420)
                 ),
                 conditionalPanel(
                   condition = "!(input.x_vars && input.x_vars.length == 1)",
                   div(class = "alert alert-info",
                       "ℹ Line chart appears when exactly one predictor is selected.")
                 )
        ),
        tabPanel("Summary",
                 h3("Least Median of Squares Summary"),
                 verbatimTextOutput("lms_summary"),
                 h3("Ordinary Least Squares Summary"),
                 verbatimTextOutput("ols_summary")),
        tabPanel("Residuals Comparison",
                 plotOutput("residuals_plot", height = 320),
                 plotOutput("residuals_histogram", height = 320)),
        tabPanel("Compare Models",
                 DTOutput("metrics_table"))   # <-- NEW TAB
      )
    )
  )
)

server <- function(input, output, session) {
  # Load data
  data <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })
  
  # Update variable choices
  observe({
    req(data())
    updateSelectizeInput(session, "x_vars", choices = names(data()))
    updateSelectInput(session, "y_var", choices = names(data()))
  })
  
  # Run analysis
  results <- eventReactive(input$run, {
    req(input$x_vars, input$y_var)
    df <- na.omit(data()[, c(input$x_vars, input$y_var), drop = FALSE])
    y <- input$y_var
    form <- as.formula(paste(y, "~", paste(setdiff(input$x_vars, y), collapse = "+")))
    
    ols <- lm(form, data = df)
    lms <- lqs(form, data = df, method = "lms", nsamp = min(input$sample_size, nrow(df)))
    
    list(data = df, y = y, form = form, ols = ols, lms = lms)
  })
  
  # ---- helper for predictions/residuals ----
  preds <- reactive({
    r <- results(); df <- r$data; yv <- df[[r$y]]
    yhat_ols <- fitted(r$ols)
    yhat_lms <- as.vector(model.matrix(r$form, df) %*% coef(r$lms))
    list(y = yv, yhat_ols = yhat_ols, yhat_lms = yhat_lms,
         res_ols = yv - yhat_ols, res_lms = yv - yhat_lms)
  })
  
  # Regression plot (Observed vs Fitted)
  output$regression_plot <- renderPlot({
    p <- preds()
    df_plot <- data.frame(
      observed   = p$y,
      fitted_ols = p$yhat_ols,
      fitted_lms = p$yhat_lms
    )
    
    ggplot(df_plot) +
      geom_point(aes(x = fitted_ols, y = observed, color = "OLS")) +
      geom_point(aes(x = fitted_lms, y = observed, color = "LMedS")) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      scale_color_manual(values = c("OLS" = "blue", "LMedS" = "red"), name = "Model") +
      labs(title = "Observed vs Fitted (OLS & LMedS)",
           x = "Fitted values", y = "Observed y") +
      theme_minimal()
  })
  
  # Single predictor overlay (only when exactly one X is selected)
  output$single_fit_plot <- renderPlot({
    res <- results()
    # only render if exactly one X
    xs <- setdiff(colnames(res$data), res$y)
    validate(need(length(xs) == 1, "Select exactly one predictor."))
    
    xname <- xs[1]
    df <- data.frame(x = res$data[[xname]], y = res$data[[res$y]])
    
    # coefficients
    b_ols <- coef(res$ols)
    b_lms <- coef(res$lms)
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.8) +
      geom_abline(intercept = b_ols[1], slope = b_ols[2], colour = "blue", size = 1) +
      geom_abline(intercept = b_lms[1], slope = b_lms[2], colour = "red", size = 1) +
      labs(title = paste0("Fit on ", xname, " (OLS vs LMedS)"),
           x = xname, y = res$y) +
      theme_minimal()
  })
  
  # Residuals vs Fitted (merged)
  output$residuals_plot <- renderPlot({
    p <- preds()
    df_plot <- data.frame(
      fitted_ols = p$yhat_ols, resid_ols = p$res_ols,
      fitted_lms = p$yhat_lms, resid_lms = p$res_lms
    )
    
    ggplot(df_plot) +
      geom_point(aes(x = fitted_ols, y = resid_ols, color = "OLS")) +
      geom_point(aes(x = fitted_lms, y = resid_lms, color = "LMedS")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("OLS" = "blue", "LMedS" = "red"), name = "Model") +
      labs(title = "Residuals vs Fitted (OLS & LMedS)",
           x = "Fitted values", y = "Residuals") +
      theme_minimal()
  })
  
  # Residuals histogram
  output$residuals_histogram <- renderPlot({
    res <- results()
    df <- res$data
    yhat_ols <- fitted(res$ols)
    yhat_lms <- as.vector(model.matrix(res$form, df) %*% coef(res$lms))
    df$ols_resid <- df[[res$y]] - yhat_ols
    df$lms_resid <- df[[res$y]] - yhat_lms
    
    ggplot(df) +
      geom_histogram(aes(x = ols_resid, fill = "OLS"), alpha = 0.5, bins = input$hist_bins) +
      geom_histogram(aes(x = lms_resid, fill = "LMedS"), alpha = 0.5, bins = input$hist_bins) +
      scale_fill_manual(
        name = "Method",
        values = c("OLS" = "blue", "LMedS" = "red")
      ) +
      labs(title = "Residuals Distribution",
           x = "Residuals", y = "Count") +
      theme_minimal()
  })
  
  # ---- Compare Models table (for the new tab) ----
  output$metrics_table <- renderDT({
    p <- preds()
    y <- p$y
    
    rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
    mae  <- function(e) mean(abs(e), na.rm = TRUE)
    r2   <- function(y, yhat) {
      ss_res <- sum((y - yhat)^2, na.rm = TRUE)
      ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
      1 - ss_res/ss_tot
    }
    
    OLS_RMSE <- rmse(p$res_ols); OLS_MAE <- mae(p$res_ols)
    OLS_R2   <- if (all(is.na(p$yhat_ols))) NA_real_ else r2(y, p$yhat_ols)
    AIC_OLS  <- tryCatch(AIC(results()$ols), error = function(e) NA_real_)
    
    LMS_RMSE <- rmse(p$res_lms); LMS_MAE <- mae(p$res_lms)
    LMS_R2   <- if (all(is.na(p$yhat_lms))) NA_real_ else r2(y, p$yhat_lms)
    
    datatable(
      data.frame(
        Model = c("OLS", "LMedS"),
        RMSE  = c(OLS_RMSE, LMS_RMSE),
        MAE   = c(OLS_MAE,  LMS_MAE),
        R2    = c(OLS_R2,   LMS_R2),
        AIC   = c(AIC_OLS,  NA)
      ),
      options = list(dom = "tip", paging = FALSE),
      rownames = FALSE
    )
  })
  
  # Summary outputs
  output$lms_summary <- renderPrint({
    res <- results()
    cat("Coefficients:\n")
    print(coef(res$lms))
    cat("\nScale estimate (residual median):", res$lms$scale, "\n")
  })
  
  output$ols_summary <- renderPrint({
    res <- results()
    summary(res$ols)
  })
}

shinyApp(ui = ui, server = server)
