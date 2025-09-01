library(shiny)
library(ggplot2)
library(MASS)
library(DT)  # for the Compare Models table

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
                 downloadButton("dl_regplot", "Download Observed vs Fitted (PNG)"),
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
                 verbatimTextOutput("ols_summary"),
                 h3("Huber (BFGS) Summary"),
                 verbatimTextOutput("huber_summary")
        ),
        tabPanel("Residuals Comparison",
                 plotOutput("residuals_plot", height = 320),
                 plotOutput("residuals_histogram", height = 320)
        ),
        tabPanel("Boxplot",
                 plotOutput("residuals_boxplot", height = 380)
        ),
        tabPanel("Compare Models",
                 DTOutput("metrics_table"),
                 br(),
                 downloadButton("dl_metrics", "Download metrics (CSV)")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ---- Defensive CSV handling ----
  data <- reactive({
    req(input$file)
    df <- tryCatch({
      read.csv(input$file$datapath, stringsAsFactors = FALSE)
    }, error = function(e) {
      showNotification(
        "Couldn't read the CSV. Ensure it's a comma-separated .csv with a header row.",
        type = "error", duration = 7
      )
      return(NULL)
    })
    validate(
      need(!is.null(df), "Uploaded file could not be read."),
      need(is.data.frame(df) && nrow(df) > 0 && ncol(df) > 0,
           "Uploaded file appears empty or invalid.")
    )
    df
  })
  
  # Update variable choices
  observe({
    df <- data()
    req(df)
    updateSelectizeInput(session, "x_vars", choices = names(df))
    updateSelectInput(session, "y_var", choices = names(df))
  })
  
  # ---- Huber (BFGS) helpers (summary only) ----
  huber_rho <- function(r, delta) {
    ifelse(abs(r) <= delta, 0.5 * r^2, delta * (abs(r) - 0.5 * delta))
  }
  fit_huber <- function(form, data) {
    mf <- model.frame(form, data)
    y  <- model.response(mf)
    X  <- model.matrix(form, mf)
    
    # start at OLS
    b0 <- tryCatch(coef(lm(form, data = data)), error = function(e) rep(0, ncol(X)))
    if (is.null(names(b0))) names(b0) <- colnames(X)[seq_along(b0)]
    if (length(b0) != ncol(X)) {                 # align to X columns if needed
      b_tmp <- rep(0, ncol(X)); names(b_tmp) <- colnames(X)
      cmn <- intersect(names(b0), colnames(X)); b_tmp[cmn] <- b0[cmn]; b0 <- b_tmp
    }
    
    # robust scale for delta (1.345 * MAD of initial residuals)
    r0 <- as.vector(y - X %*% b0)
    s  <- suppressWarnings(mad(r0, center = 0, constant = 1, na.rm = TRUE))
    if (!is.finite(s) || s <= 0) s <- 1
    delta_eff <- 1.345 * s
    
    obj <- function(b) {
      r <- as.vector(y - X %*% b)
      sum(huber_rho(r, delta_eff))
    }
    
    opt <- optim(par = b0, fn = obj, method = "BFGS", control = list(maxit = 1000))
    list(coef = setNames(opt$par, colnames(X)),
         delta = delta_eff,
         converged = (opt$convergence == 0))
  }
  
  # Run analysis
  results <- eventReactive(input$run, {
    req(input$x_vars, input$y_var)
    df0 <- data()
    req(df0)
    
    # Subset to selected columns; handle missing columns gracefully
    df <- tryCatch({
      na.omit(df0[, c(input$x_vars, input$y_var), drop = FALSE])
    }, error = function(e) {
      validate(need(FALSE, "Selected variables not found in the uploaded data."))
      NULL
    })
    validate(
      need(!is.null(df) && nrow(df) > 0, "No usable rows after subsetting/removing NAs.")
    )
    
    y <- input$y_var
    form <- as.formula(paste(y, "~", paste(setdiff(input$x_vars, y), collapse = "+")))
    
    ols <- lm(form, data = df)
    lms <- lqs(form, data = df, method = "lms",
               nsamp = min(input$sample_size, nrow(df)))
    huber <- fit_huber(form, df)
    
    list(data = df, y = y, form = form, ols = ols, lms = lms, huber = huber)
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
      geom_point(aes(x = fitted_lms, y = observed, color = "LMS")) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      scale_color_manual(values = c("OLS" = "blue", "LMS" = "red"), name = "Model") +
      labs(title = "Observed vs Fitted (OLS & LMS)",
           x = "Fitted values", y = "Observed y") +
      theme_minimal()
  })
  
  # Single predictor overlay (only when exactly one X is selected)
  output$single_fit_plot <- renderPlot({
    res <- results()
    xs <- setdiff(colnames(res$data), res$y)
    validate(need(length(xs) == 1, "Select exactly one predictor."))
    
    xname <- xs[1]
    df <- data.frame(x = res$data[[xname]], y = res$data[[res$y]])
    
    b_ols <- coef(res$ols)
    b_lms <- coef(res$lms)
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.8) +
      geom_abline(intercept = b_ols[1], slope = b_ols[2], colour = "blue", size = 1) +
      geom_abline(intercept = b_lms[1], slope = b_lms[2], colour = "red", size = 1) +
      labs(title = paste0("Fit on ", xname, " (OLS vs LMS)"),
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
      geom_point(aes(x = fitted_lms, y = resid_lms, color = "LMS")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("OLS" = "blue", "LMS" = "red"), name = "Model") +
      labs(title = "Residuals vs Fitted (OLS & LMS)",
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
      geom_histogram(aes(x = lms_resid, fill = "LMS"), alpha = 0.5, bins = input$hist_bins) +
      scale_fill_manual(
        name = "Method",
        values = c("OLS" = "blue", "LMS" = "red")
      ) +
      labs(title = "Residuals Distribution",
           x = "Residuals", y = "Count") +
      theme_minimal()
  })
  
  # ---- Box-and-whisker plot of residuals ----
  output$residuals_boxplot <- renderPlot({
    p <- preds()
    df_box <- data.frame(
      Method   = factor(c(rep("OLS", length(p$res_ols)),
                          rep("LMS", length(p$res_lms))),
                        levels = c("OLS", "LMS")),
      Residual = c(p$res_ols, p$res_lms)
    )
    ggplot(df_box, aes(x = Method, y = Residual, fill = Method)) +
      geom_boxplot(outlier.alpha = 0.4) +
      scale_fill_manual(values = c("OLS" = "blue", "LMS" = "red")) +
      labs(title = "Residuals: Box-and-Whisker (OLS vs LMS)",
           x = NULL, y = "Residual") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  # Compare Models table
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
    
    LMS_RMSE <- rmse(p$res_lms); LMS_MAE <- mae(p$res_lms)
    LMS_R2   <- if (all(is.na(p$yhat_lms))) NA_real_ else r2(y, p$yhat_lms)
    
    datatable(
      data.frame(
        Model = c("OLS", "LMS"),
        RMSE  = c(OLS_RMSE, LMS_RMSE),
        MAE   = c(OLS_MAE,  LMS_MAE),
        R2    = c(OLS_R2,   LMS_R2)
      ),
      options = list(dom = "tip", paging = FALSE),
      rownames = FALSE
    )
  })
  
  # ---- Download handler for metrics CSV ----
  output$dl_metrics <- downloadHandler(
    filename = function() paste0("model_metrics_", Sys.Date(), ".csv"),
    content = function(file) {
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
      
      LMS_RMSE <- rmse(p$res_lms); LMS_MAE <- mae(p$res_lms)
      LMS_R2   <- if (all(is.na(p$yhat_lms))) NA_real_ else r2(y, p$yhat_lms)
      
      metrics_df <- data.frame(
        Model = c("OLS", "LMS"),
        RMSE  = c(OLS_RMSE, LMS_RMSE),
        MAE   = c(OLS_MAE,  LMS_MAE),
        R2    = c(OLS_R2,   LMS_R2)
      )
      
      write.csv(metrics_df, file, row.names = FALSE)
    }
  )
  
  # ---- Download handler for Observed vs Fitted PNG ----
  output$dl_regplot <- downloadHandler(
    filename = function() paste0("observed_vs_fitted_", Sys.Date(), ".png"),
    content = function(file) {
      p <- preds()
      df_plot <- data.frame(
        observed   = p$y,
        fitted_ols = p$yhat_ols,
        fitted_lms = p$yhat_lms
      )
      
      plt <- ggplot(df_plot) +
        geom_point(aes(x = fitted_ols, y = observed, color = "OLS")) +
        geom_point(aes(x = fitted_lms, y = observed, color = "LMS")) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        scale_color_manual(values = c("OLS" = "blue", "LMS" = "red"), name = "Model") +
        labs(title = "Observed vs Fitted (OLS & LMS)",
             x = "Fitted values", y = "Observed y") +
        theme_minimal()
      
      ggsave(filename = file, plot = plt, width = 8, height = 6, dpi = 300)
    },
    contentType = "image/png"
  )
  
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
  
  # Huber (BFGS) summary
  output$huber_summary <- renderPrint({
    res <- results()
    cat("Coefficients (Huber, BFGS):\n")
    print(res$huber$coef)
    cat("\nEffective delta:", res$huber$delta,
        "\nConverged:", res$huber$converged, "\n")
  })
}

shinyApp(ui = ui, server = server)
