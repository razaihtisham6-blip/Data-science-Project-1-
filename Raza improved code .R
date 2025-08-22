# Shiny App: Robust Regression (OLS, LMS, Huber) with Extras
# Features:
#  - Upload CSV, select X/Y
#  - Fit OLS, LMS (Least Median of Squares), and Huber regression
#  - Diagnostic plots and residuals
#  - Model comparison (RMSE, MAE, R2)
#  - Bootstrap resampling (coefficients) and k-fold CV
#  - User-friendly alerts, validation, and downloads

library(shiny)
library(shinythemes)
library(ggplot2)
library(MASS)
library(DT)             # for the Compare Models table
library(shinycssloaders) # for loading spinners

# User Interface (UI) definition
ui <- fluidPage(
  theme = shinytheme("flatly"), # Apply a theme 
  titlePanel("Least Median of Squares Regression"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      # Explanatory variables
      selectizeInput("x_vars", "X Variables", choices = NULL, multiple = TRUE),
      # Response variable (only one)
      selectInput("y_var", "Y Variable", choices = NULL, multiple = FALSE),
      sliderInput("sample_size", "Sample size for LMS (if < n)",
                  min = 1, max = 100, value = 50),
      sliderInput("hist_bins", "Histogram bins", min = 5, max = 50, value = 10),
      actionButton("run", "Run Analysis"), br(), br(), # br() used to give line space
      uiOutput("xvar_warning")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions",
                 h4(class = "alert alert-success","ℹ Each group member will put their own instructions in their individual submission!!!"),
                 p("Please proceed to next tab to use the application.")
        ),
        tabPanel("Plot",
                 withSpinner(plotOutput("regression_plot", height = 450), type = 5, size = 0.8),
                 downloadButton("dl_regplot", "Download Observed vs Fitted (PNG)"), br(),br(),
                 conditionalPanel(
                   condition = "input.x_vars && input.x_vars.length == 1",
                   tags$hr(),
                   h4("Single Predictor Fit (LMS vs OLS)"),
                   withSpinner(plotOutput("single_fit_plot", height = 420), type = 5, size = 0.8)
                 ),
                 conditionalPanel(
                   condition = "!(input.x_vars && input.x_vars.length == 1)",
                   div(class = "alert alert-info",
                       "ℹ Line chart appears when exactly one predictor is selected.")
                 )
        ),
        tabPanel("Summary",
                 h3("Least Median of Squares Summary"),  
                 withSpinner(verbatimTextOutput("lms_summary"), type = 5, size = 0.8),
                 h3("Ordinary Least Squares Summary"),   
                 withSpinner(verbatimTextOutput("ols_summary"), type = 5, size = 0.8),
                 h3("Huber (BFGS) Summary"),             
                 withSpinner(verbatimTextOutput("huber_summary"), type = 5, size = 0.8)
        ),
        tabPanel("Residuals Comparison",
                 withSpinner(plotOutput("residuals_plot", height = 360), type = 5, size = 0.8),
                 withSpinner(plotOutput("residuals_histogram", height = 360), type = 5, size = 0.8)
        ),
        tabPanel("Compare Models",
                 withSpinner(DTOutput("metrics_table"), type = 5, size = 0.8),  br(),
                 downloadButton("dl_metrics", "Download metrics (CSV)")    
        ),
        tabPanel("Bootstrap",            
                 numericInput("boot_B", "Bootstrap resamples (B)", value = 1000, min = 200, step = 100),
                 actionButton("boot_run", "Run Bootstrap"),
                 tags$hr(),
                 h4("Bootstrap Coefficients (Point Estimates)"),
                 withSpinner(DTOutput("boot_ci_table"), type = 5, size = 0.8),
                 br(),
                 downloadButton("dl_bootci", "Download bootstrap coefficients (CSV)")
        ),
        tabPanel("Cross-Validation",               
                 numericInput("cv_k", "k-folds", value = 5, min = 3, max = 10, step = 1),
                 actionButton("cv_run", "Run k-fold CV"),
                 tags$hr(),
                 withSpinner(DTOutput("cv_table"), type = 5, size = 0.8)
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  data <- reactive({                             #CSV handling
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
  # Show a visible (Bootstrap) warning near the button if X is empty when Run is clicked
  output$xvar_warning <- renderUI({
    req(input$run)  # only trigger after button click
    if (length(input$x_vars) == 0) {
      div(
        class = "alert alert-info",
        "ℹ Please select at least one X variable before running the analysis.")}
  })
  # Update variable choices
  observe({
    df <- data()
    req(df)
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]     #Only keep numerical data 
    updateSelectizeInput(session, "x_vars", choices = num_cols)
    updateSelectInput(session, "y_var", choices = num_cols)
  })
  #Huber (BFGS) helpers
  huber_rho <- function(r, delta) {
    ifelse(abs(r) <= delta, 0.5 * r^2, delta * (abs(r) - 0.5 * delta))
  }
  fit_huber <- function(form, data) {
    mf <- model.frame(form, data)
    y  <- model.response(mf)
    X  <- model.matrix(form, mf)
    
    b0 <- tryCatch(coef(lm(form, data = data)), error = function(e) rep(0, ncol(X)))    # start at OLS
    if (is.null(names(b0))) names(b0) <- colnames(X)[seq_along(b0)]
    if (length(b0) != ncol(X)) {                 # align to X columns if needed
      b_tmp <- rep(0, ncol(X)); names(b_tmp) <- colnames(X)
      cmn <- intersect(names(b0), colnames(X)); b_tmp[cmn] <- b0[cmn]; b0 <- b_tmp
    }
    
    # Set Huber threshold delta using 1.345 × MAD of initial residuals (robust 95% efficient scale).
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
  
  # Core analysis pipeline: builds formula, fits OLS/LMS/Huber after input validation
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
    validate(                     # Ensure sample size is large enough (at least predictors + 1 rows)
      need(nrow(df) >= length(input$x_vars) + 1,
           "Not enough rows for the selected number of predictors.")
    )
    y <- input$y_var
    form <- as.formula(paste(y, "~", paste(setdiff(input$x_vars, y), collapse = "+")))
    
    ols <- lm(form, data = df)
    lms <- lqs(form, data = df, method = "lms",
               nsamp = min(input$sample_size, nrow(df)))
    huber <- fit_huber(form, df)                  #fit Huber (BFGS)
    list(data = df, y = y, form = form, ols = ols, lms = lms, huber = huber)
  })
  preds <- reactive({                       #helper for predictions/residuals
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
  output$dl_metrics <- downloadHandler(        #Download handler for metrics CSV 
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
  
  output$dl_regplot <- downloadHandler(                 #Download handler for Observed vs Fitted PNG
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
  output$lms_summary <- renderPrint({                         # Summary outputs
    res <- results()
    cat("Coefficients:\n")
    print(coef(res$lms))
    cat("\nScale estimate (residual median):", res$lms$scale, "\n")
  })
  output$ols_summary <- renderPrint({
    res <- results()
    summary(res$ols)
  })
  output$huber_summary <- renderPrint({               #Huber (BFGS) summary
    res <- results()
    cat("Coefficients (Huber, BFGS):\n")
    print(res$huber$coef)
    cat("\nEffective delta:", res$huber$delta,
        "\nConverged:", res$huber$converged, "\n")
  }) 
  #Bootstrap
  .fit_ols <- function(form, data) coef(lm(form, data = data))
  .fit_lms <- function(form, data) coef(MASS::lqs(form, data = data, method = "lms"))
  
  boot_pairs <- function(form, data, fit_fun, B = 1000) {
    n <- nrow(data)
    pnames <- colnames(model.matrix(form, data))
    out <- matrix(NA_real_, nrow = B, ncol = length(pnames))
    colnames(out) <- pnames
    for (b in seq_len(B)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      cf  <- suppressWarnings(try(fit_fun(form, data[idx, , drop = FALSE]), silent = TRUE))
      if (inherits(cf, "try-error")) next
      out[b, names(cf)] <- unname(cf)
    }
    as.data.frame(out, check.names = FALSE)
  }
  
  boot_res <- eventReactive(input$boot_run, {
    r <- results()
    validate(need(!is.null(r), "Run Analysis first."))
    draws_ols <- boot_pairs(r$form, r$data, .fit_ols, B = input$boot_B)
    draws_lms <- boot_pairs(r$form, r$data, .fit_lms, B = input$boot_B)
    list(ols = draws_ols, lms = draws_lms, ref = r)
  })
  
  boot_coef_table_data <- reactive({                # Table with only Model, Term, Estimate
    br <- boot_res()
    validate(need(!is.null(br), "Click 'Run Bootstrap' first."))
    r  <- br$ref
    block_df <- function(label, ref_coef) {
      terms <- names(ref_coef)
      data.frame(
        Model    = label,
        Term     = terms,
        Estimate = as.numeric(ref_coef[terms]),
        check.names = FALSE
      )
    }
    ols_block <- block_df("OLS",   coef(r$ols))
    lms_block <- block_df("LMS", coef(r$lms))
    tbl <- rbind(ols_block, lms_block)
    tbl$Estimate <- round(tbl$Estimate, 5)
    tbl
  })
  
  output$boot_ci_table <- renderDT({
    datatable(boot_coef_table_data(), options = list(dom = "tip", paging = FALSE), rownames = FALSE)
  })
  
  output$dl_bootci <- downloadHandler(
    filename = function() paste0("bootstrap_coefficients_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(boot_coef_table_data(), file, row.names = FALSE)
  )
  
  #k-fold Cross-Validation
  kfold_cv <- function(form, data, k = 5, fit_fun, pred_fun = NULL) {
    n   <- nrow(data)
    idx <- sample(rep(1:k, length.out = n))
    rmses <- numeric(k)
    
    for (i in 1:k) {
      train <- data[idx != i, , drop = FALSE]
      test  <- data[idx == i, , drop = FALSE]
      fit_try <- try(fit_fun(form, train), silent = TRUE)
      if (inherits(fit_try, "try-error")) { rmses[i] <- NA_real_; next }
      fit <- fit_try
      yhat <- try(
        if (is.null(pred_fun)) predict(fit, newdata = test) else pred_fun(fit, form, test),
        silent = TRUE
      )
      if (inherits(yhat, "try-error")) { rmses[i] <- NA_real_; next }
      rmses[i] <- sqrt(mean((test[[all.vars(form)[1]]] - yhat)^2))
    }
    mean(rmses, na.rm = TRUE)
  }
  
  pred_lqs <- function(fit, form, newdata) {
    as.vector(model.matrix(form, newdata) %*% coef(fit))
  }
  
  cv_res <- eventReactive(input$cv_run, {
    r <- results()
    validate(need(!is.null(r), "Run Analysis first."))
    k <- input$cv_k
    
    ols_rmse <- kfold_cv(r$form, r$data, k, fit_fun = function(f, d) lm(f, d))
    lms_rmse <- kfold_cv(r$form, r$data, k,
                         fit_fun = function(f, d) MASS::lqs(f, d, method = "lms"),
                         pred_fun = pred_lqs)
    
    data.frame(Model = c("OLS", "LMS"), CV_RMSE = c(ols_rmse, lms_rmse))
  })
  
  output$cv_table <- renderDT({
    datatable(cv_res(), options = list(dom = "tip", paging = FALSE), rownames = FALSE)
  })
}

shinyApp(ui = ui, server = server)