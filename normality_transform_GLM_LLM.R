# Load packages
packages <- c('readxl', 'rstatix', 'FSA', 'DescTools', 'lme4', 'ggeffects', 'Matrix', 'nlme', 'emmeans')
lapply(packages, library, character.only=TRUE, quietly = TRUE)

# Read in the Excel file
df <- as.data.frame(read_excel('sample_data.xlsx', sheet = 'analysis', col_names=TRUE))

# Convert necessary columns that are IVs to factors 
df$Diet <- as.factor(df$Diet)
df$Drug <- as.factor(df$Drug)
df$Sex <- as.factor(df$Sex)
df$Condition <- as.factor(df$Condition)

# Make variable that holds all DVs which will be analyzed
dep_vars <- names(df)
exclude_vars <- c('Dam','ID', 'Sex', 'Condition', 'Diet', 'Drug', 'Boli')
DVs <- dep_vars[!dep_vars %in% exclude_vars]
print(DVs)

# Indicate which DVs should be assessed by LLM instead of GLM
llm_dvs <- c('Center_point_velocity', 'body_weight')

# Function: remove extreme outliers (3x IQR)
remove_outliers <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 3 * IQR(x, na.rm = TRUE)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

# Function: normality tests + transformations (for N <30 Shaprio Wilk to be used, otherwise Kolmogorov-Smirnov)
normality_test <- function(data, dv) {
  column_data <- as.numeric(data[[dv]])
  transformation <- "none"
  
  # Determine which test to use based on number of entries
  if (length(column_data) < 30) {
    normality <- suppressWarnings(shapiro.test(column_data))
  } else {
    normality <- suppressWarnings(ks.test(column_data, "pnorm", mean = mean(column_data, na.rm = TRUE), sd = sd(column_data, na.rm = TRUE)))
  }
  
  # If p-value is below 0.05, remove outliers and re-test
  if (normality$p.value < 0.05) {
    # Remove extreme outliers
    column_data_sansoutliers <- remove_outliers(column_data)
    if (length(column_data_sansoutliers) < 30) {
      normality <- suppressWarnings(shapiro.test(column_data_sansoutliers))
    } else {
      normality <- suppressWarnings(ks.test(column_data_sansoutliers, "pnorm", mean = mean(column_data_sansoutliers, na.rm = TRUE), sd = sd(column_data_sansoutliers, na.rm = TRUE)))
    }
    transformation <- "sans_outliers"
    
    # If still not normal, do log transformation (with all entries) & redo normality test
    if (normality$p.value < 0.05) {
      column_data_logged <- log(column_data + 1)  # Avoid log(0)
      if (length(column_data_logged) < 30) {
        normality <- suppressWarnings(shapiro.test(column_data_logged))
      } else {
        normality <- suppressWarnings(ks.test(column_data_logged, "pnorm", mean = mean(column_data_logged, na.rm = TRUE), sd = sd(column_data_logged, na.rm = TRUE)))
      }
      transformation <- "Log"
      
      # If still not normal, do sqrt transformation (with all entries) & redo normality test
      if (normality$p.value < 0.05) {
        column_data_sqrt <- sqrt(column_data)
        if (length(column_data_sqrt) < 30) {
          normality <- suppressWarnings(shapiro.test(column_data_sqrt))
        } else {
          normality <- suppressWarnings(ks.test(column_data_sqrt, "pnorm", mean = mean(column_data_sqrt, na.rm = TRUE), sd = sd(column_data_sqrt, na.rm = TRUE)))
        }
        transformation <- "Sqrt"
        
        # If still not normal, do log transformation (sans extreme outliers) & redo normality test
        if (normality$p.value < 0.05) {
          column_data_sansoutliers_logged <- log(remove_outliers(column_data) + 1)  # Avoid log(0)
          if (length(column_data_sansoutliers_logged) < 30) {
            normality <- suppressWarnings(shapiro.test(column_data_sansoutliers_logged))
          } else {
            normality <- suppressWarnings(ks.test(column_data_sansoutliers_logged, "pnorm", mean = mean(column_data_sansoutliers_logged, na.rm = TRUE), sd = sd(column_data_sansoutliers_logged, na.rm = TRUE)))
          }
          transformation <- "log_sans_outliers"
          
          # If still not normal, do sqrt transformation (sans extreme outliers) & redo normality test
          if (normality$p.value < 0.05) {
            column_data_sansoutliers_sqrt <- sqrt(remove_outliers(column_data))
            if (length(column_data_sansoutliers_sqrt) < 30) {
              normality <- suppressWarnings(shapiro.test(column_data_sansoutliers_sqrt))
            } else {
              normality <- suppressWarnings(ks.test(column_data_sansoutliers_sqrt, "pnorm", mean = mean(column_data_sansoutliers_sqrt, na.rm = TRUE), sd = sd(column_data_sansoutliers_sqrt, na.rm = TRUE)))
            }
            transformation <- "sqrt_sans_outliers"
            
            # If still not normal, put nonparametric flag
            if (normality$p.value < 0.05) {
              transformation <- "nonparametric"
            }
          }
        }
      }
    }
  }
  list(transformation = transformation)
}

# Function: to analyze data based on normality, conduct the transformations, and conduct appropriate GLM/LLM
analyze_data <- function(df, dv) {
  # Apply normality test and transformations
  result <- normality_test(df, dv)
  
  # Transform the DV if necessary
  transformed_dv <- switch(result$transformation,
                           "none" = df[[dv]],
                           "sans_outliers" = remove_outliers(df[[dv]]),
                           "Log" = log(df[[dv]] + 1),
                           "Sqrt" = sqrt(df[[dv]]),
                           "log_sans_outliers" = log(remove_outliers(df[[dv]]) + 1),
                           "sqrt_sans_outliers" = sqrt(remove_outliers(df[[dv]])),
                           "nonparametric" = df[[dv]]
  )

  # Print the transformation used
  print(paste("Transformation used for", dv, ":", result$transformation))
  
  # Check if the DV should be assessed with LLM
  if (dv %in% llm_dvs) {
    # Fit the linear mixed model (control for dam (litter effects))
    formula <- as.formula(paste(dv, "~ Condition + (1 | Dam)"))
    mixed.lmer <- lmer(formula, data = df)
    print(paste("Results for DV:", dv))
    print(summary(mixed.lmer))
    # 1-way ANOVA for LLM
    print(anova(mixed.lmer))
    # Get EMM
    emmeans <- ggpredict(mixed.lmer, terms = "Condition")
    print(emmeans)
    # Conduct Tukey posthoc
    posthoc <- emmeans(mixed.lmer, list(pairwise ~ Condition), adjust = "tukey")
    print(posthoc)
    
  } else {
    # Factorial GLM or Kruskal-Wallis H Test
    if (result$transformation != "nonparametric") {
      # Perform Factorial GLM
      glm_result <- glm(transformed_dv ~ Diet * Drug * Sex, data = df)
      print(paste("Results for DV:", dv))
      print(summary(glm_result))
      # 1-way ANOVA
      anova_result <- aov(transformed_dv ~ Condition, data = df)
      print(summary(anova_result))
      # Conduct Tukey posthoc
      tukey_result <- TukeyHSD(anova_result)
      print(tukey_result)
      
    } else {
      # Kruskal-Wallis H Test
      kruskal_result <- kruskal_test(transformed_dv ~ Condition, data = df)
      print(paste("Results for DV:", dv))
      print(kruskal_result)
      # Conduct Tukey posthoc
      dunn_result <- dunnTest(transformed_dv ~ Condition, data = df)
      print(dunn_result)
    }
  }
}

# Make list of dependent variables to iterate thru
dep_vars_list <- DVs

# Apply the function to each dependent variable
for (dv in dep_vars_list) {
  analyze_data(df, dv)
}