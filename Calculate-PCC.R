library(tidyverse)
library(mice)
library(Hmisc)
library(imputeTS)

## ----- Load m6A and APA data  -----
df_APA <- read.csv("DEapa_list.csv")
df_m6A <- read.csv("DEm6a_list.csv")

## ----- Initialize result matrix  -----
results <- matrix(NA, nrow = nrow(df_APA), ncol = 4)
colnames(results) <- c("Gene", "Correlation", "P.Value", "P.Value.Adj")

## ----- Traverse each row  -----
for (i in 1:nrow(df_APA)) {
  if (i %in% 1967) { # Skip rows with insufficient observations
    next
  }
  
  # Extract APA and m6A data for the current row
  apa <- as.numeric(df_APA[i, ])
  m6a <- as.numeric(df_m6A[i, ])
  
  # Find indices with complete cases
  complete_cases <- complete.cases(apa, m6a)
  
  # Extract valid data using indices
  apa <- apa[complete_cases]
  m6a <- m6a[complete_cases]
  
  # Validate sample size
  if (length(apa) < 3 || length(m6a) < 3) {
    cat(sprintf("Row %d (%s): Not enough finite observations (n=%d)\n", 
                i, rownames(df_APA)[i], length(apa)))
    next
  }
  
  # Compute Pearson correlation
  result <- cor.test(apa, m6a, method = "pearson")
  
  # Save results
  results[i, ] <- c(rownames(df_APA)[i], result$estimate, result$p.value)
}

## ----- Convert to data frame -----
results <- as.data.frame(results)

## ----- Extract non-NA p-values and apply FDR correction (BH) -----
valid_p <- !is.na(results$P.Value)
results$P.Value.Adj[valid_p] <- p.adjust(results$P.Value[valid_p], method = "BH")

## ----- Convert to data frame -----
write.csv(results, "cor_results.csv", row.names = TRUE)
