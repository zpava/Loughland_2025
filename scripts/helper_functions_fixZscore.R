library(pathfindR)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# # Function to calculate z-scores from adjusted p-values
# calc_zscore_from_adjp_logfc <- function(adj_pval, logFC) {
#   # Convert adjusted p-values to absolute z-scores
#   abs_z_scores <- qnorm(adj_pval / 2, lower.tail = FALSE)
#   
#   # Handle any infinite values (very small p-values)
#   max_finite_z <- max(abs_z_scores[is.finite(abs_z_scores)])
#   abs_z_scores[is.infinite(abs_z_scores)] <- max_finite_z
#   
#   # Assign direction based on logFC
#   z_scores <- abs_z_scores * sign(logFC)
#   
#   return(z_scores)
# }
# 
# library(tidyverse)
# library(ComplexHeatmap)
# library(circlize)
# 
# # Function to calculate average Z-score for each enriched term
# # note: low_p and high_p are both adjusted p values calculated via bonferroni.
# 
# # Function to create the heatmap
# create_pathfindr_heatmap2 <- function(enriched_terms, sample_groups, low_color, mid_color, high_color) {
#   #Create direction
#   enriched_terms = enriched_terms%>%
#     mutate(
#       dir= case_when(Down_regulated != "" & Up_regulated == "" ~ -1,
#                      Down_regulated == "" & Up_regulated != "" ~  1,
#                      TRUE ~ 0),
#       Down_Reg_No = str_count(Down_regulated,  ","),
#       Up_Reg_No = str_count(Up_regulated, ","),
#       final_dir = case_when(dir == 0 & Down_Reg_No > Up_Reg_No ~ -1,
#                             dir == 0 & Down_Reg_No < Up_Reg_No ~ 1, TRUE ~ dir)) 
#   
#   # Calculate Z-scores
#   enriched_terms$z.score <- calc_zscore_from_adjp_logfc(enriched_terms$lowest_p, pathway_outcome$final_dir)
#   
#   #Filtering significant z.scores
#   enriched_terms = enriched_terms %>%
#                     filter(z.score >= 2 | z.score <=-2)
#   
#    
#   # Create a matrix of z-scores
#   z_score_matrix <- enriched_terms %>%
#     pivot_wider(names_from = Cluster_name, values_from = z.score, values_fn = list) %>%
#     mutate(across(-Term_Description, ~sapply(., function(x) if(length(x) > 0) mean(unlist(x)) else NA))) %>%
#     column_to_rownames("Term_Description")
#   
#   # Ensure all sample groups are present
#   missing_samples <- setdiff(names(sample_groups), colnames(z_score_matrix))
#   if (length(missing_samples) > 0) {
#     z_score_matrix[, missing_samples] <- NA
#   }
#   
#   # Order columns according to sample_groups
#   z_score_matrix <- z_score_matrix[, names(sample_groups)]
#   
#   # Order rows based on average z-score
#   row_means <- rowMeans(z_score_matrix, na.rm = TRUE)
#   z_score_matrix <- z_score_matrix[order(row_means, decreasing = TRUE), ]
#   
#   # Create the color function
#   col_fun <- colorRamp2(
#     breaks = c(min(z_score_matrix, na.rm = TRUE), 0, max(z_score_matrix, na.rm = TRUE)),
#     colors = c(low_color, mid_color, high_color)
#   )
# 
#   # Create the heatmap
#   heatmap <- Heatmap(as.matrix(z_score_matrix),
#                      name = "Z-score",
#                      col = col_fun,
#                      width = unit(50, "mm"),
#                      cluster_rows = FALSE,
#                      cluster_columns = FALSE,
#                      show_row_names = TRUE,
#                      show_column_names = TRUE,
#                      row_names_gp = gpar(fontsize = 8),
#                      column_names_gp = gpar(fontsize = 10),
#                      column_names_rot = 45,
#                      column_title = "Cluster Day",
#                      row_title = "Enriched Terms",
#                      row_title_gp = gpar(fontsize = 12),
#                      column_title_gp = gpar(fontsize = 12))
#   
#   return(heatmap)
# }
# Function to calculate z-scores from adjusted p-values
calc_zscore_from_adjp_logfc <- function(adj_pval, logFC) {
  # Convert adjusted p-values to absolute z-scores
  abs_z_scores <- qnorm(adj_pval / 2, lower.tail = FALSE)
  
  # Handle any infinite values (very small p-values)
  max_finite_z <- max(abs_z_scores[is.finite(abs_z_scores)])
  abs_z_scores[is.infinite(abs_z_scores)] <- max_finite_z
  
  # Assign direction based on logFC
  z_scores <- abs_z_scores * sign(logFC)
  
  return(z_scores)
}

# Function to create the heatmap
create_pathfindr_heatmap2 <- function(enriched_terms, sample_groups, low_color, mid_color, high_color) {
  # Create direction
  enriched_terms <- enriched_terms %>%
    mutate(
      dir = case_when(
        Down_regulated != "" & Up_regulated == "" ~ -1,
        Down_regulated == "" & Up_regulated != "" ~ 1,
        TRUE ~ 0
      ),
      Down_Reg_No = ifelse(Down_regulated == "", 0, str_count(Down_regulated, ",") + 1),
      Up_Reg_No = ifelse(Up_regulated == "", 0, str_count(Up_regulated, ",") + 1),
      final_dir = case_when(
        dir == 0 & Down_Reg_No > Up_Reg_No ~ -1,
        dir == 0 & Down_Reg_No < Up_Reg_No ~ 1,
        TRUE ~ dir
      )
    )
  
  # Calculate Z-scores
  enriched_terms$z.score <- calc_zscore_from_adjp_logfc(enriched_terms$lowest_p, enriched_terms$final_dir)
  
  # Filtering significant z.scores
  enriched_terms <- enriched_terms %>%
    filter(z.score >= 2 | z.score <= -2)
  
  # Create a matrix of z-scores
  # First get all unique terms
  unique_terms <- unique(enriched_terms$Term_Description)
  
  # Create an empty matrix with all possible combinations
  all_clusters <- unique(enriched_terms$Cluster_name)
  z_score_matrix <- matrix(NA, 
                           nrow = length(unique_terms), 
                           ncol = length(all_clusters),
                           dimnames = list(unique_terms, all_clusters))
  
  # Fill in the z-scores
  for(i in 1:nrow(enriched_terms)) {
    term <- enriched_terms$Term_Description[i]
    cluster <- enriched_terms$Cluster_name[i]
    z_score_matrix[term, cluster] <- enriched_terms$z.score[i]
  }
  
  # Convert to data frame
  z_score_matrix <- as.data.frame(z_score_matrix)
  
  # Ensure all sample groups are present
  missing_samples <- setdiff(names(sample_groups), colnames(z_score_matrix))
  if (length(missing_samples) > 0) {
    z_score_matrix[, missing_samples] <- NA
  }
  
  # Order columns according to sample_groups
  z_score_matrix <- z_score_matrix[, names(sample_groups)]
  
  # Order rows based on average absolute z-score
  row_means <- rowMeans(as.matrix(z_score_matrix), na.rm = TRUE)
  z_score_matrix <- z_score_matrix[order(row_means, decreasing = TRUE), ]
  
  # Create the color function
  col_fun <- colorRamp2(
    breaks = c(min(z_score_matrix, na.rm = TRUE), 0, max(z_score_matrix, na.rm = TRUE)),
    colors = c(low_color, mid_color, high_color)
  )
  
  # Create the heatmap
  heatmap <- Heatmap(
    as.matrix(z_score_matrix),
    name = "Z-score",
    col = col_fun,
    width = unit(50, "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    column_title = "Cluster Day",
    row_title = "Enriched Terms",
    row_title_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 12)
  )
  
  return(heatmap)
}

##Plotting enrichment

library(pathfindR)
library(ggplot2)

custom_enrichment_bubble_chart <- function(result_df, top_terms = 10, plot_by_cluster = FALSE, 
                                           num_bubbles = 4, even_breaks = TRUE,
                                           low_color = "#f5efef", high_color = "red") {
  message("Plotting the enrichment bubble chart")
  necessary <- c("Term_Description", "Fold_Enrichment", "lowest_p", 
                 "Up_regulated", "Down_regulated")
  if (!all(necessary %in% colnames(result_df))) {
    stop("The input data frame must have the columns:\n", 
         paste(necessary, collapse = ", "))
  }
  if (!is.logical(plot_by_cluster)) {
    stop("`plot_by_cluster` must be either TRUE or FALSE")
  }
  if (!is.numeric(top_terms) & !is.null(top_terms)) {
    stop("`top_terms` must be either numeric or NULL")
  }
  if (!is.null(top_terms)) {
    if (top_terms < 1) {
      stop("`top_terms` must be > 1")
    }
  }
  result_df <- result_df[order(result_df$lowest_p), ]
  if (!is.null(top_terms)) {
    if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
      keep_ids <- tapply(result_df$ID, result_df$Cluster, 
                         function(x) {
                           x[seq_len(min(top_terms, length(x)))]
                         })
      keep_ids <- unlist(keep_ids)
      result_df <- result_df[result_df$ID %in% keep_ids, ]
    }
    else if (top_terms < nrow(result_df)) {
      result_df <- result_df[seq_len(top_terms), ]
    }
  }
  num_genes <- vapply(result_df$Up_regulated, function(x) length(unlist(strsplit(x, 
                                                                                 ", "))), 1)
  num_genes <- num_genes + vapply(result_df$Down_regulated, 
                                  function(x) length(unlist(strsplit(x, ", "))), 1)
  result_df$Term_Description <- factor(result_df$Term_Description, 
                                       levels = rev(unique(result_df$Term_Description)))
  log_p <- -log10(result_df$lowest_p)
  g <- ggplot2::ggplot(result_df, ggplot2::aes(.data$Fold_Enrichment, 
                                               .data$Term_Description))
  g <- g + ggplot2::geom_point(ggplot2::aes(color = log_p, 
                                            size = num_genes), na.rm = TRUE)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10), 
                          axis.text.y = ggplot2::element_text(size = 10), plot.title = ggplot2::element_blank())
  g <- g + ggplot2::xlab("Fold Enrichment")
  g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  g <- g + ggplot2::labs(size = "# genes", color = expression(-log[10](p)))
  if (max(num_genes) < num_bubbles) {
    g <- g + ggplot2::scale_size_continuous(breaks = seq(0, 
                                                         max(num_genes)))
  }
  else {
    if (even_breaks) {
      brks <- base::seq(0, max(num_genes), round(max(num_genes)/(num_bubbles + 
                                                                   1)))
    }
    else {
      brks <- base::round(base::seq(0, max(num_genes), 
                                    length.out = num_bubbles + 1))
    }
    g <- g + ggplot2::scale_size_continuous(breaks = brks)
  }
  g <- g + ggplot2::scale_color_gradient(low = low_color, high = high_color)
  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster ~ ., scales = "free_y", 
                                 space = "free", drop = TRUE)
  }
  else if (plot_by_cluster) {
    message("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
  }
  return(g)
}




