## Function to plot principal components visually in R - can be used for any dataset
## with columns structured roughly like:
##
## - ID
## - PC1
## - PC2
## - PC3
## - ...etc

plot_principal_components <- function(
    data,
    pc_x = "1",
    pc_y = "2",
    color_by = NULL,
    title = "Principal Component Plot",
    point_size = 3
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required. Install it using ",
      "install.packages('ggplot2')."
    )
  }

  # Check that data is a data frame
  if (!is.data.frame(data)) {
    stop("The supplied 'data' object must be a data frame.")
  }

  # Construct the complete PC column names
  x_column <- paste0("PC", pc_x)
  y_column <- paste0("PC", pc_y)

  # Check that the specified PC columns exist
  missing_pc_columns <- setdiff(
    c(x_column, y_column),
    colnames(data)
  )

  if (length(missing_pc_columns) > 0) {
    stop(
      "The following principal component columns do not exist: ",
      paste(missing_pc_columns, collapse = ", ")
    )
  }

  # Check that color_by exists, if supplied
  if (!is.null(color_by) && !color_by %in% colnames(data)) {
    stop(
      "The specified colour column does not exist: ",
      color_by
    )
  }

  # Create plot without a colour variable
  if (is.null(color_by)) {

    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[x_column]],
        y = .data[[y_column]]
      )
    ) +
      ggplot2::geom_point(size = point_size)

  # Create plot with a colour variable
  } else {

    # Use canonical colours for common superpopulation labels. Values that
    # are not recognised retain their labels and receive stable fallback
    # colours based on their sorted order.
    superpopulation_colours <- c(
      EUR = "#377EB8",
      AFR = "#E41A1C",
      AMR = "#E6AB02",
      CSA = "#238B45",
      SAS = "#E66101",
      EAS = "#56B4E9",
      MID = "#984EA3"
    )

    colour_values <- as.character(data[[color_by]])
    colour_keys <- trimws(toupper(colour_values))
    canonical_colour_values <- ifelse(
      colour_keys == "CAS",
      "CSA",
      colour_keys
    )
    canonical_colour_keys <- unique(canonical_colour_values)
    colour_mapping <- superpopulation_colours[
      intersect(names(superpopulation_colours), canonical_colour_keys)
    ]

    unknown_colour_values <- sort(unique(
      colour_values[
        !is.na(colour_values) &
          !(canonical_colour_values %in% names(superpopulation_colours))
      ]
    ))
    fallback_colours <- c(
      "#F781BF", "#A65628", "#999999", "#66A61E", "#666666", "#A6761D"
    )

    if (length(unknown_colour_values) > length(fallback_colours)) {
      fallback_colours <- grDevices::hcl.colors(
        length(unknown_colour_values),
        palette = "Dark 3"
      )
    }

    if (length(unknown_colour_values) > 0) {
      colour_mapping <- c(
        colour_mapping,
        stats::setNames(
          fallback_colours[seq_along(unknown_colour_values)],
          unknown_colour_values
        )
      )
    }

    mapped_colour_values <- colour_values
    known_values <- !is.na(colour_values) &
      canonical_colour_values %in% names(superpopulation_colours)
    mapped_colour_values[known_values] <- canonical_colour_values[known_values]
    data_for_plot <- data
    data_for_plot[[".plot_colour"]] <- mapped_colour_values

    p <- ggplot2::ggplot(
      data_for_plot,
      ggplot2::aes(
        x = .data[[x_column]],
        y = .data[[y_column]],
        color = .data[[".plot_colour"]]
      )
    ) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::scale_color_manual(values = colour_mapping, na.value = "grey50") +
      ggplot2::labs(color = color_by)
  }

  p <- p +
    ggplot2::labs(
      title = title,
      x = x_column,
      y = y_column
    ) +
    ggplot2::theme_minimal()

  return(p)
}




### Example usage: HGDP + 1KG data (set up as a function so you can easily uncomment and run it if you have the data files downloaded)
### this requires two files: 
### - https://console.cloud.google.com/storage/browser/_details/covid19-hg-public/pca_projection/gnomad_meta_hgdp_tgp_v1.txt
### - https://console.cloud.google.com/storage/browser/_details/covid19-hg-public/pca_projection/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz (unzip this file first)

run_example_pc_plot <- function(save_as_png = FALSE) {
    
    HGDP_1KG_filepath <- "C:/Users/cw1245/Documents/PhD/LocalCoding/HGDP_1KG/"

    metadata <- read.delim(paste(HGDP_1KG_filepath, "pca_projection_gnomad_meta_hgdp_tgp_v1.txt", sep = ""), header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
    scores <- read.table(paste(HGDP_1KG_filepath, "pca_projection_hgdp_tgp_pca_covid19hgi_snps_scores.txt", sep = ""), header = TRUE, sep = "\t")

    ## add labels to scores dataframe (s is the ID column)

    scores_w_labels <- merge(scores, metadata, by = "s", all.x = TRUE)

    ## run pc plot

    plot <- plot_principal_components(scores_w_labels, pc_x = "1", pc_y = "2", color_by = "hgdp_tgp_meta.Genetic.region", title = "HGDP + 1KG PCA Plot")

    print(plot)

    ## save plot as a png
    if(save_as_png) {
      ggplot2::ggsave(filename = paste(HGDP_1KG_filepath, "HGDP_1KG_PCA_plot.png", sep = ""), plot = plot, width = 8, height = 6, dpi = 300)
    }
}


#run_example_pc_plot()
