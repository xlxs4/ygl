# Make sure to launch the R process from the
# project top-level directory for here() to work.
# If you can't do that, feel free to play with rprojroot:
# https://github.com/jennybc/here_here#tldr

suppressWarnings(suppressPackageStartupMessages(library(docopt)))


"R script to produce interactive YGL localization plots.

Usage:
  ygl_analysis_sc_localization.R [--plots]
  ygl_analysis_sc_localization.R --plots [-q | --no-color]
  ygl_analysis_sc_localization.R [-q | --no-color]
  ygl_analysis_sc_localization.R (-h | --help)
  ygl_analysis_sc_localization.R (-v | --version)

Options:
  -h --help     Show this screen.
  -v --version     Show version.

" -> doc
arguments <- docopt(doc, version = "YGL-Loc 0.1")

suppressWarnings(suppressPackageStartupMessages({
    library(biomaRt)
    library(dplyr)
    library(htmlwidgets)
    library(ggplot2)
    library(glue)
    library(here)
    library(logger)
    library(plotly)
    library(R.matlab)
    library(stringr)
    library(tibble)
    # For saving the tibbles.
    library(arrow)
    # For the flexboard plots.
    library(rmarkdown)
    library(flexdashboard)
}))


if (arguments$q) {
    log_threshold(SUCCESS)
} else {
    if (!arguments$no_color) {
        log_layout(layout_glue_colors)
    }
}

datafiles_pathname <- "singlecell_datafiles"

extract_loc_means <- function(sgd_systematic) {
    scfilename <- paste(sgd_systematic, "_MMS_high.mat", sep = "")
    scfilepath <- here::here(
        datafiles_pathname,
        scfilename
    )

    data <- R.matlab::readMat(scfilepath)
    locmean <- data$scData[3][[1]]

    return(locmean)
}

extract_scfilenames <- function() {
    names <- list()
    i <- 1
    for (filename in list.files(here::here(datafiles_pathname))) {
        names[i] <- (unlist(strsplit(filename, split = "_", fixed = TRUE))[1])
        i <- i + 1
    }
    return(unlist(names))
}

convert_to_genename <- function(sgd_systematic) {
    sgd_path <- here::here("SGD_features.tab")
    sgd <- readLines(sgd_path)
    return(strsplit(grep(sgd_systematic, sgd, value = TRUE), "\t")[[1]][5])
}

loc_means_to_df <- function(loc_means, shapes, n_tpoints = 40) {
    names <- rep(shapes, n_tpoints)
    timepoints <- rep(1:n_tpoints, each = length(shapes))
    vals <- unlist(lapply(loc_means, function(i) unlist(i, recursive = TRUE)))

    df <- tibble::tibble(names, timepoints)
    df$loc_probabilities <- vals

    return(df)
}

plot_df <- function(df) {
    p <- df %>%
        ggplot(aes(
            x = timepoints,
            y = loc_probabilities,
            group = names,
            color = names
        )) +
        geom_line() +
        labs(
            x = "Observation Timepoints",
            y = "Localization Probability",
            color = "Shapes"
        ) +
        theme(
            axis.title = element_text(
                size = 15,
                color = "firebrick"
            ),
            axis.text = element_text(
                color = "dodgerblue",
                size = 12
            ),
            legend.title = element_text(
                size = 15,
                color = "firebrick",
                face = "italic"
            )
        ) +
        ggtitle(str_glue("{df$genename} Localization"))
    l <- plotly::ggplotly(p)
    return(l)
}

save_plot <- function(plot, path) {
    htmlwidgets::saveWidget(plot, path)
}

main <- function() {
    plot_dir <- here::here("plots")
    tibble_dir <- here::here("tibbles")

    ggplot2::theme_set(ggplot2::theme_bw())

    for (sgd_systematic in extract_scfilenames()) {
        logger::log_info("{sgd_systematic} found!")
        genename <- convert_to_genename(sgd_systematic)
        logger::log_info("Gene name: {genename}")

        logger::log_info("Extracting localization means...")
        loc_means <- extract_loc_means(sgd_systematic)

        shapes <- c(
            "Periphery",
            "Structure",
            "Punctate",
            "Disk",
            "Corona",
            "Homogeneous"
        )
        logger::log_info("Converting to tibble...")
        loc_means <- loc_means_to_df(loc_means, shapes)
        loc_means$genename <- genename

        logger::log_info("Saving tibble as feather binary...")
        tibble_path <- paste(genename, ".feather", sep = "")
        arrow::write_feather(loc_means, here::here(tibble_dir, tibble_path))

        if (arguments$plots) {
            logger::log_info("Plotting...")
            p <- plot_df(loc_means)

            plot_path <- paste(genename, ".html", sep = "")
            logger::log_info("Saving plot...")
            save_plot(p, here::here(plot_dir, plot_path))
        }
    }
}

logger::log_info("Executing...")
main()
log_success("The script finished running successfully!")