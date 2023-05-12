#' @title Color palette names
#' @description Returns all currently valid color panels or -spectra.
#' @return A named list.
#'
#' @details Discrete color panels derive from the ggsci-package. Continuous
#' colorspectra derive from the colorspace-package.
#'
#' @export

all_color_palettes <- confuns::all_color_palettes

#' @rdname all_color_palettes
#' @export
#' @keywords internal
all_color_spectra <- confuns::all_color_spectra

#' @keywords internal
#' @inherit confuns::color_vector
#' @export
color_vector <- confuns::color_vector

#' @keywords internal
#' @inherit confuns::scale_color_add_on
#' @export
scale_color_add_on <- confuns::scale_color_add_on



#' @importFrom magrittr %>%
#'
NULL



# confuns -----------------------------------------------------------------

# a -----------------------------------------------------------------------

#' @title Valid colorpanels & -spectra
#' @export
#' @keywords internal
allColorpalettes <- confuns::all_color_palettes

#' @rdname allColorpalettes
#' @export
#' @keywords internal
allColorspectra <- confuns::all_color_spectra

#' @keywords internal
adapt_reference <- confuns::adapt_reference

#' @keywords internal
adjust_ggplot_params <- confuns::adjust_ggplot_params

#' @keywords internal
are_values <- confuns::are_values

# c -----------------------------------------------------------------------

#' @keywords internal
check_across_subset <- confuns::check_across_subset

#' @keywords internal
check_across_subset2 <- confuns::check_across_subset2

#' @keywords internal
check_across_subset_negate <- confuns::check_across_subset_negate

#' @keywords internal
check_data_frame <- confuns::check_data_frame

#' @keywords internal
check_none_of <- confuns::check_none_of

#' @keywords internal
check_one_of <- confuns::check_one_of

#' @keywords internal
create_progress_bar <- confuns::create_progress_bar



# d -----------------------------------------------------------------------

# imported from grid - see p.R
#drawDetails <- grid::drawDetails


# g -----------------------------------------------------------------------

#' @keywords internal
give_feedback <- confuns::give_feedback

#' @keywords internal
pull_var <- confuns::pull_var



# i -----------------------------------------------------------------------

#' @keywords internal
is_value <- confuns::is_value

#' @keywords internal
is_vec <- confuns::is_vec


# l -----------------------------------------------------------------------
#' @keywords internal
lrename <- confuns::lrename

# m -----------------------------------------------------------------------

#' @keywords internal
make_capital_letters <- confuns::make_capital_letters

#' @keywords internal
make_pretty_df <- confuns::make_pretty_df

#' @keywords internal
make_pretty_name <- confuns::make_pretty_name



# p -----------------------------------------------------------------------

#' @keywords internal
plot_barchart <- confuns::plot_barplot

#' @keywords internal
plot_boxplot <- confuns::plot_boxplot

#' @keywords internal
plot_dot_plot_1d <- confuns::plot_dot_plot_1d

#' @keywords internal
plot_dot_plot_2d <- confuns::plot_dot_plot_2d

#' @keywords internal
plot_density <- confuns::plot_density

#' @keywords internal
plot_histogram <- confuns::plot_histogram

#' @keywords internal
plot_ridgeplot <- confuns::plot_ridgeplot

#' @keywords internal
plot_scatterplot <- confuns::plot_scatterplot

#' @keywords internal
plot_violinplot <- confuns::plot_violinplot

# imported from grid - see p.R

#postDrawDetails <- grid::postDrawDetails

#preDrawDtails <- grid::preDrawDetails

# s -----------------------------------------------------------------------

#' @keywords internal
scollapse <- confuns::scollapse


# v -----------------------------------------------------------------------

#' @keywords internal
vselect <- confuns::vselect





