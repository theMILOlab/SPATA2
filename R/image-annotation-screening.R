









# b -----------------------------------------------------------------------

#' @title Bin barcode-spots by angle
#'
#' @description Bins barcode-spots according to their angle towards the position
#' specified with argument \code{center}.
#'
#' @param center Numeric vector of length two that is named. Value named \emph{x}
#' provides position on the x-axis. Value named \emph{y} provides position on
#' the y-axis.
#' @param min_bins_circle Numeric value or NULL. Indiates the minimum
#' number of circle bins the angle bin groups must have in order not
#' to be renamed or removed. Ignored if NULL.
#' @param remove Logical value. If TRUE, barcode-spots that fall into
#' angle bins that do feature less circle bins than input for
#' argument \code{min_bins_circle} are removed. Ignored if
#' \code{min_bins_circle} is NULL.
#' @param rename Logical value. If TRUE, barcode-spots that fall into
#' angle bins that feature less circle bins than input for argument
#' \code{min_bins_circle} are renamed to \emph{'Outside'}. Ignored if \code{min_bins_circle} is NULL.
#' Set \code{remove} to FALSE in order not to remove the renamed
#' barcode-spots.
#'
#' @inherit bin_by_area params
#'
#' @export
bin_by_angle <- function(coords_df,
                         center,
                         n_bins_angle = 1,
                         angle_span = c(0,360),
                         min_bins_circle = NULL,
                         rename = FALSE,
                         remove = FALSE,
                         drop = TRUE){

  confuns::is_vec(x = angle_span, mode = "numeric", of.length = 2)

  angle_span <- c(from = angle_span[1], to = angle_span[2])

  range_span <- base::range(angle_span)

  if(angle_span[1] == angle_span[2]){

    stop("Invalid input for argument `angle_span`. Must contain to different values.")

  } else if(base::min(angle_span) < 0 | base::max(angle_span) > 360){

    stop("Input for argument `angle_span` must range from 0 to 360.")

  }

  # compute angle
  prel_angle_df <-
    dplyr::group_by(.data = coords_df, barcodes) %>%
    dplyr::mutate(
      angle = compute_angle_between_two_points(
        p1 = c(x = x, y = y),
        p2 = center
      )
    )

  # create angle bins

  if(angle_span[["from"]] > angle_span[["to"]]){

    range_vec <- c(
      angle_span[["from"]]:360,
      0:angle_span[["to"]]
    )

    nth <- base::floor(base::length(range_vec)/n_bins_angle)

    bin_list <- base::vector(mode = "list", length = n_bins_angle)

    for(i in 1:n_bins_angle){

      if(i == 1){

        sub <- 1:nth

      } else {

        sub <- ((nth*(i-1))+1):(nth*i)

      }

      bin_list[[i]] <- range_vec[sub]

    }

    if(base::any(base::is.na(bin_list[[n_bins_angle]]))){

      bin_list[[(n_bins_angle)-1]] <-
        c(bin_list[[(n_bins_angle-1)]], bin_list[[n_bins_angle]]) %>%
        rm_na()

      bin_list[[n_bins_angle]] <- NULL

    }

    all_vals <- purrr::flatten_dbl(bin_list)

    bin_list[[n_bins_angle]] <-
      c(bin_list[[n_bins_angle]], range_vec[!range_vec %in% all_vals])

    prel_angle_bin_df <-
      dplyr::ungroup(prel_angle_df) %>%
      dplyr::filter(base::round(angle) %in% range_vec) %>%
      dplyr::mutate(
        angle_round = base::round(angle),
        bins_angle = ""
      )

    bin_names <- base::character(n_bins_angle)

    for(i in base::seq_along(bin_list)){

      angles <- bin_list[[i]]

      bin_names[i] <-
        stringr::str_c(
          "(", angles[1], ",", utils::tail(angles,1), "]"
        )

      prel_angle_bin_df[prel_angle_bin_df$angle_round %in% angles, "bins_angle"] <-
        bin_names[i]

    }

    prel_angle_bin_df$angle_round <- NULL

    prel_angle_bin_df$bins_angle <-
      base::factor(
        x = prel_angle_bin_df$bins_angle,
        levels = bin_names
      )

  } else {

    range_vec <- range_span[1]:range_span[2]

    sub <-
      base::seq(
        from = 1,
        to = base::length(range_vec),
        length.out = n_bins_angle+1
      ) %>%
      base::round()

    breaks <- range_vec[sub]

    prel_angle_bin_df <-
      dplyr::ungroup(prel_angle_df) %>%
      dplyr::filter(base::round(angle) %in% range_vec) %>%
      dplyr::mutate(
        bins_angle = base::cut(x = base::abs(angle), breaks = breaks)
      )

  }

  bins_angle_levels <-
    base::levels(prel_angle_bin_df[["bins_angle"]]) %>%
    c("Core", ., "Outside") %>%
    base::unique()

  bins_circle_levels <-
    base::levels(prel_angle_bin_df[["bins_circle"]]) %>%
    c("Core", ., "Outside") %>%
    base::unique()

  # add to input data.frame
  angle_df <-
    dplyr::left_join(
      x = coords_df,
      y = prel_angle_df[,c("barcodes", "angle")],
      by = "barcodes"
    ) %>%
    dplyr::left_join(
      x = .,
      y = prel_angle_bin_df[,c("barcodes", "bins_angle")],
      by = "barcodes"
    ) %>%
    dplyr::mutate(
      bins_angle = base::as.character(bins_angle),
      bins_angle = dplyr::case_when(
        !base::is.na(bins_angle) ~ bins_angle,
        TRUE ~ "Outside"
      ),
      bins_angle = base::factor(bins_angle, levels = bins_angle_levels)
    )

  # denote barcodes to remove due to insufficient number of circles in angle_bin
  if(!base::is.numeric(min_bins_circle)){

    min_bins_circle <- 0

  }

  bins_to_keep <-
    dplyr::select(angle_df ,bins_circle, bins_angle) %>%
    dplyr::distinct() %>%
    dplyr::group_by(bins_angle) %>%
    dplyr::tally() %>%
    dplyr::filter(n >= {{min_bins_circle}}) %>%
    dplyr::pull(bins_angle) %>%
    base::as.character()

  if(base::isTRUE(rename)){

    angle_df <-
      dplyr::mutate(
        .data = angle_df,
        bins_angle = base::as.character(bins_angle),
        bins_angle = dplyr::case_when(
          bins_circle == "Core" ~ "Core",
          bins_circle == "Outside" ~ "Outside",
          !bins_angle %in% {{bins_to_keep}} ~ "Outside",
          TRUE ~ bins_angle
        ),
        bins_circle = base::as.character(bins_circle),
        bins_circle = dplyr::case_when(
          bins_angle == "Outside" ~ "Outside",
          TRUE ~ bins_circle
        ),
        bins_circle = base::factor(bins_circle, levels = bins_circle_levels),
        bins_angle = base::factor(bins_angle, levels = bins_angle_levels)
      )

  }

  if(base::isTRUE(remove)){

    angle_df <- dplyr::filter(angle_df, bins_angle %in% {{bins_to_keep}})

  }

  if(base::isTRUE(drop)){

    angle_df <-
      dplyr::mutate(
        .data = angle_df,
        bins_angle = base::droplevels(bins_angle),
        bins_circle = base::droplevels(bins_circle)
      )

  }

  return(angle_df)

}


#' @title Bin barcode-spots by area extension
#'
#' @description Bins barcode-spots by consecutively expanding a polygon.
#'
#' @param coords_df The coordinates data.frame whose barcode-spots are supposed
#' to be binned.
#' @param area_df Data.frame with variables \emph{x} and \emph{y} describing the
#' vertices of the polygon that encircles the area based on which the barcode-spots
#' are binned. E.g. slot @@area of \code{ImageAnnotation}-objects.
#' @param remove Character or logical. If character, denotes circle bins that
#' are removed. If TRUE, bins \emph{'Core' and 'Outside'} are removed. If FALSE,
#' ignored.
#' @param drop Logical value. If TRUE, unused levels of the \emph{bins_circle}
#' variables are dropped.
#' @inherit imageAnnotationScreening params
#'
#' @export
bin_by_area <- function(coords_df,
                        area_df,
                        binwidth,
                        n_bins_circle,
                        remove = FALSE,
                        drop = TRUE){

  n_bins_circle <- base::max(n_bins_circle)

  circle_names <- stringr::str_c("Circle", 1:n_bins_circle, sep = " ")

  circles <-
    purrr::set_names(
      x = c((1:n_bins_circle)*binwidth),
      nm = circle_names
    )

  binwidth_vec <- c("Core" = 0, circles)

  areas <-
    purrr::imap(
      .x = binwidth_vec,
      .f = ~ buffer_area(df = area_df, buffer = .x)
    )

  # create new variable. Default is 'Outside'.
  # values will be overwritten with every additional loop
  coords_df$bins_circle <- "Outside"

  for(area in base::names(areas)){

    area_df <- areas[[area]]

    coords_df$pt_in_plg <-
      sp::point.in.polygon(
        point.x = coords_df$x,
        point.y = coords_df$y,
        pol.x = area_df$x,
        pol.y = area_df$y
      )

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        bins_circle = dplyr::case_when(
          # if bins_circle is NOT 'Outside' it has already bin binned
          bins_circle == "Outside" & pt_in_plg %in% c(1,2) ~ {{area}},
          TRUE ~ bins_circle
        )
      )
  }

  bin_levels <- c(base::names(binwidth_vec), "Outside")

  out_df <-
    dplyr::mutate(
      .data = coords_df,
      bins_circle = base::factor(x = bins_circle, levels = bin_levels),
      # as.numeric uses level order, -1 cause 'Core' is first, should be Circle 1
      bins_order = (base::as.numeric(bins_circle) - 1),
      pt_in_plg = NULL
    )

  if(base::is.character(remove)){

    out_df <- dplyr::filter(out_df, !bins_circle %in% {{remove}})

  } else if(base::isTRUE(remove)){

    out_df <- dplyr::filter(out_df, !bins_circle %in% c("Core", "Outside"))

  }

  if(base::isTRUE(drop)){

    out_df <- dplyr::mutate(out_df, bins_circle = base::droplevels(bins_circle))

  }

  return(out_df)

}


bin_by_expansion <- bin_by_area




#' @title Buffer area
#'
#' @description Buffers the area of a polygon.
#'
#' @param df Data.frame with variables \emph{x} and \emph{y} describing the
#' vertices of the polygon that encircles the area based on which the barcode-spots
#' are binned. E.g. slot @@area of \code{ImageAnnotation}-objects.
#'
#' @param buffer The distance by which to consecutively expand the
#' area that covers the image annotation screening. Given to argument
#' \code{dist} of function \code{sf::st_buffer()}. (See details for more.)
#'
#' @export
buffer_area <- function(df, buffer){

  frow <- df[1, c("x", "y")] %>% base::as.numeric()
  lrow <- df[base::nrow(df), c("x", "y")] %>% base::as.numeric()

  if(!base::identical(frow, lrow)){

    df <- close_area_df(df)

  }

  area_grown <-
    sf::st_polygon(x = list(base::as.matrix(df[,c("x", "y")]))) %>%
    sf::st_buffer(dist = buffer) %>%
    base::as.matrix() %>%
    base::as.data.frame() %>%
    magrittr::set_colnames(value = c("x", "y")) %>%
    tibble::as_tibble()

  return(area_grown)

}



# c -----------------------------------------------------------------------

#' @title Close area encircling
#'
#' @description "Closes" the area described by the vertices of \code{df} by
#' adding the starting point (first row) to the end of the data.frame.
#'
#' @export
close_area_df <- function(df){

  fr <- base::as.numeric(df[1,])
  lr <- base::as.numeric(df[base::nrow(df), ])

  if(!base::identical(x = fr, y = lr)){

    df[base::nrow(df) + 1, ] <- df[1,]

  }

  return(df)

}


#' @title Compute angle between two points
#'
#' @description Computes the angle between two points. 0° is aligned
#' with the y-axis.
#'
#' @param p1,p2 Numeric vectors of length two, named \emph{x} and \emph{y}.
#'
#' @export
compute_angle_between_two_points <- function(p1, p2){

  angle <- base::atan2(y = (p2["y"] - p1["y"]), x = (p2["x"] - p1["x"])) * 180/pi

  if(angle >= 0){

    angle <- 360 - angle

  } else {

    angle <- base::abs(angle)

  }

  angle <- angle + 90

  if(angle >= 360){

    angle <- angle - 360

  }

  angle <- angle + 180

  if(angle > 360){

    angle <- angle - 360

  }


  return(angle)

}



# e -----------------------------------------------------------------------



# g -----------------------------------------------------------------------


#' @title Obtain IAS results (data.frame)
#'
#' @description Extracts (and filters) the summarized IAS results in form
#' of a data.frame
#'
#' @param var_pval,var_eval Character value. Specifies the p-value- and the
#' evaluation variable based on which the thresholds are applied.
#' @param threshold_pval,threshold_eval Numeric value. Used to filter
#' the output accordingly.
#'
#' @inherit object_dummy params
#' @inherit add_models params
#'
#' @return Data.frame.
#' @export
#'

setGeneric(name = "getSmrdResultsDf", def = function(object, ...){

  standardGeneric(f = "getSmrdResultsDf")

})

#' @rdname getSmrdResultsDf
#' @export
setMethod(
  f = "getSmrdResultsDf",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        var_pval = "pvalue_median",
                        var_eval = "corr_median",
                        threshold_pval = 1,
                        threshold_eval = 0,
                        model_subset = NULL,
                        model_remove = NULL){

    rdf <-
      dplyr::filter(
        .data = object@results_smrd,
        !!rlang::sym(var_pval) <= {{threshold_pval}} &
        !!rlang::sym(var_eval) >= {{threshold_eval}}
        )

    if(base::is.character(model_subset)){

      keep <- confuns::vselect(base::unique(rdf[["models"]]), dplyr::contains(model_subset))

      rdf <- dplyr::filter(rdf, models %in% {{keep}})

    }

    if(base::is.character(model_remove)){

      keep <- confuns::vselect(base::unique(rdf[["models"]]), -dplyr::contains(model_subset))

      rdf <- dplyr::filter(rdf, models %in% {{keep}})

    }

    return(rdf)

  }
)

#' @title Obtain IAS results (variable names)
#'
#' @description Extracts (and filters) the summarized IAS results in form
#' of a character vector of variable names.
#'
#' @param var_arrange The evaulation variable based on which the variables
#' are arranged before being returned. Defaults to input of argument \code{var_eval}.
#' @inherit getVarDf params
#'
#' @details After the filtering the variable names are arranged according
#' to their values of \code{var_eval}.
#'
#' @return Data.frame.
#' @export
#'

setGeneric(name = "getVarNames", def = function(object, ...){

  standardGeneric(f = "getVarNames")

})

#' @rdname getVarNames
#' @export
setMethod(
  f = "getVarNames",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        var_pval = "pvalue_median",
                        var_eval = "corr_median",
                        var_arrange = var_eval,
                        threshold_pval = 1,
                        threshold_eval = 0,
                        model_subset = NULL,
                        model_remove = NULL){

    getSmrdResultsDf(
      object = object,
      var_pval = var_pval,
      var_eval = var_eval,
      threshold_pval = threshold_pval,
      threshold_eval = threshold_eval,
      model_subset = model_subset,
      model_remove = model_remove
      ) %>%
      dplyr::group_by(variables) %>%
      dplyr::slice_max(order_by = !!rlang::sym(var_eval)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(!!rlang::sym(var_arrange))) %>%
      dplyr::pull(variables) %>%
      base::unique()

  }
)


#' @title Obtain IAS screending data.frame
#'
#' @description Bins and annotates barcode-spots in the same way that
#' \code{imageAnnotationScreening()} does.
#'
#' @param remove_circle_bins,remove_angle_bins,remame_angle_bins Logical values.
#' Given to the corresponding arguments of \code{bin_by_area()}
#' and \code{bin_by_angle()}.
#'
#' @inherit imageAnnotationScreening params
#' @inherit joinWith params
#'
#' @export

getImageAnnotationScreeningDf <- function(object,
                                          id,
                                          distance = NA_integer_,
                                          binwidth = NA_integer_,
                                          n_bins_circle = NA_integer_,
                                          angle_span = c(0,360),
                                          n_bins_angle = 1,
                                          variables = NULL,
                                          normalize = TRUE,
                                          smooth = FALSE,
                                          smooth_span = 0.2,
                                          remove_circle_bins = FALSE,
                                          remove_angle_bins = FALSE,
                                          rename_angle_bins = FALSE,
                                          drop = TRUE,
                                          summarize_by_circles = TRUE,
                                          summarize_with = "mean",
                                          verbose = TRUE){

  hlpr_assign_arguments(object)

  input_list <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle
    )

  distance <- input_list$distance
  n_bins_circle <- input_list$n_bins_circle
  binwidth  <- input_list$binwidth

  max_circles <- base::max(n_bins_circle)
  min_circles <- base::min(n_bins_circle)

  img_ann <- getImageAnnotation(object = object, id = id, add_image = FALSE)

  img_ann_center <- getImageAnnotationCenter(object, id = id)

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::select(barcodes, x, y)

  if(base::length(drop) == 1){ drop <- base::rep(drop, 2)}

  ias_df <-
    bin_by_area(
      coords_df = coords_df,
      area_df = img_ann@area,
      binwidth = binwidth,
      n_bins_circle = max_circles,
      remove = remove_circle_bins,
      drop = drop[1]
    ) %>%
    bin_by_angle(
      center = img_ann_center,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      min_bins_circle = min_circles,
      rename = rename_angle_bins,
      remove = remove_angle_bins,
      drop = drop[2]
    )

  if(base::is.character(variables)){

    if(base::length(normalize) == 1){

      normalize <- base::rep(normalize, 2)

    }

    var_df <-
      joinWithVariables(
        object = object,
        spata_df = getSpataDf(object),
        variables = variables,
        smooth = smooth,
        smooth_span = smooth_span,
        normalize = normalize[1]
      )

    ias_df <-
      dplyr::left_join(
        x = ias_df,
        y = var_df,
        by = "barcodes"
      )

    if(base::isTRUE(normalize[2])){

      ias_df <-
        dplyr::mutate(
          .data = ias_df,
          dplyr::across(
            .cols = dplyr::all_of(variables),
            .fns = confuns::normalize
          )
        )

    }

  }

  if(base::isTRUE(summarize_by_circles)){

    confuns::give_feedback(
      msg = glue::glue("Summarizing by bins."),
      verbose = verbose
    )

    ias_df <-
      dplyr::group_by(.data = ias_df, bins_circle, bins_order) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::any_of(variables),
          .fns = summarize_formulas[[summarize_with]]
        )
      )

  }

  return(ias_df)

}



# i -----------------------------------------------------------------------

#' @title Implementation of the IAS-algorithm
#'
#' @description Screens the sample for numeric variables that stand
#' in meaningful, spatial relation to annotated structures/areas.
#' For a detailed explanation on how to define the parameters \code{distance},
#' \code{n_bins_circle}, \code{binwidth}, \code{angle_span} and \code{n_bins_angle}
#' see details section.
#'
#' @inherit getImageAnnotatation params
#' @param variables Character vector. All numeric variables (meaning genes,
#' gene-sets and numeric features) that are supposed to be included in
#' the screening process.
#' @param distance Numeric value. Specifies the distance from the border of the
#' image annotation to the \emph{horizon} in the periphery up to which the screening
#' is conducted. (See details for more.)
#' @param binwidth Numeric value. The width of the circular bins to which
#' the barcode-spots are assigned. We recommend to set it equal to the distance
#' between the barcode-spots: \code{binwidth = getBarcodeSpotsDistance(object)}.
#' (See details for more.)
#' @param n_bins_circle Numeric value or vector of length 2. Specifies how many times the area is buffered with the value
#' denoted in \code{binwidth}.
#'  (See details for more.)
#' @param angle_span Numeric vector of length 2. Confines the area screened by
#' an angle span relative to the center of the image annotation.
#'  (See details fore more.)
#' @param n_bins_angle Numeric value. Number of bins that are created by angle.
#' (See details for more.)
#'
#' @param summarize_with Character value. Either \emph{'mean'} or \emph{'median'}.
#' Specifies the function with which the bins are summarized.
#'
#' @inherit add_models params
#' @inherit argument_dummy params
#' @inherit buffer_area params
#'
#' @return An object of class \code{ImageAnnotationScreening}. See documentation
#' with \code{?ImageAnnotationScreening} for more information.
#'
#' @seealso createImageAnnotations()
#'
#' @details In conjunction with argument \code{id} which provides the
#' ID of the image annotation of interest the arguments \code{distance},
#' \code{binwidth}, \code{n_bins_circle}, \code{angle_span} and \code{n_bins_angle} can be used
#' to specify the exact area that is screened as well as the resolution of the screening.
#'
#' \bold{How the algorithm works:} During the IAS-algorithm the barcode spots are
#' binned according to their localisation to the image annotation. Every bin's mean
#' expression of a given gene is then aligned in an ascending order - mean expression
#' of bin 1, mean expression of bin 2, ... up to the last bin, the bin with the
#' barcode-spots that lie farest away from the image annotation. This allows to infer
#' the gene expression changes in relation to the image annotation and
#' to screen for genes whose expression changes resemble specific biological
#' behaviors. E.g. linear ascending: gene expression increases linearly with
#' the distance to the image annotation. E.g. immediate descending: gene expression
#' is high in close proximity to the image annotation and declines logarithmically
#' with the distance to the image annotation.
#'
#' \bold{How circular binning works:}
#' To bin barcode-spots according to their localisation to the image annotation
#' three parameters are required:
#'
#'  \itemize{
#'    \item{\code{distance}: The distance from the border of the image annotation to
#'     the \emph{horizon} in the periphery up to which the screening is conducted. Unit
#'     of the distance is pixel as is the unit of the image.
#'     }
#'     \item{\code{binwidth}: The width of every bin. Unit is pixel.}
#'     \item{\code{n_bins_circle}: The number of bins that are created.}
#'  }
#'
#' Regarding parameter \code{n_bins_circle}: The suffix \code{_circle} is used for one
#' thing to emphasize that bins are created in a circular fashion around the image
#' annotation (although the shape of the polygon that was created to encircle the
#' image annotation is maintained). Additionally, the suffix is needed to delineate
#' it from argument \code{n_bins_angle} which can be used to increase the
#' resolution of the screening.
#'
#' These three parameters stand in the following relation to each other:
#'
#'  \enumerate{
#'   \item{\code{n_bins_circle} = \code{distance} / \code{binwidth}}
#'   \item{\code{distance} = \code{n_bins_circle} * \code{binwidth}}
#'   \item{\code{binwidth} = \code{distance} / \code{n_bins_circle}}
#'  }
#'
#' Therefore, only two of the three arguments must be specified as the remaining
#' one is calculated. We recommend to stick to the first option: Specifying
#' \code{distance} and \code{binwidth} and letting the function calculate
#' \code{n_bins_circle}.
#'
#' Once the parameters are set and calculated the polygon that is used to
#' define the borders of the image annotation (the one you draw with
#' \code{createImageAnnotation()}) is repeatedly expanded by the distance indicated
#' by parameter \code{binwidth}. The number of times this expansion is
#' repeated is equal to the parameter \code{n_bins_circle}. Every time the
#' polygon is expanded, the newly enclosed barcode-spots are binned (grouped)
#' and the bin is given a number that is equal to the number of the expansion.
#' Thus, barcode-spots that are adjacent to the image annotation are binned into
#' bin 1, barcode spots that lie a distance of \code{binwidth} away are binned into
#' bin 2, etc.
#'
#' Note that the function \code{plotSurfaceIas()} allows to visually check
#' if your input results in the desired screening.
#'
#' \bold{How the screening works:} For every gene that is included in the
#' screening process every bin's mean expression is calculated and then
#' aligned in an ascending order - mean expression of bin 1, mean expression
#' of bin 2, ... up to the last bin, namely the bin with the barcode-spots that lie
#' farest away from the image annotation. This allows to infer
#' the gene expression changes in relation to the image annotation and
#' to screen for genes whose expression changes resemble specific biological
#' behaviors. The gene expression change is fitted to every model that is included.
#' (Use \code{showModels()} to visualize the predefined models of \code{SPATA2}).
#' A gene-model-fit is evaluated twofold:
#'
#'  \itemize{
#'    \item{Residuals area over the curve}: The area under the curve (AUC) of the
#'    residuals between the inferred expression changes and the model is calculated,
#'    normalized against the number of bins and then subtracted from 1.
#'    \item{Pearson correlation}: The inferred expression changes is correlated
#'    with the model. (Correlation as well as the corresponding p-value depend
#'    on the number of bins!)
#'   }
#'
#' Eventually, the mean of the RAOC and the Correlation for every gene-model-fit
#' is calculated and stored as the IAS-Score.
#'
#' @export
imageAnnotationScreening <- function(object,
                                     id,
                                     variables,
                                     distance = NA_integer_,
                                     binwidth = NA_integer_,
                                     n_bins_circle = NA_integer_,
                                     angle_span = c(0,360),
                                     n_bins_angle = 1,
                                     summarize_with = "mean",
                                     method_padj = "fdr",
                                     model_subset = NULL,
                                     model_remove = NULL,
                                     model_add = NULL,
                                     verbose = NULL,
                                     ...){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Starting image annotation screening.",
    verbose = verbose
  )

  img_ann <- getImageAnnotation(object, id = id)

  input_list <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle
    )

  distance <- input_list$distance
  n_bins_circle <- input_list$n_bins_circle
  binwidth  <- input_list$binwidth

  ias_df <-
    getImageAnnotationScreeningDf(
      object = object,
      id = id,
      variables = variables,
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      remove_circle_bins = TRUE,
      remove_angle_bins = TRUE,
      drop = FALSE,
      summarize_by = FALSE
    )

  bins_angle <- base::levels(ias_df$bins_angle)

  ias_df <- dplyr::mutate(ias_df, bins_angle = base::droplevels(bins_angle))

  bins_angle_remaining <- base::levels(ias_df$bins_angle)

  min_bins_circle <- base::min(n_bins_circle)
  max_bins_circle <- base::max(n_bins_circle)

  # test model input
  model_df <-
    create_model_df(
      input = max_bins_circle,
      var_order = "bins_order",
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add,
      verbose = verbose
    )

  # model fitting
  n_total <- base::length(bins_angle_remaining)

  time_start <- base::Sys.time()
  bin_duration <- NULL
  fn_envir <- base::environment()

  confuns::give_feedback(
    msg = "Fitting models by bin.",
    verbose = verbose
  )

  results <-
    purrr::map_df(
      .x = bins_angle_remaining,
      .f = function(bin){

        start_bin <- base::Sys.time()

        nth <- base::which(bins_angle_remaining == bin)

        confuns::give_feedback(
          msg = glue::glue("Working on bin {bin}. ({nth}/{n_total})"),
          verbose = TRUE
        )

        bin_dur <- base::get(x = "bin_duration", envir = fn_envir)

        if(!base::is.null(bin_dur)){

          # -1 cause nth bin is yet to be screened
          n_remaining <- n_total - (nth-1)

          dur_sec <-
            base::as.numeric(bin_dur * n_remaining) %>%
            base::round(digits = 2)

          dur_min <- base::round(dur_sec/60, digits = 2)
          dur_hours <- base::round(dur_sec/3600, digits = 2)

          est_end <- base::Sys.time() + dur_sec

          msg <- glue::glue("Estimated end of screening: {est_end}.")

          confuns::give_feedback(msg = msg, verbose = verbose)

        }

        bin_angle_df <-
          dplyr::filter(ias_df, bins_angle == {{bin}}) %>%
          dplyr::group_by(bins_order) %>%
          dplyr::summarise(
            dplyr::across(
              .cols = dplyr::all_of(variables),
              .fns = summarize_formulas[[summarize_with]]
            )
          ) %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(variables),
            names_to = "variables",
            values_to = "values"
          )

        shifted_df_with_models <-
          dplyr::left_join(
            x = bin_angle_df,
            y = model_df,
            by = "bins_order"
          ) %>%
          dplyr::arrange(variables) %>%
          shift_for_evaluation(var_order = "bins_order")

        results <-
          base::suppressWarnings({

            evaluate_model_fits(
              input_df = shifted_df_with_models,
              var_order = "bins_order",
              with_corr = TRUE,
              with_raoc = TRUE
            )

          }) %>%
          dplyr::mutate(bins_angle = {{bin}})

        end_bin <- base::Sys.time()

        base::assign(
          x = "bin_duration",
          value = base::difftime(end_bin, start_bin, units = "secs"),
          envir = fn_envir
        )

        return(results)

      }
    )

  confuns::give_feedback(
    msg = "Finished model fitting.",
    verbose = verbose
  )


  # assemble output and summarize
  confuns::give_feedback(
    msg = "Summarizing output.",
    verbose = verbose
    )

  ias_out <-
    ImageAnnotationScreening(
      angle_span = angle_span,
      binwidth = binwidth,
      coords = getCoordsDf(object),
      distance = distance,
      img_annotation = getImageAnnotation(object, id = id),
      n_bins_angle = n_bins_angle,
      n_bins_circle = n_bins_circle,
      results = results,
      sample = object@samples
    ) %>%
    summarizeIAS(method_padj = method_padj)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(ias_out)

}

#' @title Convert image annotation to segmentation
#'
#' @description Converts one or more image annotations to a binary
#' segmentation variable in the feature data.frame.
#' @param ids Character vector. Specifies the image annotation(s) of interest.
#' Barcode-spots that fall into the area of these annotations are labeled
#' with the input for argument \code{inside}.
#' @param segmentation_name Character value. The name of the new segmentation variable.
#' @param inside Character value. The group name for the barcode-spots that
#' are located inside the area of the image annotation(s).
#' @param outside Character value. The group name for the barcode-spots that
#' are located outside the area of the image annotation(s).
#' @param overwrite Logical. Set to TRUE to overwrite existing variables with
#' the same name.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata object.
#' @export
#'
imageAnnotationToSegmentation <- function(object,
                                          ids,
                                          segmentation_name,
                                          inside = "inside",
                                          outside = "outside",
                                          overwrite = FALSE){

  confuns::are_values("inside", "outside", mode = "character")

  confuns::check_none_of(
    input = segmentation_name,
    against = getFeatureNames(object),
    ref.against = "names of the feature data",
    overwrite = overwrite
  )

  bcsp_inside <- getImageAnnotationBarcodes(object, ids = ids)

  fdata <-
    getFeatureDf(object) %>%
    dplyr::mutate(
      {{segmentation_name}} := dplyr::case_when(
        condition = barcodes %in% {{bcsp_inside}} ~ {{inside}},
        TRUE ~ {{outside}}
      ),
      {{segmentation_name}} := base::factor(
        x = !!rlang::sym(segmentation_name),
        levels = c(inside, outside)
      )
    )

  object <- setFeatureDf(object, feature_df = fdata)

  return(object)

}

# m -----------------------------------------------------------------------


make_angle_bins <- function(n){

  n <- base::as.integer(n)

  mltp <- 360/n

  breaks <- 0:n * mltp

  base::cut(x = 0:360, breaks = breaks) %>%
    base::levels()

}

# p -----------------------------------------------------------------------

pick_vars <- function(df, input, order_by, neg_log){

  if(base::is.list(input)){

    var_names <-
      purrr::keep(.x = input, .p = is.character) %>%
      purrr::flatten_chr() %>%
      base::unique()

    out_df_names <-
      dplyr::filter(df, variables %in% {{var_names}})

    n <-
      purrr::keep(.x = input, .p = is.numeric) %>%
      purrr::flatten() %>%
      base::as.numeric()

    if(base::length(n) == 0){

      out_df <- out_df_names

    } else if(base::length(n) == 1){

      if(base::isTRUE(neg_log)){

        out_df <-
          dplyr::group_by(df, models) %>%
          dplyr::slice_max(order_by = !!rlang::sym(order_by), n = n) %>%
          dplyr::ungroup()

      } else {

        out_df <-
          dplyr::group_by(df, models) %>%
          dplyr::slice_min(order_by = !!rlang::sym(order_by), n = n) %>%
          dplyr::ungroup()

      }

      out_df <-
        base::rbind(out_df, out_df_names) %>%
        dplyr::distinct()

    } else {

      stop("Numeric input for argument `label_vars` must be of length 1.")

    }


  } else if(base::is.character(input)){

    confuns::check_one_of(
      input = input,
      against = df$variables
    )

    out_df <- dplyr::filter(df, variables %in% {{input}})

  } else if(base::is.numeric(input)){

    confuns::is_value(x = input, mode = "numeric")

    if(base::isTRUE(neg_log)){

      out_df <-
        dplyr::group_by(df, models) %>%
        dplyr::slice_max(order_by = !!rlang::sym(order_by), n = input) %>%
        dplyr::ungroup()

    } else {

      out_df <-
        dplyr::group_by(df, models) %>%
        dplyr::slice_min(order_by = !!rlang::sym(order_by), n = input) %>%
        dplyr::ungroup()

    }

  } else {

    out_df <- df[base::rep(FALSE, base::nrow(df))]

  }

  return(out_df)

}

#' @title Plota clockplot
#'
#' @description Visualize the evaluation of the fit of a numeric variable
#' against models around the area of an image annotation.
#'
#' @param fill Character value. The color with which the columns are filled.
#'
#' @inherit object_dummy params
#' @inherit variables_num params
#' @inherit imageAnnotationScreening params
#' @inherit ggplot2::facet_wrap params
#' @inherit ggplot2::facet_grid params
#' @inherit argument_dummy params
#'
#' @export
#'
setGeneric(name = "plotClockplot", def = function(object, ...){

  standardGeneric(f = "plotClockplot")

})

#' @rdname plotClockplot
#' @export
setMethod(
  f = "plotClockplot",
  signature = "spata2",
  definition = function(object,
                        id,
                        variables,
                        distance = NA_integer_,
                        binwidth = NA_integer_,
                        n_bins_circle = NA_integer_,
                        angle_span = c(0,360),
                        n_angle_bins = 12,
                        summarize_with = "mean",
                        model_subset = NULL,
                        model_remove = NULL,
                        model_add = NULL,
                        layout = 1,
                        switch = NULL,
                        fill = "steelblue",
                        ...){

    input_list <-
      check_ias_input(
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle
      )

    distance <- input_list$distance
    n_bins_circle <- input_list$n_bins_circle
    binwidth  <- input_list$binwidth

    temp_ias <-
      imageAnnotationScreening(
        object = object,
        id = id,
        variables = variables,
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        summarize_with = summarize_with,
        model_subset = model_subset,
        model_remove = model_remove,
        model_add = model_add
      )

    plotClockplot(
      object = temp_ias,
      layout = layout,
      switch = switch,
      fill = fill,
      ...
    )

  })



#' @rdname plotClockplot
#' @export
setMethod(
  f = "plotClockplot",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        variables,
                        model_subset = NULL,
                        model_remove = NULL,
                        layout = 1,
                        switch = NULL,
                        fill = "steelblue",
                        ...){

    ias_results_df <-
      dplyr::filter(object@results, variables %in% {{variables}}) %>%
      dplyr::mutate(bins_angle = base::factor(bins_angle, levels = make_angle_bins(object@n_bins_angle)))

    bins_angle <- base::levels(ias_results_df$bins_angle)
    models <- base::unique(ias_results_df$models)

    plot_df <-
      tidyr::expand_grid(
        variables = variables,
        models = models,
        bins_angle = base::factor(bins_angle, levels = bins_angle)
      ) %>%
      dplyr::left_join(y = ias_results_df, by = c("variables", "models", "bins_angle")) %>%
      dplyr::mutate(
        corr = tidyr::replace_na(corr, replace = 0),
      )

    if(base::is.character(model_subset)){

      plot_df <-
        dplyr::filter(plot_df, stringr::str_detect(pattern, pattern = model_subset))

    }

    if(base::is.character(model_remove)){

      plot_df <-
        dplyr::filter(plot_df, !stringr::str_detect(pattern, pattern = model_subset))
    }

    if(base::length(variables) == 1){

      facet_add_on <-
        ggplot2::facet_wrap(
          facets = . ~ models,
          nrow = nrow,
          ncol = ncol
          )

    } else if(layout == 1){

      facet_add_on <-
        ggplot2::facet_grid(
          rows = ggplot2::vars(variables),
          cols = ggplot2::vars(models),
          switch = switch
        )

    } else {

      facet_add_on <-
        ggplot2::facet_grid(
          rows = ggplot2::vars(models),
          cols = ggplot2::vars(variables),
          switch = switch
        )

    }

    plot_df$models <- make_pretty_model_names(plot_df$models)

    background_df <-
      dplyr::mutate(
        .data = plot_df,
        screened = !base::is.na(p_value),
        corr = dplyr::if_else(screened, true = 1, false = NaN)
      )

    ggplot2::ggplot(data = plot_df) +
      ggplot2::coord_polar() +
      ggplot2::theme_bw() +
      facet_add_on +
      ggplot2::geom_col(
        data = background_df,
        mapping = ggplot2::aes(x = bins_angle, y = corr),
        width = 1, color = "black", fill = "white"
      ) +
      ggplot2::geom_col(
        data = plot_df,
        mapping = ggplot2::aes(x = bins_angle, y = corr),
        width = 1, color = "black", fill = fill
      ) +
      ggplot2::scale_x_discrete(breaks = bins_angle, labels = bins_angle) +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::labs(x = NULL, y = NULL)

  }
)

#' @title Plot summary of spatial fitting
#'
#' @description Assigns every numeric variable to the model it fitted best
#' against and plots the p-value of the fit against the fit evaluation.
#'
#' @param x,y Character value. Specifies what to map to the x- and what to
#' map to the y-axis.
#'
#' @inherit argument_dummy params
#'
#' @export

setGeneric(name = "plotSummary", def = function(object, ...){

  standardGeneric(f = "plotSummary")

})

#' @rdname plotSummary
#' @export
setMethod(
  f = "plotSummary",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        x = "corr_mean",
                        y = "pvalue_mean",
                        pt_alpha = 0.75,
                        pt_color = "black",
                        pt_size = 1,
                        display_labels = FALSE,
                        model_subset = NULL,
                        model_remove = NULL,
                        pretty_names = FALSE,
                        ...){

    plot_df <-
      dplyr::group_by(object@results_smrd, variables) %>%
      dplyr::slice_max(order_by = !!rlang::sym(x), n = 1) %>%
      dplyr::ungroup()

    if(!base::is.null(model_subset)){

      keep <-
        confuns::vselect(
          input = base::unique(plot_df[["models"]]),
          dplyr::contains(match = model_subset)
        )

      plot_df <- dplyr::filter(plot_df, models %in% {{keep}})

    }

    if(!base::is.null(model_remove)){

      remove <-
        confuns::vselect(
          input = base::unique(plot_df[["models"]]),
          dplyr::contains(match = model_remove)
        )

      plot_df <- dplyr::filter(plot_df, !models %in% {{remove}})

    }

    if(base::isTRUE(pretty_names)){

      plot_df[["models"]] <- confuns::make_pretty_names(plot_df[["models"]])

    }

    ggplot2::ggplot(
      data = plot_df,
      mapping = ggplot2::aes(x = .data[[x]], y = -log10(.data[[y]]))
      ) +
      ggplot2::geom_point(
        alpha = pt_alpha,
        color = pt_color,
        size = pt_size
        ) +
      #ggplot2::scale_y_reverse(limits = c(1,0)) +
      ggplot2::facet_wrap(facets = . ~ models) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.grid = ggplot2::element_line(color = "lightgrey")
      )

  }
)



#' @title Visualize screening areaof IAS-algorithm
#'
#' @description Plots the surface of the sample three times with different
#' coloring to visualize how \code{imageAnnotationScreening()} screens
#' the sample depending on the input of arguments \code{binwidth}, \code{n_bins_circle},
#' \code{n_bins_angle}.
#'
#' @inherit getImageAnnotation params
#' @inherit imageAnnotationScreening params
#' @param color_core,color_outside Character value. Denotes
#' the colors with which the area of image annotation (\code{color_core})
#' and the area that is not included in the screening (\code{color_outside})
#' is displayed.
#' @param show_plots Logical value. If TRUE, the plots are immediately
#' plotted. If FALSE, only a list of plots is returned (invisibly).
#' @param display_angle,display_bins_angle,display_circle Logical value.
#' If TRUE, the plot is included. If FALSE, plotting is skipped.
#' @inherit argument_dummy params
#'
#' @return An invisible list of ggplots.
#'
#' @details The method for class \code{ImageAnnotationScreening} (the output of
#' the function \code{imageAnnotationScreening()}) can be used
#' to show the area on which the results base. Therefore, it does not have
#' arguments \code{binwidth}, \code{n_bins_circle} and \code{n_bins_angle}.
#'
#' @export

setGeneric(name = "plotSurfaceIAS", def = function(object, ...){

  standardGeneric(f = "plotSurfaceIAS")

})

#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "spata2",
  definition = function(object,
                        id,
                        distance = NA_integer_,
                        binwidth = NA_integer_,
                        n_bins_circle = NA_integer_,
                        angle_span = c(0,360),
                        n_bins_angle = 1,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = NULL,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        ...){

    hlpr_assign_arguments(object)

    if(base::length(pt_clrp) != 2){ pt_clrp <- base::rep(pt_clrp, 2)}

    ias_df <-
      getImageAnnotationScreeningDf(
        object = object,
        id = id,
        variables = NULL,
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        remove_circle_bins = "Core",
        rename_angle_bins = TRUE,
        drop = c(FALSE, TRUE),
        summarize_by = FALSE
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_circle",
            pt_clrp = pt_clrp[1],
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_angle",
            pt_clrp = pt_clrp[2],
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface2(
          coords_df = ias_df,
          color_by = "angle",
          pt_clrsp = pt_clrsp,
          pt_size = pt_size
        ) + ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

  }
)


#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = 2.25,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        ...){

    max_circles <- base::max(object@n_bins_circle)
    min_circles <- base::min(object@n_bins_circle)

    img_ann <- object@img_annotation
    img_ann_center <- getImageAnnotationCenter(img_ann)

    coords_df <- object@coords

    binwidth <- object@binwidth
    n_bins_angle <- object@n_bins_angle

    ias_df <-
      bin_by_area(
        coords_df = coords_df,
        area_df = img_ann@area,
        binwidth = binwidth,
        n_bins_circle = max_circles,
        remove = "Core"
      ) %>%
      bin_by_angle(
        center = img_ann_center,
        angle_span = object@angle_span,
        n_bins_angle = n_bins_angle,
        min_bins_circle = min_circles,
        rename = TRUE,
        remove = FALSE
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_circle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_angle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface2(
          coords_df = ias_df,
          color_by = "angle",
          pt_size = pt_size,
          pt_clrsp = pt_clrsp,
          pt_alpha = pt_alpha
        ) + ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

  }
)


#' @title Compare evaluatio of spatially opposing fits
#'
#' @description Plots a volcano plot by using the model evaluation
#' of spatial fitting as implemented by \code{imageAnnotationScreening()}
#' and \code{spatialTrajectoryScreening()}.
#'
#' @param eval Character value. The variable to use for the x-axis.
#' @param pval Character value. The variable to use for the y-axis.
#' @param left,right Character value. The name of the model whose best-fit variables
#' go to the left or to the right, respectively. Defaults to \code{left} = \emph{'linear_ascending'}
#' and \code{right} = \emph{'linear_descending'}.
#' @param display_threshold Logical value. If TRUE, the thresholds set by
#' \code{treshold_pval} and \code{threshold_eval} are used to color the points
#' of the plot.
#' @param threshold_pval,threshold_eval Numeric values that set the thresholds below/above
#' which the points are highlighted.
#' @param threshold_colors Character vector of length two. First denotes
#' the color of the significant variables, second denotes the color
#' of the not-significant variables.
#' @param label_vars Character value, numeric value or NULL. Useful to highlight
#' the exact position/evaulation of variables.
#'
#' If character, specifies the variables that are labeled. If numeric, specifies
#' the top n of variables that are labeled. If NULL, ignored.
#'
#' @param hstep,vstep Adjust the position of the two labels that show the
#' model names on the left and on the right.
#'
#' @param best_only Logical value. If TRUE, only variables are included in
#' the plot that have their best model fit in either the left or the right
#' model.
#'
#' @inherit argument_dummy params
#'
#'
#' @export

setGeneric(name = "plotVolcano", def = function(object, ...){

  standardGeneric(f = "plotVolcano")

})

#' @rdname plotVolcano
#' @export
setMethod(
  f = "plotVolcano",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        eval = "corr_mean",
                        pval = "pvalue_mean",
                        left = "linear_ascending",
                        right = "linear_descending",
                        display_thresholds = TRUE,
                        threshold_eval = 0.5,
                        threshold_pval = 0.05,
                        threshold_colors = c("tomato", "lightgrey"),
                        label_vars = NULL,
                        label_alpha = 0.9,
                        label_color = "black",
                        label_size = 2,
                        negative_log = TRUE,
                        pt_alpha = 0.9,
                        pt_size = 1,
                        display_names = TRUE,
                        hstep = 1.5,
                        vstep = 1.2,
                        best_only = FALSE,
                        ...){

    confuns::is_vec(x = threshold_colors, mode = "character", of.length = 2)

    ias_df_smrd <- object@results_smrd

    # if TRUE, the subsequent filtering will remove all variables that did not have
    # their best fit with the left or right model
    if(base::isTRUE(best_only)){

      ias_df_smrd <-
        dplyr::group_by(ias_df_smrd, variables) %>%
        dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
        dplyr::ungroup()

    }

    # subsequent filtering^^
    prel_plot_df <-
      dplyr::filter(
        .data = ias_df_smrd,
        stringr::str_detect(string = models, pattern = stringr::str_c(left, right, sep = "|"))
      )

    # if TRUE slice_max has already been applied above
    if(!base::isTRUE(best_only)){

      prel_plot_df <-
        dplyr::group_by(prel_plot_df, variables) %>%
        dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
        dplyr::ungroup()

    }

    prel_plot_df <-
      dplyr::mutate(
        .data = prel_plot_df,
        status = dplyr::case_when(
          !!rlang::sym(eval) >= {{threshold_eval}} & !!rlang::sym(pval) <= {{threshold_pval}} ~ "signif",
          TRUE ~ "not_signif"
        )
      )

    left_df <-
      dplyr::filter(prel_plot_df, stringr::str_detect(models, pattern = {{left}}))

    right_df <-
      dplyr::filter(prel_plot_df, stringr::str_detect(models, pattern = {{right}}))

    left_df[[eval]] <- left_df[[eval]] * -1

    plot_df <-
      base::rbind(left_df, right_df) %>%
      dplyr::ungroup()

    breaks_x <- base::seq(-1, 1, by = 0.2)

    labels_x <- stringr::str_remove(breaks_x, pattern = "^-")

    if(base::isTRUE(negative_log)){

      y_label <- stringr::str_c(pval, "(-log10)", sep = " ")

      plot_df[[pval]] <- -base::log10(x = plot_df[[pval]])

      threshold_pval <- -base::log10(threshold_pval)

    } else {

      y_label <- pval

    }

    if(!base::is.null(label_vars)){

      label_df <-
        pick_vars(
          df = dplyr::filter(plot_df, status == "signif"),
          input = label_vars,
          order_by = pval,
          neg_log = negative_log
        )

      label_add_on <-
        ggrepel::geom_text_repel(
          data = label_df,
          mapping = ggplot2::aes(x = .data[[eval]], y = .data[[pval]], label = variables),
          alpha = label_alpha,
          color = label_color,
          size = label_size,

          ...
        )

    } else {

      label_add_on <- NULL

    }

    max_y <- base::max(plot_df[[pval]])

    if(display_thresholds){

      tc <- threshold_eval
      tp <- threshold_pval

      hline_add_on <- ggplot2::geom_hline(yintercept = tp, linetype = "dashed", color = "grey")
      vline_add_on <- ggplot2::geom_vline(xintercept = c(-tc, tc), linetype = "dashed", color = "grey")

      mapping <- ggplot2::aes(x = .data[[eval]], y = .data[[pval]], color = .data[["status"]])

      color_add_on <-
        confuns::scale_color_add_on(
          variable = plot_df[["status"]],
          clrp = "milo",
          clrp.adjust = c("not_signif" = threshold_colors[2], "signif" = threshold_colors[1])
        )

      threshold_add_ons <-
        list(
          vline_add_on,
          hline_add_on,
          color_add_on
        )

    } else {

      mapping <- ggplot2::aes(x = .data[[eval]], y = .data[[pval]])
      threshold_add_ons <- NULL

    }

    if(base::is.character(display_names) | base::isTRUE(display_names)){

      if(base::is.character(display_names)){

        left <- display_names[1]
        right <- display_names[2]

      }

      annotation_df <-
        tibble::tibble(
          labels = confuns::make_pretty_names(c(left, right)),
          pos_x = c(-0.5, 0.5) * hstep,
          pos_y = max_y * vstep
        )

      text_add_on <-
        ggplot2::geom_text(
          data = annotation_df,
          mapping = ggplot2::aes(x = pos_x, y = pos_y, label = labels)
        )

    } else {

      text_add_on <- NULL

    }

    ggplot2::ggplot(data = plot_df) +
      threshold_add_ons +
      ggplot2::geom_point(
        data = plot_df,
        mapping = mapping,
        alpha = pt_alpha, size = pt_size) +
      label_add_on +
      text_add_on +
      #ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(
        limits = c(-1,1),
        breaks = breaks_x,
        labels = labels_x
      ) +
      ggplot2::labs(
        x = confuns::make_pretty_name(eval),
        y = confuns::make_pretty_name(y_label)
      ) +
      legendNone()

  }
)

#' @rdname plotVolcano
#' @export
setMethod(
  f = "plotVolcano",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        ...){


  }
)





# r -----------------------------------------------------------------------


# s -----------------------------------------------------------------------


#' @rdname export
subsetIAS <- function(ias, angle_span = NULL, angle_bins = NULL, variables = NULL, verbose = TRUE){

  if(purrr::map_lgl(c(angle_span, angle_bins, variables), .f = base::is.null)){

    stop("Please provide at least one subset input.")

  }

  stopifnot(base::min(angle_span) >= 0)
  stopifnot(base::max(angle_span) <= 360)

  amin_input <- base::min(angle_span)
  amax_input <- base::max(angle_span)

  n_bins <- ias@n_angle_bins

  confuns::give_feedback(
    msg = "Subsetting object of class ImageAnnotationScreening.",
    verbose = verbose
  )

  if(base::is.numeric(angle_span)){

    ias@results <-
      dplyr::mutate(
        .data = ias@results,
        temp = stringr::str_remove_all(base::as.character(bins_angle), pattern = "\\(|\\]")
      ) %>%
      tidyr::separate(col = temp, into = c("amin", "amax"), sep = ",") %>%
      dplyr::mutate(
        amin = base::as.numeric(amin),
        amax = base::as.numeric(amax)
      ) %>%
      dplyr::filter(amin >= {{amin_input}} & amax <= {{amax_input}}) %>%
      dplyr::select(-amin, -amax)

  }

  if(base::is.character(angle_bins)){

    ias@results <- dplyr::filter(ias@results, bins_angle %in% {{angle_bins}})

  }

  if(base::is.character(variables)){

    ias@results <- dplyr::filter(ias@results, variables %in% {{variables}})

  }

  ias@results$bins_angle <- base::droplevels(ias@results$bins_angle)

  confuns::give_feedback(
    msg = "Summarizing.",
    verbose = verbose
  )

  ias@results_smrd <- summarize_ias_df(df = ias@results)

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(ias)

}

#' @export
summarize_and_shift_variable_df <- function(grouped_df, variables){

  dplyr::summarise(
    .data = grouped_df,
    dplyr::across(
      .cols = dplyr::any_of(variables),
      .fns = ~ base::mean(.x, na.rm = TRUE)
    )
  ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::any_of(variables),
        .fns = confuns::normalize
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(variables),
      values_to = "values",
      names_to = "variables"
    ) %>%
    dplyr::mutate(
      bins_order = stringr::str_remove(bins_circle, pattern = "Circle ") %>% base::as.numeric()
    ) %>%
    # remove NA
    dplyr::group_by(variables) %>%
    dplyr::filter(!base::any(base::is.na(values)))


}


#' @title Summarize IAS-results
#'
#' @description Summarizes the results of the IAS-algorithm. Creates
#' the content of slot @@results_smrd of the \code{ImageAnnotationScreening}-class.
#'
#' @details Model fitting and evaluation happens within every angle-bin.
#' To get a single evaulation for every gene the results of every
#' angle-bin must be summarized.
#'
#' @export
summarizeIAS <- function(ias, method_padj = "fdr"){

  smrd_df <-
    dplyr::mutate(
      .data  = ias@results,
      p_value = tidyr::replace_na(data = p_value, replace = 1),
      corr = tidyr::replace_na(data = corr, replace = 0)
    ) %>%
    dplyr::group_by(variables, models) %>%
    dplyr::summarise(
      n_bins_angle = dplyr::n_distinct(bins_angle),
      corr_mean = base::mean(corr),
      corr_median = stats::median(corr),
      corr_min = base::min(corr),
      corr_max = base::max(corr),
      corr_sd = stats::sd(corr),
      raoc_mean = base::mean(raoc),
      pvalue_mean = base::mean(p_value),
      pvalue_median = stats::median(p_value),
      pvalue_combined = base::prod(p_value)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ias_score = (raoc_mean + corr_mean) / 2,
      pvalue_mean_adjusted = stats::p.adjust(p = pvalue_mean, method = method_padj),
      pvalue_median_adjusted = stats::p.adjust(p = pvalue_median, method = method_padj),
      pvalue_combined_adjusted = stats::p.adjust(p = pvalue_combined, method = method_padj)
    ) %>%
    dplyr::select(variables, models, ias_score, dplyr::everything())

  ias@method_padj <- method_padj

  ias@results_smrd <- smrd_df

  return(ias)

}
