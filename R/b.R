

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
#' @inherit bin_by_expansion params
#'
#' @export
bin_by_angle <- function(coords_df,
                         center,
                         n_bins_angle = 1,
                         angle_span = c(0,360),
                         min_bins_circle = NULL,
                         rename = FALSE,
                         remove = FALSE,
                         drop = TRUE,
                         var_to_bin = "barcodes",
                         verbose = TRUE){

  confuns::is_vec(x = angle_span, mode = "numeric", of.length = 2)

  angle_span <- c(from = angle_span[1], to = angle_span[2])

  range_span <- base::range(angle_span)

  if(angle_span[1] == angle_span[2]){

    stop("Invalid input for argument `angle_span`. Must contain to different values.")

  } else if(base::min(angle_span) < 0 | base::max(angle_span) > 360){

    stop("Input for argument `angle_span` must range from 0 to 360.")

  }

  from <- angle_span[1]
  to <- angle_span[2]

  confuns::give_feedback(
    msg = glue::glue("Including area between {from}° and {to}°."),
    verbose = verbose
  )

  if(n_bins_angle > 1){

    confuns::give_feedback(
      msg = glue::glue("Binning included area in {n_bins_angle} angle-bins."),
      verbose = verbose
    )

  }

  if(confuns::is_list(center)) {

    base::stopifnot("border" %in% base::colnames(coords_df))

    prel_angle_df <-
      purrr::imap_dfr(
        .x = center,
        .f = function(c, b){

          dplyr::filter(coords_df, border == {{b}}) %>%
            dplyr::group_by(!!rlang::sym(var_to_bin)) %>%
            dplyr::mutate(
              angle = compute_angle_between_two_points(
                p1 = c(x = x, y = y),
                p2 = {{c}}
              )
            )

        }
      )

  } else {

    # compute angle
    prel_angle_df <-
      dplyr::group_by(.data = coords_df, !!rlang::sym(var_to_bin)) %>%
      dplyr::mutate(
        angle = compute_angle_between_two_points(
          p1 = c(x = x, y = y),
          p2 = center
        )
      )

  }



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
      y = prel_angle_df[,c(var_to_bin, "angle")],
      by = var_to_bin
    ) %>%
    dplyr::left_join(
      x = .,
      y = prel_angle_bin_df[,c(var_to_bin, "bins_angle")],
      by = var_to_bin
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
    dplyr::select(angle_df, dplyr::any_of(c("bins_circle")), bins_angle) %>%
    dplyr::distinct() %>%
    dplyr::filter(bins_angle != "Outside") %>%
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
#' @export
#'
bin_by_expansion <- function(coords_df,
                             area_df,
                             binwidth,
                             n_bins_circle,
                             remove = FALSE,
                             bcsp_exclude = NULL,
                             drop = TRUE,
                             arrange = TRUE){

  n_bins_circle <- base::max(n_bins_circle)

  circle_names <- stringr::str_c("Circle", 1:n_bins_circle, sep = " ")

  circles <-
    purrr::set_names(
      x = c((1:n_bins_circle)*binwidth),
      nm = circle_names
    )

  binwidth_vec <- c("Core" = 0, circles)

  # create new variable. Default is 'Outside'.
  # values will be overwritten with every additional loop
  coords_df$bins_circle <- "Outside"
  coords_df$border <- NA_character_

  # if outer circles exist
  if("outer" %in% area_df[["border"]]){

    outer_df <- dplyr::filter(area_df, border == "outer")

    expansions_pos <-
      purrr::imap(
        .x = binwidth_vec,
        .f = ~ buffer_area(df = outer_df[c("x", "y")], buffer = .x)
      )

    for(exp in base::names(expansions_pos)){

      exp_df <- expansions_pos[[exp]]

      coords_df$pt_in_plg <-
        sp::point.in.polygon(
          point.x = coords_df$x,
          point.y = coords_df$y,
          pol.x = exp_df$x,
          pol.y = exp_df$y
        )

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          bins_circle = dplyr::case_when(
            # if bins_circle is NOT 'Outside' it has already bin binned
            bins_circle == "Outside" & pt_in_plg %in% c(1,2) ~ {{exp}},
            TRUE ~ bins_circle
          ),
          border = dplyr::case_when(
            pt_in_plg %in% c(1,2) ~ "outer",
            TRUE ~ border
          )
        )
    }

  }


  # if inner circles exist
  if("inner1" %in% area_df[["border"]]){

    holes <-
      stringr::str_subset(area_df[["border"]], pattern = "inner") %>%
      base::unique()

    # correct binwidth vec for screening towards the inside
    binwidth_vec_neg <-
      purrr::set_names(
        x = (-binwidth_vec[1:(n_bins_circle-1)]),
        nm = base::names(binwidth_vec)[2:n_bins_circle]
      )

    for(h in base::seq_along(holes)){

      hole <- holes[h]

      hole_df <- dplyr::filter(area_df, border == {{hole}})

      expansions_neg <-
        purrr::imap(
          .x = binwidth_vec_neg,
          .f = ~ buffer_area(df = hole_df[c("x", "y")], buffer = .x)
        ) %>%
        purrr::discard(.p = base::is.null) %>%
        purrr::map(.f = tibble::as_tibble)

      for(exp in base::rev(base::names(expansions_neg))){

        exp_df <- expansions_neg[[exp]]

        coords_df$pt_in_plg <-
          sp::point.in.polygon(
            point.x = coords_df$x,
            point.y = coords_df$y,
            pol.x = exp_df$x,
            pol.y = exp_df$y
          )

        coords_df <-
          dplyr::mutate(
            .data = coords_df,
            bins_circle = dplyr::case_when(
              # if bins_circle is NOT 'Outside' it has already bin binned
              # if bins_circle is 'Core' it can be overwritten as 'Core'
              # means everything inside the outer border
              bins_circle %in% c("Outside", "Core") & pt_in_plg %in% c(1,2) ~ {{exp}},
              TRUE ~ bins_circle
            ),
            border = dplyr::case_when(
              pt_in_plg %in% c(1,2) ~ {{hole}},
              TRUE ~ border
            )
          )

      }

    }

  }

  # Option bcsp_exclude: Rename manually selected spots to "Outside"
  if(base::is.character(bcsp_exclude)){
    if(any(!bcsp_exclude %in% coords_df$barcodes)){
      warning("Barcode(s) given in `bcsp_exclude` not found in spata object. Is the format correct?")
    }
    coords_df[coords_df$barcodes %in% bcsp_exclude,]$bins_circle <- "Outside"
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

  if(base::isTRUE(arrange)){

    out_df <- dplyr::arrange(out_df, bins_order)

  }

  return(out_df)

}


#' @export
bin_projection_df <- function(projection_df, n_bins = NULL, binwidth = NULL){

  # prioritize n_bins cause binwidth is defined by default with getCCD()
  # if n_bins is valid numeric input it has been set manually
  if(base::is.numeric(n_bins) & !base::is.na(n_bins)){

    binned_projection_df <-
      dplyr::mutate(
        .data = projection_df,
        proj_length_binned = base::cut(projection_length, breaks = n_bins),
        order_numeric = base::as.numeric(proj_length_binned)
      )

  } else {

    is_dist_pixel(input = binwidth, error = TRUE)

    binned_projection_df <-
      dplyr::mutate(
        .data = projection_df,
        proj_length_binned = plyr::round_any(x = projection_length, accuracy = {{binwidth}}, f = base::ceiling),
        order_numeric = base::as.factor(proj_length_binned) %>% base::as.numeric()
      )

  }

  return(binned_projection_df)

}


br_add <- function(height, break_add = NULL){

  if(base::is.null(break_add)){

    x <- base::ceiling((height - 400)/100)

    out <- x*5

  } else {

    out <- break_add

  }

  return(out)

}

breaks <- function(n){

  base::rep("<br>", n) %>%
    stringr::str_c(collapse = "") %>%
    shiny::HTML()

}


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
#' \code{dist} of function \code{sf::st_buffer()}.
#'
#' @export
buffer_area <- function(df, buffer, close_plg = TRUE){

  frow <- df[1, c("x", "y")] %>% base::as.numeric()
  lrow <- df[base::nrow(df), c("x", "y")] %>% base::as.numeric()

  if(!base::identical(frow, lrow) & base::isTRUE(close_plg)){

    df <- close_area_df(df)

  }

  area_grown <-
    sf::st_polygon(x = list(base::as.matrix(df[,c("x", "y")]))) %>%
    sf::st_buffer(dist = buffer) %>%
    base::as.matrix() %>%
    base::as.data.frame()

  if(purrr::is_empty(area_grown)){

    out <- NULL

  } else {

    out <-
      magrittr::set_colnames(area_grown, value = c("x", "y")) %>%
      tibble::as_tibble()

  }

  return(out)

}

