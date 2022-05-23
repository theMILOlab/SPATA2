



# b -----------------------------------------------------------------------


bin_by_area <- function(coords_df, area_df, buffer, n_circles){

  circle_names <- stringr::str_c("Circle", 1:n_circles, sep = " ")

  circles <-
    purrr::set_names(
      x = c((1:n_circles)*buffer),
      nm = circle_names
    )

  buffer_vec <- c("Core" = 0, circles)

  areas <-
    purrr::imap(
      .x = buffer_vec,
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

  bin_levels <-
    c(base::names(buffer_vec)) %>%
    confuns::vselect(-dplyr::any_of(c("Core", "Outside")))

  out_df <-
    dplyr::filter(coords_df, !bins_circle %in% c("Core", "Outside")) %>%
    dplyr::mutate(
      bins_circle = base::factor(x = bins_circle, levels = bin_levels),
      bins_order = stringr::str_remove(bins_circle, pattern = "Circle") %>% base::as.numeric(),
      pt_in_plg = NULL
    )

  return(out_df)

}


bin_by_angle <- function(coords_df, center, n_bins = 12){

  base::stopifnot("bins_circle" %in% base::names(coords_df))

  mltp <- 360/n_bins
  breaks <- 0:n_bins * mltp

  out_df <-
    dplyr::group_by(.data = coords_df, barcodes) %>%
    dplyr::mutate(
      angle = compute_angle_between_two_points(p1 = c(x = x, y = y), p2 = center)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      bins_angle = base::cut(x = base::abs(angle), breaks = breaks)
    )

  return(out_df)

}


buffer_area <- function(df, buffer){

  area_grown <- Ternary::GrowPolygon(x = df$x, y = df$y, buffer = buffer)

  new_df <-
    base::data.frame(
      x = area_grown$x,
      y = area_grown$y
    ) %>%
    tibble::as_tibble()

  return(new_df)

}



# c -----------------------------------------------------------------------

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

  return(angle)

}

# e -----------------------------------------------------------------------

evaluate_spatial_patterns <- function(input_df){

  smrd_df <- input_df

  pb1 <-
    confuns::create_progress_bar(
      total = dplyr::n_distinct(smrd_df$variables),
      format = "Fitting models: [:bar] :percent eta: :eta"
    )

  pb2 <-
    confuns::create_progress_bar(
      total = dplyr::n_distinct(smrd_df$variables),
      format = "Evaluating AUC: [:bar] :percent eta: :eta"
    )

  pb3 <-
    confuns::create_progress_bar(
      total = dplyr::n_distinct(smrd_df$variables),
      format = "Evaluating Correlation: [:bar] :percent eta: :eta"
    )

  pb4 <-
    confuns::create_progress_bar(
      total = dplyr::n_distinct(smrd_df$variables),
      format = "Joining results: [:bar] :percent eta: :eta"
    )

  pdf <- create_model_df(input = dplyr::n_distinct(smrd_df$bins_circle))

  all_models <- base::names(pdf)

  nested_df <-
    dplyr::group_by(smrd_df, variables) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      residuals = purrr::map(
        .x = data,
        .f = ~
          hlpr_add_residuals2(
            df = base::cbind(.x, pdf) %>% dplyr::select_if(.predicate = base::is.numeric),
            pb = pb1,
            column_order = "bins_order"
          )
      )
    ) %>%
    dplyr::mutate(
      auc = purrr::map(
        .x = residuals,
        .f = ~
          hlpr_summarize_residuals(
            df = .x,
            pb = pb2,
            column_order = "bins_order",
            shift_longer = TRUE
          )
      )
    ) %>%
    dplyr::mutate(
      corr = purrr::map(
        .x = data,
        pb = pb3,
        .f = function(df, pb = NULL){

          if(!base::is.null(pb)){ pb$tick() }

          base::cbind(df, pdf) %>%
            tibble::as_tibble() %>%
            tidyr::pivot_longer(
              cols = dplyr::all_of(all_models),
              names_to = "pattern_names",
              values_to = "pattern_values"
            ) %>%
            dplyr::group_by(pattern_names) %>%
            tidyr::nest() %>%
            dplyr::summarize(
              corr_res = purrr::map(
                .x = data,
                .f = ~
                  stats::cor.test(x = .x$pattern_values, y = .x$values) %>%
                  broom::tidy() %>%
                  dplyr::select(corr_value = estimate, p_value = p.value)
              )
            ) %>%
            tidyr::unnest(cols = corr_res)

        }
      )
    ) %>%
    dplyr::select(-data, -residuals) %>%
    dplyr::mutate(
      eval = purrr::map2(
        .x = auc,
        .y = corr,
        pb = pb4,
        .f = function(xdf, ydf, pb){

          if(!base::is.null(pb)){ pb$tick() }

          dplyr::left_join(x = xdf, y = ydf, by = "pattern_names")

        })
    ) %>%
    dplyr::select(-auc, -corr) %>%
    tidyr::unnest(eval)

  out_df <-
    dplyr::mutate(
      nested_df,
      pattern = hlpr_name_models(pattern_names),
      auc_residuals = auc,
      auc_residuals_scaled = auc / n_circles
    ) %>%
    dplyr::group_by(variables) %>%
    dplyr::ungroup()

  return(out_df)

}


# g -----------------------------------------------------------------------


getPatternEvaluationDf <- function(object,
                                   id,
                                   variables,
                                   n_circles,
                                   buffer){

  img_ann <- getImageAnnotation(object = object, id = id, add_image = FALSE)

  circle_names <- stringr::str_c("Circle", 1:n_circles, sep = " ")

  circles <-
    purrr::set_names(
      x = c((1:n_circles)*buffer),
      nm = circle_names
    )

  buffer_vec <- c("Core" = 0, circles)

  areas <-
    purrr::imap(
      .x = buffer_vec,
      .f = ~ buffer_area(df = img_ann@area, buffer = .x)
    )

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::select(barcodes, x, y)

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
          bins_circle == "Outside" & pt_in_plg %in% c(1,2) ~ {{area}},
          TRUE ~ bins_circle
        )
      )

  }

  center <- getImageAnnotationCenter(object, id = id)

  ymax <- base::max(getImageRange(object)[["y"]])

  bin_levels <-
    c(base::names(buffer_vec)) %>%
    confuns::vselect(-dplyr::any_of(c("Core", "Outside")))

  n_angle_bins <- 12
  breaks <- 0:12 * 30

  binned_df <-
    dplyr::filter(coords_df, !bins_circle  %in% c("Core", "Outside")) %>%
    dplyr::mutate(
      bins_circle = base::factor(bins_circle, levels = bin_levels),
      pt_in_plg = NULL
    ) %>%
    dplyr::group_by(barcodes) %>%
    dplyr::mutate(
      angle = compute_angle_between_two_points(p1 = center, p2 = c(x = x, y = y))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      bins_angle = base::cut(
        x = base::abs(angle),
        breaks = breaks
      )
    )

  bins_angle_levels <- base::levels(binned_df$bins_angle)

  keep_bins_angle <-
    dplyr::distinct(binned_df, bins_circle, bins_angle) %>%
    dplyr::group_by(bins_angle) %>%
    dplyr::tally() %>%
    dplyr::filter(n == {{n_circles}}) %>%
    dplyr::pull(bins_angle) %>%
    base::as.character()

  binned_df <- dplyr::filter(binned_df, bins_angle %in% {{keep_bins_angle}})

  variable_df <-
    joinWithVariables(
      object = object,
      spata_df = binned_df,
      variables = variables,
      smooth = FALSE
    )


  bins_angle_to_screen <- base::as.character(base::unique(variable_df$bins_angle))

  n_bins <- base::length(bins_angle_to_screen)

  eval_df <-
    purrr::map(
      .x = bins_angle_to_screen,
      .f = purrr::safely(.f = function(bin){

        ref_bin <- base::which(bins_angle_to_screen == bin)

        confuns::give_feedback(msg = glue::glue("Working on angle bin {ref_bin}/{n_bins}."))

        grouped_df <-
          dplyr::filter(variable_df, bins_angle == {{bin}}) %>%
          dplyr::group_by(bins_circle, bins_angle)

        smrd_df <-
          summarize_and_shift_variable_df(
            grouped_df = grouped_df,
            variables = variables
          )

        eval_df <- evaluate_spatial_patterns(input_df = smrd_df)

        return(eval_df)

      }, otherwise = NULL)
    ) %>%
    purrr::set_names(nm = bins_angle_to_screen) %>%
    purrr::keep(.p = ~ base::is.null(.x$error)) %>%
    purrr::imap_dfr(
      .f = ~ dplyr::mutate(.x$result, bins_angle = base::factor(x = .y, levels = bins_angle_levels))
    )

  return(eval_df)

}


getResultsRAPI <- function(object, id){

  out <- object@spatial[[1]][["rapi"]][[id]]

  check_availability(
    test = !base::is.null(out),
    ref_x = "region associated pattern identifcation",
    ref_fns = "runRAPI()"
  )

  return(out)

}




# p -----------------------------------------------------------------------


plotEvaluationRAPI <- function(object,
                               id,
                               variables,
                               n_circles,
                               buffer,
                               pattern = NULL,
                               layout = 1,
                               switch = NULL,
                               ...){

  eval_df <-
    getPatternEvaluationDf(
      object = object,
      id = id,
      variables = variables,
      n_circles = n_circles,
      buffer = buffer
    ) %>% dplyr::mutate(pattern = pattern_names)

  bins_angle <- base::levels(eval_df$bins_angle)

  pattern <- base::unique(eval_df$pattern)

  plot_df <-
    tidyr::expand_grid(
      variables = variables,
      pattern = pattern,
      bins_angle = base::factor(bins_angle, levels = bins_angle)
    ) %>%
    dplyr::left_join(y = eval_df, by = c("variables", "pattern", "bins_angle")) %>%
    dplyr::mutate(
      corr_value = tidyr::replace_na(corr_value, replace = 0)
    )


  if(base::length(variables) == 1){

    facet_add_on <- ggplot2::facet_wrap(facets = . ~ pattern, ...)

  } else if(layout == 1){

    facet_add_on <-
      ggplot2::facet_grid(
        rows = ggplot2::vars(variables),
        cols = ggplot2::vars(pattern),
        switch = switch
      )

  } else {

    facet_add_on <-
      ggplot2::facet_grid(
        rows = ggplot2::vars(pattern),
        cols = ggplot2::vars(variables),
        switch = switch
      )

  }

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_col(
      mapping = ggplot2::aes(x = bins_angle, y = corr_value),
      width = 1, color = "black", fill = "steelblue",
    ) +
    facet_add_on +
    ggplot2::scale_y_continuous(limits = c(0,1)) +
    ggplot2::theme_bw() +
    ggplot2::coord_polar() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL)

}


plotSummaryRAPI <- function(object,
                            id,
                            x = "corr_median",
                            y = "pval_median_adjusted",
                            pt_alpha = 0.9,
                            pt_color = "black",
                            pt_size = 1,
                            display_labels = FALSE,
                            n_labels = 10,
                            var_labels = "x",
                            pattern_subset = NULL,
                            ...){

  rapi_df <-
    getResultsRAPI(object, id)$eval_smrd %>%
    dplyr::mutate(
      pattern_names = hlpr_name_models(pattern_names) %>% stringr::str_remove(pattern = "^p_")
    ) %>%
    confuns::check_across_subset(
      df = .,
      across = "pattern_names",
      across.subset = pattern_subset,
      relevel = TRUE
    )

  if(base::isTRUE(display_labels)){

    label_df <- rapi_df

    if(base::length(n_labels) == 1){

      label_df <- dplyr::group_by(label_df, pattern_names)

      if(var_labels == "x"){

        label_df <- dplyr::slice_max(label_df, order_by = !!rlang::sym(x), n = n_labels)

      } else {

        label_df <- dplyr::slice_min(label_df, order_by = !!rlang::sym(y), n = n_labels)

      }

    }

    label_add_on <-
      ggrepel::geom_label_repel(
        data = label_df,
        mapping = ggplot2::aes(label = variables),
        ...
      )

  } else {

    label_add_on <- NULL

  }

  ggplot2::ggplot(data = rapi_df, mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(alpha = pt_alpha, color = pt_color, size = pt_size) +
    ggplot2::scale_y_reverse(limits = c(1,0)) +
    ggplot2::facet_wrap(facets = . ~ pattern_names) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = ggplot2::element_line(color = "lightgrey")
    ) +
    label_add_on

}

plotSurfaceGradientBins <- function(object,
                                    id,
                                    n_circles,
                                    buffer,
                                    pt_alpha = NA_integer_,
                                    pt_clrp = "inferno",
                                    pt_size = NULL,
                                    color_core = ggplot2::alpha("grey", 0),
                                    color_outside = ggplot2::alpha("lightgrey", 0.25),
                                    direction = -1,
                                    display_image = FALSE,
                                    ...){

  hlpr_assign_arguments(object)

  img_ann <- getImageAnnotation(object = object, id = id, add_image = FALSE)

  circle_names <- stringr::str_c("Circle", 1:n_circles, sep = " ")

  circles <-
    purrr::set_names(
      x = c((1:n_circles)*buffer),
      nm = circle_names
    )

  buffer_vec <- c("Core" = 0, circles)

  areas <-
    purrr::imap(
      .x = buffer_vec,
      .f =
        ~ buffer_area(df = img_ann@area, buffer = .x) %>%
        dplyr::mutate(., circle = .y)
    )

  coords_df <- getCoordsDf(object)

  coords_df$bins <- "Outside"

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
        bins = dplyr::case_when(
          bins == "Outside" & pt_in_plg %in% c(1,2) ~ {{area}},
          TRUE ~ bins
        )
      )

  }

  bin_levels <- c(base::names(buffer_vec), "Outside")

  clrp_adjust <-
    c(
      color_core, viridis::inferno(n = length(buffer_vec)-1, direction = -1), color_outside
    ) %>% set_names(bin_levels)

  clrp_adjust <-
    c(
      color_core,
      viridis::viridis_pal(option = pt_clrp)(n = base::length(buffer_vec)-1),
      color_outside
    ) %>%
    purrr::set_names(nm = bin_levels)

  if(!base::is.na(pt_alpha)){

    clrp_adjust[2:(base::length(clrp_adjust)-1)] <-
      ggplot2::alpha(clrp_adjust[2:(base::length(clrp_adjust)-1)], alpha = pt_alpha)

  }

  if(base::isTRUE(display_image)){

    img <- getImage(object)

  } else {

    img <- NULL

  }

  plotSurface2(
    coords_df = coords_df,
    color_by = "bins",
    pt_alpha = NA_integer_,
    pt_clrp = "milo",
    pt_size = pt_size,
    clrp_adjust = clrp_adjust,
    image = img
  ) +
    ggplot2::labs(color = "Bins")

}



# r -----------------------------------------------------------------------


runRAPI <- function(object,
                    id,
                    variables,
                    buffer,
                    n_circles){

  out_eval_df <-
    getPatternEvaluationDf(
      object = object,
      id = id,
      variables = variables,
      n_circles = n_circles,
      buffer = buffer
    )

  out_smrd_eval_df <-
    dplyr::group_by(out_eval_df, variables, pattern_names) %>%
    dplyr::summarise(
      n_bins_angle = dplyr::n_distinct(bins_angle),
      corr_mean = base::mean(corr_value),
      corr_median = stats::median(corr_value),
      corr_min = base::min(corr_value),
      corr_max = base::max(corr_value),
      corr_sd = stats::sd(corr_value),
      pvalue_mean = base::mean(p_value),
      pvalue_median = stats::median(p_value),
      pvalue_combined = base::prod(p_value),
      auc_mean = base::mean(auc),
      auc_median = stats::median(auc),
      auc_min = base::min(auc),
      auc_max = base::max(auc)
    ) %>%
    dplyr::mutate(
      pval_mean_adjusted = stats::p.adjust(p = pvalue_mean, method = "fdr"),
      pval_median_adjusted = stats::p.adjust(p = pvalue_median, method = "fdr"),
      pval_combined_adjusted = stats::p.adjust(p = pvalue_combined, method = "fdr")
    ) %>%
    dplyr::group_by(variables) %>%
    dplyr::slice_max(order_by = corr_median, n = 1) %>%
    dplyr::arrange(pval_median_adjusted)

  res <-
    list(
      bins_angle = base::sort(bins_angle_to_screen),
      buffer = buffer,
      id = id,
      n_circles = n_circles,
      eval = out_eval_df,
      eval_smrd = out_smrd_eval_df,
      variables = variables
    )

  object@spatial[[1]][["rapi"]][[img_ann@id]] <- res

  return(object)

}


# s -----------------------------------------------------------------------


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




