continuous_ggplot2_aes <- c("x", "y", "alpha", "color", "fill", "shape", "size")

transformation_fns <-
  list(
    "exp" = base::exp,
    "log" = base::log,
    "log2" = base::log2,
    "log10" = base::log10,
    "sqrt" = base::sqrt
  )


#' @export
transform_df <- function(df, transform.with, sep = "_"){

  if(is_list(transform.with) & !purrr::is_empty(transform.with)){

    names_tw <- base::names(transform.with)

    check_one_of(
      input = names_tw,
      against = dplyr::select_if(df, base::is.numeric) %>% names(),
      ref.opt.2 = "numeric variable names",
      fdb.opt = 2
    )

    for(name_tw in names_tw){

      var_to_transform <- name_tw

      transform_info <- transform.with[[name_tw]]

      if(base::is.character(transform_info)){

        for(tf in transform_info){

          df <-
            dplyr::mutate(
              .data = df,
              dplyr::across(
                .cols = {{var_to_transform}},
                .fns = transformation_fns[[tf]],
                .names = NULL
              )
            )

        }

      } else if(is_list(transform_info)){

        for(tf in transform_info){

          if(base::is.function(tf)){

            df <-
              dplyr::mutate(
                .data = df,
                dplyr::across(
                  .cols = {{var_to_transform}},
                  .fns = tf,
                  .names = NULL
                )
              )

          } else if(base::is.character(tf)){

            df <-
              dplyr::mutate(
                .data = df,
                dplyr::across(
                  .cols = {{var_to_transform}},
                  .fns = transformation_fns[[tf]],
                  .names = NULL
                )
              )

          }

        }

      }

    }

  }

  return(df)

}
