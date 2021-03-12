
# get-functions -----------------------------------------------------------


#' @title Obtain information about the optimal neural network set up
#'
#' @description Extracts the results from \code{assessAutoencoderOptions()}.
#'
#' @inherit check_object params
#'
#' @return A data.frame containing the total variance measured by \code{irlba::prcomp_irlba()} after each
#' combination of activations/bottlenecks.
#' @export

getAutoencoderAssessment <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  assessment <- object@autoencoder[[of_sample]]$assessment

  test <- !(base::identical(assessment, list()) & base::is.null(assessment))
  ref_x <- "autoencoder assessment information"
  ref_fns <- "function runAutoencoderAssessment() first"

  check_availability(
    test = test,
    ref_x = ref_x,
    ref_fns = ref_fns
  )

  base::return(assessment)

}


#' @title Obtain information on neural network
#'
#' @description Returns the argument input that was chosen to construct the
#' neural network that generated the matrix denoted in \code{mtr_name}.
#'
#' @inherit getExpressionMatrix params
#'
#' @return A named list.
#' @export

getAutoencoderSetUp <- function(object, mtr_name, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  nn_set_up <-
    object@autoencoder[[of_sample]][["nn_set_ups"]][[mtr_name]]

  if(base::is.null(nn_set_up)){

    base::stop(glue::glue("Could not find any autoencoder information for matrix '{mtr_name}' of sample '{of_sample}'"))

  }

  base::return(nn_set_up)

}

# -----



# main-functions ----------------------------------------------------------

#' @title Assessment of Neural Network Set Up
#'
#' @description Assesses different neural network set ups regarding
#' the activation function and the number of bottleneck neurons.
#'
#' @inherit check_object params
#' @param expr_mtr The expression matrix that is to be used as input for the neural network.
#' @param activations Character vector. Denotes the activation functions to be assessed.
#' @param bottlenecks Numeric vector. Denotes the different numbers of bottleneck neurons to be assessed.
#' @inherit runAutoencoderDenoising params
#'
#' @return
#'
#' \itemize{
#' \item{\code{runAutoencoderAssessment()}: The spata object containing the list that holds the total variance measured by \code{irlba::prcomp_irlba()} after each
#' combination of activations/bottlenecks as well as the additional set up.}
#' \item{\code{assessAutoencoderOptions()}:
#' The list that holds the total variance measured by \code{irlba::prcomp_irlba()} after each combination
#' of activations/bottlenecks as well as the additional set up.}
#' }
#'
#' @export

runAutoencoderAssessment <- function(object,
                                     activations,
                                     bottlenecks,
                                     layers = c(128, 64, 32),
                                     dropout = 0.1,
                                     epochs = 20,
                                     verbose = TRUE){

  check_object(object)

  assessment_list <-
    assessAutoencoderOptions(expr_mtr = getExpressionMatrix(object, of_sample = "", mtr_name = "scaled"),
                              activations = activations,
                              bottlenecks = bottlenecks,
                              layers = layers,
                              dropout = dropout,
                              epochs = epochs,
                              verbose = verbose)

  object <- setAutoencoderAssessment(object = object, assessment_list = assessment_list)

  base::return(object)

}

#' @rdname runAutoencoderAssessment
#' @export
assessAutoencoderOptions <- function(expr_mtr,
                                     activations,
                                     bottlenecks,
                                     layers = c(128, 64, 32),
                                     dropout = 0.1,
                                     epochs = 20,
                                     verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  confuns::check_one_of(input = activations, against = activation_fns)

  confuns::are_values(c("dropout", "epochs"), mode = "numeric")

  confuns::is_vec(x = layers, mode = "numeric", of.length = 3)
  confuns::is_vec(x = bottlenecks, mode = "numeric")

  # 2. Assess all combinations in for loop ----------------------------------

  activations_list <-
    base::vector(mode = "list", length = base::length(activations)) %>%
    purrr::set_names(nm = activations)

  for(a in base::seq_along(activations)){

    activation <- activations[a]

    bottlenecks_list <-
      base::vector(mode = "list", length = base::length(bottlenecks)) %>%
      purrr::set_names(nm = stringr::str_c("bn", bottlenecks, sep = "_"))

    for(b in base::seq_along(bottlenecks)){

      bottleneck <- bottlenecks[b]

      base::message(Sys.time())
      base::message(glue::glue("Assessing activation option {a}/{base::length(activations)}:'{activation}' and bottleneck option {b}/{base::length(bottlenecks)}: {bottleneck}"))

      # Neural network ----------------------------------------------------------

      input_layer <-
        keras::layer_input(shape = c(base::ncol(expr_mtr)))

      encoder <-
        input_layer %>%
        keras::layer_dense(units = layers[1], activation = activation) %>%
        keras::layer_batch_normalization() %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[2], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[3], activation = activation) %>%
        keras::layer_dense(units = bottleneck)

      decoder <-
        encoder %>%
        keras::layer_dense(units = layers[3], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[2], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[1], activation = activation) %>%
        keras::layer_dense(units = c(ncol(expr_mtr)))

      autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

      autoencoder_model %>% keras::compile(
        loss = 'mean_squared_error',
        optimizer = 'adam',
        metrics = c('accuracy')
      )

      history <-
        autoencoder_model %>%
        keras::fit(expr_mtr, expr_mtr, epochs = epochs, shuffle = TRUE,
                   validation_data = list(expr_mtr, expr_mtr), verbose = verbose)

      reconstructed_points <-
        autoencoder_model %>%
        keras::predict_on_batch(x = expr_mtr)

      base::rownames(reconstructed_points) <- base::rownames(expr_mtr)
      base::colnames(reconstructed_points) <- base::colnames(expr_mtr)


      # PCA afterwards ----------------------------------------------------------

      bottlenecks_list[[b]] <- irlba::prcomp_irlba(base::t(reconstructed_points), n = 30)

    }

    activations_list[[a]] <- bottlenecks_list

  }

  # 3. Summarize in data.frame ----------------------------------------------

  res_df <-
    purrr::imap_dfr(.x = activations_list, .f = function(.list, .name){

      data.frame(
        activation = .name,
        bottleneck = stringr::str_remove(string = base::names(.list), pattern = "^bn_"),
        total_var = purrr::map_dbl(.x = .list, .f = "totalvar")
      )

    }) %>% tibble::remove_rownames()

  res_df$bottleneck <- base::factor(res_df$bottleneck, levels = base::unique(res_df$bottleneck))

  pca_scaled <- irlba::prcomp_irlba(x = base::t(expr_mtr), n = 30)

  assessment_list <- list("df" = res_df,
                          "set_up" = list("epochs" = epochs, "dropout" = dropout, "layers" = layers),
                          "scaled_var" = pca_scaled$totalvar)

  base::return(assessment_list)

}


#' @title Denoise expression matrix
#'
#' @description This function constructs and uses a neural network to denoise
#' expression levels spatially.
#'
#' @inherit check_object params
#' @param layers Numeric vector of length 3. Denotes the number of neurons in the three hidden layers.
#'  (default = c(128, 64, 32))
#' @param bottleneck Numeric value. Denotes the number of bottleneck neurons.
#' @param mtr_name_output Character value. Denotes the name under which the denoised matrix is stored
#' in the data slot.
#' @param dropout Numeric value. Denotes the dropout. (defaults to 0.1)
#' @param activation Character value. Denotes the activation function. (defaults to \emph{'relu'})
#' @param epochs Numeric value. Denotes the epochs of the neural network. (defaults to 20)
#' @param display_plot Logical. If set to TRUE a scatter plot of the result is displayed in the viewer pane.
#' See documentation for \code{plotAutoencoderResults()} for more information.
#' @param genes Character vector of length two. Denotes the genes to be used for the validation plot.
#' @param set_as_active Logical. If set to TRUE the denoised matrix is set as the active matrix via
#' \code{setActiveExpressionMatrix()}.
#'
#' @return A spata-object containing the denoised expression matrix in slot @@data$denoised. This matrix
#' is then denoted as the active matrix.
#'
#' @importFrom Seurat ScaleData
#'
#' @export

runAutoencoderDenoising <- function(object,
                                    activation,
                                    bottleneck,
                                    mtr_name_output = "denoised",
                                    layers = c(128, 64, 32),
                                    dropout = 0.1,
                                    epochs = 20,
                                    display_plot = FALSE,
                                    genes,
                                    set_as_active = FALSE,
                                    verbose = TRUE,
                                    of_sample = NA){

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(
    msg = base::message("Checking input for validity."),
    verbose = verbose
  )

  check_object(object)

  confuns::are_values(c("dropout", "epochs"), mode = "numeric")
  confuns::are_vectors(c("activation"), mode = "character")
  confuns::are_values(c("display_plot", "set_as_active", "verbose"), mode = "logical")

  confuns::is_vec(x = layers, mode = "numeric", of.length = 3)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::isTRUE(display_plot)){

    # check validation genes
    val_genes <- check_genes(object, genes = genes, max_length = 2, of.length = 2)

    base::stopifnot(base::length(val_genes) == 2)

  }

  x_train <- getExpressionMatrix(object, mtr_name = "scaled" , of_sample = of_sample)

  # assess optimum if any of the two inputs are vectors
  if(base::any(purrr::map_int(.x = list(bottleneck, activation), .f = base::length) > 1)){

    assessment_list <-
      assessAutoencoderOptions2(expr_mtr = x_train,
                                dropout = dropout,
                                epochs = epochs,
                                layers = layers,
                                bottlenecks = bottleneck,
                                activations = activation,
                                verbose = FALSE)

    assessment_df <- assessment_list$df

    max_df <- dplyr::filter(assessment_df, total_var == base::max(total_var))

    activation <- max_df$activation[1]
    bottleneck <- base::as.character(max_df$bottleneck[1]) %>% base::as.numeric()

    msg <- glue::glue("Assessment done. Running autoencoder with: \nActivation function: '{activation}'\nBottleneck neurons: {bottleneck} ")

    confuns::give_feedback(
      msg = msg,
      verbose = verbose
    )

  } else {

    assessment_list <- base::tryCatch({

      getAutoencoderAssessment(object, of_sample = of_sample)

    }, error = function(error){

      return(list())

    })

  }

  # -----

  # 2. Create network --------------------------------------------------------

  input_layer <-
    keras::layer_input(shape = c(ncol(x_train)))

  encoder <-
    input_layer %>%
    keras::layer_dense(units = layers[1], activation = activation) %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[2], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[3], activation = activation) %>%
    keras::layer_dense(units = bottleneck)

  decoder <-
    encoder %>%
    keras::layer_dense(units = layers[3], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[2], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[1], activation = activation) %>%
    keras::layer_dense(units = c(ncol(x_train)))

  autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

  autoencoder_model %>% keras::compile(
    loss = 'mean_squared_error',
    optimizer = 'adam',
    metrics = c('accuracy')
  )

  history <-
    autoencoder_model %>%
    keras::fit(x_train, x_train, epochs = epochs, shuffle = TRUE,
               validation_data = list(x_train, x_train), verbose = verbose)

  reconstructed_points <-
    autoencoder_model %>%
    keras::predict_on_batch(x = x_train)

  base::rownames(reconstructed_points) <- base::rownames(x_train)
  base::colnames(reconstructed_points) <- base::colnames(x_train)

  if(base::isTRUE(display_plot)){

    plot_df <-
      base::rbind(
        data.frame(base::t(reconstructed_points[val_genes, ]), type = "Denoised"),
        data.frame(base::t(x_train[val_genes, ]), type = "Scaled")
      ) %>%
      dplyr::mutate(type = base::factor(x = type, levels = c("Scaled", "Denoised")))

    val_plot <-
      ggplot2::ggplot(data = plot_df, ggplot2::aes(x = .data[[val_genes[1]]], y = .data[[val_genes[2]]], color = type)) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
      ggplot2::facet_wrap(. ~ type, scales = "free") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      scale_color_add_on(variable = "discrete", clrp = "milo")

    plot(val_plot)

  }

  # 3. Return updated object ------------------------------------------------

  set_up <- list("activation" = activation,
                 "bottleneck" = bottleneck,
                 "dropout" = dropout,
                 "epochs" = epochs,
                 "input_mtr" = "scaled",
                 "output_mtr" = mtr_name_output,
                 "layers" = layers)

  object <- addAutoencoderSetUp(object = object,
                                mtr_name = mtr_name_output,
                                set_up_list = set_up,
                                of_sample = of_sample)

  object <- addExpressionMatrix(object = object,
                                mtr_name = mtr_name_output,
                                expr_mtr = reconstructed_points,
                                of_sample = of_sample)

  object <- computeGeneMetaData(object,
                                mtr_name = mtr_name_output,
                                of_sample = of_sample)

  object <-
    setActiveExpressionMatrix(object = object, mtr_name = mtr_name_output)


  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )



  return(object)

}

# -----

# print-functions ---------------------------------------------------------


#' @title Print autoencoder summary
#'
#' @description Prints a human readable summary about the set up of the last neural network that
#' was constructed to generate a denoised expression matrix.
#'
#' @inherit check_sample params
#'
#' @inherit print_family return
#' @export

printAutoencoderSummary <- function(object, mtr_name = "denoised", of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  info_list <- object@information$autoencoder[[of_sample]][["nn_set_ups"]]

  info_list <- getAutoencoderSetUp(object = object, of_sample = of_sample, mtr_name = mtr_name)

  if(base::is.null(info_list)){

    base::stop("Could not find any information. It seems as if function 'runAutoEncoderDenoising()' has not been run yet.")

  }

  feedback <- glue::glue("{introduction}: \n\nActivation function: {activation}\nBottleneck neurons: {bn}\nDropout: {do}\nEpochs: {epochs}\nLayers: {layers}",
                         introduction = glue::glue("The neural network that generated matrix '{mtr_name}' was constructed with the following adjustments"),
                         activation = info_list$activation,
                         bn = info_list$bottleneck,
                         do = info_list$dropout,
                         epochs = info_list$epochs,
                         layers = glue::glue_collapse(x = info_list$layers, sep = ", ", last = " and "))

  base::return(feedback)

}


# -----



# set-functions -----------------------------------------------------------

#' @title Set results of autoencoder assessment
#'
#' @inherit check_object params
#' @param assessment_list Named list with slots \code{$df} and \code{$set_up}.
#'
#' @return A spata-object.

setAutoencoderAssessment <- function(object, assessment_list, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_data_frame(
    df = assessment_list$df,
    var.class = list("activation" = c("character", "factor"),
                     "bottleneck" = c("character", "factor"),
                     "total_var" = c("numeric", "integer", "double")),
    ref = "assessment_list$df"
  )

  object@autoencoder[[of_sample]][["assessment"]] <- assessment_list

  base::return(object)

}


# -----

