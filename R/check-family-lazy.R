#' @include S4-generic-functions.R
NULL

#' @title lazy_check_dummy
#'
#' @description A member of the \code{lazy-check_*()}-family.
#'
#' @details Members of the \code{lazy-check_*()}-family take the arguments
#' of their kind - that are used in the function they are called in - and
#' checks whether these arguments input fit the requirements. They stop and return an
#' error immediately once they stumble upon something invalid. They do not alter or adjust input
#' and return TRUE if the whole function has been executed without anything
#' invalid being found.
#'
#' @return A logical value TRUE if nothing invalid has been detected or an informative
#' error message.

lazy_check_dummy <- function(){}



#################################################################################################


# Miscellaneous -----------------------------------------------------------

#' @title Gives feedback if object was not found
#'
#' @param test Logical test that checks if something was found.
#' @param ref_x Object reference.
#' @param ref_fns Function(s) that needs


check_availability <- function(test, ref_x, ref_fns){

  if(!base::isTRUE(test)){

    msg <- glue::glue("Could not find {ref_x}. Make sure to use {ref_fns}.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

  }

  base::invisible(TRUE)

}

#' @title Gives feedback about input validity
#'
#' @param to Character value. Denotes the object to which
#' the directory leads. Must be one of \emph{'cell_data_set', 'seurat_object'}
#' or \emph{'spata_object'}.
#'
#' @inherit lazy_check_dummy return description details
#' @export

check_to <- function(object, to){

  confuns::check_one_of(
    input = to,
    against = validDirectoryInstructionSlots(),
    ref.input = "input for argument 'to'"
  )

  not_defined_directories <-
    purrr::keep(.x = to, .p = ~ object@information$instructions$directories[[.x]]  == "not defined")

  if(base::length(not_defined_directories) >= 1){

    msg <- glue::glue("The {ref_dir} for '{not_defined}' {ref_have} not been defined yet.",
                      ref_dir = confuns::adapt_reference(not_defined_directories, sg = "directory", pl = "directories"),
                      ref_have = confuns::adapt_reference(not_defined_directories, sg = "has", pl = "have"),
                      not_defined = glue::glue_collapse(x = not_defined_directories, sep = "', '", last = "' and '"))

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "warning",
      with.time = FALSE
    )

  }


  base::invisible(TRUE)

}

# Recurring data.frame input ----------------------------------------------

#' @title Check assessed trajectory data.frame
#'
#' @param atdf A data.frame containing the results of trajectory-modelling. Must contain the variables:
#'
#'  \describe{
#'   \item{\emph{variables}}{Character. The genes, gene-sets and features that have been assessed.}
#'   \item{\emph{pattern}}{Character. The respective trajectory pattern.}
#'   \item{\emph{auc}}{Numeric. The assessment value which is the residual's area under the curve.}
#'  }
#'
#' @inherit lazy_check_dummy description details return


check_atdf <- function(atdf){

  confuns::check_data_frame(
    df = atdf,
    var.class = list(
      variables = c("character"),
      pattern = c("character", "factor"),
      auc = c("numeric", "integer", "double")
    ),
    ref = "atdf")

}

#' @rdname check_atdf
check_rtdf <- function(rtdf, variable = NULL){

  # check classes
  confuns::check_data_frame(df = rtdf,
                            var.class =
                              list(
                                variables = "character",
                                data = "list",
                                residuals = "list",
                                auc = "list"),
                            ref = "rtdf")

  base::return(base::invisible(TRUE))

}

#' @title Check coordinate data.frame
#'
#' @param spata_df A data.frame containing information of barcode-spots. Must contain the variables.
#'
#' \describe{
#'  \item{\emph{barcodes}}{Character. The barcode-sequences (+ the sample belonging) of every barcode spot.}
#'  \item{\emph{sample}}{Character. The sample belonging of every barcode-spot.}
#'  }
#'
#' @inherit lazy_check_dummy description details return

check_spata_df <- function(spata_df){

  confuns::check_data_frame(
    df = spata_df,
    var.class = list(
      barcodes = c("character"),
      sample = c("character")
    ),
    ref = "spata_df"
  )

}


#' @title Check summarized trajectory data.frame
#'
#' @param stdf A data.frame containing information about a spatial trajectory. Must contain the variables:
#'
#'  \describe{
#'   \item{\emph{trajectory_part}}{Character. The trajectory's subparts.}
#'   \item{\emph{trajectory_part_order}}{Numeric. Denotes every trajectory-bin's position along the trajectory's subpart.}
#'   \item{\emph{trajectory_order}}{Numeric. Denotes every trajectory-bin's position along the whole trajectory.}
#'   \item{\emph{variables}}{Character. The genes, gene-sets and features that have been summarized along the trajectory.}
#'   \item{\emph{values}}{Numeric. The binned values of every gene, gene-set and feature that has been summarized along the trajectory.}
#'   }
#'
#' @inherit lazy_check_dummy description details return

check_stdf <- function(stdf, shift = NULL){

  if(!base::is.null(shift)){ confuns::check_one_of(input = shift, against = c("wider", "longer"))}

  if(base::is.null(shift) || shift == "wider"){

    confuns::check_data_frame(
      df = stdf,
      var.class = list(
        trajectory_part = c("character"),
        trajectory_part_order = c("integer", "numeric", "double"),
        trajectory_order = c("integer", "numeric", "double"),
        variables = c("character"),
        values = c("integer", "numeric", "double")
      ),
      ref = "stdf"
    )

  } else if(!base::is.null(shift) && shift == "longer"){

    cnames <- base::colnames(stdf)

    stdf_variables <- cnames[!cnames %in% trajectory_df_colnames]

    all_variables <- c(stdf_variables, "trajectory_order", "trajectory_part_order")

    var.class <-
      purrr::map(.x = all_variables,
                 .f = ~ c("integer", "numeric", "double")) %>%
      purrr::set_names(nm = all_variables)

    confuns::check_data_frame(
      df = stdf,
      var.class = var.class,
      ref = "stdf"
    )

  }

}


#' @title Check de data.frame
#' @param dea_df A data.frame containing information about differentially expressed genes. Must contain the variables:
#'
#'  \describe{
#'   \item{\emph{gene}}{Character. The differentially expressed genes.}
#'   \item{\emph{cluster}}{Character. The clusters (or experimental groups) across which the analysis was performed.}
#'   \item{\emph{avg_logFC}}{Numeric. The average log-fold change to which the belonging gene was differentially expressed..}
#'   \item{\emph{p_val}}{Numeric. The p-values.}
#'   \item{\emph{p_val_adj}}{Numeric. The adjusted p-values.}
#'  }
#'
#'Hint: Use the resulting data.frame of \code{SPATA::findDE()} or it's descendants as input.
#'
#' @inherit lazy_check_dummy description details return

check_dea_df <- function(dea_df){

  confuns::check_data_frame(
    df = dea_df,
    var.class = list(
      p_val = "numeric",
      avg_logFC = "numeric",
      p_val_adj = "numeric",
      gene = "character"
    ),
    ref = "de_df"
  )

}


#' @title Check coords data.frame
#'
#' @param coords_df A data.frame containing information about every barcode-spot. Must contain the variables:
#'  \itemize{
#'   \item{\emph{barcodes} Character. The barcode-sequences (+ the sample belonging) of every barcode spot.}
#'   \item{\emph{sample} Character. The sample belonging of every barcode-spot.}
#'   \item{\emph{x} Numeric. The x-coordinates of every barcode.}
#'   \item{\emph{y} Numeric. The y-coordinates of every barcode.}
#'  }
#'
#' @inherit lazy_check_dummy description details return

check_coords_df <- function(coords_df){

  confuns::check_data_frame(
    df = coords_df,
    var.class = list(
      barcodes = "character",
      x = c("numeric", "integer", "double"),
      y = c("numeric", "integer", "double")
    ),
    ref = "coords_df"
  )

}



#' @title Check customized trend list
#'
#' @param length_trajectory Numeric value. Length of trajectory according to which the trends have been
#' customized.
#' @param customized_trends A data.frame or a named list. All numeric variables are considered to correspond to customized trends
#' the trajectory of interest might adopt. The names of the respective variables will correspond to the name
#' with which you want to refer to the trend later on.

check_customized_trends <- function(length_trajectory,
                                    customized_trends){

  if(!base::is.list(customized_trends)){

    base::stop("Input for argument 'customized_trends' must be a named list or a data.frame.")

  }

  # keep only numeric slots
  all_numerics <-
    purrr::keep(.x = customized_trends, .p = ~ base::is.numeric(.x))

  # check names
  trend_names <- base::names(all_numerics)

  if(base::is.null(trend_names) | base::length(trend_names) != base::length(all_numerics)){

    base::stop("Please make sure that all numeric slots of the list 'customized_trends' are named.")

  }

  # check lengths
  all_lengths <-
    purrr::map_int(.x = all_numerics, .f = ~ base::length(.x))

  if(dplyr::n_distinct(all_lengths) != 1){

    base::stop("Please make sure that all numeric slots of the list 'customized_trends' are of the same length.")

  }

  # compare length of trajectory with length of customized trends

  if(base::is.numeric(length_trajectory)){

    length_customized_trends <- base::unique(all_lengths)

    if(length_trajectory != length_customized_trends){

      base::stop(glue::glue("Please make sure that the lengths of the customized trends are equal to the length of the trajectory (= {length_trajectory})."))

    }

  }

  # check for nas
  has_nas <-
    purrr::map(.x = all_numerics, .f = ~ base::is.na(.x) %>% base::sum()) %>%
    purrr::keep(.x = ., .p = ~ .x > 0)

  if(base::length(has_nas) >= 1){

    slots_with_nas <- stringr::str_c(base::names(has_nas), collapse = "', '")

    base::warning(glue::glue("Ignoring slots '{slots_with_nas}' as they contain NAs."))

  }

  no_nas <-
    purrr::keep(.x = all_numerics, .p = ~ base::is.na(.x) %>% base::sum() == 0) %>%
    purrr::map_df(.x = . , .f = ~ confuns::normalize(x = .x))

  base::return(no_nas)

}

# -----

#' @title Check assign input
#'
#' @param assign Logical. If set to TRUE a named list will be assigned to the global
#' environment. This list contains data and information to rebuild or additionally
#' customize the output plot of this function.
#' @param assign_name The name the assigned list is supposed to have specified as
#' a single character value.
#'
#' @inherit lazy_check_dummy description details return

check_assign <- function(assign = FALSE,
                         assign_name = character(1)){


  if(!base::is.logical(assign)){

    base::stop("Argument 'assign' needs to be logical.")

  }

  if(base::isTRUE(assign)){

    if(!base::is.character(assign_name) | !base::length(assign_name) == 1){

      base::stop("Argument 'assign_name' needs to be a single character value.")

    }

    if(assign_name == ""){

      base::stop("Argument 'assign_name' must not be ''.")

    }

    if(base::exists(x = assign_name, where = .GlobalEnv)){

      base::stop(stringr::str_c("It already exists an object named '",
                                assign_name, "' in the global environment.",
                                sep = ""))

    }


  }

  base::return(TRUE)

}




#' @title Check compiled trajectory data.frame
#'
#' @param ctdf A compiled trajectory data.frame containing the variables
#' \emph{'barcodes', 'sample', 'x', 'y', 'projection_length', 'trajectory_part'}.
#'
#' @inherit lazy_check_dummy description details return

check_compiled_trajectory_df <- function(ctdf){

  check_spata_df(spata_df = ctdf)
  check_coordinate_variables(data = ctdf, x = "x", y = "y")

  vc <- confuns::variable_classes2(data = ctdf)

  if(!base::all(c("projection_length", "trajectory_part") %in% base::names(vc))){
    base::stop("Variables must contain 'projection_length' and 'trajectory_part'.")
  }

  if(vc["projection_length"] != "numeric"){
    base::stop("Variable 'projection_length' needs to be of class numeric.")
  }

  if(vc["trajectory_part"] != "character"){
    base::stop("Variable 'projection_length' needs to be of class character.")
  }

}



#' @title Check coordinate variables
#'
#' @param data A data.frame containing the variables of interest as well
#' as the variables needed to map onto the x- and y axis of the plot.
#' @param x The name of the numeric variable to be plotted on the x axis.
#' @param y The name of the numeric variable to be plotted on the y axis.
#'
#' @inherit lazy_check_dummy description details return

check_coordinate_variables <- function(data, x = "x", y = "y"){

  if(!base::all(c(x, y) %in% base::colnames(data))){

    base::stop(glue::glue("Provided data.frame needs to have numeric variables '{x}' and '{y}'."))

  }

  if(base::any(!base::sapply(X = data[,c(x,y)], FUN = base::is.numeric))){

    base::stop(glue::glue("Both variables '{x}' and '{y}' of 'data' need to be numeric."))

  }

  base::return(TRUE)

}




#' @title Check display input
#'
#' @param display_image Logical. If set to TRUE the histology image of the specified sample
#' is displayed underneath the plot.
#'
#' @param display_title Logical. If set to TRUE an informative title is displayed.
#'
#' @inherit lazy_check_dummy description details return

check_display <- function(display_title = FALSE,
                          display_image = FALSE){

  if(!base::is.logical(display_title)){

    base::stop("Argument 'display_title' needs to be logical.")

  }

  if(!base::is.logical(display_image)){

    base::stop("Argument 'display_image' needs to be logical.")

  }

}



#' @title Check feature data.frame
#'
#' @param feature_name Character value. The name of the feature that is to be added
#'  to the object.
#' @param feature_df A data.frame that contains the feature and the key-variables.
#'
#' @inherit lazy_check_dummy description details return

check_feature_df <- function(feature_name,
                             feature_df){

  if(!base::length(feature_name) == 1 | !base::is.character(feature_name)){

    base::stop("Argument 'feature_name' needs to be a single character value.")

  }

  if(!base::is.data.frame(feature_df)){

    base::stop("Argument 'feature_df' needs to be a data.frame.")

  } else if(!base::all(c("barcodes", "sample", feature_name) %in% base::colnames(feature_df))){

    base::stop(glue::glue("Data.frame 'feature_df' needs to have the variables 'barcodes', 'sample' and '{{feature_name}}'."))

  } else {

    classes <- base::sapply(X = feature_df[,c("barcodes", "sample")],
                            FUN = base::class)

    if(!base::all(classes == "character")){

      base::stop("Variables 'barcodes' and 'sample' need to be of class character.")

    } else {

      base::return(TRUE)

    }

  }

}



#' @title Check method input
#'
#' @param method_dr Character value. The dimensional reduction method of
#' interest specified as a single character value. (Currently
#' either \emph{'pca'}, \emph{'umap'} or \emph{'tsne'}).
#' @param method_gs Character value. The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{'mean'} or one
#' of \emph{'gsva', 'ssgsea', 'zscore', or 'plage'}. The latter four will be given to
#' \code{gsva::GSVA()}.
#' @param method_hclust Character value. Denotes the hierarchical clustering method  according
#' to which the clustering is performed. Valid input options are \emph{'ward.D', 'ward.D2', 'single',
#'  'complete', 'average', 'mcquitty', 'median'} and \emph{'centroid'}.
#' @param method_padj Character value. The method according to which the adjusted p-values will
#' be calculated. Given to \code{stats::p.adjust()}. Run \code{stats::p.adjust.methods} to obtain
#' all valid input options.
#' @param method_de Character value. Denotes the method to according to which the de-analysis is performed.
#' Given to argument \code{test.use} of the \code{Seurat::FindAllMarkers()}-function. Run \code{SPATA::dea_methods}
#' to obtain all valid input options.
#' @param method_ovl Character value. One of \emph{'classic', 'bayesian'}. Decides
#' according to which method the spatial overlap is calculated.
#' @inherit lazy_check_dummy description details return
#' @export

check_method <- function(method_csr = NULL,
                         method_de = NULL,
                         method_dr = NULL,
                         method_gs = NULL,
                         method_hclust = NULL,
                         method_ovl = NULL,
                         method_padj = NULL
                         ){


  # complete spatial randomness ---------------------------------------------

  if(!base::is.null(method_csr)){

    if(confuns::is_value(x = method_csr, mode = "character")){

      confuns::check_one_of(
        input = method_csr,
        against = csr_methods,
        ref.input = "for argument 'method_csr'"
      )

    }


  }



  # differential expression methods -----------------------------------------

  if(!base::is.null(method_de)){

    if(confuns::is_value(x = method_de, mode = "character")){

      confuns::check_one_of(
        input = method_de,
        against = de_methods,
        ref.input = "for argument 'method_de'"
      )

    }

  }

  # -----

  # dimensional reduction methods -------------------------------------------

  if(!base::is.null(method_dr)){

    if(confuns::is_value(x = method_dr, mode = "character")){

      confuns::check_one_of(
        input = method_dr,
        against = dim_red_methods,
        ref.input = "argument 'method_dr'"
      )

    }

  }

  # -----


  # gene set methods --------------------------------------------------------

  if(!base::is.null(method_gs)){

    if(confuns::is_value(x = method_gs, mode = "character")){

      confuns::check_one_of(
        input = method_gs,
        against = gene_set_methods,
        ref.input = "for argument 'method_gs'"
      )

    }

  }

  # -----


  # hierarchical clustering methods -----------------------------------------

  if(!base::is.null(method_hclust)){

    if(confuns::is_value(x = method_hclust, mode = "character")){

      confuns::check_one_of(
        input = method_hclust,
        against = hclust_methods,
        ref.input = "for argumnet 'method_hclust'"
      )

    }


  }

  # find overlap methods ----------------------------------------------------

  if(!base::is.null(method_ovl)){

    if(confuns::is_value(x = method_ovl, mode = "character")) {

      confuns::check_one_of(
        input = method_ovl,
        against = c("classic", "bayesian"),
        ref.input = "for argument 'method_ovl'"
      )

    }

  }

  # adjusted pvalue methods -------------------------------------------------

  if(!base::is.null(method_padj)){

    if(confuns::is_value(x = method_padj, mode = "character")){

      confuns::check_one_of(
        input = method_padj,
        against = stats::p.adjust.methods,
        ref.input = "for argument 'method_padj'"
      )

    }

  }

  # -----

}



#' @title Check monocle input
#'
#' @param preprocess_method Monocle3 - description:
#'
#' A string specifying the initial dimension method to use,
#' currently either PCA or LSI. For LSI (latent semantic
#' indexing), it converts the (sparse) expression matrix into
#' tf-idf matrix and then performs SVD to decompose the gene
#' expression / cells into certain modules / topics. Default
#' is "PCA".
#'
#' @param reduction_method Monocle3 - description:
#'
#' A character string specifying the algorithm to use for
#' dimensionality reduction. Currently "UMAP", "tSNE", "PCA"
#' and "LSI" are supported.
#'
#' @param cluster_method Monocle3 - description:
#'
#' String indicating the clustering method to use. Options are
#' "louvain" or "leiden". Default is "leiden". Resolution parameter
#' is ignored if set to "louvain".
#'
#' @param k Monocle3 - description:
#'
#' Integer number of nearest neighbors to use when creating the k
#' nearest neighbor graph for Louvain/Leiden clustering. k is related
#' to the resolution of the clustering result, a bigger k will result
#' in lower resolution and vice versa. Default is 20.
#'
#' @param num_iter Monocle3 - description:
#'
#' Integer number of iterations used for Louvain/Leiden clustering.
#' The clustering result giving the largest modularity score will be
#' used as the final clustering result. Default is 1. Note that if
#' num_iter is greater than 1, the random_seed argument will be ignored
#' for the louvain method.
#'
#' @return
#'
#' @details With respect to the arguments \code{preprocess_method},
#' \code{reduction_method} and \code{cluster_method}:
#'
#' If you provide a vector instead of a single character value (string)
#' the function will iterate over all inputs via a for-loop to compute all
#' valid combinations.
#'

check_monocle_input <- function(preprocess_method,
                                reduction_method,
                                cluster_method,
                                k = base::numeric(1),
                                num_iter = base::numeric(1)){

  confuns::is_value(k, "numeric", "k")
  confuns::is_value(num_iter, "numeric", "num_iter")

  msg <- NULL

  if(!base::is.null(preprocess_method) &&
     base::any(!preprocess_method %in% c("PCA", "LSI"))){

    msg <- "Invalid input for argument 'preprocess_method'. Valid input options are: 'PCA', 'LSI'."

  }

  if(!base::is.null(reduction_method) &&
     base::any(!reduction_method %in% c("UMAP", "tSNE", "PCA", "LSI"))){

    msg <- "Invalid input for argument 'reduction_method'. Valid input options are: 'UMAP', 'tSNE', 'PCA', 'LSI'."

  }

  if(!base::is.null(cluster_method) &&
     base::any(!cluster_method %in% c("louvain", "leiden"))){

    msg <- "Invalid input for argument 'cluster_method'. Valid inputs are: 'louvain', 'leiden'."

  }

  confuns::give_feedback(
    msg = msg,
    fdb.fn = "stop"
  )

  base::return(TRUE)

}

#' Makes sure that both packages are installed
#'
check_monocle_packages <- function(){

  pkgs <-
    utils::installed.packages() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column("pkgs") %>%
    dplyr::pull(pkgs)

  missing <- base::setdiff(x = c("leidenbase", "monocle3"), y = pkgs)

  if(base::length(missing) >= 1){

    msg <- glue::glue("The {ref1} '{ref2}' must be installed in order to use this function.",
                      ref1 = confuns::adapt_reference(input = missing, sg = "package"),
                      ref2 = glue::glue_collapse(x = missing, sep = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

}


#' Check spata object input
#'
#' @param object A valid spata-object.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_object <- function(object){

  validation(object)

  base::invisible(TRUE)

}



#' @title Check pt input

#' @param pt_alpha Numeric value. Specifies the degree of transparency of all points.
#' @param pt_size Numeric value. Specifies the size of all points.
#' @param pt_clr Character value. Specifies the color of all points.
#' @param pt_clrp The color palette to be used if the specified variable displayed by
#' color is categorical/discrete. Run \code{validColorPalettes()} to see valid input.
#' @param pt_clrsp The color spectrum to be used if the specified variable displayed by
#' color is continuous. Run \code{validColorSpectra()} to see valid input.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_pt <- function(pt_size = NULL,
                     pt_alpha = NULL,
                     pt_clrsp = NULL,
                     pt_clrp = NULL,
                     pt_clr = NULL){

  msg <- NULL

  confuns::are_values(c("pt_clrp", "pt_clrsp"), mode = "character", skip.allow = TRUE, skip.val = NULL)
  confuns::are_values(c("pt_size", "pt_alpha"), mode = "numeric", skip.allow = TRUE, skip.val = NULL)

  if(!base::is.null(pt_clrsp) && !pt_clrsp %in% base::unlist(validColorSpectra(), use.names = FALSE)){

    msg <- "Invalid input for argument 'pt_clrsp'. Run validColorSpectra() to see all valid input options."

  }

  if(!base::is.null(pt_clrp) && !pt_clrp %in% base::unlist(validColorPalettes(), use.names = FALSE)){

    msg <- "Invalid input for argument 'pt_clrp'. Run validColorPalettes() to see all valid input options."

  }

  confuns::give_feedback(
    msg = msg,
    fdb.fn = "stop"
  )

  base::return(TRUE)

}


#' @title Check smooth input
#'
#' @param df A data.frame that is to be smoothed spatially. That data frame must have
#' numeric \emph{x}- and \emph{y}-variables.
#' @param smooth Logical. If set to TRUE values will be smoothed according to the
#' \code{smoooth_}-parameters.
#' @param smooth_clr Character value. The color with which to display the smoothed
#' model.
#' @param smooth_method The smoothing method that will be used specified as a
#' single character value (e.g. \emph{"lm", "glm", "gam", "loess"}).
#' @param smooth_se Logical. If set to TRUE the confidence interval will be
#' @param smooth_span The amount of smoothing specified as a single numeric value.
#' displayed.
#'
#' @inherit lazy_check_dummy description details return

check_smooth <- function(df = NULL,
                         smooth = NULL,
                         smooth_span = NULL,
                         smooth_method = NULL,
                         smooth_se = NULL){

  if(!base::is.null(smooth) &&
     !base::isTRUE(smooth) &
     !base::isFALSE(smooth)){

    base::stop("Argument 'smooth' needs to be TRUE or FALSE.")

  }

  if(!base::is.null(smooth) && base::isTRUE(smooth)){

    if(!base::is.null(df) &&
       !base::all(c("x", "y") %in% base::colnames(df))){

      base::stop("Input data.frame doesn't contain x and y variables." )

    }

  }

  if(!base::is.null(smooth_span) && !base::is.numeric(smooth_span)){

    base::stop("Argument 'smooth_span' needs to be numeric.")

  }

  if(!base::is.null(smooth_method)){

    if(!base::is.character(smooth_method) |
       !base::length(smooth_method) == 1){

      base::stop("Argument 'smooth_method' needs to be a single character value.")

    }

  }

  if(!base::is.null(smooth_se) &&
     !base::isTRUE(smooth_se) &
     !base::isFALSE(smooth_se)){

    base::stop("Argument 'smooth_se' needs to be TRUE or FALSE.")

  }

  base::return(TRUE)

}



#' @title Check stdf-input
#'
#' @param stdf A summarized trajectory data.frame
#'
#' @inherit lazy_check_dummy description details return

check_summarized_trajectory_df <- function(stdf){

  warning("check_summarized_trajectory_df is deprecated. Use check_stdf()")

  confuns::check_data_frame(
    df = stdf, var.class = list(
      trajectory_part = "character",
      trajectory_part_order = c("numeric", "integer"),
      trajectory_order = c("numeric", "integer"),
      values = c("numeric", "integer"),
      variables = "character"
    ),
    ref = "stdf")

}



#' Check trajectory name input
#'
#' @inherit check_sample params
#' @param trajectory_name The trajectory of interest specified
#' as a single character value.
#'
#'
#' @inherit lazy_check_dummy description details return

check_trajectory <- function(object,
                             trajectory_name,
                             of_sample = NA){

  confuns::is_value(x = trajectory_name, mode = "character")

  if(!trajectory_name %in% getTrajectoryNames(object, of_sample = of_sample)){

    base::stop(stringr::str_c("There is no trajectory of name '", trajectory_name,
                              "' in sample '", of_sample, "'.", sep = ""))

  }

}


#' @title Check trajectory binwdith input
#'
#' @param binwidth Numeric value. Denotes the binwidth with which to sort all
#' relevant barcode spots into groups that are then aligned with respect to the
#' chosen trajectory's direction.#'
#'

check_trajectory_binwidth <- function(binwidth){

  confuns::is_value(x = binwidth, mode = "numeric")

}


#' @title Check uniform genes input
#'
#' @param uniform_genes Character value. If set to \emph{'discard'} genes that are
#' uniformly expressed across all barcode-spots of the specified coordinates
#' data.frame are discarded. If set to \emph{'keep'} they are kept.
#'
#' @return

check_uniform_genes <- function(uniform_genes){

  confuns::is_value(uniform_genes, "character", "uniform_genes")

  if(!uniform_genes %in% c("keep", "discard")){

    base::stop("Argument 'uniform genes' must be set to 'keep' or 'discard'.")

  } else {

    base::return(base::invisible(TRUE))
  }

}








