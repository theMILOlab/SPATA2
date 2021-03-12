

#' @title Wrapper around Seurat processing functions
#'
#' @inherit argument_dummy params
#' @inherit transformSpataToSeurat params
#' @param seurat_object A valid seurat-object.
#'
#'
#' @return A processed seurat-object.
#'

process_seurat_object <- function(seurat_object,
                                  assay_name = NULL,
                                  calculate_rb_and_mt = TRUE,
                                  remove_stress_and_mt = TRUE,
                                  SCTransform = FALSE,
                                  NormalizeData = TRUE,
                                  FindVariableFeatures = TRUE,
                                  ScaleData = TRUE,
                                  RunPCA = TRUE,
                                  FindNeighbors = TRUE,
                                  FindClusters = TRUE,
                                  RunTSNE = TRUE,
                                  RunUMAP = TRUE,
                                  verbose = TRUE){

# 1. Control --------------------------------------------------------------

  base::stopifnot(methods::is(object = seurat_object, class2 = "Seurat"))

  confuns::is_value(x = assay_name, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.null(assay_name)){

    assay_name <- base::names(seurat_object@assays)

    if(base::length(assay_name) != 1){

      msg <- glue::glue("Found more than one assay in provided seurat-object. Please specify one of the options '{ref_assays}' using argument 'assay_name'.",
                        ref_assays = glue::glue_collapse(x = assay_name, sep = "', '", last = "' or '"))

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  }

  for(fn in seurat_process_fns){

    input <- base::parse(text = fn) %>% base::eval()

    if(base::is.data.frame(input) | (!base::isTRUE(input) && !base::is.list(input) &&!base::isFALSE(input))){

      base::stop(glue::glue("Invalid input for argument '{fn}'. Must either be TRUE, FALSE or a named list."))

    }

  }

  # calculate ribosomal and mitochondrial percentage
  if(base::isTRUE(calculate_rb_and_mt)){

    msg <- "Calculating percentage of ribosomal and mitochondrial genes."

    confuns::give_feedback(msg = msg, verbose = verbose)

    seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT.")
    seurat_object[["percent.RB"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^RPS")

  }

  # remove stress and mitochondrial genes
  if(base::isTRUE(remove_stress_and_mt)){

    msg <- "Removing stress genes and mitochondrial genes."

    confuns::give_feedback(msg = msg, verbose = verbose)

    exclude <- c(base::rownames(seurat_object@assays[[assay_name]])[base::grepl("^RPL", base::rownames(seurat_object@assays[[assay_name]]))],
                 base::rownames(seurat_object@assays[[assay_name]])[base::grepl("^RPS", base::rownames(seurat_object@assays[[assay_name]]))],
                 base::rownames(seurat_object@assays[[assay_name]])[base::grepl("^MT-", base::rownames(seurat_object@assays[[assay_name]]))],
                 c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1'))

    feat_keep <- base::rownames(seurat_object@assays[[assay_name]][!(base::rownames(seurat_object@assays[[assay_name]]) %in% exclude), ])

    seurat_object <- base::subset(x = seurat_object, features = feat_keep)

  }



# 2. Process seurat object ------------------------------------------------

  functions_to_call <- seurat_process_fns

  for(fn in functions_to_call){

    input <-
      base::parse(text = fn) %>%
      base::eval()

    if(base::isTRUE(input)){

      msg <- glue::glue("Running 'Seurat::{fn}()' with default parameters.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      args <- base::list("object" = seurat_object)

      if(fn == "ScaleData"){

        args <- base::append(x = args,
                             values = list("features" = base::rownames(seurat_object)))

      }

      # ensure that function is called from Seurat-namespace
      fn <- stringr::str_c("Seurat::", fn, sep = "")

      seurat_object <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = fn)), args),

        error = function(error){

          msg <- glue::glue("Running'Seurat::{fn}()' resulted in the following error: {error$message}. Abort and continue with next function.")

          confuns::give_feedback(msg = msg, verbose = TRUE)

          base::return(seurat_object)

        })

    } else if(base::is.list(input) &
              !base::is.data.frame(input)){

      msg <- glue::glue("Running 'Seurat::{fn}()' with specified parameters.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      args <- purrr::prepend(x = input, values = seurat_object)

      if(fn == "ScaleData" && !"features" %in% base::names(args)){

        args <- base::append(x = args,
                             values = list("features" = base::rownames(seurat_object)))

      }

      # ensure that function is called from Seurat-namespace
      fn <- stringr::str_c("Seurat::", fn, sep = "")

      seurat_object <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = fn)), args),

        error = function(error){

          msg <- glue::glue("Running'Seurat::{fn}()' resulted in the following error: {error$message}. Abort and continue with next function.")

          confuns::give_feedback(msg = msg, verbose = TRUE)

          base::return(seurat_object)

        }

      )

    } else {

      msg <- glue::glue("Skip running '{fn}()' as it's argument input is neither TRUE nor a list.")

      confuns::give_feedback(msg = msg, verbose = verbose)

    }

  }

  base::return(seurat_object)

}









