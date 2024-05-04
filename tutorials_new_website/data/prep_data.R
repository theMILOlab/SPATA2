

obj_old <- readRDS("E:/Lab/Manuscripts/SPATA_Kueckelhaus_et_al_2022/data/processed/269_T.RDS")

cnv_res <- obj_old@cnv[[1]]

obj_old <- readSpataObject("275_T")
cnv_res <- getCnvResults(obj_old)

meta_df <- obj_old@meta_obs

meta_df <- left_join(meta_df, y = bayes_space_df[[1]][,c("barcodes", "bayes_space")], by = "barcodes")

object <-
  initiateSpataObjectVisium(
    sample_name = "UKF275T",
    directory_visium = visium_output_dir("275_T_P")
  )

object <- identifyTissueOutline(object, method = "obs")
object <- identifySpatialOutliers(object, method = "dbscan")
object <- removeSpatialOutliers(object)
object <- removeGenesStress(object)
object <- removeGenesZeroCounts(object)

object <- setCnvResults(object, cnv_list = cnv_res)

object <- processImage(object, img_name = "lowres")

object <- addFeatures(object, feature_df = meta_df)

object_t313 <- saveSpataObject(object_t313, "tutorials_new_website/data/object_UKF313T.RDS")


getLogfileDf(object)

test_fn <- function(object, x = 3, y = 44, z = SPATA2::cnv_ref, list_in = 33, df_in = 22){

  returnObject(object)

}






x <- test_fn(object = yolo, x = 1, list_in = list(x = iris), df_in = mtcars)

ce <- rlang::caller_env()

fn <- rlang::caller_fn()

fn_frame <- base::sys.parent()

init_call <- base::sys.call(which = fn_frame)

fn_name <- base::as.character(init_call)[1]

init_args <-
  rlang::fn_fmls_names(fn = fn)

init_args <- init_args[init_args != "..."]

init_args_input <-
  purrr::map(
    .x = init_args,
    .f = function(arg){

      base::parse(text = arg) %>%
        base::eval(envir = ce)

    }
  ) %>%
  purrr::set_names(nm = init_args)

init_args_input <-
  init_args_input[!base::names(init_args_input) %in% c("cds",  "coords_df", "count_mtr", "expr_mtr", "object", "seurat_object")]

init_args_input <-
  c(init_args_input, additional_input)

initiation_list <- list(
  fn = fn_name,
  input = init_args_input,
  time = base::Sys.time()
)

object@obj_info$initiation <- initiation_list

return(object)
