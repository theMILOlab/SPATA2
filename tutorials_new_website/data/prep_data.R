

obj_old <- readRDS("review/data/object_T313.RDS")

obj_old <- readSpataObject("313_T")

object <-
  initiateSpataObjectVisium(
    sample_name = "T313",
    directory_visium = visium_output_dir("313_T_P")
  )

object <- identifyTissueOutline(object, method = "obs")
object <- identifySpatialOutliers(object, method = "dbscan")
object <- removeSpatialOutliers(object)
object <- removeGenesStress(object)
object <- removeGenesZeroCounts(object)

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
