
library(tidyverse)
library(devtools)
load_all()

meta_df <- readRDS(file = "E:/Lab/Data/Kuppe_et_al_2022/meta_df_processed.RDS")

for(i in nrow(meta_df):nrow(meta_df)){

  print(glue::glue("Working on {i}/{nrow(meta_df)}. Sample: {meta_df$hca_sample_id[i]}."))

  object <-
    initiateSpataObjectVisium(
      directory_visium = str_c(meta_df$sample_folder[i], "/outs"),
      sample_name = meta_df$hca_sample_id[i]
    )

  object <- runImagePipeline(object)

  object <- processWithSeurat(object)

  object <- setActiveMatrix(object, "scaled")

  object <- runBayesSpaceClustering(object,name = "bayes_space", nrep = c(1001, 5000), overwrite = T)

  object <- runDEA(object, across = "bayes_space")

  #object <- runGSEA(object, across = "bayes_space")

  object <- setSpataDir(object, dir = stringr::str_c(meta_df$sample_folder[i], "/spata2_object.RDS"))

  saveSpataObject(object)

}
