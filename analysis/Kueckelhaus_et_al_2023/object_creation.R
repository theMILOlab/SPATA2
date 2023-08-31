
library(tidyverse)
library(devtools)
load_all()

dir <- "E:/Lab/Data/Milo Lab/10XVisium"

dir_sample <- str_c(dir, "/#UKF265_C/outs")

if(FALSE){

  object <- initiateSpataObjectVisium(directory_visium = dir_sample, "C265")

  object <- runImagePipeline(object)

  object <- processWithSeurat(object)

  object <- runBayesSpaceClustering(object)

  object <- setSpataDir(object, dir = str_c(dir, "/spata_object.RDS"))

  saveSpataObject(object)

}

object <- loadSpataObject(str_c(dir, "/spata_object.RDS"))

object <- runBayesSpaceClustering(object, name = "bspace_7", overwrite = T, qs = 7)

object <- setSpataDir(object, dir = str_c(dir_sample, "/spata_object.RDS"))

saveSpataObject(object)

saveSpataObject(object)


