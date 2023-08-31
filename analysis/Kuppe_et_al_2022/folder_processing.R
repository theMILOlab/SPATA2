
library(tidyverse)

data_folder <- "E:/Lab/Data/Kuppe_et_al_2022/"

meta_df <- read_delim(file = "E:/Lab/Data/Kuppe_et_al_2022/metadata-Visium.csv")

files <-
  list.files("E:/Lab/Data/Kuppe_et_al_2022", full.names = T) %>%
  str_subset(pattern = ".tar.gz$")

meta_df$dir <- character(1)

# unzipping tar.gz. files
for(sample_name in meta_df$hca_sample_id){


  file_ending <- stringr::str_c(sample_name, ".tar.gz$")

  file_pos <- which(str_detect(files, file_ending))

  file <- files[file_pos]

  sample_pos <- which(meta_df$hca_sample_id == sample_name)

  meta_df$dir[sample_pos] <- file

  if(FALSE){

    input_archive <- file
    output_directory <- data_folder

    #dir.create(output_directory)

    # Construct and execute the tar command
    command <- paste("tar -xvzf", shQuote(input_archive), "-C", shQuote(output_directory))
    result <- system(command, intern = TRUE, invisible = FALSE)

  }

}

meta_df$sample_folder <- str_c(data_folder, meta_df$hca_sample_id)

saveRDS(meta_df, file = str_c(data_folder, "meta_df_processed.RDS"))

# unzipping spatial folder
spatial_zip_files <-
  list.files("E:/Lab/Data/Kuppe_et_al_2022", full.names = T, recursive = T) %>%
  str_subset(pattern = "spatial\\.zip$")

spatial_files <-
  c('aligned_fiducials.jpg',
    'detected_tissue_image.jpg',
    'scalefactors_json.json',
    'tissue_hires_image.png',
    'tissue_lowres_image.png',
    'tissue_positions_list.csv'
    )

pb <- create_progress_bar(length(spatial_zip_files))

for(szf in spatial_zip_files){

  pb$tick()

  exdir <-
    confuns::str_extract_before(szf, pattern = "outs\\/") %>%
    stringr::str_c(., "outs")

  zip::unzip(zipfile = szf, exdir = exdir, overwrite = TRUE)

  spatial_folder <- str_c(exdir, "/", "spatial")

  dir.create(spatial_folder)

  for(sf in spatial_files){

    from <-
      list.files(exdir, recursive = T, full.names = T) %>%
      str_subset(pattern = sf)

    dir <- dirname(from)

    filename <-
      stringr::str_remove(from, pattern = dir) %>%
      str_remove(pattern = "\\/")

    file.copy(from = from, to = str_c(spatial_folder, "/", filename))

  }

  print(szf)

}

# unzipping filtered_feature_bc_matrix

feature_bc_files <-
  c("barcodes.tsv.gz",
    "features.tsv.gz",
    "matrix.mtx.gz"
  )

fbc_zip_files <-
  list.files("E:/Lab/Data/Kuppe_et_al_2022", full.names = T, recursive = T) %>%
  str_subset(pattern = "filtered_feature_bc_matrix\\.zip$")

pb <- create_progress_bar(length(fbc_zip_files))

for(file in fbc_zip_files){

  pb$tick()

  exdir <-
    confuns::str_extract_before(file, pattern = "outs\\/") %>%
    stringr::str_c(., "outs")

  zip::unzip(zipfile = file, exdir = exdir, overwrite = TRUE)

  fbc_folder <- str_c(exdir, "/", "filtered_feature_bc_matrix")

  dir.create(fbc_folder)

  for(ff in feature_bc_files){

    from <-
      list.files(exdir, recursive = T, full.names = T) %>%
      str_subset(pattern = ff)

    dir <- dirname(from)

    filename <-
      stringr::str_remove(from, pattern = dir) %>%
      str_remove(pattern = "\\/")

    file.copy(from = from, to = str_c(fbc_folder, "/", filename))

  }

  print(file)

}


# renaming matrices
for(i in 2:nrow(meta_df)){

  sample_id <- meta_df[i,][["hca_sample_id"]]

  dir <-
    stringr::str_c("E:/Lab/Data/Kuppe_et_al_2022/", sample_id, "/outs")

  file.rename(
    from = stringr::str_c(dir, "/spatial_", sample_id, "_filtered_feature_bc_matrix.h5"),
    to = stringr::str_c(dir, "/filtered_feature_bc_matrix.h5")
  )

}


