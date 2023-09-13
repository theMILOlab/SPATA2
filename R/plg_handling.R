

increase_polygon_vertices <- function(polygon_df, avg_dist) {

  polygon_df <- base::as.data.frame(polygon_df)

  # Ensure the polygon is closed (first and last point are the same)
  if (!identical(polygon_df[1, ], polygon_df[nrow(polygon_df), ])) {
    polygon_df <- rbind(polygon_df, polygon_df[1, ])
  }

  # Initialize a new data frame to store interpolated vertices
  interpolated_df <- data.frame(x = numeric(0), y = numeric(0))

  # Loop through each pair of consecutive vertices
  for (i in 1:(nrow(polygon_df) - 1)) {
    x1 <- polygon_df[i, "x"]
    y1 <- polygon_df[i, "y"]
    x2 <- polygon_df[i + 1, "x"]
    y2 <- polygon_df[i + 1, "y"]

    # Calculate the distance between the consecutive vertices
    dist_between_vertices <- sqrt((x2 - x1)^2 + (y2 - y1)^2)

    # Calculate the number of interpolated vertices needed
    num_interpolated <- max(1, floor(dist_between_vertices / avg_dist))

    # Calculate the step size for interpolation
    step_x <- (x2 - x1) / (num_interpolated + 1)
    step_y <- (y2 - y1) / (num_interpolated + 1)

    # Add the original vertex to the interpolated data frame
    interpolated_df <- rbind(interpolated_df, data.frame(x = x1, y = y1))

    # Interpolate new vertices between the consecutive vertices
    for (j in 1:num_interpolated) {
      new_x <- x1 + j * step_x
      new_y <- y1 + j * step_y
      interpolated_df <- rbind(interpolated_df, data.frame(x = new_x, y = new_y))
    }
  }

  # Combine the original and interpolated vertices
  new_polygon_df <- rbind(polygon_df, interpolated_df)

  return(new_polygon_df)
}

compute_avg_vertex_distance <- function(polygon_df) {

  # Ensure the polygon is closed (first and last point are the same)

  if (!identical(polygon_df[1, ], polygon_df[nrow(polygon_df), ])) {
    polygon_df <- rbind(polygon_df, polygon_df[1, ])
  }

  # Initialize a vector to store vertex distances
  vertex_distances <- numeric(nrow(polygon_df))

  # Loop through each vertex
  for (i in 1:nrow(polygon_df)) {

    x1 <- polygon_df[i, "x"]
    y1 <- polygon_df[i, "y"]

    # Calculate the distances to all other vertices
    distances <- sqrt((polygon_df$x - x1)^2 + (polygon_df$y - y1)^2)

    # Set the distance to itself to infinity
    distances[i] <- Inf

    # Find the minimum distance
    min_distance <- min(distances)

    # Store the minimum distance in the vector
    vertex_distances[i] <- min_distance
  }

  # Compute the average vertex distance
  avg_distance <- mean(vertex_distances)

  return(avg_distance)
}

compute_avg_dp_distance <- function(object, vars = c("x_orig", "y_orig")){

  getCoordsDf(object) %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    base::as.matrix() %>%
    FNN::knn.dist(data = ., k = 1) %>%
    base::mean()

}

make_bins <- function(numeric_vector, binwidth, neg = FALSE) {

  numeric_vector <- base::abs(numeric_vector)

  # Calculate the minimum and maximum values of the numeric vector
  min_value <- min(numeric_vector)
  max_value <- max(numeric_vector)

  # Create a sequence of breaks (bin edges)
  breaks <- seq(min_value, max_value, by = binwidth)

  # Bin the numeric vector
  bin_indices <- floor((abs(numeric_vector) - min_value) / binwidth) + 1

  if(base::isTRUE(neg)){

    prefix <- "-"

    ranges <-
      stringr::str_c(
        "[",
        prefix,
        bin_indices*base::round(binwidth,2),
        ",",
        prefix,
        (bin_indices-1)*base::round(binwidth,2),
        "]"
      )

  } else {

    prefix <- ""

    ranges <-
      stringr::str_c(
        "[",
        prefix,
        (bin_indices-1)*base::round(binwidth,2),
        ",",
        prefix,
        bin_indices*base::round(binwidth,2),
        "]"
      )

  }

  ranges <- stringr::str_replace(ranges, pattern = "-0", replacement = "0")

  levels_out <-
    ranges[base::order(bin_indices)] %>%
    base::unique()

  if(base::isTRUE(neg)){

    levels_out <- base::rev(levels_out)

  }

  out <- base::factor(ranges, levels = levels_out)

  return(out)

}





