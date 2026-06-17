library(geometry)

#' Preprocess hexagons by finding the convex hull of required nodes
#'
#' @param hexa_coordinates Matrix of centroid coordinates for all hexagons
#' @param all_required_nodes Indices of start and target nodes
#' @param hexa The sf/sp object containing hexagon spatial geometries
#' @param nb The neighborhood object (spdep nb)
#' @return A data frame of coordinates of nodes that lie within the convex hull
convex_hull_preprocess <- function(hexa_coordinates, all_required_nodes, hexa, nb) {
  # Convex hull with geometry library
  start_target_points <- as.matrix(hexa_coordinates[all_required_nodes, ])
  convex_hull_hexagon <- convhulln(start_target_points, options = "FA")
  
  # Plot map topology and neighborhood connections
  plot(hexa$geometry, border = grey(0.8), axes = TRUE)
  plot(nb, coords = hexa_coordinates, pch = 20, cex = 0.5, col = grey(0.5), add = TRUE)

  # Find which points stay inside the convex hull
  hexa_matrix <- as.matrix(hexa_coordinates)
  in_ch <- inhulln(convex_hull_hexagon, hexa_matrix)

  reduced_solution_surface <- 100 - round(nrow(hexa_matrix[in_ch, , drop = FALSE]) * 100 / nrow(hexa_matrix), 2)
  cat("There were ", nrow(hexa_matrix), " nodes.\n", sep = "")
  cat("Preprocessing with convex hull reduced the total number of nodes to ", 
      nrow(hexa_matrix[in_ch, , drop = FALSE]), " nodes (", reduced_solution_surface, "% reduction)\n", sep = "")

  # Plot reduced nodes and convex hull boundary
  points(hexa_matrix[in_ch, , drop = FALSE], col = "red", pch = 21, cex = 2)
  plot(convex_hull_hexagon, col = "blue", lwd = 2, add = TRUE)

  plot_title <- paste0("From: ", nrow(hexa_matrix), ", To: ",
                       nrow(hexa_matrix[in_ch, , drop = FALSE]), " nodes\n",
                       reduced_solution_surface, "% less nodes")
  title(plot_title)

  return(as.data.frame(hexa_matrix[in_ch, , drop = FALSE]))
}

