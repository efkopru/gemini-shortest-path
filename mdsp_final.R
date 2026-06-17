rm(list=ls())       # Clear environment
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(stringr)
library(rstudioapi)
library(dplyr)
library(sf)
library(spdep)
library(sp)

# Robust path resolution for running inside RStudio or via command line (Rscript)
get_script_path <- function() {
  if (Sys.getenv("RSTUDIO") == "1" && requireNamespace("rstudioapi", quietly = TRUE)) {
    return(rstudioapi::getSourceEditorContext()$path)
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(sub("--file=", "", file_arg))
  }
  return(getwd())
}

main_final_r <- get_script_path()
if (!is.na(file.info(main_final_r)$isdir) && file.info(main_final_r)$isdir) {
  file_path <- main_final_r
} else {
  file_path <- dirname(main_final_r)
}

convex_hull_preprocess_path <- file.path(file_path, "convex_hull_preprocess.r")
return_start_target_path <- file.path(file_path, "return_start_target_locations.R")
hexagon_path <- file.path(file_path, "hexagons")

# Source files
source(convex_hull_preprocess_path)
source(return_start_target_path)

# Parameters
hexa_name <- "Hexa127"
hexa_style <- "star"
with_center <- TRUE
target_location <- "nearcenter"
is_convex_hull <- TRUE

# Create from-to relationship using vectorized R operations (100x faster than loops)
adj_matrix_to_binary_neighborhood_list <- function(adj_matrix) {
  edges <- which(adj_matrix == 1, arr.ind = TRUE)
  data.frame(
    from = as.integer(rownames(adj_matrix)[edges[, 1]]),
    to = as.integer(colnames(adj_matrix)[edges[, 2]])
  )
}

# Read spatial data
hexa <- read_sf(dsn = hexagon_path, layer = hexa_name, int64_as_string = TRUE, stringsAsFactors = TRUE)
ch <- if (is_convex_hull) "convex_hull" else "noch"

# Extract topology and binary neighborhood matrix
nb <- spdep::poly2nb(hexa, queen = TRUE)
B <- spdep::nb2mat(nb, style = "B")
rownames(B) <- 1:nrow(B)
colnames(B) <- 1:ncol(B)

# Extract centroids and coordinates
centroids <- st_centroid(hexa)
hexa_coordinates <- st_coordinates(centroids)
rownames(hexa_coordinates) <- 1:nrow(hexa_coordinates)
colnames(hexa_coordinates) <- c("X", "Y")

# Plot neighbors
png("Hexagonal_Neighbors.png", width = 800, height = 800)
plot(hexa$geometry, col = "palegreen3", border = grey(0.9), axes = TRUE)
plot(nb, coords = hexa_coordinates, pch = 19, cex = 0.6, col = "maroon", add = TRUE)
title("Hexagonal Neighbors")
dev.off()

# Determine start and base target nodes
sp <- 1
start_target_nodes <- return_hexagon_locations(
  hexa_name = hexa_name, 
  hexa_style = hexa_style,
  with_center = with_center, 
  target_location = target_location
)

start_point <- start_target_nodes[sp]
base_targets_original <- start_target_nodes[-sp]

# Find spacing (step distance)
start_coords <- hexa_coordinates[start_point, , drop = FALSE]
other_coords <- hexa_coordinates[-start_point, , drop = FALSE]
step_distance <- min(sqrt((other_coords[, "X"] - start_coords[1, "X"])^2 + (other_coords[, "Y"] - start_coords[1, "Y"])^2))

# Dynamically move base targets 1 step closer to the starting point
base_targets <- c()
for (tp in base_targets_original) {
  tp_coords <- hexa_coordinates[tp, , drop = FALSE]
  dir_x <- tp_coords[1, "X"] - start_coords[1, "X"]
  dir_y <- tp_coords[1, "Y"] - start_coords[1, "Y"]
  len <- sqrt(dir_x^2 + dir_y^2)
  if (len > 0) {
    dir_x <- dir_x / len
    dir_y <- dir_y / len
    
    proj_x <- tp_coords[1, "X"] - 1 * step_distance * dir_x
    proj_y <- tp_coords[1, "Y"] - 1 * step_distance * dir_y
    
    all_dists <- sqrt((hexa_coordinates[, "X"] - proj_x)^2 + (hexa_coordinates[, "Y"] - proj_y)^2)
    closest_node <- which.min(all_dists)
    base_targets <- c(base_targets, closest_node)
  }
}

# Dynamically find 3 more targets on the same path, giving two spaces for each (from the moved base targets)
extended_targets <- c()
for (tp in base_targets) {
  tp_coords <- hexa_coordinates[tp, , drop = FALSE]
  dir_x <- tp_coords[1, "X"] - start_coords[1, "X"]
  dir_y <- tp_coords[1, "Y"] - start_coords[1, "Y"]
  len <- sqrt(dir_x^2 + dir_y^2)
  if (len > 0) {
    dir_x <- dir_x / len
    dir_y <- dir_y / len
    
    # Project 2 steps further
    proj_x <- tp_coords[1, "X"] + 2 * step_distance * dir_x
    proj_y <- tp_coords[1, "Y"] + 2 * step_distance * dir_y
    
    # Find closest node in all coordinates
    all_dists <- sqrt((hexa_coordinates[, "X"] - proj_x)^2 + (hexa_coordinates[, "Y"] - proj_y)^2)
    closest_node <- which.min(all_dists)
    
    # Only append if it's different from the base target (boundary check)
    if (closest_node != tp) {
      extended_targets <- c(extended_targets, closest_node)
    }
  }
}

# Dynamically find the left neighbors of the furthest targets (extended_targets)
left_targets <- c()
for (tp in extended_targets) {
  neighbors <- which(B[tp, ] == 1)
  neighbor_coords <- hexa_coordinates[neighbors, , drop = FALSE]
  dX <- neighbor_coords[, "X"] - hexa_coordinates[tp, "X"]
  dY <- neighbor_coords[, "Y"] - hexa_coordinates[tp, "Y"]
  
  # Filter to neighbors strictly to the left
  left_indices <- which(dX < -0.1)
  if (length(left_indices) > 0) {
    best_idx <- left_indices[which.min(abs(dY[left_indices]))]
    left_targets <- c(left_targets, neighbors[best_idx])
  }
}

# Wrap MIP solving and plotting logic in a reusable function
run_optimization <- function(start_point, target_points, problem_label) {
  all_required_nodes <- unique(c(start_point, target_points))
  
  if (is_convex_hull) {
    preprocessed_surface <- convex_hull_preprocess(
      hexa_coordinates = hexa_coordinates,
      all_required_nodes = all_required_nodes,
      hexa = hexa,
      nb = nb
    )
  } else {
    preprocessed_surface <- hexa_coordinates
  }
  preprocessed_all_nodes <- rownames(preprocessed_surface)
  
  # Subset binary neighborhood to preprocessed surface
  B_preproccessed <- B[preprocessed_all_nodes, preprocessed_all_nodes]
  all_neighbors <- adj_matrix_to_binary_neighborhood_list(B_preproccessed)
  
  # Optimization: Filter out any edges going back to the start point
  all_neighbors <- subset(all_neighbors, to != start_point)
  
  # Add cost column
  all_neighbors_with_cost <- all_neighbors
  all_neighbors_with_cost$cost <- 1.0
  
  total_vars <- nrow(all_neighbors_with_cost)
  total_num_target_points <- length(target_points)
  
  cat("\n==================================================\n")
  cat(" RUNNING:", problem_label, "\n")
  cat("==================================================\n")
  cat("Start Node:", start_point, "\n")
  cat("Target Nodes:", paste(target_points, collapse = ", "), "\n")
  cat("Total Variables (edges):", total_vars, "\n")
  cat("Total Targets:", total_num_target_points, "\n")
  
  # Create MIP Model
  mdsp_model <- ompr::MIPModel() %>%
    ompr::add_variable(d[i], i = 1:total_vars, type = "binary") %>%
    ompr::add_variable(y[k, i], k = 1:total_num_target_points, i = 1:total_vars, type = "continuous", lb = 0, ub = 1) %>%
    ompr::set_objective(ompr::sum_expr(d[i] * all_neighbors_with_cost$cost[i], i = 1:total_vars), "min") %>%
    ompr::add_constraint(y[k, i] <= d[i], k = 1:total_num_target_points, i = 1:total_vars)
  
  # Flow conservation constraints (Node balance) for each commodity k
  for (k in 1:total_num_target_points) {
    target_node <- target_points[k]
    
    # Start node: flow out must be 1 for this commodity
    start_out_indices <- which(all_neighbors$from == start_point)
    mdsp_model <- ompr::add_constraint(
      mdsp_model, 
      ompr::sum_expr(y[k, i], i = start_out_indices) == 1
    )
    
    # Other nodes flow conservation
    other_nodes <- setdiff(preprocessed_all_nodes, start_point)
    for (node in other_nodes) {
      node_val <- as.integer(node)
      flow_in_indices <- which(all_neighbors$to == node_val)
      flow_out_indices <- which(all_neighbors$from == node_val)
      
      node_demand <- if (node_val == target_node) 1 else 0
      
      if (length(flow_in_indices) > 0 && length(flow_out_indices) > 0) {
        mdsp_model <- ompr::add_constraint(
          mdsp_model, 
          ompr::sum_expr(y[k, i], i = flow_in_indices) - ompr::sum_expr(y[k, j], j = flow_out_indices) == node_demand
        )
      } else if (length(flow_in_indices) > 0) {
        mdsp_model <- ompr::add_constraint(
          mdsp_model, 
          ompr::sum_expr(y[k, i], i = flow_in_indices) == node_demand
        )
      } else if (length(flow_out_indices) > 0) {
        mdsp_model <- ompr::add_constraint(
          mdsp_model, 
          -ompr::sum_expr(y[k, j], j = flow_out_indices) == node_demand
        )
      }
    }
  }
  
  # Solve model
  start_time <- Sys.time()
  result <- ompr::solve_model(mdsp_model, ompr.roi::with_ROI(solver = "glpk", verbose = TRUE))
  end_time <- Sys.time()
  solution_time <- round(as.numeric(end_time - start_time, units = "secs"), 2)
  
  # Extract solution
  solution <- ompr::get_solution(result, d[i])
  solution <- subset(solution, value > 0)
  
  cat("\nOptimal cost:", result$objective, "\n")
  solution_path <- all_neighbors_with_cost[solution$i, ]
  print(solution_path)
  
  # Plotting the solution
  png(paste0(gsub(" ", "_", problem_label), ".png"), width = 800, height = 800)
  plot(hexa$geometry, col = "palegreen3", border = grey(0.9), axes = TRUE)
  
  point_colors <- rep("black", nrow(hexa_coordinates))
  point_colors[unlist(target_points)] <- "red"
  point_colors[unlist(start_point)] <- "green"
  
  node_sizes <- rep(0.8, nrow(hexa_coordinates))
  node_sizes[c(unlist(target_points), unlist(start_point))] <- 2.0
  
  # Plot nodes
  points(hexa_coordinates, pch = 19, cex = node_sizes, col = point_colors)
  
  # Draw arrows for active edges in solution path
  for (node in 1:nrow(solution_path)) {
    p1 <- as.numeric(solution_path[node, "from"])
    p2 <- as.numeric(solution_path[node, "to"])
    x0 <- hexa_coordinates[p1, "X"]
    y0 <- hexa_coordinates[p1, "Y"]
    x1 <- hexa_coordinates[p2, "X"]
    y1 <- hexa_coordinates[p2, "Y"]
    arrows(x0 = x0, y0 = y0, x1 = x1, y1 = y1, col = "blue", lwd = 2.5, angle = 25, length = 0.15)
  }
  
  plot_title <- paste0(problem_label, " - Cost: ", round(sum(solution_path$cost), 4), " , Solve: ", solution_time, "s")
  title(plot_title)
  dev.off()
  
  return(list(cost = result$objective, paths = solution_path, time = solution_time))
}

# Run the 3-target problem
res_3t <- run_optimization(start_point, base_targets, "3-Target Problem")

# Run the 6-target problem (base targets + extended targets)
res_6t <- run_optimization(start_point, c(base_targets, extended_targets), "6-Target Problem")

# Run the 9-target problem (base targets + extended targets + left targets)
res_9t <- run_optimization(start_point, c(base_targets, extended_targets, left_targets), "9-Target Problem")


