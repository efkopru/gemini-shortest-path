# Systemmatic tests configuration lookup for hexagon start/target locations.
# Supports Hexa061, 091, 127, 169, 217, 271 for various size experiments.

#' Return start and target node locations based on hexagon settings
#'
#' @param hexa_name Character name of the hexagon shape file (e.g. "Hexa061")
#' @param hexa_style Character style of configuration (e.g. "star")
#' @param with_center Logical indicating if center node is included
#' @param target_location Character indicating target proximity ("nearcenter" or "nearcorner")
#' @return Integer vector representing start (first element) and target nodes
return_hexagon_locations <- function(hexa_name, hexa_style, with_center, target_location) {
  # Define configurations in a structured lookup list
  configs <- list(
    "star" = list(
      "TRUE" = list(
        "nearcenter" = list(
          "Hexa061" = c(31, 27, 61, 5),
          "Hexa091" = c(46, 42, 84, 12),
          "Hexa127" = c(64, 60, 110, 22),
          "Hexa169" = c(85, 81, 139, 35),
          "Hexa217" = c(109, 105, 171, 51),
          "Hexa271" = c(136, 132, 206, 70)
        ),
        "nearcorner" = list(
          "Hexa061" = c(31, 27, 61, 5),
          "Hexa091" = c(46, 42, 84, 12),
          "Hexa127" = c(64, 59, 119, 14),
          "Hexa169" = c(85, 79, 160, 16),
          "Hexa217" = c(109, 102, 207, 18),
          "Hexa271" = c(136, 128, 260, 20)
        )
      )
    )
  )

  # Fetch the matching configuration
  style_cfg <- configs[[hexa_style]]
  center_cfg <- style_cfg[[as.character(with_center)]]
  loc_cfg <- center_cfg[[target_location]]
  start_target_nodes <- loc_cfg[[hexa_name]]

  if (is.null(start_target_nodes)) {
    warning(sprintf("No start/target configuration found for: name='%s', style='%s', with_center=%s, target_location='%s'",
                    hexa_name, hexa_style, with_center, target_location))
    return(integer(0))
  }

  return(start_target_nodes)
}

