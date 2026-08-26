#' Plot Persistence Diagram
#'
#' @param df Dataframe from plot_persistence.
#' @return A ggplot2 object representing the persistence diagram.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_minimal labs scale_x_continuous scale_y_continuous coord_fixed theme .data
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_filtration(points, method = "VR", eps_max = 1.2)
#' res <- boundary_info(filtration)
#' pairs <- extract_persistence_pairs(filtration, res$last_1, res$pivot_owner)
#' plot_persistence(pairs)
plot_persistence <- function(df) {
  # add liitle padding for infinite death or it will cut the point
  finite_deaths <- df$death[is.finite(df$death)]
  death_max <- max(finite_deaths, na.rm = TRUE)
  death_pad <- 0.2 * death_max
  death_plot <- ifelse(is.infinite(df$death), death_max + death_pad, df$death)
  birth_min <- min(df$birth, na.rm = TRUE)
  x_lower <- min(0, birth_min - 0.05 * death_max)

  # "H0"/"H1"/"H2"... as a discrete factor, not the raw integer dim column,
  # so each dimension gets its own qualitative colour and legend entry
  # instead of shading indistinguishably along a continuous gradient.
  df$dim_label <- factor(paste0("H", df$dim), levels = paste0("H", sort(unique(df$dim))))

  plot <- ggplot(df, aes(x = .data$birth, y = death_plot, color = .data$dim_label)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = "Persistence Diagram", x = "Birth", y = "Death", color = "Dimension") +
    scale_x_continuous(limits = c(x_lower, death_max + death_pad), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, death_max + death_pad + 0.1), expand = c(0,0)) +
    coord_fixed() + # same x, y scale
    theme(legend.position = "right")

  return(plot)
}

#' Plot a Persistence Landscape
#'
#' @param landscape_df Data frame from \code{persistence_landscape()}, with
#'   columns \code{t}, \code{k}, \code{value}.
#' @return A ggplot2 object with one line per landscape level
#'   \eqn{\lambda(k, \cdot)}.
#'
#' @importFrom ggplot2 ggplot aes geom_line theme_minimal labs theme .data
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_filtration(points, method = "VR", eps_max = 1.2)
#' res <- boundary_info(filtration)
#' pairs <- extract_persistence_pairs(filtration, res$last_1, res$pivot_owner)
#' landscape <- persistence_landscape(pairs, dimension = 0)
#' plot_landscape(landscape)
plot_landscape <- function(landscape_df) {
  plot <- ggplot(landscape_df, aes(x = .data$t, y = .data$value,
                                    color = factor(.data$k), group = .data$k)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Persistence Landscape", x = "t", y = "Landscape value", color = "k") +
    theme(legend.position = "right")

  return(plot)
}
