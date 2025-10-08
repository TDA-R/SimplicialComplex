#' Plot Persistence Diagram
#'
#' @param df Dataframe from plot_persistence
#'
#' @return A ggplot2 object representing the persistence diagram
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_minimal labs scale_x_continuous scale_y_continuous coord_fixed theme
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_vr_filtration(points, eps_max=1.2)
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

  plot <- ggplot(df, aes(x = birth, y = death_plot, color = dim)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = "Persistence Diagram", x = "Birth", y = "Death") +
    scale_x_continuous(limits = c(x_lower, death_max + death_pad), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, death_max + death_pad + 0.1), expand = c(0,0)) +
    coord_fixed() + # same x, y scale
    theme(legend.position="none")

  return(plot)
}
