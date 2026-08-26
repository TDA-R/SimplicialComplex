# ---------------------------------------------------------------------------
# Visualizing the matching behind a Wasserstein/bottleneck distance
#
# plot_matching() draws both persistence diagrams on one set of axes and
# connects every matched pair with a dashed line: a real-to-real match gets
# a line straight between the two points, and a point matched to the
# diagonal (i.e. effectively unmatched, the standard device described in
# R/Distance.R) gets a dashed line down to its own perpendicular projection
# onto birth = death. The underlying correspondence comes from
# diagram_matching() (R/Distance.R), which reduces to the exact same
# assignment problem wasserstein_distance()/bottleneck_distance() solve.
# ---------------------------------------------------------------------------

#' Plot the optimal matching between two persistence diagrams
#'
#' Visualizes the point correspondence that realizes the Wasserstein or
#' bottleneck distance between two persistence diagrams in a given
#' homological dimension: both diagrams' points, overlaid on the same axes,
#' joined by dashed lines to their matched partner - either a point of the
#' other diagram, or (for a point left unmatched) its own projection onto
#' the diagonal \code{birth = death}.
#'
#' @param df1,df2 Persistence diagram data frames with columns \code{dim},
#'   \code{birth}, \code{death}.
#' @param dimension Homological dimension to compare.
#' @param distance Which distance's optimal matching to visualize,
#'   \code{"wasserstein"} (default) or \code{"bottleneck"}.
#' @param p Wasserstein order, used only when \code{distance = "wasserstein"}.
#' @param ground Ground metric, \code{"L2"} or \code{"Linf"}. Defaults to
#'   each distance's own default (\code{"L2"} for Wasserstein, \code{"Linf"}
#'   for bottleneck) when left \code{NULL} - see
#'   \code{\link{wasserstein_distance}}/\code{\link{bottleneck_distance}}.
#' @param labels Legend labels identifying \code{df1} and \code{df2},
#'   e.g. \code{c("Clean", "Noisy")}.
#' @return A ggplot2 object; the title reports the resulting distance.
#'
#' @details
#' Essential (death = \code{Inf}) points are capped exactly the way
#' \code{\link{wasserstein_distance}} does (see its Details) so they take
#' part in the matching instead of being dropped, and are drawn as triangles
#' at their capped height - marked "essential" in the legend - so they stay
#' visually distinguishable from genuine finite points that happen to reach
#' that height. Matches to the diagonal are always drawn to the point's
#' perpendicular projection \code{((birth+death)/2, (birth+death)/2)}
#' regardless of \code{ground}, since that is the standard, readable way to
#' depict an unmatched point whichever ground metric produced the matching.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_abline theme_minimal labs scale_shape_manual coord_fixed theme .data
#' @export
#' @examples
#' df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 2))
#' df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 3))
#' plot_matching(df1, df2, dimension = 0, distance = "wasserstein", p = 2)
plot_matching <- function(df1, df2, dimension, distance = c("wasserstein", "bottleneck"),
                           p = 2, ground = NULL, labels = c("Diagram 1", "Diagram 2")) {
  distance <- match.arg(distance)
  res <- diagram_matching(df1, df2, dimension, distance = distance, p = p, ground = ground)

  proj <- function(b, d) (b + d) / 2 # perpendicular projection onto birth = death

  seg_list <- list()

  rr <- res$matches[res$matches$type == "real-real", , drop = FALSE]
  if (nrow(rr) > 0) {
    seg_list[[length(seg_list) + 1]] <- data.frame(
      x = rr$x_birth, y = rr$x_death, xend = rr$y_birth, yend = rr$y_death
    )
  }
  xd <- res$matches[res$matches$type == "x-diagonal", , drop = FALSE]
  if (nrow(xd) > 0) {
    pr <- proj(xd$x_birth, xd$x_death)
    seg_list[[length(seg_list) + 1]] <- data.frame(x = xd$x_birth, y = xd$x_death, xend = pr, yend = pr)
  }
  yd <- res$matches[res$matches$type == "y-diagonal", , drop = FALSE]
  if (nrow(yd) > 0) {
    pr <- proj(yd$y_birth, yd$y_death)
    seg_list[[length(seg_list) + 1]] <- data.frame(x = yd$y_birth, y = yd$y_death, xend = pr, yend = pr)
  }
  seg_df <- if (length(seg_list) > 0) {
    do.call(rbind, seg_list)
  } else {
    data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0))
  }

  pts_df <- rbind(
    data.frame(birth = res$X[, "birth"], death = res$X[, "death"],
               essential = res$essential_X, diagram = rep(labels[1], nrow(res$X))),
    data.frame(birth = res$Y[, "birth"], death = res$Y[, "death"],
               essential = res$essential_Y, diagram = rep(labels[2], nrow(res$Y)))
  )
  pts_df$point_type <- ifelse(pts_df$essential, "essential", "finite")
  pts_df$diagram <- factor(pts_df$diagram, levels = labels)

  vmax <- max(pts_df$birth, pts_df$death)
  pad <- max(0.08 * vmax, 1e-6)

  distance_label <- paste0(toupper(substring(distance, 1, 1)), substring(distance, 2))

  ggplot() +
    geom_segment(
      data = seg_df, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      linetype = "dashed", color = "grey50", alpha = 0.8
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(
      data = pts_df,
      aes(x = .data$birth, y = .data$death, color = .data$diagram, shape = .data$point_type),
      size = 3, alpha = 0.85
    ) +
    scale_shape_manual(values = c(finite = 16, essential = 17), name = NULL) +
    theme_minimal() +
    labs(
      title = sprintf("%s matching (H%d), distance = %.4g", distance_label, dimension, res$distance),
      x = "Birth", y = "Death", color = NULL
    ) +
    coord_fixed(xlim = c(min(0, -pad), vmax + pad), ylim = c(min(0, -pad), vmax + pad)) +
    theme(legend.position = "right")
}
