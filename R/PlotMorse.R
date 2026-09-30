#' Plot the graph reconstructed by \code{morse_recon}
#'
#' Draws the reconstructed graph \eqn{\hat G} (the union of 1-unstable manifolds of the surviving high-persistence critical edges)
#' over the input vertices, optionally coloured by the input density \code{rho}.
#'
#' @param mr A \code{"morse_recon"} object returned by \code{\link{morse_recon}}.
#' @param show_density Logical; if \code{TRUE} (default) colour vertices by \code{mr$rho} using a continuous viridis scale.
#' @param point_size,edge_size Point/line sizes passed to ggplot2.
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_color_viridis_c theme_minimal labs coord_fixed theme .data
#' @export
plot_morse_recon <- function(mr, show_density = TRUE, point_size = 0.6, edge_size = 0.9) {
  stopifnot(inherits(mr, "morse_recon"))
  if (is.null(mr$coords)) {
    stop("morse_recon() was run without `vertex_coords`; cannot plot without vertex positions.")
  }

  coords <- mr$coords
  vdf <- data.frame(x = coords[, 1], y = coords[, 2], rho = mr$rho)

  ge <- mr$graph_edges
  if (nrow(ge) > 0) {
    edf <- data.frame(
      x = coords[ge$u, 1], y = coords[ge$u, 2],
      xend = coords[ge$v, 1], yend = coords[ge$v, 2]
    )
  } else {
    edf <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0))
  }

  p <- ggplot()
  if (show_density) {
    p <- p + geom_point(data = vdf, aes(x = .data$x, y = .data$y, color = .data$rho), size = point_size, alpha = 0.5) + scale_color_viridis_c(name = "density")
  } else {
    p <- p + geom_point(data = vdf, aes(x = .data$x, y = .data$y), size = point_size, color = "grey70", alpha = 0.5)
  }

  p + geom_segment(data = edf, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), color = "firebrick", linewidth = edge_size) +
    coord_fixed() + theme_minimal() +
    labs(title = "Discrete Morse graph reconstruction", subtitle = sprintf("delta = %.3g, %d reconstructed edges", mr$delta, nrow(ge)), x = NULL, y = NULL) +
    theme(legend.position = "right")
}

#' Plot a local patch of the triangulation with its discrete gradient field
#'
#' Figure zooms into a small neighbourhood of the mesh (as most DMT papers illustrate the
#' vector field, since drawing it over the whole domain is unreadable) and
#' draws the actual computed structure on top of the local triangulation:
#' \itemize{
#'   \item every triangulation edge with both endpoints in the window, in light grey, for context;
#'   \item the surviving critical edges of \code{mr$dmvf$critical_edges} that fall inside the window, as thick red segments.
#'     These are the "saddles" that did not get matched to a vertex;
#'   \item for every non-critical vertex \eqn{v} in the window whose DMVF match \code{mr$dmvf$parent[[v]]} is also in the window, a black arrow
#'     from \eqn{v} to the midpoint of its matched edge \eqn{(v, \mathrm{parent}(v))}, the vertex-edge V-path arrows themselves.
#' }
#'
#' @param mr A \code{"morse_recon"} object returned by \code{\link{morse_recon}}.
#' @param center Length-2 numeric, an \code{(x, y)} point in the same units as \code{mr$coords} to centre the window on.
#' @param radius Euclidean radius (same units as \code{mr$coords}) of the window around \code{center}.
#' @param vertex_size,arrow_size,critical_size Point size, arrow-head size, and line width for the critical edges, respectively.
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_fixed theme_void labs
#' @export
plot_morse_vpath <- function(mr, center, radius, vertex_size = 1.2, arrow_size = 0.12, critical_size = 1.3) {
  stopifnot(inherits(mr, "morse_recon"))
  if (is.null(mr$coords)) {
    stop("morse_recon() was run without `vertex_coords`; cannot plot without vertex positions.")
  }

  coords <- mr$coords
  dist_to_center <- sqrt((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2)
  win <- which(dist_to_center <= radius)
  if (length(win) == 0) {
    stop("No vertices found within `radius` of `center`; try a larger radius.")
  }
  win_set <- win

  # local triangulation edges (both endpoints in the window), for context
  edges1 <- Filter(function(x) length(x$simplex) == 2L, mr$filtration)
  e_mat <- do.call(rbind, lapply(edges1, `[[`, "simplex"))
  in_win <- e_mat[, 1] %in% win_set & e_mat[, 2] %in% win_set
  e_mat <- e_mat[in_win, , drop = FALSE]
  mesh_df <- data.frame(
    x = coords[e_mat[, 1], 1], y = coords[e_mat[, 1], 2],
    xend = coords[e_mat[, 2], 1], yend = coords[e_mat[, 2], 2]
  )

  # critical edges inside the window
  ce <- mr$dmvf$critical_edges
  ce_win <- ce[ce$u %in% win_set & ce$v %in% win_set, , drop = FALSE]
  crit_df <- data.frame(
    x = coords[ce_win$u, 1], y = coords[ce_win$u, 2],
    xend = coords[ce_win$v, 1], yend = coords[ce_win$v, 2]
  )

  # vertex-edge V-path arrows: v -> midpoint(v, parent(v)), for matched
  # vertices whose partner is also inside the window
  parent <- mr$dmvf$parent
  par_of_win <- parent[as.character(win_set)]
  has_par <- !is.na(par_of_win)
  par_in_win <- has_par & (as.integer(par_of_win) %in% win_set)
  arrow_v <- win_set[par_in_win]
  arrow_p <- as.integer(par_of_win[par_in_win])
  arrow_df <- data.frame(
    x = coords[arrow_v, 1], y = coords[arrow_v, 2],
    xend = (coords[arrow_v, 1] + coords[arrow_p, 1]) / 2,
    yend = (coords[arrow_v, 2] + coords[arrow_p, 2]) / 2
  )

  vdf <- data.frame(x = coords[win_set, 1], y = coords[win_set, 2])

  p <- ggplot() +
    geom_segment(data = mesh_df, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                 color = "grey75", linewidth = 0.3)

  if (nrow(arrow_df) > 0) {
    p <- p + geom_segment(
      data = arrow_df, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "black", linewidth = 0.45,
      arrow = grid::arrow(length = grid::unit(arrow_size, "cm"), type = "closed")
    )
  }

  if (nrow(crit_df) > 0) {
    p <- p + geom_segment(
      data = crit_df, aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "firebrick", linewidth = critical_size
    )
  }

  p +
    geom_point(data = vdf, aes(x = .data$x, y = .data$y), size = vertex_size, color = "grey40") +
    coord_fixed() +
    theme_void() +
    labs(
      title = "Local discrete gradient vector field",
      subtitle = sprintf("%d vertices, %d critical edges in window", length(win_set), nrow(ce_win))
    )
}

# interactive 3D visualisation (rgl)

#' Plot the density landscape and critical structure in interactive 3D
#'
#' A 3D companion to \code{\link{plot_morse_recon}}: renders the input density \code{rho} as an actual terrain surface over the triangulated
#' domain (using the triangles stored in \code{mr$filtration}), and overlays the discrete Morse critical structure on it with the \pkg{rgl} package,
#' producing a mouse-rotatable 3D scene rather than a flat projection:
#' \itemize{
#'   \item the root (critical) vertices of \code{mr$dmvf}, the surviving density peaks after persistence simplification - as labelled spheres ("M1", "M2", ...);
#'   \item the endpoints of \code{mr$dmvf$critical_edges}, the surviving "saddle-like" edges of Section 10.4.1 as smaller spheres;
#'   \item \code{mr$graph_edges} itself, i.e. the collected 1-unstable manifolds from \code{\link{collect_g}}, literally the vertex-to-
#'     saddle connections along the reconstructed ridge lines as raised line segments on the surface.
#' }
#'
#' @param mr A \code{"morse_recon"} object returned by \code{\link{morse_recon}} (must have been run with \code{vertex_coords}).
#' @param z_scale Vertical exaggeration applied to \code{mr$rho} before plotting. Default \code{NULL} auto-scales so the density relief spans
#'   about 35\% of the horizontal extent of \code{mr$coords}.
#' @param show_critical_edges Logical; if \code{TRUE} (default) also mark the endpoints of \code{mr$dmvf$critical_edges} (the saddle-like edges).
#' @param label_roots Logical; if \code{TRUE} (default) label root vertices "M1", "M2", ... in the scene.
#' @param surface_col,root_col,saddle_col,edge_col Colour for the density surface, the root markers, the saddle markers, and the \code{graph_edges} ridge lines, respectively.
#' @param point_radius,edge_lwd Marker radius, auto-scaled from the domain size when \code{NULL}) and line width for \code{graph_edges}.
#' @param window_size Length-2 \code{c(width, height)} in pixels for the new \pkg{rgl} device (only used when \code{new_window = TRUE}).
#' @param title_cex Character expansion for the \code{rgl::title3d()} title.
#'
#' @return A interactive 3D \pkg{rgl} scene.
#'
#' @section Display:
#' On macOS, opening a native \pkg{rgl} window requires XQuartz; Windows and Linux do not need it.
#' To avoid a native window altogether, render the scene as a WebGL widget in the RStudio Viewer or a browser:
#' \preformatted{
#' options(rgl.useNULL = TRUE)
#' plot_morse_landscape(mr)
#' rgl::rglwidget()
#' }
#'
#' @export
plot_morse_landscape <- function(
    mr, z_scale = NULL, show_critical_edges = TRUE, label_roots = TRUE, surface_col = c("#f7fbff", "#6baed6", "#08306b"), root_col = "red",
    saddle_col = "blue", edge_col = "black", point_radius = NULL, edge_lwd = 3, window_size = c(1400, 1000), title_cex = 1.4) {

  stopifnot(inherits(mr, "morse_recon"))

  if (is.null(mr$coords)) {
    stop("morse_recon() was run without `vertex_coords`; cannot plot without vertex positions.")
  }
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("plot_morse_landscape() needs the 'rgl' package; install it with install.packages(\"rgl\").")
  }

  coords <- mr$coords
  rho <- as.numeric(mr$rho)

  tris <- Filter(function(x) length(x$simplex) == 3L, mr$filtration)

  if (length(tris) == 0) {
    stop("mr$filtration has no triangles; plot_morse_landscape() needs morse_recon() to have been run on a 2-complex.")
  }

  tri_idx <- do.call(rbind, lapply(tris, `[[`, "simplex"))

  xy_span <- max(diff(range(coords[, 1])), diff(range(coords[, 2])))
  rho_span <- diff(range(rho))
  if (is.null(z_scale)) {
    z_scale <- if (rho_span > 0) 0.35 * xy_span / rho_span else 1
  }
  z <- rho * z_scale

  if (is.null(point_radius)) point_radius <- 0.01 * xy_span

  edge_offset <- 0.02 * max(diff(range(z)), 1e-8)

  # per-vertex colour ramp over density (Gouraud-shaded by rgl::triangles3d when `col` has one entry per vertex)
  if (length(surface_col) >= 2) {
    ramp <- grDevices::colorRampPalette(surface_col)(256)
    rr <- if (rho_span > 0) (rho - min(rho)) / rho_span else rep(0, length(rho))
    vcols <- ramp[pmin(256L, pmax(1L, round(rr * 255) + 1L))]
  } else {
    vcols <- rep(surface_col, length(rho))
  }

  flat <- as.vector(t(tri_idx)) # triangle-major order: v1,v2,v3, v1,v2,v3, ...

  rgl::triangles3d(coords[flat, 1], coords[flat, 2], z[flat], col = vcols[flat], specular = "black")

  roots <- mr$dmvf$roots

  if (length(roots) > 0) {
    rgl::spheres3d(
      coords[roots, 1], coords[roots, 2], z[roots] + edge_offset,
      radius = point_radius, color = root_col
    )
    if (label_roots) {
      rgl::text3d(
        coords[roots, 1], coords[roots, 2], z[roots] + edge_offset + 3 * point_radius,
        texts = paste0("M", seq_along(roots)), color = root_col, cex = 1.2
      )
    }
  }

  if (show_critical_edges) {

    ce <- mr$dmvf$critical_edges

    if (nrow(ce) > 0) {

      sv <- setdiff(unique(c(ce$u, ce$v)), roots)

      if (length(sv) > 0) {

        rgl::spheres3d(
          coords[sv, 1], coords[sv, 2], z[sv] + edge_offset,
          radius = 0.7 * point_radius, color = saddle_col
        )
      }
    }
  }

  ge <- mr$graph_edges

  if (nrow(ge) > 0) {

    idx_pairs <- as.vector(rbind(ge$u, ge$v))
    rgl::segments3d(coords[idx_pairs, 1], coords[idx_pairs, 2], z[idx_pairs] + edge_offset, color = edge_col, lwd = edge_lwd)
  }

  rgl::title3d(main = sprintf("Discrete Morse landscape (delta = %.3g)", mr$delta), col = "black", cex = title_cex)
  rgl::aspect3d(1, diff(range(coords[, 2])) / max(diff(range(coords[, 1])), 1e-8), 0.5)
  rgl::view3d(theta = -35, phi = 25, zoom = 0.8)

  invisible(rgl::cur3d())
}
