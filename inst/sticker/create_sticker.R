# Creates the bsimms hex sticker (inst/sticker/logo.svg).
#
# The figure is an isospace plot (plot_isospace()) of five simulated mixtures
# and 3 sources, drawn in front of the density of each source's proportion:
# the mixtures' proportions are Dirichlet draws around a global proportion
# vector, and the curves are the marginal densities of that Dirichlet.
#
# Run from the package root with the development version of bsimms loaded
# (e.g. after `devtools::load_all()`), and with the `hexSticker` and `svglite`
# packages installed.

library(bsimms)
library(ggplot2)
library(hexSticker)

seed <- 2
p_global <- c(0.2, 0.3, 0.5) # global source proportions
kappa <- 14 # Dirichlet concentration: higher = less noise around p_global
n_mixtures <- 5
cols <- c(a = "#DC267F", b = "#FFB000", c = "#648FFF") # IBM colour-blind safe
out <- "inst/sticker/logo.svg"

set.seed(seed)

# sources (means/SDs, chosen for an irregular triangle that fits well the
# sticker)
src <- data.frame(
  Source = c("a", "b", "c"),
  d13C_mean = c(0.5, 0.9, 0),
  d13C_sd = c(0.24, 0.16, 0.10),
  d15N_mean = c(0, 0.6, 0.5),
  d15N_sd = c(0.10, 0.22, 0.13)
)
tdf <- data.frame(
  Source = src$Source,
  d13C_mean = 0,
  d13C_sd = 0.01,
  d15N_mean = 0,
  d15N_sd = 0.01
)

# mixtures: each sample's proportions vary randomly around the global ones
g <- matrix(
  rgamma(n_mixtures * 3, shape = kappa * p_global),
  ncol = 3,
  byrow = TRUE
)
w <- g / rowSums(g)
mix <- data.frame(
  d13C = as.vector(w %*% src$d13C_mean) + rnorm(n_mixtures, 0, 0.02),
  d15N = as.vector(w %*% src$d15N_mean) + rnorm(n_mixtures, 0, 0.02)
)

# isospace plot, stripped down to the bare geoms
p <- plot_isospace(mix, src, tdf, c("d13C", "d15N"), source_means_sds = TRUE) +
  theme_void() +
  theme(legend.position = "none")
mix_layer <- 3
p$layers[[mix_layer]]$aes_params[c("fill", "colour", "size", "alpha")] <-
  list("white", "white", 1.2, 1)
for (i in 1:2) {
  p$layers[[i]]$aes_params$linewidth <- 0.6
}
p$layers[[4]]$aes_params$size <- 2.4

# background: densities of the source proportions, as their own plot (with its
# own axes, so independent of the isospace scales) placed behind the isospace
xs <- seq(0, 1, length.out = 300)
dens <- do.call(
  rbind,
  lapply(1:3, function(k) {
    data.frame(
      x = xs,
      y = dbeta(xs, kappa * p_global[k], kappa * (1 - p_global[k])),
      source = c("a", "b", "c")[k]
    )
  })
)
dens$y <- dens$y / max(dens$y)
dens <- dens[dens$y > 0.01, ]
bg_plot <- ggplot(dens, aes(x, y, colour = source)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.08, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none") +
  theme_transparent()
inset <- grid::grobTree(
  grid::editGrob(
    ggplotGrob(bg_plot),
    vp = grid::viewport(x = 0.56, y = 0.58, width = 0.66, height = 0.68)
  )
)
bg_layer <- annotation_custom(
  inset,
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf
)
p$layers <- c(list(bg_layer), p$layers)
p <- p +
  scale_colour_manual(values = cols) +
  coord_cartesian(xlim = c(-0.2, 1.15), ylim = c(-0.2, 0.9)) +
  theme_transparent()

sticker(
  p,
  package = "bsimms",
  p_size = 8,
  p_y = 0.6,
  p_color = "white",
  s_x = 1,
  s_y = 1.17,
  s_width = 1.45,
  s_height = 1.1,
  h_fill = "#000000",
  h_color = "#C8C8C8",
  filename = out
)
