library(hexSticker)
library(adoptr)
library(tidyverse)
library(cowplot)

order <- 7L

design <- TwoStageDesign(
    n1 = 25,
    c1f = 0,
    c1e = 2,
    n2 = 40,
    c2 = 2,
    order = order
)

null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)
datadist    <- Normal(two_armed = FALSE)

ess  <- ExpectedSampleSize(datadist, alternative)
pow  <- Power(datadist, alternative)
toer <- Power(datadist, null)

optimal_design <- minimize(
    ess,
    subject_to(
        pow  >= 0.8,
        toer <= .05
    ),
    design,
    opts = list(
        algorithm   = "NLOPT_LN_COBYLA",
        xtol_rel    = 1e-5,
        maxeval     = 50000, # super precise but super slow
        print_level = 1
    )
)

optimal_design@rounded <- TRUE


dat   <- data.frame(t=seq(0, 2*pi, by=0.1) )
xhrt  <- function(t) 16*sin(t)^3
yhrt  <- function(t) 13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t)
dat$y <- yhrt(dat$t)
dat$x <- xhrt(dat$t)
with(dat, plot(x,y, type = "l"))


df <- tibble(
    x1 = seq(-1, 3, by = .01),
    n  = adoptr::n(optimal_design$design, x1),
    c2 = c2(optimal_design$design, x1)
) %>%
    mutate(
        tmp1 = abs(x1 - optimal_design$design@c1f),
        tmp2 = abs(x1 - optimal_design$design@c1e),
        n = ifelse(tmp1 == min(tmp1) | tmp2 == min(tmp2), NA, n)
    ) %>%
    gather(variable, value, n, c2) %>%
    filter(is.finite(value) | is.na(value))

p1 <- dat %>%
    ggplot(aes(x, y)) +
        geom_polygon(aes(x, y), data = dat, fill = "white") +
        xlim(c(-100, 90)) +
        ylim(c(-70, 20)) +
        theme_void() +
        theme_transparent()

p2 <- df %>% ggplot(aes(x = x1, y = value)) +
    geom_line(color = "white", size = .5) +
    facet_wrap(~factor(variable, levels = c("n", "c2")), scales = "free_y") +
    theme_void() +
    ylim(c(0, NA)) +
    theme_transparent() +
    coord_cartesian() +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "white", fill = NA, size = 1)
    )

pp <- cowplot::plot_grid(p1, p2, nrow = 2)

sticker(pp, package="adoptr", p_size = 120, s_x=1, s_y=1.125, s_width=1.5, s_height=1.15,
        filename="man/figures/logo.png", h_fill = "#DE29F5", h_color = "#8B1999", dpi = 1200)
sticker(pp, package="adoptr", p_size = 10, s_x=1, s_y=1.125, s_width=1.5, s_height=1.15,
        filename="man/figures/logo.pdf", h_fill = "#DE29F5", h_color = "#8B1999")
