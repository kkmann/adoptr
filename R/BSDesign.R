setClass("BSDesign",
        representation(
                n1         = "numeric",
                early_stopping_for_futility = "numeric",
                early_stopping_for_efficacy = "numeric",
                weights_n2 = "numeric",
                weights_c2 = "numeric",
                knots      = "numeric",
                degree     = "integer"
        ),
        contains = "Design"
)

BSDesign <- function(
        n1,
        early_stopping_for_futility,
        early_stopping_for_efficacy,
        knots,
        degree = 3L,
        weights_n2,
        weights_c2
) {

        new("BSDesign",
                n1 = n1,
                early_stopping_for_futility = early_stopping_for_futility,
                early_stopping_for_efficacy = early_stopping_for_efficacy,
                knots = knots,
                degree = degree,
                weights_n2 = weights_n2,
                weights_c2 = weights_c2
        )

}

BSDesign_ <- function(
        tunable_parameters, knots, degree = 3L
) {
        if (!(length(tunable_parameters) == 3 + 2*(length(knots))))
                stop("parameter length does not fit")
        # n weights
        names(tunable_parameters) <- NULL
        n_weights <- length(knots)
        new("BSDesign", # need to
            n1 = tunable_parameters[1],
            early_stopping_for_futility = tunable_parameters[2],
            early_stopping_for_efficacy = tunable_parameters[3],
            knots = knots,
            degree = degree,
            weights_n2 = tunable_parameters[4:(3 + n_weights)],
            weights_c2 = tunable_parameters[(4 + n_weights):length(tunable_parameters)]
        )

}

setMethod("n1",
          signature("BSDesign"),
          function(d, ...) d@n1
)

setMethod("n2",
        signature("BSDesign", "numeric"),
        function(d, z1, ...) {
                # tmp = splines::bs(
                #         x = z1,
                #         knots = d@knots,
                #         degree = d@degree,
                #         Boundary.knots = range(d@knots)
                # ) %*% d@weights_n2
                tmp <- approx(d@knots, d@weights_n2, xout = z1, yleft = 0, yright = 0, rule = 2, method = "linear")$y
                tmp <- ifelse(z1 <= d@early_stopping_for_futility | z1 >= d@early_stopping_for_efficacy, 0, tmp)
                return(tmp)
        }
)

setMethod("c2",
        signature("BSDesign", "numeric"),
        function(d, z1, ...) {
                # tmp = splines::bs(
                #         x = z1,
                #         knots = d@knots,
                #         degree = d@degree,
                #         Boundary.knots = range(d@knots)
                # ) %*% d@weights_c2
                tmp <- approx(d@knots, d@weights_c2, xout = z1, yleft = Inf, yright = -Inf, rule = 2, method = "linear")$y
                tmp <- ifelse(z1 <= d@early_stopping_for_futility, Inf, tmp)
                tmp <- ifelse(z1 >= d@early_stopping_for_efficacy, -Inf, tmp)
                return(tmp)
        }
)

setMethod("conditional_power",
          signature("BSDesign", "numeric", "numeric"),
          function(d, z1, delta, ...) {
                  n2 <- n2(d, z1)
                  c2 <- c2(d, z1)
                  return(1 - pnorm(c2 - sqrt(n2)*delta))
          }
)

setMethod("get_tunable_parameters",
          signature("BSDesign"),
          function(design, simplify = FALSE, ...) {
                  res <- list(
                          n1 = n1(design),
                          early_stopping_for_futility = design@early_stopping_for_futility,
                          early_stopping_for_efficacy = design@early_stopping_for_efficacy,
                          weights_n2 = design@weights_n2,
                          weights_c2 = design@weights_c2
                        )
                  if (simplify) {
                          return(as.numeric(do.call(c, res)))
                  } else {
                          return(res)
                  }
          }
)
