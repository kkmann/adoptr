#' Boundary designs
#'
#' Compute boundary designs. These can be used for optimization and are
#' implemented by default in \link{minimize}.
#'
#' @param objective An object of class \code{\link{UnconditionalScore-class}} to be minimized
#' @param subject_to A \code{\link{ConstraintsCollection-class}} defining the constraints
#' @param initial_design The initial design
#' @param ... optimal arguments
#'
#' @rdname boundary-designs
#' @export
setGeneric("get_lower_boundary_design",
           function(objective, subject_to, initial_design, ...) standardGeneric("get_lower_boundary_design"))



#' @rdname boundary-designs
#' @export
setGeneric("get_upper_boundary_design",
           function(objective, subject_to, initial_design, ...) standardGeneric("get_upper_boundary_design"))



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "OneStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              OneStageDesign(5, 0.0)
          })


#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "GroupSequentialDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

              os_design <- minimize(
                  objective      = objective,
                  subject_to     = cnstrs,
                  initial_design = OneStageDesign(200, 2.0)
              )$design


              lb_design <- GroupSequentialDesign(5,
                                                 -1,
                                                 os_design@c1f + .05,
                                                 2,
                                                 rep(-2, length(initial_design@c2_pivots)))

              return(lb_design)
          })



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "TwoStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

              os_design <- minimize(
                  objective      = objective,
                  subject_to     = cnstrs,
                  initial_design = OneStageDesign(200, 2.0)
              )$design


              lb_design <- TwoStageDesign(5,
                                          -2,
                                          os_design@c1f + .05,
                                          rep(2, length(initial_design@c2_pivots)),
                                          rep(-2, length(initial_design@c2_pivots)))

              return(lb_design)
          })






#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "OneStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              OneStageDesign(10000, 5.0)
          })


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "GroupSequentialDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

              os_design <- minimize(
                  objective      = objective,
                  subject_to     = cnstrs,
                  initial_design = OneStageDesign(200, 2.0)
              )$design


              ub_design <- GroupSequentialDesign(os_design@n1,
                                                 os_design@c1f - .05,
                                                 5,
                                                 2*os_design@n1,
                                                 rep(5, length(initial_design@c2_pivots)))

              return(ub_design)
          })


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "TwoStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

              os_design <- minimize(
                  objective      = objective,
                  subject_to     = cnstrs,
                  initial_design = OneStageDesign(200, 2.0)
              )$design


              ub_design <- TwoStageDesign(os_design@n1,
                                          os_design@c1f - .05,
                                          5,
                                          rep(2*os_design@n1, length(initial_design@c2_pivots)),
                                          rep(5, length(initial_design@c2_pivots)))

              return(ub_design)
          })
