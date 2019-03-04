#' Boundary designs
#'
#' Compute boundary designs. These can be used for optimization and are
#' implemented by default in \link{minimize}.
#'
#' @param initial_design The initial design
#' @param ... optimal arguments
#'
#' @rdname boundary-designs
#' @export
setGeneric("get_lower_boundary_design",
           function(initial_design, ...) standardGeneric("get_lower_boundary_design"))



#' @rdname boundary-designs
#' @export
setGeneric("get_upper_boundary_design",
           function(initial_design, ...) standardGeneric("get_upper_boundary_design"))



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("OneStageDesign"),
          function(initial_design, ...) {
              OneStageDesign(5, 0.0)
          })


#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("GroupSequentialDesign"),
          function(initial_design, ...) {
              GroupSequentialDesign(5,
                                    -1,
                                    1.5,
                                    2,
                                    rep(-2, length(initial_design@c2_pivots)))
})



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("TwoStageDesign"),
          function(initial_design, ...) {
              TwoStageDesign(5,
                             -2,
                             1.5,
                             rep(2, length(initial_design@c2_pivots)),
                             rep(-2, length(initial_design@c2_pivots)))
})






#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("OneStageDesign"),
          function(initial_design, ...) {
              OneStageDesign(5 * initial_design@n1, 5.0)
})


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("GroupSequentialDesign"),
          function(initial_design, ...) {
              GroupSequentialDesign(5 * initial_design@n1,
                                    2.0,
                                    5.0,
                                    5 * initial_design@n2_pivots,
                                    rep(5.0, length(initial_design@c2_pivots)))
})


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("TwoStageDesign"),
          function(initial_design, ...) {
              TwoStageDesign(5 * initial_design@n1,
                             2.0,
                             5.0,
                             5 * initial_design@n2_pivots,
                             rep(5.0, length(initial_design@c2_pivots)))

})
