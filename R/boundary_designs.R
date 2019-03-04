#' Boundary designs
#'
#' Compute boundary designs. These can be used for optimization and are
#' implemented by default in \link{minimize}.
#'
#' @param initial_design The initial design
#' @param n1 bound for n1
#' @param c1f bound for c1f
#' @param c1e bound for c1e
#' @param n2_pivots bound for n2_pivots
#' @param c2_pivots bound for c2_pivots
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
          function(initial_design, n1 = 1, c1f = 0.0, ...) {
              OneStageDesign(n1, c1f)
})


#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("GroupSequentialDesign"),
          function(initial_design, n1 = 1, c1f = -2, c1e = 1.5, n2_pivots = 1, c2_pivots = -2, ...) {
              GroupSequentialDesign(n1,
                                    c1f,
                                    c1e,
                                    n2_pivots,
                                    c2_pivots,
                                    order = length(initial_design@c2_pivots))
})



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("TwoStageDesign"),
          function(initial_design, n1 = 1, c1f = -2, c1e = 1.5, n2_pivots = 1, c2_pivots = -2, ...) {
              TwoStageDesign(n1,
                             c1f,
                             c1e,
                             n2_pivots,
                             c2_pivots,
                             order = length(initial_design@c2_pivots))
})





#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("OneStageDesign"),
          function(initial_design, n1 = NULL, c1f = 5.0, ...) {
              if(is.null(n1))
                  n1 = 5 * initial_design@n1

              return(OneStageDesign(n1, c1f))
})


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("GroupSequentialDesign"),
          function(initial_design, n1 = NULL, c1f = 2.0, c1e = 5.0, n2_pivots = NULL, c2_pivots = 5.0, ...) {
              if(is.null(n1))
                  n1 = 5 * initial_design@n1

              if(is.null(n2_pivots))
                  n2_pivots = 5 * initial_design@n2_pivots

              return(GroupSequentialDesign(n1,
                                           c1f,
                                           c1e,
                                           n2_pivots,
                                           rep(c2_pivots, length(initial_design@c2_pivots))))
})


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("TwoStageDesign"),
          function(initial_design, n1 = NULL, c1f = 2.0, c1e = 5.0, n2_pivots = NULL, c2_pivots = 5.0, ...) {
              if(is.null(n1))
                  n1 = 5 * initial_design@n1

              if(is.null(n2_pivots))
                  n2_pivots = 5 * initial_design@n2_pivots

              return(TwoStageDesign(n1,
                                    c1f,
                                    c1e,
                                    n2_pivots,
                                    rep(c2_pivots, length(initial_design@c2_pivots))))

})
