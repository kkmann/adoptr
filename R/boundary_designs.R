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
          function(initial_design, n1 = 1, c_buffer = 2, ...) {
              OneStageDesign(n1, min(0, initial_design@c1f - c_buffer))
})


#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("GroupSequentialDesign"),
          function(
              initial_design,
              n1        = 1,
              n2_pivots = 1,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              GroupSequentialDesign(
                  n1,
                  initial_design@c1f - c1_buffer,
                  initial_design@c1e - c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots - c2_buffer,
                  order = length(initial_design@c2_pivots)
              )
})



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("TwoStageDesign"),
          function(
              initial_design,
              n1        = 1,
              n2_pivots = 1,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              TwoStageDesign(
                  n1,
                  initial_design@c1f - c1_buffer,
                  initial_design@c1e - c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots - c2_buffer,
                  order = length(initial_design@c2_pivots)
              )
})





#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("OneStageDesign"),
          function(initial_design, n1_fctr = 5, c_buffer = 2, ...) {
              return(OneStageDesign(n1_fctr * n1(initial_design), initial_design@c1f + c_buffer))
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
