setClass("IsoperimetricConstraint", representation(score = "Score", orientation = "character", RHS = "numeric"))

IsoperimetricConstraint <- function(score, orientation, RHS) {
        new("IsoperimetricConstraint", score = score, orientation = orientation, RHS = RHS)
}
