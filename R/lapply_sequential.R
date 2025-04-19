# lapply_sequential.r

# allow switching between lapply and future_lapply without future.seed issues


#' @title Sequential lapply
#'
#' @description This function allows switching between `lapply` and `future_lapply`
#'     without `future.seed` issues by using `lapply` sequentially.
#'
#' @param X A vector (atomic or list) or an expression object. Other objects
#'     (including classed objects) will be coerced by `as.list`.
#' @param FUN The function to be applied to each element of `X`.
#' @param ... Optional arguments to `FUN`.
#' @param future.seed Ignored in this function, included for compatibility with
#'     `future_lapply`.
#'
#' @return A list of the same length as `X` and named by `X`.
#' @export



lapply.sequential <- function( X, FUN, ..., future.seed = NULL ) {
  lapply( X, FUN, ... )
}
