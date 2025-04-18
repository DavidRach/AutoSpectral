# lapply_sequential.r
# allow switching between lapply and future_lapply without future.seed issues

lapply.sequential <- function( X, FUN, ..., future.seed = NULL ) {
  lapply( X, FUN, ... )
}